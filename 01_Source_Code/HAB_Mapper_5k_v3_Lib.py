"""
   Copyright 2023 Ian Housman, RedCastle Resources Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
# Template script to map algal blooms in Google Earth Engine using field observations for training 
# Intended to work within the geeViz package
# python -m pip install geeViz
# Also requires numpy, and pandas
####################################################################################################
import os,sys,threading,json,pandas,time,glob,datetime
import geeViz.getImagesLib as getImagesLib
import geeViz.assetManagerLib as aml
import geeViz.taskManagerLib as tml
import geeViz.gee2Pandas as g2p
import geeViz.foliumView as fv
import numpy as np
ee = getImagesLib.ee
Map = getImagesLib.Map
# Map = fv.foliumMapper()
####################################################################################################
# Function to check for the date of the most recent Sentinel 2 data
def getMostRecentS2Date():
    today = ee.Date(datetime.datetime.today())
    s2s = getImagesLib.getS2(ee.Geometry.Point([-111,41]),
            today.advance(-25,'day'),
            today.advance(2,'day'),
            1,
            365,
            resampleMethod = 'near',
            toaOrSR = 'TOA',
            convertToDailyMosaics = False,
            addCloudProbability = False).sort('system:time_start',False)
    return  int(s2s.first().date().format('DDD').getInfo())
####################################################################################################
# Function to save a Random Forest model from GEE
def saveModel(rfModel,filename,delimiter='split_trees_here'):
    trees = rfModel.explain().getInfo()
    # print('Tree keys:',trees.keys())
    # print(rfModel.getInfo())
    trees=trees['trees']
    o = open(filename,'w')
    o.write(delimiter.join(trees))
    o.close()
# Function to read in saved Random Forest model from GEE
def openModel(filename,delimiter='split_trees_here'):
    o = open(filename,'r')
    trees = o.read().split(delimiter)
    o.close()
    return ee.Classifier.decisionTreeEnsemble(trees)
####################################################################################################    
# Function to get a predictor stack for a supervised classification of water
def getWaterPredictors(studyArea,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],addToMap=False):

    # Set up a standard name ending based on the date params
    nameEnd = 'yr{}-{} jd{}-{}'.format(water_startYear,water_endYear,water_startJulian,water_endJulian)

    # Get the elevation data
    ned = ee.Image("USGS/3DEP/10m").resample('bicubic').rename(['elevation'])
    # ned = ee.Image("USGS/NED").resample('bicubic')
    slope = ee.Terrain.slope(ned).rename(['slope'])
    aspect = ee.Terrain.aspect(ned).rename(['aspect'])

    # Compute some topographic position indices at different scales
    tpi_29 = ned.subtract(ned.focal_mean(14.5)).rename(['tpi_29'])
    tpi_59 = ned.subtract(ned.focal_mean(29.5)).rename(['tpi_59'])

    # Handle the need for a good number of years for TDOM (temporal dark outlier mask cloud shadow masking) needs to work well
    if water_endYear- water_startYear < 3:
        water_startYearT,water_endYearT = water_startYear-2,water_endYear+2
    else:
        water_startYearT,water_endYearT = water_startYear,water_endYear
    
    # Get the S2 imagery
    water_comp =  getImagesLib.getProcessedSentinel2Scenes(studyArea,water_startYearT,water_endYearT,water_startJulian,water_endJulian,resampleMethod='bicubic',convertToDailyMosaics = False,cloudProbThresh=20)

    # Filter the imagery back to the specified years if it had to be extended for TDOM to work
    water_comp = water_comp.filter(ee.Filter.calendarRange(water_startYear,water_endYear,'year'))

    # Get a hillshade
    mean_azimuth = water_comp.reduceColumns(ee.Reducer.mean(),['MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A']).get('mean')
    mean_zenith = water_comp.reduceColumns(ee.Reducer.mean(),['MEAN_INCIDENCE_ZENITH_ANGLE_B8A']).get('mean')
    hillshade = ee.Terrain.hillshade(ned, mean_azimuth, mean_zenith).rename(['hillshade'])

    # Add relevant indices
    water_comp = water_comp.map(getImagesLib.addTCAngles)
    water_comp = water_comp.map(getImagesLib.HoCalcAlgorithm2).median()
    water_mask = getImagesLib.simpleWaterMask(water_comp,0).rename(['water_mask'])

    water_comp = ee.Image.cat([water_comp,ned,slope,aspect,water_mask,hillshade,tpi_29,tpi_59])
    water_comp = water_comp.select(pred_bands)

    # Add layer to map if specified
    if addToMap:
        Map.addLayer(slope.reproject(crs,transform),{'palette':'00F,F00'},'Slope',False)
        Map.addLayer(aspect.reproject(crs,transform),{'palette':'00F,F00'},'Aspect',False)
        Map.addLayer(tpi_29.reproject(crs,transform),{'min':-50,'max':50,'palette':'00F,888,F00'},'TPI 29',False)
        Map.addLayer(tpi_59.reproject(crs,transform),{'min':-50,'max':50,'palette':'00F,888,F00'},'TPI 59',False)

        Map.addLayer(water_comp.reproject(crs,transform),getImagesLib.vizParamsFalse,'Water Comp {}'.format(nameEnd),False)
        Map.addLayer(hillshade,{},'Hillshade',False)
        Map.addLayer(water_mask.selfMask().reproject(crs,transform),{'palette':'00F'},'Simple Water Mask {}'.format(nameEnd),False)

    return water_comp
####################################################################################################
# Function go get a RF supervised classification model of snow/ice, water, and other
# This was developed to avoid commiting water over water-saturated snow and ice that is common in early summer in higher elevations when using the getImagesLib.simpleWaterMask algorithm
# Uses training data manually created in the GEE Playground: https://code.earthengine.google.com/?scriptPath=users%2Faaronkamoske%2FAlgalBlooms%3AWater_Training.js
def getSupervisedWaterModel(water_training_points,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs' ,nTrees=250):
    if not os.path.exists(output_dir):os.makedirs(output_dir)
    # Set up some output filenames
    out_training = os.path.join(output_dir,'Water_Training_yr{}-{}_jd{}-{}.geojson'.format(water_startYear,water_endYear,water_startJulian,water_endJulian))
    out_model = os.path.join(output_dir,'Water_Model_nTrees{}_yr{}-{}_jd{}-{}.txt'.format(nTrees,water_startYear,water_endYear,water_startJulian,water_endJulian))

    # Get the model if it doesn't already exist
    if not os.path.exists(out_model):

        # Read in water training points from geojson to GEE featureCollection
        o = open(water_training_points)
        water_training = json.load(o)
        o.close()
        water_training = ee.FeatureCollection(water_training)

        # Get a predictor stack for the area the points cover
        water_comp = getWaterPredictors(water_training.geometry().bounds(),water_startYear,water_endYear,water_startJulian,water_endJulian,pred_bands)

        # Extract the training data predictor values if they don't already exist
        if not os.path.exists(out_training):
            Map.addLayer(water_training,{'layerType':'geeVector'},'Water Training',True)

            training = water_comp.reduceRegions(water_training, ee.Reducer.first(), None, crs, transform, 4)
            training=training.filter(ee.Filter.notNull(pred_bands)).getInfo()
            o = open(out_training,'w')
            o.write(json.dumps(training))
            o.close()

        # Open the saved water training data
        o = open(out_training)
        water_training = json.load(o)
        o.close()

        # Train the RF model
        rf= ee.Classifier.smileRandomForest(nTrees)
        rf=rf.train(water_training, 'Cls', pred_bands)

        # Get some model info
        getRFModelInfoClassification(rf,'Water_Classification',nTrees,['Snow/Ice','Water','Other'],output_dir)

        # Save the model to a text file
        saveModel(rf,out_model,delimiter='split_trees_here')
        

    # Read in the saved water model
    rf2 = openModel(out_model,delimiter='split_trees_here')
    return rf2
####################################################################################################
# Function to get a saved water RF model
def getCachedWaterModel(output_dir,model_name=None):
  if model_name == None:
    model_name = glob.glob(os.path.join(output_dir),'Water_Model_nTrees*.txt')[0]
  print('Opening cached water model:',model_name)
  o = open(model_name,'r')
  trees = o.read().split('split_trees_here')
  o.close()
  rf2 = ee.Classifier.decisionTreeEnsemble(trees)
  return rf2
####################################################################################################
# Function to apply a water classification model
def applyWaterClassification(output_dir,model,studyArea,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],addToMap=False,water_contract_pixels=1):

    # Get a saved model if one is not provided
    if model == None:
        model = getCachedWaterModel(output_dir)

    # Set up a common ending name based on the date params
    nameEnd = 'yr{}-{} jd{}-{}'.format(water_startYear,water_endYear,water_startJulian,water_endJulian)

    # Get the water predictors for the specified dates
    water_comp = getWaterPredictors(studyArea,water_startYear,water_endYear,water_startJulian,water_endJulian,pred_bands,addToMap)

    # Apply the model
    predicted = water_comp.classify(model)
    predicted = predicted.rename(['SnowIce_Water'])

    # Pull out just the water class and contract it by n pixels
    water = predicted.eq(2).focal_min(water_contract_pixels).rename(['Water'])
    waterLegendDict = {'Snow/Ice':'0FF','Water':'00F'}

    # Add output to the map if specified
    if addToMap:
        Map.addLayer(predicted.updateMask(predicted.lte(2)).reproject(crs,transform),{'min':1,'max':2,'palette':'0FF,00F','layerType':'geeImage','classLegendDict':waterLegendDict},'Water Classification {}'.format(nameEnd))
    
    # Return the composite along with the classified output and just the water mask
    return water_comp.addBands(predicted).addBands(water)
####################################################################################################
# Function to get a clean waterbody sample
# Takes clean waterbody centroid points provided by Paul Marone, filters out ones that are likely not actually clean, gets a late summer composite, finds the areas covered in water, and derives a sample over those water bodies
# Current workflow does convert the shapefile provided by Paul to geojson using ogr2ogr (not included in this script)
def prepareCleanTrainingData(clean_points,clean_startYear = 2018,clean_endYear = 2022,clean_startJulian = 213,clean_endJulian = 229,nSamples = 150,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI','elevation'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',water_model=None):
    if not os.path.exists(output_dir):os.makedirs(output_dir)
    # Load the clean waterbody centroid points
    o = open(clean_points)
    clean_points = ee.FeatureCollection(json.load(o))
    o.close()

    # Filter out points likely to not be clean
    clean_points = ee.FeatureCollection(clean_points).filter(ee.Filter.stringContains('Control','Bloom Recorded').Not())
    clean_points = ee.FeatureCollection(clean_points).filter(ee.Filter.stringContains('Control','Bloom Whitnessed').Not())
    clean_points = ee.FeatureCollection(clean_points).filter(ee.Filter.stringContains('Control','No Bloom Seen').Not())
    clean_points = ee.FeatureCollection(clean_points).filter(ee.Filter.stringContains('Control','Sampled No Bloom').Not())

    # Add the centroid points to the map
    Map.addLayer(clean_points.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'Clean Water Body Centroids')

    # Get the bounds of all the points for a study area
    clean_studyArea = clean_points.geometry().bounds(500,'EPSG:5070').buffer(2000)

    # Get the water mask and composite 
    water_comp = applyWaterClassification(output_dir,water_model,clean_studyArea,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,crs,transform,water_pred_bands,addToMap=False,water_contract_pixels=2)

    # Extract the wtaer mask and mask out the composite to its extent
    water = water_comp.select(['Water']).selfMask()
    compWater = water_comp.select(pred_bands).updateMask(water)
 
    # Add the composite to the map
    Map.addLayer(water_comp.reproject(crs,transform),getImagesLib.vizParamsFalse,'Clean Comp {}-{}'.format(clean_startYear,clean_endYear),False)
    Map.addLayer(water.reproject(crs,transform),{'palette':'00F','classLegendDict':{'Water':'00F'}},'Clean Water Mask {}-{}'.format(clean_startYear,clean_endYear),False)

    oldWaterMask = getImagesLib.simpleWaterMask(water_comp).selfMask()
                                                
    Map.addLayer(oldWaterMask,{'palette':'0FF','classLegendDict':{'Water':'0FF'}},'Old Water Mask',False)

    # Get the ids of each point to then extract
    ids  = clean_points.aggregate_histogram('system:index').keys().getInfo()
    
    # Extract each water body by centroid id
    for id in ids:
        # Set up the output filename
        of = os.path.join(output_dir,'Clean_Training_Samples_id{}_n{}_yr{}-{}_jd{}-{}.geojson'.format(id,nSamples,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian))

        if not os.path.exists(of):
            print('Extracting clean data for point id:',id)

            # Get the point for the id
            clean_pt = clean_points.filter(ee.Filter.eq('system:index',id))

            # Get a generous study area aroudn the point
            clean_point_b = clean_pt.geometry().bounds(500,'EPSG:5070').buffer(2000)
            
            # Convert the water mask to vector and then filter it down to the one the point falls within
            water_v = water.clip(clean_point_b).reduceToVectors(scale = 30).filterBounds(clean_pt)

            # Create a sample within that waterbody outline
            pts = ee.FeatureCollection.randomPoints(water_v, nSamples, 1, 50)

            # Extract the predictor values
            training = compWater.reduceRegions(pts, ee.Reducer.first(), None, crs, transform, 4)
            training=training.filter(ee.Filter.notNull(pred_bands)).getInfo()

            # Write the output
            o = open(of,'w')
            o.write(json.dumps(training))
            o.close()
####################################################################################################
# Function to extract the HCB predictor values for each unique sample collection date
# Gets a composite for a +- 2 week window and extracts those values
def extractHCBPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,water_model,d,output_dir,focal_radius=None,weeks_before=2,weeks_after=2):
    # Get date window 
    startDate = ee.Date(d).advance(-weeks_before,'week')
    endDate = ee.Date(d).advance(weeks_after,'week')
    startYear = int(startDate.get('year').getInfo())
    endYear = int(endDate.get('year').getInfo())
    startJulian = int(startDate.format('DDD').getInfo())
    endJulian = int(endDate.format('DDD').getInfo())

    # Filter HCB data to the provided date
    hcb_d = hcb_fc.filter(ee.Filter.eq(dateProp,d))

    # Get a water mask and composite for those dates
    water_comp = applyWaterClassification(output_dir,water_model,hcb_d,startYear,endYear,startJulian,endJulian,crs,transform,addToMap=False)
    water = water_comp.select(['Water']).selfMask()
    compWater = water_comp.select(pred_bands).updateMask(water)
    
    # Extract the values for a 5x5 pixel window mean around the point (this is necessary since many samples are right on the shorline)
    if focal_radius != None:
        compWater = compWater.focalMean(focal_radius)
    t = compWater.reduceRegions(hcb_d, ee.Reducer.first(), None, crs, transform, 4)
    return t
####################################################################################################
# Function to write out the predictor table from GEE to geojson
def writePredictorTables(featureCollection,filename):
    out_f = ee.FeatureCollection(featureCollection).getInfo()
    o = open(filename,'w')
    o.write(json.dumps(out_f))
    o.close()
####################################################################################################
# Function to extract HCB samples by sets
# Originally written to extract using multiple threads/processes, but that approach doesn't work. This now extracts by groups of samples
def batchExtractHCBPredictorTables(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,output_dir,water_model,nameStart='HCB',focal_radius=None):
    if not os.path.exists(output_dir):os.makedirs(output_dir)

    # Filter training data dates to when Sentinel 2 is available (roughly 2017 although there are data for 2016 as well)
    hcb_fc = hcb_fc.map(lambda f:f.set('system:time_start',ee.Date(f.get(dateProp)).millis()))
    hcb_fc = hcb_fc.filterDate('2017-01-01','2200-01-01')

    # Get the dates and group them into 20 groups
    dates = hcb_fc.aggregate_histogram(dateProp).keys().getInfo()

    date_sets = new_set_maker(dates,20)

    # Iterate across each group and extract the predictor tables
    date_set_i = 1
    for date_set in date_sets:
        of = os.path.join(output_dir,'{}_Training_{}_foc{}.geojson'.format(nameStart,date_set_i,focal_radius))
        if not os.path.exists(of):
            print('Extracting dates:',date_set)
            out_f = []
            for d in date_set:
                t = extractHCBPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,water_model,d,output_dir,focal_radius)
                out_f.append(t)
            out_f = ee.FeatureCollection(out_f).flatten()
            writePredictorTables(out_f,of)
        date_set_i +=1
####################################################################################################
# Function to get model info for a RF regression model
# Includes the variable importance and out of bag error (bottom row)
def getRFModelInfoRegression(rf,name,nTrees,table_dir):
    out_importance_filename = os.path.join(table_dir,'{}_nTrees{}_Var_Importance.csv'.format(name,nTrees))
    if not os.path.exists(out_importance_filename):
        # Get the importance and OOB error
        exp = ee.Dictionary({'importance':rf.explain().get('importance'),\
                                'oobError':rf.explain().get('outOfBagErrorEstimate')}).getInfo()

        obError = exp['oobError']

        # Convert the importance into a comma delimited set of rows
        importance = list(zip(list(exp['importance'].keys()),list(exp['importance'].values())))
        importance = sorted(importance, key=lambda row: row[1]*-1)
        importance_lines = 'Name,Importance\n'
        for line in importance:
            importance_lines+=','.join([str(i) for i in line])+'\n'

        # Add in the OOB Error
        importance_lines+='OOB Error,{}\n'.format(obError)
        
        # Write output csv
        o = open(out_importance_filename,'w')
        o.write(importance_lines)
        o.close()
####################################################################################################
# Function to get RF classification model error matrix, oob accuracy, and variable importance
def getRFModelInfoClassification(rf,name,nTrees,class_names,table_dir):

    # Set up output filenames
    out_table_filename = os.path.join(table_dir,'{}_nTrees{}_OOB_Error_Matrix.csv'.format(name,nTrees))
    out_importance_filename = os.path.join(table_dir,'{}_nTrees{}_Var_Importance.csv'.format(name,nTrees))

    if not os.path.exists(out_table_filename) or not os.path.exists(out_importance_filename):

        # Get the confusion matrix, kappa, and OOB error
        cm = rf.confusionMatrix()
        exp = ee.Dictionary({'importance':rf.explain().get('importance'),\
                                'kappa':cm.kappa(),\
                                'oobError':rf.explain().get('outOfBagErrorEstimate'),
                                'cm':cm}).getInfo()
        oobAcc = 1-exp['oobError']
        kappa = exp['kappa']
        cm = exp['cm']

        # Convert the importance to comma delimited lines
        importance = list(zip(list(exp['importance'].keys()),list(exp['importance'].values())))
        importance = sorted(importance, key=lambda row: row[1]*-1)
        importance_lines = 'Name,Importance\n'
        for line in importance:
            importance_lines+=','.join([str(i) for i in line])+'\n'
        
        # Convert the confusion matrix to comma delimtied and compute commission and omission error
        cm = cm[1:]
        cm = [row[1:] for row in cm]

        allSum = np.sum(cm)
        diag = []
        rowAccs = []
        colAccs = []

        for ri in range(0,len(cm)):
            row =cm[ri]
            col = [cm[i][ri] for i in range(0,len(cm))]
            v = row[ri]
            diag.append(v)
            rowSum = np.sum(row)
            colSum = np.sum(col)
            rowAcc = v/rowSum
            colAcc = v/colSum
            rowAccs.append(rowAcc)
            colAccs.append(colAcc)

            
        overall_acc = np.sum(diag)/allSum

        out_lines = ',{},Commission Error (1-n)\n'.format(','.join(class_names))
        for ri in range(0,len(cm)):
            out_lines+='{},{},{}\n'.format(class_names[ri],','.join([str(i) for i in cm[ri]]),rowAccs[ri])
        out_lines+='Omission Error (1-n),{},,\n'.format(','.join([str(i) for i in colAccs]))
        out_lines+='Overall Accuracy,{}\n'.format(overall_acc)
        out_lines+='OOB Accuracy,{}\n'.format(oobAcc)
        out_lines+='Kappa,{}\n'.format(overall_acc)
       

        # Write output files
        o = open(out_table_filename,'w')
        o.write(out_lines)
        o.close()

        o = open(out_importance_filename,'w')
        o.write(importance_lines)
        o.close()
####################################################################################################
# Load multiple geojson files as a single GEE featureCollection
def loadGeoJSON(files):
    out_fc = []
    for t in files:
        o = open(t,'r')
        out_fc.append(ee.FeatureCollection(json.load(o)))
        o.close()
    return ee.FeatureCollection(out_fc).flatten()
####################################################################################################
# Function to get the RF algal models
def getAlgalModels(table_dir,pred_bands,nTrees=100,getError=False,hcb_focal_radius=2.5,wdeq_focal_radius=None,apply_HCB_Model=True,apply_WDEQ_Model=True):

    # Name of field for algal classes that must be added to the training data
    training_field = 'HAB'

    # Saved model names
    class_model = os.path.join(table_dir,'HAB_Classification_RF_Model.txt')
    count_model = os.path.join(table_dir,'Cell_Count_RF_Model.txt')
    biovolume_model = os.path.join(table_dir,'BioVolume_RF_Model.txt')

    # Get the RF models for classification, count, and biovolume
    def getModel(training,training_field,outputMode):
        rf= ee.Classifier.smileRandomForest(nTrees)
        rf = rf.setOutputMode(outputMode)
        rf=rf.train(training, training_field, pred_bands)
        return rf

    # Set up some class names and numbers for classification
    class_names=['Not HAB','HAB']
    class_numbers = [1,2]

    # Read in HCB and clean training geojson files as GEE featureCollections
    if apply_HCB_Model:
        hcb_training = glob.glob(os.path.join(table_dir,'HCB*_foc{}*.geojson'.format(hcb_focal_radius)))
        clean_training = glob.glob(os.path.join(table_dir,'Clean_*.geojson'))
        hcb_training_fc = loadGeoJSON(hcb_training)
        clean_training_fc = loadGeoJSON(clean_training)
        # Set the clean training dependent variable fields to the respective clean values
        # 1 for classification, and 0 for both regression variables (Count and Biovolume)
        clean_training_fc = clean_training_fc.map(lambda f:ee.Feature(f).set(training_field,1))
        clean_training_fc = clean_training_fc.map(lambda f:ee.Feature(f).set('Cyanobacteria Count (cells/mL)',0))
        clean_training_fc = clean_training_fc.map(lambda f:ee.Feature(f).set('Cyanobacteria Biovolume (um^3)',0))
        # print(clean_training_fc.first().toDictionary().keys().getInfo())
        # print('N HCB Points:',hcb_training_fc.size().getInfo())
        # print('N Clean Points:',clean_training_fc.size().getInfo())
        # Cyanobacteria Count (cells/mL) > 25000 for threshold
        # Cyanobacteria Biovolume (um^3)
        # If failing bin taxa

        # Filter out HCB and WDEQ data to ensure there are observations and the count is > the WY threshold of 25000
        hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(pred_bands))
        
        hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(['Cyanobacteria Count (cells/mL)']))
        hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(['Cyanobacteria Biovolume (um^3)']))
        hcb_training_fc = hcb_training_fc.filter(ee.Filter.hasType('Cyanobacteria Count (cells/mL)','float'))
        hcb_training_fc = hcb_training_fc.filter(ee.Filter.hasType('Cyanobacteria Biovolume (um^3)','float'))
        hcb_training_fc = hcb_training_fc.filter(ee.Filter.gte('Cyanobacteria Count (cells/mL)',25000))
        
        # Set the HCB training dependent variable field to the respective algal value
        hcb_training_fc = hcb_training_fc.map(lambda f:ee.Feature(f).set(training_field,2))

        Map.addLayer(hcb_training_fc.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector','strokeColor':'F00'},'HCB Training')

        
        Map.addLayer(clean_training_fc.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'Clean Training')

         # Combine the HCB and clean training data
        training = hcb_training_fc.merge(clean_training_fc)
        Map.addLayer(training.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'All Training')
        
        rf_hab = getModel(training,training_field,'CLASSIFICATION')
        rf_cells=getModel(training,'Cyanobacteria Count (cells/mL)','REGRESSION')
        rf_vol=getModel(training,'Cyanobacteria Biovolume (um^3)','REGRESSION')

        # Save models (can only save classification model so all models are retrained each time currently)
        saveModel(rf_hab,class_model)
        # saveModel(rf_cells,count_model)
        # saveModel(rf_vol,biovolume_model)
        # Get error info
        if getError:
            getRFModelInfoClassification(rf_hab,'HAB_Classification',nTrees,class_names,table_dir)
            getRFModelInfoRegression(rf_cells,'Cell_Count',nTrees,table_dir)
            getRFModelInfoRegression(rf_vol,'Biovolume',nTrees,table_dir)
    if apply_WDEQ_Model:
        wdeq_training = glob.glob(os.path.join(table_dir,'WDEQ*_foc{}*.geojson'.format(wdeq_focal_radius)))
        wdeq_training_fc = loadGeoJSON(wdeq_training)
    
        wdeq_training_fc=wdeq_training_fc.filter(ee.Filter.notNull(pred_bands))
        # wdeq_training_fc  = wdeq_training_fc.filter(ee.Filter.eq('Class','Cyanophyceae'))
        Map.addLayer(wdeq_training_fc.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector','strokeColor':'FF0'},'WDEQ Training')
    
        rf_wdeq_density=getModel(wdeq_training_fc,'Density (cells/L)','REGRESSION')
        rf_wdeq_count=getModel(wdeq_training_fc,'Individuals (Raw Cnt)','REGRESSION')

        # Get error info
        if getError:
            getRFModelInfoRegression(rf_wdeq_count,'WDEQ_Count',nTrees,table_dir)
            getRFModelInfoRegression(rf_wdeq_density,'WDEQ_Density',nTrees,table_dir)
    
    # # Return the models
    if not apply_WDEQ_Model:rf_wdeq_density,rf_wdeq_count=None,None
    if not apply_HCB_Model:rf_hab,rf_cells,rf_vol=None,None,None

    return rf_hab,rf_cells,rf_vol,rf_wdeq_density,rf_wdeq_count,class_names,class_numbers
####################################################################################################
# Function to apply an algal model created in the functions above
def applyAlgalModels(rf_hab,rf_cells,rf_vol,rf_wdeq_density,rf_wdeq_count,water_model,class_names,class_numbers,applyStudyArea,applyYear = 2022,applyStartJulian = 190, applyEndJulian = 210,pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI','elevation'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',export_outputs=False,output_asset_collection='projects/gtac-algal-blooms/assets/outputs/HAB-RF-Images',studyAreaName='WY',apply_HCB_Model=True,apply_WDEQ_Model=True):

    # Get a common name ending for the dates provided
    applyStartDate = ee.Date.fromYMD(applyYear,1,1).advance(applyStartJulian,'day')
    
    applyEndDate = ee.Date.fromYMD(applyYear,1,1).advance(applyEndJulian,'day')
    
    nameEnd = 'yr{}_{}_{}'.format(applyYear,applyStartDate.format('MM-dd').getInfo(),applyEndDate.format('MM-dd').getInfo())

    outputAssetName = 'RF-Algal-Stack_{}_{}'.format(studyAreaName,nameEnd)
    
    # Get the water mask and composite for the dates
    water_comp = applyWaterClassification(output_dir,water_model,applyStudyArea,applyYear,applyYear,applyStartJulian,applyEndJulian,crs,transform,water_pred_bands,addToMap=False)
    water = water_comp.select(['Water']).selfMask()
    comp = water_comp.select(pred_bands)

    # Mask composite to extent of water
    compWater = comp.updateMask(water)

    # Add layers to map
    Map.addLayer(applyStudyArea,{},'Apply Study Area',False)
    Map.addLayer(comp.reproject(crs,transform),getImagesLib.vizParamsFalse,'Median S2 {}'.format(nameEnd),False)
    Map.addLayer(water.selfMask().reproject(crs,transform),{'palette':'00D'},'Median S2 Water {}'.format(nameEnd),False)

    # Apply classification model
    predicted = compWater.classify(rf_hab)
    algalLegendDict={'Algal Negative':'00D','Algal Positive':'D00'}
    # Map.addLayer(predicted.reproject(crs,transform),{'min':class_numbers[0],'max':class_numbers[-1],'palette':'00D,D00','layerType':'geeImage','classLegendDict':algalLegendDict},'AB Classified {}'.format(nameEnd),False)

    def getMapStretch(predicted,nameStart,map_stretch_pctls = [5,95]):
        stats_filename = os.path.join(output_dir,'Algal_{}_Stats_{}-{}_{}.txt'.format(nameStart,map_stretch_pctls[0],map_stretch_pctls[1],nameEnd))
        if not os.path.exists(stats_filename):
            print('Computing output image stats:',stats_filename)
            predicted_stats= [str(int(i)) for i in predicted.reduceRegion(ee.Reducer.percentile(map_stretch_pctls), applyStudyArea, 3000,tileScale=4).values().getInfo()]
            o = open(stats_filename,'w')
            o.write(','.join(predicted_stats))
            o.close()
        o = open(stats_filename,'r')
        stats = o.read().split(',')
        o.close()
        return [int(i) for i in stats]
    
    if apply_HCB_Model:
        # Apply count regression model
        cell_predicted = compWater.classify(rf_cells)
        # cell_predicted_stats = getMapStretch(cell_predicted,'Count')
        Map.addLayer(cell_predicted,{'min':1000000,'max':5000000,'palette':'00D,D00','layerType':'geeImage'},'Cyanobacteria Count  {}'.format(nameEnd),False)

        # Apply biovolume regression model
        biovolume_predicted = compWater.classify(rf_vol)
        # biovolume_predicted_stats = getMapStretch(biovolume_predicted,'Biovolume')
        Map.addLayer(biovolume_predicted,{'min':200000000,'max':1000000000,'palette':'00D,D00','layerType':'geeImage'},'Cyanobacteria Biovolume  {}'.format(nameEnd),False)

    if apply_WDEQ_Model:
        # Apply WDEQ density regression model
        wdeq_density_predicted = compWater.classify(rf_wdeq_density)
        # wdeq_density_stats = getMapStretch(wdeq_density_predicted,'WDEQ_Density')
        Map.addLayer(wdeq_density_predicted.divide(1000),{'min':0,'max':1000000000,'palette':'00D,D00','layerType':'geeImage'},'WDEQ Density  {}'.format(nameEnd),False)

        # Apply WDEQ count regression model
        wdeq_count_predicted = compWater.classify(rf_wdeq_count)
        # wdeq_count_stats = getMapStretch(wdeq_count_predicted,'WDEQ_Count')
        Map.addLayer(wdeq_count_predicted,{'min':0,'max':1000,'palette':'00D,D00','layerType':'geeImage'},'WDEQ Count {}'.format(nameEnd),False)
    if export_outputs:
        if not aml.ee_asset_exists(output_asset_collection + '/'+ outputAssetName):
            # Set up outputs for export
            out_pred_stack = []
            out_pred_stack_names = []
            if apply_HCB_Model:
                out_pred_stack.extend([cell_predicted.uint32(),biovolume_predicted.uint32()])
                out_pred_stack_names.extend(['HCB_Count','HCB_Biovolume'])
            if apply_WDEQ_Model:
                out_pred_stack.extend([wdeq_count_predicted.uint32(),wdeq_density_predicted.uint32()])
                out_pred_stack_names.extend(['WDEQ_Count','WDEQ_Density'])
            out_pred_stack = ee.Image.cat(out_pred_stack).rename(out_pred_stack_names).set({'system:time_start':applyStartDate.millis(),
                                        'system:time_end':applyEndDate.millis(),
                                        'year':applyYear,
                                        'startJulian':applyStartJulian,
                                        'endJulian':applyEndJulian,
                                        'studyAreaName':studyAreaName})
            # task_list.append(outputAssetName)
            t = ee.batch.Export.image.toAsset(out_pred_stack.clip(applyStudyArea), 
                            description = outputAssetName, 
                            assetId = output_asset_collection + '/'+ outputAssetName, 
                            pyramidingPolicy = {'AB':'mode','HCB_Count':'mean','HCB_Count':'mean','WDEQ_Count':'mean','WDEQ_Density':'mean'}, 
                            dimensions = None, 
                            region = None, 
                            scale = None, 
                            crs = crs, 
                            crsTransform = transform, 
                            maxPixels = 1e13)
            # print('Exporting:',task_list)
            print(t)
            t.start()
####################################################################################################
# Divide list into a list of lists length n threads
def new_set_maker(in_list,threads):
    out_sets =[]
    for t in range(threads):
        out_sets.append([])
    i =0
    for il in in_list:

        out_sets[i].append(il)
        i += 1
        if i >= threads:
            i = 0
    return out_sets
####################################################################################################
def limitThreads(limit):
  while threading.active_count() > limit:
    time.sleep(2)
    print(threading.active_count(),'threads running')
####################################################################################################
# Wrapper function for it all
def supervised_algal_mapper(water_training_points,hcb_data,wdeq_data,clean_points,\
clean_startYear = 2018,clean_endYear = 2022,clean_startJulian = 213,clean_endJulian = 229,clean_nSamples = 150,\
hcb_lat='Sampling Latitude', hcb_lon='Sampling Longitude',hcb_dateProp = 'Sample Date',hcb_properties=[],\
    wdeq_lat='Sampling Latitude', wdeq_lon='Sampling Longitude',wdeq_dateProp = 'Sample Date',wdeq_properties=[],\
        crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI','elevation'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',nTrees=250,getError=False,applyStudyArea=None,applyYears = [2022],applyStartJulians = range(150,200,14),applyNDayWindow = 28,export_outputs=False,output_asset_collection='projects/gtac-algal-blooms/assets/outputs/HAB-RF-Images',studyAreaName='WY',hcb_focal_radius = 2.5,wdeq_focal_radius=None,apply_HCB_Model=True,apply_WDEQ_Model=True):

    # Get the water RF model
    water_model = getSupervisedWaterModel(water_training_points,water_startYear=clean_startYear,water_endYear=clean_endYear,water_startJulian=clean_startJulian,water_endJulian=clean_endJulian,crs = crs,transform = transform,pred_bands=water_pred_bands,output_dir =output_dir,nTrees=nTrees)
    # applyWaterClassification(output_dir,water_model,applyStudyArea,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,crs,transform,water_pred_bands,addToMap=True)
    
    if apply_HCB_Model:
        # Get clean training data
        prepareCleanTrainingData(clean_points,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,clean_nSamples,crs,transform ,pred_bands,water_pred_bands,output_dir,water_model)

        # Get HCB and WDEQ training data
        hcb_fc =  g2p.tableToFeatureCollection(hcb_data,hcb_lat,hcb_lon,hcb_properties) 

        batchExtractHCBPredictorTables(hcb_fc,hcb_dateProp,crs,transform,pred_bands,water_pred_bands,output_dir,water_model,'HCB',hcb_focal_radius)

    if apply_WDEQ_Model:
        wdeq_fc = g2p.tableToFeatureCollection(wdeq_data,wdeq_lat,wdeq_lon,wdeq_properties,wdeq_dateProp,[wdeq_lat,wdeq_lon,wdeq_dateProp,'Class'])
        wdeq_fc  = wdeq_fc.filter(ee.Filter.eq('Class','Cyanophyceae'))
        batchExtractHCBPredictorTables(wdeq_fc,wdeq_dateProp,crs,transform,pred_bands,water_pred_bands,output_dir,water_model,'WDEQ',wdeq_focal_radius)
    
    

    # Get algal RF models
    rf_hab,rf_cells,rf_vol,rf_wdeq_density,rf_wdeq_count,class_names,class_numbers = getAlgalModels(output_dir,pred_bands,nTrees,getError,hcb_focal_radius,wdeq_focal_radius,apply_HCB_Model,apply_WDEQ_Model)

    # Apply algal models
    # if export_outputs:
    for applyYear in applyYears:
        for applyStartJulian in applyStartJulians:
            applyEndJulian = applyStartJulian+applyNDayWindow
            applyAlgalModels(rf_hab,rf_cells,rf_vol,rf_wdeq_density,rf_wdeq_count,water_model,class_names,class_numbers,applyStudyArea,applyYear,applyStartJulian, applyEndJulian,pred_bands,water_pred_bands,crs,transform,output_dir,export_outputs,output_asset_collection,studyAreaName,apply_HCB_Model,apply_WDEQ_Model)
####################################################################################################
# Function to view outputs (this function is deprecated by the Bloom-Mapper)
# Dev: https://dev.wrk.fs.usda.gov/forest-atlas/lcms-viewer/bloom-mapper.html 
# Prod: https://apps.fs.usda.gov/lcms-viewer/bloom-mapper.html
def viewExportedOutputs(output_dir,output_asset_collection):
    ab = ee.ImageCollection(output_asset_collection)
    algalLegendDict={'Algal Negative':'00D','Algal Positive':'D00'}
    Map.addTimeLapse(ab.select([0]),{'min':1,'max':2,'palette':'00D,D00','classLegendDict':algalLegendDict,'dateFormat':'YYMMdd','advanceInterval':'day'},'AB Classified')

    def getStats(name_ending):
        stats_files = glob.glob(os.path.join(output_dir,name_ending))
        stats = []
        for f in stats_files:
            o = open(f,'r')
            l = [int(n) for n in o.read().split(',')]
            stats.append(l)
            o.close()
        stats = np.array(stats)
        return np.mean(stats,0)
    
    count_stats = getStats('Algal_Count*')
    Map.addTimeLapse(ab.select([1]),{'min':count_stats[0],'max':count_stats[1],'palette':'00D,D00','dateFormat':'YYMMdd','advanceInterval':'day'},'Cyanobacteria Count (cells/mL)')

    biovolume_stats = getStats('Algal_Biovolume*')
    Map.addTimeLapse(ab.select([2]),{'min':biovolume_stats[0],'max':biovolume_stats[1],'palette':'00D,D00','dateFormat':'YYMMdd','advanceInterval':'day'},'Cyanobacteria Biovolume (um3)')
    
####################################################################################################