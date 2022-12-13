"""
   Copyright 2022 Ian Housman, RedCastle Resources Inc.

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

# Script to map algal blooms in Google Earth Engine using field observations for training 
# Intended to work within the geeViz package
# python -m pip install geeViz
# Also requires numpy
####################################################################################################
import os,sys,threading,json,pandas,time,glob
import geeViz.getImagesLib as getImagesLib
import numpy as np
ee = getImagesLib.ee
Map = getImagesLib.Map
####################################################################################################
# Function to convert a Pandas dataframe to geojson
# Function taken from: https://notebook.community/captainsafia/nteract/applications/desktop/example-notebooks/pandas-to-geojson
def df_to_geojson(df, properties, lat='latitude', lon='longitude'):
    # create a new python dict to contain our geojson data, using geojson format
    geojson = {'type':'FeatureCollection', 'features':[]}

    # loop through each row in the dataframe and convert each row to geojson format
    for _, row in df.iterrows():
      if not pandas.isnull(row[lon]) and not pandas.isnull(row[lat]):
        # create a feature template to fill in
        feature = {'type':'Feature',
                    'properties':{},
                    'geometry':{'type':'Point',
                                'coordinates':[]}}

        # fill in the coordinates
        feature['geometry']['coordinates'] = [row[lon],row[lat]]

        # for each column, get the value and add it as a new feature property
        for prop in properties:
          p = row[prop]
          if pandas.isnull(p):p= 'NA'
          feature['properties'][prop] = p
        
        # add this feature (aka, converted dataframe row) to the list of features inside our dict
        
        geojson['features'].append(feature)
    
    return geojson
####################################################################################################
# Function to take the Excel HCB spreadsheet and convert it to a GEE featureCollection
def prepareTrainingData(hcb_data,lat='Sampling Latitude', lon='Sampling Longitude',properties=[]):
    # Read in the Excel table as a Pandas dataframe
    hcb_df = pandas.read_excel(hcb_data)

    # Convert the time to a user-friendly format of yyyy-mm-DD
    for c in hcb_df.columns[hcb_df.dtypes=='datetime64[ns]']:
        hcb_df[c]= hcb_df[c].dt.strftime('%Y-%m-%d')

    # Convert the dataframe to geojson
    hfb_json = df_to_geojson(hcb_df, properties, lat, lon)

    # Read in the geojson as a GEE featureCollection
    hfb_json = ee.FeatureCollection(hfb_json)
    return hfb_json
####################################################################################################
# Function to get a predictor stack for a supervised classification of water
def getWaterPredictors(studyArea,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],addToMap=False):

    # Set up a standard name ending based on the date params
    nameEnd = 'yr{}-{} jd{}-{}'.format(water_startYear,water_endYear,water_startJulian,water_endJulian)

    # Get the elevation data
    ned = ee.Image("USGS/NED").resample('bicubic')
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

    water_comp = ee.Image.cat([water_comp,slope,aspect,water_mask,hillshade,tpi_29,tpi_59])
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
def getSupervisedWaterModel(water_training_points,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs' ,nTrees=250):

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
        trees = rf.explain().getInfo()['trees']
        o = open(out_model,'w')
        o.write('split_trees_here'.join(trees))
        o.close()

    # Read in the saved water model
    o = open(out_model,'r')
    trees = o.read().split('split_trees_here')
    o.close()
    rf2 = ee.Classifier.decisionTreeEnsemble(trees)
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
def applyWaterClassification(output_dir,model,studyArea,water_startYear=2018,water_endYear=2022,water_startJulian=213,water_endJulian=229,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59'],addToMap=False,water_contract_pixels=1):

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
def prepareCleanTrainingData(clean_points,clean_startYear = 2018,clean_endYear = 2022,clean_startJulian = 213,clean_endJulian = 229,nSamples = 150,crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',water_model=None):

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
    Map.addLayer(clean_points.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'Clean Training')

    # Get the bounds of all the points for a study area
    clean_studyArea = clean_points.geometry().bounds(500,'EPSG:5070').buffer(2000)

    # Get the water mask and composite 
    water_comp = applyWaterClassification(output_dir,water_model,clean_studyArea,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,crs,transform,water_pred_bands,addToMap=False,water_contract_pixels=3)

    # Extract the wtaer mask and mask out the composite to its extent
    water = water_comp.select(['Water']).selfMask()
    compWater = water_comp.select(pred_bands).updateMask(water)
 
    # Add the composite to the map
    Map.addLayer(water_comp,getImagesLib.vizParamsFalse,'Clean Comp {}-{}'.format(clean_startYear,clean_endYear),False)

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
def extractHCBPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,water_model,d):
    # Get date window 
    startDate = ee.Date(d).advance(-2,'week')
    endDate = ee.Date(d).advance(2,'week')
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
    
    # Extract the values for a 9x9 pixel window mean around the point (this is necessary since many samples are right on the shorline)
    t = compWater.focalMean(4.5).reduceRegions(hcb_d, ee.Reducer.first(), None, crs, transform, 4)
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
def batchExtractHCBPredictorTables(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,output_dir,water_model):
    if not os.path.exists(output_dir):os.makedirs(output_dir)

    # Get the dates and group them into 20 groups
    dates = hcb_fc.aggregate_histogram(dateProp).keys().getInfo()
    date_sets = new_set_maker(dates,20)

    # Iterate across each group and extract the predictor tables
    date_set_i = 1
    for date_set in date_sets:
        of = os.path.join(output_dir,'HCB_Training_{}.geojson'.format(date_set_i))
        if not os.path.exists(of):
            out_f = []
            for d in date_set:
                t = extractHCBPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,water_model,d)
                out_f.append(t)
                out_f = ee.FeatureCollection(out_f).flatten()
            writePredictorTables(out_f,of)
        date_set_i +=1
####################################################################################################
# Function to get model info for a RF regression model
# Includes the variable importance and out of bag error (bottom row)
def getRFModelInfoRegression(rf,name,nTrees,table_dir):
    out_importance_filename = os.path.join(table_dir,'{}_{}_Var_Importance.csv'.format(name,nTrees))
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
def getAlgalModels(table_dir,training_field,pred_bands,nTrees=100,outputMode='CLASSIFICATION',getError=False):

    training_field = 'HAB'

    # Read in HCB and clean training geojson files as GEE featureCollections
    hcb_training = glob.glob(os.path.join(table_dir,'HCB*.geojson'))
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

    # Filter out HCB data to ensure there are observations and the count is > the WY threshold of 25000
    hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(pred_bands))
    hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(['Cyanobacteria Count (cells/mL)']))
    hcb_training_fc=hcb_training_fc.filter(ee.Filter.notNull(['Cyanobacteria Biovolume (um^3)']))
    hcb_training_fc = hcb_training_fc.filter(ee.Filter.hasType('Cyanobacteria Count (cells/mL)','float'))
    hcb_training_fc = hcb_training_fc.filter(ee.Filter.hasType('Cyanobacteria Biovolume (um^3)','float'))
    hcb_training_fc = hcb_training_fc.filter(ee.Filter.gte('Cyanobacteria Count (cells/mL)',25000))

    # Set the HCB training dependent variable field to the respective algal value
    hcb_training_fc = hcb_training_fc.map(lambda f:ee.Feature(f).set(training_field,2))
    
    # Set up some class names and numbers for classification
    class_names=['Not HAB','HAB']
    class_numbers = [1,2]

    Map.addLayer(hcb_training_fc.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'HCB Training')
    Map.addLayer(clean_training_fc.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'Clean Training')

    # Combine the HCB and clean training data
    training = hcb_training_fc.merge(clean_training_fc)
    Map.addLayer(training.map(lambda f:ee.Feature(f).buffer(15).bounds()),{'layerType':'geeVector'},'All Training')
    
    # Get the RF models for classification, count, and biovolume
    def getModel(training,training_field,outputMode):
        rf= ee.Classifier.smileRandomForest(nTrees)
        rf = rf.setOutputMode(outputMode)
        rf=rf.train(training, training_field, pred_bands)
        return rf

    rf_hab = getModel(training,training_field,'CLASSIFICATION')
    rf_cells=getModel(training,'Cyanobacteria Count (cells/mL)','REGRESSION')
    rf_vol=getModel(training,'Cyanobacteria Biovolume (um^3)','REGRESSION')

    # Get error info
    if getError:
        getRFModelInfoClassification(rf_hab,'HAB_Classification',nTrees,class_names,table_dir)
        getRFModelInfoRegression(rf_cells,'Cell_Count',nTrees,table_dir)
        getRFModelInfoRegression(rf_vol,'Biovolume',nTrees,table_dir)
    
    # Return the models
    return rf_hab,rf_cells,rf_vol,class_names,class_numbers
####################################################################################################
# Function to apply an algal model created in the functions above
def applyAlgalModels(rf_hab,rf_cells,rf_vol,water_model,class_names,class_numbers,applyStudyArea,applyYear = 2022,applyStartJulian = 190, applyEndJulian = 210,pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59'],crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs'):

    # Get a common name ending for the dates provided
    nameEnd = 'yr{} jd{}-{}'.format(applyYear,applyStartJulian,applyEndJulian)

    # Get the water mask and composite for the dates
    water_comp = applyWaterClassification(output_dir,water_model,applyStudyArea,applyYear,applyYear,applyStartJulian,applyEndJulian,crs,transform,water_pred_bands,addToMap=False)
    water = water_comp.select(['Water']).selfMask()
    comp = water_comp.select(pred_bands)

    # Mask composite to extent of water
    compWater = comp.updateMask(water)

    # Add layers to map
    Map.addLayer(applyStudyArea,{},'Apply Study Area')
    Map.addLayer(comp.reproject(crs,transform),getImagesLib.vizParamsFalse,'Median S2 {}'.format(nameEnd),False)
    Map.addLayer(water.selfMask().reproject(crs,transform),{'palette':'00D'},'Median S2 Water {}'.format(nameEnd),False)

    # Apply classification model
    predicted = compWater.classify(rf_hab)
    algalLegendDict={'Algal Negative':'00D','Algal Positive':'D00'}
    Map.addLayer(predicted,{'min':class_numbers[0],'max':class_numbers[-1],'palette':'00D,D00','layerType':'geeImage','classLegendDict':algalLegendDict},'AB Classified {}'.format(nameEnd))

    def getMapStretch(predicted,nameStart,map_stretch_pctls = [5,95]):
        stats_filename = os.path.join(output_dir,'Algal_{}_Stats_{}-{}_{}.txt'.format(nameStart,map_stretch_pctls[0],map_stretch_pctls[1],nameEnd))
        if not os.path.exists(stats_filename):
            predicted_stats= [str(int(i)) for i in predicted.reduceRegion(ee.Reducer.percentile(map_stretch_pctls), applyStudyArea, 3000,tileScale=4).values().getInfo()]
            o = open(stats_filename,'w')
            o.write(','.join(predicted_stats))
            o.close()
        o = open(stats_filename,'r')
        return [int(i) for i in o.read().split(',')]

    # Apply count regression model
    cell_predicted = compWater.classify(rf_cells)
    cell_predicted_stats = getMapStretch(cell_predicted,'Count')
    Map.addLayer(cell_predicted,{'min':cell_predicted_stats[0],'max':cell_predicted_stats[1],'palette':'00D,D00','layerType':'geeImage'},'Cyanobacteria Count (cells/mL) {}'.format(nameEnd))

    # Apply biovolume regression model
    biovolume_predicted = compWater.classify(rf_vol)
    biovolume_predicted_stats = getMapStretch(biovolume_predicted,'Biovolume')
    Map.addLayer(biovolume_predicted,{'min':biovolume_predicted_stats[0],'max':biovolume_predicted_stats[1],'palette':'00D,D00','layerType':'geeImage'},'Cyanobacteria Biovolume (um3) {}'.format(nameEnd))
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
def supervised_algal_mapper(water_training_points,hcb_data,clean_points,\
clean_startYear = 2018,clean_endYear = 2022,clean_startJulian = 213,clean_endJulian = 229,clean_nSamples = 150,\
lat='Sampling Latitude', lon='Sampling Longitude',dateProp = 'Sample Date',properties=[],crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform = [10,0,-2361915.0,0,-10,3177735.0],pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI'],water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',training_field = 'Bloom Type',nTrees=250,outputMode='CLASSIFICATION',getError=False,applyStudyArea=None,applyYear = 2022,applyStartJulian = 190, applyEndJulian = 210):

    # Get the water RF model
    water_model = getSupervisedWaterModel(water_training_points,water_startYear=clean_startYear,water_endYear=clean_endYear,water_startJulian=clean_startJulian,water_endJulian=clean_endJulian,crs = crs,transform = transform,pred_bands=water_pred_bands,output_dir =output_dir)

  
    # Get clean training data
    prepareCleanTrainingData(clean_points,clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,clean_nSamples,crs,transform ,pred_bands,water_pred_bands,output_dir,water_model)

    # Get HCB training data
    hcb_fc = prepareTrainingData(hcb_data,lat,lon,properties)
    batchExtractHCBPredictorTables(hcb_fc,dateProp,crs,transform,pred_bands,water_pred_bands,output_dir,water_model)

    # Get algal RF models
    rf_hab,rf_cells,rf_vol,class_names,class_numbers = getAlgalModels(output_dir,training_field,pred_bands,nTrees,outputMode,getError)

    # Apply algal models
    applyAlgalModels(rf_hab,rf_cells,rf_vol,water_model,class_names,class_numbers,applyStudyArea,applyYear,applyStartJulian, applyEndJulian,pred_bands,water_pred_bands,crs,transform,output_dir)

####################################################################################################
####################################################################################################
####################################################################################################
# Params
water_training_points = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\waterTraining.geojson"
hcb_data = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\HCB Data Sharing.xlsx"
clean_points = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\HABControlPoints\HABControlPoints.geojson"

clean_startYear=2018
clean_endYear=2022
clean_startJulian = 213
clean_endJulian = 229
clean_nSamples = 150

lat='Sampling Latitude'
lon='Sampling Longitude'
dateProp = 'Sample Date'
properties=['Sample Date','Wind Conditions','Cloud Cover','Bloom Description','Bloom Type','Sample Collection Method','Cyanobacteria Count (cells/mL)','Cyanobacteria Biovolume (um^3)']
pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI']
water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'slope', 'aspect', 'hillshade', 'tpi_29', 'tpi_59']
training_field = 'Bloom Type'
outputMode = 'CLASSIFICATION'
nTrees = 150
getError = True
crs = getImagesLib.common_projections['NLCD_CONUS']['crs']
transform = [10,0,-2361915.0,0,-10,3177735.0]
output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs'

applyYear = 2022
applyStartJulian = 213
applyEndJulian = 229
applyStudyArea = ee.Geometry.Polygon(
        [[[-111.19609375, 45.04612495487054],
          [-111.19609375, 41.00037953566339],
          [-104.076953125, 41.00037953566339],
          [-104.076953125, 45.04612495487054]]], None, False)
#####################################################
if __name__ == '__main__':
  
  supervised_algal_mapper(water_training_points,hcb_data,clean_points,\
    clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,clean_nSamples,\
    lat,lon,dateProp,properties,crs,transform,pred_bands,water_pred_bands,output_dir,training_field,nTrees,outputMode,getError,applyStudyArea,applyYear,applyStartJulian,applyEndJulian)

#####################################################
# if  not exportZAndTables:

Map.turnOnInspector()
Map.view()  

#####################################################
#Make table outputs public if needed (for public viewer)
# makeTablesPublic(hab_summary_table_folder)

