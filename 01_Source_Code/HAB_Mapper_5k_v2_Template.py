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

# Script to map algal blooms in Google Earth Engine
# Intended to work within the geeViz package
# python -m pip install geeViz
# Also requires gdal, pandas, openpyxl
####################################################################################################
import os,sys,threading,json,pandas,time,glob
import geeViz.getImagesLib as getImagesLib
import numpy as np

ee = getImagesLib.ee
Map = getImagesLib.Map
##################################################
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
#####################################################
def prepareTrainingData(hcb_data,lat='Sampling Latitude', lon='Sampling Longitude',properties=[]):
  print(hcb_data)
  hcb_df = pandas.read_excel(hcb_data)
  
  for c in hcb_df.columns[hcb_df.dtypes=='datetime64[ns]']:
    hcb_df[c]= hcb_df[c].dt.strftime('%Y-%m-%d')
  
  hfb_json = df_to_geojson(hcb_df, properties, lat, lon)
  # print(hfb_json)
  hfb_json = ee.FeatureCollection(hfb_json)
  return hfb_json
#####################################################

def extractPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,d):

  startDate = ee.Date(d).advance(-2,'week')
  endDate = ee.Date(d).advance(2,'week')
  startYear = int(startDate.get('year').getInfo())
  endYear = int(endDate.get('year').getInfo())
  startJulian = int(startDate.format('DDD').getInfo())
  endJulian = int(endDate.format('DDD').getInfo())

  hcb_d = hcb_fc.filter(ee.Filter.eq(dateProp,d))

  comp = getImagesLib.getProcessedSentinel2Scenes(hcb_d,startYear-1,endYear+1,startJulian,endJulian,resampleMethod='bicubic',convertToDailyMosaics = False)
  comp = comp.filter(ee.Filter.calendarRange(startYear,endYear,'year'))
  comp = comp.map(getImagesLib.HoCalcAlgorithm2).select(pred_bands).median()
  water = getImagesLib.simpleWaterMask(comp)
  compWater = comp.updateMask(water)
  # k =ee.Kernel.square(2.5,'pixels',False)

  
  # print(d, hcb_d.size().getInfo(),startDate.format('YYYY-MM-dd').getInfo(),images.size().getInfo())
  # Map.addLayer(comp,getImagesLib.vizParamsFalse,'Median S2 {}'.format(d),False)
  # Map.addLayer(water.selfMask(),{'palette':'00D'},'Median S2 Water {}'.format(d),False)
  compWater = compWater.focalMean(2.5);
  # Map.addLayer(compWater,getImagesLib.vizParamsFalse,'NB Median S2 {}'.format(d),False)

  # Map.addLayer(hcb_d.map(lambda f:ee.Feature(f).buffer(250).bounds()),{'layerType':'geeVector'},'HCB {}'.format(d),False)
  
  
  

  t = compWater.reduceRegions(hcb_d, ee.Reducer.first(), None, crs, transform, 4)
  return t
 

def writePredictorTables(featureCollection,filename):
  out_f = ee.FeatureCollection(featureCollection).getInfo()
  o = open(filename,'w')
  o.write(json.dumps(out_f))
  o.close()
    
def batchExtractPredictorTables(hcb_fc,dateProp,crs,transform,pred_bands,output_dir):
  if not os.path.exists(output_dir):os.makedirs(output_dir)
  dates = hcb_fc.aggregate_histogram(dateProp).keys().getInfo()
  date_sets = new_set_maker(dates,20)
  date_set_i = 1
  for date_set in date_sets:
    # print(date_set)
    of = os.path.join(output_dir,'HCB_Training_{}.geojson'.format(date_set_i))
    if not os.path.exists(of):
      out_f = []
      for d in date_set:
        t = extractPredictorValues(hcb_fc,dateProp,crs,transform,pred_bands,d)
        out_f.append(t)
      out_f = ee.FeatureCollection(out_f).flatten()
      writePredictorTables(out_f,of)
    date_set_i +=1

    
#####################################################
def getAlgalModel(table_dir,training_field,pred_bands,nTrees=100,outputMode='CLASSIFICATION',getError=False):
  tables = glob.glob(os.path.join(table_dir,'*.geojson'))
  training = []
  for t in tables:
    o = open(t,'r')
    training.append(ee.FeatureCollection(json.load(o)))
    o.close()
  training = ee.FeatureCollection(training).flatten()

  class_names= training.aggregate_histogram(training_field).keys().getInfo()
  class_numbers = ee.List.sequence(1,len(class_names))
  training = training.remap(class_names, class_numbers, training_field)
  training=training.filter(ee.Filter.notNull(pred_bands))
  # print(class_names,class_numbers.getInfo())

  rf= ee.Classifier.smileRandomForest(nTrees)
  rf = rf.setOutputMode(outputMode)
  rf=rf.train(training, training_field, pred_bands)

  if getError:
    cm = rf.confusionMatrix()
    exp = ee.Dictionary({'importance':rf.explain().get('importance'),\
                          'kappa':cm.kappa(),\
                          'oobError':rf.explain().get('outOfBagErrorEstimate'),
                          'cm':cm}).getInfo()

    oobAcc = 1-exp['oobError']
    kappa = exp['kappa']
    cm = exp['cm']

    importance = list(zip(list(exp['importance'].keys()),list(exp['importance'].values())))
    importance = sorted(importance, key=lambda row: row[1]*-1)
    importance_lines = 'Name,Importance\n'
    for line in importance:
      importance_lines+=','.join([str(i) for i in line])+'\n'
    print(importance_lines)
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

      
      # print(ri,row,col,v,rowSum,colSum,rowAcc,colAcc)
    overall_acc = np.sum(diag)/allSum

    # out_table_lines+='{},{},{},{},{}\n'.format(run,field,oobAcc,overall_acc,kappa)

    out_lines = ',{},Commission Error (1-n)\n'.format(','.join(class_names))
    for ri in range(0,len(cm)):
      out_lines+='{},{},{}\n'.format(class_names[ri],','.join([str(i) for i in cm[ri]]),rowAccs[ri])
    out_lines+='Omission Error (1-n),{},,\n'.format(','.join([str(i) for i in colAccs]))
    out_lines+='Overall Accuracy,{}\n'.format(overall_acc)
    out_lines+='OOB Accuracy,{}\n'.format(oobAcc)
    out_lines+='Kappa,{}\n'.format(overall_acc)
    print(out_lines,rowAccs,colAccs,kappa)

    out_table_filename = os.path.join(table_dir,'HAB_{}_OOB_Error_Matrix.csv'.format(nTrees))
    o = open(out_table_filename,'w')
    o.write(out_lines)
    o.close()

    out_importance_filename = os.path.join(table_dir,'HAB_{}_Var_Importance.csv'.format(nTrees))
    o = open(out_importance_filename,'w')
    o.write(importance_lines)
    o.close()
  return rf
#####################################################
def applyRFModel(rfModel,applyStudyArea,applyYear = 2022,applyStartJulian = 190, applyEndJulian = 210):
  
  comp = getImagesLib.getProcessedSentinel2Scenes(applyStudyArea,applyYear-1,applyYear+1,applyStartJulian,applyEndJulian,resampleMethod='bicubic',convertToDailyMosaics = False)
  comp = comp.filter(ee.Filter.calendarRange(applyYear,applyYear,'year'))
  comp = comp.map(getImagesLib.HoCalcAlgorithm2).select(pred_bands).median()
  water = getImagesLib.simpleWaterMask(comp)
  compWater = comp.updateMask(water)
  Map.addLayer(applyStudyArea,{},'Apply Study Area')
  Map.addLayer(comp,getImagesLib.vizParamsFalse,'Median S2',False)
  Map.addLayer(water.selfMask(),{'palette':'00D'},'Median S2 Water',False)

  predicted = compWater.classify(rfModel)
  Map.addLayer(predicted.randomVisualizer(),{},'Classified')
#####################################################
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
#####################################################
def limitThreads(limit):
  while threading.active_count() > limit:
    time.sleep(2)
    print(threading.active_count(),'threads running')
#####################################################
def supervised_algal_mapper(hcb_data,lat='Sampling Latitude', lon='Sampling Longitude',dateProp = 'Sample Date',properties=[],crs = getImagesLib.common_projections['NLCD_CONUS']['crs'],transform= getImagesLib.common_projections['NLCD_CONUS']['transform'],pred_bands = ['blue','green','red','nir','nir2','swir1','swir2','NBR','NDVI','NDMI','NDSI','brightness','greenness','wetness','tcAngleBG','NDGI'],output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs',training_field = 'Bloom Type',applyStudyArea=None,applyYear = 2022,applyStartJulian = 190, applyEndJulian = 210):
  hcb_fc = prepareTrainingData(hcb_data,lat,lon,properties)
  batchExtractPredictorTables(hcb_fc,dateProp,crs,transform,pred_bands,output_dir)
  rfModel = getAlgalModel(output_dir,training_field,pred_bands)
  applyRFModel(rfModel,applyStudyArea,applyYear,applyStartJulian, applyEndJulian)
  # print(hcb_df.info())
  # print(hfb_json)
#####################################################
hcb_data = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\HCB Data Sharing.xlsx"
lat='Sampling Latitude'
lon='Sampling Longitude'
dateProp = 'Sample Date'
properties=['Sample Date','Wind Conditions','Cloud Cover','Bloom Description','Bloom Type','Sample Collection Method','Cyanobacteria Count (cells/mL)','Cyanobacteria Biovolume (um^3)']
pred_bands = ['blue','green','red','nir','nir2','swir1','swir2','NBR','NDVI','NDMI','NDSI','brightness','greenness','wetness','tcAngleBG','NDGI']
training_field = 'Bloom Type'
crs = getImagesLib.common_projections['NLCD_CONUS']['crs']
transform = getImagesLib.common_projections['NLCD_CONUS']['transform']
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
  supervised_algal_mapper(hcb_data,lat,lon,dateProp,properties,crs,transform,pred_bands,output_dir,training_field,applyStudyArea,applyYear,applyStartJulian,applyEndJulian)

#####################################################
# if  not exportZAndTables:
Map.turnOnInspector()
Map.view()  

#####################################################
#Make table outputs public if needed (for public viewer)
# makeTablesPublic(hab_summary_table_folder)

