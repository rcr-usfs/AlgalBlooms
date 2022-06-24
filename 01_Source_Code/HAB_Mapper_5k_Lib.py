"""
   Copyright 2021 Ian Housman, RedCastle Resources Inc.

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
#Script to map algal blooms in Google Earth Engine and then summarize outputs and provide
#final deliverables for local viewing
#Intended to work within the geeViz package
# python -m pip install geeViz
#Also requires gdal
####################################################################################################
import geeViz.getImagesLib as getImagesLib
import geeViz.taskManagerLib as tml
ee = getImagesLib.ee
Map = getImagesLib.Map

import make_esri_viewer as mew
import json,os,threading,time,glob
from osgeo import gdal
##################################################
#Bring in some existing assets to use for cloud masking
preComputedCloudScoreOffset = getImagesLib.getPrecomputedCloudScoreOffsets(10)
preComputedLandsatCloudScoreOffset = preComputedCloudScoreOffset['landsat']
preComputedSentinel2CloudScoreOffset = preComputedCloudScoreOffset['sentinel2']

#The TDOM stats are the mean and standard deviations of the two bands used in TDOM
#By default, TDOM uses the nir and swir1 bands
preComputedTDOMStats = getImagesLib.getPrecomputedTDOMStats()
preComputedLandsatTDOMIRMean = preComputedTDOMStats['landsat']['mean']
preComputedLandsatTDOMIRStdDev = preComputedTDOMStats['landsat']['stdDev']

preComputedSentinel2TDOMIRMean = preComputedTDOMStats['sentinel2']['mean']
preComputedSentinel2TDOMIRStdDev = preComputedTDOMStats['sentinel2']['stdDev']
############################################################################
############################################################################
#On-the-fly basic water masking method
#This method is used to provide a time-sensitive water mask
#This method tends to work well if there is no wet snow present
#Wet snow over flat areas can result in false positives
def simpleWaterMask(img,contractPixels = 1.5):
  img = getImagesLib.addTCAngles(img);
  ned = ee.Image("USGS/NED").resample('bicubic')
  slope = ee.Terrain.slope(ned.focal_mean(5.5))
  flat = slope.lte(10)
  
  waterMask = img.select(['tcAngleBW']).gte(-0.05)\
    .And(img.select(['tcAngleBG']).lte(0.05))\
    .And(img.select(['brightness']).lt(0.3))\
    .And(flat).focal_min(contractPixels)
  
  return waterMask

#Function to get water mask from JRC water data
#As of August 2021, these data stop in 2018, so they cannot represent dates after
#It can be challenging to get a timely water mask with this method 
def getWaterMask(startYear,endYear,startMonth,endMonth,contractPixels = 2.5):
  startWaterYear = startYear
  endWaterYear = endYear
  if startYear > 2018: startWaterYear = 2018
  if endYear > 2018: endWaterYear = 2018
  permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
  tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
                .filter(ee.Filter.calendarRange(startWaterYear,endWaterYear,'year'))\
                .filter(ee.Filter.calendarRange(startMonth,endMonth,'month')).mode().eq(2).unmask(0).focal_min(contractPixels)

  water_mask = permWater.Or(tempWater)
  return water_mask
############################################################################
#Function to get stats for clean water bodies
#This serves as the control population to compare water bodies against
#Stats are computed for each study area for each month
def getStats(studyAreas,startYear,endYear,startMonth,endMonth,stats_json,bands = ['bloom2','NDGI'],maxTries = 20,crs = 'EPSG:5070',transform = [30,0,-2361915.0,0,-30,3177735.0]):

  #Set up stats outputs
  out_stats = {}
  if os.path.exists(stats_json):
    o = open(stats_json,'r')
    out_stats = json.loads(o.read())
  
  #Iterate across each study area and compute stats for each month
  for studyArea in list(studyAreas.keys()):
    
    sa = studyAreas[studyArea]
    try:
      saBounds = sa.geometry().bounds()
    except:
      saBounds = sa.bounds()
    Map.addLayer(sa,{},studyArea + ' outline',False)
    Map.addLayer(saBounds,{},studyArea + ' bounds',False)
    for month in range(startMonth,endMonth+1):

      sa_month_key = '{}_{}'.format(studyArea,month)
      if sa_month_key not in out_stats.keys():
        print(studyArea,month)
     

        startJulian = int(ee.Date.fromYMD(2003,month,1).format('DD').getInfo())
        endJulian = int(ee.Date.fromYMD(2003,month,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
        
        #Get Landsat and Sentinel 2 images
        # clean_imgs = getImagesLib.getProcessedLandsatAndSentinel2Scenes(
        # 	saBounds,
        # 	startYear,
        # 	endYear,
        # 	startJulian,
        # 	endJulian,
        #   includeSLCOffL7 = True,
        #   convertToDailyMosaics = True,
        #   landsatResampleMethod = 'bicubic',
        #   sentinel2ResampleMethod = 'bicubic',
        # 	preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,
        # 	preComputedSentinel2CloudScoreOffset=preComputedSentinel2CloudScoreOffset,
        # 	preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,
        # 	preComputedLandsatTDOMIRStdDev=preComputedLandsatTDOMIRStdDev,
        # 	preComputedSentinel2TDOMIRMean=preComputedSentinel2TDOMIRMean,
        # 	preComputedSentinel2TDOMIRStdDev=preComputedSentinel2TDOMIRStdDev)
        
        clean_imgs = getImagesLib.getProcessedSentinel2Scenes(\
            saBounds,
            startYear,
            endYear,
            startJulian,
            endJulian,
            applyQABand = False,
            applyCloudScore = False,
            applyShadowShift = False,
            applyTDOM = True,
            cloudScoreThresh = 20,
            performCloudScoreOffset = True,
            cloudScorePctl = 10,
            cloudHeights = ee.List.sequence(500,10000,500),
            zScoreThresh = -1,
            shadowSumThresh = 0.35,
            contractPixels = 1.5,
            dilatePixels = 3.5,
            shadowSumBands = ['nir','swir1'],
            resampleMethod = 'aggregate',
            toaOrSR = 'TOA',
            convertToDailyMosaics = False,
            applyCloudProbability = True,
            preComputedCloudScoreOffset = preComputedSentinel2CloudScoreOffset,
            preComputedTDOMIRMean = preComputedSentinel2TDOMIRMean,
            preComputedTDOMIRStdDev = preComputedSentinel2TDOMIRStdDev,
            cloudProbThresh = 40)

        #Reduce to a median composite
        clean_imgs = clean_imgs.map(getImagesLib.HoCalcAlgorithm2)
        clean_composite = clean_imgs.median()



           #Get water bodies and convert them to a vector to compute stats within
        water_mask = getWaterMask(startYear,endYear,startMonth,endMonth).selfMask()
        # water_mask = simpleWaterMask(clean_composite).selfMask()
        Map.addLayer(clean_composite,getImagesLib.vizParamsFalse,'Clean Composite {} {}'.format(studyArea,month),False)
        
        clean_composite = clean_composite.updateMask(water_mask).select(bands)

        Map.addLayer(clean_composite,{'min':-0.2,'max':0.2,'palette':'00D,DDD,0D0'},'Selected Bands {} {}'.format(studyArea,month),False)
        # clean_lakes = gtacWaterMask.clip(sa).reduceToVectors(scale = 30)
        Map.addLayer(water_mask.clip(sa),{'palette':'00D'},'Clean Lakes {} {}'.format(studyArea,month))
        
        #Compute zonal stats (mean and stdDev) within water areas
        clean_stats = clean_composite.clip(sa).reduceRegion(ee.Reducer.mean().combine(ee.Reducer.stdDev(),'',True),saBounds,None,crs,transform,True,1e13)
  
        #Convert to local json
        stats = None
        tryCount = 1
        def getStatsTryer():return clean_stats.getInfo()
        while stats == None and tryCount < maxTries:
          try:
            print('Computing stats. Try number: ',tryCount)
            stats = getStatsTryer()
            print(stats)
            out_stats[sa_month_key] = stats
            print(out_stats)
          except Exception as e:
            print(e)
        
          tryCount+=1
      else:
        print('Already computed stats for:',studyArea,month)
      o = open(stats_json,'w')
      o.write(json.dumps(out_stats))
      o.close()
  return out_stats

	
############################################################################
def batchMapHABs(summary_areas_dict,study_area_stats_keys,analysisStartYear,analysisEndYear,startMonth,endMonth,clean_stats,reducer = ee.Reducer.percentile([50]),bands = ['bloom2','NDGI'],z_threshs = [1],exportZAndTables = False,hab_summary_table_folder = 'projects/gtac-algal-blooms/assets/outputs/HAB-Summary-Tables',hab_z_imageCollection = 'projects/gtac-algal-blooms/assets/outputs/HAB-Z-Images',crs = 'EPSG:5070',transform = [30,0,-2361915.0,0,-30,3177735.0]):
  task_list = []
  for study_area_stats_key in study_area_stats_keys:
    summary_areas = summary_areas_dict[study_area_stats_key]
    tasks = mapHABs(summary_areas,study_area_stats_key,analysisStartYear,analysisEndYear,startMonth,endMonth,clean_stats,study_area_stats_key,reducer,bands,z_threshs,exportZAndTables,hab_summary_table_folder,hab_z_imageCollection,crs,transform)
    task_list.extend(tasks)
  return task_list
############################################################################
#Function to map algal blooms
#These aren't always harmful
#This function produces a raster z score of algal bloom and a summary table
def mapHABs(summaryAreas,studyAreaName,analysisStartYear,analysisEndYear,startMonth,endMonth,clean_stats,stats_sa,reducer = ee.Reducer.percentile([50]),bands = ['bloom2','NDGI'],z_threshs = [1],exportZAndTables = False,hab_summary_table_folder = 'projects/gtac-algal-blooms/assets/outputs/HAB-Summary-Tables',hab_z_imageCollection = 'projects/gtac-algal-blooms/assets/outputs/HAB-Z-Images',crs = 'EPSG:5070',transform = [30,0,-2361915.0,0,-30,3177735.0]):

  #Set up output areas and dates
  try:
    saBounds = summaryAreas.geometry().bounds()
  except:
    saBounds = summaryAreas.bounds()

  Map.addLayer(summaryAreas,{},'Study Area',False)
  Map.addLayer(saBounds,{},'Clip Area',False)


  startJulian = int(ee.Date.fromYMD(analysisStartYear,startMonth,1).format('DD').getInfo())
  endJulian = int(ee.Date.fromYMD(analysisEndYear,endMonth,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
  print(analysisStartYear,analysisEndYear,startMonth,endMonth, startJulian,endJulian)
  
  #Get Landsat and Sentinel 2 images
  # dirty_imgs = getImagesLib.getProcessedLandsatAndSentinel2Scenes(
  #   saBounds,
  #   analysisStartYear,
  #   analysisEndYear,
  #   startJulian,
  #   endJulian,
  #   includeSLCOffL7 = True,
  #   convertToDailyMosaics = True,
  #   landsatResampleMethod = 'bicubic',
  #   sentinel2ResampleMethod = 'bicubic',
  #   applyTDOMLandsat = True,
  #   applyTDOMSentinel2 = True,
  #   preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,
  #   preComputedSentinel2CloudScoreOffset=preComputedSentinel2CloudScoreOffset,
  #   preComputedLandsatTDOMIRMean = preComputedSentinel2TDOMIRMean,
  #   preComputedLandsatTDOMIRStdDev=preComputedSentinel2TDOMIRStdDev,
  #   preComputedSentinel2TDOMIRMean=preComputedSentinel2TDOMIRMean,
  #   preComputedSentinel2TDOMIRStdDev=preComputedSentinel2TDOMIRStdDev)

  dirty_imgs = getImagesLib.getProcessedSentinel2Scenes(\
    saBounds,
    analysisStartYear,
    analysisEndYear,
    startJulian,
    endJulian,
    applyQABand = False,
    applyCloudScore = False,
    applyShadowShift = False,
    applyTDOM = True,
    cloudScoreThresh = 20,
    performCloudScoreOffset = True,
    cloudScorePctl = 10,
    cloudHeights = ee.List.sequence(500,10000,500),
    zScoreThresh = -1,
    shadowSumThresh = 0.35,
    contractPixels = 1.5,
    dilatePixels = 3.5,
    shadowSumBands = ['nir','swir1'],
    resampleMethod = 'bicubic',
    toaOrSR = 'TOA',
    convertToDailyMosaics = False,
    applyCloudProbability = True,
    preComputedCloudScoreOffset = preComputedSentinel2CloudScoreOffset,
    preComputedTDOMIRMean = preComputedSentinel2TDOMIRMean,
    preComputedTDOMIRStdDev = preComputedSentinel2TDOMIRStdDev,
    cloudProbThresh = 40)

  dirty_imgs = dirty_imgs.map(getImagesLib.HoCalcAlgorithm2)

  #Specify which tasks for tracking
  task_list = []
  #Iterate across each year and month
  for analysisYear in range(analysisStartYear,analysisEndYear+1):
    for month in range(startMonth,endMonth+1):

      #Filter out images for the year-month
      startJulianT = int(ee.Date.fromYMD(analysisYear,month,1).format('DD').getInfo())
      endJulianT = int(ee.Date.fromYMD(analysisYear,month,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
      dirty_imgsT = dirty_imgs.filter(ee.Filter.calendarRange(analysisYear,analysisYear,'year'))\
                              .filter(ee.Filter.calendarRange(startJulianT,endJulianT))

      #Water mask can be computed as the mode of the water mask for each image
      # gtacWaterMask = dirty_imgsT.map(simpleWaterMask).mode()

      #Or simply by computing it from the median of the images
      gtacWaterMask = simpleWaterMask(dirty_imgsT.median())

      #JRC water mask isn't currently used since it commits large areas at edges of water bodies as water which
      #results in false positives
      # jrcWaterMask = getWaterMask(analysisYear,analysisYear,month,month)
      # print(dirty_imgsT.size().getInfo(),analysisYear,month,startJulianT,endJulianT)
  
  
  # # Map.addLayer(waterMask,{'min':1,'max':1,'palette':'00F'},'Water Mask')

  
 
      Map.addLayer(dirty_imgsT.median(),getImagesLib.vizParamsFalse,'Composite yr{} m{}'.format(analysisYear,month),False)
      Map.addLayer(dirty_imgsT.median().select(bands).updateMask(gtacWaterMask),{'min':-0.2,'max':0.2,'palette':'00D,DDD,0D0'},'Selected Bands yr{} m{}'.format(analysisYear,month),False)
      Map.addLayer(gtacWaterMask.selfMask(),{'min':1,'max':1,'palette':'00F','classLegendDict':{'Water':'00F'}},'GTAC Water Mask yr{} m{}'.format(analysisYear,month),False)
      # Map.addLayer(jrcWaterMask.selfMask(),{'min':1,'max':1,'palette':'00F','classLegendDict':{'Water':'00F'}},'JRC Water Mask yr{} m{}'.format(analysisYear,month),False)
      dirty_imgsT = dirty_imgsT.select(bands)

      #Specify which water mask to use
      waterMask = gtacWaterMask

      #Bring in the stats from the clean lakes
      means = ee.Image([clean_stats['{}_{}'.format(stats_sa,month)][band+'_mean'] for band in bands])
      stdDevs = ee.Image([clean_stats['{}_{}'.format(stats_sa,month)][band+'_stdDev'] for band in bands])

      #Compute the z score based on the stats from the clean lakes
      dirty_zs = dirty_imgsT.map(lambda img: img.subtract(means).divide(stdDevs).updateMask(waterMask))
      
      #Summarize based on the provided reducer, and then take the max of that if multiple bands were used
      dirty_z = dirty_zs.reduce(reducer).reduce(ee.Reducer.max())

      #Threshold hab/not hab
      hab = dirty_z.gte(z_threshs)
      bns = hab.bandNames()
      # print('here')
      bns = list(map(lambda z:'AB_Pct_Z_' + f'{int(z*100):03}' ,z_threshs))
      # print('bns:',bns)
      hab = hab.rename(bns)
      # print(hab.bandNames().getInfo())

      Map.addLayer(dirty_z,{'min':0,'max':2,'palette':'00F,0F0,F00'},'Dirty Z yr{} m{}'.format(analysisYear,month),False)
      Map.addLayer(hab,{'min':0,'max':1,'classLegendDict':{'Not HAB':'00F','HAB':'F00'}},'AB yr{} m{}'.format(analysisYear,month),False)
      
      #Function to convert to pct of an area as hab
      def toPct(f,inField,outField):
        a = f.get(inField);
        a = ee.Array(a).slice(1,1,2).project([0])
        sum = ee.Number(a.reduce(ee.Reducer.sum(),[0]).get([0]))
        a = ee.Number(a.toList().get(1))
        pct = ee.Number(a.divide(sum).multiply(100))
        return f.set({inField:pct,outField:pct,'year':analysisYear,'month':month})

      #Set up output table and raster for export
      if exportZAndTables:

        #Summarize the pct of each area mapped as algal bloom
        summary_table = hab.reduceRegions(summaryAreas, ee.Reducer.fixedHistogram(0, 2, 2), None, crs, transform, 4)
        if len(z_threshs) == 1:
          summary_table = summary_table.filter(ee.Filter.notNull(['histogram']))
          summary_table = summary_table.map(lambda i:toPct(i,'histogram',bns[0]))
        
        else:
          for bn in bns:
            summary_table = summary_table.filter(ee.Filter.notNull([bn]))
            summary_table = summary_table.map(lambda i:toPct(i,bn,bn))

        #Summarize the pct of each area mapped as water
        summary_table = waterMask.reduceRegions(summary_table, ee.Reducer.fixedHistogram(0, 2, 2), None, crs, transform, 4)
        summary_table = summary_table.map(lambda i:toPct(i,'histogram','Pct_Water'))
        
        summary_table = dirty_z.reduceRegions(summary_table, ee.Reducer.mean(), None, crs, transform, 4)
        summary_table = summary_table.map(lambda f: f.set('AB_Mean_Z',f.get('mean')))
    
        #Export output table
        output_table_name = '{}_HAB_Summary_Table_yr{}_m{}'.format(studyAreaName,analysisYear,month)
        output_z_name = '{}_HAB_Z_yr{}_m{}'.format(studyAreaName,analysisYear,month)
        task_list.append(output_table_name)
        t = ee.batch.Export.table.toAsset(summary_table,\
                  description=output_table_name,\
                  assetId=hab_summary_table_folder + '/'+output_table_name)
        t.start()
       
        #Set up the z score raster for export and export it
        dirty_z_for_export = dirty_z.multiply(1000).int16().set({'system:time_start':ee.Date.fromYMD(analysisYear,month,1).millis(),
                                                                  'year':analysisYear,
                                                                  'month':month,
                                                                  'studyAreaName':studyAreaName,
                                                                  'zThresh':z_thresh})
        # task_list.append(output_z_name)
        t = ee.batch.Export.image.toAsset(dirty_z_for_export.clip(saBounds), 
                      description = output_z_name, 
                      assetId = hab_z_imageCollection + '/'+ output_z_name, 
                      pyramidingPolicy = {'.default':'mean'}, 
                      dimensions = None, 
                      region = None, 
                      scale = None, 
                      crs = crs, 
                      crsTransform = transform, 
                      maxPixels = 1e13)
        print('Exporting:',output_z_name)
        print(t)
        t.start()
  return task_list
     
############################################################################
#Make all assets public
def makeTablesPublic(table_dir):
  tables = ee.data.getList({'id':table_dir})
  for table in tables:
      table = table['id']
      print('Making public: ',table)
      ee.data.setAssetAcl(table, json.dumps({u'writers': [], u'all_users_can_read': True, u'readers': []}))
###############################################################
def limitThreads(limit):
  while threading.activeCount() > limit:
    time.sleep(1)
    print(threading.activeCount(),'threads running')

#Function to summarize output algal bloom tables over multiple stats and z thresholds simultaneously
def batchSummarizeTables(input_table_dir,output_table_dir,startYear,endYear,startMonth,endMonth,summary_areas_dict,states = ['WA','OR'],z_threshs = [1],pct_HAB_threshold = 5,unique_area_field = 'UNID',viewer_template = 'HAB_Mapper_Viewer_Template.html'):
  for z_thresh in z_threshs:
    for state in states:
      tt = threading.Thread(target = summarizeTables, args = (input_table_dir,output_table_dir,startYear,endYear,startMonth,endMonth,summary_areas_dict,[state],z_thresh,pct_HAB_threshold,unique_area_field,viewer_template))
      tt.start()
      time.sleep(0.1)
      limitThreads(4)
  limitThreads(1)
###############################################################
#Simplify tables
def summarizeTables(input_table_dir,output_table_dir,startYear,endYear,startMonth,endMonth,summary_areas_dict,states = ['WA','OR'],z_thresh = 1,pct_HAB_threshold = 5,unique_area_field = 'UNID',viewer_template = 'HAB_Mapper_Viewer_Template.html'):
  #Potential fields available
  # keep_fields = ['ADMINFORES', 'AreaSqKm', 'DISTRICTNA', 'DISTRICTNU', 'DISTRICTOR', 'Elevation', 'FCode', 'FDate', 'FID_FS_Bou', 'FID_S_USA_', 'FID_WA_FS_', 'FID_WA_Nam', 'FORESTNAME', 'FORESTNUMB', 'FORESTORGC', 'FType', 'GIS_ACRES', 'GIS_ACRE_1', 'GNIS_ID', 'GNIS_Name', 'OBJECTID', 'Permanent_', 'RANGERDIST', 'REGION', 'REGION_1', 'ReachCode', 'Resolution', 'SHAPE_AR_1', 'SHAPE_LEN', 'SHAPE__Are', 'SHAPE__Len', 'Shape_Area', 'Shape_Leng', 'Visibility', 'state']
  #Specify which fields are kept in the summarized table
  keep_fields = ['UNID','state','REGION','FORESTNAME','GNIS_Name']

  state_dict = {'WA':'Washington','OR':'Oregon','WY':'Wyoming'}

  month_dict = {1:'January',2:'February',3:'March',4:'April',5:'May',6:'June',7:'July',8:'August',9:'September',10:'October',11:'November',12:'December'}

  #Set up output dir
  if not os.path.exists(output_table_dir):os.makedirs(output_table_dir)

  #Get date ranges
  years = range(startYear,endYear+1)
  months = range(startMonth,endMonth+1)

  #Filter available tables
  print('Reading in tables')
  tables = ee.data.getList({'id':input_table_dir})
  tables = [i['id'] for i in tables]
  tables = [i for i in tables if i.split('/')[-1].split('_')[0] in states]
  tables = [i for i in tables if int(i.split('/')[-1].split('_yr')[-1].split('_')[0]) in years]
  tables = [i for i in tables if int(i.split('/')[-1].split('_m')[-1]) in months]
  
  z_thresh_padded =  f'{int(z_thresh*100):03}' 
  ab_field = 'AB_Pct_Z_' +z_thresh_padded
  #Handle creating a unid to summarize across
  #Currently Wyoming has a different pair of fields needed than WA and
  def getUNID(t,state):
    if state == 'WY':
      return t.set('UNID',ee.String(t.get('RANGERDIST')).cat('_').cat(ee.String(t.get('ReachCode'))))
    else:
      return t.set('UNID',ee.Number(t.get('FID_S_USA_')).format().cat('_').cat(ee.Number(t.get('FID_{}_FS_'.format(state))).format()))
    
  #Read in tables and set up a UNID
  def getTable(table_name):
    f = ee.FeatureCollection(table_name)
    
    state = table_name.split('/')[-1].split('_')[0]
    f = f.map(lambda ft: ft.set({'state':state}))
    f = f.map(lambda t:getUNID(t,state))
    f = f.map(lambda t: t.set('is_hab',ee.Number(t.get(ab_field)).gt(pct_HAB_threshold)))
    f = f.map(lambda t: t.set('is_not_hab',ee.Number(t.get(ab_field)).lte(pct_HAB_threshold)))

    # n_featues = f.size().getInfo()
    # n_unids = f.aggregate_histogram('UNID').keys().length().getInfo()
    # print('N features/UNIDs{}:'.format(table_name),n_featues,n_unids)

    return f

  tables = list(map(getTable,tables))
  tables = ee.FeatureCollection(tables).flatten()#.sort('date')

  
  
  #Read in state study areas for getting geometry from and add UNID
  print('Reading in summary area spatial info')
  study_areas = list(map(lambda s:ee.FeatureCollection(summary_areas_dict[s]).map(lambda f:getUNID(f,s)),states))
  study_areas = ee.FeatureCollection(study_areas).flatten()

  #Select fields to keep
  study_areas = study_areas.select(keep_fields)

  #Convert to geojson and get the features
  study_areas_geojson = study_areas.getInfo()
  study_areas = study_areas_geojson['features']

  #Summarize all the tables using several fields and reducers
  print('Summarizing tables')
  stats = tables.reduceColumns(**{
    'selectors': ['UNID','is_hab','is_not_hab',ab_field,'AB_Mean_Z','Pct_Water'],
    'reducer':ee.Reducer.max()
    .combine(ee.Reducer.sum(),None,True)
    .combine(ee.Reducer.mean(),None,True)
    .repeat(5)
    .group(**{
      'groupField': 0,
      'groupName': 'UNID',
    })}).get('groups')

  #Extract the proper fields/stats and give it a descriptive field name
  def extractSummary(stat):
    stat = ee.Dictionary(stat)
    max = ee.List(stat.get('max'))
    sum = ee.List(stat.get('sum'))
    mean = ee.List(stat.get('mean'))
    is_hab = ee.Number(max.get(0)).byte()
    pos_ct = ee.Number(sum.get(0)).byte()
    neg_ct = ee.Number(sum.get(1)).byte()
    pct_water = ee.Number(mean.get(4)).format('%.2f')
    max_pct_hab = ee.Number(max.get(2)).format('%.2f')
    max_mean_z = ee.Number(max.get(3)).format('%.2f')
    out_dict = ee.Dictionary({
      'UNID':stat.get('UNID'),
      'A_Any_Pos':is_hab,
      'A_Pos_Ct':pos_ct,
      'A_Neg_Ct':neg_ct,
      'Pct_Water':pct_water,
      'A_Max_Pct':max_pct_hab,
      'A_Max_Zmean':max_mean_z
      })
    return out_dict

  out = ee.List(stats).map(extractSummary).getInfo()

  #Function to join the summary dictionary to the feature from the original summary areas
  def joinGeo(stat):
    f = [i for i in study_areas if i['properties']['UNID'] == stat['UNID']][0]
    f['properties'].update(stat)
    return f
  
  #Join the geometry and keep fields from original study areas to the summary stats
  print('Joining summary stats to spatial info')
  out_geo = list(map(joinGeo,out))

  #Finalize the geojson formatting
  tables_json = study_areas_geojson
  tables_json['features'] = out_geo


  out_name = os.path.join(output_table_dir,'Algal_Bloom_Summaries_{}_yrs{}-{}_mths{}-{}_zthresh{}'.format('-'.join(states),startYear,endYear,startMonth,endMonth,z_thresh_padded))

  print('Writing final outputs')
  #Write out combined geojson
  out_name_geojson = out_name+ '.geojson'
  o = open(out_name_geojson,'w')
  o.write(json.dumps(tables_json))
  o.close()


  #Convert to csv
  out_name_csv = out_name + '.csv'
  tables = [i['properties'] for i in tables_json['features']]
  header = ','.join([str(i) for i in list(tables[0].keys())])
  out = [header]
  for table in tables:out.append(','.join([str(i) for i in list(table.values())]))
  out = '\n'.join(out)
  o = open(out_name_csv,'w')
  o.write(out)
  o.close()

  #Set up self-contained viewer
  out_name_html = out_name+'_viewer.html'
  
  #Set up different fields that need replaced in the template
  title = ' '.join(os.path.basename(out_name).split('_'))
 
  geojson_title = title
  title = title + ' Viewer'
  # print('Theyre equal:',startYear==endYear)
  if startYear == endYear:
    year_range_string = '{}'.format(startYear)
  elif endYear - startYear == 1:
    year_range_string = '{} and {}'.format(startYear,endYear)
  else:
    year_range_string = 'years {} to {}'.format(startYear,endYear)

  if len(states) == 1:
    states_list_string = state_dict[states[0]]

  elif len(states) == 2:
    states_list_string = state_dict[states[0]] + ' and ' + state_dict[states[1]]

  else:
    states_list_string = ', '.join([state_dict[i] for i in states[:-1]]) + ', and '+ state_dict[states[-1]]
   

  replace_dict = {'{TITLE}':title,
                  '{RAW_GEOJSON}':json.dumps(tables_json),
                  '{START_MONTH}':month_dict[startMonth],
                  '{END_MONTH}':month_dict[endMonth],
                  '{YEAR_RANGE}':year_range_string,
                  '{startYear}':str(startYear),
                  '{endYear}':str(endYear),
                  '{startMonth}':str(startMonth),
                  '{endMonth}':str(endMonth),
                  '{STATES}':str(states),
                  '{PCT_THRESH}':str(pct_HAB_threshold),
                  '{GEOJSON_TITLE}':'geojson_title',
                  '{Z_THRESH}':str(z_thresh),
                  '{Z_PCTL}':'80',
                  '{STATES_LIST}':states_list_string
                  }

  mew.setup_viewer(out_name_html, replace_dict, viewer_template)

  #Convert geojson to shapefile
  out_name_shp = out_name+'.shp'
  srcDS = gdal.OpenEx(out_name_geojson)
  ds = gdal.VectorTranslate(out_name_shp, srcDS, format='ESRI Shapefile' ,options = '-skipfailures')

############################################################################

dirty_lakes = {}
dirty_lakes['billy_chinook'] = ee.Geometry.Polygon(\
        [[[-121.48239392202757, 44.62035796573114],\
          [-121.48239392202757, 44.509308410535326],\
          [-121.2262751476135, 44.509308410535326],\
          [-121.2262751476135, 44.62035796573114]]], None, False)
        
dirty_lakes['odell_lake'] = ee.Geometry.Polygon(\
        [[[-122.0549625378963, 43.59540724921347],\
          [-122.0549625378963, 43.547151047382215],\
          [-121.95436897100177, 43.547151047382215],\
          [-121.95436897100177, 43.59540724921347]]], None, False)
dirty_lakes['wy_dirty_lake_combo'] = ee.Geometry.MultiPolygon(\
       [[[[-108.22911901517031, 43.41087661396606],\
          [-108.22911901517031, 43.18224132965033],\
          [-108.13642187161562, 43.18224132965033],\
          [-108.13642187161562, 43.41087661396606]]],\
        [[[-109.31414671598651, 44.50635067370683],\
          [-109.31414671598651, 44.44044923506513],\
          [-109.15587492643573, 44.44044923506513],\
          [-109.15587492643573, 44.50635067370683]]],\
        [[[-105.24916665039864, 41.18222447048324],\
          [-105.24916665039864, 41.17072504799633],\
          [-105.21994130097237, 41.17072504799633],\
          [-105.21994130097237, 41.18222447048324]]],\
        [[[-110.00372861158245, 43.69250055304031],\
          [-110.00372861158245, 43.68871474128848],\
          [-109.99746297132366, 43.68871474128848],\
          [-109.99746297132366, 43.69250055304031]]],\
        [[[-110.68328269422415, 42.011520615832595],\
          [-110.68328269422415, 41.96418338823834],\
          [-110.64277060926321, 41.96418338823834],\
          [-110.64277060926321, 42.011520615832595]]],\
        [[[-105.73231371755158, 41.900976006257636],\
          [-105.73231371755158, 41.86340128602617],\
          [-105.70553454274689, 41.86340128602617],\
          [-105.70553454274689, 41.900976006257636]]]], None, False)
dirty_lakes['wy_turbid_lake_combo'] = ee.Geometry.MultiPolygon(\
       [[[[-108.64299538824541, 43.21152592379036],\
          [-108.64299538824541, 43.155449940915105],\
          [-108.56197121832354, 43.155449940915105],\
          [-108.56197121832354, 43.21152592379036]]],\
        [[[-111.03544902486978, 41.505363985440525],\
          [-111.03544902486978, 41.44465706685474],\
          [-111.00249004049478, 41.44465706685474],\
          [-111.00249004049478, 41.505363985440525]]],
        [[[-104.92643828230109, 44.39755383373756],\
          [-104.92643828230109, 44.321459875907976],\
          [-104.7458505137464, 44.321459875907976],\
          [-104.7458505137464, 44.39755383373756]]],\
        [[[-109.39235326037969, 42.234835938460854],\
          [-109.39235326037969, 42.21411566114201],\
          [-109.34772130237188, 42.21411566114201],\
          [-109.34772130237188, 42.234835938460854]]],\
        [[[-109.47472392655484, 42.29657504761792],\
          [-109.47472392655484, 42.245002911853064],\
          [-109.40829097367399, 42.245002911853064],\
          [-109.40829097367399, 42.29657504761792]]]], None, False)