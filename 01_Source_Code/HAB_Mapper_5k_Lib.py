import geeViz.getImagesLib as getImagesLib
ee = getImagesLib.ee
Map = getImagesLib.Map
import json,os
##################################################
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
def getWaterMask(startYear,endYear,startMonth,endMonth,contractPixels = 2.5):
  startWaterYear = startYear
  endWaterYear = endYear
  if startYear > 2018: startWaterYear = 2018
  if endYear > 2018: endWaterYear = 2018
  permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
  tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
                .filter(ee.Filter.calendarRange(startWaterYear,endWaterYear,'year'))\
                .filter(ee.Filter.calendarRange(startMonth,endMonth,'month')).mode().eq(2).unmask(0).focal_min(contractPixels)

  water_mask = permWater.Or(tempWater).selfMask()
  return water_mask
############################################################################
def getStats(studyAreas,startYear,endYear,startMonth,endMonth,stats_json,bands = ['bloom2','NDGI'],maxTries = 10,scale = 30,crs = 'EPSG:5070'):
  out_stats = {}
  if os.path.exists(stats_json):
    o = open(stats_json,'r')
    out_stats = json.loads(o.read())
  
  for studyArea in list(studyAreas.keys()):
    
    sa = studyAreas[studyArea]
    for month in range(startMonth,endMonth+1):

      sa_month_key = '{}_{}'.format(studyArea,month)
      if sa_month_key not in out_stats.keys():
        print(studyArea,month)
        water_mask = getWaterMask(startYear,endYear,startMonth,endMonth)
        clean_lakes = water_mask.clip(sa).reduceToVectors(scale = 30)
        Map.addLayer(clean_lakes,{},'Clean Lakes {} {}'.format(studyArea,month))
        try:
          saBounds = sa.geometry().bounds()
        except:
          saBounds = sa.bounds()

        startJulian = int(ee.Date.fromYMD(2003,month,1).format('DD').getInfo())
        endJulian = int(ee.Date.fromYMD(2003,month,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
        
        clean_imgs = getImagesLib.getProcessedLandsatAndSentinel2Scenes(
        	saBounds,
        	startYear,
        	endYear,
        	startJulian,
        	endJulian,
        	preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,
        	preComputedSentinel2CloudScoreOffset=preComputedSentinel2CloudScoreOffset,
        	preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,
        	preComputedLandsatTDOMIRStdDev=preComputedLandsatTDOMIRStdDev,
        	preComputedSentinel2TDOMIRMean=preComputedSentinel2TDOMIRMean,
        	preComputedSentinel2TDOMIRStdDev=preComputedSentinel2TDOMIRStdDev)
        
        clean_composite = clean_imgs.map(getImagesLib.HoCalcAlgorithm2).median()

        clean_stats = clean_composite.select(bands).reduceRegion(ee.Reducer.mean().combine(ee.Reducer.stdDev(),'',True),sa,scale,crs,None,True,1e13)
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
  # Map.addLayer(clean_composite,getImagesLib.vizParamsFalse,'Clean Lakes Composite')
	
############################################################################
############################################################################
def mapHABs(studyArea,analysisStartYear,analysisEndYear,startMonth,endMonth,clean_stats,stats_sa,reducer = ee.Reducer.percentile([50]),bands = ['bloom2','NDGI'],z_thresh = 2):
  try:
    saBounds = studyArea.geometry().bounds()
  except:
    saBounds = studyArea.bounds()

  startJulian = int(ee.Date.fromYMD(analysisStartYear,startMonth,1).format('DD').getInfo())
  endJulian = int(ee.Date.fromYMD(analysisEndYear,endMonth,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
  print(analysisStartYear,analysisEndYear,startMonth,endMonth, startJulian,endJulian)
  
  dirty_imgs = getImagesLib.getProcessedLandsatAndSentinel2Scenes(
    saBounds,
    analysisStartYear,
    analysisEndYear,
    startJulian,
    endJulian,
    convertToDailyMosaics = True,
    preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,
    preComputedSentinel2CloudScoreOffset=preComputedSentinel2CloudScoreOffset,
    preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,
    preComputedLandsatTDOMIRStdDev=preComputedLandsatTDOMIRStdDev,
    preComputedSentinel2TDOMIRMean=preComputedSentinel2TDOMIRMean,
    preComputedSentinel2TDOMIRStdDev=preComputedSentinel2TDOMIRStdDev)
  dirty_imgs = dirty_imgs.map(getImagesLib.HoCalcAlgorithm2)
  for analysisYear in range(analysisStartYear,analysisEndYear+1):
    for month in range(startMonth,endMonth+1):
      startJulianT = int(ee.Date.fromYMD(analysisYear,month,1).format('DD').getInfo())
      endJulianT = int(ee.Date.fromYMD(analysisYear,month,1).advance(1,'month').advance(-1,'day').format('DD').getInfo())
      dirty_imgsT = dirty_imgs.filter(ee.Filter.calendarRange(analysisYear,analysisYear,'year'))\
                              .filter(ee.Filter.calendarRange(startJulianT,endJulianT))
      # print(dirty_imgsT.size().getInfo(),analysisYear,month,startJulianT,endJulianT)
  
  # # Map.addLayer(studyArea,{},'Study Area')
  # # Map.addLayer(waterMask,{'min':1,'max':1,'palette':'00F'},'Water Mask')

  
 
      Map.addLayer(dirty_imgsT.median(),getImagesLib.vizParamsFalse,'Composite yr{} m{}'.format(analysisYear,month),False)
      dirty_imgsT = dirty_imgsT.select(bands)

      waterMask = getWaterMask(analysisYear,analysisYear,month,month)
      means = ee.Image([clean_stats['{}_{}'.format(stats_sa,month)][band+'_mean'] for band in bands])
      stdDevs = ee.Image([clean_stats['{}_{}'.format(stats_sa,month)][band+'_stdDev'] for band in bands])
      dirty_zs = dirty_imgsT.map(lambda img: img.subtract(means).divide(stdDevs).updateMask(waterMask))
    
      dirty_z = dirty_zs.reduce(reducer).reduce(ee.Reducer.max())
      hab = dirty_z.gte(z_thresh).selfMask()

      Map.addLayer(dirty_z,{'min':0,'max':2,'palette':'00F,0F0'},'Dirty Z yr{} m{}'.format(analysisYear,month),False)
      Map.addLayer(hab,{'min':1,'max':1,'palette':'0FF'},'HAB yr{} m{}'.format(analysisYear,month),False)


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