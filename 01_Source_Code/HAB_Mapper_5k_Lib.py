import geeViz.getImagesLib as getImagesLib
ee = getImagesLib.ee
Map = getImagesLib.Map
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
def getWaterMask(startYear,endYear,startJulian,endJulian):
  startWaterYear = startYear
  endWaterYear = endYear
  if startYear > 2018: startWaterYear = 2018
  if endYear > 2018: endWaterYear = 2018
  permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
  tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
                .filter(ee.Filter.calendarRange(startWaterYear,endWaterYear,'year'))\
                .filter(ee.Filter.calendarRange(startJulian,endJulian)).mode().eq(2).unmask(0)

  water_mask = permWater.Or(tempWater).selfMask()
  return water_mask
############################################################################
def getStats(studyArea,startYear,endYear,startJulian,endJulian,bands,maxTries = 10):
  try:
    saBounds = studyArea.geometry().bounds()
  except:
    saBounds = studyArea.bounds()

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

  clean_stats = clean_composite.select(bands).reduceRegion(ee.Reducer.mean().combine(ee.Reducer.stdDev(),'',True),studyArea,30,'EPSG:5070',None,True,1e13)
  stats = None
  tryCount = 1
  def getStatsTryer():return clean_stats.getInfo()
  while stats == None and tryCount < maxTries:
    try:
      print('Computing stats. Try number: ',tryCount)
      stats = getStatsTryer()
      print(stats)
    except Exception as e:
      print(e)
  		
    tryCount+=1
  Map.addLayer(clean_composite,getImagesLib.vizParamsFalse,'Clean Lakes Composite')
	
############################################################################
precomputed_stats = {'WA':{'NDGI_mean': -0.2892009068943467, 'NDGI_stdDev': 0.04485800075759838, 'bloom2_mean': 0.5522615758740659, 'bloom2_stdDev': 0.05982976634202523},
'WY':{'NDGI_mean': -0.21767293738121518, 'NDGI_stdDev': 0.058372681366962116, 'bloom2_mean': 0.6454459387948787, 'bloom2_stdDev': 0.08495050203230246},
'OR':{'NDGI_mean': -0.3231765959110561, 'NDGI_stdDev': 0.06533652226615914, 'bloom2_mean': 0.5149265225223071, 'bloom2_stdDev': 0.08665382662586683}}
############################################################################
def mapHABs(studyArea,analysisYear,startJulian,endJulian,reducer = ee.Reducer.percentile([50]),bands = ['bloom2','NDGI'],clean_stats = precomputed_stats['WA']):
  try:
    saBounds = studyArea.geometry().bounds()
  except:
    saBounds = studyArea.bounds()

  waterMask = getWaterMask(analysisYear,analysisYear,startJulian,endJulian)
  Map.addLayer(studyArea,{},'Study Area')
  Map.addLayer(waterMask,{'min':1,'max':1,'palette':'00F'},'Water Mask')

  dirty_imgs = getImagesLib.getProcessedLandsatAndSentinel2Scenes(
    saBounds,
    analysisYear,
    analysisYear,
    startJulian,
    endJulian,
    preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,
    preComputedSentinel2CloudScoreOffset=preComputedSentinel2CloudScoreOffset,
    preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,
    preComputedLandsatTDOMIRStdDev=preComputedLandsatTDOMIRStdDev,
    preComputedSentinel2TDOMIRMean=preComputedSentinel2TDOMIRMean,
    preComputedSentinel2TDOMIRStdDev=preComputedSentinel2TDOMIRStdDev)
  dirty_imgs = dirty_imgs.map(getImagesLib.HoCalcAlgorithm2)
  Map.addLayer(dirty_imgs.median(),getImagesLib.vizParamsFalse,'Dirty Lakes Composite')
  dirty_imgs = dirty_imgs.select(bands)

  means = ee.Image([clean_stats[band+'_mean'] for band in bands])
  stdDevs = ee.Image([clean_stats[band+'_stdDev'] for band in bands])
  dirty_zs = dirty_imgs.map(lambda img: img.subtract(means).divide(stdDevs).updateMask(waterMask))
 
  Map.addLayer(dirty_zs.reduce(reducer).select([0]),{'min':0,'max':3,'palette':'00F,0F0'},'Dirty Z')
  # clean_stats_img = ee.Image(clean_stats.keys())