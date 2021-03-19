
from  algal_lib import *
#####################################################################################
##Define user parameters:
output_table_dir = r'Q:\Algal_detection_GEE_work\Viz_Outputs'
output_table_name = os.path.join(output_table_dir,'test4.json')
baselineStartYear = 2018
baselineEndYear = 2019
analysisStartYear = 2020
analysisEndYear = 2020
zThresh = 1.5
startJulian = 1
endJulian = 365
nSamples = 500
compositePeriod = 16
exportBands = ['bloom2','NDSI','NBR','NDMI','NDVI']
mn_test =ee.Geometry.Polygon(\
        [[[-92.40898254274096, 48.17134997592254],\
          [-92.40898254274096, 47.855327704615014],\
          [-91.66465881227221, 47.855327704615014],\
          [-91.66465881227221, 48.17134997592254]]], None, False)

dirty_billy_chinook = ee.Geometry.Polygon(\
        [[[-121.48239392202757, 44.62035796573114],\
          [-121.48239392202757, 44.509308410535326],\
          [-121.2262751476135, 44.509308410535326],\
          [-121.2262751476135, 44.62035796573114]]], None, False)

studyArea = dirty_billy_chinook
#If available, bring in preComputed cloudScore offsets and TDOM stats
#Set to null if computing on-the-fly is wanted
#These have been pre-computed for all CONUS for Landsat and Setinel 2 (separately)
#and are appropriate to use for any time period within the growing season
#The cloudScore offset is generally some lower percentile of cloudScores on a pixel-wise basis
preComputedCloudScoreOffset = getPrecomputedCloudScoreOffsets(10)
preComputedLandsatCloudScoreOffset = preComputedCloudScoreOffset['landsat']
preComputedSentinel2CloudScoreOffset = preComputedCloudScoreOffset['sentinel2']

#The TDOM stats are the mean and standard deviations of the two bands used in TDOM
#By default, TDOM uses the nir and swir1 bands
preComputedTDOMStats = getPrecomputedTDOMStats()
preComputedLandsatTDOMIRMean = preComputedTDOMStats['landsat']['mean']
preComputedLandsatTDOMIRStdDev = preComputedTDOMStats['landsat']['stdDev']

preComputedSentinel2TDOMIRMean = preComputedTDOMStats['sentinel2']['mean']
preComputedSentinel2TDOMIRStdDev = preComputedTDOMStats['sentinel2']['stdDev']

#####################################################################################
#Function Calls
#Get all images
images = getProcessedLandsatAndSentinel2Scenes(studyArea,\
      baselineStartYear,\
      analysisEndYear,\
      startJulian,\
      endJulian,\
      preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,\
      preComputedSentinel2CloudScoreOffset = preComputedSentinel2CloudScoreOffset,\
      preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,\
      preComputedLandsatTDOMIRStdDev = preComputedLandsatTDOMIRStdDev,\
      preComputedSentinel2TDOMIRMean = preComputedSentinel2TDOMIRMean,\
      preComputedSentinel2TDOMIRStdDev = preComputedSentinel2TDOMIRStdDev
      )
Map.addLayer(images.median(),vizParamsFalse,'Median',False)

#Get water mask
permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
              .filter(ee.Filter.calendarRange(baselineStartYear,analysisEndYear,'year'))\
              .filter(ee.Filter.calendarRange(startJulian,endJulian)).mode().eq(2).unmask(0)

water_mask = permWater.Or(tempWater).selfMask()
Map.addLayer(water_mask, {'min':1, 'max':1, 'palette': "00F"}, 'Water Mask', False)

usfs = ee.FeatureCollection('projects/lcms-292214/assets/CONUS-Ancillary-Data/FS_Boundaries')\
          .filterBounds(studyArea)
usfs = usfs.map(lambda f:f.set('constant',1))
usfsMask = usfs.reduceToImage(['constant'], ee.Reducer.first()).focal_max(2000,'circle','meters')

Map.addLayer(usfsMask,{},'USFS Lands',False)


images = images.map(HoCalcAlgorithm2)
Map.addLayer(images.select(exportBands),{},'ts')

composites = nDayComposites(images,baselineStartYear,analysisEndYear,startJulian,endJulian,compositePeriod)
Map.addLayer(composites.select(exportBands),{},'16 day composites')

stack = composites.select(exportBands).toBands();

#fix band names to be yyyy_startDay_endDay
bns = stack.bandNames()
bns = bns.map(lambda bn: ee.String(bn).split('_').slice(1,None).join('_'))
stack = stack.rename(bns);



waterVector = water_mask.clip(studyArea).reduceToVectors()
Map.addLayer(waterVector,{},'Water Vector')

randomSample = ee.FeatureCollection.randomPoints(waterVector, nSamples, 0, 50)
Map.addLayer(randomSample,{},'Samples',False)

if not os.path.exists(output_table_name):
	tt = threading.Thread(target = getTableWrapper, args = (stack,randomSample,output_table_name))
	tt.start()
	limitThreads(1)



Map.addLayer(studyArea,{},'Study Area')
# Map.centerObject(studyArea)
Map.view()
