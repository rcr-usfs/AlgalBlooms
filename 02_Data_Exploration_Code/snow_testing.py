from geeViz.getImagesLib import *
#####################################################################################
AOI = ee.FeatureCollection("projects/gtac-cheatgrass/assets/Mangum_Fire_1kmBuffer")
AOIBounds = AOI.geometry().bounds()

#Set up dates
startYear = 2020
endYear = 2020
startJulian = 1
endJulian = 30

studyArea = AOI


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
      startYear,\
      endYear,\
      startJulian,\
      endJulian,\
      toaOrSR = 'TOA',
      includeSLCOffL7 = False,
      defringeL5 = False,
      applyQABand = False,
      applyCloudProbability = True,
      applyShadowShift = False,
      applyCloudScoreLandsat = False,
      applyCloudScoreSentinel2 = False,
      applyTDOMLandsat = True,
      applyTDOMSentinel2 = True,
      applyFmaskCloudMask = True,
      applyFmaskCloudShadowMask = True,
      applyFmaskSnowMask = False,
      cloudHeights = ee.List.sequence(500,10000,500),
      cloudScoreThresh = 20,
      performCloudScoreOffset = True,
      cloudScorePctl = 10,
      zScoreThresh = -1,
      shadowSumThresh = 0.35,
      contractPixels = 1.5,
      dilatePixels = 3.5,
      shadowSumBands = ['nir','swir1'],
      landsatResampleMethod = 'near',
      sentinel2ResampleMethod = 'aggregate',
      convertToDailyMosaics = True,
      runChastainHarmonization = True,
      preComputedLandsatCloudScoreOffset = preComputedLandsatCloudScoreOffset,\
      preComputedSentinel2CloudScoreOffset = preComputedSentinel2CloudScoreOffset,\
      preComputedLandsatTDOMIRMean = preComputedLandsatTDOMIRMean,\
      preComputedLandsatTDOMIRStdDev = preComputedLandsatTDOMIRStdDev,\
      preComputedSentinel2TDOMIRMean = preComputedSentinel2TDOMIRMean,\
      preComputedSentinel2TDOMIRStdDev = preComputedSentinel2TDOMIRStdDev
      )
#Vizualize the median of the images
print(images.aggregate_histogram('sensor').getInfo())
lImages = images.filter(ee.Filter.stringContains('sensor','LANDS'))
s2Images = images.filter(ee.Filter.stringContains('sensor','Sentinel'))
 
Map.addLayer(s2Images.median(),vizParamsFalse,'S2 Median Before Snow Masking',False)
s2Images = s2Images.map(sentinel2SnowMask)
Map.addLayer(s2Images.median(),vizParamsFalse,'S2 Median After Snow Masking',False)


Map.addLayer(lImages.median(),vizParamsFalse,'Landsat Median Before Snow Masking',False)
lImages = lImages.map(sentinel2SnowMask)
Map.addLayer(lImages.median(),vizParamsFalse,'Landsat Median After Snow Masking',False)

# #Vizualize outputs
Map.addLayer(studyArea,{},'Study Area')
Map.centerObject(studyArea)
Map.view()
