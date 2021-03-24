from geeViz.getImagesLib import *
import threading,time
###############################################################
#Apply bloom detection algorithm
def HoCalcAlgorithm2(image):
  # Algorithm 2 based on: 
  # Matthews, M. (2011) A current review of empirical procedures 
  #  of remote sensing in inland and near-coastal transitional 
  #  waters, International Journal of Remote Sensing, 32:21, 
  #  6855-6899, DOI: 10.1080/01431161.2010.512947
  
  # Apply algorithm 2: B2/B1
  bloom2 = image.select('green')\
              .divide(image.select('blue'))\
              .rename(['bloom2'])
  ndgi = image.normalizedDifference(['green','blue']).rename(['NDGI'])
  return image.addBands(bloom2).addBands(ndgi)

###############################################################
#Function to take images and create a median composite every n days
def nDayComposites(images,startYear,endYear,startJulian,endJulian,compositePeriod):
  
  #create dummy image for with no values
  dummyImage = ee.Image(images.first())

  #convert to composites as defined above
  def getYrImages(yr):
    #take the year of the image
    yr = ee.Number(yr).int16()
    #filter out images for the year
    yrImages = images.filter(ee.Filter.calendarRange(yr,yr,'year'))
  
    #use dummy image to fill in gaps for GEE processing
    yrImages = fillEmptyCollections(yrImages,dummyImage)
    return yrImages

  #Get images for a specified start day
  def getJdImages(yr,yrImages,start):
    yr = ee.Number(yr).int16()
    start = ee.Number(start).int16()
    index = ee.Date.fromYMD(yr,1,1).advance(start.subtract(1),'day').format('yyyy-MM-dd')
    end = start.add(compositePeriod-1).int16()
    jdImages = yrImages.filter(ee.Filter.calendarRange(start,end))
    jdImages = fillEmptyCollections(jdImages,dummyImage)
    composite = jdImages.median()
    return composite.set('system:index',index)

  #Set up wrappers
  def jdWrapper(yr,yrImages):
    return ee.FeatureCollection(ee.List.sequence(startJulian,endJulian,compositePeriod).map(lambda start: getJdImages(yr,yrImages,start)))
  def yrWrapper(yr):
    yrImages = getYrImages(yr)
    return jdWrapper(yr,yrImages)

  composites = ee.FeatureCollection(ee.List.sequence(startYear,endYear).map(lambda yr:yrWrapper(yr)))
  #return the composites as an image collection
  composites = ee.ImageCollection(composites.flatten());

  return composites
###############################################################
#Function to extract a json table of
def getTableWrapper(image,fc, outputName,reducer = ee.Reducer.first(),tryNumber = 1,maxTries = 15):
  #get zonal stats and export
  table = image.reduceRegions(\
      fc,\
      reducer, \
      30, \
      'EPSG:5070',\
      None, \
      4)
  try:
    print('Exporting table:',outputName)
    t = table.getInfo()
    o = open(outputName,'w')
    o.write(json.dumps(t))
    o.close()
  except Exception as e:
    print('Error encountered:',e)
    tryNumber += 1
    if tryNumber < maxTries:
      print('Trying to convert table again. Try number:',tryNumber)
      getTableWrapper(image,fc, outputName,reducer,tryNumber,maxTries)

###############################################################
def limitThreads(limit):
  while threading.activeCount() > limit:
    time.sleep(1)
    print(threading.activeCount(),'threads running')
###############################################################
#Adapted from: https://earth.esa.int/documents/10174/3166008/ESA_Training_Vilnius_07072017_SAR_Optical_Snow_Ice_Exercises.pdf
def sentinel2SnowMask(img,dilatePixels = 3.5):
  ndsi =  img.normalizedDifference(['green', 'swir1'])

  #IF NDSI > 0.40 AND ρ(NIR) > 0.11 THEN snow in open land
  #IF 0.1 < NDSI < 0.4 THEN snow in forest
  snowOpenLand = ndsi.gt(0.4).And(img.select(['nir']).gt(0.11))
  snowForest = ndsi.gt(0.1).And(ndsi.lt(0.4))

  # Map.addLayer(snowOpenLand.selfMask(),{'min':1,'max':1,'palette':'88F'},'Snow Open Land')
  # Map.addLayer(snowForest.selfMask(),{'min':1,'max':1,'palette':'00F'},'Snow Forest')
  #Fractional snow cover (FSC, 0 % - 100% snow) can be detected by the approach of Salomonson
  #and Appel (2004, 2006), which was originally developed for MODIS data:
  #FSC = –0.01 + 1.45 * NDSI
  fsc = ndsi.multiply(1.45).subtract(0.01)
  # Map.addLayer(fsc,{'min':0,'max':1,'palette':'080,008'},'Fractional Snow Cover')
  snowMask = ((snowOpenLand.Or(snowForest)).Not()).focal_min(dilatePixels)
  return img.updateMask(snowMask)
#Wrapper to get a sample of locations for a given area
def getTimeSeriesSample(startYear,endYear,startJulian,endJulian,compositePeriod,exportBands,studyArea,nSamples,output_table_name):
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
  images = getProcessedLandsatAndSentinel2Scenes(studyArea.geometry().bounds(),\
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
  Map.addLayer(images.median(),vizParamsFalse,'Median Composite',False)

  #Mask snow
  images = images.map(sentinel2SnowMask)


  #Add greenness ratio/indices
  images = images.map(HoCalcAlgorithm2)
  Map.addLayer(images.select(exportBands),{},'Raw Time Series')

  #Convert to n day composites
  composites = nDayComposites(images,startYear,endYear,startJulian,endJulian,compositePeriod)
  Map.addLayer(composites.select(exportBands),{},str(compositePeriod) +' day composites')

  #Convert to a stack
  stack = composites.select(exportBands).toBands();

  #Fix band names to be yyyy_mm_dd
  bns = stack.bandNames()
  bns = bns.map(lambda bn: ee.String(bn).split('_').slice(1,None).join('_'))
  stack = stack.rename(bns);


  #Sample the water polygons
  randomSample = ee.FeatureCollection.randomPoints(studyArea, nSamples, 0, 50)
  Map.addLayer(randomSample,{},'Samples',False)

  #Export table
  if not os.path.exists(output_table_name):
    tt = threading.Thread(target = getTableWrapper, args = (stack,randomSample,output_table_name))
    tt.start()
    limitThreads(1)


  #Vizualize outputs
  Map.addLayer(studyArea,{},'Study Area')
  Map.centerObject(studyArea)
  Map.view()
