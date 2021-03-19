from geeViz.getImagesLib import *
import threading,time
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

  return image.addBands(bloom2)
###############################################################
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
  def getJdImages(yr,yrImages,start):
    yr = ee.Number(yr).int16()
    start = ee.Number(start).int16()
    end = start.add(compositePeriod-1).int16()
    jdImages = yrImages.filter(ee.Filter.calendarRange(start,end))
    jdImages = fillEmptyCollections(jdImages,dummyImage)
    composite = jdImages.median()
    return composite.set('system:index',ee.Number(yr).format().cat('_').cat(ee.Number(start).format('%03d')).cat('_').cat(ee.Number(end).format('%03d')))
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
def getTableWrapper(image,fc, outputName,reducer = ee.Reducer.first()):
  #get zonal stats and export
  table = image.reduceRegions(\
      fc,\
      reducer, \
      30, \
      'EPSG:5070',\
      None, \
      4)
  t = table.getInfo()
  o = open(outputName,'w')
  o.write(json.dumps(t))
  o.close()

###############################################################
def limitThreads(limit):
  while threading.activeCount() > limit:
    time.sleep(1)
    print(threading.activeCount(),'threads running')