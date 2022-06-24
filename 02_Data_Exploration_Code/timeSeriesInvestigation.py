import os,sys
sys.path.append(os.getcwd())
import study_areas as sas
#Module imports
import geeViz.getImagesLib as getImagesLib
import geeViz.taskManagerLib as taskManagerLib
ee = getImagesLib.ee
Map = getImagesLib.Map
Map.clearMap()
##############################################################
study_areas = sas.study_areas[list(sas.study_areas.keys())[0]]
for sa in sas.study_areas.keys():
	study_areas.extend(sas.study_areas[sa])
study_areas = ee.Geometry.MultiPolygon(study_areas, None, False)
study_area = study_areas
startYear = 2021
endYear = 2021
startJulian = int(ee.Date.fromYMD(2001,7,15).format('DDD').getInfo())
endJulian = int(ee.Date.fromYMD(2001,9,15).format('DDD').getInfo())

####################################################################################################
#Harmonic regression parameters

#Which harmonics to include
#Is a list of numbers of the n PI per year
#Typical assumption of 1 cycle/yr would be [2] (2*pi)
#If trying to overfit, or expected bimodal phenology try adding a higher frequency as well
#ex. [2,4]
whichHarmonics = [2]

#Which bands/indices to run harmonic regression across
indexNames =['NDCI','NDGI'];#,'NBR','NDMI','nir','swir1','swir2','tcAngleBG'];//['nir','swir1','swir2','NDMI','NDVI','NBR','tcAngleBG'];//['blue','green','red','nir','swir1','swir2','NDMI','NDVI','NBR','tcAngleBG'];

#Choose which band/index to use for visualizing seasonality in hue, saturation, value color space (generally NDVI works best)
seasonalityVizIndexName = 'NDCI'


#Whether to apply a linear detrending of data.  Can be useful if long-term change is not of interest
detrend = True
####################################################################################################
#The TDOM stats are the mean and standard deviations of the two bands used in TDOM
#By default, TDOM uses the nir and swir1 bands
# preComputedTDOMStats = getImagesLib.getPrecomputedTDOMStats()
# preComputedLandsatTDOMIRMean = preComputedTDOMStats['landsat']['mean']
# preComputedLandsatTDOMIRStdDev = preComputedTDOMStats['landsat']['stdDev']

# preComputedSentinel2TDOMIRMean = preComputedTDOMStats['sentinel2']['mean']
# preComputedSentinel2TDOMIRStdDev = preComputedTDOMStats['sentinel2']['stdDev']


# allScenes = getImagesLib.getProcessedSentinel2Scenes(\
#   study_area,
#   startYear,
#   endYear,
#   startJulian,
#   endJulian,
#   applyQABand = False,
#   applyCloudScore = False,
#   applyShadowShift = False,
#   applyTDOM = True,
#   cloudScoreThresh = 20,
#   performCloudScoreOffset = True,
#   cloudScorePctl = 10,
#   cloudHeights = ee.List.sequence(500,10000,500),
#   zScoreThresh = -1,
#   shadowSumThresh = 0.35,
#   contractPixels = 1.5,
#   dilatePixels = 3.5,
#   shadowSumBands = ['nir','swir1'],
#   resampleMethod = 'near',
#   toaOrSR = 'SR',
#   convertToDailyMosaics = False,
#   applyCloudProbability = True,
#   preComputedCloudScoreOffset = None,
#   preComputedTDOMIRMean = preComputedSentinel2TDOMIRMean,
#   preComputedTDOMIRStdDev = preComputedSentinel2TDOMIRStdDev,
#   cloudProbThresh = 40)
# allScenes = allScenes.map(getImagesLib.HoCalcAlgorithm2)
# composite =allScenes.median()
# waterMask = getImagesLib.simpleWaterMask(composite).selfMask()
# Map.addLayer(waterMask,{'min':1,'max':1,'palette':'0000FF'},'Water Mask',False)
# seasonalityMedian = composite.select([seasonalityVizIndexName])
# nameStart = str(startYear) + '_'+str(endYear)
# Map.addLayer(composite,getImagesLib.vizParamsFalse,'Median Composite')
# nDayComps = getImagesLib.nDayComposites(allScenes,startYear,endYear,startJulian,endJulian,30)

# vizParams ={'min':-0.1,'max':0.1,'dateFormat':'YYYYMMdd','advanceInterval':'day'}
# # print(nDayComps.aggregate_histogram('system:index').keys().getInfo())
# Map.addTimeLapse(nDayComps.select(indexNames),vizParams,'Weekly Composites')
#Fit harmonic model
# coeffsPredicted =getImagesLib.getHarmonicCoefficientsAndFit(allScenes,indexNames,whichHarmonics,detrend)
# coeffs = coeffsPredicted[0]
# #Get predicted values for visualization
# predicted = coeffsPredicted[1]
# Map.addLayer(predicted,{},nameStart+ '_predicted',False);

# #Optionally simplify coeffs to phase, amplitude, and date of peak
# if 2 in whichHarmonics :
# 	pap = ee.Image(getImagesLib.getPhaseAmplitudePeak(coeffs))


# 	vals = coeffs.select(['.*_intercept'])
# 	amplitudes = pap.select(['.*_amplitude'])
# 	phases = pap.select(['.*_phase'])
# 	peakJulians = pap.select(['.*peakJulianDay'])
# 	AUCs = pap.select(['.*AUC'])

# 	Map.addLayer(phases,{},nameStart+ '_phases',False)
# 	Map.addLayer(amplitudes,{'min':0,'max':0.6},nameStart+ '_amplitudes',False)
# 	Map.addLayer(AUCs,{'min':0,'max':0.3},nameStart+ '_AUCs',False)
# 	Map.addLayer(peakJulians,{'min':0,'max':365},nameStart+ '_peakJulians',False)

# 	#Create synthetic image for peak julian day according the the seasonalityVizIndexName band
# 	dateImage = ee.Image(startYear).add(peakJulians.select([seasonalityVizIndexName + '_peakJulianDay']).divide(365))
# 	synth = getImagesLib.synthImage(coeffs,dateImage,indexNames,whichHarmonics,detrend);
# 	Map.addLayer(synth,{'min':0.1,'max':0.4},nameStart + '_Date_of_Max_'+seasonalityVizIndexName+'_Synth_Image',False);

# 	#Turn the HSV data into an RGB image and add it to the map.
# 	seasonality = ee.Image.cat(phases.select([seasonalityVizIndexName+'.*']).clamp(0,1),amplitudes.select([seasonalityVizIndexName+'.*']).unitScale(0,0.5).clamp(0,1),seasonalityMedian.unitScale(0,0.8).clamp(0,1)).hsvToRgb()

# 	Map.addLayer(seasonality, {'min':0,'max':1}, nameStart+ '_'+seasonalityVizIndexName+'_Seasonality',True);

composite = ee.ImageCollection('projects/lcms-tcc-shared/assets/Composites/Composite-Collection-yesL7-1984-2020')\
				.filter(ee.Filter.calendarRange(startYear,endYear,'year'))\
				.mosaic()\
				.select(['blue','green','red','nir','swir1','swir2'])\
				.divide(10000)
composite = getImagesLib.simpleAddIndices(composite)
composite = getImagesLib.getTasseledCap(composite)
composite = getImagesLib.simpleAddTCAngles(composite)
composite = getImagesLib.HoCalcAlgorithm2(composite)

Map.addLayer(composite,getImagesLib.vizParamsFalse,'Composite',False)
waterMask = getImagesLib.simpleWaterMask(composite).selfMask()
Map.addLayer(waterMask,{'min':1,'max':1,'palette':'0000FF'},'Water Mask',False)
cluster_bands = ['NDGI','greenness','wetness','green','blue']
cluster_sample_n = 10000
cluster_n = 5

cluster_img = composite.select(cluster_bands).updateMask(waterMask)
training = cluster_img.sample(region= study_area,scale=20,numPixels= cluster_sample_n);
# Map.addLayer(training,{'layerType':'geeVector'},'Cluster Training')
# print(training.limit(5).getInfo());
kmeansClust = ee.Clusterer.wekaKMeans(cluster_n+1).train(training)
# print(kmeansClust.getInfo())
kmeansResult = cluster_img.cluster(kmeansClust)
clusterColors = '0FF,880,008,0F8,008,00F'.split(',')
queryDict = {'0':'0-Glacial','1':'1-Murky','2':'2-Clean','3':'3-Glacial/Murky','4':'4-Clean','5':'5-Clean'}
legendDict =  {queryDict[str(i)]: clusterColors[i] for i in range(len(queryDict.keys()))}
print(legendDict)
clusterViz = {'min':0,'max':cluster_n,'palette':clusterColors	,'classLegendDict':legendDict,'queryDict':queryDict}

#Add it to the map

Map.addLayer(kmeansResult,clusterViz,'Clusters')
Map.addLayer(cluster_img.select(['NDGI']),{'min':-0.2,'max':0.2,'palette':'00D,DDD,0D0'},'NDGI',False)

Map.addLayer(study_area,{},'Study Area')

Map.turnOnInspector()
Map.view()