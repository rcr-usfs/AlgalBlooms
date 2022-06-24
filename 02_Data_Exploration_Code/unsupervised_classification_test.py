# Originally written by: Ian Housman
# Goal: To explore unsupervised classification methods of water bodies for potential improvements to algal detection
##############################################################
import os,sys
sys.path.append(os.getcwd())
import study_areas as sas
#Module imports
import geeViz.getImagesLib as getImagesLib
ee = getImagesLib.ee
Map = getImagesLib.Map
Map.clearMap()
##############################################################
# Specify parameters

# Bring in study areas from study_areas.py script
study_areas = sas.study_areas[list(sas.study_areas.keys())[0]]
for sa in list(sas.study_areas.keys())[1:]:
	study_areas.extend(sas.study_areas[sa])
study_areas = ee.Geometry.MultiPolygon(study_areas, None, False)
study_area = study_areas

# Specify which year to look at
year = 2021

# Define cluster parameters
cluster_bands = ['NDGI','blue','green','brightness','greenness','wetness']
cluster_sample_n = 100000
cluster_n = 5
cluster_scale = 30
cluster_projection = getImagesLib.common_projections['NLCD_CONUS']['crs']

# Populate names and colors of each cluster
cluster_lookup = {	0:{'name':'0-Clean',
					'color':'00A'},
					1:{'name':'1-Murky/Silty',
					'color':'880'},
					2:{'name':'2-Shallow',
					'color':'F0F'},
					3:{'name':'3-Murky',
					'color':'DD0'},
					4:{'name':'4-Clean',
					'color':'00D'}
					}

####################################################################################################
# Get the composite (precomputed from LCMS)
composite = ee.ImageCollection('projects/lcms-tcc-shared/assets/Composites/Composite-Collection-yesL7-1984-2020')\
				.filter(ee.Filter.calendarRange(year,year,'year'))\
				.mosaic()\
				.select(['blue','green','red','nir','swir1','swir2'])\
				.divide(10000)

# Add indices
composite = getImagesLib.simpleAddIndices(composite)
composite = getImagesLib.getTasseledCap(composite)
composite = getImagesLib.simpleAddTCAngles(composite)
composite = getImagesLib.HoCalcAlgorithm2(composite)

Map.addLayer(composite,getImagesLib.vizParamsFalse,'Composite',False)

# Get the water mask
waterMask = getImagesLib.simpleWaterMask(composite,2.5).selfMask()
# Map.addLayer(waterMask,{'min':1,'max':1,'palette':'0000FF'},'Water Mask',False)

# Set up clustering
cluster_img = composite.select(cluster_bands).updateMask(waterMask)
training = cluster_img.sample(region= study_area,scale=cluster_scale,projection=cluster_projection,numPixels= cluster_sample_n)

kmeansClust = ee.Clusterer.wekaKMeans(cluster_n).train(training)
kmeansResult = cluster_img.cluster(kmeansClust)

# Set up how to visualize clusters
clusterColors = [cluster_lookup[i]['color'] for i in list(range(cluster_n))]
queryDict = {k: cluster_lookup[k]['name'] for k in cluster_lookup.keys()}
legendDict =  {cluster_lookup[k]['name']: cluster_lookup[k]['color'] for k in cluster_lookup.keys()}
clusterViz = {'min':0,'max':cluster_n-1,'palette':clusterColors	,'classLegendDict':legendDict,'queryDict':queryDict}

Map.addLayer(kmeansResult,clusterViz,'Clusters')


Map.addLayer(cluster_img.select(['NDGI']),{'min':-0.1,'max':0.1,'palette':'00F,888,0D0'},'NDGI',False)
##############################################################
Map.addLayer(study_area,{},'Study Area')
##############################################################
Map.turnOnInspector()
Map.view()