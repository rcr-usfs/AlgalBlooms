# Originally written by: Ian Housman
# Goal: To explore unsupervised classification methods of water bodies for potential improvements to algal detection
##############################################################
import os,sys,math
import pandas as pd
sys.path.append(os.getcwd())
import study_areas as sas
#Module imports
import geeViz.getImagesLib as getImagesLib
ee = getImagesLib.ee
Map = getImagesLib.Map
Map.clearMap()
##############################################################
# Specify parameters
obs_csv = r"Q:\Algal_detection_GEE_work\2021-1025_Wyoming-HCB-Investivations_CellDensity-Microcystins.csv"
lat_column = 'Sampling Latitude'
lng_column = 'Sampling Longitude'
filter_column = 'HCB Confirmation'
filter_condition = 'yes'
buffer_size = 30
##############################################################
def csv_to_gee_featureCollection(csv,lat_column,lng_column):
	df = pd.read_csv(csv).dropna(axis=1)
	return ee.FeatureCollection([ee.Feature(ee.Geometry.Point([row[lng_column], row[lat_column]]).buffer(buffer_size,math.ceil(buffer_size/10)).bounds(),row.to_dict()) for i,row in df.iterrows()])

obs = csv_to_gee_featureCollection(obs_csv,lat_column,lng_column)
print(obs.size().getInfo())
Map.addLayer(obs.filter(ee.Filter.eq(filter_column,filter_condition)),{'layerType':'geeVector','strokeColor':'0FF'},'Observations')
##############################################################
Map.turnOnInspector()
Map.view()