"""
   Copyright 2023 Ian Housman, RedCastle Resources Inc.

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

# Template script to map algal blooms in Google Earth Engine using field observations for training 
# Intended to work within the geeViz package
# python -m pip install geeViz
# Also requires numpy and pandas
from HAB_Mapper_5k_v3_Lib import *
Map.port=1233
####################################################################################################
####################################################################################################
# Params
# Paths for input training data

# Water training. Manually interpreted here: https://code.earthengine.google.com/?scriptPath=users%2Faaronkamoske%2FAlgalBlooms%3AWater_Training.js
water_training_points = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\waterTraining.geojson"

# HCB Data
# Accessed at: https://app.smartsheet.com/sheets/PjFmf547rwChWWwPGx4jfFhwmrgjq8Jp46ccHW61
hcb_data = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\HCB Data Sharing.xlsx"

# WDEQ Data
wdeq_data = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\WDEQ_Lake_Phytoplankton_Data_2002-2021.xlsx"

# Centroids of clean water bodies
clean_points = r"Q:\Algal_detection_GEE_work\Supervised_Method\Training_Data\HABControlPoints\HABControlPoints.geojson"

# Time period for extracting clean training data
# These dates are also being used to train the water model
clean_startYear=2018
clean_endYear=2022
clean_startJulian = 213
clean_endJulian = 229

# Number of randomly located samples to draw from each clean water body
clean_nSamples = 150

# HCB spreadsheet params
# Field names for coordinates and dates
hcb_lat='Sampling Latitude'
hcb_lon='Sampling Longitude'
hcb_dateProp = 'Sample Date'

# WDEQ spreadsheet params
# Field names for coordinates and dates
wdeq_lat='Lat_Dec'
wdeq_lon='Long_Dec'
wdeq_dateProp = 'CollDate'

# Properties to keep from the HCB table in the GEE featureCollection
hcb_properties=['Sample Date','Wind Conditions','Cloud Cover','Bloom Description','Bloom Type','Sample Collection Method','Cyanobacteria Count (cells/mL)','Cyanobacteria Biovolume (um^3)']

# Properties to keep from the WDEQ table in the GEE featureCollection
wdeq_properties=['CollDate','Individuals (Raw Cnt)','Density (cells/L)','Class']#,'Genus/Species/Variety','Family','Genus','Class','Order','Division','CollMeth']

# Bands to be used to train/apply the algal models
pred_bands = ['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI','elevation']

# Bands to be used to train/apply the water model
water_pred_bands=['blue', 'green', 'red', 're1', 're2', 're3', 'nir', 'nir2', 'waterVapor', 'cirrus', 'swir1', 'swir2', 'NDVI', 'NBR', 'NDMI', 'NDSI', 'brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth', 'tcAngleBG', 'NDCI', 'tcAngleGW', 'tcAngleBW', 'tcDistBG', 'tcDistGW', 'tcDistBW', 'bloom2', 'NDGI', 'elevation','slope', 'aspect',  'hillshade', 'tpi_29', 'tpi_59']

# Number of trees in RF model 
nTrees = 150

# Whether or not to get variable importance and error info about models (can take some time to run so only run when needed)
getError = True

# Projection info
# Currently using the USGS NLCD grid at 10m spatial res
crs = getImagesLib.common_projections['NLCD_CONUS']['crs']
transform = [10,0,-2361915.0,0,-10,3177735.0]

# Directory to store intermediate model inputs and outputs
output_dir = r'Q:\Algal_detection_GEE_work\Supervised_Method\Outputs'

# Applying model and exporting to asset info
# GEE image collection to export to
export_outputs=False
output_asset_collection = 'projects/gtac-algal-blooms/assets/outputs/HAB-RF-Images'

# What years to apply the model to
applyYears = [2022]#[int(ee.Date(datetime.datetime.today()).format('YYYY').getInfo())]#[2020,2021,2022]

# Automatically find the dates that can be run
mostRecentS2Date = getMostRecentS2Date()
earliestJulian = 150
applyFrequency = 14
applyNDayWindow = 28
applyStartJulians = [210]#list(range(earliestJulian,mostRecentS2Date-applyNDayWindow,applyFrequency))

# Which models to apply
apply_HCB_Model=False
apply_WDEQ_Model=True

# Whether to extract points from a focal mean of a given radius kernel
# Set to None if no focal mean should be computed
hcb_focal_radius = 2.5
wdeq_focal_radius = None

# Give some unique study area name to the outputs
studyAreaName='GYE-WY'

# Specify a study area to model across
# Currently using the union of WY and the GYE
# applyStudyAreaWY = ee.Geometry.Polygon(
#         [[[-111.19609375, 45.04612495487054],
#           [-111.19609375, 41.00037953566339],
#           [-104.076953125, 41.00037953566339],
#           [-104.076953125, 45.04612495487054]]], None, False)
# applyStudyAreaLarge = ee.Geometry.Polygon(
#         [[[-112.57620954084283,46.01879051083097],[-112.4369301748728,42.83567849910429],[-112.22232561927319,41.60142657357238],[-111.94479109005457,41.12094139364396],[-111.85037331522834,40.4124045531744],[-111.71992429071048,39.96012600761873],[-109.41392169381751,40.25547306772276],[-107.20268212854833,40.40797952136478],[-105.08741729424692,40.55629989931333],[-104.86876131220322,41.238072102013746],[-105.04135172293292,42.544038977968626],[-105.96583234921064,43.323851744452966],[-106.11418240351291,44.142756413130016],[-106.47606538791294,44.91707272795046],[-108.4712502176446,45.787698126347934],[-109.78713579601984,46.47829894398842],[-111.44315471144157,46.32029673542297],[-112.57620954084283,46.01879051083097]]], None, False)

states = ee.FeatureCollection("TIGER/2018/States")
wy = states.filter(ee.Filter.eq('STUSPS','WY'))
gya = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/gyaoutline')
applyStudyAreaLarge = wy.geometry().union(gya,500)    

#####################################################
# Wrapper to run it all
if __name__ == '__main__':
    Map.addLayer(applyStudyAreaLarge,{'layerType':'geeVectorImage','strokeColor':'F00'},'Bloom Mapper Study Area')
    supervised_algal_mapper(water_training_points,hcb_data,wdeq_data,clean_points,\
    clean_startYear,clean_endYear,clean_startJulian,clean_endJulian,clean_nSamples,\
    hcb_lat,hcb_lon,hcb_dateProp,hcb_properties,\
    wdeq_lat,wdeq_lon,wdeq_dateProp,wdeq_properties,\
    crs,transform,pred_bands,water_pred_bands,output_dir,nTrees,getError,applyStudyAreaLarge,applyYears,applyStartJulians,applyNDayWindow,export_outputs,output_asset_collection,studyAreaName,hcb_focal_radius,wdeq_focal_radius,apply_HCB_Model,apply_WDEQ_Model)

    # tml.trackTasks2()
#     # tml.batchCancel()
#     # viewExportedOutputs(output_dir,output_asset_collection)
# #####################################################
# # Uncomment this out if map viewing is needed
Map.turnOnInspector()
Map.view()  

#####################################################

