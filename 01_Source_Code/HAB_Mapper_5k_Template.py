"""
   Copyright 2021 Ian Housman, RedCastle Resources Inc.

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
#Script to map algal blooms in Google Earth Engine
#Intended to work within the geeViz package
# python -m pip install geeViz
#Also requires gdal
####################################################################################################
from HAB_Mapper_5k_Lib import *
##################################################
#Specify years to compute statistis (mean and stdDev) for clean lakes
cleanStartYear = 2019
cleanEndYear = 2020

#Specify years to map HABs
analysisStartYear = 2021
analysisEndYear = 2021

#Specify which months to map HABs for
startMonth = 5
endMonth = 9

#Location of clean lakes statistics - local path
stats_json = r'Q:\Algal_detection_GEE_work\HAB_Mapper_5k_Outputs\clean_stats.json'


#Export Params
exportZAndTables = True

#GEE Asset folder location for HAB percent summary tables
hab_summary_table_folder = 'projects/gtac-algal-blooms/assets/outputs/HAB-Summary-Tables'
local_hab_summary_table_folder = r'Q:\Algal_detection_GEE_work\HAB_Mapper_5k_Outputs\summary_tables'
#GEE Asset imageCollection location for HAB z score images
hab_z_imageCollection = 'projects/gtac-algal-blooms/assets/outputs/HAB-Z-Images'

#Export and summary projection info
crs = 'EPSG:5070'
transform = [30,0,-2361915.0,0,-30,3177735.0]

#Z score reducer - how to summarize the multiple z scores for a given month
#If a more conservative depiction of algal blooms is needed, put a lower percentile and visa versa
reducer = ee.Reducer.percentile([50])

#Z score threshold for identifying HABs (generally 1-3 or so works well) - anything above this z-score is mapped as HAB
z_thresh = 1

#Indices to use 
#Good choices are NDGI or bloom2
hab_indices = ['NDGI']

#Specify study areas to map
run_study_area_keys = ['OR','WA','WY']

#Define clean control stats
clean_lakes = {}

clean_lakes['OR'] = ee.Geometry.MultiPolygon(\
        [[[[-122.16247277538403, 42.98060571118538],\
           [-122.16247277538403, 42.9016878744337],\
           [-122.0512362031184, 42.9016878744337],\
           [-122.0512362031184, 42.98060571118538]]],\
         [[[-122.03353176573546, 43.50210081409495],\
           [-122.03353176573546, 43.45925265588565],\
           [-121.95216427305968, 43.45925265588565],\
           [-121.95216427305968, 43.50210081409495]]],\
         [[[-122.06863020193826, 43.76538659547769],\
           [-122.06863020193826, 43.6812758137654],\
           [-122.00648878348123, 43.6812758137654],\
           [-122.00648878348123, 43.76538659547769]]]], None, False) 
clean_lakes['WY'] = ee.Geometry.MultiPolygon(\
        [[[[-109.56790010304167, 42.96292917782831],\
           [-109.56790010304167, 42.953161377243354],\
           [-109.55584089131071, 42.953161377243354],\
           [-109.55584089131071, 42.96292917782831]]],\
         [[[-109.70411340565397, 42.94700466781123],\
           [-109.70411340565397, 42.92224596795357],\
           [-109.66875116200163, 42.92224596795357],\
           [-109.66875116200163, 42.94700466781123]]],\
         [[[-109.90752207668908, 43.02675227419442],\
           [-109.90752207668908, 42.98683332211115],\
           [-109.82924448879845, 42.98683332211115],\
           [-109.82924448879845, 43.02675227419442]]]], None, False)

clean_lakes['WA'] = ee.Geometry.MultiPolygon(\
        [[[[-121.32593536613808, 47.586467669618045],\
           [-121.32593536613808, 47.56342205752472],\
           [-121.27993011711465, 47.56342205752472],\
           [-121.27993011711465, 47.586467669618045]]],\
         [[[-121.39322662590371, 47.59920209760001],\
           [-121.39322662590371, 47.57685691359207],\
           [-121.37297058342324, 47.57685691359207],\
           [-121.37297058342324, 47.59920209760001]]],\
         [[[-121.33949661491738, 47.60695700338792],\
           [-121.33949661491738, 47.59723425208618],\
           [-121.32748031853066, 47.59723425208618],\
           [-121.32748031853066, 47.60695700338792]]]], None, False)

#Bring in summary areas
summary_areas_dict = {}
summary_areas_dict['WA'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WA_FS_Named_Recreation_Lakes_v2')
summary_areas_dict['OR'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/OR_FS_Named_Recreation_Lakes_v2')
summary_areas_dict['WY'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WY_Named_FS_Lakes')

#Set up study area
# states = ee.FeatureCollection("TIGER/2018/States")
# study_area = states.filter(ee.Filter.eq('NAME',study_area_state)).first().geometry()
# study_area = ee.Geometry.Polygon(
#         [[[-110.0178203886775, 43.26704900177293],
#           [-110.0178203886775, 42.793226687150714],
#           [-109.419065505865, 42.793226687150714],
#           [-109.419065505865, 43.26704900177293]]], None, False)
# summary_areas = summary_areas_dict[study_area_stats_key]
############################################################################
#Produce clean lake stats, HAB tables, and z score images
#Get clean lake stats
#Will compute them if they do not exist
# clean_stats = getStats(clean_lakes,cleanStartYear,cleanEndYear,startMonth,endMonth,stats_json,hab_indices)

#Map habs
# tasks = batchMapHABs(summary_areas_dict,run_study_area_keys,analysisStartYear,analysisEndYear,startMonth,endMonth,clean_stats,reducer,hab_indices,z_thresh,exportZAndTables,hab_summary_table_folder,hab_z_imageCollection,crs,transform)

#Wait until tasks are finished
# tml.trackTasks2(id_list = tasks)

#####################################################
#Summarize tables to produce deliverables
#Use this to create individual state summaries
batchSummarizeTables(hab_summary_table_folder,local_hab_summary_table_folder,analysisStartYear,analysisEndYear,startMonth,endMonth,summary_areas_dict,run_study_area_keys)

#Use this to summarize all specified states at once
summarizeTables(hab_summary_table_folder,local_hab_summary_table_folder,analysisStartYear,analysisEndYear,startMonth,endMonth,summary_areas_dict,run_study_area_keys)

#####################################################
if  not exportZAndTables:
  Map.view()  

#####################################################
#Make table outputs public if needed (for public viewer)
# makeTablesPublic(hab_summary_table_folder)

