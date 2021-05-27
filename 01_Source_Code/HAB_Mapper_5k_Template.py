from HAB_Mapper_5k_Lib import *
##################################################
cleanStartYear = 2020
cleanEndYear = 2020
startMonth = 5
endMonth =9

stats_json = r'Q:\Algal_detection_GEE_work\HAB_Mapper_5k_Outputs\clean_stats.json'
analysisYears = [2020]

reducer = ee.Reducer.percentile([50])


hab_indices = ['bloom2','NDGI']

study_area_state = 'Oregon'
study_area_stats_key = 'OR'

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
         [[[-108.95214927727804, 42.591368975972856],\
           [-108.95214927727804, 42.57948841139656],\
           [-108.92897499138937, 42.57948841139656],\
           [-108.92897499138937, 42.591368975972856]]]], None, False)

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
summary_areas = {}
summary_areas['WA'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WA_FS_Named_Recreation_Lakes_v2')
summary_areas['OR'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/OR_FS_Named_Recreation_Lakes_v2')
summary_areas['WY'] = ee.FeatureCollection('projects/gtac-algal-blooms/assets/ancillary/WY_FS_Named_Recreation_Lakes_v2')


states = ee.FeatureCollection("TIGER/2018/States")
study_area = states.filter(ee.Filter.eq('NAME',study_area_state))







############################################################################
clean_stats = getStats(clean_lakes,cleanStartYear,cleanEndYear,startMonth,endMonth,stats_json,hab_indices)
for analysisYear in analysisYears:
  for month in range(startMonth,endMonth+1):
    mapHABs(study_area,analysisYear,month,clean_stats,study_area_stats_key,reducer,hab_indices)

Map.addLayer(summary_areas[study_area_stats_key],{},'Rec Lakes',False)

Map.view()