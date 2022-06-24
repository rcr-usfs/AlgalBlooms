"""
   Copyright 2022 Ian Housman

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

import os,glob
import study_areas as sas
import geeViz.phEEnoViz  as phEEnoViz
ee = phEEnoViz.ee
Map = phEEnoViz.Map
#####################################################################################
#Exports table of n day median composite values for a sample of point locations
#It then creates a hovmuller-like plot of the time series
#This is a useful tool for exploring the variability of a given band/index across space and time
#####################################################################################
##Define user parameters:
#Define location of outputs
output_table_dir = r'C:\PheenoViz_Outputs'

#Define output table name (no extension needed)
# output_table_name ='WY_Dirty'
#Set up dates
#Years can range from 1984-present
#Julian days can range from 1-365
startYear = 2017
endYear = 2021
startJulian = int(ee.Date.fromYMD(2001,7,15).format('DDD').getInfo())
endJulian = int(ee.Date.fromYMD(2001,9,15).format('DDD').getInfo())


#Number of samples to pull (Generally < 3000 or so will work. 5000 seems to be the max)
nSamples = 5000

#Number of days in each median composite
compositePeriod = 16

#Which programs to include
#Options are Landsat and Sentinel 2
#E.g. ['Landsat','Sentinel2'] will include both Landsat 4-8 and Sentinel 2a and 2b TOA data
#If choosing both Landsat and Sentinel 2, available bands/indices are limited to those that can be computed
#Using bands from each (e.g. Sentinel 2 red edge bands would not be available if using only Landsat or Landsat and Sentinel 2)
programs = ['Sentinel2']

#Bands to export. 
#Landsat and Sentinel2 combined band options include: blue, green, red, nir, swir1, swir2, NDVI, NBR, NDMI, brightness, greenness, wetness, bloom2, NDGI
#Sentinel 2 only band options include: 'cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2', NDVI, NBR, NDMI, brightness, greenness, wetness, bloom2, NDGI, NDCI
#Landsat only band options include: blue, green, red, nir, swir1, swir2, temp, NDVI, NBR, NDMI, brightness, greenness, wetness, bloom2, NDG
exportBands = ['NDGI','blue','green','brightness','greenness','wetness','tcAngleBG']#['NBR','NDVI','brightness','greenness','wetness']#['bloom2','NDGI','NDSI','NBR','NDVI','brightness','greenness','wetness']

#Whether to apply snow mask
maskSnow = True

# How many harmonics to include in fit
# Choices are 1, 2, or 3
# 1 will only include the first harmonic (1 htz) - This is best for fitting traditional seasonality
# 2 will include the first (1 htz) and second (2 htz) harmonics
# 3 will include the first (1 htz), second (2 htz), and third (3 htz) harmonics - This is best to fit to all regular cycles in the time series
howManyHarmonics = 1

#Whether to annotate dates of peaks and troughs of harmonic curve on chart
annotate_harmonic_peaks = True

#Whether to show the study area and samples
showGEEViz = False

#Whether to show charts in an interactive, non-png version as they are being created
showChart = False

#Whether to overwrite already produced png chargs
overwriteCharts = False


##############################################################
#Set up study areas
#Generally copy and paste from the Playground

# {'WY_Clean':wy_clean_lake_combo2}
                    
                    # 'WY_Merky':wy_merky_ocean_big_sandy_eden_lakes,
                    # 'WY_Dark_Bloom':wy_dark_bloom_flaming_gorge_rob_roy_half_moon,
                    # 'WY_Glacial_Till':wy_glacial_till_12}
# study_area_dict = {  'WY_Dark_Bloom':wy_dark_bloom_flaming_gorge_rob_roy_half_moon}
#Can set up any mask or really any polygon you'd like to sample
#Here are some examples

#If you want to sample water, using the JRC water layers works well
startWaterYear = startYear
endWaterYear = endYear
if startYear > 2018: startWaterYear = 2018
if endYear > 2018: endWaterYear = 2018
permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
              .filter(ee.Filter.calendarRange(startWaterYear,endWaterYear,'year'))\
              .filter(ee.Filter.calendarRange(startJulian,endJulian)).mode().eq(2).unmask(0)

water_mask = permWater.Or(tempWater).focal_min(3).selfMask()

# studyArea = dirty_odell_lake
#If you would like to visualize phenology of trees, the LCMS tree layer works well
#LCMS land cover classes are as follows:
# 1: Trees
# 2: Tall Shrubs & Trees Mix (SEAK Only)
# 3: Shrubs & Trees Mix
# 4: Grass/Forb/Herb & Trees Mix
# 5: Barren & Trees Mix
# 6: Tall Shrubs (SEAK Only)
# 7: Shrubs
# 8: Grass/Forb/Herb & Shrubs Mix
# 9: Barren & Shrubs Mix
# 10: Grass/Forb/Herb
# 11: Barren & Grass/Forb/Herb Mix
# 12: Barren or Impervious
# 13: Snow or Ice
# 14: Water
# 15: Non-Processing Area Mask
# lcmsLC = ee.ImageCollection("USFS/GTAC/LCMS/v2020-5").select(['Land_Cover']).mode()
# lcmsTreeMask = lcmsLC.eq(1).selfMask()
# studyArea = lcmsTreeMask.clip(uinta_tree).reduceToVectors(scale = 30)

# #Can also look at other land cover classes
# #This will look at grasses and shrubs
# lcmsShrubMask = lcmsLC.eq(7).Or(lcmsLC.eq(10)).selfMask()
# studyArea = lcmsShrubMask.clip(wy_shrub).reduceToVectors(scale = 90)

# #Or you could use the NLCD Tree Canopy Cover layer to get a tree mask
# nlcdTCC =  ee.ImageCollection("USGS/NLCD_RELEASES/2016_REL")\
#               .filter(ee.Filter.calendarRange(2016,2016,'year'))\
#               .select(['percent_tree_cover']).mosaic()
# tccTreeMask = nlcdTCC.gte(30).selfMask()
# studyArea = tccTreeMask.clip(ga_test).reduceToVectors(scale = 30)


# #You can also just pass a feature collection (keeping the area relatively small helps ensure it will run)
# states = ee.FeatureCollection("TIGER/2018/States")
# studyArea = states.filter(ee.Filter.eq('STUSPS','VI'))
#Or also apply a mask to it as well
# studyArea = tccTreeMask.clip(states.filter(ee.Filter.inList('STUSPS',['PR','VI'])).geometry()).reduceToVectors(scale = 300)

# studyArea = ga_test#states.filter(ee.Filter.inList('STUSPS',['PR','VI']))
######################################################################################
#Main function calls
if __name__ == '__main__':
  for output_table_name in list(sas.study_areas.keys()):
    study_area_outline = ee.Geometry.MultiPolygon(sas.study_areas[output_table_name],None,False)
    studyArea = water_mask.clip(study_area_outline).reduceToVectors(scale = 30)

    #Set up output directories
    table_dir = os.path.join(output_table_dir,output_table_name,'tables')
    chart_dir = os.path.join(output_table_dir,output_table_name)
    
    #Get raw json table of samples of time series across provided area
    if not os.path.exists(output_table_name):
      phEEnoViz.getTimeSeriesSample(startYear,endYear,startJulian,endJulian,compositePeriod,exportBands,studyArea,nSamples,os.path.join(table_dir,output_table_name+'.json',),showGEEViz = showGEEViz,maskSnow = maskSnow,programs = programs)


    #Create plots for each band
    #Get a list of csvs for the specified parameters
    csvs = glob.glob(os.path.join(table_dir,'*{}*_{}-{}_{}_{}*.csv'.format('-'.join(programs),startJulian,endJulian,compositePeriod,nSamples)))
    csvs = [i for i in csvs if int(os.path.splitext(os.path.basename(i))[0].split('_')[-5]) in range(startYear,endYear+1)]
    # print(csvs)
    #Create plots
    phEEnoViz.chartTimeSeriesDistributions(csvs,chart_dir,output_table_name + '_{}_{}-{}_{}-{}_{}_{}'.format('-'.join(programs),startYear,endYear,startJulian,endJulian,compositePeriod,nSamples),overwrite = overwriteCharts,howManyHarmonics = howManyHarmonics,showChart =showChart,annotate_harmonic_peaks = annotate_harmonic_peaks,min_pctl = 0.5,max_pctl = 99.5)