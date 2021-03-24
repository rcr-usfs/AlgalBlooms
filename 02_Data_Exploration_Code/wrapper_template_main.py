from pheno_lib import *
from hovmuller_lib import *
#####################################################################################
#Exports table of n day median composite values for a sample of point locations
#####################################################################################
##Define user parameters:
#Define location of outputs
output_table_dir = r'Q:\Algal_detection_GEE_work\Viz_Outputs'

#Define output table name
output_table_name = os.path.join(output_table_dir,'Chinook_Test.json')

#Set up dates
startYear = 2017
endYear = 2020
startJulian = 1
endJulian = 365

#Number of samples to pull (Generally < 3000 or so will work)
nSamples = 1000

#Number of days in each median composite
compositePeriod = 8

#Bands to export. Can include: blue, green, red, nir, swir1, swir2, NDVI, NBR, NDMI, brightness, greenness, wetness, bloom2, NDGI
exportBands = ['bloom2','NDGI','NDSI','NBR','NDVI']


#Set up study areas
#Generally copy and paste from the Playground
mn_test =ee.Geometry.Polygon(\
        [[[-92.40898254274096, 48.17134997592254],\
          [-92.40898254274096, 47.855327704615014],\
          [-91.66465881227221, 47.855327704615014],\
          [-91.66465881227221, 48.17134997592254]]], None, False)

dirty_billy_chinook = ee.Geometry.Polygon(\
        [[[-121.48239392202757, 44.62035796573114],\
          [-121.48239392202757, 44.509308410535326],\
          [-121.2262751476135, 44.509308410535326],\
          [-121.2262751476135, 44.62035796573114]]], None, False)

studyArea = dirty_billy_chinook


#Get water mask
#Can set up any mask or really any polygon however you'd like to sample
permWater = ee.Image("JRC/GSW1_1/GlobalSurfaceWater").select([0]).gte(90).unmask(0)
tempWater =ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")\
              .filter(ee.Filter.calendarRange(startYear,endYear,'year'))\
              .filter(ee.Filter.calendarRange(startJulian,endJulian)).mode().eq(2).unmask(0)

water_mask = permWater.Or(tempWater).selfMask()
# Map.addLayer(water_mask, {'min':1, 'max':1, 'palette': "00F"}, 'Water Mask', False)


#Pull out water polygons and set to new study area
studyArea = water_mask.clip(studyArea).reduceToVectors()
# Map.addLayer(studyArea,{},'Water Vector')

######################################################################################
if __name__ == '__main__':
  if not os.path.exists(output_table_name):
    getTimeSeriesSample(startYear,endYear,startJulian,endJulian,compositePeriod,exportBands,studyArea,nSamples,output_table_name)

  #Convert json to csvs by band
  convert_to_csv(output_table_name)

  #Create plots for each table
  csvs = glob.glob(os.path.join(output_table_dir,'*.csv'))
  for csv in csvs:
    hovmuller(csv)