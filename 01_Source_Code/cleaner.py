import ee
ee.Initialize()


hab_summary_table_folder = 'projects/gtac-algal-blooms/assets/outputs/HAB-Summary-Tables'

hab_z_imageCollection = 'projects/gtac-algal-blooms/assets/outputs/HAB-Z-Images'



tables = ee.data.listAssets({'parent':hab_summary_table_folder})

for table in tables['assets']:
  
#   if table['id'].find('HAB_Summary_Table_yr2021_m8')>-1 and table['id'].find('WY_') == -1:
	print('Deleting: ',table['id'])
	ee.data.deleteAsset(table['id'])


# z_imgs = ee.ImageCollection(hab_z_imageCollection)
# for img in z_imgs.aggregate_array('system:index').getInfo():
# 	if img.find('HAB_Z_yr2021_m8')>-1 and img.find('WY_') == -1:
# 		print('Deleting: ',hab_z_imageCollection+'/'+img)
# 		ee.data.deleteAsset(hab_z_imageCollection+'/'+img)


# tasks = ee.data.getTaskList()
# for task in tasks:
#   # print(task)
#   if task['state'] in ['RUNNING','READY']:
#     print('Cancelling: ',task['id'])

#     ee.data.cancelOperation(task['name'])
#   # if task['state'] in ['COMPLETED']:
#     # print('Completed: ',task['description'])