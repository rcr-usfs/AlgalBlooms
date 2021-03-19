import os,json, pdb,glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

######################################################################
output_table_dir = r'Q:\Algal_detection_GEE_work\Viz_Outputs'
output_table_name = os.path.join(output_table_dir,'test4.json')
######################################################################
def convert_to_csv(output_table_name):
  with open(output_table_name) as jf:
    table = json.load(jf)
  years = []
  jds = []
  bands = []
  dates = []
  for feature in table['features']:
    props = feature['properties']
    for prop in list(props.keys()):
      value = props[prop]
      # print(prop,value)
      band = prop.split('_')[-1]
      if band not in bands:bands.append(band)
      jd = '_'.join(prop.split('_')[1:3])
      if jd not in jds:jds.append(jd)

      date = '_'.join(prop.split('_')[0:3])
      if date not in dates:dates.append(date)
      year = prop.split('_')[0]
      if year not in years:years.append(year)

  for band in bands:
    out_table = 'ID,{}\n'.format(','.join(dates))
    print('Parsing:',band)
    for feature in table['features']:
      id = feature['id']
      values = []
      props = feature['properties']
      for prop in list(props.keys()):
        value = str(props[prop])
        if value == 'None':value = ''
        values.append(value)
      out_line = '{},{}\n'.format(id,','.join(values))
      out_table += out_line
    o = open(output_table_name + '_{}.csv'.format(band),'w')
    o.write(out_table)
    o.close()
######################################################################
def hovmuller(table,scaler = 1000,min= -400,max = 1000,bin_step = 200):
  data = pd.read_csv(table)
  columns =data.columns[1:]

  dates = [int(''.join(i.split('_')[0:2])) for i in columns]

  allpercentiles = []
  alldatatables = []

  for column in columns:
    c = data[column]*scaler
    # print(c)
    allpercentiles.append(c.quantile([0.05,0.25,0.5,0.75,0.95]))
    alldatatables.append(pd.cut(c,np.arange(min,max +bin_step,bin_step)).value_counts(sort=False))

  percentiles = pd.concat(allpercentiles, axis=1).transpose()
  allbins = pd.concat(alldatatables, axis =1).transpose()
  allbins_np = allbins.to_numpy()
  bins = np.arange(min,max,bin_step)
  counts = allbins.to_numpy().transpose()
  sumCounts = np.sum(counts, axis = 0)
  normalized = np.multiply(np.divide(counts, sumCounts), 100)

  # ------------------------------------------------
  #     plot - statewide hovmuller by gde
  # ------------------------------------------------
  # Plot
  fig = plt.figure(figsize=(6.75,3),frameon=True,facecolor='w')

  ax = fig.add_axes([0.1, 0.14, 0.76, 0.85])
  cmap = plt.get_cmap('viridis')
  cf = plt.pcolormesh(columns, bins, normalized, cmap = cmap)#, vmin = 500)
  q25 = plt.plot(columns, percentiles[0.05], linestyle = ':', color = 'w', linewidth = 2)
  q50 = plt.plot(columns, percentiles[0.50], linestyle = '-', color = 'k', linewidth = 2)
  q75 = plt.plot(columns, percentiles[0.95], linestyle = ':', color = 'w', linewidth = 2)
  ax.set_ylim([min, max])
  ax.set_ylabel('Index Value', fontsize = 10)
  ax.set_xlabel('Year', fontsize = 10)
  ax.xaxis.set_major_locator(plt.MultipleLocator(5))
  ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
  ax.yaxis.set_major_locator(plt.MultipleLocator(2))
  ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
  ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
  ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
  # ax.set_yticklabels([str(abs(i)) for i in range(min, max, int((max-min)/10))])
  # ax.tick_params(axis='both', which='major', labelsize=10)
  # ax.grid(True, which='major', axis='y', linestyle='--', color='k')
  # ax.text(0.99, 0.97, 'n = '+str(sumCounts[0]), transform=ax.transAxes, 
  #   fontsize = 10, color = 'w', verticalalignment='top', horizontalalignment = 'right')

  # cbax = fig.add_axes([0.88, 0.14, 0.03, 0.85]) 
  # cb = plt.colorbar(cf, cax = cbax, orientation = 'vertical')
  # cb.ax.tick_params(labelsize=10)
  # cb.set_label('Percent of Obs (%)', rotation = 270, labelpad = 15, fontsize = 10) 

  # cbar = fig.colorbar(cf)
  # #cbar.ax.set_yticklabels(['0','1','2','>3'])
  # cbar.set_label('Percent of GDEs (%)', rotation=270, labelpad=20)
  fig.savefig(os.path.splitext(table)[0] + '.png')
  plt.close()
######################################################################
if __name__ == '__main__':
  # convert_to_csv(output_table_name)
  csvs = glob.glob(os.path.join(output_table_dir,'*.csv'))
  for csv in csvs[:1]:
    hovmuller(csv)
######################################################################