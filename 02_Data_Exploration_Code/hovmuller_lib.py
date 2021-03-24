import os,json, pdb,glob,math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
# from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
# from scipy import stats

# from scipy.interpolate import UnivariateSpline
###############################################################

###############################################################
#Function to convert json gee table into csvs
#Assumes id format is bandName_yyyy-dd-mm
def convert_to_csv(output_table_name):
  with open(output_table_name) as jf:
    table = json.load(jf)
  
  #First find the bands and dates in the table
  bands = []
  dates = []
  print('Finding dates and bands in json')
  for feature in table['features'][:1]:
    props = feature['properties']
    for prop in list(props.keys()):
      value = props[prop]
      band = prop.split('_')[-1]
      if band not in bands:bands.append(band)

      date = prop.split('_')[0]
      if date not in dates:dates.append(date)
      
  #For each of the bands, create a csv
  for band in bands:
    output_csv = output_table_name + '_{}.csv'.format(band)
    out_table = 'ID,{}\n'.format(','.join(dates))
    if not os.path.exists(output_csv):
      print('Parsing:',band)

      #Iterate across each feature and pull all properties and values that have that band
      for feature in table['features']:
        id = feature['id']
        values = []
        props = feature['properties']
        prop_keys = list(props.keys())
        prop_keys = [i for i in prop_keys if i.split('_')[-1] == band]
        
        for prop in prop_keys:
          value = str(props[prop])
          if value == 'None':value = ''
          values.append(value)
        out_line = '{},{}\n'.format(id,','.join(values))
        out_table += out_line

      #Write out table
      o = open(output_csv,'w')
      o.write(out_table)
      o.close()

###############################################################
def hovmuller(table,bin_step = 0.01,background_color = '#1B1716',font_color = '#D6D1CA'):
  print('Creating hovmuller for:',table)
  #Read in data
  data = pd.read_csv(table)
  # data = data.head(50)

  #Find the name of the band/index
  index_name = os.path.splitext(os.path.basename(table))[0].split('_')[-1]

  #Get rid of the id column
  columns =data.columns[1:]

  #Convert values to numpy array
  values = data[columns].to_numpy()

  #Get min and max
  flat = values.flatten()
  flat = flat[~(np.isnan(flat))]
  min = np.percentile(flat,0.05)#flat.min()
  max = np.percentile(flat,99.95)#flat.max()
  values = values.clip(min,max)
  #Get dates (yyyy-mm-dd)
  dates = [i.split('_')[0] for i in columns]

  #Set up bins for histogram
  bins = np.arange(min,max+bin_step,bin_step)

  #Get histograms for each date and clip out outlier frequency values
  hist = np.array([np.histogram(data[column], bins=bins,density = True)[0] for column in columns]).transpose()
  hist = np.nan_to_num(hist, nan=0)
  hist = hist.clip(np.percentile(hist,10),np.percentile(hist,99))
  


  #Still under construction
  #Trying to fit a line to the time series
  table_xs = np.array([])
  table_ys = np.array([])
  table_all_xs = np.array([])
  d0 = date(1970, 1, 1)

  for i,column in enumerate(columns):
    d = dates[i]
    d1 = date(int(d.split('-')[0]),int(d.split('-')[1]),int(d.split('-')[2]))
    delta = d1 - d0
    delta_fraction = math.modf(delta.days/365.25)[0]
    decimal_date = int(d.split('-')[0])+delta_fraction

    ys = data[column].to_numpy()
    ys = ys[~(np.isnan(ys))]
  #   # print(ys)
    xs = np.repeat(decimal_date,len(ys))
    table_ys = np.append(table_ys,ys)
    table_xs = np.append(table_xs,xs)
    table_all_xs = np.append(table_all_xs,decimal_date)
  table_all_xs = np.array(table_all_xs).flatten()
  table_xs = np.array(table_xs).flatten()
  table_ys = np.array(table_ys).flatten()
  # print(len(table_xs),len(table_ys),table_xs,table_ys)

  # slope, intercept, r_value, p_value, std_err = stats.linregress(table_xs,table_ys)

  xs = np.array([table_xs]).T
  sin1Term = np.sin(xs*2*math.pi)
  cos1Term = np.cos(xs*2*math.pi)
  sin2Term = np.sin(xs*4*math.pi)
  cos2Term = np.cos(xs*4*math.pi)
  sin3Term = np.sin(xs*6*math.pi)
  cos3Term = np.cos(xs*6*math.pi)
  intTerm = np.ones(xs.shape[0])
  harm_1 = np.c_[sin1Term,cos1Term,xs, intTerm] 
  harm_1_2 = np.c_[sin1Term,cos1Term,sin2Term,cos2Term,xs, intTerm]
  harm_1_2_3 = np.c_[sin1Term,cos1Term,sin2Term,cos2Term,sin3Term,cos3Term,xs, intTerm]

  harm_1_model = np.linalg.lstsq(harm_1, table_ys, rcond=None)
  harm_1_2_model = np.linalg.lstsq(harm_1_2, table_ys, rcond=None)
  harm_1_2_3_model = np.linalg.lstsq(harm_1_2_3, table_ys, rcond=None)
  # print(beta_hat)
  xs = np.array([table_all_xs]).T
  sin1Term = np.sin(xs*2*math.pi)
  cos1Term = np.cos(xs*2*math.pi)
  sin2Term = np.sin(xs*4*math.pi)
  cos2Term = np.cos(xs*4*math.pi)
  sin3Term = np.sin(xs*6*math.pi)
  cos3Term = np.cos(xs*6*math.pi)
  intTerm = np.ones(xs.shape[0])
  harm_1 = np.c_[sin1Term,cos1Term,xs, intTerm] 
  harm_1_2 = np.c_[sin1Term,cos1Term,sin2Term,cos2Term,xs, intTerm]
  harm_1_2_3 = np.c_[sin1Term,cos1Term,sin2Term,cos2Term,sin3Term,cos3Term,xs, intTerm]
  # print('beta hat:',beta_hat[0],xs)
  pred_1 = np.dot(harm_1,harm_1_model[0])
  pred_1_2 = np.dot(harm_1_2,harm_1_2_model[0])
  pred_1_2_3 = np.dot(harm_1_2_3,harm_1_2_3_model[0])
  
  # pred = np.sum(xs*beta_hat[0],axis = 1)


  # predicted = (table_all_dates_i*slope) + intercept

  # spl = UnivariateSpline(table_xs, table_ys,s  =  len(table_ys) * np.var(table_ys))
  # spl.set_smoothing_factor(0.1)
  # predicted_spline = spl(table_all_dates_i)
  # print(predicted_spline)
  # print(predicted,len(predicted),len(dates))

  #Plot
  fig = plt.figure(figsize=(9.75,6),frameon=True,facecolor='w')

  fig.patch.set_facecolor(background_color)
  
  params = {"ytick.color" : font_color,
          "xtick.color" : font_color,
          "axes.labelcolor" : font_color,
          "axes.edgecolor" : font_color}
  plt.rcParams.update(params)

  ax = fig.add_axes([0.1, 0.14, 0.76, 0.85])

  cmap = plt.get_cmap('viridis')
  cf = plt.pcolormesh(dates, bins, hist, cmap = cmap,shading='flat')#, vmin = 500)
  degrees = 45
  plt.xticks(rotation=degrees, fontsize = 8)


  # q25 = plt.plot(dates, percentiles[0.05], linestyle = ':', color = 'w', linewidth = 2)
  # q50 = plt.plot(dates, pred_1, linestyle = '--', color = '#FF0000', linewidth = 2)
  # q50 = plt.plot(dates, pred_1_2, linestyle = '--', color = '#FF8800', linewidth = 2)
  q50 = plt.plot(dates, pred_1_2_3, linestyle = '--', color = '#FFFFFF', linewidth = 1)
  # q75 = plt.plot(dates, spl(table_all_dates_i), linestyle = ':', color = 'w', linewidth = 2)
  # for i in np.arange(0,1,0.1):
  #   spl.set_smoothing_factor(i)
  #   q75 = plt.plot(dates, spl(table_all_dates_i), linestyle = ':', color = 'w', linewidth = 2)
  ax.set_ylim([min, max])
  ax.set_ylabel('{} Value'.format(index_name), fontsize = 10)
  ax.set_xlabel('Date', fontsize = 10)
  ax.xaxis.set_major_locator(plt.MultipleLocator(5))
  ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
  ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
  ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
  
  ax.grid(True, which='major', axis='y', linestyle='--', color='k')
  ax.grid(True, which='major', axis='x', linestyle='--', color='k')
  ax.text(0.99, 0.99, index_name, transform=ax.transAxes, 
    fontsize = 15, color = font_color, verticalalignment='top', horizontalalignment = 'right')

  cbax = fig.add_axes([0.88, 0.14, 0.03, 0.85]) 
  cb = plt.colorbar(cf, cax = cbax, orientation = 'vertical')
  cb.ax.tick_params(labelsize=10)
  cb.set_label('Percent of Obs (%)', rotation = 270, labelpad = 15, fontsize = 10) 


  fig.savefig(os.path.splitext(table)[0] + '.png')
  # plt.show()
  plt.close()