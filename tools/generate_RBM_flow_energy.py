#!/usr/local/anaconda/bin/python

# This script generates a flow and an energy file (in the format of RBM input)

import numpy as np
import sys
import datetime as dt
import pandas as pd
import xray
import my_functions

cfg = my_functions.read_config(sys.argv[1])  # Read config file

#====================================================#
# Parameter loading from config file
#====================================================#
# [INPUT]
# Routing station file - output from 'prepare_RBM_param.py'
route_station_file = cfg['INPUT']['route_station_file']
# RVIC output nc file (must be "grid" format)
RVIC_output_nc = cfg['INPUT']['RVIC_output_nc']
# VIC output nc file - energy
vic_output_energy_nc = cfg['INPUT']['vic_output_energy_nc']

# [RBM_OPTIONS]
# RBM start and end date (dt.date)
start_date = dt.date(cfg['RBM_OPTIONS']['start_date'][0], \
                         cfg['RBM_OPTIONS']['start_date'][1], \
                         cfg['RBM_OPTIONS']['start_date'][2])
end_date = dt.date(cfg['RBM_OPTIONS']['end_date'][0], \
                         cfg['RBM_OPTIONS']['end_date'][1], \
                         cfg['RBM_OPTIONS']['end_date'][2])
# Same, but string 'YYYY-MM-DD'
start_date_str = '{:4d}-{:02d}-{:02d}'.format(cfg['RBM_OPTIONS']['start_date'][0], \
                         cfg['RBM_OPTIONS']['start_date'][1], \
                         cfg['RBM_OPTIONS']['start_date'][2])
end_date_str = '{:4d}-{:02d}-{:02d}'.format(cfg['RBM_OPTIONS']['end_date'][0], \
                         cfg['RBM_OPTIONS']['end_date'][1], \
                         cfg['RBM_OPTIONS']['end_date'][2])

# [HYDRAULIC_OPTIONS]
# Leopold coefficients (constant throughout basin)
# <a_d> and <b_d>: a and b coefficients in flow depth estimation: z=aQ^b, where Q is flow discharge [cfs], z is flow depth [ft]
a_d = cfg['HYDRAULIC_OPTIONS']['a_d']
b_d = cfg['HYDRAULIC_OPTIONS']['b_d']
# <a_w> and <b_w>: a and b coefficients in channel width estimation: w=aQ^b, where Q is flow discharge [cfs], w is channel width [ft]
a_w = cfg['HYDRAULIC_OPTIONS']['a_w']
b_w = cfg['HYDRAULIC_OPTIONS']['b_w']

#====================================================#
# Load latlon list for flow and energy grid cells
#====================================================#
print 'Loading routing station file...'
f = open(route_station_file, 'r')
#=== Read first line (# flow cells; # energy cells) ===#
line = f.readline().rstrip("\n")
n_flow = int(line.split()[0])
n_energy = int(line.split()[1])
#=== Loop over each flow/energy grid cell ==#
list_flow_lat_lon = []  # list of flow cells ([lat, lon])
list_flow_cell_number = []  # list of cell number of flow cells (corresponding to list_flow_lat_lon)
list_energy_lat_lon = []  # list of energy cells ([lat, lon])
cell_number = 1
while 1:
    line = f.readline().rstrip("\n")
    if line=="":
        break
    if line.split()[0]=="1":  # if this is a flow cell (i.e., not the end of a reach)
        lat_lon = line.split()[2]
        list_flow_lat_lon.append(lat_lon)
        list_flow_cell_number.append(cell_number)
        list_energy_lat_lon.append(lat_lon)
        line = f.readline().rstrip("\n")  # Read the next line - None or path for uh_s file
    elif line.split()[0]=="0":  # if this is NOT a flow cell (i.e., the end of a reach)
                                # Only energy grid cell
        lat_lon = line.split()[2]
        list_energy_lat_lon.append(lat_lon)
    cell_number = cell_number + 1
#=== Check whether the number of flow/energy cells are correct ===#
if len(list_flow_lat_lon)!=n_flow or len(list_energy_lat_lon)!=n_energy:
    print 'Error: incorrect number of flow/energy cells!'
    exit()

#====================================================#
# Load and process VIC output data - flow
#====================================================#
print 'Loading and processing VIC output flow data...'

#=== Load data ===#
ds_flow = xray.open_dataset(RVIC_output_nc)
da_flow = ds_flow['streamflow'][:-1,:,:]  # Remove the last junk time step

#=== Select time range ===#
da_flow = da_flow.sel(time=slice(start_date_str, end_date_str))

#=== Convert units ===#
da_flow = da_flow * pow(1000.0/25.4/12, 3)  # convert m3/s to cfs

#=== Set zero flow to 5.0 cfs ====#
da_flow.values[da_flow.values<5.0] = 5.0

#=== Calculate flow depth, width and velocity ===#
da_depth = a_d * pow(da_flow, b_d)  # flow depth [ft]
da_width = a_w * pow(da_flow, b_w)  # flow width [ft]
da_velocity = da_flow / da_depth / da_width  # flow velocoty [ft/s]

#====================================================#
# Rearrange data - flow
#====================================================#
#=== Put data for each grid cell into a df ===#
list_df = []
for i, lat_lon in enumerate(list_flow_lat_lon):
    print i+1
    lat = float(lat_lon.split('_')[0])
    lon = float(lat_lon.split('_')[1])
    df = pd.DataFrame(index=da_flow.coords['time'].values) # create df
    df['day'] = range(1, len(df)+1)  # day number (1,2,3,...)
    df['cell'] = list_flow_cell_number[i]  # cell number
    df['Q_in'] = da_flow.loc[:,lat,lon].values  # inflow discharge [cfs]
    df['Q_out'] = df['Q_in']  # outflow discharge [cfs]
    df['Q_diff'] = 0.0  # lateral flow [cfs]
    df['depth'] = da_depth.loc[:,lat,lon].values  # flow depth [ft]
    df['width'] = da_width.loc[:,lat,lon].values  # flow width [ft]
    df['velocity'] = da_velocity.loc[:,lat,lon].values  # flow velocity [ft/s]
    list_df.append(df)
    
#=== Combine df of all grid cells together, in a single multiindex df (indice: cell; date) ===#
df_flow = pd.concat(list_df, keys=list_flow_cell_number)

#=== Switch order of indice (to: date; cell), then sort ===#
df_flow = df_flow.reorder_levels([1,0], axis=0)
df_flow = df_flow.sortlevel(0)

#====================================================#
# Load and process VIC output data - energy
#====================================================#
print 'Loading and processing VIC output energy data...'

vic_energy_daily = {}  # keys: energy veriables; content: xray.DataArray of daily data

#=== Load data ===#
ds_vic_energy = xray.open_dataset(vic_output_energy_nc)

#=== If sub-daily data, convert to daily ===#
print 'Converting sub-daily to daily data...'
for var in ['Tair', 'vp', 'Shortwave', 'Longwave', 'Density', 'Pressure', 'Wind']:
    vic_energy_daily[var] = ds_vic_energy[cfg['INPUT'][var]]\
                                .groupby('time.date').mean(dim='time')

#=== Select time range for RBM ===#
print 'Selecting time range...'
for var in ['Tair', 'vp', 'Shortwave', 'Longwave', 'Density', 'Pressure', 'Wind']:
    vic_energy_daily[var] = vic_energy_daily[var].loc[start_date:end_date, :, :]

#=== Converting units ===#
print 'Converting units...'
vic_energy_daily['vp'] = vic_energy_daily['vp'] * 10.0 # convert [kPa] to [mb]
vic_energy_daily['Shortwave'] = vic_energy_daily['Shortwave'] \
                                * 2.388 * pow(10, -4) # convert [W/m2] to [mm*K/s]
vic_energy_daily['Longwave'] = vic_energy_daily['Longwave'] \
                               * 2.388 * pow(10, -4) # convert [W/m2] to [mm*K/s]
vic_energy_daily['Pressure'] = vic_energy_daily['Pressure'] * 10.0 # convert [kPa] to [mb]

#=== Close dataset ===#
ds_vic_energy.close()

#====================================================#
# Rearrange data - energy
#====================================================#
#=== Put data for each grid cell into a df ===#
list_df = []
for i, lat_lon in enumerate(list_energy_lat_lon):
    print i+1
    lat = float(lat_lon.split('_')[0])
    lon = float(lat_lon.split('_')[1])
    df = pd.DataFrame(index=vic_energy_daily['Tair'].coords['date'].values) # create df
    df['cell'] = i+1  # cell number
    df['Tair'] = vic_energy_daily['Tair'].loc[:,lat,lon].values # Tair
    df['vp'] = vic_energy_daily['vp'].loc[:,lat,lon].values # vp
    df['Shortwave'] = vic_energy_daily['Shortwave'].loc[:,lat,lon].values # Shortwave
    df['Longwave'] = vic_energy_daily['Longwave'].loc[:,lat,lon].values # Longwave
    df['Density'] = vic_energy_daily['Density'].loc[:,lat,lon].values # Density
    df['Pressure'] = vic_energy_daily['Pressure'].loc[:,lat,lon].values # Pressure
    df['Wind'] = vic_energy_daily['Wind'].loc[:,lat,lon].values # Wind
    list_df.append(df)

#=== Combine df of all grid cells together, in a single multiindex df (indice: cell; date) ===#
list_cell_index = []
for i in range(len(list_energy_lat_lon)):
    list_cell_index.append('cell{}'.format(i+1))
df_energy = pd.concat(list_df, keys=list_cell_index)

#=== Switch order of indice (to: date; cell), then sort ===#
df_energy = df_energy.reorder_levels([1,0], axis=0)
df_energy = df_energy.sortlevel(0)

#====================================================#
# Write data to file
#====================================================#
print 'Writing data to files...'
np.savetxt(cfg['OUTPUT']['rbm_flow_file'], df_flow.values, fmt='%d %d %.1f %.1f %.1f %.1f %.1f %.2f')
np.savetxt(cfg['OUTPUT']['rbm_energy_file'], df_energy.values, fmt='%d %.1f %.1f %.4f %.4f %.3f %.1f %.1f')
