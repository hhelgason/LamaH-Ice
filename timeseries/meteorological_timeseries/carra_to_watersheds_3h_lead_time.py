from pathlib import Path
import xagg as xa
import xarray as xr
import pandas as pds
import geopandas as gpds
import numpy as np
import shapefile as shp
import xesmf as xe
from shapely.geometry import Polygon, MultiPolygon
import glob
import datetime
import matplotlib.pyplot as plt

# We use the xesmf package to calculate spatial averages

# Folder to save a concatenated netcdf file for each year:
folder = '/home/xhordurbhe/Verkefni/carra/carra_watershed_averages/3h_lead_time'

# Read the watershed shapefile
path = '/home/xhordurbhe/Verkefni/lamah_ice/data/compiled_with_glaciers_wgs84.shp'
gdf = gpds.read_file(path)

#Simplify the geometries to a 0.001 deg tolerance
gdf["geometry"] = gdf.simplify(tolerance=0.001, preserve_topology=True)
gdf = gdf.set_index('id')

# First we read a carra file and calculate the xesmf averager
ds = xr.open_dataset('/data/carra/carra_island/CARRA_WestDomain_daily_Iceland_1990_10_17.nc')

# We need to add the grid-cell bounds to the dataset as dimensions
# So we calculate the bounds as follows:

# This is the difference between center-points of grid cells
dx = ds.longitude[1:,:] - ds.longitude[:-1,:]
dy = ds.latitude[:,1:] - ds.latitude[:,:-1]

# We divide by two and insert the missing edge line
a = dx.values/2
a = np.insert(a,0,a[0],axis=0) # Insert first edge line 
# a = np.insert(a,-1,a[-1],axis=0) # Insert last edge line 

# Now we calculate the bounds by shifting the original longitude values by the distance between grid centers divided by 2
new_x = ds.longitude.values - 360 - a
# Append last edge column
new_x = np.insert(new_x,-1,new_x[:,-1]+a[:,-1],axis=1) 
# Also insert the last edge column to a before we can create the last edge row of the new bounds
a = np.insert(a,-1,a[:,-1],axis=1) # Insert last edge column to a 
new_xx = np.insert(new_x,280,new_x[-1]+2*a[-1],axis=0) # Append last edge row to new_x
# Eða: ds.longitude[-1] + a[-1]

b = dy.values/2
b = np.insert(b,0,b[:,0],axis=1)
new_y = ds.latitude.values - b

# Append last edge line to new_y
new_y = np.insert(new_y,-1,new_y[-1],axis=0) # Insert last edge line
# Insert last edge line to b
b = np.insert(b,-1,b[-1],axis=0) # Insert last edge line
# Insert last edge column to new_y
new_yy = np.insert(new_y,330,new_y[:,-1]-2*b[:,-1],axis=1) # Insert last edge column

ds['lon_b'] = (['x_b','y_b'],new_xx)
ds['lat_b'] = (['x_b','y_b'],new_yy)

# Calculate the spatial averager:
savg = xe.SpatialAverager(ds.sel(time=ds.time[0]), gdf.geometry, geom_dim_name="wshed")

# Now we process the 3h lead time CARRA data (the precipitation is correct in this data)
# We calculate the averages of the weather parameters on the watersheds

ds_dict = dict()
for year in np.arange(1991,2010):
    listi = []
    filenames=glob.glob('/data/carra/carra_island_3h_lead/CARRA_WestDomain_daily_Iceland_'+str(year)+'*')
    for filename in filenames:
        ds_c = xr.open_dataset(filename)
        # Calculate spatial averages
        out = savg(ds_c)
        out = out.assign_coords(wshed=xr.DataArray(gdf.index, dims=("wshed",)))
        out=out.drop_vars(['heightAboveGround','longitude','latitude'])
        # Calculate daily means
        ds_wsheds = out[['wdir10', 'si10', 'rsn', 'sd', 'sp', 'sro', 'ssrd', 'strd']].resample(time='1D').mean()
        ds_wsheds['tp'] = out['tp'].resample(time='1D').sum()
        # We now have an xarray dataset with a single day´s average over all watersheds (dimension: "wshed")
        # We append this to the list of xarray datasets
        listi.append(ds_wsheds)
        
    # Create a concatted xarray dataset with the full year
    concatted = xr.concat(listi,dim='time') 
    # Sort by time and save the dataset
    concatted = concatted.sortby('time')
    path = folder + '/%s.nc' % year
    concatted.to_netcdf(path)
    print('xesmf spatial averaging complete for %s' % year)

# Combine the netCDF datasets and export .csv files:
filenames = glob.glob(folder+'/*.nc')
ds_all_years = xr.open_mfdataset(filenames) 
ds_all_years = ds_all_years.sortby('time')

for wshed in ds_all_years.wshed.values:
    dff = ds_all_years.sel(wshed=wshed).to_dataframe() 
    dff.to_csv(folder + '/csv/ID_%s.csv' % wshed)

# Create csv files that only contains precipitation:
for wshed in ds_all_years.wshed.values:
    df_h = pds.read_csv(folder + '/csv/ID_%s.csv' % wshed)
    df_h = df_h.set_index('time')
    df_h.index=pds.to_datetime(df_h.index)
    df_h = pds.DataFrame(df_h['tp'])
    df_h.to_csv(folder + '/csv/precip_daily/ID_%s.csv' % wshed)
