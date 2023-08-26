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
pds.set_option('display.max_rows', 200)
import datetime
import matplotlib.pyplot as plt
import regionmask
import multiprocessing

# Define what variables are to be extracted 
variables = ['GLW','GRAUPELNC','GRDFLX','HFX','LH','OLR','Q2','RAINNC','SNOW','SR','SWDOWN','T2','TSK','U10','V10'] #WIND=SQRT(U102+V102) (http://www.meteo.unican.es/wiki/cordexwrf/OutputVariables)
var_names = ['DOWNWARD LONG WAVE FLUX AT GROUND SURFACE','ACCUMULATED TOTAL GRID SCALE GRAUPEL','GROUND HEAT FLUX','UPWARD HEAT FLUX AT THE SURFACE','LATENT HEAT FLUX AT THE SURFACE','TOA OUTGOING LONG WAVE','Specific humidity','ACCUMULATED TOTAL GRID SCALE PRECIPITATION','SNOW WATER EQUIVALENT','fraction of frozen precipitation','DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE','TEMP at 2 M','SURFACE SKIN TEMPERATURE','U at 10 M','V at 10 M']

# First we read a RAV file and create an xarray dataset to calculate the xesmf averager
year = 1980
path = '/data/rav/data/%s-%s' % (year,str(year+1)[-2:])
filenames = glob.glob(path + '/*')
ds_coords = xr.open_dataset(filenames[0])

times = []
for x in range(len(ds_coords['Times'].values)): 
    times.append(pds.to_datetime((ds_coords['Times'].values[x]).decode("utf-8").replace('_',' ')))

ds = xr.open_dataset(filenames[0])
ds = ds.assign_coords({"time":times ,"latitude": ds_coords['XLAT'].sel(Time=1),"longitude": ds_coords['XLONG'].sel(Time=1)})

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
new_x = ds.longitude.values - a # In CARRA we had to subtract 360
# Append last edge column
new_x = np.insert(new_x,-1,new_x[:,-1]+a[:,-1],axis=1) 
# Also insert the last edge column to a before we can create the last edge row of the new bounds
a = np.insert(a,-1,a[:,-1],axis=1) # Insert last edge column to a 
new_xx = np.insert(new_x,len(ds.south_north),new_x[-1]+2*a[-1],axis=0) # Append last edge row to new_x
# Eða: ds.longitude[-1] + a[-1]

b = dy.values/2
b = np.insert(b,0,b[:,0],axis=1)
new_y = ds.latitude.values - b

# Append last edge line to new_y
new_y = np.insert(new_y,-1,new_y[-1],axis=0) # Insert last edge line
# Insert last edge line to b
b = np.insert(b,-1,b[-1],axis=0) # Insert last edge line
# Insert last edge column to new_y
new_yy = np.insert(new_y,len(ds.west_east),new_y[:,-1]-2*b[:,-1],axis=1) # Insert last edge column

ds['lon_b'] = (['x_b','y_b'],new_xx)
ds['lat_b'] = (['x_b','y_b'],new_yy)

# Read the watershed shapefile
path = '/home/xhordurbhe/Verkefni/lamah_ice/data/Basins_A_with_glaciers_wgs84.shp'
gdf = gpds.read_file(path)

#Simplify the geometries to a 0.001 deg tolerance
gdf["geometry"] = gdf.simplify(tolerance=0.001, preserve_topology=True)
gdf = gdf.set_index('id')

# Calculate the spatial averager:
savg = xe.SpatialAverager(ds.sel(time=ds.time[0]), gdf.geometry, geom_dim_name="wshed")

cumulated_variables=['RAINNC','GRAUPELNC']

# Now we calculate the averages of the weather parameters on the watersheds
# Folder to save a concatenated netcdf file for each year:
folder = '/home/xhordurbhe/Verkefni/lamah_ice/rav2_on_watersheds_Basins_A'


def process_file(filename):
    ds_c = xr.open_dataset(filename)
    # Calculate spatial averages
    out = savg(ds_c)
    # Extract datetime values:
    times = []
    for x in range(len(ds_c['Times'].values)):
        times.append(pds.to_datetime((ds_c['Times'].values[x]).decode("utf-8").replace('_', ' ')))

    out = out.assign_coords(wshed=xr.DataArray(gdf.index, dims=("wshed",)), Time=xr.DataArray(times, dims=("Time",)))
    out = out.drop_vars(['longitude', 'latitude'])

    # Calculate daily means
    ds_wsheds = out[variables]  # .resample(time='1D').mean()
    # We now have an xarray dataset with hourly averages over all watersheds (dimension: "wshed")
    return ds_wsheds

def process_year_range(year):
    path = '/data/rav/data/%s-%s' % (year, str(year + 1)[-2:])
    filenames = glob.glob(path + '/*')

    listi = []
    for filename in filenames:
        print(filename)
        ds_wsheds = process_file(filename)
        listi.append(ds_wsheds)

    # Create a concatted xarray dataset with the full year
    print('finished all files for year %s. now concatenating:' % year)
    concatted = xr.concat(listi, dim='Time')

    # Sort by time and save the dataset
    concatted = concatted.sortby('Time')
    path = folder + '/%s.nc' % year
    concatted.to_netcdf(path)
    print('xesmf spatial averaging complete for %s' % year)

def main():
    # ... (Existing code for creating the spatial averager, loading the watershed shapefile, etc.)

    # Define the range of years to process
    start_year = 1958
    end_year = 2019
    years_to_process = range(start_year, end_year)

    # Define the number of processes (you can adjust this based on the number of CPU cores available)
    num_processes = 16 #multiprocessing.cpu_count()

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=num_processes)

    # Process the years in parallel using the pool
    pool.map(process_year_range, years_to_process)

    # Close the pool to release resources
    pool.close()
    pool.join()

    # Combine the netCDF datasets and export .csv files:
    filenames = glob.glob(folder + '/*.nc')
    ds_all_years = xr.open_mfdataset(filenames)
    ds_all_years = ds_all_years.sortby('Time')

    for wshed in ds_all_years.wshed.values:
        dff = ds_all_years.sel(wshed=wshed).to_dataframe()

        # Handle variables that WRF accumulates:
        cumulated_variables = ['RAINNC', 'GRAUPELNC']
        for cumvar in cumulated_variables:
            raindiff = dff[cumvar].diff()
            raindiff[raindiff < 0] = 0
            raindiff = raindiff.fillna(0)
            dff[cumvar] = raindiff
        dff.to_csv(folder + '/csv/ID_%s.csv' % wshed)

if __name__ == "__main__":
    main()

