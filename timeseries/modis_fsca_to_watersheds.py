import xarray as xr
from pathlib import Path
import pyproj
import numpy as np
import pandas as pds
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import geopandas as gpds
import regionmask
import pickle
import os
import glob
import datetime
from multiprocessing import Pool

# Define a function that creates an xarray dataset containing modis data for one day

def modis_to_masked_xr(ds, geo, date, mask, variable_name):
   
    """" This function creates a masked xarray dataset containing modis data for one day
         Inputs are:
         ds: Raw Modis netcdf file
         geo: The modis geo file (values already rotated)
         date: The day of the modis file
         mask: The regionmask mask (watersheds)
         variable_name: The name of the variable we're extracting from the modis file. Ex.: 'R_median_filter_fsca_GFD' """

    # First we rotate the values in the MODIS file
    m = ds[variable_name].values
    rotated = np.rot90(m, 1)
    ds[variable_name].values = rotated    
    
    # Now we create an xarray dataset
    ds_new = xr.Dataset(
        data_vars=dict(R_median_filter_fsca_GFD=(["x", "y", "time"], ds[variable_name].values[..., np.newaxis])), #'R_median_fsca'
        coords=dict(
            lon=(["x", "y"], geo['longitude'].values),
            lat=(["x", "y"], geo['latitude'].values),
            time=pds.date_range(date, periods=1),
            reference_time=pds.Timestamp(date)),
        attrs=dict(description='MODIS fractional snow cover'),)
    
    # Now we mask the dataset
    ds_masked = ds_new['R_median_filter_fsca_GFD'].where(mask)
    return(ds_masked)

# Define a function to process a single file
def process_file(filename):
    # extract the date string from the filename (assuming the date always appears as the last 10 characters)
    date_str = filename[-13:-3]
    print(date_str)
    # convert the date string to a datetime object
    date_obj = datetime.datetime.strptime(date_str, '%d_%m_%Y')
    # format the datetime object as a string in the desired format
    date = date_obj.strftime('%Y-%m-%d')
    # open the modis dataset
    dataset = xr.open_dataset(filename)
    ds = dataset.copy()
    dataset.close()
    
    # Run the modis_to_masked_xr function for a single date
    masked_modis_fsca = modis_to_masked_xr(ds, geo, date, mask, variable_name)

    # Calculate the mean snow cover on watersheds
    mean = masked_modis_fsca.mean(dim=['x', 'y'])

    # convert to a pandas DataFrame
    df = mean.to_dataframe()

    # reset the index to get the region number as a column
    df = df.reset_index()

    # Set the region as index
    df = df.set_index('region')
    df = df.sort_index()

    # drop the unwanted columns
    df = df.drop(columns=['time']) #, 'reference_time'])
    try:
        df = df.drop(columns=['reference_time'])
    except:
        pass
    df.columns = [date]    
    df = df.transpose()
    
    save_path = os.path.join(save_folder,date+'.csv')
    df.round(3).to_csv(save_path)

# Read the reprojected watersheds dataframe:
path = '/home/xhordurbhe/Verkefni/lamah_ice/data/final_shapefiles/Basins_A_with_glaciers_wgs84.shp'
gdf = gpds.read_file(path)
gdf = gdf.set_index('id')
gdf = gdf.sort_index()

# Simplify the geometries to a 0.01 deg tolerance
gdf["geometry"] = gdf.simplify(tolerance=0.01, preserve_topology=True)

# Read a sample modis file
file_folder = '/data/MODreSWE/reSWE_MOD/'
file_to_open = file_folder + 'reSWE_mod_fsca_alb_01_04_2016.nc'
ds = xr.open_dataset(file_to_open)
variable_name = 'R_median_filter_fsca_GFD'
date = '2016-04-01'

# Read the geo file
file_folder = '/home/xhordurbhe/Verkefni/lamah_ice/data/'
file_to_open = file_folder + "MOD_h17v02_geo.nc"
geo = xr.open_dataset(file_to_open)

# Rotate the values
m = ds[variable_name].values
rotated = np.rot90(m, 1)
ds[variable_name].values = rotated

m = geo['latitude'].values
rotated = np.rot90(m, 1)
geo['latitude'].values = rotated

m = geo['longitude'].values
rotated = np.rot90(m, 1)
geo['longitude'].values = rotated

mask_path = '/home/xhordurbhe/Verkefni/lamah_ice/modis_on_watersheds/regionmask_mask_final.p'
if os.path.isfile(mask_path):
    # The mask already exists. We load the mask:
    mask = pickle.load( open( mask_path, "rb" ) )
else:
    # We create a mask for the modis data on watersheds, using regionmask.       
    # First we create an xarray dataset from the rotated data:
    ds_for_mask_gen = xr.Dataset(
        data_vars=dict(R_median_filter_fsca_GFD=(["x", "y", "time"], ds[variable_name].values[..., np.newaxis])), #'R_median_fsca'
        coords=dict(
            lon=(["x", "y"], geo['longitude'].values),
            lat=(["x", "y"], geo['latitude'].values),
            time=pds.date_range(date, periods=1),
            reference_time=pds.Timestamp(date)),
        attrs=dict(description='MODIS fractional snow cover'),)

    mask = regionmask.mask_3D_geopandas(gdf, ds_for_mask_gen.lon, ds_for_mask_gen.lat,overlap=True)

    # Save the mask (because it takes a while to generate)
    pickle.dump( mask, open( mask_path, "wb" ) )



if __name__ == '__main__':
    filenames=glob.glob('/data/MODreSWE/reSWE_MOD/reSWE_mod_fsca_alb_*')
    save_folder='/home/xhordurbhe/Verkefni/lamah_ice/modis_on_watersheds/data_fcsa_full_watershed_final/'
    now = datetime.datetime.now()
    print(now)
    pool = Pool(processes=8)  # Use 4 processors
    pool.map(process_file, filenames)  # Process all files
    pool.close()
    pool.join()
    now = datetime.datetime.now()
    print(now)

# Read the .csv files back in and combine them
filenames=glob.glob(save_folder+'/*')

read_back_in = dict()
for filename in filenames:
    date = filename[-14:-4]
    df = pds.read_csv(filename)
    df = df.set_index('Unnamed: 0')
    df.index = df.index.rename('date_time')
    read_back_in[date] = df
    
df_combined = pds.concat(read_back_in.values())
df_combined.index=pds.to_datetime(df_combined.index)
df_combined = df_combined.sort_index()

# Resample to daily frequency
daily_df_combined = df_combined.resample('D').mean()

# Check for missing values
missing_dates = daily_df_combined[daily_df_combined.isna().any(axis=1)].index
print(missing_dates)

save_path = '/home/xhordurbhe/Verkefni/lamah_ice/modis_on_watersheds/timeseries/modis_fsca_full_watershed_final.csv'
daily_df_combined.to_csv(save_path)
