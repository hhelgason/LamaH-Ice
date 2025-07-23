import os
import urllib.request
from datetime import datetime
import shutil
from pathlib import Path
import sys
import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
import xesmf as xe
import pickle
import hashlib
from multiprocessing import Pool

def print_progress(message):
    """Print a message with timestamp."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def check_free_space_gb(path):
    """Check available disk space in GB."""
    total, used, free = shutil.disk_usage(path)
    return free / (1024 ** 3)  # Convert bytes to GB

def generate_monthly_dates(start_date=199101, end_date=202503):
    """Generate list of dates from start_date to end_date (format: YYYYMM)."""
    dates = []
    current = start_date
    while current <= end_date:
        dates.append(current)
        year = current // 100
        month = current % 100
        if month == 12:
            current = (year + 1) * 100 + 1
        else:
            current = year * 100 + month + 1
    return dates

def download_file(url, filepath):
    """Download a file and return True if successful."""
    try:
        urllib.request.urlretrieve(url, filepath)
        return True
    except Exception as e:
        print_progress(f"❌ Download failed for {url}: {e}")
        return False

def get_grid_hash(ds, gdf):
    """Generate a hash based on grid and watershed properties."""
    grid_info = {
        'grid_shape': (ds.dims['y'], ds.dims['x']),
        'lon_bounds': (float(ds.lon.min()), float(ds.lon.max())),
        'lat_bounds': (float(ds.lat.min()), float(ds.lat.max())),
        'n_watersheds': len(gdf),
        'watershed_ids': sorted(gdf.index.tolist())
    }
    grid_str = str(grid_info)
    return hashlib.md5(grid_str.encode()).hexdigest()

def load_or_create_averager(ds, gdf, cache_dir):
    """Load spatial averager from cache if it exists, otherwise create and cache it."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True)
    
    grid_hash = get_grid_hash(ds, gdf)
    cache_file = cache_dir / f"spatial_averager_{grid_hash}.pkl"
    
    if cache_file.exists():
        print_progress(f"Loading cached spatial averager...")
        with open(cache_file, 'rb') as f:
            return pickle.load(f)
    
    print_progress("Creating xESMF spatial averager (this may take a while)...")
    savg = xe.SpatialAverager(
        ds.isel(Time=0),  # Use first timestep for grid setup
        gdf.geometry,
        geom_dim_name="id"
    )
    
    print_progress(f"Caching spatial averager...")
    with open(cache_file, 'wb') as f:
        pickle.dump(savg, f)
    
    return savg

def setup_coordinates(ds):
    """Prepare dataset coordinates for xESMF."""
    ds = ds.rename({'south_north': 'y', 'west_east': 'x'})
    
    lon = ds.XLONG.data
    lat = ds.XLAT.data
    
    # Calculate cell boundaries
    lon_b = np.zeros((lon.shape[0] + 1, lon.shape[1] + 1))
    lat_b = np.zeros((lat.shape[0] + 1, lat.shape[1] + 1))
    
    # Interior points
    lon_b[1:-1, 1:-1] = 0.25 * (lon[:-1, :-1] + lon[1:, :-1] + lon[:-1, 1:] + lon[1:, 1:])
    lat_b[1:-1, 1:-1] = 0.25 * (lat[:-1, :-1] + lat[1:, :-1] + lat[:-1, 1:] + lat[1:, 1:])
    
    # Edges and corners
    for arr in [lon_b, lat_b]:
        arr[0, 1:-1] = arr[1, 1:-1] + (arr[1, 1:-1] - arr[2, 1:-1])
        arr[-1, 1:-1] = arr[-2, 1:-1] + (arr[-2, 1:-1] - arr[-3, 1:-1])
        arr[1:-1, 0] = arr[1:-1, 1] + (arr[1:-1, 1] - arr[1:-1, 2])
        arr[1:-1, -1] = arr[1:-1, -2] + (arr[1:-1, -2] - arr[1:-1, -3])
        arr[0, 0] = arr[1, 1] + (arr[1, 1] - arr[2, 2])
        arr[0, -1] = arr[1, -2] + (arr[1, -2] - arr[2, -3])
        arr[-1, 0] = arr[-2, 1] + (arr[-2, 1] - arr[-3, 2])
        arr[-1, -1] = arr[-2, -2] + (arr[-2, -2] - arr[-3, -3])
    
    return ds.assign_coords({
        "lon": (("y", "x"), lon),
        "lat": (("y", "x"), lat),
        "lon_b": (("y_b", "x_b"), lon_b),
        "lat_b": (("y_b", "x_b"), lat_b)
    })

# List of variables to process
VARIABLES = [
    'wind_speed_at_10m_agl',
    'air_temperature_at_2m_agl',
    'relative_humidity_at_2m_agl',
    'lwe_precipitation_rate',
    'lwe_snow_depth',
    'net_downward_shortwave_flux',
    'surface_net_thermal_radiation',
    'lwe_solid_precipitation_rate',
    'surface_downward_sensible_heat_flux',
    'surface_downward_latent_heat_flux',
    'downward_shortwave_flux',
    'downward_longwave_flux'
]

def create_initial_averager(carra_dir, shapefile, cache_dir):
    """Create the spatial averager once at startup using the first CARRA file."""
    print_progress("Setting up spatial averager...")
    
    # Find first available CARRA file
    for date in generate_monthly_dates():
        carra_file = carra_dir / f"{date:06d}-carra-sfc_wod.nc"
        if carra_file.exists():
            print_progress(f"Using {carra_file.name} to create spatial averager")
            
            # Read data and prepare coordinates
            ds = xr.open_dataset(carra_file)
            ds = setup_coordinates(ds)
            
            # Read and prepare watersheds
            gdf = gpd.read_file(shapefile)
            gdf = gdf.to_crs("EPSG:4326")
            gdf["geometry"] = gdf.simplify(tolerance=0.001, preserve_topology=True)
            
            # Create a mapping between original IDs and sequential indices
            original_ids = sorted(gdf['id'].unique())
            id_to_index = {id_: idx for idx, id_ in enumerate(original_ids)}
            index_to_id = {idx: id_ for id_, idx in id_to_index.items()}
            
            # Create a new sequential index column and set it as the index
            gdf['sequential_index'] = gdf['id'].map(id_to_index)
            gdf = gdf.set_index('sequential_index')
            
            # Create averager
            savg = load_or_create_averager(ds, gdf, cache_dir)
            
            ds.close()
            return savg, gdf, id_to_index, index_to_id
            
    raise FileNotFoundError("No CARRA files found to create spatial averager")

def process_month(args):
    """Process one month of CARRA data."""
    date, carra_dir, savg, gdf, id_to_index, index_to_id, output_dir = args
    
    try:
        # Setup paths
        filename = f"{date:06d}-carra-sfc_wod.nc"
        carra_file = carra_dir / filename
        month_dir = output_dir / f"{date:06d}"
        month_dir.mkdir(exist_ok=True, parents=True)
        (month_dir / "csv").mkdir(exist_ok=True)
        
        if not carra_file.exists():
            print_progress(f"❌ CARRA file not found: {carra_file}")
            return date, False
        
        # Read data
        print_progress(f"Reading data: {filename}")
        ds = xr.open_dataset(carra_file)
        ds = setup_coordinates(ds)
        
        # Process variables sequentially
        results = {}
        for var in VARIABLES:
            if var not in ds:
                continue
                
            if 'Time' in ds[var].dims:
                result = savg(ds[var])
                times = pd.to_datetime([t.decode('utf-8').replace('_', ' ') for t in ds.Times.values])
                # Map sequential indices back to original IDs
                result = result.assign_coords(id=('id', [index_to_id[i] for i in result.id.values]))
                results[var] = result.assign_coords(Time=times)
            else:
                result = savg(ds[var])
                # Map sequential indices back to original IDs
                result = result.assign_coords(id=('id', [index_to_id[i] for i in result.id.values]))
                results[var] = result
        
        # Create output dataset
        ds_wsheds = xr.Dataset(results)
        
        # Save as netCDF
        ds_wsheds.to_netcdf(month_dir / "carra_watershed_averages.nc")
        
        # Save individual CSV files
        original_ids = sorted(index_to_id.values())
        for wshed_id in original_ids:
            wshed_data = {var: ds_wsheds[var].sel(id=wshed_id).to_pandas() 
                         for var in VARIABLES if var in ds_wsheds}
            df = pd.DataFrame(wshed_data)
            df = df.round(3)
            df.to_csv(month_dir / "csv" / f"watershed_{wshed_id}.csv")
        
        # Cleanup
        ds.close()
        return date, True
        
    except Exception as e:
        print_progress(f"❌ Error processing {date:06d}: {str(e)}")
        return date, False

def main():
    # Configuration - Windows paths
    carra_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\data")
    shapefile = Path(r"C:\Users\hordurbhe\OneDrive - Landsvirkjun\Documents\Vinna\lamah\lamah_ice\lamah_ice\A_basins_total_upstrm\3_shapefiles\Basins_A.shp")
    output_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds")
    cache_dir = output_dir / "cache"  # Central cache location
    
    # Set the number of cores to use.  
    physical_cores = 14  
    num_workers = max(1, physical_cores - 2)  
    
    print_progress(f"System has {physical_cores} physical cores")
    print_progress(f"Using {num_workers} workers for processing")
    
    # Create directories
    for d in [output_dir, cache_dir]:
        d.mkdir(exist_ok=True, parents=True)
    
    # Verify CARRA directory exists
    if not carra_dir.exists():
        print_progress(f"❌ CARRA directory not found: {carra_dir}")
        sys.exit(1)
    
    # Create spatial averager once at startup
    try:
        savg, gdf, id_to_index, index_to_id = create_initial_averager(carra_dir, shapefile, cache_dir)
        print_progress("✔ Spatial averager created successfully")
    except Exception as e:
        print_progress(f"❌ Failed to create spatial averager: {str(e)}")
        sys.exit(1)
    
    # Generate list of dates to process
    dates = generate_monthly_dates()
    
    # Create arguments for each process
    process_args = [(date, carra_dir, savg, gdf, id_to_index, index_to_id, output_dir) for date in dates]
    
    # Process months in parallel
    with Pool(num_workers) as pool:
        results = pool.map(process_month, process_args)
    
    # Report results
    successful = [date for date, success in results if success]
    failed = [date for date, success in results if not success]
    
    print_progress("\nProcessing complete!")
    print_progress(f"Successfully processed {len(successful)} months")
    if failed:
        print_progress(f"Failed to process {len(failed)} months: {failed}")
    
if __name__ == "__main__":
    main() 