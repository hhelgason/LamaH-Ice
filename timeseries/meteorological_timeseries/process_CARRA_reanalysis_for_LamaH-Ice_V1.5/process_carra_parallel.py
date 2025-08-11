import os
import urllib.request
from datetime import datetime
import shutil
from pathlib import Path
import sys
import argparse
from multiprocessing import Pool

import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
import xesmf as xe
import pickle  
import hashlib
import warnings

# -----------------------------
# Utilities & logging
# -----------------------------
def print_progress(message):
    """Print a message with timestamp (flush immediately)."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def check_free_space_gb(path):
    """Check available disk space in GB."""
    total, used, free = shutil.disk_usage(path)
    return free / (1024 ** 3)  # Convert bytes to GB

def ensure_proj_env():
    """
    Ensure pyproj/GDAL can locate PROJ data in this conda env.
    Safe to call multiple times.
    """
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        proj_dir = Path(conda_prefix) / "share" / "proj"
        if proj_dir.exists():
            os.environ.setdefault("PROJ_DATA", str(proj_dir))
            os.environ.setdefault("PROJ_LIB", str(proj_dir))

# -----------------------------
# Date helpers
# -----------------------------
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

# -----------------------------
# Download helper (unused here but kept for completeness)
# -----------------------------
def download_file(url, filepath):
    """Download a file and return True if successful."""
    try:
        urllib.request.urlretrieve(url, filepath)
        return True
    except Exception as e:
        print_progress(f"❌ Download failed for {url}: {e}")
        return False

# -----------------------------
# Grid hashing & coordinates
# -----------------------------
def get_grid_hash(ds, gdf):
    """
    Generate a hash based on grid and watershed properties to key the weight file.
    """
    # Robust min/max extraction (works if lon/lat are DataArrays or numpy arrays)
    lon_values = np.asarray(ds["lon"].values)
    lat_values = np.asarray(ds["lat"].values)

    grid_info = {
        'grid_shape': (int(ds.dims['y']), int(ds.dims['x'])),
        'lon_bounds': (float(np.nanmin(lon_values)), float(np.nanmax(lon_values))),
        'lat_bounds': (float(np.nanmin(lat_values)), float(np.nanmax(lat_values))),
        'n_watersheds': int(len(gdf)),
        # Use list of IDs; ensure deterministic order
        'watershed_ids': tuple(sorted(map(int, gdf.index.tolist()))) if gdf.index.dtype.kind in "iu" else tuple(sorted(gdf.index.tolist()))
    }
    grid_str = str(grid_info)
    return hashlib.md5(grid_str.encode()).hexdigest()

def setup_coordinates(ds):
    """
    Prepare dataset coordinates for xESMF SpatialAverager:
    - rename dims to ('y','x')
    - add lon/lat & lon_b/lat_b cell boundaries
    Assumes WRF-like fields: XLONG, XLAT, and dims south_north/west_east.
    """
    # rename dims if present
    if 'south_north' in ds.dims and 'west_east' in ds.dims:
        ds = ds.rename({'south_north': 'y', 'west_east': 'x'})
    elif 'y' not in ds.dims or 'x' not in ds.dims:
        raise ValueError("Dataset must have dims ('y','x') or ('south_north','west_east').")

    if 'XLONG' not in ds or 'XLAT' not in ds:
        raise KeyError("Dataset must contain 'XLONG' and 'XLAT' variables for grid coordinates.")

    lon = np.asarray(ds['XLONG'])
    lat = np.asarray(ds['XLAT'])

    if lon.shape != lat.shape:
        raise ValueError("XLONG and XLAT must have the same shape.")

    # Build cell boundary arrays (simple finite-difference extrapolation)
    ny, nx = lon.shape
    lon_b = np.zeros((ny + 1, nx + 1), dtype=float)
    lat_b = np.zeros((ny + 1, nx + 1), dtype=float)

    # Interior averages of surrounding 4 points
    lon_b[1:-1, 1:-1] = 0.25 * (lon[:-1, :-1] + lon[1:, :-1] + lon[:-1, 1:] + lon[1:, 1:])
    lat_b[1:-1, 1:-1] = 0.25 * (lat[:-1, :-1] + lat[1:, :-1] + lat[:-1, 1:] + lat[1:, 1:])

    # Edges (linear extrapolation)
    # top/bottom rows
    lon_b[0, 1:-1]  = lon_b[1, 1:-1]  + (lon_b[1, 1:-1]  - lon_b[2, 1:-1])
    lon_b[-1, 1:-1] = lon_b[-2, 1:-1] + (lon_b[-2, 1:-1] - lon_b[-3, 1:-1])
    lat_b[0, 1:-1]  = lat_b[1, 1:-1]  + (lat_b[1, 1:-1]  - lat_b[2, 1:-1])
    lat_b[-1, 1:-1] = lat_b[-2, 1:-1] + (lat_b[-2, 1:-1] - lat_b[-3, 1:-1])

    # left/right cols
    lon_b[1:-1, 0]  = lon_b[1:-1, 1]  + (lon_b[1:-1, 1]  - lon_b[1:-1, 2])
    lon_b[1:-1, -1] = lon_b[1:-1, -2] + (lon_b[1:-1, -2] - lon_b[1:-1, -3])
    lat_b[1:-1, 0]  = lat_b[1:-1, 1]  + (lat_b[1:-1, 1]  - lat_b[1:-1, 2])
    lat_b[1:-1, -1] = lat_b[1:-1, -2] + (lat_b[1:-1, -2] - lat_b[1:-1, -3])

    # corners
    lon_b[0, 0]   = lon_b[1, 1]   + (lon_b[1, 1]   - lon_b[2, 2])
    lon_b[0, -1]  = lon_b[1, -2]  + (lon_b[1, -2]  - lon_b[2, -3])
    lon_b[-1, 0]  = lon_b[-2, 1]  + (lon_b[-2, 1]  - lon_b[-3, 2])
    lon_b[-1, -1] = lon_b[-2, -2] + (lon_b[-2, -2] - lon_b[-3, -3])

    # attach as coords
    ds = ds.assign_coords({
        "lon": (("y", "x"), lon),
        "lat": (("y", "x"), lat),
        "lon_b": (("y_b", "x_b"), lon_b),
        "lat_b": (("y_b", "x_b"), lat_b)
    })
    return ds

# -----------------------------
# xESMF weights caching
# -----------------------------
def load_or_create_averager_weights(ds, gdf, cache_dir):
    """
    Ensure a weights file exists for this grid+watershed combo and return its path.
    We cache the expensive step (building weights) to a NetCDF file.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)

    grid_hash = get_grid_hash(ds, gdf)
    weights_path = cache_dir / f"spatial_averager_weights_{grid_hash}.nc"

    if weights_path.exists():
        print_progress(f"Using cached spatial-averager weights: {weights_path.name}")
        return weights_path

    print_progress("Creating xESMF spatial averager weights (this may take a while)...")
    # Build once; with filename=..., xESMF writes weights to disk.
    _ = xe.SpatialAverager(
        ds.isel(Time=0) if 'Time' in ds.dims else ds,
        gdf.geometry,
        geom_dim_name="id",
        filename=str(weights_path)
    )
    print_progress(f"Saved weights to {weights_path}")
    return weights_path

def build_averager_from_weights(ds, gdf, weights_path):
    """
    Reconstruct a SpatialAverager cheaply from on-disk weights in the worker process.
    """
    savg = xe.SpatialAverager(
        ds.isel(Time=0) if 'Time' in ds.dims else ds,
        gdf.geometry,
        geom_dim_name="id",
        filename=str(weights_path),
        reuse_weights=True
    )
    return savg

# -----------------------------
# Variables to extract
# -----------------------------
VARIABLES = [
    'lwe_precipitation_rate',                           # prec_carra
    'lwe_solid_precipitation_rate',                     # solid_prec_carra
    'air_temperature_at_2m_agl',                        # 2m_temp_carra
    'wind_speed_at_10m_agl',                            # 10m_wind_speed_carra
    'wind_from_direction_at_10m_agl',                   # 10m_wind_dir_carra
    'relative_humidity_at_2m_agl',                      # 2m_rel_hum_carra
    'specific_humidity_at_2m_agl',                      # 2m_spec_hum_carra
    'net_downward_shortwave_flux',                      # surf_net_solar_rad_carra
    'surface_net_thermal_radiation',                    # surf_net_therm_rad
    'downward_shortwave_flux',                          # surf_dwn_solar_rad_carra
    'downward_longwave_flux',                           # surf_dwn_therm_rad_carra
    'lwe_snow_depth',                                   # swe_carra
    'surface_upward_latent_heat_flux_due_to_sublimation',  # snow_sublimation_carra
    'surface_downward_sensible_heat_flux',              # surf_dwn_sens_heat_flux_carra
    'surface_downward_latent_heat_flux',                # surf_dwn_lat_heat_flux_carra
    'water_evaporation_amount',                         # total_et_carra
    'lwe_percolation_rate',                             # percolation_carra
    'lwe_runoff_surface_rate'                           # runoff_carra
]

# -----------------------------
# Setup & precompute
# -----------------------------
def create_initial_weights_and_metadata(carra_dir, shapefile, cache_dir):
    """
    Use the first available CARRA file to:
    - build lon/lat coordinates
    - read/simplify watersheds
    - compute and cache averager weights (once)
    Returns: (weights_path, gdf, id_to_index, index_to_id)
    """
    print_progress("Setting up spatial averager...")

    # Find first available CARRA file
    for date in generate_monthly_dates():
        carra_file = carra_dir / f"{date:06d}-carra-sfc_wod.nc"
        if carra_file.exists():
            print_progress(f"Using {carra_file.name} to create spatial averager")

            ds = xr.open_dataset(carra_file)
            ds = setup_coordinates(ds)

            # Read and prepare watersheds
            ensure_proj_env()
            gdf = gpd.read_file(shapefile)

            if gdf.crs is None:
                # If your shapefile lacks CRS, set it to the known source CRS before reprojecting
                warnings.warn("Shapefile has no CRS. Assuming EPSG:3057 (Iceland TM) — adjust if incorrect.")
                gdf.set_crs(epsg=3057, inplace=True)

            gdf = gdf.to_crs("EPSG:4326")
            gdf["geometry"] = gdf.simplify(tolerance=0.001, preserve_topology=True)

            # Create a mapping between original IDs and sequential indices
            if 'id' not in gdf.columns:
                raise KeyError("Shapefile must contain an 'id' column with unique watershed IDs.")

            original_ids = sorted(gdf['id'].unique())
            id_to_index = {id_: idx for idx, id_ in enumerate(original_ids)}
            index_to_id = {idx: id_ for id_, idx in id_to_index.items()}

            # Create a new sequential index column and set it as the index
            gdf['sequential_index'] = gdf['id'].map(id_to_index)
            gdf = gdf.set_index('sequential_index')

            # Create/locate weights file (expensive step done once)
            weights_path = load_or_create_averager_weights(ds, gdf, cache_dir)

            ds.close()
            return weights_path, gdf, id_to_index, index_to_id

    raise FileNotFoundError("No CARRA files found in carra_dir to create spatial averager")

# -----------------------------
# Worker
# -----------------------------
def _decode_times(ds):
    """
    Try to decode times from either WRF-like byte strings in ds.Times
    or a 'Time' coordinate; fall back safely.
    """
    if 'Times' in ds.variables:
        try:
            # WRF-style bytes -> strings like 'YYYY-MM-DD_HH:MM:SS'
            times_raw = ds['Times'].values
            if hasattr(times_raw, 'dtype') and times_raw.dtype.kind in ('S', 'O'):
                strings = [t.decode('utf-8').replace('_', ' ') if isinstance(t, (bytes, bytearray)) else str(t) for t in times_raw]
                return pd.to_datetime(strings)
        except Exception:
            pass

    if 'Time' in ds.coords:
        try:
            return pd.to_datetime(ds['Time'].values)
        except Exception:
            pass

    # Last resort: make a RangeIndex of length of the first time-like dim
    if 'Time' in ds.dims:
        return pd.date_range("1970-01-01", periods=ds.dims['Time'], freq='H')

    return None

def process_month(args):
    """
    Process one month of CARRA data.
    args: (date, carra_dir, weights_path, gdf, id_to_index, index_to_id, output_dir)
    """
    (date, carra_dir, weights_path, gdf, id_to_index, index_to_id, output_dir) = args

    try:
        ensure_proj_env()

        # Paths
        filename = f"{date:06d}-carra-sfc_wod.nc"
        carra_file = carra_dir / filename
        month_dir = output_dir / f"{date:06d}"
        (month_dir / "csv").mkdir(exist_ok=True, parents=True)

        if not carra_file.exists():
            print_progress(f"❌ CARRA file not found: {carra_file}")
            return date, False

        # Read data & ensure coords exist
        print_progress(f"Reading data: {filename}")
        ds = xr.open_dataset(carra_file)
        ds = setup_coordinates(ds)

        # Rebuild SpatialAverager locally (fast) from cached weights
        savg = build_averager_from_weights(ds, gdf, weights_path)

        # Decode times once (if present)
        times = _decode_times(ds)

        # Process variables sequentially
        results = {}
        for var in VARIABLES:
            if var not in ds:
                continue

            da = ds[var]
            # Apply spatial averaging
            out = savg(da)

            # Map sequential indices back to original IDs
            out = out.assign_coords(id=('id', [index_to_id[i] for i in out['id'].values]))

            # Attach times if present
            if 'Time' in out.dims and times is not None:
                out = out.assign_coords(Time=times)

            results[var] = out

        if not results:
            print_progress(f"❌ No expected variables found for {filename}")
            ds.close()
            return date, False

        # Create output dataset
        ds_wsheds = xr.Dataset(results)

        # Save as netCDF
        out_nc = month_dir / "carra_watershed_averages.nc"
        ds_wsheds.to_netcdf(out_nc)

        # Save individual CSV files
        original_ids_sorted = sorted(index_to_id.values())
        for wshed_id in original_ids_sorted:
            series_dict = {}
            for var in VARIABLES:
                if var in ds_wsheds:
                    try:
                        series = ds_wsheds[var].sel(id=wshed_id).to_pandas()
                        series_dict[var] = series
                    except Exception:
                        pass
            if series_dict:
                df = pd.DataFrame(series_dict).round(3)
                df.to_csv(month_dir / "csv" / f"watershed_{wshed_id}.csv")

        ds.close()
        return date, True

    except Exception as e:
        print_progress(f"❌ Error processing {date:06d}: {str(e)}")
        return date, False

# -----------------------------
# Main
# -----------------------------
def main():

    carra_dir = Path(r"/home/LV/xhordurbhe/verkefni/CARRA")
    shapefile = Path(r"/home/LV/xhordurbhe/verkefni/CARRA/data/shapefiles/Basins_A.shp")
    output_dir = Path(r"/home/LV/xhordurbhe/verkefni/CARRA/data/output")
    cache_dir = output_dir / "cache"  # Central cache location

    # Choose workers
    num_workers = 1 

    # Create directories
    for d in [output_dir, cache_dir]:
        d.mkdir(exist_ok=True, parents=True)

    # Verify CARRA directory exists
    if not carra_dir.exists():
        print_progress(f"❌ CARRA directory not found: {carra_dir}")
        sys.exit(1)

    # Precompute weights and metadata (does NOT pickle ESMF objects)
    try:
        weights_path, gdf, id_to_index, index_to_id = create_initial_weights_and_metadata(
            carra_dir, shapefile, cache_dir
        )
        print_progress("✔ Spatial-averager weights ready")
    except Exception as e:
        print_progress(f"❌ Failed to create/load spatial-averager weights: {str(e)}")
        sys.exit(1)

    # Generate list of dates to process - using default values from generate_monthly_dates()
    dates = generate_monthly_dates()  # This will use default start_date=199101, end_date=202503

    # Prepare args for each month
    process_args = [
        (date, carra_dir, weights_path, gdf, id_to_index, index_to_id, output_dir)
        for date in dates
    ]

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
    # Make sure PROJ vars are present in this process
    ensure_proj_env()
    main()
