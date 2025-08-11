import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

def print_progress(message):
    """Print a message with timestamp."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def load_carra_data(carra_dir, watershed_id):
    """Load CARRA daily timeseries for a watershed."""
    carra_file = Path(carra_dir) / f"watershed_{watershed_id}.csv"
    
    if not carra_file.exists():
        print_progress(f"❌ CARRA file not found: {carra_file}")
        return None
    
    try:
        df = pd.read_csv(carra_file, index_col=0, parse_dates=True)
        print_progress(f"  Loaded CARRA data: {df.index.min()} to {df.index.max()}")
        return df
    except Exception as e:
        print_progress(f"❌ Error reading CARRA file: {str(e)}")
        return None

def load_lamah_data(lamah_dir, watershed_id):
    """Load LAMAH-ICE daily timeseries for a watershed."""
    lamah_file = Path(lamah_dir) / f"ID_{watershed_id}.csv"
    
    if not lamah_file.exists():
        print_progress(f"❌ LAMAH file not found: {lamah_file}")
        return None
    
    try:
        # Read the semicolon-separated file, keeping empty strings as empty strings
        df = pd.read_csv(lamah_file, sep=';', keep_default_na=False)
        
        # Create date index for merging (but don't save it)
        df['_date'] = pd.to_datetime(df[['YYYY', 'MM', 'DD']].astype(str).agg('-'.join, axis=1))
        df.set_index('_date', inplace=True)
        
        print_progress(f"  Loaded LAMAH data: {df.index.min()} to {df.index.max()}")
        return df
        
    except Exception as e:
        print_progress(f"❌ Error reading LAMAH file: {str(e)}")
        return None

def merge_data(lamah_df, carra_df):
    """Merge CARRA data into LAMAH dataframe with appropriate column renaming."""
    
    # Create a copy of LAMAH data to avoid modifying the original
    merged_df = lamah_df.copy()
    
    # Define mapping of CARRA variable names to new column names
    carra_columns = {
        'lwe_precipitation_rate': 'prec_carra',
        'lwe_solid_precipitation_rate': 'solid_prec_carra',
        'air_temperature_at_2m_agl': '2m_temp_carra',
        'air_temperature_at_2m_agl_min': '2m_temp_min_carra',
        'air_temperature_at_2m_agl_max': '2m_temp_max_carra',
        'wind_speed_at_10m_agl': '10m_wind_speed_carra',
        'wind_from_direction_at_10m_agl': '10m_wind_dir_carra',
        'relative_humidity_at_2m_agl': '2m_rel_hum_carra',
        'specific_humidity_at_2m_agl': '2m_spec_hum_carra',
        'net_downward_shortwave_flux': 'surf_net_solar_rad_carra',
        'surface_net_thermal_radiation': 'surf_net_therm_rad_carra',
        'downward_shortwave_flux': 'surf_dwn_solar_rad_carra',
        'downward_longwave_flux': 'surf_dwn_therm_rad_carra',
        'lwe_snow_depth': 'swe_carra',
        'surface_upward_latent_heat_flux_due_to_sublimation': 'snow_sublimation_carra',
        'surface_downward_sensible_heat_flux': 'surf_dwn_sens_heat_flux_carra',
        'surface_downward_latent_heat_flux': 'surf_dwn_lat_heat_flux_carra',
        'water_evaporation_amount': 'total_et_carra',
        'lwe_percolation_rate': 'percolation_carra',
        'lwe_runoff_surface_rate': 'runoff_carra'
    }
    
    # Print info about the dataframes before merging
    print_progress("\nData ranges:")
    print_progress(f"LAMAH: {lamah_df.index.min()} to {lamah_df.index.max()}")
    print_progress(f"CARRA: {carra_df.index.min()} to {carra_df.index.max()}")
    
    # Limit CARRA data to LAMAH period
    carra_df = carra_df[carra_df.index <= lamah_df.index.max()]
    print_progress(f"CARRA (limited to LAMAH period): {carra_df.index.min()} to {carra_df.index.max()}")
    
    # Add CARRA columns to merged dataframe
    n_added = 0
    for carra_name, new_name in carra_columns.items():
        if carra_name in carra_df.columns:
            # Initialize the column with empty string
            merged_df[new_name] = ''
            
            # Fill in values only where CARRA data exists
            merged_df.loc[carra_df.index, new_name] = carra_df[carra_name].round(3)
            
            # Count non-empty values
            n_values = (merged_df[new_name] != '').sum()
            n_total = len(merged_df)
            print_progress(f"  Added {new_name}: {n_values} values ({n_values/n_total*100:.1f}% coverage)")
            n_added += 1
    
    print_progress(f"\nSummary:")
    print_progress(f"  Added {n_added} CARRA variables")
    print_progress(f"  CARRA data available from {carra_df.index.min()} to {carra_df.index.max()}")
    print_progress(f"  Full time range (with LAMAH): {merged_df.index.min()} to {merged_df.index.max()}")
    
    # Reset index to remove the date column
    merged_df.reset_index(drop=True, inplace=True)
    
    return merged_df

def main():
    # Paths
    carra_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds_all_vars\daily_series")
    lamah_dir = Path(r"C:\Users\hordurbhe\OneDrive - Landsvirkjun\Documents\Vinna\lamah\lamah_ice\lamah_ice\A_basins_total_upstrm\2_timeseries\daily\meteorological_data")
    output_dir = lamah_dir 
    
    # Get list of all watershed IDs
    watershed_ids = []
    for file in lamah_dir.glob("ID_*.csv"):
        try:
            wshed_id = int(file.stem.split('_')[1])
            watershed_ids.append(wshed_id)
        except (IndexError, ValueError):
            continue
    
    watershed_ids = sorted(watershed_ids)
    total_watersheds = len(watershed_ids)
    
    if not watershed_ids:
        print_progress("❌ No watershed files found in LAMAH directory")
        return
    
    print_progress(f"Found {total_watersheds} watersheds to process")
    print_progress(f"Reading CARRA data from: {carra_dir}")
    print_progress(f"Reading LAMAH data from: {lamah_dir}")
    print_progress(f"Saving merged data to: {output_dir}")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Process each watershed
    success = 0
    failed = 0
    start_time = datetime.now()
    
    for i, watershed_id in enumerate(watershed_ids, 1):
        print_progress(f"\nProcessing watershed {watershed_id} ({i}/{total_watersheds})")
        try:
            # Load data
            carra_df = load_carra_data(carra_dir, watershed_id)
            if carra_df is None:
                failed += 1
                continue
            
            lamah_df = load_lamah_data(lamah_dir, watershed_id)
            if lamah_df is None:
                failed += 1
                continue
            
            # Merge data
            merged_df = merge_data(lamah_df, carra_df)
            
            # Save merged data
            output_file = output_dir / f"ID_{watershed_id}.csv"
            merged_df.to_csv(output_file, sep=';', index=False)
            print_progress(f"  ✔ Saved merged data to {output_file}")
            
            # Print summary
            n_carra_cols = sum(1 for col in merged_df.columns if col.endswith('_carra'))
            print_progress(f"  Added {n_carra_cols} CARRA variables")
            print_progress(f"  Total records: {len(merged_df)}")
            
            success += 1
            
            # Estimate remaining time
            elapsed = datetime.now() - start_time
            avg_time = elapsed / i
            remaining = avg_time * (total_watersheds - i)
            print_progress(f"  Progress: {i}/{total_watersheds} ({i/total_watersheds*100:.1f}%)")
            print_progress(f"  Estimated time remaining: {str(remaining).split('.')[0]}")
            
        except Exception as e:
            print_progress(f"❌ Error processing watershed {watershed_id}: {str(e)}")
            failed += 1
            continue
    
    # Print final summary
    elapsed = datetime.now() - start_time
    print_progress("\nProcessing complete!")
    print_progress(f"Successfully processed {success} watersheds")
    if failed > 0:
        print_progress(f"Failed to process {failed} watersheds")
    print_progress(f"Total time: {str(elapsed).split('.')[0]}")

if __name__ == "__main__":
    from datetime import datetime
    main() 