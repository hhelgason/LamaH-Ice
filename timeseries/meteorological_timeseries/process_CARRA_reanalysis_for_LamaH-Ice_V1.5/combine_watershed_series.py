import pandas as pd
from pathlib import Path
import os
import glob
from datetime import datetime
import numpy as np

def print_progress(message):
    """Print a message with timestamp."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def handle_extreme_values(df, column='surface_net_thermal_radiation', threshold=10000):
    """Handle extreme values in radiation and sublimation columns by replacing them with previous valid values."""
    if column not in df.columns:
        return df, 0
    
    # Create a copy to avoid modifying the original
    df = df.copy()
    
    # Set thresholds for different variables
    if column == 'surface_net_thermal_radiation':
        value_threshold = 10000  # Original threshold for thermal radiation
    elif column == 'surface_upward_latent_heat_flux_due_to_sublimation':
        value_threshold = 300   # Threshold for sublimation in W/m² (about 3 mm/day)
    else:
        value_threshold = threshold
    
    # Find extreme values (using absolute value)
    extreme_mask = abs(df[column]) > value_threshold
    n_extremes = extreme_mask.sum()
    
    if n_extremes > 0:
        # Replace extreme values with NaN first
        df.loc[extreme_mask, column] = np.nan
        # Forward fill NaN values (use previous valid value)
        df[column] = df[column].fillna(method='ffill')
        # Backward fill any remaining NaN values at the start
        df[column] = df[column].fillna(method='bfill')
        return df, n_extremes
    
    return df, 0

def combine_watershed_csvs(monthly_dir, output_dir):
    """Combine monthly CSV files for each watershed into continuous time series."""
    monthly_dir = Path(monthly_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Get list of all watershed IDs by scanning the first month's directory
    first_month = sorted(monthly_dir.glob("[0-9]" * 6))[0]
    watershed_files = list((first_month / "csv").glob("watershed_*.csv"))
    watershed_ids = [int(f.stem.split('_')[1]) for f in watershed_files]
    
    print_progress(f"Found {len(watershed_ids)} watersheds to process")
    
    # Process each watershed
    for wshed_id in sorted(watershed_ids):
        print_progress(f"Processing watershed {wshed_id}")
        
        # Collect all monthly files for this watershed
        monthly_data = []
        for month_dir in sorted(monthly_dir.glob("[0-9]" * 6)):
            csv_file = month_dir / "csv" / f"watershed_{wshed_id}.csv"
            if csv_file.exists():
                try:
                    df = pd.read_csv(csv_file)
                    df['Time'] = pd.to_datetime(df['Time'])
                    df = df.set_index('Time')
                    monthly_data.append(df)
                except Exception as e:
                    print_progress(f"❌ Error reading {csv_file}: {str(e)}")
        
        if not monthly_data:
            print_progress(f"❌ No data found for watershed {wshed_id}")
            continue
        
        # Combine all months
        combined_df = pd.concat(monthly_data, axis=0)
        
        # Sort by time and remove any duplicates
        combined_df = combined_df.sort_index()
        combined_df = combined_df[~combined_df.index.duplicated(keep='first')]
        
        # Handle extreme values in thermal radiation
        combined_df, n_extremes_therm = handle_extreme_values(combined_df, 'surface_net_thermal_radiation')
        if n_extremes_therm > 0:
            print_progress(f"  Fixed {n_extremes_therm} extreme thermal radiation values")
        
        # Handle extreme values in sublimation
        combined_df, n_extremes_subl = handle_extreme_values(combined_df, 'surface_upward_latent_heat_flux_due_to_sublimation')
        if n_extremes_subl > 0:
            print_progress(f"  Fixed {n_extremes_subl} extreme sublimation values")
        
        # Basic quality checks
        print_progress(f"Watershed {wshed_id} summary:")
        print_progress(f"  Time range: {combined_df.index.min()} to {combined_df.index.max()}")
        print_progress(f"  Total records: {len(combined_df)}")
        
        # Check for gaps
        time_diff = combined_df.index.to_series().diff()
        expected_diff = pd.Timedelta(hours=3)  # CARRA data is 3-hourly
        gaps = time_diff[time_diff > expected_diff]
        if not gaps.empty:
            print_progress("  ⚠ Found time gaps:")
            for start, diff in gaps.items():
                print_progress(f"    Gap from {start - diff} to {start} ({diff})")
        
        # Save combined series
        output_file = output_dir / f"watershed_{wshed_id}.csv"
        combined_df.to_csv(output_file)
        print_progress(f"  ✔ Saved to {output_file}")

def main():
    # Paths
    monthly_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds_all_vars")
    output_dir = monthly_dir / "combined_series"
    
    print_progress("Starting to combine watershed time series")
    print_progress(f"Reading from: {monthly_dir}")
    print_progress(f"Saving to: {output_dir}")
    
    combine_watershed_csvs(monthly_dir, output_dir)
    
    print_progress("\nProcessing complete!")

if __name__ == "__main__":
    main() 