import pandas as pd
from pathlib import Path
from datetime import datetime
import numpy as np

def print_progress(message):
    """Print a message with timestamp."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

# Variables that should be summed over the day (precipitation rates)
SUM_VARIABLES = [
    'lwe_precipitation_rate',
    'lwe_solid_precipitation_rate'
]

# Temperature variables that should include min/max
TEMPERATURE_VARIABLES = [
    'air_temperature_at_2m_agl'
]

def process_watershed_to_daily(input_file, output_file):
    """Convert 3-hourly watershed data to daily values."""
    print_progress(f"Processing {input_file.name}")
    
    # Read the data
    df = pd.read_csv(input_file)
    df['Time'] = pd.to_datetime(df['Time'])
    df = df.set_index('Time')
    
    # Create empty daily dataframe
    daily_df = pd.DataFrame()
    
    # Process each column with appropriate aggregation
    for column in df.columns:
        if column in SUM_VARIABLES:
            # For precipitation rates, first convert to amounts by multiplying by seconds in 3 hours
            # Then sum up the daily amounts
            amounts = df[column] * (3 * 3600)  # 3 hours in seconds
            daily_amounts = amounts.resample('D').sum()
            # Convert back to average daily rate (divide by seconds in a day)
            daily_df[column] = daily_amounts / (24 * 3600)
        
        elif column in TEMPERATURE_VARIABLES:
            # For temperature, calculate mean, min, and max
            daily_df[column] = df[column].resample('D').mean()
            daily_df[f"{column}_min"] = df[column].resample('D').min()
            daily_df[f"{column}_max"] = df[column].resample('D').max()
        
        else:
            # For all other variables (including radiation), take the daily mean
            daily_df[column] = df[column].resample('D').mean()
    
    # Basic quality checks
    print_progress(f"  Time range: {daily_df.index.min()} to {daily_df.index.max()}")
    print_progress(f"  Total days: {len(daily_df)}")
    
    # Check for gaps
    time_diff = daily_df.index.to_series().diff()
    expected_diff = pd.Timedelta(days=1)
    gaps = time_diff[time_diff > expected_diff]
    if not gaps.empty:
        print_progress("  ⚠ Found time gaps:")
        for start, diff in gaps.items():
            print_progress(f"    Gap from {start - diff} to {start} ({diff})")
    
    # Save daily series
    daily_df.to_csv(output_file)
    print_progress(f"  ✔ Saved to {output_file}")

def main():
    # Paths
    base_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds")
    input_dir = base_dir / "combined_series"
    output_dir = base_dir / "daily_series"
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print_progress("Starting to create daily watershed time series")
    print_progress(f"Reading from: {input_dir}")
    print_progress(f"Saving to: {output_dir}")
    
    # Process each watershed file
    for input_file in sorted(input_dir.glob("watershed_*.csv")):
        output_file = output_dir / input_file.name
        process_watershed_to_daily(input_file, output_file)
    
    print_progress("\nProcessing complete!")

if __name__ == "__main__":
    main() 