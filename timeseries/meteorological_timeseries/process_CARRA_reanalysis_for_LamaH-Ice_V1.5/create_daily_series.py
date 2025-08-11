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
    'lwe_solid_precipitation_rate',
    'lwe_percolation_rate',
    'lwe_runoff_surface_rate'
]

# Temperature variables that should include min/max
TEMPERATURE_VARIABLES = [
    'air_temperature_at_2m_agl'
]

# Energy flux variables that need conversion to daily water equivalent
ENERGY_FLUX_VARIABLES = [
    'surface_upward_latent_heat_flux_due_to_sublimation'
]

# Evaporation variables (already in mm per 3-hour period)
EVAPORATION_VARIABLES = [
    'water_evaporation_amount'
]

# Radiation variables that need sign flipping (to make downward fluxes positive)
FLIP_SIGN_VARIABLES = [
    'surface_net_thermal_radiation'
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
            # For precipitation rates (mm/hour), multiply by 3 to get mm per 3-hour period
            # Then sum up the daily amounts to get mm/day
            amounts = df[column] * 3  # Convert from mm/hour to mm/3hours
            daily_df[column] = amounts.resample('D').sum()
        
        elif column in TEMPERATURE_VARIABLES:
            # For temperature, calculate mean, min, and max
            daily_df[column] = df[column].resample('D').mean()
            daily_df[f"{column}_min"] = df[column].resample('D').min()
            daily_df[f"{column}_max"] = df[column].resample('D').max()
        
        elif column in ENERGY_FLUX_VARIABLES:
            # For energy flux in W/m²:
            # 1. Convert power (W/m²) to energy (J/m²) for each 3-hour period
            # 2. Sum up to get daily total energy (J/m²/day)
            # 3. Convert daily energy to water mass using latent heat of sublimation
            # Note: In CARRA, negative flux means energy loss (sublimation)
            # so we need to flip the sign to get positive values for sublimation
            
            # First convert W/m² to J/m² for each 3-hour period
            seconds_per_period = 3 * 3600  # 3 hours in seconds
            energy_3h = df[column] * seconds_per_period  # Convert from W/m² to J/m² per 3-hour period
            
            # Sum to daily totals (J/m²/day)
            # resample('D').sum() adds up all eight 3-hour periods in the day
            daily_energy = energy_3h.resample('D').sum()
            
            # Convert energy to water equivalent
            # Flip sign so positive = sublimation (mass loss)
            latent_heat_sublimation = 2.834e6  # J/kg
            daily_df[column] = -1 * daily_energy / latent_heat_sublimation  # kg/m² = mm/day
            
            # Add some quality checks for this variable
            if abs(daily_df[column]).max() > 5:  # mm/day seems very high
                print_progress(f"  ⚠ Warning: Very high sublimation values detected")
                print_progress(f"    Max absolute value: {abs(daily_df[column]).max():.2f} mm/day")
                print_progress(f"    This might indicate a unit conversion issue in the source data")
        
        elif column in EVAPORATION_VARIABLES:
            # For evaporation, values are already in mm per 3-hour period
            # Just sum up to get daily totals
            # Flip sign so positive means water loss (evaporation)
            daily_df[column] = -1 * df[column].resample('D').sum()
            
            # Add some quality checks
            if abs(daily_df[column]).max() > 10:  # mm/day seems very high
                print_progress(f"  ⚠ Warning: Very high evaporation values detected")
                print_progress(f"    Max absolute value: {abs(daily_df[column]).max():.2f} mm/day")
                print_progress(f"    This might indicate a unit conversion issue in the source data")
        
        elif column in FLIP_SIGN_VARIABLES:
            # For radiation variables where we want downward fluxes to be positive
            daily_df[column] = -1 * df[column].resample('D').mean()
        
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
    base_dir = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds_all_vars")
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