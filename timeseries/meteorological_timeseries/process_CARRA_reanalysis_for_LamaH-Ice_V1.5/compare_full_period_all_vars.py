import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from datetime import datetime
import os
import geopandas as gpd

def print_progress(message):
    """Print a message with timestamp."""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {message}", flush=True)

def read_new_carra_data(watershed_id):
    """Read the newly generated CARRA data for a specific watershed."""
    file_path = Path(r"C:\Users\hordurbhe\Not_backed_up\CARRA\carra_on_watersheds_all_vars\daily_series") / f"watershed_{watershed_id}.csv"
    print_progress(f"\nReading new CARRA data from: {file_path}")
    
    try:
        df = pd.read_csv(file_path)
        df['Time'] = pd.to_datetime(df['Time'])
        df = df.set_index('Time')
        
        print_progress(f"New CARRA data shape: {df.shape}")
        print_progress(f"Date range: {df.index.min()} to {df.index.max()}")
        print_progress(f"Available variables: {', '.join(df.columns)}")
        
        return df
        
    except Exception as e:
        print_progress(f"Error reading new CARRA data: {str(e)}")
        raise

def read_comparison_data(watershed_id):
    """Read the comparison data for a specific watershed."""
    file_path = Path(r"C:\Users\hordurbhe\OneDrive - Landsvirkjun\Documents\Vinna\lamah\lamah_ice\lamah_ice\A_basins_total_upstrm\2_timeseries\daily\meteorological_data") / f"ID_{watershed_id}.csv"
    print_progress(f"\nReading comparison data from: {file_path}")
    
    try:
        # Read the semicolon-separated file
        df = pd.read_csv(file_path, sep=';')
        
        # First copy the columns we want
        result_df = df[['2m_temp_rav', '2m_temp_mean', 'prec_rav', 'prec_carra',
                       '2m_qv_rav',  # specific humidity from RAV-II
                       '10m_wind_u_rav', '10m_wind_v_rav',  # wind components from RAV-II
                       '10m_wind_u', '10m_wind_v',  # wind components from ERA5-Land
                       'surf_net_solar_rad_mean', 'surf_net_therm_rad_mean',  # ERA5-Land radiation
                       'surf_outg_therm_rad_rav', 'surf_dwn_therm_rad_rav',  # RAV radiation
                       'total_et_rav', 'total_et'  # ET variables
                       ]].copy()
        
        # Calculate wind speed and direction from components
        # RAV-II wind
        result_df['wind_speed_rav'] = np.sqrt(result_df['10m_wind_u_rav']**2 + result_df['10m_wind_v_rav']**2)
        result_df['wind_dir_rav'] = (270 - np.degrees(np.arctan2(result_df['10m_wind_v_rav'], 
                                                                result_df['10m_wind_u_rav']))) % 360
        # ERA5-Land wind
        result_df['wind_speed_era5'] = np.sqrt(result_df['10m_wind_u']**2 + result_df['10m_wind_v']**2)
        result_df['wind_dir_era5'] = (270 - np.degrees(np.arctan2(result_df['10m_wind_v'], 
                                                                 result_df['10m_wind_u']))) % 360
        
        # Calculate net thermal radiation for RAV
        result_df['surf_net_therm_rad_rav'] = result_df['surf_dwn_therm_rad_rav'] - result_df['surf_outg_therm_rad_rav']
        
        # Create date index
        dates = pd.to_datetime(df[['YYYY', 'MM', 'DD']].astype(str).agg('-'.join, axis=1))
        result_df.index = dates
        
        return result_df
        
    except Exception as e:
        print_progress(f"Error reading comparison data: {str(e)}")
        raise

def plot_detailed_et_comparison(new_data, old_data, watershed_id, output_dir):
    """Create a detailed daily ET comparison for June-July 2013."""
    plt.figure(figsize=(15, 8))
    
    # Define the period of interest
    start_date = '2013-06-01'
    end_date = '2013-07-31'
    
    # Filter data to the period
    mask_new = (new_data.index >= start_date) & (new_data.index <= end_date)
    mask_old = (old_data.index >= start_date) & (old_data.index <= end_date)
    
    # Get ET data
    et_new = new_data.loc[mask_new, 'water_evaporation_amount']
    et_rav = old_data.loc[mask_old, 'total_et_rav']
    et_era5 = old_data.loc[mask_old, 'total_et']
    
    # Plot daily values
    plt.plot(et_new.index, et_new, label='New CARRA', linewidth=2)
    plt.plot(et_rav.index, et_rav, label='RAV II', linewidth=2)
    plt.plot(et_era5.index, et_era5, label='ERA5-Land', linewidth=2)
    
    # Calculate statistics for the period
    stats_text = (
        f"Mean daily ET (mm/day):\n"
        f"New CARRA: {et_new.mean():.2f}\n"
        f"RAV II: {et_rav.mean():.2f}\n"
        f"ERA5-Land: {et_era5.mean():.2f}"
    )
    
    plt.text(0.02, 0.98, stats_text,
            transform=plt.gca().transAxes,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8),
            fontsize=10)
    
    plt.title(f'Daily Evapotranspiration Comparison (June-July 2013) - Watershed {watershed_id}')
    plt.xlabel('Date')
    plt.ylabel('Evapotranspiration (mm/day)')
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True, facecolor='white', framealpha=1)
    
    # Format x-axis to show dates nicely
    plt.gcf().autofmt_xdate()
    
    plt.savefig(output_dir / f'detailed_et_watershed_{watershed_id}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_variable_groups(new_data, old_data, watershed_id, output_dir):
    """Create comparison plots for all variable groups."""
    
    # Define variable groups with their properties
    variable_groups = {
        'temperature': {
            'vars': ['air_temperature_at_2m_agl'],
            'comparison_vars': ['2m_temp_rav', '2m_temp_mean'],
            'title': 'Temperature at 2m',
            'ylabel': 'Temperature (°C)',
            'include_minmax': True
        },
        'precipitation': {
            'vars': ['lwe_precipitation_rate', 'lwe_solid_precipitation_rate'],
            'comparison_vars': ['prec_rav', 'prec_carra'],
            'title': 'Precipitation',
            'ylabel': 'Precipitation (mm/day)',
            'cumulative': True
        },
        'humidity': {
            'vars': ['specific_humidity_at_2m_agl'],
            'comparison_vars': ['2m_qv_rav'],
            'title': 'Specific Humidity',
            'ylabel': 'kg/kg'
        },
        'wind': {
            'vars': ['wind_speed_at_10m_agl'],
            'comparison_vars': ['wind_speed_rav', 'wind_speed_era5'],
            'title': 'Wind Speed at 10m',
            'ylabel': 'm/s'
        },
        'wind_direction': {
            'vars': ['wind_from_direction_at_10m_agl'],
            'comparison_vars': ['wind_dir_rav', 'wind_dir_era5'],
            'title': 'Wind Direction',
            'ylabel': 'Degrees'
        },
        'new_carra_radiation_components': {
            'vars': [
                'net_downward_shortwave_flux',
                'surface_net_thermal_radiation',
                'downward_shortwave_flux',
                'downward_longwave_flux'
            ],
            'title': 'CARRA Radiation Components',
            'ylabel': 'W/m²'
        },
        'heat_fluxes': {
            'vars': [
                'surface_upward_latent_heat_flux_due_to_sublimation',
                'surface_downward_sensible_heat_flux',
                'surface_downward_latent_heat_flux'
            ],
            'title': 'Heat Fluxes',
            'ylabel': 'W/m²'
        },
        'snow': {
            'vars': ['lwe_snow_depth'],
            'title': 'Snow Water Equivalent',
            'ylabel': 'mm w.e.'
        },
        'evapotranspiration': {
            'vars': ['water_evaporation_amount'],
            'comparison_vars': ['total_et_rav', 'total_et'],
            'title': 'Evapotranspiration',
            'ylabel': 'mm/day',
            'cumulative': True
        },
        'water_balance': {
            'vars': [
                'lwe_precipitation_rate',  # Input
                'surface_upward_latent_heat_flux_due_to_sublimation',  # Output (as mm/day)
                'water_evaporation_amount',  # Output
                'lwe_percolation_rate',  # Output
                'lwe_runoff_surface_rate'  # Output
            ],
            'title': 'Water Balance Components',
            'ylabel': 'mm/day',
            'cumulative': True
        }
    }
    
    # Process each variable group
    for group_name, group_info in variable_groups.items():
        print_progress(f"\nProcessing {group_name} plots...")
        
        # Create output directory for this group
        group_dir = output_dir / group_name
        group_dir.mkdir(exist_ok=True, parents=True)
        
        if group_name == 'precipitation':
            plot_precipitation_comparison(new_data, old_data, watershed_id, group_dir)
        elif group_name == 'temperature':
            plot_temperature_comparison(new_data, old_data, watershed_id, group_dir)
        elif group_name == 'wind_direction':
            plot_wind_direction_comparison(new_data, old_data, watershed_id, group_dir)
        elif group_name == 'evapotranspiration':
            # Create both cumulative and detailed daily plots
            plot_comparison_variables(new_data, old_data, group_info, watershed_id, group_dir)
            plot_detailed_et_comparison(new_data, old_data, watershed_id, group_dir)
        elif group_name == 'water_balance':
            water_balance_dir = output_dir / 'water_balance'
            water_balance_dir.mkdir(exist_ok=True, parents=True)
            plot_water_balance(new_data, watershed_id, water_balance_dir)
        elif 'comparison_vars' in group_info:
            plot_comparison_variables(new_data, old_data, group_info, watershed_id, group_dir)
        else:
            plot_standard_variables(new_data, group_info, watershed_id, group_dir)
    
    # Create radiation comparison plots
    radiation_dir = output_dir / 'radiation_comparison'
    radiation_dir.mkdir(exist_ok=True, parents=True)
    plot_radiation_comparison(new_data, old_data, watershed_id, radiation_dir)

def plot_precipitation_comparison(new_data, old_data, watershed_id, output_dir):
    """Create cumulative precipitation comparison plot."""
    plt.figure(figsize=(15, 10))
    
    # Find common time period
    start_date = max(new_data.index.min(), old_data.index.min())
    end_date = min(new_data.index.max(), old_data.index.max())
    
    # Filter data to common period
    mask_new = (new_data.index >= start_date) & (new_data.index <= end_date)
    mask_old = (old_data.index >= start_date) & (old_data.index <= end_date)
    
    # Get precipitation data for common period
    new_prec = new_data.loc[mask_new, 'lwe_precipitation_rate']
    old_prec_rav = old_data.loc[mask_old, 'prec_rav']
    old_prec_carra = old_data.loc[mask_old, 'prec_carra']
    
    # Calculate cumulative precipitation starting from zero
    new_cum = new_prec.cumsum()
    old_cum_rav = old_prec_rav.cumsum()
    old_cum_carra = old_prec_carra.cumsum()
    
    # Plot cumulative lines
    plt.plot(new_cum.index, new_cum, label='New CARRA', linewidth=2)
    plt.plot(old_cum_rav.index, old_cum_rav, label='RAV II', linewidth=2)
    plt.plot(old_cum_carra.index, old_cum_carra, label='Previous CARRA', linewidth=2)
    
    # Calculate total precipitation
    total_new = new_cum.iloc[-1]
    total_rav = old_cum_rav.iloc[-1]
    total_carra = old_cum_carra.iloc[-1]
    
    # Add statistics to plot
    stats_text = (
        f"Total Precipitation for period {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}:\n"
        f"New CARRA: {total_new:.1f} mm\n"
        f"RAV II: {total_rav:.1f} mm (Diff: {((total_new-total_rav)/total_rav*100):.1f}%)\n"
        f"Prev CARRA: {total_carra:.1f} mm (Diff: {((total_new-total_carra)/total_carra*100):.1f}%)"
    )
    
    plt.text(0.02, 0.98, stats_text,
            transform=plt.gca().transAxes,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8),
            fontsize=10)
    
    plt.title(f'Cumulative Precipitation Comparison - Watershed {watershed_id}', pad=20)
    plt.xlabel('Date')
    plt.ylabel('Cumulative Precipitation (mm)')
    plt.grid(True, alpha=0.3)
    plt.legend(frameon=True, facecolor='white', framealpha=1)
    
    plt.savefig(output_dir / f'cumulative_precip_watershed_{watershed_id}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # If solid precipitation is available, create a separate cumulative plot
    if 'lwe_solid_precipitation_rate' in new_data.columns:
        plt.figure(figsize=(15, 10))
        
        # Get solid precipitation for the same period
        new_solid = new_data.loc[mask_new, 'lwe_solid_precipitation_rate']
        cum_solid = new_solid.cumsum()
        
        # Plot both total and solid precipitation
        plt.plot(new_cum.index, new_cum, label='Total Precipitation', linewidth=2)
        plt.plot(cum_solid.index, cum_solid, label='Solid Precipitation', linewidth=2)
        
        # Calculate percentages
        total_solid = cum_solid.iloc[-1]
        percent_solid = (total_solid / total_new * 100) if total_new > 0 else 0
        
        # Add statistics
        stats_text = (
            f"Precipitation Totals ({start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}):\n"
            f"Total: {total_new:.1f} mm\n"
            f"Solid: {total_solid:.1f} mm ({percent_solid:.1f}% of total)"
        )
        
        plt.text(0.02, 0.98, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8),
                fontsize=10)
        
        plt.title(f'Total vs Solid Precipitation - Watershed {watershed_id}', pad=20)
        plt.xlabel('Date')
        plt.ylabel('Cumulative Precipitation (mm)')
        plt.grid(True, alpha=0.3)
        plt.legend(frameon=True, facecolor='white', framealpha=1)
        
        plt.savefig(output_dir / f'solid_precip_watershed_{watershed_id}.png', 
                    dpi=300, bbox_inches='tight')
        plt.close()

def plot_temperature_comparison(new_data, old_data, watershed_id, output_dir):
    """Create detailed temperature comparison plots."""
    plt.figure(figsize=(15, 10))
    
    # Calculate weekly means for smoother plotting
    temp_new = new_data['air_temperature_at_2m_agl'].resample('W').mean()
    temp_rav = old_data['2m_temp_rav'].resample('W').mean()
    temp_mean = old_data['2m_temp_mean'].resample('W').mean()
    
    # Plot temperature comparison
    plt.plot(temp_new.index, temp_new, label='New CARRA')
    plt.plot(temp_rav.index, temp_rav, label='RAV II')
    plt.plot(temp_mean.index, temp_mean, label='ERA5-Land')
    
    # Add min/max if available
    if 'air_temperature_at_2m_agl_min' in new_data.columns:
        temp_min = new_data['air_temperature_at_2m_agl_min'].resample('W').min()
        temp_max = new_data['air_temperature_at_2m_agl_max'].resample('W').max()
        plt.fill_between(temp_min.index, temp_min, temp_max, 
                        alpha=0.2, label='New CARRA Min-Max Range')
    
    plt.title(f'Weekly Mean Temperature at 2m - Watershed {watershed_id}')
    plt.ylabel('Temperature (°C)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(output_dir / f'temperature_comparison_watershed_{watershed_id}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_comparison_variables(new_data, old_data, group_info, watershed_id, output_dir):
    """Create plots for variables that have comparison data."""
    plt.figure(figsize=(15, 10))
    
    # Find common time period
    start_date = max(new_data.index.min(), old_data.index.min())
    end_date = min(new_data.index.max(), old_data.index.max())
    
    # Filter data to common period
    mask_new = (new_data.index >= start_date) & (new_data.index <= end_date)
    mask_old = (old_data.index >= start_date) & (old_data.index <= end_date)
    
    # Plot new CARRA data
    for var in group_info['vars']:
        if var in new_data.columns:
            if group_info.get('cumulative', False):
                # For cumulative variables like ET
                values = new_data.loc[mask_new, var].cumsum()
                plt.plot(values.index, values, label='New CARRA', alpha=0.7)
            else:
                # For regular variables, use weekly means
                weekly_mean = new_data.loc[mask_new, var].resample('W').mean()
                plt.plot(weekly_mean.index, weekly_mean, label='New CARRA', alpha=0.7)
    
    # Plot comparison data
    for var in group_info['comparison_vars']:
        if var in old_data.columns:
            if group_info.get('cumulative', False):
                # For cumulative variables like ET
                values = old_data.loc[mask_old, var].cumsum()
                label = 'RAV II' if 'rav' in var else 'ERA5-Land'
                plt.plot(values.index, values, label=label, alpha=0.7)
            else:
                # For regular variables, use weekly means
                weekly_mean = old_data.loc[mask_old, var].resample('W').mean()
                label = 'RAV II' if 'rav' in var else 'ERA5-Land'
                plt.plot(weekly_mean.index, weekly_mean, label=label, alpha=0.7)
    
    plt.title(f"{group_info['title']} - Watershed {watershed_id}")
    plt.ylabel(group_info['ylabel'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add statistics for cumulative variables
    if group_info.get('cumulative', False):
        # Calculate totals
        new_total = new_data.loc[mask_new, group_info['vars'][0]].sum()
        stats_text = f"Total for period {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}:\n"
        stats_text += f"New CARRA: {new_total:.1f} mm\n"
        
        for var in group_info['comparison_vars']:
            if var in old_data.columns:
                old_total = old_data.loc[mask_old, var].sum()
                label = 'RAV II' if 'rav' in var else 'ERA5-Land'
                diff_percent = ((new_total - old_total) / old_total * 100) if old_total != 0 else float('inf')
                stats_text += f"{label}: {old_total:.1f} mm (Diff: {diff_percent:.1f}%)\n"
        
        plt.text(0.02, 0.98, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8),
                fontsize=10)
    
    plt.savefig(output_dir / f'{group_info["title"].lower().replace(" ", "_")}_watershed_{watershed_id}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_standard_variables(new_data, group_info, watershed_id, output_dir):
    """Create plots for standard variable groups."""
    vars_to_plot = [var for var in group_info['vars'] if var in new_data.columns]
    
    if not vars_to_plot:
        return
    
    plt.figure(figsize=(15, 10))
    
    for var in vars_to_plot:
        # Calculate weekly means for smoother plotting
        weekly_mean = new_data[var].resample('W').mean()
        plt.plot(weekly_mean.index, weekly_mean, 
                 label=var.replace('_', ' ').title(), alpha=0.7)
    
    plt.title(f"{group_info['title']} - Watershed {watershed_id}")
    plt.ylabel(group_info['ylabel'])
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(output_dir / f'{group_info["title"].lower()}_watershed_{watershed_id}.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_wind_direction_comparison(new_data, old_data, watershed_id, output_dir):
    """Create enhanced wind direction comparison plots."""
    plt.figure(figsize=(20, 10))
    
    # Find common time period
    start_date = max(new_data.index.min(), old_data.index.min())
    end_date = min(new_data.index.max(), old_data.index.max())
    
    # Filter data to common period
    mask_new = (new_data.index >= start_date) & (new_data.index <= end_date)
    mask_old = (old_data.index >= start_date) & (old_data.index <= end_date)
    
    # Get wind data
    carra_dir = new_data.loc[mask_new, 'wind_from_direction_at_10m_agl']
    carra_speed = new_data.loc[mask_new, 'wind_speed_at_10m_agl']
    
    # RAV-II comparison
    plt.subplot(121)
    rav_dir = old_data.loc[mask_old, 'wind_dir_rav']
    rav_speed = old_data.loc[mask_old, 'wind_speed_rav']
    
    # Create scatter plot colored by wind speed (minimum of both datasets)
    min_speeds = np.minimum(carra_speed, rav_speed)
    scatter = plt.scatter(carra_dir, rav_dir, 
                         c=min_speeds, cmap='viridis', 
                         alpha=0.3, s=10,
                         vmin=0, vmax=15)
    
    plt.colorbar(scatter, label='Wind Speed (m/s)')
    plt.plot([0, 360], [0, 360], 'k--', alpha=0.5)  # 1:1 line
    
    # Calculate and display statistics
    stats = calculate_wind_stats(carra_dir, carra_speed, rav_dir, rav_speed)
    if stats:
        stats_text = (
            f"Statistics (wind speed ≥ 2 m/s):\n"
            f"Mean angular diff: {stats['mean_diff']:.1f}°\n"
            f"Median angular diff: {stats['median_diff']:.1f}°\n"
            f"RMSE: {stats['rmse']:.1f}°\n"
            f"Within 45°: {stats['percent_agreement_45']:.1f}%\n"
            f"N = {stats['num_filtered']} of {stats['num_total']}"
        )
        plt.text(0.02, 0.98, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8),
                fontsize=10)
    
    plt.xlabel('CARRA Wind Direction (degrees)')
    plt.ylabel('RAV II Wind Direction (degrees)')
    plt.title('Wind Direction: CARRA vs RAV II')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    
    # ERA5-Land comparison
    plt.subplot(122)
    era5_dir = old_data.loc[mask_old, 'wind_dir_era5']
    era5_speed = old_data.loc[mask_old, 'wind_speed_era5']
    
    # Create scatter plot colored by wind speed
    min_speeds = np.minimum(carra_speed, era5_speed)
    scatter = plt.scatter(carra_dir, era5_dir, 
                         c=min_speeds, cmap='viridis', 
                         alpha=0.3, s=10,
                         vmin=0, vmax=15)
    
    plt.colorbar(scatter, label='Wind Speed (m/s)')
    plt.plot([0, 360], [0, 360], 'k--', alpha=0.5)  # 1:1 line
    
    # Calculate and display statistics
    stats = calculate_wind_stats(carra_dir, carra_speed, era5_dir, era5_speed)
    if stats:
        stats_text = (
            f"Statistics (wind speed ≥ 2 m/s):\n"
            f"Mean angular diff: {stats['mean_diff']:.1f}°\n"
            f"Median angular diff: {stats['median_diff']:.1f}°\n"
            f"RMSE: {stats['rmse']:.1f}°\n"
            f"Within 45°: {stats['percent_agreement_45']:.1f}%\n"
            f"N = {stats['num_filtered']} of {stats['num_total']}"
        )
        plt.text(0.02, 0.98, stats_text,
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8),
                fontsize=10)
    
    plt.xlabel('CARRA Wind Direction (degrees)')
    plt.ylabel('ERA5-Land Wind Direction (degrees)')
    plt.title('Wind Direction: CARRA vs ERA5-Land')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'wind_direction_comparison_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()

def plot_radiation_comparison(new_data, old_data, watershed_id, output_dir):
    """Create radiation comparison plots."""
    
    # Create subdirectories for different radiation components
    solar_dir = output_dir / 'solar'
    thermal_dir = output_dir / 'thermal'
    solar_dir.mkdir(exist_ok=True, parents=True)
    thermal_dir.mkdir(exist_ok=True, parents=True)
    
    # Find common time period
    start_date = max(new_data.index.min(), old_data.index.min())
    end_date = min(new_data.index.max(), old_data.index.max())
    
    # Filter data to common period
    mask_new = (new_data.index >= start_date) & (new_data.index <= end_date)
    mask_old = (old_data.index >= start_date) & (old_data.index <= end_date)
    
    # Net Solar Radiation Comparison
    plt.figure(figsize=(15, 10))
    plt.plot(new_data.loc[mask_new].index, 
             new_data.loc[mask_new, 'net_downward_shortwave_flux'].rolling('7D').mean(),
             label='CARRA', alpha=0.7)
    plt.plot(old_data.loc[mask_old].index, 
             old_data.loc[mask_old, 'surf_net_solar_rad_mean'].rolling('7D').mean(),
             label='ERA5-Land', alpha=0.7)
    
    plt.title(f'Net Solar Radiation Comparison - Watershed {watershed_id}')
    plt.ylabel('Net Solar Radiation (W/m²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(solar_dir / f'net_solar_radiation_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Net Thermal Radiation Comparison
    plt.figure(figsize=(15, 10))
    plt.plot(new_data.loc[mask_new].index, 
             new_data.loc[mask_new, 'surface_net_thermal_radiation'].rolling('7D').mean(),
             label='CARRA', alpha=0.7)
    plt.plot(old_data.loc[mask_old].index, 
             old_data.loc[mask_old, 'surf_net_therm_rad_mean'].rolling('7D').mean(),
             label='ERA5-Land', alpha=0.7)
    plt.plot(old_data.loc[mask_old].index, 
             old_data.loc[mask_old, 'surf_net_therm_rad_rav'].rolling('7D').mean(),
             label='RAV II', alpha=0.7)
    
    plt.title(f'Net Thermal Radiation Comparison - Watershed {watershed_id}')
    plt.ylabel('Net Thermal Radiation (W/m²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(thermal_dir / f'net_thermal_radiation_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Downward Thermal Radiation Comparison
    plt.figure(figsize=(15, 10))
    plt.plot(new_data.loc[mask_new].index, 
             new_data.loc[mask_new, 'downward_longwave_flux'].rolling('7D').mean(),
             label='CARRA', alpha=0.7)
    plt.plot(old_data.loc[mask_old].index, 
             old_data.loc[mask_old, 'surf_dwn_therm_rad_rav'].rolling('7D').mean(),
             label='RAV II', alpha=0.7)
    
    plt.title(f'Downward Thermal Radiation Comparison - Watershed {watershed_id}')
    plt.ylabel('Downward Thermal Radiation (W/m²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(thermal_dir / f'downward_thermal_radiation_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate and plot outgoing thermal radiation for CARRA
    plt.figure(figsize=(15, 10))
    
    # CARRA: outgoing = net + downward (since net = downward - outgoing)
    carra_outgoing = (new_data.loc[mask_new, 'surface_net_thermal_radiation'] + 
                     new_data.loc[mask_new, 'downward_longwave_flux'])
    
    plt.plot(new_data.loc[mask_new].index, 
             carra_outgoing.rolling('7D').mean(),
             label='CARRA', alpha=0.7)
    plt.plot(old_data.loc[mask_old].index, 
             old_data.loc[mask_old, 'surf_outg_therm_rad_rav'].rolling('7D').mean(),
             label='RAV II', alpha=0.7)
    
    # Calculate statistics
    carra_mean = carra_outgoing.mean()
    rav_mean = old_data.loc[mask_old, 'surf_outg_therm_rad_rav'].mean()
    diff_percent = ((carra_mean - rav_mean) / rav_mean * 100)
    
    stats_text = (
        f"Mean Values (W/m²):\n"
        f"CARRA: {carra_mean:.1f}\n"
        f"RAV II: {rav_mean:.1f}\n"
        f"Difference: {diff_percent:.1f}%"
    )
    
    plt.text(0.02, 0.98, stats_text,
            transform=plt.gca().transAxes,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8),
            fontsize=10)
    
    plt.title(f'Outgoing Thermal Radiation Comparison - Watershed {watershed_id}')
    plt.ylabel('Outgoing Thermal Radiation (W/m²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(thermal_dir / f'outgoing_thermal_radiation_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()

    # Add new plot: All CARRA longwave components
    plt.figure(figsize=(15, 10))
    
    # Get the components
    downward_lw = new_data.loc[mask_new, 'downward_longwave_flux']
    net_lw = new_data.loc[mask_new, 'surface_net_thermal_radiation']
    outgoing_lw = downward_lw + net_lw  # since net = downward - outgoing
    
    # Plot with 7-day rolling means
    plt.plot(downward_lw.index, downward_lw.rolling('7D').mean(), 
             label='Downward LW', color='blue', alpha=0.7)
    plt.plot(outgoing_lw.index, outgoing_lw.rolling('7D').mean(), 
             label='Outgoing LW', color='red', alpha=0.7)
    plt.plot(net_lw.index, net_lw.rolling('7D').mean(), 
             label='Net LW', color='black', alpha=0.7)
    
    # Calculate statistics
    stats_text = (
        f"Mean Values (W/m²):\n"
        f"Downward LW: {downward_lw.mean():.1f}\n"
        f"Outgoing LW: {outgoing_lw.mean():.1f}\n"
        f"Net LW: {net_lw.mean():.1f}"
    )
    
    plt.text(0.02, 0.98, stats_text,
            transform=plt.gca().transAxes,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8),
            fontsize=10)
    
    plt.title(f'CARRA Longwave Radiation Components - Watershed {watershed_id}')
    plt.ylabel('Radiation (W/m²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add a horizontal line at y=0
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    plt.savefig(thermal_dir / f'carra_longwave_components_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()

def calculate_angular_difference(dir1, dir2):
    """Calculate the smallest angular difference between two directions (0-360)."""
    diff = np.abs(dir1 - dir2)
    return np.minimum(diff, 360 - diff)

def calculate_wind_stats(carra_dir, carra_speed, other_dir, other_speed, min_speed=2.0):
    """Calculate wind direction statistics considering wind speed."""
    # Filter for significant wind speeds
    mask = (carra_speed >= min_speed) & (other_speed >= min_speed)
    
    if not mask.any():
        return None
    
    carra_dir_filtered = carra_dir[mask]
    other_dir_filtered = other_dir[mask]
    
    # Calculate angular differences
    differences = calculate_angular_difference(carra_dir_filtered, other_dir_filtered)
    
    stats = {
        'mean_diff': np.mean(differences),
        'median_diff': np.median(differences),
        'rmse': np.sqrt(np.mean(differences**2)),
        'num_total': len(carra_dir),
        'num_filtered': mask.sum(),
        'percent_agreement_45': (differences <= 45).mean() * 100,  # % within 45 degrees
    }
    
    return stats

def calculate_statistics(new_data, old_data):
    """Calculate comparison statistics for all variables."""
    stats = {}
    
    # Temperature statistics
    if 'air_temperature_at_2m_agl' in new_data.columns and '2m_temp_rav' in old_data.columns:
        common_dates = new_data.index.intersection(old_data.index)
        if len(common_dates) > 0:
            new_temp = new_data.loc[common_dates, 'air_temperature_at_2m_agl']
            old_temp = old_data.loc[common_dates, '2m_temp_rav']
            
            stats['temperature'] = {
                'correlation': new_temp.corr(old_temp),
                'mean_difference': (new_temp - old_temp).mean(),
                'rmse': np.sqrt(((new_temp - old_temp) ** 2).mean())
            }
    
    # Precipitation statistics
    if 'lwe_precipitation_rate' in new_data.columns:
        for old_col in ['prec_rav', 'prec_carra']:
            if old_col in old_data.columns:
                common_dates = new_data.index.intersection(old_data.index)
                if len(common_dates) > 0:
                    new_prec = new_data.loc[common_dates, 'lwe_precipitation_rate']
                    old_prec = old_data.loc[common_dates, old_col]
                    
                    stats[f'precipitation_{old_col}'] = {
                        'correlation': new_prec.corr(old_prec),
                        'total_new': new_prec.sum(),
                        'total_old': old_prec.sum(),
                        'difference_percent': ((new_prec.sum() - old_prec.sum()) / old_prec.sum() * 100)
                    }
    
    return stats

def get_water_year(date):
    """Convert date to water year (Oct-Sep)."""
    if date.month >= 10:
        return date.year + 1
    return date.year

def get_glacier_fraction(watershed_id):
    """Get glacier fraction for a watershed from attributes file."""
    try:
        attributes_file = Path(r"C:\Users\hordurbhe\OneDrive - Landsvirkjun\Documents\Vinna\lamah\lamah_ice\lamah_ice\A_basins_total_upstrm\1_attributes\Catchment_attributes.csv")
        df = pd.read_csv(attributes_file, sep=';')
        watershed_info = df[df['id'] == watershed_id].iloc[0]
        return watershed_info['glac_fra']
    except Exception as e:
        print_progress(f"Could not determine glacier coverage: {str(e)}")
        return None

def plot_water_balance(new_data, watershed_id, output_dir):
    """Create water balance analysis plots."""
    
    # Check glacier coverage
    glacier_fraction = get_glacier_fraction(watershed_id)
    glacier_text = ""
    balance_text = ""
    if glacier_fraction is not None:
        glacier_percent = glacier_fraction * 100
        glacier_text = f" (Glaciated: {glacier_percent:.1f}%)"
        if glacier_percent > 1:
            print_progress(f"⚠ Warning: Watershed {watershed_id} is {glacier_percent:.1f}% glaciated.")
            print_progress("  Water balance analysis may not be meaningful due to glacier storage/melt.")
            if glacier_percent > 10:
                print_progress("  Skipping water balance analysis due to high glacier coverage.")
                return
    
    # Get components
    components = {
        'Precipitation': {'var': 'lwe_precipitation_rate', 'color': 'blue'},
        'Snow Sublimation': {'var': 'surface_upward_latent_heat_flux_due_to_sublimation', 'color': 'purple'},
        'Evaporation': {'var': 'water_evaporation_amount', 'color': 'red'},
        'Percolation': {'var': 'lwe_percolation_rate', 'color': 'brown'},
        'Surface Runoff': {'var': 'lwe_runoff_surface_rate', 'color': 'green'}
    }
    
    # Ensure all components are present
    available_components = {name: info for name, info in components.items() 
                          if info['var'] in new_data.columns}
    
    if len(available_components) < len(components):
        missing = set(components.keys()) - set(available_components.keys())
        print_progress(f"Warning: Missing water balance components: {missing}")
    
    # Add water year to the data
    data_with_wy = new_data.copy()
    data_with_wy['water_year'] = [get_water_year(date) for date in data_with_wy.index]
    
    # Check which water years are complete
    # A complete water year should have data from October through September
    complete_years = []
    for year in data_with_wy['water_year'].unique():
        year_data = data_with_wy[data_with_wy['water_year'] == year]
        # Get unique months in this water year
        months_present = set()
        for date in year_data.index:
            # For October-December, we need previous year
            if date.month >= 10:
                if date.year == year - 1:
                    months_present.add(date.month)
            # For January-September, we need the water year itself
            else:
                if date.year == year:
                    months_present.add(date.month)
        
        # Check if we have all months (1-12)
        if len(months_present) == 12:
            complete_years.append(year)
    
    if not complete_years:
        print_progress(f"Warning: No complete water years found for watershed {watershed_id}")
        return
    
    # Filter data to only complete water years
    data_with_wy = data_with_wy[data_with_wy['water_year'].isin(complete_years)]
    
    # Calculate water year totals
    plt.figure(figsize=(15, 8))
    
    water_years = sorted(complete_years)
    x_positions = np.arange(len(water_years))
    
    # Plot annual totals for each component
    for name, info in available_components.items():
        yearly_totals = data_with_wy.groupby('water_year')[info['var']].sum()
        plt.plot(x_positions, yearly_totals, 
                label=name, color=info['color'],
                marker='o', markersize=6, linewidth=2)
    
    # Calculate and plot water balance
    if 'Precipitation' in available_components:
        yearly_p = data_with_wy.groupby('water_year')[components['Precipitation']['var']].sum()
        outputs = []
        for name, info in available_components.items():
            if name != 'Precipitation':
                outputs.append(data_with_wy.groupby('water_year')[info['var']].sum())
        
        total_output = sum(outputs)
        balance = yearly_p - total_output
        
        # Calculate average balance as percentage of precipitation
        avg_balance = balance.mean()
        avg_precip = yearly_p.mean()
        balance_percent = (avg_balance / avg_precip) * 100
        balance_text = f", Avg. Net Balance: {balance_percent:.1f}% of P"
        
        # Plot balance line
        plt.plot(x_positions, balance, label='Net Balance',
                color='black', marker='s', markersize=6,
                linewidth=2, linestyle='-')
    
    # Customize plot
    plt.title(f'Annual Water Balance Components - Watershed {watershed_id}{glacier_text}{balance_text}')
    plt.xlabel('Water Year')
    plt.ylabel('Annual Total (mm/year)')
    
    # Set x-axis labels as water years
    plt.xticks(x_positions, [f'{year-1}/{year}' for year in water_years], rotation=45)
    
    # Add horizontal line at y=0
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Adjust legend
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust layout
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_dir / f'annual_water_balance_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also create the monthly pattern plot
    plt.figure(figsize=(15, 8))
    
    # Calculate and plot monthly means for each component
    for name, info in available_components.items():
        # Group by month but calculate means for each water year
        monthly_means = data_with_wy.groupby([data_with_wy.index.month, 'water_year'])[info['var']].mean()
        # Reorder to start with October
        months = list(range(10, 13)) + list(range(1, 10))
        monthly_means = monthly_means.unstack()  # Convert to months x years
        monthly_means = monthly_means.reindex(months)  # Reorder months
        
        # Plot mean and range
        mean_values = monthly_means.mean(axis=1)
        std_values = monthly_means.std(axis=1)
        
        plt.plot(range(len(months)), mean_values, 
                label=name, color=info['color'], alpha=0.7)
        plt.fill_between(range(len(months)), 
                        mean_values - std_values,
                        mean_values + std_values,
                        color=info['color'], alpha=0.2)
    
    plt.title(f'Monthly Mean Water Balance Components - Watershed {watershed_id}{glacier_text}{balance_text}')
    plt.xlabel('Month')
    plt.ylabel('Mean Rate (mm/day)')
    
    # Set x-axis labels to month names starting from October
    month_labels = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep']
    plt.xticks(range(len(month_labels)), month_labels)
    
    # Add horizontal line at y=0
    plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    
    # Adjust legend and layout
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save monthly pattern plot
    plt.savefig(output_dir / f'monthly_pattern_water_balance_watershed_{watershed_id}.png',
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Configuration
    output_dir = Path("full_period_comparison_all_vars")
    output_dir.mkdir(exist_ok=True)
    
    # Get list of first 5 watershed IDs
    comparison_dir = Path(r"C:\Users\hordurbhe\OneDrive - Landsvirkjun\Documents\Vinna\lamah\lamah_ice\lamah_ice\A_basins_total_upstrm\2_timeseries\daily\meteorological_data")
    watershed_ids = sorted([int(f.stem.split('_')[1]) for f in comparison_dir.glob("ID_*.csv")])[:5]
    
    print_progress(f"Processing first {len(watershed_ids)} watersheds")
    
    # Process each watershed
    for watershed_id in watershed_ids:
        print_progress(f"\nProcessing watershed {watershed_id}")
        try:
            # Read data
            new_data = read_new_carra_data(watershed_id)
            old_data = read_comparison_data(watershed_id)
            
            # Create plots
            plot_variable_groups(new_data, old_data, watershed_id, output_dir)
            
            # Calculate statistics
            stats = calculate_statistics(new_data, old_data)
            
            # Save statistics
            stats_file = output_dir / f"statistics_watershed_{watershed_id}.txt"
            with open(stats_file, 'w') as f:
                f.write(f"Statistics for Watershed {watershed_id}\n")
                f.write("=" * 40 + "\n\n")
                for var_name, var_stats in stats.items():
                    f.write(f"{var_name}:\n")
                    for stat_name, value in var_stats.items():
                        f.write(f"  {stat_name}: {value:.3f}\n")
                    f.write("\n")
            
            print_progress(f"✓ Successfully processed watershed {watershed_id}")
            
        except Exception as e:
            print_progress(f"✗ Error processing watershed {watershed_id}: {str(e)}")
            continue

if __name__ == '__main__':
    main() 