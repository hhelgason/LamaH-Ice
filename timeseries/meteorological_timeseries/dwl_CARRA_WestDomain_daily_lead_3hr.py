# This script downloads CARRA data from the CDS
import cdsapi

c = cdsapi.Client()
for year in range(1992,2010):

    print(year)

    for month in range(1,12+1):
        for day in range(1,31+1):
            print(year,month)
            filename = "CARRA_WestDomain_daily_lead_3hr_" + str(year) + "_" + str(month) + "_" + str(day) + ".nc"
            print(filename)
            try:

                c.retrieve(
                    'reanalysis-carra-single-levels',
                    {
                        'format': 'netcdf',
                        'domain': 'west_domain',
                        'level_type': 'surface_or_atmosphere',
                        'variable': [
                            '10m_wind_direction', '10m_wind_speed', '2m_relative_humidity',
                            '2m_temperature', 'snow_density', 'snow_depth_water_equivalent',
                            'surface_pressure', 'surface_runoff', 'surface_solar_radiation_downwards',
                            'surface_thermal_radiation_downwards', 'total_precipitation',
                        ],
                        'product_type': 'forecast',
                        'time': [
                            '00:00', '03:00', '06:00',
                            '09:00', '12:00', '15:00',
                            '18:00', '21:00',
                        ],
                        'year': str(year),
                        'leadtime_hour': [
                            '3',
                        ],
                    'month': str(month),
                    'day': str(day),
                        },
                        filename)

            except:
                pass
