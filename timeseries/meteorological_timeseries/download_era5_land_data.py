import cdsapi
import numpy as np
c = cdsapi.Client()

variables = ["10m_u_component_of_wind", "10m_v_component_of_wind", "surface_net_solar_radiation", "surface_net_thermal_radiation", "snow_depth_water_equivalent","volumetric_soil_water_layer_1",
             "volumetric_soil_water_layer_2", "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4","forecast_albedo", "snow_albedo", "snow_cover", "snow_density", "snowmelt"] # "2m_dewpoint_temperature", "2m_temperature", "potential_evaporation", "total_precipitation", "total_evaporation", "surface_pressure", "leaf_area_index_low_vegetation", "leaf_area_index_high_vegetation",

years = np.arange(1980,2022)

for var in variables:
    for year in years:
        c.retrieve(
    'reanalysis-era5-land',
    {
        'variable': var,
        'year': str(year),
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            66.75, -24.8, 63.3,
            -13.3,
        ],
        'format': 'netcdf',
    },
    'data/ERA5-Land_%s_%s.nc' % (var,str(year)))
