# %% [markdown]


# %%
# If not running from the same folder, uncomment these two lines and adapt the path
import os   
os.chdir('/home/dkivits/STILT/WRF-ARL-converter/metprog/CDSAPItools')
from CDSAPItools import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray
from windrose import WindroseAxes
import matplotlib.cm as cm

# %% Concatenate all monthly data into a time-continuous dataset
# https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14

# Read and concatenate all nc files
# ds = xarray.open_mfdataset(f'{out_path}/*.nc',combine = 'by_coords', concat_dim="time")

# # Save as a new netCDF file
# ds.to_netcdf(f'{out_path}/{year_start}{month_start:02d}-{year_end}{month_end:02d}.nc')

# %% PARIS - STILT Model level data

year_start  = 2021
year_end    = 2021
month_start = 1
month_end   = 12
area = [72, -15, 33, 35]
# level_list = [
#                     '1', '3', '50',
#                     '100', '250', '350',
#                     '500', '650', '700',
#                     '750', '800', '825',
#                     '850', '875', '900',
#                     '925', '950', '975',
#                     '1000',
#                 ]
level_list = [
                    '350',
                    '500', '650', '700',
                    '750', '800', '825',
                    '850', '875', '900',
                    '925', '950', '975',
                    '1000',
                ]
variable_list = [
                    'relative_humidity', 'temperature', 'geopotential',
                    'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',
                ]
time_list = [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
                ]
prepend = False
dataset = 'reanalysis-era5-pressure-levels'
out_path = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/ML'

rDict = {
    'product_type': 'reanalysis',
    'area': area,
    'pressure_level': level_list,
    'variable': variable_list,
    'stream':'oper',
    'time': time_list,
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grb')
#checkERA(out_path, prepend)

# %% PARIS - STILT Surface data

dataset = 'reanalysis-era5-single-levels'

out_path = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/SFC'
variable_list = [
                '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
                'surface_pressure', 'total_precipitation', 'geopotential',
                'boundary_layer_height', 'friction_velocity', 'surface_net_solar_radiation',
                'surface_sensible_heat_flux',
            ]

rDict = {
    'product_type': 'reanalysis',
    'area': area,
    'variable': variable_list,
    'stream':'oper',
    'time': time_list,
    'grid': "0.25/0.25",
}

# %%
submitERA(out_path, year_start, year_end, month_start, month_end, dataset, rDict, format='grb')
#checkERA(out_path, prepend)

# %%
