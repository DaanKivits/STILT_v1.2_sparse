#!/usr/bin/env python
"""Script to download ERA data for STILT processing.
Written by Auke van der Woude,
last modified 17-05-2023"""

# Import modules.
# We download data from ERA, so import the climate data store api
# For this, a key is required!
import cdsapi
import datetime as dtm
import subprocess
from dateutil import rrule
from calendar import monthrange
import os
#from cdo import Cdo
#cdo = Cdo()

# Set up the connection to the server
c = cdsapi.Client()
# Download data for these days
starttime = dtm.datetime(2021,1,1)
endtime = dtm.datetime(2022,2,1)
for date in rrule.rrule(rrule.MONTHLY, dtstart=starttime, until=endtime):
    # We'll download 2 different files per day:
    # One for pressure levels and one at surface level
    # We then combine these files
    files = []

    ndays_month = monthrange(date.year, date.month)[1]
    # This is the outfile to write to: pressure levels!
    outfile = '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-{}.pressure.grb'.format(date.strftime('%Y%m'))
    if not os.path.exists(outfile):
        # Download the data
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'format': 'grb',
                'variable': [
                    'relative_humidity', 'temperature', 'geopotential',
                    'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',
                ],
                'pressure_level': [
                    '1', '3', '50',
                    '100', '250', '350',
                    '500', '650', '700',
                    '750', '800', '825',
                    '850', '875', '900',
                    '925', '950', '975',
                    '1000',
                ],
                'year': '{}'.format(date.year),
                'month': '{0:02}'.format(date.month),
                'day': ['%02d'%(day+1) for day in range(ndays_month)],
                'area': [
                    72, -15, 33, 35,
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
            },
        outfile,
        )

    # We save this filename to combine later
    files.append(outfile)

    # Now retrieve the single-level data:
    outfile = '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-{}.single.grb'.format(date.strftime('%Y%m'))
    if not os.path.exists(outfile):
        c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grb',
            'variable': [
                '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
                'surface_pressure', 'total_precipitation', 'geopotential',
                'boundary_layer_height', 'friction_velocity', 'surface_net_solar_radiation',
                'surface_sensible_heat_flux',
            ],
                'year': '{}'.format(date.year),
                'month': '{0:02}'.format(date.month),
                'day': ['%02d'%(day+1) for day in range(ndays_month)],
                'area': [
                    72, -15, 33, 35,
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
            },
        outfile,
        )
    # Append this out file to combine
    files.append(outfile)

    # we'll combine these files:
    files = ' '.join(files)
    # Into this file:
    newfile = f'/projects/0/ctdas/PARIS/DATA/meteo/STILT/scripts/combined_auke/ERA5-{date:%Y%m}.grb'
    # With grb, we can just concatenate them
    command = f'cdo merge {files} {newfile}'
    p = subprocess.Popen(command, shell=True)
    p.communicate()

    # Now, we have to make single days out of the .grb file
    # We can do that by cdo seltime,1/24 ..
    for i in range(ndays_month):
        d = i + 1
        starthour = i * 24 + 1 # indexing starts at 1
        endhour = starthour + 23 # inclusive
        dailyfile = f'/projects/0/ctdas/PARIS/DATA/meteo/STILT/scripts/combined_auke/ERA5-{date:%Y%m}{d:02}.grb'
        command = f'cdo seltimestep,{starthour}/{endhour} {newfile} {dailyfile}'
        p = subprocess.Popen(command, shell=True)
        p.communicate()
        #cdo.seltimestep(f'{starthour}/{endhour}', input=outfile, output=f'./out/ERA5-{date:%Y%m}{d:02}.grb')
    assert False
