#!/usr/bin/env python

import cdsapi
import datetime as dtm
import subprocess
from dateutil import rrule
from calendar import monthrange
import os 

c = cdsapi.Client()
starttime = dtm.datetime(2022,1,1)
endtime = dtm.datetime(2022,2,1)
for date in rrule.rrule(rrule.MONTHLY, dtstart=starttime, until=endtime):
    files = []
    ndays_month = monthrange(date.year, date.month)[1]
    outfile = '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-{}.pressure.nc'.format(date.strftime('%Y%m%d'))
    
    if not os.path.exists(outfile):
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'format': 'nc',
                'variable': [
                    'relative_humidity', 'temperature',
                    'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',
                    'geopotential'
                    
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
    else:
        print('File {} already exists'.format(outfile))
        files.append(outfile)

#    command = 'cdo -f nc setgridtype,regular {} {}'.format(outfile, outfile.replace('.grib', '.nc'))
#    p = subprocess.Popen(command, shell=True)
#    p.communicate()
    #files.append(outfile)


    outfile = '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-{}.single.nc'.format(date.strftime('%Y%m%d'))
    if not os.path.exists(outfile):
        c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'nc',
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
    else:
        print('File {} already exists'.format(outfile))

    files.append(outfile)

    files = ' '.join(files)
    newfile =  '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-{}.nc'.format(date.strftime('%Y%m%d'))
    command = 'cdo merge {} {}'.format(files, newfile)
    p = subprocess.Popen(command, shell=True)
    p.communicate()

    command = 'cdo -v -f grb -copy {} {}'.format(newfile, newfile.replace('.nc', '.grb')) 
    p = subprocess.Popen(command, shell=True)
    p.communicate()
