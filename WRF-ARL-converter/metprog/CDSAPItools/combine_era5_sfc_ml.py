#!/usr/bin/env python
"""Script to download ERA data for STILT processing.
Written by Auke van der Woude,
last modified 17-05-2023"""

# Import modules.
# We download data from ERA, so import the climate data store api
# For this, a key is required!
import datetime as dtm
import subprocess
from dateutil import rrule
from calendar import monthrange
import os

# Download data for these days
starttime = dtm.datetime(2021,1,1)
endtime = dtm.datetime(2021,12,31)
for date in rrule.rrule(rrule.MONTHLY, dtstart=starttime, until=endtime):
    # We'll download 2 different files per day:
    # One for pressure levels and one at surface level
    # We then combine these files
    files = []

    ndays_month = monthrange(date.year, date.month)[1]
    # This is the outfile to write to: pressure levels!
    ml = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/ML/{}.grb'.format(date.strftime('%Y_%m'))

    # Now retrieve the single-level data:
    sfc = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/SFC/{}.grb'.format(date.strftime('%Y_%m'))
    
    # Append this out file to combine
    files.extend([ml, sfc])

    # we'll combine these files:
    files = ' '.join(files)

    # Into this new intermediate nc file:
    intermediatefile = f'/projects/0/ctdas/PARIS/DATA/meteo/STILT/combined/ERA5intermediate-{date:%Y%m}.grb'

    # Merge nc files into one
    command = f'cdo merge {files} {intermediatefile}'
    p = subprocess.Popen(command, shell=True)
    p.communicate()
    
    """
    # If not already done, extract variable names from intermediate file
    if not os.path.exists(variabledef_file):
        command = f'cdo vardes {intermediatefile} > {variabledef_file}'
        p = subprocess.Popen(command, shell=True)
        p.communicate()
        
    # Into this new GRIB file:
    newfile = f'/projects/0/ctdas/PARIS/DATA/meteo/STILT/combined/ERA5-{date:%Y%m}.grb'
    
    # Check if the directory exists; if not, create it
    if not os.path.exists(newfile):
        # Create the directory
        os.makedirs(os.path.dirname(newfile), exist_ok=True)

    command = f'cdo -v -f --history grb -copy {intermediatefile} {newfile}'
    p = subprocess.Popen(command, shell=True)
    p.communicate()

    # Change variable names in grib file with ones extracted earlier
    command = f'cdo chname,{variabledef_file} {newfile} {newfile}'
"""

    # Now, we have to make single days out of the .grb file
    # We can do that by cdo seltime,1/24 ..
    for i in range(ndays_month):
        d = i + 1
        starthour = i * 24 + 1 # indexing starts at 1
        endhour = starthour + 23 # inclusive
        dailyfile = f'/projects/0/ctdas/PARIS/DATA/meteo/STILT/combined/ERA5-{date:%Y%m}{d:02}.grb'
        command = f'cdo seltimestep,{starthour}/{endhour} {intermediatefile} {dailyfile}'
        p = subprocess.Popen(command, shell=True)
        p.communicate()

    # Remove the intermediate files
    #os.remove([newfile, catfile, intermediatefile])
