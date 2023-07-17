import subprocess
import glob
from shutil import move
import datetime as dtm
import os

filedir = '/home/dkivits/STILT/WRF-ARL-converter/NWP/arl-format/Metfiles/intermediatedata/ERA5-2022'

files = glob.glob(filedir + '*.grb')

for f in files:
    datestr = f.split('/')[-1][5:-4]

    outfile = './arl/ecmw.{}.arl'.format(datestr)
    
    command = ['./grib2arl'] + [r'-i' + f ]
    command = ' '.join(command)
    print(command)
    p = subprocess.Popen(command, shell=True)
    p.communicate()
    move('DATA.ARL', outfile)
        
