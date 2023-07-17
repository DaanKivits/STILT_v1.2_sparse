import subprocess
import glob
import shutil
import datetime as dtm
import os

os.chdir('/home/dkivits/STILT/WRF-ARL-converter/metprog/grib2arl/')

filedir = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/combined_20lvls/'
#filedir = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/testfiles/'
files = sorted(glob.glob(filedir + '*.grb2'))

for f in files:
    datestr = f.split('/')[-1][5:-4]

    #outfile = '/projects/0/ctdas/PARIS/DATA/meteo/STILT/combined/arl/ecmw.{}.arl'.format(datestr)
    outfile = './arl/ecmw.{}.arl'.format(datestr)
    #outfile = './arl/ecmw.test.arl'.format(datestr)

    #command = ['/home/dkivits/STILT/WRF-ARL-converter/metprog/grib2arl/grib2arl'] + [r'-i' + f ]
    command = ['./grib2arl'] + [r'-i' + f ]
    command = ' '.join(command)
    print(command)
    p = subprocess.Popen(command, shell=True)
    p.communicate()
    shutil.move('DATA.ARL', outfile)
        
