#bash_run.py

""" Run wrfnc2arl to create ARL from wrfout for multiple WRF output files and times

Author: Ingrid Super
Last revisions: 20-1-2017"""

import os
import subprocess
import shutil
import glob
import sys
import datetime as dtm
import netCDF4 as nc
####################################specify the following##########################################

ndom=2
execfile='./wrfnc2arl'
wrfpath='/projects/0/ctdas/RINGO/wrfout/run04/'
arlpath='./'

###################################################################################################

print("Processing", ndom, "domains")
for k in range(ndom):
    domain = str(k + 1)
    dstilt = ndom - k
    wrffiles = sorted(glob.glob(wrfpath + 'wrfout_d0' + domain + '*'))
    for file in wrffiles:
        print("processing", file)
        mf = nc.Dataset(file)
        Times = mf.variables['Times'][:]
        mf.close()
        for it, tim_ in enumerate(Times):
            tim = b''.join(tim_).decode()
            tim = dtm.datetime.strptime(tim, '%Y-%m-%d_%H:%M:%S')
            dat = tim.date()
            dat2 = dat.strftime('%Y%m%d')
            newdate = (tim + dtm.timedelta(hours=1)).strftime('%Y%m%d_%H')
            args = [execfile,'-T',str(it+1),'-P','var_sample',file]
            err_flag = 1
            while err_flag == 1:
                a = subprocess.call(args,shell=False)
                outfile1 = glob.glob('CFG*')[0]
                outfile2 = glob.glob('DATA*')[0]
                with open(outfile2, 'rb') as f:
                    first_line = f.readline()
                dum = first_line[14: 18]
                if dum == b'INDX': 
                    err_flag = 0
            newfile1 = 'wrfout_d0' + str(dstilt) + '_' + newdate + '.'+outfile1
            newfile2 = 'wrfout_d0' + str(dstilt) + '_' + newdate + '.'+outfile2
            if not os.path.exists('arl/' + dat2):
                os.makedirs('arl/' + dat2)
            dummy = shutil.move(outfile1, os.path.join('arl/' + dat2, newfile1))
            dummy = shutil.move(outfile2, os.path.join('arl/' + dat2, newfile2))
