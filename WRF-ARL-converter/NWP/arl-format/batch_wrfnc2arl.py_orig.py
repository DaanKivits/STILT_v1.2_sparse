#sh_run.py

""" Run wrfnc2arl to create ARL from wrfout for multiple WRF output files and times

Author: Ingrid Super
Last revisions: 20-1-2017"""

import datetime as dtm
import os
import netCDF4 as nc
import subprocess
import shutil
import glob
####################################specify the following#############################################################
ndom=1
execfile='./wrfnc2arl'
wrfpath='/home/awoude/WRF41/WRF/run/'

#######################################################################################################################

for k in range(ndom):
    domain=k+1
    dstilt=ndom-k
    wrffiles = glob.glob(wrfpath+'wrfout_d0'+str(domain)+"*")
    print(wrffiles)
#    wrffiles=[os.path.join(filename) for filename in os.listdir(wrfpath) if filename.startswith('wrfout_d%02d'%(domain))]
    for file in wrffiles:
        mf=nc.Dataset(file)
        Times=mf.variables['Times'][:]
        mf.close()
        nlen=len(Times)
        for it in range(nlen):
            tim=dtm.datetime(int(''.join(Times[it][:4])),int(''.join(Times[it][5:7])),int(''.join(Times[it][8:10])),int(''.join(Times[it][11:13])),int(''.join(Times[it][14:16])))
            dat=dtm.date(int(''.join(Times[it][:4])),int(''.join(Times[it][5:7])),int(''.join(Times[it][8:10])))
            dat2='%04d%02d%02d'%(int(''.join(Times[it][:4])),int(''.join(Times[it][5:7])),int(''.join(Times[it][8:10])))
            id=str(tim.hour+1)
            id2=id.zfill(2)
            args=[execfile,'-T',str(it+1),'-P','var_sample',file]
            err_flag=1
            while err_flag==1:
                subprocess.call(args,shell=False)
                outfile1 = [os.path.join(filename) for filename in os.listdir(wrfpath) if filename.startswith('CFG')]
                outfile2 = [os.path.join(filename) for filename in os.listdir(wrfpath) if filename.startswith('DATA')]
                with open(outfile2[0],'r') as f:
                    first_line=f.readline()
                dum=first_line[14:18]
                if dum=='INDX':
                    err_flag=0
                newfile1 = 'wrfout_d0'+str(dstilt)+'_'+str(dat)+'_'+id2+'.'+outfile1[0]
                newfile2 = 'wrfout_d0'+str(dstilt)+'_'+str(dat)+'_'+id2+'.'+outfile2[0]
                if not os.path.exists('arl/'+dat2):
                    os.makedirs('arl/'+dat2)
                dummy = shutil.move(outfile1[0],os.path.join('arl/'+dat2,newfile1))
                dummy = shutil.move(outfile2[0],os.path.join('arl/'+dat2,newfile2))

