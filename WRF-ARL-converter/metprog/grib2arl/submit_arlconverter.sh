#!/bin/bash
#SBATCH -p thin
#SBATCH -t 00:30:00 
#SBATCH -n 64
#SBATCH -N 1
#SBATCH --mail-user=daan.kivits@wur.nl
#SBATCH --mail-type=FAIL,END

source /home/dkivits/.bashrc 

# Only if using virtual environment!
module purge
module load 2021
module load CDO/1.9.10-gompi-2021a
module load NCO/5.0.1-foss-2021a
module load Anaconda3/2021.05
module load netCDF/4.8.0-gompi-2021a 

source activate cte-hr-env

python /home/dkivits/STILT/WRF-ARL-converter/metprog/grib2arl/process_multiple.py