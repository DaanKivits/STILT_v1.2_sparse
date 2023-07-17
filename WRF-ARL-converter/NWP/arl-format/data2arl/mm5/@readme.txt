./data2arl/mm5				Last Revised: 09 Mar 2004  
----------------------------------------------------------------------

The conversion program for MM5 Version 3 output files is simply run 
with the command mm5toarl.exe.  The MM5 input file can be defined on
the command line:

	mm5toarl MMOUT_DOMAIN3

or it will default to file name "fort.10".  The output file will be
created according to the domain number, so that for the above 
example the converted output file will be called:

	DATA3.ARL

Note that the conversion program rounds the output times to the 
nearest 15 minutes.  This is required because when MM5 is run on
a multiprocessor system, output times may differ be a few minutes 
with each cycle.
