#***************************************************************************************************
# source all required R functions
#***************************************************************************************************
# $Id: sourceall.r,v 1.22 2009-11-25 08:51:03 gerbig Exp $
#---------------------------------------------------------------------------------------------------
sourceall <- function(mode='multi', sparse=TRUE) {
    go <- graphics.off

    # Default sourcepath is current directory
    if (!is.element("sourcepath",objects())) sourcepath <- "./"

    source(paste(sourcepath, "assignr.r", sep=""))
    source(paste(sourcepath, "col.grads.r", sep=""))
    source(paste(sourcepath, "combine.met.r", sep=""))
    source(paste(sourcepath, "create_sparse_footprints.r", sep=""))
    source(paste(sourcepath, "day.of.week.r", sep=""))
    source(paste(sourcepath, "distance.r", sep=""))
    source(paste(sourcepath, "existsr.r", sep=""))
    source(paste(sourcepath, "flttrack.r", sep=""))
    source(paste(sourcepath, "footplot.r", sep=""))
    source(paste(sourcepath, "get.CarbonTracker.netcdf.r", sep=""))
    source(paste(sourcepath, "get.fireBARCA.netcdf.r", sep=""))
    source(paste(sourcepath, "get.fossEU.netcdf.r", sep=""))
    source(paste(sourcepath, "get.edgar.time.r", sep=""))
    source(paste(sourcepath, "get.Emis.netcdf.r", sep=""))
    source(paste(sourcepath, "get.Emis.annual.netcdf.r", sep=""))
    source(paste(sourcepath, "get.GEMS_CO.netcdf.r", sep=""))
    source(paste(sourcepath, "get.LMDZ.netcdf.r", sep=""))
    source(paste(sourcepath, "get.MACC_CO.netcdf.r", sep=""))
    source(paste(sourcepath, "get.MACC_CO2.netcdf.r", sep=""))
    source(paste(sourcepath, "get.mean.traj.r", sep=""))
    source(paste(sourcepath, "get.modis.netcdf.r", sep=""))
    source(paste(sourcepath, "get.vegfrac.netcdf.r", sep=""))
    source(paste(sourcepath, "get.vprm.dim.r", sep=""))
    source(paste(sourcepath, "get.TM3.bin.r", sep=""))
    source(paste(sourcepath, "get.TM3.netcdf.r", sep=""))
    source(paste(sourcepath, "get.TM3CH4.bin.r", sep=""))
    source(paste(sourcepath, "getgridp.r", sep=""))
    source(paste(sourcepath, "getmetfile.r", sep=""))
    source(paste(sourcepath, "getr.r", sep=""))
    source(paste(sourcepath, "id2pos.r", sep=""))
    source(paste(sourcepath, "image.plot.fix.r", sep=""))
    source(paste(sourcepath, "image.plot.plt.fix.r", sep=""))
    source(paste(sourcepath, "imagell.r", sep=""))
    source(paste(sourcepath, "is.inf.r", sep=""))
    source(paste(sourcepath, "julian.r", sep=""))
    source(paste(sourcepath, "leap.year.r", sep=""))
    source(paste(sourcepath, "lsr.r", sep=""))
    source(paste(sourcepath, "measured.co2.r", sep=""))
    source(paste(sourcepath, "memory.size.r", sep=""))
    source(paste(sourcepath, "month.day.year.r", sep=""))
    source(paste(sourcepath, "motif.r", sep=""))
    source(paste(sourcepath, "my.file.exists.r", sep=""))
    source(paste(sourcepath, "pos2id.r", sep=""))
    source(paste(sourcepath, "project.nice.r", sep=""))
    source(paste(sourcepath, "read.asc.r", sep=""))
    source(paste(sourcepath, "read.bground.r", sep=""))
    source(paste(sourcepath, "rmr.r", sep=""))
    source(paste(sourcepath, "rp2ll.r", sep=""))
    source(paste(sourcepath, "stdev.r", sep=""))
    source(paste(sourcepath, "stilt3D.r", sep=""))
    source(paste(sourcepath, "Trajec.r", sep=""))
    source(paste(sourcepath, "Trajecflux.r", sep=""))
    source(paste(sourcepath, "Trajecmulti.r", sep=""))
    source(paste(sourcepath, "Trajecvprm.r", sep=""))
    source(paste(sourcepath, "trajwind.r", sep=""))
    source(paste(sourcepath, "unix.r", sep=""))
    source(paste(sourcepath, "unix.shell.r", sep=""))
    source(paste(sourcepath, "weekdayhr.r", sep=""))

    if (mode == 'single' & sparse == FALSE) {
    # Switch these on if doing a single-site run:
    source(paste(sourcepath, "Trajecmod.r", sep=""))
    source(paste(sourcepath, "Trajecfoot.r", sep=""))
    source("setStiltparam.r")

    } else if (mode == 'single' & sparse == TRUE) {
    # Switch these on if doing a single-site run:
    source(paste(sourcepath, "Trajecmod_sparse.r", sep=""))
    source(paste(sourcepath, "Trajecfoot_sparse.r", sep=""))
    source("setStiltparam.r")

    } else if (mode == 'multi' & sparse == FALSE){
    # Switch these on if doing a multi-site run:
    source(paste(sourcepath, "Trajecmod_loop.r", sep=""))
    source(paste(sourcepath, "Trajecfoot.r", sep=""))
    source("setStiltparam_loop.r")

    } else if (mode == 'multi' & sparse == TRUE){
    # Switch these on if running for a multi-site sparse footprints run:
    source(paste(sourcepath, "Trajecmod_loop_sparse.r", sep=""))
    source(paste(sourcepath, "Trajecfoot_sparse.r", sep=""))
    source("setStiltparam_loop.r")
    
    }
}