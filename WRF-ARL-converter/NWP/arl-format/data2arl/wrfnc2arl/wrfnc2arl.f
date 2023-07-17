! Use emacs f90-mode:   -*- mode:f90  -*-

! $Id: wrfnc2arl.f,v 1.3 2008/01/17 18:41:09 trn Exp $

!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: WRFNC2ARL      DECODE netcdf WRF MODEL FIELD FOR HYSPLIT
!
! Written by: Thomas Nehrkorn, AER, Inc. 
! Based on:
! WRFGRIB2ARL
! MAIN PROGRAM: RSMS2ARL      DECODE GRIB RSM MODEL FIELD FOR HYSPLIT
!   PRGMMR: DRAXLER          ORG: R/E/AR      DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE written at AER ...
!     converts a WRF model netcdf file to arl packed format.
!     processes only one time period per execution.  Multiple input files
!     may be specified on the command line.  Input data are
!     assumed to be already on a conformal map projection on model levels.
!     Only packing to ARL format and units conversion is required.
!
! PROGRAM HISTORY LOG:
! $Log: wrfnc2arl.f,v $
! Revision 1.3  2008/01/17 18:41:09  trn
! Support use of WPS geolocation routines
!
! Revision 1.3  2008/01/17 18:41:09  trn
! Support use of WPS geolocation routines
!
! Revision 1.2  2007/11/20 22:37:49  trn
! Support ncdf files with multiple file periods
!
! Revision 1.1  2007/02/07 19:46:38  trn
! Initial version of netcdf_to_arl converter
!
!
! USAGE:  WRFNC2ARL -P par_fname [-V iprint] [-C G] [-S G] [-C L] [-S L] [-T itime] 
!                   [-M mdlid] [-G geoloc_type] NETCDF_FILE_NAME [NETCDF_FILE_NAME2] [...]
!   INPUT ARGUMENT LIST:
!     par_fname - filename of text file containing parameters (variables) info
!                 (see below for format)
!     -C G/L     - combine grids and/or levels (default is to combine both)
!     -S G/L     - keep grids and/or levels separate
!                  grids: if combined, all fields are output to the unstaggered (mass) file;
!                         all fields will be dimensioned by the dimensions of the mass grid, but
!                         augmented by an additional row and column
!                  levels: if combined, only full levels are output, but half-level variables
!                          are output as if they were located at the next higher full level
!     -M mdlid   - model ID to output to cfg, data files
!     -V iprint - printout control (default:2, summary info; less:0-1; 3,4: more detail; 8,9: full detail )
!     -T itime  - time period to extract (default:1)
!     -G geoloc_type - type of mapping routines to use (one of: cmap, wps)
!     NETCDF_FILE_NAME - multiple netcdf files may be defined (only one at a time is supported for now)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     unit 10 - used for par_fname
!     NETCDF_DATA - netcdf input data files use netcdf library calls
!   OUTPUT FILES:
!     unit 20 DATA_MASS.WRF - ARL packed data output file for mass points
!     unit 21 DATA_STAX.WRF - ARL packed data output file for x-staggered points
!     unit 22 DATA_STAY.WRF - ARL packed data output file for y-staggered points
!     unit 30 - used for CFG_MASS.WRF, CFG_STAX.WRF, CFG_STAY.WRF (temp file area)
!     unit 40 WRFTIME - time indicator file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

PROGRAM WRFNC2ARL

      use read_header_module
      use xtract_module

      implicit none       

      include 'netcdf.inc'

      LOGICAL   :: ftest, lhelp
      integer, parameter :: lchar=256
      CHARACTER (len=lchar) :: option, FNAME
      CHARACTER (len=lchar) :: par_fname, vlabel
      character (len=4) :: mdlid, geoloc_type, saved_mdlid
      INTEGER :: handle, status

      integer :: iarg, narg, iargc
      integer :: iutext=10, iumass=20, iustax=21, iustay=22, iucfg=30, iutime=40
      integer :: ierr, idatlin, ndatlin, iprint, itime

      integer, parameter :: nlevtyp=3, p_sfc=1, p_half=2, p_full=3
      integer, parameter :: ngridtyp=3, p_mass=1, p_stax=2, p_stay=3
      integer :: nvars(nlevtyp,ngridtyp), ndim1v
      character (len=lchar), allocatable :: vlabels(:,:,:)
      real, allocatable :: vcnvrt(:,:,:),vadd(:,:,:)
      character (len=4), allocatable :: vnames(:,:,:)
      character (len=4) :: vname
      character (len=1) :: sta_xy, sta_z
      logical :: have_half, have_full, combine_levs, combine_grids, use_16
      integer :: i, iltyp, igtyp, iladd, i_23d, iladd_soil, iladd_save
      real :: cnvrt,add

!=>check for command line arguments
      NARG=IARGC()

      lhelp = narg .eq. 0

!     process options
      iarg = 0
      par_fname=' '
      mdlid = 'AWRF'  !'DWRF' for 16-bit files, '9WRF/NWRF' for WRF with hybrid vertical coordinate
      geoloc_type = 'WPS' !Use WPS geolocation routines by default, 'CMAP' for old hymodelc built-ins
      combine_levs=.TRUE.
      combine_grids=.TRUE.
      iprint=2
      itime=1
      do while (iarg .lt. narg .and. .not. lhelp)
         iarg=iarg+1
         CALL GETARG(iarg,option)
         if (option(1:1) .ne. '-') exit
         iarg=iarg+1
         CALL GETARG(iarg,fname)
         SELECT CASE (option)
         CASE ('-V','-v')
            read(fname,*,iostat=ierr) iprint
            if (ierr .ne. 0) then
               write (*,*) 'Error decoding iprint from -V option:',trim(fname)
               lhelp = .TRUE.
            endif
         CASE ('-C','-c')
            SELECT CASE (fname)
               CASE ('G','g')
                  combine_grids=.TRUE.
               CASE ('L','l')
                  combine_levs=.TRUE.
               CASE DEFAULT
                  write (*,*) 'Invalid -C option:',trim(fname)
                  lhelp=.TRUE.
               END SELECT
         CASE ('-S','-s')
            SELECT CASE (fname)
               CASE ('G','g')
                  combine_grids=.FALSE.
               CASE ('L','l')
                  combine_levs=.FALSE.
               CASE DEFAULT
                  write (*,*) 'Invalid -S option:',trim(fname)
                  lhelp=.TRUE.
               END SELECT
         CASE ('-P','-p')
            par_fname=fname
         CASE ('-M','-m')
            mdlid=fname(1:4)
         CASE ('-T','-t')
            read(fname,*,iostat=ierr) itime
            if (ierr .ne. 0) then
               write (*,*) 'Error decoding itime from -T option:',trim(fname)
               lhelp = .TRUE.
            endif
         CASE ('-G','-g')
            geoloc_type=fname
         CASE DEFAULT
            write (*,*) 'Invalid option:',trim(option)
            lhelp=.TRUE.
         END SELECT
      end do
      lhelp = lhelp .or. par_fname .eq. ' '
      lhelp = lhelp .or. iarg .gt. narg
      use_16 = .TRUE.
      if (mdlid(1:1) .ne. 'D' .and. mdlid(1:1) .ne. 'N') use_16 = .FALSE.
      write (*,'(a)') 'Running wrfnc2arl $Revision: 1.3 $ with the following options:'
      write (*,'(a,t30,a)') 'Parameter codes file [-P]:',trim(par_fname)
      write (*,'(a,t30,i10)') 'Extracting time period [-T]:',itime
      write (*,'(a,t30,l1)') 'Combining grids [-C/-S G]:',combine_grids
      write (*,'(a,t30,l1)') 'Combining levels [-C/-S L]:',combine_levs
      write (*,'(a,t30,a)') 'Model ID [-M mdlid]:',trim(mdlid)
      write (*,'(a)') '  NOTE: '
      write (*,'(a)') '   use AWRF/9WRF for 8 bit, DWRF/NWRF for 16 bit precision'
      write (*,'(a)') '   AWRF/DWRF is for sigma-p (terrain-following, TF) vertical coordinate output, '
      write (*,'(a)') '   9WRF/NWRF is for hybrid vertical coordinate output '
      write (*,'(a)') '   for WRF output with hybrid_opt=-1, 9WRF/NWRF will get changed to AWRF/DWRF'
      write (*,'(a)') '    (not compiled for hybrid or pre-v39, no hybrid constants arrays)'
      write (*,'(a)') '   for WRF output with hybrid_opt=0, mdlid will remain unchanged'
      write (*,'(a)') '    (compiled for hybrid, but run with TF, will have hybrid constants arrays that'
      write (*,'(a)') '     are functionally equivalent to TF, so either option will work)'
      write (*,'(a)') '   for WRF output with hybrid_opt=2, AWRF/DWRF will get changed to 9WRF/NWRF'
      write (*,'(a)') '    (compiled and run for hybrid, need hybrid constants arrays)'
      write (*,'(a,t30,a)') 'Geolocation type [-G type]:',trim(geoloc_type)
      write (*,'(a,t30,a)') '   - can be one of: cmap, wps'
      write (*,'(a,t50,l1)') 'Using 16 bits [mdlid(1:1) .eq. "D" or "N"]:',use_16
      write (*,'(a,t30,i10)') 'Printout control [-V]:',iprint

      if (lhelp) then
         write (*,'(a)') 'Usage: '
         write (*,'(2a)') '  WRFNC2ARL -P par_fname [-V iprint]', &
              &' [-C G] [-C L] [-S G] [-S L] [-T itime] [-G geoloc_type]', &
              &' NETCDF_FILE_NAME [NETCDF_FILE_NAME2] [...]'
         stop
      endif

! Read in variable info from text file:
! Format: 
!   Header:
!         1st line: contains integer (list-directed I/O) ndatlin: 
!                   number of data lines to follow the header
!         subsequent lines: any number of header lines (ignored), followed
!         by line with the string ENDHEADER in columns 1-9
!   Data lines: (ndatlin lines, each containing in list-directed I/O format)
!     4-character variable name, enclosed in ''
!     1-character xy-staggering indicator, one of ('x', 'y', or ' '), case-insensitive
!     1-character z-staggering indicator, one of ('z', or ' '), case-insensitive
!     integer dimensionality indicator, either 2 or 3 (for 2d/3d fields)
!     character variable label in netcdf file, enclosed in ''
!     real conversion factor to be applied to data before output
      call read_header(iutext,par_fname,ndatlin,ierr,iprint)
      if (ierr .ne. 0) stop 'Error reading header from par_fname'
      ndim1v=ndatlin
      allocate(vlabels(ndatlin,nlevtyp,ngridtyp), vnames(ndatlin,nlevtyp,ngridtyp), &
           & vcnvrt(ndatlin,nlevtyp,ngridtyp), vadd(ndatlin,nlevtyp,ngridtyp), &
           & stat=ierr)
      if (ierr .ne. 0) stop 'Allocate vlabel... failure'
      nvars(:,:) = 0
      have_half = .false.
      have_full = .false.
      do idatlin=1,ndatlin
         read (iutext,*,iostat=ierr) vname,sta_xy,sta_z,i_23d,vlabel,cnvrt,add
         if (ierr .ne. 0) stop 'Error reading datlin from par_fname'
         igtyp=0
         if (sta_xy .eq. ' ') then
            igtyp=p_mass
         elseif (sta_xy .eq. 'X' .or. sta_xy .eq. 'x') then
            igtyp=p_stax
         elseif (sta_xy .eq. 'Y' .or. sta_xy .eq. 'y') then
            igtyp=p_stay
         else
            stop 'invalid value for sta_xy in par_fname'
         endif
         iltyp=0
         if (sta_z .eq. ' ') then
            iltyp=p_half
         elseif (sta_z .eq. 'Z' .or. sta_z .eq. 'z') then
            iltyp=p_full
         else
            stop 'invalid value for sta_z in par_fname'
         endif
         if (i_23d .eq. 2) then
            iltyp = p_sfc
         elseif (i_23d .eq. 3) then
            have_half = have_half .or. iltyp .eq. p_half
            have_full = have_full .or. iltyp .eq. p_full
         else
            stop 'invalid value for i_23d in par_fname'
         endif
         nvars(iltyp,igtyp) = nvars(iltyp,igtyp) + 1
         vnames(nvars(iltyp,igtyp),iltyp,igtyp) = vname
         vlabels(nvars(iltyp,igtyp),iltyp,igtyp) = vlabel
         vcnvrt(nvars(iltyp,igtyp),iltyp,igtyp) = cnvrt
         vadd(nvars(iltyp,igtyp),iltyp,igtyp) = add
      end do
      close (iutext)

      if (.not. (have_half .or. have_full)) stop 'need to specify non-sfc variables'
      
! Read in and process GRIB filenames (at present, xtract can only handle one file)
      do while (iarg .le. narg)
         CALL GETARG(iarg,fname)
         INQUIRE(file=fname,exist=ftest)
         if(ftest)then
            status = NF_OPEN(fname, nf_nowrite, handle)
            if (status .ne. nf_noerr) then
               write (*,*) 'Error opening ',trim(fname),' : ',nf_strerror(status)
            else
               write(*,'(2a)') 'Started processing: ',trim(fname)
               call xtract(handle, iumass, iustax, iustay, iucfg, iutime, iprint, &
                    & nlevtyp, p_sfc, p_half, p_full, ngridtyp, p_mass, p_stax, p_stay, &
                    & nvars, ndim1v, &
                    & vlabels, lchar, vnames, vcnvrt, vadd, &
                    & combine_grids, combine_levs, mdlid, use_16, itime, geoloc_type )
               status=nf_close(handle)
               if (status .ne. nf_noerr) write (*,*) 'Error closing ',trim(fname), &
                    & ' : ',nf_strerror(status)
            endif
         else
            write(*,*)'File not found:',fname
         end if
         iarg = iarg + 1
      end do

!     close out time period and write index record(s)
      INQUIRE(unit=iumass,opened=ftest)
      if (ftest) then
         call pakndx_16(iumass,use_16)
         CLOSE (iumass)
      endif
      INQUIRE(unit=iustax,opened=ftest)
      if (ftest) then
         call pakndx_16(iustax,use_16)
         CLOSE (iustax)
      endif
      INQUIRE(unit=iustay,opened=ftest)
      if (ftest) then
         call pakndx_16(iustay,use_16)
         CLOSE (iustay)
      endif

    end PROGRAM WRFNC2ARL
