! Use emacs f90-mode:   -*- mode:f90  -*-

! $Id: geos2arl.f,v 1.4 2013/04/01 19:07:41 trn Exp $

!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: GEOS2ARL      DECODE netcdf GEOS MODEL FIELD FOR HYSPLIT
!
! Written by: Thomas Nehrkorn, AER, Inc. 
! Based on:
! WRFNC2ARL
! MAIN PROGRAM: RSMS2ARL      DECODE GRIB RSM MODEL FIELD FOR HYSPLIT
!   PRGMMR: DRAXLER          ORG: R/E/AR      DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE written at AER ...
!     converts a GEOS model netcdf file to arl packed format.
!     processes only one time period per execution.  
!
! PROGRAM HISTORY LOG:
! $Log: geos2arl.f,v $
! Revision 1.4  2013/04/01 19:07:41  trn
! Cosmetic changes
!
! Revision 1.4  2013/04/01 19:07:41  trn
! Cosmetic changes
!
! Revision 1.3  2013/04/01 16:45:26  trn
! Changes for first successful run through using GEOS572 data
!
! Revision 1.2  2013/03/29 17:21:40  trn
! Intermediate version, still debugging
!
! Revision 1.1  2013/03/28 19:02:23  trn
! Added initial version of geos2arl converter: compiles and links
!
!
! USAGE:  GEOS2ARL -P par_fname [-V iprint]
!                   [-M mdlid] [-T top_sigma] [-GP geos_prefix] [-GV geos_sysver] [-GS geos_suffix] valid_time
!   INPUT ARGUMENT LIST:
!     par_fname - filename of text file containing parameters (variables) info
!                 (see below for format)
!     -M mdlid   - model ID to output to cfg, data files
!     -V iprint - printout control (default:3, grid/level/summary info; less:0-2; 4: slightly more detail; >= 8: full detail )
!     -T top_sigma - include levels only as high up as top_sigma
!     -GP geos_prefix - DAS.config.mode: defaults to DAS.ops.asm
!     -GV geos_sysver - e.g., "GEOS572", "GEOS591"
!     -GS geos_suffix - file_ver.nc4: defaults to V01.nc4
!     valid_time  - time stamp for end of averaging time period (yyyymmdd_hhmn)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     unit 10 - used for par_fname
!     NETCDF_DATA - netcdf input data files use netcdf library calls
!   OUTPUT FILES:
!     unit 20 DATA_GEOS - ARL packed data output file for mass points
!     unit 30 - used for CFG_GEOS (temp file area)
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

PROGRAM GEOS2ARL

  use read_header_module
  use process_onetime_module

  implicit none       

  include 'netcdf.inc'

  LOGICAL   :: ftest, lhelp
  integer, parameter :: lchar=256
  CHARACTER (len=lchar) :: option, FNAME
  CHARACTER (len=lchar) :: par_fname, vlabel
  character (len=4) :: mdlid
  INTEGER :: handle, status

  integer :: iarg, narg, iargc
  integer :: iutext=10, iudata=20, iucfg=30
  integer :: ierr, idatlin, ndatlin, iprint, itime

  integer, parameter :: nlevtyp=2, p_sfc=1, p_lay=2
  integer :: nvars(nlevtyp), ndim1v
  character (len=lchar), allocatable :: vlabels(:,:)
  real, allocatable :: vcnvrt(:,:),vadd(:,:)
  character (len=4), allocatable :: vnames(:,:)
  character (len=lchar), allocatable :: vfnames(:,:)
  integer, allocatable :: vi_tavg(:,:)
  character (len=4) :: vname
  character (len=1) :: sta_xy, sta_z
  logical :: have_lay, use_16
  integer :: i, iltyp, iladd, i_23d, i_tavg, iladd_soil, iladd_save
  real :: cnvrt,add
  character (len=lchar) :: geos_prefix, geos_sysver, valid_time, geos_suffix
  integer :: ifname, iloop, ivar
  real :: top_sigma

  !=>check for command line arguments
  NARG=IARGC()

  lhelp = narg .eq. 0

  !     process options
  iarg = 1
  par_fname=' '
  mdlid = 'AGEO'  !'DGEO' for 16-bit files
  geos_prefix='DAS.ops.asm'
  geos_sysver='GEOS572'
  geos_suffix='V01.nc4'
  iprint=3
  top_sigma=0.010  ! 10 mb
  do while (iarg .le. narg .and. .not. lhelp)
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
     CASE ('-T','-t')
        read(fname,*,iostat=ierr) top_sigma
        if (ierr .ne. 0) then
           write (*,*) 'Error decoding top_sigma from -T option:',trim(fname)
           lhelp = .TRUE.
        endif
     CASE ('-P','-p')
        par_fname=fname
     CASE ('-M','-m')
        mdlid=fname(1:4)
     CASE ('-GP','-gp')
        geos_prefix=fname
     CASE ('-GV','-gv')
        geos_sysver=fname
     CASE ('-GS','-gs')
        geos_suffix=fname
     CASE DEFAULT
        write (*,*) 'Invalid option:',trim(option)
        lhelp=.TRUE.
     END SELECT
     iarg=iarg+1
  end do
  lhelp = lhelp .or. par_fname .eq. ' '
  use_16 = .TRUE.
  if (mdlid(1:1) .ne. 'D') use_16 = .FALSE.
  write (*,'(a)') 'Running geos2arl $Revision: 1.4 $ with the following options:'
  write (*,'(a,t40,a)') 'Parameter codes file [-P]:',trim(par_fname)
  write (*,'(a,t40,a)') 'Model ID [-M mdlid]:',trim(mdlid)
  write (*,'(a,t40,l1)') 'Using 16 bits [mdlid(1:1) .eq. "D"]:',use_16
  write (*,'(a,t40,g15.6)') 'Go up to sigma= [-T top_sigma]:',top_sigma
  write (*,'(a,t40,a)') 'File name prefix [-GP geos_prefix]:',trim(geos_prefix)
  write (*,'(a,t40,a)') 'File name suffix [-GP geos_suffix]:',trim(geos_suffix)
  write (*,'(a,t40,a)') 'GEOS version [-GV geos_sysver]:',trim(geos_sysver)
  write (*,'(a,t40,i10)') 'Verbosity (0 - 8, def: 3) [-V iprint]:',iprint

  valid_time=' '
  ifname=1
  do while (iarg .le. narg)
     CALL GETARG(iarg,fname)
     SELECT CASE (ifname)
     CASE (1)
        valid_time=trim(fname)
     CASE DEFAULT
        write (*,*) 'Too many non-option args:',ifname
        lhelp=.TRUE.
     END SELECT
     iarg = iarg + 1
     ifname = ifname + 1
  end do
  lhelp = lhelp .or. len_trim(valid_time) .lt. 13
  write (*,'(a,t40,a)') 'valid_time:', trim(valid_time)

  if (lhelp) then
     write (*,'(a)') 'Aborting because of missing or invalid arguments.  Usage: '
     write (*,'(2a)') '  GEOS2ARL -P par_fname [-M mdlid] [-V iprint]  [-T top_sigma]', &
          &' [-GP geos_prefix] [-GV geos_sysver] [-GS geos_suffix] valid_time (as yyyymmdd_hhmn)'
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
  if (iprint .ge. 8) write (*,*) 'Allocating vlabels etc with ndatlin, ndim1v=',ndatlin, ndim1v
  allocate(vlabels(ndatlin,nlevtyp), vnames(ndatlin,nlevtyp), vi_tavg (ndatlin,nlevtyp), &
       & vcnvrt(ndatlin,nlevtyp), vadd(ndatlin,nlevtyp), vfnames (ndatlin,nlevtyp), &
       & stat=ierr)
  if (ierr .ne. 0) stop 'Allocate vlabel... failure'
  nvars(:) = 0
  have_lay = .false.
  do idatlin=1,ndatlin
     read (iutext,*,iostat=ierr) vname,sta_xy,sta_z,i_23d,i_tavg,vlabel,fname,cnvrt,add
     if (ierr .ne. 0) then
        write (*,*) 'Error reading data line ',idatlin,' from par_fname=',trim(par_fname)
        stop 'Error reading datlin from par_fname'
     end if
     if (iprint .ge. 8) write (*,'(3(1x,a,1x),2i4,2(1x,a,1x),2g15.6)') &
          vname,sta_xy,sta_z,i_23d,i_tavg,trim(vlabel),trim(fname),cnvrt,add
     ! ignore sta_xy, sta_z
     if (i_23d .eq. 2) then
        iltyp = p_sfc
     elseif (i_23d .eq. 3) then
        iltyp = p_lay
        have_lay = .true.
     else
        stop 'invalid value for i_23d in par_fname'
     endif
     nvars(iltyp) = nvars(iltyp) + 1
     vnames(nvars(iltyp),iltyp) = vname
     vfnames(nvars(iltyp),iltyp) = fname
     vlabels(nvars(iltyp),iltyp) = vlabel
     vcnvrt(nvars(iltyp),iltyp) = cnvrt
     vadd(nvars(iltyp),iltyp) = add
     vi_tavg(nvars(iltyp),iltyp) = i_tavg
  end do
  close (iutext)

  if (iprint .ge. 3) then
     do iloop=1,nlevtyp
        if (iloop .eq. 1) then
           iltyp=p_sfc
        elseif (iloop .eq. 2) then
           iltyp=p_lay
        else
           write (*,*) 'Internal logic error: nlevtyp>2'
           cycle
        end if
        write (*,'(a,2i4,/,3(1x,a4),2(1x,a20),2A15,/,(i5,a5,i5,2(1x,a20),2g15.6))') &
             'iltyp, nvars=',iltyp,nvars(iltyp), &
             'ivar','name','tavg','label','fname','cnvrt','add', &
             (ivar,vnames(ivar,iltyp),vi_tavg(ivar,iltyp),&
              trim(vlabels(ivar,iltyp)),trim(vfnames(ivar,iltyp)), &
              vcnvrt(ivar,iltyp),vadd(ivar,iltyp), ivar=1,nvars(iltyp))
     end do
  end if
  
  if (.not. have_lay) stop 'need to specify non-sfc variables'

  ! Read in and process GEOS files
  call process_onetime(iudata, iucfg, iprint, &
       nlevtyp, p_sfc, p_lay, top_sigma, &
       nvars, ndim1v, &
       vlabels, lchar, vnames, vcnvrt, vadd, vfnames, vi_tavg, &
       mdlid, use_16, geos_prefix, geos_sysver, valid_time, geos_suffix)

  !     close out time period and write index record(s)
  INQUIRE(unit=iudata,opened=ftest)
  if (ftest) then
     call pakndx_16(iudata,use_16)
     CLOSE (iudata)
  endif

end PROGRAM GEOS2ARL
