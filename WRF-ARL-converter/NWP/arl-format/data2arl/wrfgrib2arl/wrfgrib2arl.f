! Use emacs f90-mode:   -*- mode:f90  -*-
!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: WRFGRIB2ARL      DECODE GRIB WRF MODEL FIELD FOR HYSPLIT
!
! Written by: Thomas Nehrkorn, AER, Inc. 
! Based on:
! MAIN PROGRAM: RSMS2ARL      DECODE GRIB RSM MODEL FIELD FOR HYSPLIT
!   PRGMMR: DRAXLER          ORG: R/E/AR      DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE written at AER ...
!     converts a WRF model grib file to arl packed format.
!     processes only one time period per execution.  Multiple input files
!     may be specified on the command line.  Input data are
!     assumed to be already on a conformal map projection on model levels.
!     Only packing to ARL format and units conversion is required.
!
! PROGRAM HISTORY LOG:
! $Log: wrfgrib2arl.f,v $
! Revision 1.16  2016/04/18 14:48:21  trn
! fixed typo introduced in last change
!
! Revision 1.16  2016/04/18 14:48:21  trn
! fixed typo introduced in last change
!
! Revision 1.15  2016/04/18 14:26:14  trn
! Fixed pakset call diagnostics: nzp is never used
!
! Revision 1.14  2016/04/18 12:47:14  trn
! Fix call to pakset_16: nxp,nyp,nzp are output
!
! Revision 1.13  2016/04/12 14:54:44  trn
! Bug fix for new pakset_16: provide explicit interface block for subroutine with optional argument
!
! Revision 1.12  2008/01/17 18:41:21  trn
! Support use of WPS geolocation routines
!
! Revision 1.11  2008/01/15 19:58:27  trn
! Bug fix
!
! Revision 1.10  2008/01/15 19:56:46  trn
! Add approximate support for secant cone projection
!
! Revision 1.9  2007/10/08 15:25:55  trn
! Do not use recl from inquire function for file open
!
! Revision 1.8  2006/10/06 14:58:20  trn
! Added support for 16-bit files (off by default)
!
! Revision 1.7  2006/05/05 15:13:06  trn
! Fixed bug in call to datini
!
! Revision 1.6  2006/04/28 18:57:56  trn
! Added code for Lambert conformal projection
!
! Revision 1.5  2005/09/11 19:20:48  trn
! Increased buffer size.  Changed padding of extra row/column
!
! Revision 1.4  2005/07/25 14:08:25  trn
! Added -M option, fixed some comments
!
! Revision 1.3  2005/07/19 15:32:44  trn
! Fixed revision printout
!
! Revision 1.2  2005/07/19 15:29:35  trn
! Added options for combining grids and/or levels, fixed output longitude
!
! Revision 1.1  2005/07/05 19:26:03  trn
! Initial revision of WRF GRIB to ARL converter
!
!
! USAGE:  WRFGRIB2ARL -P par_fname -H half_fname -F full_fname [-V iprint] [-C G] [-S G] [-C L] [-S L] 
!                   [-M mdlid] [-G geoloc_type] GRIB_FILE_NAME [GRIB_FILE_NAME2] [...]
!   INPUT ARGUMENT LIST:
!     par_fname - filename of text file containing parameters (variables) info
!                 (see below for format)
!     half_fname - filename of tex file containing sigma half-layer values
!                 (see below for format)
!     full_fname - filename of tex file containing sigma full-level values
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
!     -G geoloc_type - type of mapping routines to use (one of: cmap, wps)
!     GRIB_FILE_NAME - multiple grib files may be defined (only one at a time is supported for now)
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     unit 10 - used for par_fname, half_fname, full_fname
!     GRIB_DATA - grib input data files use special directIO routines
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

      PROGRAM WRFGRIB2ARL

      implicit none       

      LOGICAL   :: ftest, lhelp
      integer, parameter :: lchar=256
      CHARACTER (len=lchar) :: option, FNAME
      CHARACTER (len=lchar) :: par_fname, half_fname, full_fname
      character (len=4) :: mdlid, geoloc_type
      INTEGER*4 :: handle, fcopen

      integer :: iarg, narg, iargc
      integer :: iutext=10, iumass=20, iustax=21, iustay=22, iucfg=30, iutime=40
      integer :: ierr, idatlin, ndatlin, iprint

      integer, parameter :: nlevtyp=3, p_sfc=1, p_half=2, p_full=3
      integer, parameter :: ngridtyp=3, p_mass=1, p_stax=2, p_stay=3
      integer :: nvars(nlevtyp,ngridtyp), nsigl(ngridtyp), ndim1l, ndim1v
      integer, allocatable :: sigl(:,:)
      integer, allocatable :: vkpds5(:,:,:), vkpds6(:,:,:), vkpds7(:,:,:)
      real, allocatable :: vcnvrt(:,:,:)
      character (len=4), allocatable :: vnames(:,:,:)
      character (len=4) :: vname
      character (len=1) :: sta_xy, sta_z
      logical :: have_half, have_full, combine_levs, combine_grids, use_16
      integer :: i, iltyp, igtyp, iladd, i_23d, kpds5, kpds6, kpds7
      real :: xsigl, sig_factor, cnvrt

!=>check for command line arguments
      NARG=IARGC()

      lhelp = narg .eq. 0

!     process options
      iarg = 0
      par_fname=' '
      full_fname=' '
      half_fname=' '
      mdlid = 'AWRF'  !'DWRF' for 16-bit files
      geoloc_type = 'WPS' !Use WPS geolocation routines by default, 'CMAP' for old hymodelc built-ins
      combine_levs=.TRUE.
      combine_grids=.TRUE.
      iprint=2
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
         CASE ('-H','-h')
            half_fname=fname
         CASE ('-F','-f')
            full_fname=fname
         CASE ('-M','-m')
            mdlid=fname(1:4)
         CASE ('-G','-g')
            geoloc_type=fname
         CASE DEFAULT
            write (*,*) 'Invalid option:',trim(option)
            lhelp=.TRUE.
         END SELECT
      end do
      lhelp = lhelp .or. par_fname .eq. ' ' .or. &
           & (full_fname .eq. ' ' .and. half_fname .eq. ' ')
      lhelp = lhelp .or. iarg .gt. narg
      use_16 = .TRUE.
      if (mdlid(1:1) .ne. 'D') use_16 = .FALSE.
      write (*,'(a)') 'Running wrfgrib2arl $Revision: 1.16 $ with the following options:'
      write (*,'(a,t30,a)') 'Parameter codes file [-P]:',trim(par_fname)
      write (*,'(a,t30,a)') 'Half-levels file [-H]:',trim(half_fname)
      write (*,'(a,t30,a)') 'Full-levels file [-F]:',trim(full_fname)
      write (*,'(a,t30,l1)') 'Combining grids [-C/-S G]:',combine_grids
      write (*,'(a,t30,l1)') 'Combining levels [-C/-S L]:',combine_levs
      write (*,'(a,t30,a)') 'Model ID [-M mdlid]:',trim(mdlid)
      write (*,'(a,t30,a)') 'Geolocation type [-G type]:',trim(geoloc_type)
      write (*,'(a,t30,a)') '   - can be one of: cmap, wps'
      write (*,'(a,t40,l1)') 'Using 16 bits [mdlid(1:1) .eq. "D"]:',use_16
      write (*,'(a,t30,i10)') 'Printout control [-V]:',iprint

      if (lhelp) then
         write (*,'(a)') 'Usage: '
         write (*,'(2a)') '  WRFGRIB2ARL -P par_fname -H half_fname -F full_fname [-V iprint]', &
              &' [-C G] [-C L] [-S G] [-S L] [-G geoloc_type]', &
              &' GRIB_FILE_NAME [GRIB_FILE_NAME2] [...]'
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
!     integer variable code (kpds5)
!     integer level type (kpds6)
!     integer level value (kpds7; ignored for 3d variables)
!     real conversion factor to be applied to data before output
      call read_header(iutext,par_fname,ndatlin,ierr,iprint)
      if (ierr .ne. 0) stop 'Error reading header from par_fname'
      ndim1v=ndatlin
      allocate(vkpds5(ndatlin,nlevtyp,ngridtyp),vkpds6(ndatlin,nlevtyp,ngridtyp), &
           & vkpds7(ndatlin,nlevtyp,ngridtyp), vnames(ndatlin,nlevtyp,ngridtyp), &
           & vcnvrt(ndatlin,nlevtyp,ngridtyp), &
           & stat=ierr)
      if (ierr .ne. 0) stop 'Allocate vkpds5... failure'
      nvars(:,:) = 0
      have_half = .false.
      have_full = .false.
      do idatlin=1,ndatlin
         read (iutext,*,iostat=ierr) vname,sta_xy,sta_z,i_23d,kpds5,kpds6,kpds7,cnvrt
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
            kpds7 = -99  !value used to enforce match for all values of level
         else
            stop 'invalid value for i_23d in par_fname'
         endif
         nvars(iltyp,igtyp) = nvars(iltyp,igtyp) + 1
         vnames(nvars(iltyp,igtyp),iltyp,igtyp) = vname
         vkpds5(nvars(iltyp,igtyp),iltyp,igtyp) = kpds5
         vkpds6(nvars(iltyp,igtyp),iltyp,igtyp) = kpds6
         vkpds7(nvars(iltyp,igtyp),iltyp,igtyp) = kpds7
         vcnvrt(nvars(iltyp,igtyp),iltyp,igtyp) = cnvrt
      end do
      close (iutext)

! Read in level info from text file(s)
! Half-layer sigma values (sta_z=' ') and full-level sigma values (sta_z='z')
! are read in from separate files
! Levels should be specified bottom-up, and include the surface (10000) for full
! levels. One or both files can be specified
! Format:
!   Header (same as for variable file):
!         1st line: contains integer (list-directed I/O) ndatlin: 
!                   number of data lines to follow the header
!         subsequent lines: any number of header lines (ignored), followed
!         by line with the string ENDHEADER in columns 1-9
!   Data lines: (ndatlin lines, each containing in list-directed I/O format)
!     real (or integer) sigma value - if the lowest level value is less than 100
!        fractional sigma values [0-1] are assumed, and converted to [0-10000]
      nsigl(:) = 0
      nsigl(p_sfc) = 1
      ndim1l=0
      do i=1,2
         sig_factor=1.
         if (i .eq. 1) then
            iltyp=p_full
            fname=full_fname
            iladd=0
         else
            iltyp=p_half
            fname=half_fname
            iladd=1
         endif
         if (fname .eq. ' ') cycle
         call read_header(iutext,fname,ndatlin,ierr,iprint)
         if (ierr .ne. 0) stop 'Error reading header from level file'
         if (ndim1l .eq. 0) then
            ndim1l=ndatlin+iladd
            allocate(sigl(ndatlin+iladd,nlevtyp),stat=ierr)
            if (ierr .ne. 0) stop 'Allocate sigl...'
            sigl(1,p_sfc)=10000
         endif
         nsigl(iltyp) = ndatlin+iladd
         sigl(1,iltyp)=10000
         do idatlin=1,ndatlin
            read (iutext,*,iostat=ierr) xsigl
            if (ierr .ne. 0) stop 'Error reading datlin from level file'
            if (idatlin .eq. 1) then
               if (xsigl .lt. 100) then
                  sig_factor=10000.
               else
                  sig_factor=1.
               endif
            endif
            sigl(idatlin+iladd,iltyp)=nint(sig_factor*xsigl)
         end do
         close(iutext)
      end do
      
! check that full/half levels are provided if needed, and that they are consistent
      if (.not. (have_half .or. have_full)) stop 'need to specify non-sfc variables'
      if (have_half .and. nsigl(p_half) .le. 1) &
           & stop 'need to provide half levels when specifying half-level variables'
      if (nsigl(p_full) .le. 1 .and. (have_full .or. combine_levs)) &
           & stop &
           & 'need to provide full levels when specifying full-level variables or combining levels'
      if (have_half .and. have_full .and. nsigl(p_half) .ne. nsigl(p_full)) &
           & stop 'inconsistent half and full levels specified'

      if (iprint .ge. 3) then
         if (nsigl(p_half) .gt. 1) write (*,'(a/,(2i10))') 'Specified half-levels:', &
              & (i,sigl(i,p_half),i=1,nsigl(p_half))
         if (nsigl(p_full) .gt. 1) write (*,'(a/,(2i10))') 'Specified full-levels:', &
              & (i,sigl(i,p_full),i=1,nsigl(p_full))
         do igtyp=1,ngridtyp
            do iltyp=1,nlevtyp
               write (*,'(a,3i10)') 'Grid type, level type, nvars=', &
                    & igtyp,iltyp,nvars(iltyp,igtyp)
               if (nvars(iltyp,igtyp) .gt. 0) write (*,'(a4,1x,a4,1x,3a10)') &
                    & 'ivar','name',' kpds5 ',' kpds6 ',' kpds7 '
               do i=1,nvars(iltyp,igtyp)
                  write (*,'(i4,1x,a4,1x,3i10)') i,vnames(i,iltyp,igtyp), &
                       & vkpds5(i,iltyp,igtyp),vkpds6(i,iltyp,igtyp),vkpds7(i,iltyp,igtyp)
               end do
            end do
         end do
      end if
      
! Read in and process GRIB filenames (at present, xtract can only handle one file)
      do while (iarg .le. narg)
         CALL GETARG(iarg,fname)
         INQUIRE(file=fname,exist=ftest)
         if(ftest)then
            HANDLE=FCOPEN(fname,'r')
            write(*,'(2a)') 'Started processing: ',trim(fname)
            call xtract(handle, iumass, iustax, iustay, iucfg, iutime, iprint, &
                 & nlevtyp, p_sfc, p_half, p_full, ngridtyp, p_mass, p_stax, p_stay, &
                 & nvars, nsigl, ndim1l, ndim1v, &
                 & sigl, vkpds5, vkpds6, vkpds7, vnames, vcnvrt, &
                 & combine_grids, combine_levs, mdlid, use_16, geoloc_type )
            CALL FCCLOS(handle,*900)
  900       continue
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

      end PROGRAM WRFGRIB2ARL
    
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EXTRACT          MAIN ROUTINE FOR DECODING WRF GRIB DATA
! Written by Thomas Nehrkorn, AER, Inc.
!   Based on:
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  EXTRACTS all RECORDs FROM INPUT GRIB FILE, finds those that match
!            the specified variables, levels, and writes them
!            TO THE OUTPUT FILE
!
! USAGE:  see below
!   INPUT FILES:
!     HANDLE defines input data file
!   OUTPUT FILES:
!     unit 20(iumass) DATA_MASS.WRF - ARL packed data output file for mass points
!     unit 21(iustax) DATA_STAX.WRF - ARL packed data output file for x-staggered points
!     unit 22(iustay) DATA_STAY.WRF - ARL packed data output file for y-staggered points
!     unit 30(iucfg) - used for CFG_MASS.WRF, CFG_STAX.WRF, CFG_STAY.WRF (temp file area)
!     unit 40(iutime) WRFTIME - time indicator file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

      SUBROUTINE xtract(handle, iumass, iustax, iustay, iucfg, iutime, iprint, &
                 & nlevtyp, p_sfc, p_half, p_full, ngridtyp, p_mass, p_stax, p_stay, &
                 & nvars, nsigl, ndim1l, ndim1v, &
                 & sigl, vkpds5, vkpds6, vkpds7, vnames, vcnvrt, &
                 & combine_grids, combine_levs, mdlid, use_16, geoloc_type )

      use delta_time_m

      implicit none 

      include 'pakset_16.INC'

!   INPUT ARGUMENT LIST:
      INTEGER*4, intent(in) :: handle
!     HANDLE - defines input file to directIO routines
      integer, intent(in) :: iumass, iustax, iustay, iucfg, iutime, iprint !units, printout ctl
      integer, intent(in) :: nlevtyp, p_sfc, p_half, p_full !level type dimension and values
      integer, intent(in) :: ngridtyp, p_mass, p_stax, p_stay !grid type dimension and values
      integer, intent(in) :: nvars(nlevtyp,ngridtyp), nsigl(ngridtyp)
!     nvars - active number of array elements of vnames, by level and grid type
!     nsigl - active number of levels, by level type
      integer, intent(in) :: ndim1l, ndim1v !leading dimensions (max values) for levels and variables
      integer, intent(in) :: sigl(ndim1l,nlevtyp) !level values (including sfc), by level type
      integer, intent(in) :: vkpds5(ndim1v,nlevtyp,ngridtyp) !integer variable code (kpds5)
      integer, intent(in) :: vkpds6(ndim1v,nlevtyp,ngridtyp) !integer level type (kpds6)
      integer, intent(in) :: vkpds7(ndim1v,nlevtyp,ngridtyp) !integer level value (kpds7; ignored for 3d variables)
      real, intent(in) :: vcnvrt(ndim1v,nlevtyp,ngridtyp) !conversion factor
      character (len=4), intent(in) :: vnames(ndim1v,nlevtyp,ngridtyp) !variable names
      logical, intent(in) :: combine_grids, combine_levs, use_16
      integer :: nxy_factor, nxy_arg
      character (len=4), intent(in) :: mdlid, geoloc_type

!     number of levels, variables, and array limits
      integer :: MAXX, MAXY, MAXB
      integer :: MAXC

!     arrays to hold grib, character, and level variable information
      INTEGER     VGRIB
      CHARACTER*4 VCHAR
!      REAL        CNVRT, CNVRT0(MVAR), CNVRT1(MVAR)  !conversions not supported

!     define arrays for product definition section, grid description,
!     and other grib information
      INTEGER   KPDS(25),KGDS(25),KPTR(25)

!     time decoding
      character (len=14) :: init_dattim,current_dattim
      integer :: iyyyy, lead_time, lead_units, iss, ierr_tim, any_ierr_tim
      

!     input data buffer and bit map section
      CHARACTER (len=1), allocatable ::  buff(:)
      LOGICAL*1, allocatable :: KBMS(:)
      LOGICAL :: FTEST

!     remap array from one to two dimensions
      integer :: nxp(ngridtyp), nyp(ngridtyp), nxy(ngridtyp)
      integer :: nxp_pakset, nyp_pakset, nzp_pakset
      REAL, allocatable :: SVAR(:)
!     unpacked output array
      REAL, allocatable :: rvar_mass(:,:),rvar_stax(:,:),rvar_stay(:,:)
      REAL, allocatable :: rvar_out(:,:)
      INTEGER :: shape(2,ngridtyp),nvar(2*ndim1l,ngridtyp),nlvl(ngridtyp), &
           & lvl(2*ndim1l,ngridtyp)

!     packed output array
      CHARACTER (len=1), allocatable ::  CVAR(:)

!     default information (grid id, output records)
      character (len=80) :: fname
      character (len=4) :: fname_grid
      integer :: iudata
      integer :: ig(ngridtyp)
      integer :: krec(ngridtyp), krec_out(ngridtyp)
      integer, save :: icall=0

      logical :: have_index(ngridtyp), have_init(ngridtyp)
      integer :: ierr
      integer :: kbyte
      integer :: klen
      integer :: koff
      integer :: kvarb
      integer :: ktype
      integer :: level
      logical :: lmatch
      integer :: igtyp, igtyp_out
      integer :: iltyp
      integer :: ivar
      integer :: ilev
      integer :: iltype
      integer :: iladd
      integer :: kret
      integer :: iyr
      integer :: imo
      integer :: ida
      integer :: ihr
      integer :: imn
      integer :: ifh
      integer :: kl

!-------------------------------------------------------------------------------
! only required when dealing with some F90 compiler
! replace ICHAR below with internally defined JCHAR function
  CHARACTER(1)                 :: mychr    
  INTEGER                      :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
!-------------------------------------------------------------------------------

  INTERFACE
  SUBROUTINE W3FI63(BUFF,KPDS,KGDS,KBMS,RVAR,KPTR,KRET)
  CHARACTER(1), INTENT(IN)  :: BUFF(:)
  LOGICAL(1),   INTENT(OUT) :: KBMS(:)
  INTEGER,      INTENT(OUT) :: KPDS(:)
  INTEGER,      INTENT(OUT) :: KGDS(:)
  REAL,         INTENT(OUT) :: RVAR(:)
  INTEGER,      INTENT(OUT) :: KPTR(:)
  INTEGER,      INTENT(OUT) :: KRET
  END SUBROUTINE w3fi63
  END INTERFACE
!-------------------------------------------------------------------

      if (icall .eq. 0) then
         ig(:)=99
         krec(:)=0
         krec_out(:)=0
         have_index(:) = .FALSE.
         have_init(:) = .FALSE.
      else
         stop 'multiple gribfiles not supported, call it one at a time'
      endif
      icall = icall + 1

!     Initial array sizes
      nxy_factor = 1
      if (use_16) nxy_factor = 2
      MAXX=300
      MAXY=300
      MAXB=150000
      MAXC=nxy_factor * 90000

      deallocate (buff, kbms, svar, cvar, stat=ierr)
      allocate (buff(maxb), kbms(maxb), cvar(maxc), stat=ierr)
      if (ierr .ne. 0) stop 'Initial allocate buff...'
      allocate (svar(maxx*maxy), stat=ierr)
      if (ierr .ne. 0) stop 'Initial allocate svar...'

!     input grib record byte counter
      kbyte=0

!     read the indicator section
  100 CONTINUE
      call fcptps(handle,kbyte,*900)
      call fcread(handle,buff,1,8,*900)
      IF(BUFF(1)//BUFF(2)//BUFF(3)//BUFF(4).NE.'GRIB')THEN
         KBYTE=KBYTE+1
         GOTO 100
      END IF

!     determine the length of the entire grib record
      klen=jchar(buff(7))+jchar(buff(6))*256+jchar(buff(5))*65536
      if(klen.gt.maxb)then
         write(*,*)'Grib record: ',klen
         write(*,*)'Exceedes buffer: ',maxb,' , reallocating'
         maxb=klen+1000
         maxc=nxy_factor * maxb
         deallocate (buff, kbms, cvar, stat=ierr)
         allocate (buff(maxb), kbms(maxb), cvar(maxc), stat=ierr)
         if (ierr .ne. 0) stop 'Reallocate buff...'
         goto 100
      end if

!     load the indicator section and pds segment into the buffer
      call fcptps(handle,kbyte,*900)
      call fcread(handle,buff,1,36,*900)

!     product definition section (+8 byte indicator section offset)
      koff=8
      kvarb=jchar(buff(koff+9))
      ktype=jchar(buff(koff+10))
      level=jchar(buff(koff+12))+jchar(buff(koff+11))*256

! look for match: determine igtyp, iltyp, ivar, ilev
      lmatch = .FALSE.
!     check if variable present in selection table
      MATCH_LOOP: do igtyp=1,ngridtyp
         do iltyp=1,nlevtyp
            do ivar=1,nvars(iltyp,igtyp)
               lmatch=kvarb .eq. vkpds5(ivar,iltyp,igtyp) .and. &
                    & ktype .eq. vkpds6(ivar,iltyp,igtyp) .and. &
                    & (level .eq. vkpds7(ivar,iltyp,igtyp) .or. &
                    &  vkpds7(ivar,iltyp,igtyp) .eq. -99)
               if (lmatch) exit MATCH_LOOP
            end do
         end do
      end do MATCH_LOOP

!     check if 3d level is present in selection table
      ilev = 1
      if (lmatch .and. iltyp .ne. p_sfc) then
         lmatch = .FALSE.
         if (iltype .eq. p_full) iladd = 0
         if (iltype .eq. p_half) iladd = 1
         do ilev=1+iladd,nsigl(iltyp)
            lmatch = level .eq. sigl(ilev,iltyp)
            if (lmatch) exit
         end do
      end if

      if (iprint .ge. 9) write (*,'(a,l10,10(a,i10))') 'lmatch=',lmatch,' kbyte=',kbyte , &
           & ' kvarb=',kvarb,' ktype=',ktype, ' level=',level, &
           & ' igtyp=',igtyp,' iltyp=',iltyp,' ivar=',ivar,' ilev=',ilev
!     if all tests fail go and read next grib record
      if (.not. lmatch) go to 800

!     skip until the mass grid info is available for makndx
      if (combine_grids .and. igtyp .ne. p_mass .and. .not. have_index(p_mass)) goto 800

!     load the entire grib data record into the buffer
      krec(igtyp) = krec(igtyp) + 1
      igtyp_out = igtyp
      if (combine_grids) igtyp_out = p_mass
      krec_out(igtyp_out) = krec_out(igtyp_out) + 1
      if (iprint .ge. 8) write (*,'(10(a,i10))') ' kbyte=',kbyte , &
           & ' kvarb=',kvarb,' ktype=',ktype, ' level=',level, &
           & ' igtyp=',igtyp,' krec(igtyp)=',krec(igtyp), &
           & ' igtyp_out=',igtyp_out,' krec_out(igtyp_out)=',krec_out(igtyp_out), &
           & ' iltyp=',iltyp,' ivar=',ivar,' ilev=',ilev

      call fcptps(handle,kbyte,*900)
      call fcread(handle,buff,1,klen,*900)

  300 continue

!     call the nmc grib unpacker
      call W3FI63(buff,kpds,kgds,kbms,svar,kptr,kret)
      if(kret.ne.0)then
         write(*,*)'Error W3FI63: ',kret
         stop
      end if

!#    generate valid date/time
      any_ierr_tim = 0
      iyyyy = kpds(8)
      imo = kpds(9)
      ida = kpds(10)
      ihr = kpds(11)
      imn = kpds(12)
      ifh = 0 !treat all data as initial time
      iss = 0
      if (iyyyy .lt. 1000) iyyyy = iyyyy + 100*(kpds(21)-1)
      iyr = mod(iyyyy,100)
      write (init_dattim,'(i4.4,5i2.2)',iostat=ierr_tim) &
           & iyyyy, imo, ida, ihr, imn, iss
      if (ierr_tim .ne. 0) any_ierr_tim = any_ierr_tim + 1
      if (any_ierr_tim .ne. 0 .or. iprint .ge. 8 .or. &
           & (krec_out(igtyp_out) .eq. 1 .and. iprint .ge. 2)) then
         write (*,'(11(a,i10))') ' iyyyy=',iyyyy,' iyr=',iyr,' kpds(8)=',kpds(8), &
              & ' imo=',imo,' kpds(9)=',kpds(9),' ida=',ida,' kpds(10)=',kpds(10), &
              & ' ihr=',ihr,' kpds(11)=',kpds(11),' imn=',imn,' kpds(12)=',kpds(12)
         write (*,'(2a)') 'init_dattim= ',init_dattim
         if (any_ierr_tim .ne. 0) stop 'Error decoding time'
      endif
      lead_units=kpds(13)
      if (lead_units .eq. 254) lead_units=-1 !delta_time uses -1 for seconds
      if (kpds(16) .le. 1) then
         lead_time=kpds(14) !valid at time range 1
      elseif (kpds(16) .eq. 10) then
         lead_time=max(kpds(14), kpds(15))  !valid at time range 1 or 2 ??
      else
         lead_time=kpds(15) !valid at time range 2
      endif

      if (lead_time .eq. 0) then	!valid at initial time
         current_dattim = init_dattim
      else	!compute valid time
         call delta_time(init_dattim, lead_units, lead_time, &
              & current_dattim, ierr_tim)
         if (ierr_tim .ne. 0) any_ierr_tim = any_ierr_tim + 1
      endif
      read (current_dattim, '(i4.4,5i2.2)',iostat=ierr_tim) iyyyy,imo,ida,ihr,imn,iss
      if (ierr_tim .ne. 0) any_ierr_tim = any_ierr_tim + 1
      iyr=mod(iyyyy,100)
      if (any_ierr_tim .ne. 0 .or. iprint .ge. 8 .or. &
           & (krec_out(igtyp_out) .eq. 1 .and. iprint .ge. 2)) then
         write (*,'(10(a,i10))') ' lead_units',lead_units,' kpds(13)=',kpds(13), &
              & ' lead_time=',lead_time,' kpds(14)=',kpds(14),' kpds(15)=',kpds(15), &
              & ' kpds(16)=',kpds(16)
         write (*,'(2a)') 'current_dattim= ',current_dattim
         if (any_ierr_tim .ne. 0) stop 'Error decoding time'
      endif

      if (igtyp_out .eq. p_mass) then
         iudata=iumass
      elseif (igtyp_out .eq. p_stax) then
         iudata=iustax
      elseif (igtyp_out .eq. p_stay) then
         iudata=iustay
      endif

!     for the first record create an index record for pakrec and perform initializations
      if (krec(igtyp) .eq. 1) then

!        decode grib information and create index record structure
         call makndx(kgds, igtyp, nlevtyp, p_sfc, p_half, p_full, &
              & ngridtyp, p_mass, p_stax, p_stay, &
              & vnames, nvars, nsigl, ndim1l, ndim1v, &
              & sigl, iprint, iucfg, ig, nxp, nyp, &
              & nvar, lvl, nlvl, combine_grids, combine_levs, have_index(igtyp_out),mdlid, geoloc_type)
         have_index(igtyp_out) = .TRUE.
            
         if(nxp(igtyp) .gt. maxx .or. nyp(igtyp) .gt. maxy) then
            write(*,*)'Real array size: ',nxp(igtyp),nyp(igtyp)
            write(*,*)'Exceeds dimensions: ',maxx,maxy,' ,reallocating'
            deallocate(svar,stat=ierr)
            maxx=nxp(igtyp) + 1
            maxy=nyp(igtyp) + 1
            allocate(svar(maxx*maxy),stat=ierr)
            if (ierr .ne. 0) stop 'Reallocating svar...'
            goto 300
         end if
         if (.not. have_init(igtyp)) then
            shape(1,igtyp)=nxp(igtyp)
            shape(2,igtyp)=nyp(igtyp)
            nxy(igtyp) = nxp(igtyp)*nyp(igtyp)
            if (combine_grids .and. .not. have_init(igtyp_out)) then
               if (igtyp .ne. igtyp_out) stop 'xtract: igtyp internal logic error'
               nxy(igtyp) = (nxp(igtyp)+1)*(nyp(igtyp)+1)
               allocate(rvar_out(nxp(igtyp)+1,nyp(igtyp)+1),stat=ierr)
               if (ierr .ne. 0) stop 'Allocate rvar_out...'
            endif
            if(nxy_factor * nxy(igtyp) .gt. maxc) then
               write(*,*)'Packed output array: ',nxy(igtyp)
               write(*,*)'Exceeds dimensions: ',maxc,' ,reallocating'
               deallocate(cvar,stat=ierr)
               maxc=(nxp(igtyp)+1)*(nyp(igtyp)+1) + 1000
               maxc = nxy_factor * maxc
               allocate(cvar(maxc),stat=ierr)
               if (ierr .ne. 0) stop 'Reallocating cvar...'
            end if

            if (combine_grids) FNAME_grid='MASS'
            if (igtyp .eq. p_mass) then
               if (.not. combine_grids) FNAME_grid='MASS'
               allocate(rvar_mass(nxp(igtyp),nyp(igtyp)),stat=ierr)
               if (ierr .ne. 0) stop 'Allocate rvar_mass...'
            elseif (igtyp .eq. p_stax) then
               if (.not. combine_grids) FNAME_grid='STAX'
               allocate(rvar_stax(nxp(igtyp),nyp(igtyp)),stat=ierr)
               if (ierr .ne. 0) stop 'Allocate rvar_stax...'
            elseif (igtyp .eq. p_stay) then
               if (.not. combine_grids) FNAME_grid='STAY'
               allocate(rvar_stay(nxp(igtyp),nyp(igtyp)),stat=ierr)
               if (ierr .ne. 0) stop 'Allocate rvar_stay...'
            endif

            if (.not. have_init(igtyp_out)) then
               fname = 'CFG_' // fname_grid // '.WRF'
               if (iprint .ge. 2) write (*,'(a,i10,1x,a)') &
                    & 'Calling pakset with unit, file=',iudata,trim(fname)
               call pakset_16(iudata,FNAME,1,nxp_pakset,nyp_pakset,nzp_pakset,use_16)
               if (combine_grids) then
                  if (nxp_pakset .ne. nxp(igtyp_out)+1 .or. nyp_pakset .ne. nyp(igtyp_out)+1) &
                     write (*,*) 'Inconsistent nxp/nyp from pakset:',&
                          nxp_pakset, nxp(igtyp_out)+1, nyp_pakset, nyp(igtyp_out)+1
               else
                  if (nxp_pakset .ne. nxp(igtyp_out) .or. nyp_pakset .ne. nyp(igtyp_out)) &
                     write (*,*) 'Inconsistent nxp/nyp from pakset:',&
                          nxp_pakset, nxp(igtyp_out), nyp_pakset, nyp(igtyp_out)
               endif
               !        standard file name for output
               fname = 'DATA_' // fname_grid // '.WRF'

!        open output data set and initialize to missing
               if (iprint .ge. 2) write (*,'(a,i10,1x,a,4i10)') &
                    & 'Calling datini with unit, file, iyr,imo,ida,ihr=', iudata,trim(fname), &
                    & IYR,IMO,IDA,IHR
               nxy_arg=nxy_factor * NXY(igtyp)
               call datini(iudata,fname,CVAR,nxy_arg,NLVL(igtyp),NVAR(1,igtyp),IG(igtyp), &
                    & IYR,IMO,IDA,IHR)
               if (iprint .ge. 1) write(*,'(2a)') 'Initialized output data set: ',trim(fname)
               if (combine_grids) then
!        for combine_grids run, start over from beginning of file after initializations
                  kbyte=0
                  krec(:)=0
                  krec_out(:)=0
                  have_init(igtyp_out) = .TRUE.
                  goto 100
               endif
            endif
            have_init(igtyp) = .TRUE.
         endif !have_init(igtyp)
      END IF !krec(igtyp)
      
! determine ARL level number
      if (iltyp .eq. p_sfc) then
         kl=1
      else
         if (.not. combine_levs .or. iltyp .eq. p_full) then
            do kl=1,nlvl(igtyp_out)
               if (level .eq. lvl(kl,igtyp_out)) exit
            enddo
         else
            kl = ilev !half-levels are assigned to next higher full level if combine_levs
         endif
      endif

!     remap input data from one- to two-dimensional array
!     perform any required units conversion
!     [if combining grids, copy into output array]
!     then pack into ARL format and continue
      if (iprint .ge. 8) write (*,'(a,i10,8i10)') &
           & 'Calling pakrec with unit, iyr,imo,ida,ihr,imn,ifh,kl,0=', iudata, &
           & IYR,IMO,IDA,IHR,IMN,IFH,KL,0
      if (igtyp .eq. p_mass) then
         RVAR_mass=RESHAPE(SVAR(1:nxy(igtyp)),SHAPE(1:2,igtyp))
         if(vcnvrt(ivar,iltyp,igtyp) .ne. 1.0) &
              & call datcnv(rvar_mass,nxp(igtyp),nyp(igtyp),vcnvrt(ivar,iltyp,igtyp))
         if (.not. combine_grids) then
            nxy_arg=nxy_factor * NXY(igtyp)
            call pakrec_16(iudata,rvar_mass,cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                 & vnames(ivar,iltyp,igtyp), &
                 & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
         else
            call copy_arr(rvar_mass,nxp(igtyp),nyp(igtyp), &
                 & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
         endif
      elseif (igtyp .eq. p_stax) then
         RVAR_stax=RESHAPE(SVAR(1:nxy(igtyp)),SHAPE(1:2,igtyp))
         if(vcnvrt(ivar,iltyp,igtyp) .ne. 1.0) &
              & call datcnv(rvar_stax,nxp(igtyp),nyp(igtyp),vcnvrt(ivar,iltyp,igtyp))
         if (.not. combine_grids) then
            nxy_arg=nxy_factor * NXY(igtyp)
            call pakrec_16(iudata,rvar_stax,cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                 & vnames(ivar,iltyp,igtyp), &
                 & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
         else
            call copy_arr(rvar_stax,nxp(igtyp),nyp(igtyp), &
                 & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
         endif
      elseif (igtyp .eq. p_stay) then
         RVAR_stay=RESHAPE(SVAR(1:nxy(igtyp)),SHAPE(1:2,igtyp))
         if(vcnvrt(ivar,iltyp,igtyp) .ne. 1.0) &
              & call datcnv(rvar_stay,nxp(igtyp),nyp(igtyp),vcnvrt(ivar,iltyp,igtyp))
         if (.not. combine_grids) then
            nxy_arg=nxy_factor * NXY(igtyp)
            call pakrec_16(iudata,rvar_stay,cvar,nxp(igtyp),nyp(igtyp),nxy_arg, &
                 & vnames(ivar,iltyp,igtyp), &
                 & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
         else
            call copy_arr(rvar_stay,nxp(igtyp),nyp(igtyp), &
                 & rvar_out,nxp(igtyp_out)+1,nyp(igtyp_out)+1)
         endif
      endif
      if (combine_grids) then
         nxy_arg=nxy_factor * NXY(igtyp_out)
         call pakrec_16(iudata,rvar_out,cvar,nxp(igtyp_out)+1,nyp(igtyp_out)+1,nxy_arg, &
              & vnames(ivar,iltyp,igtyp), &
              & IYR,IMO,IDA,IHR,IMN,IFH,KL,0,use_16)
      end if

  800 kbyte=kbyte+klen
      GO TO 100

  900 continue !error exit from fcpts, fcread ?
      do igtyp=1,ngridtyp
         write (*,'(a,i10,a,2i10)') 'For grid type ',igtyp, &
              & ' the number of processed records is:', krec(igtyp)
         write (*,'(a,i10,a,2i10)') 'For grid type ',igtyp, &
              & ' the number of processed output records is:', krec_out(igtyp)
         if (nlvl(igtyp) .gt. 0) then
            if (sum(nvar(1:nlvl(igtyp),igtyp)) .ne. krec_out(igtyp)) &
                 & write (*,'(a,i10)') &
                 & 'WARNING: The expected number of output records is:', &
                 & sum(nvar(1:nlvl(igtyp),igtyp)), &
                 & '      in the case of missing sfc variables:', &
                 & sum(nvar(2:nlvl(igtyp),igtyp))
         else
            if (krec_out(igtyp) .gt. 0) &
                 & write (*,'(a,i10)') 'WARNING: The expected number of output records is:', 0
         end if
      end do
      RETURN
      end SUBROUTINE xtract

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAKNDX           CREATE THE INDEX RECORD CONFIGURATION FILE
!
! Written by Thomas Nehrkorn, AER, Inc.
! Based on:
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  CREATES THE CONFIGURATION FILE FOR THE OUTPUT DATA WHICH
!            DEFINES THE GRID SYSTEM AS WELL AS ALL VARIABLES THAT ARE
!            TO BE WRITTEN TO THE OUTPUT FILE.
!
! USAGE:  see below
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT iucfg - one of CFG_(MASS|STAX|STAY).WRF, defines the output file grid and structure
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  portable
!
!$$$

!USAGE: 
      SUBROUTINE MAKNDX(kgds, igtyp, nlevtyp, p_sfc, p_half, p_full, &    !
              & ngridtyp, p_mass, p_stax, p_stay, &                       !
              & vnames, nvars, nsigl, ndim1l, ndim1v, &                   !
              & sigl, iprint, iucfg, ig, nxp, nyp, &                      !
              & nvar, lvl, nlvl, combine_grids, combine_levs, have_index, mdlid, geoloc_type)             !

      implicit none

!   INPUT ARGUMENT LIST:

      INTEGER, intent(in) :: kgds(25), igtyp, nlevtyp, p_sfc, p_half, p_full, &
           & ngridtyp, p_mass, p_stax, p_stay, ndim1l, ndim1v
!     KGDS - grid definitions from w3lib decoder
!     igtyp - grid type (one of the ngridtyp values: p_mass, p_stax, p_stay)
!     nlevtyp - number of level types: p_sfc, p_half, p_full
!     ndim1l, ndim1v - leading dimensions (max values) for levels and variables

      character (len=4), intent(in) :: vnames(ndim1v,nlevtyp,ngridtyp), mdlid, geoloc_type
!     Vnames - character array of variable names
!     mdlid - model ID for output in cfg file
!     geoloc_type - which mapping routines to use

      integer, intent(in) :: nvars(nlevtyp,ngridtyp), nsigl(nlevtyp), &
           & sigl(ndim1l,nlevtyp), iprint, iucfg, ig(ngridtyp)
!     sigl - level values (including sfc), by level type
!     iprint - printout control
!     iucfg - unit number to use for configuration file (reusable)
!     IG - grid identification number

      logical, intent(in) :: combine_grids, combine_levs, have_index
!     flags controlling whether to combine grids, levels, and whether to create index

!   OUTPUT ARGUMENT LIST (output for current igtyp, unchanged for others):
      integer, intent(inout) :: nxp(ngridtyp), nyp(ngridtyp)
!     NXP,NYP - output grid dimensions

      integer, intent(inout) :: nvar(2*ndim1l,ngridtyp), lvl(2*ndim1l,ngridtyp), nlvl(ngridtyp)
!     NLVL - number of data levels in output file, by grid type
!     nvar - number of variables in output file, by level and grid type
!     lvl - height of each output level, by grid type
!           (merges sfc, full, and half levels as needed)

!   Local data:
      LOGICAL FTEST
      REAL   grids(12), sigma, true1, true2, tang_lat
      real, external :: eqvlat
      integer :: i, ii, ii1, ii2, iladd, ilev, kl, kladd, lev2, nl, nv
      integer :: nxp_out, nyp_out, igtyp_out

!     arrays to hold variable selection information

      CHARACTER (len=4) :: VCHAR0(2*ndim1v), VCHAR1(2*ndim1v,2*ndim1l), fname_grid
!     VCHAR0 - character array of surface field identification
!     VCHAR1 - character array of upper level field identifications, by level

      character (len=80) :: fname
      integer :: iproj
      real, parameter :: vmiss=-999.99

!     the number of grid points
      igtyp_out=igtyp
      if (combine_grids) igtyp_out = p_mass

      nxp(igtyp)=kgds(2)
      nyp(igtyp)=kgds(3)

      if (have_index) return   !only decode nxp, nyp if index has already been created

!     grid orientation
      grids(6)=0.0

!     lat/lon of grid point 1,1
      grids(8)=1.0
      grids(9)=1.0
      grids(10)=kgds(4)/1000.0
      grids(11)=kgds(5)/1000.0
!     constrain grid point 1,1 longitude to [-180, +180] (needed ??)      
      if (grids(11) .gt.  180.0) grids(11) = grids(11) - 360.0
      if (grids(11) .lt. -180.0) grids(11) = grids(11) + 360.0
      grids(12)=0.0

!     defines a polar sterographic projection
      if(kgds(1).eq.5)then
         iproj = 2
!        set the pole position and reference lat/lon
         if(kgds(10).eq.0)then
            grids(1)=90.0
            grids(3)=60.0
         else
            grids(1)=-90.0
            grids(3)=-60.0
         end if
         grids(2)=0.0

!        reference longitude = grid orientation
         grids(4)=kgds(7)/1000.0

!        delta=x grid size
         grids(5)=kgds(8)/1000.0

!        tangent latitude
         grids(7)=90.0

!     defines a mercator projection
      elseif(kgds(1).eq.1)then

         iproj = 3

!        pole lat/lon axis through pole
         grids(1)=90.0
         grids(2)=0.0

!        reference lat = latitude of projection intersection
         grids(3)=kgds(9)/1000.0
!        reference lon =
         grids(4)=0.0

!        longitudinal direction grid length
         grids(5)=kgds(12)/1000.0

!        tangent latitude
         grids(7)=0.0

      elseif(kgds(1).eq.3)then

         iproj = 1

!     defines a lambert conformal projection (taken from eta40arl.f)

!        set the pole position and reference lat/lon
         grids(1)=90.0
         grids(2)=0.0

!        reference latitude
         true1 = kgds(12)/1000.0
         true2 = kgds(13)/1000.0
         tang_lat = eqvlat(true1,true2)
         grids(3)=tang_lat

!        reference longitude = grid orientation
         grids(4)=kgds(7)/1000.0

!        delta=x grid size in km
         grids(5)=kgds(8)/1000.0

!        tangent latitude
         grids(7)=tang_lat

      else
         write(*,*)'Undefined projection: ',kgds(1)
         stop
      end if

      if (trim(geoloc_type) .eq. 'wps' .or. trim(geoloc_type) .eq. 'WPS') then
! Redefine parts of the header for the WPS geolocation routines
         grids(1) = vmiss
         grids(2) = vmiss
         grids(3) = vmiss
         grids(6) = real(iproj)
         if (iproj == 1) then
            grids(7) = true1
            grids(12) = true2
         endif
      else
         if (trim(geoloc_type) .ne. 'cmap' .and. trim(geoloc_type) .ne. 'CMAP') &
              & write (*,*) 'Unsupported geoloc_type=',trim(geoloc_type),' using default=cmap'
      end if
      

!     write the packer configuration file
      
      if (igtyp_out .eq. p_mass) fname_grid='MASS'
      if (igtyp_out .eq. p_stax) fname_grid='STAX'
      if (igtyp_out .eq. p_stay) fname_grid='STAY'
      fname = 'CFG_' // fname_grid // '.WRF'

      if (iprint .ge. 2) write (*,'(a,i10,1x,a)') &
           & 'Writing cfg info to unit, file=',iucfg,trim(fname)
      OPEN(iucfg,file=fname)
      WRITE(iucfg,'(20X,A4)') mdlid

!     grid number 99 and 1 for sigma coordinate system
      WRITE(iucfg,'(20X,I4)') IG(igtyp_out), 1

      WRITE(iucfg,'(20X,F10.2)')(GRIDS(I),I=1,12)

! level and variable information
      if (combine_grids) then
         ii1=1
         ii2=ngridtyp
      else
         ii1=igtyp
         ii2=igtyp
      endif
      nvar(:,igtyp_out) = 0
! sfc:
! sfc variables, if any:
      do ii=ii1,ii2
         do i=1,nvars(p_sfc,ii)
            vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_sfc,ii)
         end do
         nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_sfc,ii)
! add lowest full level variables, if any
         if (nvars(p_full,igtyp) .gt. 0) then
            do i=1,nvars(p_full,ii)
               vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_full,ii)
            end do
            nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_full,ii)
         end if
      end do
! if neither sfc nor full variables there, use lowest half-level
! (will be missing-filled)
      if (nvar(1,igtyp_out) .eq. 0) then
         do ii=ii1,ii2
            do i=1,nvars(p_half,ii)
               vchar0(nvar(1,igtyp_out)+i) = vnames(i,p_half,ii)
            end do
            nvar(1,igtyp_out) = nvar(1,igtyp_out) + nvars(p_half,ii)
         enddo
      endif
      if (iprint .ge. 4) write (*,'(a,(10(a4,1x)))') 'vchar0=',(vchar0(i),i=1,nvar(1,igtyp_out))

! continue here

! First level is at surface:
      kl = 1
      lvl(kl,igtyp_out) = 10000
! Upper levels:
      lev2 = max(nsigl(p_half),nsigl(p_full))
      do ilev=1,lev2-1
         kladd=1
! For half-levels, loop over all grid types that are to be combined:
         do ii=ii1,ii2
            if (nvars(p_half,ii) .gt. 0) then
! Increment output level counter
               kl = kl + kladd
! Half-level variables:
               do i=1,nvars(p_half,ii)
                  vchar1(nvar(kl,igtyp_out)+i,kl) = vnames(i,p_half,ii)
               end do
               nvar(kl,igtyp_out) = nvar(kl,igtyp_out) + nvars(p_half,ii)
! Disable level incrementing for subsequent grid and level types
               kladd = 0
               if (.not. combine_levs) then
                  !           Keep half-levels as separate levels:
                  iladd=1
                  lvl(kl,igtyp_out) = sigl(ilev+iladd,p_half)
                  if (iprint .ge. 4) write (*,'(a,4i10)') 'half ilev,kl,lvl,nvar=', &
                       & ilev,kl,lvl(kl,igtyp_out),nvar(kl,igtyp_out)
                  if (lvl(kl,igtyp_out) .ge. lvl(kl-1,igtyp_out)) stop 'levels out of order'
               end if
            end if
         end do

! Enable level incrementing for full levels if not combining levels
         if (.not. combine_levs) kladd = 1
! For full-levels, loop over all grid types that are to be combined:
         do ii=ii1,ii2
            if (nvars(p_full,ii) .gt. 0) then
! Increment output level counter
               kl = kl + kladd
! Half-level variables:
               do i=1,nvars(p_full,ii)
                  vchar1(nvar(kl,igtyp_out)+i,kl) = vnames(i,p_full,ii)
               end do
               nvar(kl,igtyp_out) = nvar(kl,igtyp_out) + nvars(p_full,ii)
! Disable level incrementing for subsequent grid and level types
               kladd = 0
               iladd=1
               lvl(kl,igtyp_out) = sigl(ilev+iladd,p_full)
               if (iprint .ge. 4) &
                    & write (*,'(a,4i10)') 'full ilev,kl,lvl,nvar=', &
                    & ilev,kl,lvl(kl,igtyp_out),nvar(kl,igtyp_out)
               if (lvl(kl,igtyp_out) .ge. lvl(kl-1,igtyp_out)) stop 'levels out of order'
            end if
         end do

      end do
       
      nlvl(igtyp_out) = kl

      if (combine_grids) then
         if (igtyp .ne. p_mass) stop 'makndx: igtyp internal logic error'
         nxp_out=nxp(igtyp)+1
         nyp_out=nyp(igtyp)+1
      else
         nxp_out=nxp(igtyp)
         nyp_out=nyp(igtyp)
      endif
      WRITE(iucfg,'(20X,I4)')nxp_out,nyp_out,nlvl(igtyp_out)

      DO NL=1,NLVL(igtyp_out)
         sigma=lvl(nl,igtyp_out)/10000.0
         if(nl.eq.1)then
!     sfc level information
            WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
            sigma,nvar(nl,igtyp_out),(vchar0(nv),nv=1,nvar(nl,igtyp_out))
         else
!     upper level information
            WRITE(iucfg,'(20X,F6.4,I3,99(1X,A4))')                                &
            sigma,nvar(nl,igtyp_out),(vchar1(nv,nl),nv=1,nvar(nl,igtyp_out))
         end if
      END DO
      CLOSE (iucfg)
      RETURN
      end SUBROUTINE MAKNDX
    

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DATINI           INITIALIZE OUTPUT FILE FOR ONE TIME PERIOD
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            DUMMY ROUTINE TO FILL IN HEADER BYTES ON ALL RECORS OF THE
!            OUTPUT DATA FILE WITH CORRECT TIME BUT NULL DATA FIELDS.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 - RRD
!
! USAGE:  CALL DATINI(KUNIT,fname,CVAR,NXY,NLVL,NVAR,IG,IY,IM,ID,IH)
!   INPUT ARGUMENT LIST:
!     KUNIT - unit number of the output file
!     CVAR - dummy character string to represent the gridded data
!     NXY - the length of the gridded data field in bytes
!     NLVL - number of data levels in output file
!     NVAR - array defining the number of variables on each level
!     IG - grid identification number
!     IY,IM,ID,IH - current date/time
!   OUTPUT ARGUMENT LIST:
!     NONE
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     UNIT 20 - DATA.RSM is the output data file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE DATINI(KUNIT,fname,CVAR,NXY,NLVL,NVAR,IG,IY,IM,ID,IH)

      implicit none 

      integer ,intent(in) :: ig, iy, im, id, ih, kunit, nxy, nlvl
      INTEGER ,intent(in) ::   NVAR(NLVL)
      character (len=*), intent(in) :: fname

      CHARACTER (len=1), intent(out) :: CVAR(NXY)

      integer :: ic, il, k, lrec, lrectest, mrec, nexp, nl, nv
      real :: prec, var1
      character (len=50) :: LABEL

      inquire(iolength=lrectest) label,cvar
      lrec = nxy+50
      if (lrec .ne. lrectest) print *,&
           & 'Warning from datini: nxy, lrectest, lrec=',&
           & nxy,lrectest,lrec
      OPEN(kunit,FILE=FNAME,RECL=LREC,ACCESS='DIRECT',FORM='UNFORMATTED')

      IC=-1
      NEXP=0
      PREC=0.0
      VAR1=0.0
      MREC=0

!     initialize packed data array
      DO K=1,NXY
         CVAR(K)=' '
      END DO
      CVAR(NXY)=CHAR(13)

!     header label output format
  100 FORMAT(7I2,A4,I4,2E14.7)

!     index record
      WRITE(LABEL,100)IY,IM,ID,IH,IC,0,IG,'NULL',NEXP,PREC,VAR1
      MREC=MREC+1
      WRITE(KUNIT,REC=MREC)LABEL,CVAR

      DO NL=1,NLVL
      DO NV=1,NVAR(NL)

         IL=NL-1
         WRITE(LABEL,100)IY,IM,ID,IH,IC,IL,IG,'NULL',NEXP,PREC,VAR1
         MREC=MREC+1
         WRITE(KUNIT,REC=MREC)LABEL,CVAR

      END DO
      END DO

      RETURN
      end SUBROUTINE DATINI
    

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DATCNV           CONVERT UNITS FOR ALL ELEMENTS OF THE DATA
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            IS USED TO CONVERT EACH ELEMENT OF THE OUTPUT DATA GRID
!            TO CONFORM WITH UNIT CONVENTIONS USED IN OTHER APPLICATIONS
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 - RRD
!
! USAGE:  CALL DATCNV(RVAR,NXP,NYP,CNVRT)
!   INPUT ARGUMENT LIST:
!     RVAR - real data array
!     NXP,NYP - dimensions of the array
!     CNVRT - conversion factor
!   OUTPUT ARGUMENT LIST:
!     RVAR - real data array after conversion
!   INPUT FILES:
!     NONE
!   OUTPUT FILES:
!     NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  CRAY
!
!$$$

      SUBROUTINE datcnv(rvar,nxp,nyp,cnvrt)

      implicit none

      integer, intent(in) :: nxp, nyp
      real, intent(in) :: cnvrt
      REAL, intent(inout) ::   RVAR(NXP,NYP)

      integer :: i,j

      DO J=1,NYP
      DO I=1,NXP

         rvar(i,j)=rvar(i,j)*cnvrt

      END DO
      END DO
      RETURN
      end SUBROUTINE datcnv
    
      subroutine read_header(iu,fname,ndatlin,ierr,iprint)

      implicit none

      integer, intent(in) :: iu, iprint  !unit number, printout ctl
      character (len=*) :: fname !file name
      integer, intent(out) :: ndatlin, ierr !no of datalines, error code

      integer :: status, iline
      character(len=256) :: line

      ndatlin=0
      ierr = -1
      if (iprint .ge. 1) write (*,'(2a)') 'Reading header from file: ',trim(fname)
      open(unit=iu,file=fname,iostat=status)
      if (status .ne. 0) goto 900

      ierr = 1
      read(iu,*,iostat=status) ndatlin
      if (status .ne. 0) goto 900
      if (iprint .ge. 4) write (*,'(a,i10)') 'Decoded ndatlin= ',ndatlin

      LINE_LOOP: do iline=2,10000
         ierr = iline
         read (iu,'(a)',iostat=status) line
         if (status .ne. 0) goto 900
         if (iprint .ge. 4) write (*,'(2a)') 'Header line: ',trim(line)
         if (line(1:9) .eq. 'ENDHEADER') then
            ierr = 0
            exit LINE_LOOP
         end if
      end do LINE_LOOP
      
900   if (ierr .ne. 0 .and. iprint .ge. 1) then
         write (*,'(a,i10)') 'Error reading header with error code: ', &
              & ierr
         write (*,'(a,i10)') '(-1: opening file; >0 : line number of error)'
      end if
      return
      end subroutine read_header
    
      subroutine copy_arr(xin,nin1,nin2,xout,nout1,nout2)

        implicit none 
        integer, intent(in) :: nin1,nin2,nout1,nout2
        real, intent(in) :: xin(nin1,nin2)
        real, intent(out) :: xout(nout1,nout2)

        integer :: i,j

        do j=1,min(nin2,nout2)
           do i=1,min(nin1,nout1)
              xout(i,j) = xin(i,j)
           end do
        end do
! Fill extra column(s), if any
        do j=nin2+1,nout2
           do i=1,min(nin1,nout1)
              xout(i,j) = xin(i,nin2)
           end do
        end do
! Fill extra row(s), if any
        do j=1,min(nin2,nout2)
           do i=nin1+1,nout1
              xout(i,j) = xin(nin1,j)
           end do
        end do
! Fill upper right corner, if needed
        do j=nin2+1,nout2
           do i=nin1+1,nout1
              xout(i,j) = xin(nin1,nin2)
           end do
        end do
        
        return
      end subroutine copy_arr
