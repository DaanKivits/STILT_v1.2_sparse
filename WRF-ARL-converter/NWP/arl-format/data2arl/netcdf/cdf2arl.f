!===============================================================================
! Convert netCDF files that are constructed with one variable per file to 
! packed ARL format with data organized by time such that a consecutive 
! group of records contains all the variables at the same time.  The ARL format
! consists of an index record followed by data records.  One record per        
! variable per level, then followed by records for the next time period.  All
! records are of fixed length and packed one byte per variable.  Packing
! information is coded in the header portion of each record. One month of data
! is processed per execution. A full global repacking (144x73) results in an
! output file of about 122 Mb per month. 
!
! The program repacks the input data into ARL format in two different modes. 
! If no grid spacing is defined (=0) then the latlon input grid is repacked 
! and written directly to ARL format.  A latlon subgrid may be specified by
! setting the number of ouput grid points to a non-zero value.  The subgrid
! is geolocated by setting the lower left corner point.        
!
! The latlon input data may also be interpolated to a conformal map projection
! by specifying the output grid size (km), grid center latlon, and number of
! grid points in each direction.  The projection is automatically defined 
! depending upon the latitude of the grid center, where Polar Sterographic
! is used for latitudes above 60, Mercator for latitudes less than 30, and 
! Lambert Conformal for intermediate latitudes.
!-------------------------------------------------------------------------------
! Last Revised: 08 Mar 2001
!               04 Apr 2001 - added Gaussian precipitation field
!-------------------------------------------------------------------------------

PROGRAM cdf2arl

  IMPLICIT NONE

  INTEGER, PARAMETER :: nvar  = 11      ! number of variables (and files)
  INTEGER, PARAMETER :: nsfc  = 5       ! number of surface variables
  INTEGER, PARAMETER :: kunit = 50      ! output unit for ARL packed data
  INTEGER, PARAMETER :: igmax = 192     ! number of gaussian grid longitudes
  INTEGER, PARAMETER :: jgmax = 94      ! number of gaussian grid latitudes
  INTEGER, PARAMETER :: lonpt = 144     ! number of longitude points
  INTEGER, PARAMETER :: latpt = 73      ! number of latitude points

  LOGICAL            :: global          ! flag to output global latlon grid
  INTEGER            :: klvls           ! number of levels to output
  INTEGER(8)         :: ctime           ! netCDF time indicator
! CHARACTER(8)       :: ctime           ! integer*8 not valid on some PCs
  REAL               :: gridkm          ! output grid size in km
  INTEGER            :: nxp, nyp        ! number of pts in output grid
  INTEGER            :: nxy, lrec       ! product of pts and record length
  REAL               :: clat, clon      ! center position of output grid
  LOGICAL            :: ftest           ! file exist test result
  CHARACTER(80)      :: fname           ! file name holder 
  INTEGER            :: ndat            ! number of elements (nx*ny*nz)
  INTEGER            :: nz              ! maximum number of levels
  INTEGER            :: kv              ! variable with max numb of levels
  INTEGER            :: i1, j1          ! lower left index of input grid
  INTEGER            :: n, k, m, lev    ! dummy indicies
  INTEGER            :: iy,im,id,ih     ! current obs time 
  INTEGER            :: ic, np          ! hours between obs, skip records
  INTEGER            :: julh            ! hours since start of year
  INTEGER            :: year, month     ! processing month
  INTEGER            :: fcopen          ! external file IO routines
  INTEGER            :: klen            ! dummy byte counter 
  INTEGER            :: krec            ! output record counter
  REAL               :: clat1, clon1    ! corner of input grid (1,1)
  REAL               :: dlat,  dlon     ! increment between points
  INTEGER            :: nlat,  nlon     ! number of input grid points

  INTEGER(2)         :: kdata(lonpt,latpt,17)
  REAL               :: rdata(lonpt,latpt,17)
  REAL               :: sdata(lonpt,latpt,17)

!---------------------------------------------------------
! gaussian and lat-lon grid arrays

  REAL           :: glat(jgmax),coa(jgmax),sia(jgmax),gw(jgmax)

  INTEGER(2)     :: mdata (igmax,jgmax) ! gaussian packed array
  REAL           :: gdata (igmax,jgmax) ! gaussian real array
  REAL           :: gxp   (igmax)       ! gaussian x points
  REAL           :: gyp   (jgmax)       ! gaussian y points
  REAL           :: ldata (lonpt,latpt) ! lat-lon interpolated
  REAL           :: lxp   (lonpt)       ! lat points
  REAL           :: lyp   (latpt)       ! lon points

!---------------------------------------------------------
! file dependent variables

  REAL           :: level  (20,nvar)    ! level values
  CHARACTER(8)   :: varb   (nvar)       ! data variable
  INTEGER        :: handle (nvar)       ! unit number
  INTEGER        :: kdate  (3,nvar)     ! starting year month day
  INTEGER        :: kdim   (4,nvar)     ! number in x,y,z,t
  INTEGER        :: vdata  (nvar)       ! byte position vertical levels
  INTEGER        :: nbyte  (nvar)       ! bytes per record variable    
  INTEGER        :: begin  (nvar)       ! byte position variable data   
  REAL           :: offset (nvar)       ! packing offset value    
  REAL           :: scale  (nvar)       ! packing scaling value   
  CHARACTER(4)   :: vchar  (nvar)       ! ARL equivalents    
  REAL           :: cnvrt  (nvar)       ! ARL unit conversion factors 

!---------------------------------------------------------
! conformal grid variables

  REAL, ALLOCATABLE :: tlat (:,:)       ! lat position of each point 
  REAL, ALLOCATABLE :: tlon (:,:)       ! lon position of each point
  REAL, ALLOCATABLE :: xvar (:,:)       ! dummy variable #1
  REAL, ALLOCATABLE :: yvar (:,:)       ! dummy variable #2

  CHARACTER(1), ALLOCATABLE :: cvar(:)  ! packed output array

!---------------------------------------------------------
! 4 surface variables followed by 6 3D variables

! CDC netCDF variables as expressed in file name
  DATA varb/'pres&','air&','uwnd&','vwnd&','prate&',              &
            'hgt&', 'air&','uwnd&','vwnd&','omega&','rhum&'/

! ARL packed equivalent names and unit conversion factors
! precip kg/m/s to m over a 6h accumulation 0.001 * 6 * 3600.0
  DATA vchar/'PRSS','T02M','U10M','V10M','TPP6',                  &
             'HGTS','TEMP','UWND','VWND','WWND','RELH'/
  DATA cnvrt/ 0.01 , 1.0  , 1.0  , 1.0  , 21.6,                   &
               1.0 , 1.0  , 1.0  , 1.0  , 0.01 , 1.0  /

!---------------------------------------------------------
  INTERFACE
  SUBROUTINE ncscan(handle,varid,kdate,kdim,vdata,                &
                    nbyte,begin,offset,scale,diag)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: handle
  CHARACTER(8), INTENT(IN)  :: varid
  INTEGER,      INTENT(OUT) :: kdate (3)
  INTEGER,      INTENT(OUT) :: kdim  (4)
  INTEGER,      INTENT(OUT) :: vdata 
  INTEGER,      INTENT(OUT) :: nbyte 
  INTEGER,      INTENT(OUT) :: begin 
  REAL,         INTENT(OUT) :: offset
  REAL,         INTENT(OUT) :: scale
  LOGICAL,      INTENT(IN)  :: diag
  END SUBROUTINE ncscan

  SUBROUTINE ncjulh (iy,im,id,ih,ic,julh)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: iy, im, id, ih
  INTEGER, INTENT(IN)    :: ic
  INTEGER, INTENT(OUT)   :: julh
  END SUBROUTINE ncjulh 

  SUBROUTINE MKGRID(NX,NY,TLAT,TLON)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NX,NY
  REAL,    INTENT(OUT) :: TLAT(:,:)
  REAL,    INTENT(OUT) :: TLON(:,:)
  END SUBROUTINE mkgrid

  SUBROUTINE SETNDX(GLOBAL,FTEST,KLVLS,GRIDKM,NXP,NYP,CLAT,CLON, &
                    NLAT,NLON,KUNIT)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT)   :: GLOBAL
  LOGICAL, INTENT(OUT)   :: FTEST
  INTEGER, INTENT(INOUT) :: KLVLS
  REAL,    INTENT(OUT)   :: GRIDKM
  INTEGER, INTENT(OUT)   :: NXP 
  INTEGER, INTENT(OUT)   :: NYP
  REAL,    INTENT(OUT)   :: CLAT  
  REAL,    INTENT(OUT)   :: CLON 
  INTEGER, INTENT(IN)    :: NLAT
  INTEGER, INTENT(IN)    :: NLON
  INTEGER, INTENT(IN)    :: KUNIT
  END SUBROUTINE setndx

  SUBROUTINE MAKNDX (nsfc,nvar,klvls,vchar,level,gridkm,nxp,nyp,clat,clon, &
                     global,clat1,clon1,dlat,dlon,nlat,nlon,kunit)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)   :: nsfc, nvar  
  INTEGER,      INTENT(IN)   :: klvls
  CHARACTER(4), INTENT(IN)   :: vchar (nvar)    
  REAL,         INTENT(IN)   :: level (20,nvar)
  REAL,         INTENT(IN)   :: gridkm
  INTEGER,      INTENT(IN)   :: nxp 
  INTEGER,      INTENT(IN)   :: nyp 
  REAL,         INTENT(IN)   :: clat  
  REAL,         INTENT(IN)   :: clon 
  LOGICAL,      INTENT(IN)   :: global
  REAL,         INTENT(IN)   :: clat1 
  REAL,         INTENT(IN)   :: clon1
  REAL,         INTENT(IN)   :: dlat  
  REAL,         INTENT(IN)   :: dlon 
  INTEGER,      INTENT(IN)   :: nlon
  INTEGER,      INTENT(IN)   :: nlat
  INTEGER,      INTENT(IN)   :: kunit
  END SUBROUTINE makndx

  SUBROUTINE REGRID(GLOBAL,I1,J1,V1,NX1,NY1,V2,NX2,NY2,TLAT,TLON,  &
                    CLAT1,CLON1,DLON,DLAT)
  IMPLICIT NONE
  LOGICAL, INTENT(IN)    :: GLOBAL
  INTEGER, INTENT(IN)    :: I1, J1   
  INTEGER, INTENT(IN)    :: NX1, NY1
  INTEGER, INTENT(IN)    :: NX2, NY2
  REAL,    INTENT(INOUT) :: V1(NX1,NY1)
  REAL,    INTENT(OUT)   :: V2(NX2,NY2)
  REAL ,   INTENT(IN)    :: TLAT(NX2,NY2)
  REAL ,   INTENT(IN)    :: TLON(NX2,NY2)
  REAL,    INTENT(IN)    :: CLAT1, CLON1
  REAL,    INTENT(IN)    :: DLAT,  DLON
  END SUBROUTINE regrid

  SUBROUTINE U2GRID(UU,VV,NX,NY)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY
  REAL, INTENT(INOUT) :: UU(NX,NY)
  REAL, INTENT(INOUT) :: VV(NX,NY)
  END SUBROUTINE u2grid

  END INTERFACE

!---------------------------------------------------------
! diagnostic output file

  OPEN(KUNIT+1,FILE='MESSAGE')

!---------------------------------------------------------
! select the year and month for processing

  WRITE(*,*)'Enter four digit year: '
  READ(*,*)year             
  WRITE(*,*)'Enter two digit month: '
  READ(*,*)month

!---------------------------------------------------------
! construct and open files based upon the year selected

  DO K=1,nvar
     klen=INDEX(varb(k),'&')-1
     IF(K.EQ.1)THEN
!       surface pressure
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.sfc.',year,'.nc'
     ELSEIF(K.GE.2.AND.K.LE.4)THEN
!       first sigma level variables
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.sig995.',year,'.nc'
     ELSEIF(K.EQ.5)THEN
!       gaussian flux variables - precipitation
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.sfc.gauss.',year,'.nc'
     ELSE     
!       remaining 3D variables
        WRITE(fname,'(2a,i4,a3)')varb(k)(1:klen),'.',year,'.nc'
     END IF
     WRITE(KUNIT+1,*)'Opening: ', fname(1:40)
     HANDLE(k)=FCOPEN(fname,'r')
  END DO

!---------------------------------------------------------
! scan each file for required starting information

  WRITE(KUNIT+1,*)
  DO k=1,nvar 
     CALL ncscan(handle(k),varb(k),kdate(:,k),kdim(:,k),vdata(k),              &
                 nbyte(k),begin(k),offset(k),scale(k),.false.)
     klen=INDEX(varb(k),'&')-1
     write(KUNIT+1,*)'Set: ',varb(k)(1:klen),kdate(:,k),kdim(:,k)
  END DO

!---------------------------------------------------------
! construct level information

  WRITE(KUNIT+1,*)
  WRITE(KUNIT+1,*)'Vertical data -'

  nz=0
  level=0.0  

  DO k=1,nvar
!    find variable with max # levels
     IF(kdim(3,k).GT.nz)THEN
        nz=kdim(3,k)
        kv=k
     END IF

!    save level information if >1
     IF(kdim(3,k).gt.1)THEN
        CALL fcptps(handle(k),vdata(k),*900)
        CALL fcread(handle(k),level(1,k),4,kdim(3,k),*900)
     END IF
  END DO

  k=nz
  DO WHILE (k.GT.0)
     write(KUNIT+1,'(i4,a,10i6)')k,' - ',(int(level(k,m)),m=nsfc+1,nvar)
     k=k-1
  END DO

!---------------------------------------------------------
! set remaining input grid parameters

  nlon=kdim(1,kv)
  nlat=kdim(2,kv)
  WRITE(KUNIT+1,*)
  WRITE(KUNIT+1,*)'Input grid: ',nlon,nlat
! initialize max level selection variable 
  klvls=nz+1

! corner point and grid spacing of CDC global data
! note that CDC data starts at +90 but the data array is
! transformed from south to north in subroutine regrid
  clat1=-90.0
  clon1=0.0
  dlat=2.5
  dlon=2.5

! initialize gaussian grid to lat-lon conversion
  CALL SETXY(lxp,lonpt,lyp,latpt,0.0,dlon,-90.0,dlat)
  CALL SETXY(gxp,igmax,gyp,jgmax,0.0,(360.0/FLOAT(igmax)),0.,0.)
  CALL GAUSSL(sia,coa,glat,gw,jgmax)
  gyp = glat * 180.0/3.14159265

!---------------------------------------------------------
! configure the ARL packed data

! read parameters for the output grid
  CALL SETNDX(global,ftest,klvls,gridkm,nxp,nyp,clat,clon,  &
              nlat,nlon,kunit)

! create the configuration file if it doesn't exist
  IF(.NOT.FTEST) CALL MAKNDX                                &
     (nsfc,nvar,klvls,vchar,level,gridkm,nxp,nyp,clat,clon, &
      global,clat1,clon1,dlat,dlon,nlat,nlon,kunit)

! configure the packing routines
  KREC=1
  CALL PAKSET(KUNIT,'METDATA.CFG',KREC,NXP,NYP,KLVLS)

! open output file (RS - reanalysis sigma; RP - reanalysis pressure)
  NXY =NXP*NYP
  LREC=NXY+50
  IF(GLOBAL)THEN
     WRITE(FNAME,'(a2,i4,i2.2,a4)')'RP',year,month,'.gbl'
  ELSE
     WRITE(FNAME,'(a2,i4,i2.2,a4)')'RP',year,month,'.arl'
  END IF
  OPEN(KUNIT,FILE=FNAME,RECL=LREC,ACCESS='DIRECT',   &
       FORM='UNFORMATTED')
  ALLOCATE (cvar(nxp*nyp))

! set up the grid system
  ALLOCATE (tlat(nxp,nyp),tlon(nxp,nyp))
  ALLOCATE (xvar(nxp,nyp),yvar(nxp,nyp))

! determine how data will be remapped
  IF(GLOBAL)THEN
!    lower left corner of global grid
     i1=1
     j1=1
     IF(nxp.lt.nlon)i1=1+(clon-clon1)/dlon
     IF(nyp.lt.nlat)j1=1+(clat-clat1)/dlat
     WRITE(KUNIT+1,*)'Dim of latlon output grid: ',nxp,nyp
     WRITE(KUNIT+1,*)'I,J of lower left corner : ',i1,j1
  ELSE
     CALL MKGRID(nxp,nyp,tlat,tlon)
     WRITE(KUNIT+1,*)'Dim of conformal output grid: ',nxp,nyp
  END IF

!---------------------------------------------------------
! construct time loop

  iy=kdate(1,kv)
  im=kdate(2,kv)
  id=kdate(3,kv)
  ih=0
  ic=6
  
! compute number of time periods to skip 
  CALL NCJULH (iy,month,id,ih,0,julh)
  NP=julh/ic
  im=month

! initial position for all input data files
  WRITE(KUNIT+1,*)
  DO n=1,nvar 
     begin(n)=begin(n)+(nbyte(n)+8)*np
     write(KUNIT+1,*)'File: ',n,'  skipping: ',(nbyte(n)+8)*np,'  to byte: ',begin(n)
!    back up one double integer for time record variaable
     CALL fcptps(handle(n),(begin(n)-8),*900)
  END DO

  WRITE(KUNIT+1,*)
  tloop : DO n=1,kdim(4,kv)
     WRITE(KUNIT+1,*)'Processing: ',iy,im,id,ih
     WRITE( *,*)'Processing: ',iy,im,id,ih

     vloop : DO k=1,nvar 

!       load time variable
        CALL fcread(handle(k),ctime,8,1,*900)

!       load data (assume short integer)
        ndat=kdim(1,k)*kdim(2,k)*kdim(3,k)
        IF(vchar(k).EQ.'TPP6')THEN
!          gaussian grid
           CALL fcread(handle(k),mdata,2,ndat,*900)
!          unpack the data from integer to real
           gdata =(mdata*scale(k))+offset(k)
        ELSE
!          standard lat-lon grid
           CALL fcread(handle(k),kdata,2,ndat,*900)
!          unpack the data from integer to real
           rdata =(kdata*scale(k))+offset(k)
        END IF

!       diagnostic output at bottom level 
        klen=INDEX(varb(k),'&')-1
        write(KUNIT+1,*)'   ',varb(k)(1:klen),' - ',rdata(1,1,1) 

!       shift pointer
        begin(k)=begin(k)+nbyte(k)

!       check for wind vector variables that require rotation
        klen=INDEX(varb(k),'&')-1
        IF(varb(k)(1:klen).EQ.'uwnd')THEN
           sdata = rdata
!          next variable must be v component
           CYCLE vloop
        END IF

!       write one record for each data level
        mloop : DO m=1,kdim(3,k)

!          ARL packed data level count includes surface
           lev=m
           if(kdim(3,k).gt.1)lev=m+1
           if(lev.gt.klvls) CYCLE mloop

!          special processing for Gaussian precipitation field
           IF(vchar(k).EQ.'TPP6')THEN
              CALL INTP2D(gdata,igmax,jgmax,gxp,gyp,ldata,lonpt,latpt,lxp,lyp,0)
!             convert to conformal grid
              CALL REGRID(global,i1,j1,ldata(1,1),nlon,nlat,yvar,nxp,nyp,      &
                          tlat,tlon,clat1,clon1,dlon,dlat)
           ELSE
!             convert to conformal grid
              CALL REGRID(global,i1,j1,rdata(1,1,m),nlon,nlat,yvar,nxp,nyp,    &
                          tlat,tlon,clat1,clon1,dlon,dlat)
           END IF

!          convert to ARL standard units
           yvar = yvar * cnvrt(k)          

           IF(vchar(k).EQ.'V10M'.OR.vchar(k).EQ.'VWND')THEN
!             convert the previously skipped u component
              CALL REGRID(global,i1,j1,sdata(1,1,m),nlon,nlat,xvar,nxp,nyp,    &
                          tlat,tlon,clat1,clon1,dlon,dlat)
              xvar = xvar * cnvrt(k-1)          

              IF(.NOT.GLOBAL) CALL U2GRID(xvar,yvar,nxp,nyp)

              CALL PAKREC(kunit,xvar,cvar,nxp,nyp,nxy,vchar(k-1),              &
                          mod(iy,100),im,id,ih,0,0,lev,0) 
              CALL PAKREC(kunit,yvar,cvar,nxp,nyp,nxy,vchar(k),                &
                          mod(iy,100),im,id,ih,0,0,lev,0) 
           ELSE
              CALL PAKREC(kunit,yvar,cvar,nxp,nyp,nxy,vchar(k),                &
                          mod(iy,100),im,id,ih,0,0,lev,0) 
           END IF

        END DO mloop

     END DO vloop

!    write index record to output file
     CALL PAKNDX(KUNIT)
     WRITE(KUNIT+1,*)'Index record written'
     WRITE(KUNIT+1,*)

!    increment internal clock
     CALL NCJULH (iy,im,id,ih,ic,julh)

!    only process one month per execution 
     IF(im.ne.month) EXIT tloop

  END DO tloop

!---------------------------------------------------------
! close all open files

  DO K=1,nvar 
     CALL FCCLOS(handle(k),*900)
  END DO
  CLOSE(KUNIT+1)

  900 CONTINUE

END PROGRAM cdf2arl

!===============================================================================
! Read METDATA.CFG configuration file for packing if required

SUBROUTINE SETNDX(GLOBAL,FTEST,KLVLS,GRIDKM,NXP,NYP,CLAT,CLON,NLAT,NLON,KUNIT)

  IMPLICIT NONE

  LOGICAL, INTENT(OUT)   :: GLOBAL
  LOGICAL, INTENT(OUT)   :: FTEST
  INTEGER, INTENT(INOUT) :: KLVLS
  REAL,    INTENT(OUT)   :: GRIDKM
  INTEGER, INTENT(OUT)   :: NXP 
  INTEGER, INTENT(OUT)   :: NYP
  REAL,    INTENT(OUT)   :: CLAT  
  REAL,    INTENT(OUT)   :: CLON 
  INTEGER, INTENT(IN)    :: NLAT
  INTEGER, INTENT(IN)    :: NLON
  INTEGER, INTENT(IN)    :: KUNIT

  INTEGER                :: K
  CHARACTER(80)          :: LABEL
  REAL                   :: GRIDS(12), PARMAP(9)

  COMMON / SETUP / GRIDS, PARMAP

! default is to map data to conformal projection
  GLOBAL=.FALSE.

! if configuration exists exit
  INQUIRE(FILE='METDATA.CFG',EXIST=FTEST)

  IF(FTEST)THEN
     OPEN(KUNIT,FILE='METDATA.CFG')
     READ(KUNIT,'(A)')LABEL
     READ(KUNIT,'(A)')LABEL
     READ(KUNIT,'(A)')LABEL
     DO K=1,12
        READ(KUNIT,'(A20,F10.2)')LABEL,GRIDS(K)
     END DO
     READ(KUNIT,'(A20,I4)')LABEL,NXP
     READ(KUNIT,'(A20,I4)')LABEL,NYP
     READ(KUNIT,'(A20,I4)')LABEL,KLVLS
     DO K=1,KLVLS
        READ(KUNIT,'(A)')LABEL
     END DO
     CLOSE (KUNIT)
     WRITE(*,*)'Using existing METFILE.CFG ... '
     WRITE(*,*)'Delete file and rerun program if new grid required'

!    definition of a global latlon grid (spacing=0.0)
     IF(GRIDS(5).LE.0.0)THEN
        GLOBAL=.TRUE.
        CLAT=GRIDS(10)
        CLON=GRIDS(11)
     END IF

  ELSE
     WRITE(*,*)'Output grid resolution (km or 0 for latlon)'
     READ(*,*)GRIDKM

     IF(GRIDKM.LE.0.0)THEN
!       defines that the output grid should be latlon 
        GLOBAL=.TRUE.

        WRITE(*,*)'Output grid dimensions (nlon,nlat or 0,0 if global)'
        READ(*,*)NXP,NYP

        IF(NXP.EQ.0.AND.NYP.EQ.0)THEN
           CLAT=0.0
           CLON=0.0
        ELSE 
           WRITE(*,*)'Output grid lower left corner (lat,lon)'
           READ(*,*)CLAT,CLON
!          use 0-360 system  
           IF(CLON.LT.0.0)CLON=360.0+CLON
        END IF
        IF(NXP.EQ.0)NXP=NLON
        IF(NYP.EQ.0)NYP=NLAT

        WRITE(*,*)'Number of levels required (incl sfc) - ',klvls
        READ(*,*)KLVLS

     ELSE
!       defines that the latlon grid is remapped to conformal for output
        WRITE(*,*)'Output grid dimensions (nx,ny)'
        READ(*,*)NXP,NYP
        WRITE(*,*)'Output grid center (lat,lon)'
        READ(*,*)CLAT,CLON
        WRITE(*,*)'Number of levels required (incl sfc) - ',klvls
        READ(*,*)KLVLS

     END IF
  END IF

END SUBROUTINE setndx

!===============================================================================
! Determine the lat/lon at grid intersections

SUBROUTINE MKGRID(NX,NY,TLAT,TLON)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NX,NY
  REAL,    INTENT(OUT) :: TLAT(:,:)
  REAL,    INTENT(OUT) :: TLON(:,:)

  INTEGER              :: I,J
  REAL                 :: GRIDS(12), PARMAP(9)

  COMMON / SETUP / GRIDS, PARMAP

! define the tangent latitude and reference longitude
  CALL STLMBR(PARMAP,GRIDS(7),GRIDS(4))

! define the grid by a one-point specification
  CALL STCM1P(PARMAP,GRIDS(8),GRIDS(9),GRIDS(10),GRIDS(11),                &
                     GRIDS(3),GRIDS(4),GRIDS(5),GRIDS(6))

! determine the lat/lon at the grid locations
  DO J=1,NY
  DO I=1,NX

!    cxy2ll returns -180 to +180
     CALL CXY2LL(PARMAP,FLOAT(I),FLOAT(J),TLAT(I,J),TLON(I,J))
!    shift to 0 to 360 system
     IF(TLON(I,J).LT.0.0)TLON(I,J)=360.0+TLON(I,J)

  END DO
  END DO

END SUBROUTINE mkgrid

!===============================================================================
! Create METDATA.CFG configuration file for packing if required

SUBROUTINE MAKNDX (nsfc,nvar,klvls,vchar,level,gridkm,nxp,nyp,clat,clon, &
                   global,clat1,clon1,dlat,dlon,nlat,nlon,kunit)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)   :: nsfc, nvar  
  INTEGER,      INTENT(IN)   :: klvls
  CHARACTER(4), INTENT(IN)   :: vchar (nvar)    
  REAL,         INTENT(IN)   :: level (20,nvar)
  REAL,         INTENT(IN)   :: gridkm
  INTEGER,      INTENT(IN)   :: nxp 
  INTEGER,      INTENT(IN)   :: nyp 
  REAL,         INTENT(IN)   :: clat  
  REAL,         INTENT(IN)   :: clon 
  LOGICAL,      INTENT(IN)   :: global
  REAL,         INTENT(IN)   :: clat1 
  REAL,         INTENT(IN)   :: clon1
  REAL,         INTENT(IN)   :: dlat  
  REAL,         INTENT(IN)   :: dlon 
  INTEGER,      INTENT(IN)   :: nlon
  INTEGER,      INTENT(IN)   :: nlat
  INTEGER,      INTENT(IN)   :: kunit
  
  CHARACTER(20) :: LABEL(18)                ! optional field label
  INTEGER       :: i,n,ig,nl,nvb,nve,kvc,mvar 
  REAL          :: sig
  REAL          :: GRIDS(12), PARMAP(9)

  COMMON / SETUP / GRIDS, PARMAP

! optional field label string
  DATA LABEL/'Model Type:','Grid Numb:','Vert Coord:','Pole Lat:',         &
    'Pole Lon:','Ref Lat:','Ref Lon:','Grid Size:','Orientation:',         &
    'Cone Angle:','Sync X Pt:','Sync Y Pt:','Sync Lat:','Sync Lon:',       &
    'Reserved:','Numb X pt:','Numb Y pt:','Numb Levels:'/

! grid orientation
  GRIDS(6)=0.0
! delta=x grid size in km
  GRIDS(5)=GRIDKM
! variable reserved for future use
  GRIDS(12)=0.0

  IF(GLOBAL)THEN
!    sync x,y defines lower left grid point 
     GRIDS(8)=1.0 
     GRIDS(9)=1.0 
!    lat/lon of lower left point
     GRIDS(10)=CLAT1
     GRIDS(11)=CLON1
!    latlon grid not global
     IF(NYP.LT.NLAT)GRIDS(10)=CLAT
     IF(NXP.LT.NLON)GRIDS(11)=CLON
  ELSE
!    synch point in x,y coordintes
     GRIDS(8)=(NXP+1.0)/2.0
     GRIDS(9)=(NYP+1.0)/2.0
!    synch point in lat/lon coordinates
     GRIDS(10)=CLAT
     GRIDS(11)=CLON
  END IF

! defines a global latlon grid
  IF(GLOBAL)THEN

!    pole lat/lon is used to identify the 
!    latlon point of the maximum index
     GRIDS(1)=GRIDS(10)+DLAT*(NYP-1)
     GRIDS(2)=GRIDS(11)+DLON*(NXP-1)
     GRIDS(7)=0.0  

!    the reference lat/lon defines grid spacing  
     GRIDS(3)=DLAT
     GRIDS(4)=DLON

! defines a polar sterographic projection
  ELSEIF(ABS(CLAT).GT.60.0)THEN

!    set the pole position and reference lat/lon
     IF(CLAT.GT.0.0)THEN
        GRIDS(1)=90.0
     ELSE
        GRIDS(1)=-90.0
     END IF

!    pole longtitude (+180 from cut)
     GRIDS(2)=CLON

!    reference lat/lon (at which grid size specified)
     GRIDS(3)=CLAT
!    reference longitude and grid alignment
     GRIDS(4)=CLON

!    tangent latitude
     GRIDS(7)=GRIDS(1)

! defines a mercator projection
  ELSEIF(ABS(CLAT).LT.30.0)THEN

!    pole lat/lon axis through pole
     GRIDS(1)=0.0
     GRIDS(2)=CLON

!    reference lat
     GRIDS(3)=CLAT
!    reference lon
     GRIDS(4)=CLON

!    tangent latitude
     GRIDS(7)=0.0

! defines a lambert conformal projection
  ELSE

!    pole lat/lon axis through pole
     GRIDS(1)=CLAT
     GRIDS(2)=CLON

!    reference lat
     GRIDS(3)=CLAT
!    reference lon
     GRIDS(4)=CLON

!    tangent latitude
     GRIDS(7)=CLAT

  END IF

! write the packer configuration file
  OPEN(KUNIT,FILE='METDATA.CFG')
  WRITE(KUNIT,'(A20,A4)')LABEL(1),'CDC1'

! default grid number 99 (field not used)
  IG=99
  WRITE(KUNIT,'(A20,I4)') LABEL(2),IG

! coordinate- 1:sigma  2:pressure  3:terrain  4:hybrid
  KVC=2
  WRITE(KUNIT,'(A20,I4)') LABEL(3), KVC

! grid geolocation parameters and projection
  DO I=1,12
     WRITE(KUNIT,'(A20,F10.2)')LABEL(I+3),GRIDS(I)
  END DO

! grid dimensions
  WRITE(KUNIT,'(A20,I4)') LABEL(16),NXP
  WRITE(KUNIT,'(A20,I4)') LABEL(17),NYP
  WRITE(KUNIT,'(A20,I4)') LABEL(18),KLVLS

! upper level information
  DO NL=1,KLVLS

     WRITE(LABEL(1),'(A6,I2,A1)')'Level ',NL,':'
     IF(NL.EQ.1)THEN
        MVAR=NSFC
        NVB=1
        NVE=MVAR
        SIG=0.0     
     ELSE
!       all variables available
        NVB=NSFC+1
        NVE=NVAR

!       count actual number
        MVAR=0
        DO N=NVB,NVE
           IF(LEVEL(NL-1,N).GT.0.0)MVAR=MVAR+1
        END DO
        NVE=NVB+MVAR-1

!       first upper level always complete
        SIG=LEVEL(NL-1,NSFC+1)
     END IF

     IF(SIG.LT.1.0)THEN
        WRITE(KUNIT,'(A20,F6.5,I3,10(1X,A4))')LABEL(1),                 &
              SIG,MVAR,(VCHAR(N),N=NVB,NVE)

     ELSEIF(SIG.GE.1.AND.SIG.LT.10.0)THEN
        WRITE(KUNIT,'(A20,F6.4,I3,10(1X,A4))')LABEL(1),                 &
              SIG,MVAR,(VCHAR(N),N=NVB,NVE)

     ELSEIF(SIG.GE.10.AND.SIG.LT.100.0)THEN
        WRITE(KUNIT,'(A20,F6.3,I3,10(1X,A4))')LABEL(1),                 &
              SIG,MVAR,(VCHAR(N),N=NVB,NVE)

     ELSEIF(SIG.GE.100.AND.SIG.LT.1000.0)THEN
        WRITE(KUNIT,'(A20,F6.2,I3,10(1X,A4))')LABEL(1),                 &
              SIG,MVAR,(VCHAR(N),N=NVB,NVE)

     ELSEIF(SIG.GE.1000)THEN
        WRITE(KUNIT,'(A20,F6.1,I3,10(1X,A4))')LABEL(1),                 &
              SIG,MVAR,(VCHAR(N),N=NVB,NVE)
     END IF

  END DO
  CLOSE (KUNIT)

END SUBROUTINE makndx

!===============================================================================
! Interpolate data to grid

SUBROUTINE REGRID(GLOBAL,I1,J1,V1,NX1,NY1,V2,NX2,NY2,TLAT,TLON,    &
                  CLAT1,CLON1,DLON,DLAT)

   IMPLICIT NONE

   LOGICAL, INTENT(IN)    :: GLOBAL
   INTEGER, INTENT(IN)    :: I1, J1   
   INTEGER, INTENT(IN)    :: NX1, NY1
   INTEGER, INTENT(IN)    :: NX2, NY2
   REAL,    INTENT(INOUT) :: V1(NX1,NY1)
   REAL,    INTENT(OUT)   :: V2(NX2,NY2)
   REAL ,   INTENT(IN)    :: TLAT(NX2,NY2)
   REAL ,   INTENT(IN)    :: TLON(NX2,NY2)
   REAL,    INTENT(IN)    :: CLAT1, CLON1
   REAL,    INTENT(IN)    :: DLAT,  DLON

   INTEGER :: i,j,ii,jj,ilo,ihi,jlo,jhi
   REAL    :: xp,yp,fxi,fyj,top,bot,temp

!  reverse the direction of the north-south component
   DO I=1,NX1
      JLO=1
      JHI=NY1
      DO WHILE (JLO.LT.JHI)
         TEMP=V1(I,JHI)
         V1(I,JHI)=V1(I,JLO)
         V1(I,JLO)=TEMP         
         JLO=JLO+1
         JHI=JHI-1
      END DO 
   END DO

! just move data into output array if latlon grid 
  IF(GLOBAL)THEN
     DO J=1,NY2  
        JJ=J1+J-1
!       assume bad definition if above 90
        JJ=MIN(JJ,NY1) 
        DO I=1,NX2   
           II=I1+I-1
!          assume wrap if grid goes around dateline
           IF(II.GT.NX1)II=II-NX1
           V2(I,J)=V1(II,JJ)        
        END DO
     END DO
     RETURN 
  END IF

! interpolate values to new grid
  DO I=1,NX2
  DO J=1,NY2

!    compute adjacent index values on grid 1
     XP=1.0+(TLON(I,J)-CLON1)/DLON
     YP=1.0+(TLAT(I,J)-CLAT1)/DLAT

!    compute indicies
     ILO=INT(XP)
     IHI=ILO+1
     JLO=INT(YP)
     JHI=JLO+1

!    interpolation fractions (off grid extrapolated)
     FXI=XP-ILO
     FYJ=YP-JLO

!    global grids check for wrap at prime meridian
     IF(IHI.GT.NX1)IHI=1
     IF(JHI.GT.NY1)JHI=1

!    interpolate across at top and bottom
     TOP=(V1(IHI,JHI)-V1(ILO,JHI))*FXI+V1(ILO,JHI)
     BOT=(V1(IHI,JLO)-V1(ILO,JLO))*FXI+V1(ILO,JLO)

!    interpolate between top and bottom
     V2(I,J)=(TOP-BOT)*FYJ+BOT

  END DO
  END DO

END SUBROUTINE regrid

!===============================================================================
! Convert wind components from true orientation to grid

SUBROUTINE U2GRID(UU,VV,NX,NY)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NX, NY
  REAL, INTENT(INOUT) :: UU(NX,NY)
  REAL, INTENT(INOUT) :: VV(NX,NY)

  INTEGER :: i,j
  REAL    :: ug,vg 
  REAL    :: GRIDS(12), PARMAP(9)

  COMMON/ SETUP / GRIDS, PARMAP

! convert compass winds to grid-orientation
  DO I=1,NX
  DO J=1,NY

     CALL CC2GXY(PARMAP,FLOAT(I),FLOAT(J),UU(I,J),VV(I,J),UG,VG)
     UU(I,J)=UG
     VV(I,J)=VG

  END DO
  END DO

END SUBROUTINE u2grid
