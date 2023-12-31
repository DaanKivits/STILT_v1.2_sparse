!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM: AVN2ARL      DECODE GRIB AVN MODEL FIELD FOR HYSPLIT
!   PRGMMR: DRAXLER          ORG: R/ARL       DATE: 1998-08-26
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
! Convert GRIB to ARL format with data organized by time such that a consecutive 
! group of records contains all the variables at the same time.  The ARL format
! consists of an index record followed by data records.  One record per        
! variable per level, then followed by records for the next time period.  All
! records are of fixed length and packed one byte per variable.  Packing
! information is coded in the header portion of each record. 
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
! Lambert Conformal for intermediate latitudes. The conformal map is at 100km 
! resolution of 100x100 points centered at the command line lat/lon point.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD)
!                  8 Jul 1998 (RRD) - avn extra levels relh
!                  3 Feb 1999 (RRD) - fixed bug in mercator configuration
!                                   - added lambert conformal option
!                                   - prime meridian interpolation patch
!                 12 Nov 1999 (RRD) - pc argument list compatibility
!                 10 Mar 2000 (RRD) - standardized PC/UNIX IO
!                 12 Apr 2000 (RRD) - increased resolution and fields
!                 21 Dec 2000 (RRD) - fortran90 upgrade
!                 07 Mar 2001 (RRD) - GRIB id test
!                 15 Nov 2001 (RRD) - global grid option added
!                 05 Mar 2002 (RRD) - improved w3lib consistency
!                 09 Apr 2002 (RRD) - grib test
!                 10 Jul 2002 (BJS) - add fields, levels
!
! USAGE:  AVN2ARL [-options]
!   INPUT ARGUMENT LIST:
!         -i [grib input file]
!         -x [center longitude]
!         -y [center latitude]
!         -g [1:nhem 2:shem 3:global]
!         -n [number of x,y extract grid points]
!   OUTPUT ARGUMENT LIST: none 
!   INPUT FILES:          grib input data files 
!   OUTPUT FILES:
!      unit 20 DATA.AVN - ARL packed data output file
!      unit 30 CFG_AVN  - characteristics of output grid
!      unit 40 AVNTIME  - time indicator file
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

PROGRAM AVN2ARL

  IMPLICIT NONE

  LOGICAL       :: FTEST
  CHARACTER(80) :: LABEL, FNAME
  INTEGER       :: IARGC, HANDLE, FCOPEN, NARG, KGRID, NUMPTS
  REAL          :: CLAT, CLON

  INTERFACE
    SUBROUTINE XTRACT(HANDLE,CLAT,CLON,KGRID,NUMPTS)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: HANDLE       ! IO handle
    INTEGER,INTENT(IN)  :: KGRID        ! projection type
    INTEGER,INTENT(IN)  :: NUMPTS       ! extract grid size
    REAL,   INTENT(IN)  :: CLAT, CLON   ! center position
    END SUBROUTINE xtract
  END INTERFACE

! required for NCEP operational implementation
! CALL W3TAGB('AVN2ARL ',1998,0238,0067,'R/ARL  ')

! check for command line arguments
  NARG=IARGC()

  IF(NARG.EQ.0)THEN
     WRITE(*,*)'Usage: avn2arl [-options]'
     WRITE(*,*)' -i[grib input file]'
     WRITE(*,*)' -x[extract center longitude]'
     WRITE(*,*)' -y[extract center latitude]'
     WRITE(*,*)' -g[0:conformal 1:nhem 2:shem 3:global]'
     WRITE(*,*)' -n[number of x,y extract grid points]'
     STOP
  END IF

! default values
  CLAT= 60.0
  CLON=-80.0
  KGRID=0 
  NUMPTS=100

! go through each argument
  DO WHILE (NARG.GT.0)

     CALL GETARG(NARG,LABEL)
     SELECT CASE (LABEL(1:2))

!    grib input file name   
     CASE ('-i','-I')
        READ(LABEL(3:),'(A)')FNAME

!    extract grid center longitude 
     CASE ('-x','-X')
        READ(LABEL(3:),'(F10.0)')CLON  

!    extract grid center latitude    
     CASE ('-y','-Y')
        READ(LABEL(3:),'(F10.0)')CLAT  

!    extract grid selection
     CASE ('-g','-G')
        READ(LABEL(3:),'(I1)')KGRID 

!    extract grid dimension
     CASE ('-n','-N')
        READ(LABEL(3:),'(I3)')NUMPTS

     END SELECT
     NARG=NARG-1

  END DO

! open input data file
  INQUIRE(FILE=FNAME,EXIST=FTEST)
  IF(FTEST)THEN
!    direct IO subroutine
     HANDLE=FCOPEN(FNAME,'r')

!    main decoding routine (assume one time period per process)
     WRITE(*,*)'started processing: ',FNAME
     CALL XTRACT(HANDLE,CLAT,CLON,KGRID,NUMPTS)

!    direct IO subroutine
     CALL FCCLOS(HANDLE,*900)
900  CONTINUE
  ELSE
     WRITE(*,*)'File not found:',FNAME
     STOP
  END IF

! required for NCEP operational implementation
! CALL W3TAGE('AVN2ARL ')

END PROGRAM avn2arl



!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EXTRACT          MAIN ROUTINE FOR DECODING AVN GRIB DATA
!   PRGMMR:    ROLAND DRAXLER   ORG: R/E/AR      DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            EXTRACTS ONE RECORD FROM INPUT GRIB FILE AND ADDS ONE
!            RECORD TO THE OUTPUT FILE CENTERED ABOUT INDICATED LAT/LON
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 10 Jul 1997 (RRD)
!                 11 Apr 2000 (RRD) - increased vertical resolution
!                 15 Nov 2001 (RRD) - global grid options
!                 10 Jul 2002 (BJS) - add fields, levels
!                 25 Feb 2003 (RRD) - ichar function replacement
!
! USAGE:  CALL XTRACT(HANDLE,CLAT,CLON,KGRID,NUMPTS)
!   INPUT ARGUMENT LIST:  see below
!   OUTPUT ARGUMENT LIST: none
!   INPUT FILES:          none
!   OUTPUT FILES:         none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE XTRACT(HANDLE,CLAT,CLON,KGRID,NUMPTS)

  IMPLICIT NONE

! argument list variables  
  INTEGER,INTENT(IN)  :: HANDLE       ! IO handle
  INTEGER,INTENT(IN)  :: KGRID        ! projection type
  INTEGER,INTENT(IN)  :: NUMPTS       ! extract grid size
  REAL,INTENT(INOUT)  :: CLAT, CLON   ! center position

! input buffer 
  INTEGER, PARAMETER   :: NLON=360
  INTEGER, PARAMETER   :: NLAT=181
  INTEGER, PARAMETER   :: MAXB=150000

! output array limits
  INTEGER, PARAMETER   :: NLVL=23     ! number of levels
  INTEGER, PARAMETER   :: N2DV=12     ! number of 2D variables 
  INTEGER, PARAMETER   :: N3DV=6      ! number of 3D variables

! temporary holding variable
  INTEGER      :: VGRIB
  CHARACTER(4) :: VCHAR
  REAL         :: CNVRT

! sfc arrays to hold grib, character, and level information
  INTEGER      :: VGRIB0(N2DV)        ! sfc grib ID  
  INTEGER      :: STYP  (N2DV)        ! sfc variable code 
  INTEGER      :: SIG0  (N2DV)        ! sfc level code 
  CHARACTER(4) :: VCHAR0(N2DV)        ! arl variable string
  REAL         :: CNVRT0(N2DV)        ! units conversion factor

! 3d arrays to hold grib, character, and level information
  INTEGER      :: VGRIB1(N3DV)        ! 3d grib variable code
  CHARACTER(4) :: VCHAR1(N3DV)        ! arl variable string
  REAL         :: CNVRT1(N3DV)        ! units conversion factor
  INTEGER      :: SIGL  (NLVL)        ! input/output levels
  INTEGER      :: NVAR  (NLVL)        ! variables each level

! arrays for grib description
  INTEGER      :: KPDS(25)            ! product definition
  INTEGER      :: KGDS(25)            ! grid definition
  INTEGER      :: KPTR(25)

! input data buffer and bit map section
  CHARACTER(1) :: BUFF(MAXB)
  LOGICAL(1)   :: KBMS(MAXB)

! unpacked input and output array
  REAL                     :: RVAR(NLON,NLAT)
  REAL,ALLOCATABLE         :: XVAR(:,:)      

! lat/lon grid conversion equivalence table
  REAL,ALLOCATABLE         :: TLAT(:,:)
  REAL,ALLOCATABLE         :: TLON(:,:)       

! packed output array
  CHARACTER(1),ALLOCATABLE :: CVAR(:)

! file information
  CHARACTER(80):: FNAME
  LOGICAL      :: FTEST

  LOGICAL :: global 
  INTEGER :: kvarb,koff,kunit,klvls,ltype,level 
  INTEGER :: lrec,nrec,nxy,nxp,nyp
  INTEGER :: k,i1,j1,iyr,imo,ida,ihr,ifh,imn 
  INTEGER :: nx,ny,nz,kl,kv,kret,krec,klen,kbyte 
  REAL    :: clon1,clat1,dlon,dlat,gridkm

! remap array from one to two dimensions
  REAL SVAR(NLON*NLAT)
  INTEGER SHAPE(2)
  DATA SHAPE/NLON,NLAT/

!-------------------------------------------------------------------------------

! surface variable 
  DATA VGRIB0/1,  2, 61,124,125,122,204,  52,  33,  34,  11, 71 /     ! GRIB ID
  DATA SIG0/  0,  0,  0,  0,  0,  0,  0,   2,  10,  10,   2,  0 /     ! level id
  DATA STYP/  1,102,  1,  1,  1,  1,  1, 105, 105, 105, 105, 200/     ! type 
  DATA VCHAR0  &
   /'PRSS','MSLP','TPP6','UMOF','VMOF','SHTF','DSWF','RH2M', &
    'U10M','V10M','T02M','TCLD'/
  DATA CNVRT0  &
   /0.01 , 0.01 ,0.001 , 1.0  ,  1.0 , 1.0  , 1.0  , 1.0  , 1.0 , 1.0, 1.0, 1.0/

! upper level variables
  DATA VGRIB1/    7,    11,    33,    34,    39,    52 /    ! GRIB ID
  DATA VCHAR1/'HGTS','TEMP','UWND','VWND','WWND','RELH'/    ! ARL ID
  DATA CNVRT1/  1.0,   1.0,   1.0,   1.0,  0.01,   1.0 /    ! CONVERSION

! set levels for input and output
  DATA SIGL/1000,975,950,925,900,850,800,750,700,650,600,550,500,450,400, &
             350,300,250,200,150,100,50,20/

! predefine number of variables aloft 
  DATA NVAR / 21*6, 2*4 /			

! lat/lon grid corner (1,1) and increment
  DATA CLAT1/90.0/, CLON1/0.0/, DLAT/1.0/, DLON/1.0/

! output record counter
  DATA KREC/0/
  SAVE KREC

!-------------------------------------------------------------------------------
! When dealing with some F90 compilers, replace ICHAR below with JCHAR function
  CHARACTER(1)        :: mychr    
  INTEGER             :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
!-------------------------------------------------------------------------------

  INTERFACE

  SUBROUTINE REGRID(GLOBAL,V1,NX1,NY1,V2,NX2,NY2,TLAT,TLON,    &
                    CLAT1,CLON1,DLON,DLAT)
  IMPLICIT NONE
  LOGICAL, INTENT(IN)    :: GLOBAL          ! global lat/lon grid 
  INTEGER, INTENT(IN)    :: NX1, NY1        ! input dimensions
  INTEGER, INTENT(IN)    :: NX2, NY2        ! output dimensions
  REAL,    INTENT(INOUT) :: V1(:,:)         ! input array
  REAL,    INTENT(OUT)   :: V2(:,:)         ! output array
  REAL ,   INTENT(IN)    :: TLAT(:,:)       ! position for each node 
  REAL ,   INTENT(IN)    :: TLON(:,:)       ! position for each node 
  REAL,    INTENT(IN)    :: CLAT1, CLON1    ! position of corner (1,1)
  REAL,    INTENT(IN)    :: DLAT,  DLON     ! spacing between input grid
  END SUBROUTINE regrid

  SUBROUTINE MKGRID(NX,NY,TLAT,TLON)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NX,NY      ! output grid dimensions
  REAL,    INTENT(OUT) :: TLAT(:,:)  ! position of nodes
  REAL,    INTENT(OUT) :: TLON(:,:)  ! position of nodes
  END SUBROUTINE mkgrid

  SUBROUTINE MAKNDX (N2DV,NVAR,VCHAR0,VCHAR1,LEVEL,NXP,NYP,NZP,CLAT,CLON, &
                     KGRID,GLOBAL,DLAT,DLON,NLAT,NLON)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)   :: n2dv          ! number sfc variables 
  INTEGER,      INTENT(IN)   :: nvar   (:)    ! number of 3d variables per level
  CHARACTER(4), INTENT(IN)   :: vchar0 (:)    ! sfc variable id
  CHARACTER(4), INTENT(IN)   :: vchar1 (:)    ! 3d variable id
  INTEGER,      INTENT(IN)   :: level  (:)    ! level information
  INTEGER,      INTENT(IN)   :: nxp           ! x dimension
  INTEGER,      INTENT(IN)   :: nyp           ! y dimension
  INTEGER,      INTENT(IN)   :: nzp           ! z dimension
  REAL,         INTENT(IN)   :: clat          ! center grid lat
  REAL,         INTENT(IN)   :: clon          ! center grid lon
  INTEGER,      INTENT(IN)   :: kgrid         ! grid type 
  LOGICAL,      INTENT(IN)   :: global        ! global lat/lon flag
  REAL,         INTENT(IN)   :: dlat          ! lat spacing 
  REAL,         INTENT(IN)   :: dlon          ! lon spacing
  INTEGER,      INTENT(IN)   :: nlat          ! number of lats
  INTEGER,      INTENT(IN)   :: nlon          ! number of lons 
  END SUBROUTINE makndx

  SUBROUTINE ROTATE(NREC,NXP,NYP)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NREC       ! number of records per time period
  INTEGER, INTENT(IN) :: NXP,NYP    ! dimensions of real arrays
  END SUBROUTINE rotate

  SUBROUTINE PAKSET(LUNIT,FNAME,KREC1,NXP,NYP,NZP)
  IMPLICIT NONE
  INTEGER,       INTENT(IN)    :: lunit     ! output unit number
  CHARACTER(*),  INTENT(INOUT) :: fname     ! file name of METDATA.CFG
  INTEGER,       INTENT(IN)    :: krec1     ! position of index record at time-1
  INTEGER,       INTENT(OUT)   :: nxp, nyp  ! horizontal grid dimensions
  INTEGER,       INTENT(OUT)   :: nzp       ! vertical grid dimension (incl sfc)
  END SUBROUTINE pakset

  SUBROUTINE PAKREC(LUNIT,RVAR,CVAR,NX,NY,NXY,KVAR,IY,IM,ID,IH,MN,IC,LL,KINI)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: LUNIT       ! output unit number
  INTEGER,      INTENT(IN)  :: NX,NY       ! dimensions of RVAR
  INTEGER,      INTENT(IN)  :: NXY         ! dimensions of CVAR
  REAL,         INTENT(IN)  :: RVAR(NX,NY) ! input data to be packed
  CHARACTER(1), INTENT(OUT) :: CVAR(NXY)   ! packed data array
  CHARACTER(4), INTENT(IN)  :: KVAR        ! descriptor of variable written
  INTEGER,      INTENT(IN)  :: IY,IM,ID    ! date identification
  INTEGER,      INTENT(IN)  :: IH,MN       ! time identification (MN-minutes)
  INTEGER,      INTENT(IN)  :: IC          ! forecast hour, ICX hour for >99
  INTEGER,      INTENT(IN)  :: LL          ! level indicator 
  INTEGER,      INTENT(IN)  :: KINI        ! initialization (0-no 1-yes)
  END SUBROUTINE pakrec

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

!-------------------------------------------------------------------------------

! set up the grid system
  GLOBAL=.FALSE.
  IF(KGRID.EQ.0)THEN
!    conformal extract  
     NXP=NUMPTS
     NYP=NUMPTS
     GRIDKM=100.0
  ELSEIF(KGRID.EQ.1)THEN
!    northern hemisphere polar stereographic
     NXP=257 
     NYP=257
     CLAT=90.0
     CLON=0.0
     GRIDKM=95.25
  ELSEIF(KGRID.EQ.2)THEN
!    southern hemisphere polar stereographic
     NXP=257 
     NYP=257
     CLAT=-90.0
     CLON=0.0
     GRIDKM=95.25
  ELSEIF(KGRID.EQ.3)THEN
!    global latitude-longitude grid (same as input)
     NXP=NLON
     NYP=NLAT
     GRIDKM=0.0   
     GLOBAL=.TRUE.
  END IF

! reserve space for output data arrays
  NXY=NXP*NYP
  ALLOCATE (cvar(nxy))
  ALLOCATE (tlat(nxp,nyp),tlon(nxp,nyp))
  ALLOCATE (xvar(nxp,nyp))

!-------------------------------------------------------------------------------

! create the configuration file if it doesn't exist
  CALL MAKNDX(N2DV,NVAR,VCHAR0,VCHAR1,SIGL,NXP,NYP,NLVL,CLAT,CLON, &
              KGRID,GLOBAL,DLAT,DLON,NLAT,NLON)

! configure the packing routines
  FNAME='CFG_AVN'
  CALL PAKSET(20,FNAME,1,NX,NY,NZ)
  WRITE(*,*)'Set grid from pakset: ',nx,ny,nz

! determine how the input data will be remapped 
  IF(.NOT.GLOBAL) CALL MKGRID(NXP,NYP,TLAT,TLON)

! standard output name for packed data
  INQUIRE(FILE='DATA.AVN',EXIST=FTEST)
  LREC=NXY+50
  OPEN(20,FILE='DATA.AVN',RECL=LREC,ACCESS='DIRECT',FORM='UNFORMATTED')

!-------------------------------------------------------------------------------

! input grib record byte counter
  KBYTE=0

!   read the indicator section
100 CONTINUE

!   direct IO subroutines
    CALL FCPTPS(HANDLE,KBYTE,*900)
    CALL FCREAD(HANDLE,BUFF,1,8,*900)
    IF(BUFF(1)//BUFF(2)//BUFF(3)//BUFF(4).NE.'GRIB')THEN
       KBYTE=KBYTE+1
       GOTO 100
    END IF

!   determine the length of the entire grib record
    KLEN=JCHAR(BUFF(7))+JCHAR(BUFF(6))*256+JCHAR(BUFF(5))*65536
    IF(KLEN.GT.MAXB)THEN
       WRITE(*,*)'Grib record: ',KLEN
       WRITE(*,*)'Exceedes buffer: ',MAXB
       STOP
    END IF

!   load the indicator section and pds segment into the buffer
!   direct IO subroutines
    CALL FCPTPS(HANDLE,KBYTE,*900)
    CALL FCREAD(HANDLE,BUFF,1,36,*900)

!   product definition section (+8 byte indicator section offset)
    KOFF=8
    KVARB=JCHAR(BUFF(KOFF+9))
    LTYPE=JCHAR(BUFF(KOFF+10))
    LEVEL=JCHAR(BUFF(KOFF+12))+JCHAR(BUFF(KOFF+11))*256

! check if 2d variable present in selection table

    KL=0
    DO KV=1,N2DV 
       VGRIB=VGRIB0(KV)
       VCHAR=VCHAR0(KV)
       CNVRT=CNVRT0(KV)
!      matches id and special level indicator
       IF(KVARB.EQ.VGRIB.AND.LEVEL.EQ.SIG0(KV).AND.                          &
          STYP(KV).EQ.LTYPE)GO TO 300
    END DO

! then check for 3d variable

    DO KV=1,N3DV 
       VGRIB=VGRIB1(KV)
       VCHAR=VCHAR1(KV)
       CNVRT=CNVRT1(KV)
       IF(KVARB.EQ.VGRIB)GO TO 200
    END DO
    GO TO 800

! check if 3d level is present in selection table

200 DO KL=1,NLVL
       IF(LEVEL.EQ.SIGL(KL))GO TO 300
    END DO
!   if all tests fail go and read next grib record
!   WRITE(*,*)'Level match not found: ',KREC,VCHAR,LEVEL 
    GO TO 800

! load the entire grib data record into the buffer

300 KREC=KREC+1
!   direct IO subroutines
    CALL FCPTPS(HANDLE,KBYTE,*900)
    CALL FCREAD(HANDLE,BUFF,1,KLEN,*900)

!   call the nmc grib unpacker
    CALL W3FI63(BUFF,KPDS,KGDS,KBMS,SVAR,KPTR,KRET)
    IF(KRET.NE.0)THEN
       WRITE(*,*)'Error W3FI63: ',KRET
       STOP
    END IF

!   century fix
    KPDS(8)=MOD(KPDS(8),100)

!   set the current time
    IYR=KPDS(8)
    IMO=KPDS(9)
    IDA=KPDS(10)
    IHR=KPDS(11)
    IMN=KPDS(12)
    IFH=MAX(KPDS(14),KPDS(15))
    CALL TMPLUS(IYR,IMO,IDA,IHR,IFH)

!   for the the first record create an index record for pakrec
    IF(KREC.EQ.1)THEN
!      write current output time to special file
       OPEN(40,FILE='AVNTIME')
       WRITE(40,'(5I2.2)')KPDS(8),KPDS(9),KPDS(10),KPDS(11),IFH
       CLOSE (40)
    END IF

!   remap input data from one- to two-dimensional array
    RVAR=RESHAPE(SVAR,SHAPE)

!   interpolate from lat/lon to conformal grid
    CALL REGRID(GLOBAL,RVAR,NLON,NLAT,XVAR,NXP,NYP,TLAT,TLON,    &
                 CLAT1,CLON1,DLON,DLAT)
    IF(CNVRT.NE.1.0) XVAR=XVAR*CNVRT

!   then pack into ARL format and continue
!   WRITE(*,'(1X,3I5,1X,A,1X,2I5)')KREC,KV,KL,VCHAR,LTYPE,LEVEL
    CALL PAKREC(20,XVAR,CVAR,NXP,NYP,NXY,VCHAR,                        &
                IYR,IMO,IDA,IHR,IMN,IFH,(KL+1),1)

800 KBYTE=KBYTE+KLEN
    GO TO 100
900 CONTINUE

! compute number of records per time period in output file
  NREC=1
  DO K=1,NLVL
     NREC=NREC+NVAR(K)
  END DO

! rotate vector variables from true to grid orientation
  IF(.NOT.GLOBAL)CALL ROTATE(NREC,NXP,NYP)

! close out time period and write index record
  CALL PAKNDX(20)
  CLOSE (20)

  END SUBROUTINE xtract


!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MAKNDX           CREATE THE INDEX RECORD CONFIGURATION FILE
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            CREATES THE CONFIGURATION FILE FOR THE OUTPUT DATA WHICH
!            DEFINES THE GRID SYSTEM AS WELL AS ALL VARIABLES THAT ARE
!            TO BE WRITTEN TO THE OUTPUT FILE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 03 Feb 1998 (RRD) - mercator fix, added lambert
!                 15 Nov 2001 (RRD) - dynamic array allocation
!                 10 Jul 2002 (BJS) - add fields, levels
!
! USAGE:  CALL MAKNDX(N2DV,NVAR,VCHAR0,VCHAR1,LEVEL,NXP,NYP,NZP,CLAT,CLON, 
!                     KGRID,GLOBAL,DLAT,DLON,NLAT,NLON)
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            NONE
!   OUTPUT FILES:           UNIT 30 CFG_AVN output file structure
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE MAKNDX (N2DV,NVAR,VCHAR0,VCHAR1,LEVEL,NXP,NYP,NZP,CLAT,CLON, &
                   KGRID,GLOBAL,DLAT,DLON,NLAT,NLON)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)   :: n2dv          ! number sfc variables 
  INTEGER,      INTENT(IN)   :: nvar   (:)    ! number of 3d variables per level
  CHARACTER(4), INTENT(IN)   :: vchar0 (:)    ! sfc variable id
  CHARACTER(4), INTENT(IN)   :: vchar1 (:)    ! 3d variable id
  INTEGER,      INTENT(IN)   :: level  (:)    ! level information
  INTEGER,      INTENT(IN)   :: nxp           ! x dimension
  INTEGER,      INTENT(IN)   :: nyp           ! y dimension
  INTEGER,      INTENT(IN)   :: nzp           ! z dimension
  REAL,         INTENT(IN)   :: clat          ! center grid lat
  REAL,         INTENT(IN)   :: clon          ! center grid lon
  INTEGER,      INTENT(IN)   :: kgrid         ! grid type 
  LOGICAL,      INTENT(IN)   :: global        ! global lat/lon flag
  REAL,         INTENT(IN)   :: dlat          ! lat spacing 
  REAL,         INTENT(IN)   :: dlon          ! lon spacing
  INTEGER,      INTENT(IN)   :: nlat          ! number of lats
  INTEGER,      INTENT(IN)   :: nlon          ! number of lons 
  
  CHARACTER(4)  :: VCHAR(50) ! variable id
  CHARACTER(4)  :: GTYPE     ! grid identification
  CHARACTER(20) :: LABEL(18) ! optional field label

  INTEGER       :: I,N,NL,MVAR  
  REAL          :: SIG, GRIDKM 
  REAL          :: GRIDS(12), PARMAP(9)

  COMMON / SETUP / GRIDS, PARMAP

! optional field label string
  DATA LABEL/'Model Type:','Grid Numb:','Vert Coord:','Pole Lat:',         &
    'Pole Lon:','Ref Lat:','Ref Lon:','Grid Size:','Orientation:',         &
    'Cone Angle:','Sync X Pt:','Sync Y Pt:','Sync Lat:','Sync Lon:',       &
    'Reserved:','Numb X pt:','Numb Y pt:','Numb Levels:'/

! reserve space for output data arrays
  IF(KGRID.EQ.0)THEN
     GRIDKM=100.0
     GTYPE='AVNX'
  ELSEIF(KGRID.EQ.1)THEN
     GRIDKM=95.25
     GTYPE='AVNN'
  ELSEIF(KGRID.EQ.2)THEN
     GRIDKM=95.25
     GTYPE='AVNS'
  ELSEIF(KGRID.EQ.3)THEN
     GRIDKM=0.0   
     GTYPE='AVNG'
  END IF

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
     GRIDS(10)=-90.0
     GRIDS(11)=  0.0    
!    latlon grid not global
     IF(NYP.LT.NLAT)GRIDS(10)=CLAT-NYP*DLAT/2.0
     IF(NXP.LT.NLON)GRIDS(11)=CLON-NXP*DLON/2.0
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
  OPEN(30,FILE='CFG_AVN')

! default grid number 99 (field not used)
  WRITE(30,'(A20,A4)')LABEL(1),GTYPE 
  WRITE(30,'(A20,A4)') LABEL(2),'  99'

! coordinate (1:sigma 2:pressure 3:terrain 4:hybrid)
  WRITE(30,'(A20,A4)') LABEL(3), '   2'

! grid geolocation parameters and projection
  DO I=1,12
     WRITE(30,'(A20,F10.2)')LABEL(I+3),GRIDS(I)
  END DO

! grid dimensions
  WRITE(30,'(A20,I4)') LABEL(16), NXP
  WRITE(30,'(A20,I4)') LABEL(17), NYP
  WRITE(30,'(A20,I4)') LABEL(18),(NZP+1)   

! upper level information
  DO NL=0,NZP   

     WRITE(LABEL(1),'(A6,I2,A1)')'Level ',NL,':'
     IF(NL.EQ.0)THEN
        SIG=0.0     
        MVAR=N2DV
        VCHAR=VCHAR0
     ELSE
        SIG=LEVEL(NL)
        MVAR=NVAR(NL)
        VCHAR=VCHAR1
     END IF

!    assumes no more than 15 fields per level
     IF(SIG.LT.1.0)THEN
        WRITE(30,'(A20,F6.5,I3,15(1X,A4))')LABEL(1),   &
        SIG,MVAR,(VCHAR(N),N=1,MVAR)
     ELSEIF(SIG.GE.1.AND.SIG.LT.10.0)THEN
        WRITE(30,'(A20,F6.4,I3,15(1X,A4))')LABEL(1),   &
        SIG,MVAR,(VCHAR(N),N=1,MVAR)
     ELSEIF(SIG.GE.10.AND.SIG.LT.100.0)THEN
        WRITE(30,'(A20,F6.3,I3,15(1X,A4))')LABEL(1),   &
        SIG,MVAR,(VCHAR(N),N=1,MVAR)
     ELSEIF(SIG.GE.100.AND.SIG.LT.1000.0)THEN
        WRITE(30,'(A20,F6.2,I3,15(1X,A4))')LABEL(1),   &
        SIG,MVAR,(VCHAR(N),N=1,MVAR)
     ELSEIF(SIG.GE.1000)THEN
        WRITE(30,'(A20,F6.1,I3,15(1X,A4))')LABEL(1),   &
        SIG,MVAR,(VCHAR(N),N=1,MVAR)
     END IF

  END DO
  CLOSE (30) 

END SUBROUTINE makndx



!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  MKGRID           DETERMINE THE LAT/LON OF EACH GRID POINT
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            GOES THROUGH EACH NODE OF THE OUTPUT GRID AND PLACES THE
!            LAT/LON VALUE OF THAT POINT INTO AN ARRAY
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 03 Feb 1999 (RRD) - mercator definition patch
!                 15 Nov 2001 (RRD) - dynamic array allocation 
!
! USAGE:  CALL MKGRID(NX,NY,TLAT,TLON)
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           NONE
!   OUTPUT FILES:          NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE MKGRID(NX,NY,TLAT,TLON)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NX,NY      ! output grid dimensions
  REAL,    INTENT(OUT) :: TLAT(:,:)  ! position of nodes
  REAL,    INTENT(OUT) :: TLON(:,:)  ! position of nodes

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



!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  REGRID           INTERPOLATES DATA TO THE OUTPUT GRID
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            FOR A GIVEN VARIABLE WILL GO THROUGH EACH NODE OF THE OUTPUT
!            GRID AND ITERPOLATE A VALUE TO THAT POINT FROM THE INPUT
!            DATA GRID.  ONLY LINEAR INTERPOLATION METHODS ARE USED.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 03 Feb 1999 (RRD) - fixed interpolation around prime
!                 15 Nov 2001 (RRD) - dynamic array allocation
!
! USAGE:  CALL REGRID(GLOBAL,I1,J1,V1,NX1,NY1,V2,NX2,NY2,TLAT,TLON, 
!                     CLAT1,CLON1,DLON,DLAT)
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           NONE
!   OUTPUT FILES:          NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE REGRID(GLOBAL,V1,NX1,NY1,V2,NX2,NY2,TLAT,TLON,    &
                  CLAT1,CLON1,DLON,DLAT)

  IMPLICIT NONE

  LOGICAL, INTENT(IN)    :: GLOBAL          ! global lat/lon grid 
  INTEGER, INTENT(IN)    :: NX1, NY1        ! input dimensions
  INTEGER, INTENT(IN)    :: NX2, NY2        ! output dimensions
  REAL,    INTENT(INOUT) :: V1(:,:)         ! input array
  REAL,    INTENT(OUT)   :: V2(:,:)         ! output array
  REAL ,   INTENT(IN)    :: TLAT(:,:)       ! position for each node 
  REAL ,   INTENT(IN)    :: TLON(:,:)       ! position for each node 
  REAL,    INTENT(IN)    :: CLAT1, CLON1    ! position of corner (1,1)
  REAL,    INTENT(IN)    :: DLAT,  DLON     ! spacing between input grid
 
  INTEGER :: i,j,ii,jj,ilo,ihi,jlo,jhi
  REAL    :: xp,yp,fxi,fyj,top,bot,temp

! just move data into output array if latlon grid 
  IF(GLOBAL)THEN
     DO J=1,NY2  
!       flip order of vertical array element
        JJ=NY2-J+1
        DO I=1,NX2   
           V2(I,J)=V1(I,JJ)        
        END DO
     END DO
     RETURN 
  END IF

! interpolate values to new grid
  DO I=1,NX2
  DO J=1,NY2

!    compute adjacent index values on grid 1
     XP=1.0+(TLON(I,J)-CLON1)/DLON
     YP=1.0+(CLAT1-TLAT(I,J))/DLAT

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



!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  ROTATE           CONVERT WINDS FROM TRUE TO GRID ORIENTATION
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            AFTER THE OUTPUT FILE HAS BEEN WRITTEN IT IS NECESSARY TO
!            GO BACK AND CONVERT THE WINDS FROM COMPASS ORIENTATION TO
!            GRID ORIENTATION.  THIS WAS NOT POSSIBLE AS THE FILE WAS
!            INITIALLY WRITTEN BECAUSE THE U,V WIND COMPONENTS WERE NOT
!            IN SEQUENTIAL ORDER IN THE INPUT DATA FILE.  BOTH COMPONENTS
!            ARE REQUIRED SIMULTANEOUSLY TO DO THE ROTATION. THIS ROUTINE
!            READS THE OUTPUT FILE, ROTATES THE WINDS, AND WRITES THOSE
!            RECORDS BACK INTO THE FILE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD)
!                 15 Nov 2001 (RRD) - dynamic array allocation
!
! USAGE:  CALL ROTATE(NREC,NXP,NYP)
!   INPUT ARGUMENT LIST:  see below
!   OUTPUT ARGUMENT LIST: none
!   INPUT FILES:          UNIT 20 - DATA.AVN direct access
!   OUTPUT FILES:         UNIT 20 - DATA.AVN direct access
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE ROTATE(NREC,NXP,NYP)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NREC       ! number of records per time period
  INTEGER, INTENT(IN) :: NXP,NYP    ! dimensions of real arrays

! data array
  CHARACTER(1), ALLOCATABLE :: CVAR(:)
  REAL,         ALLOCATABLE :: XVAR(:,:), YVAR(:,:) 

  CHARACTER(50) :: LABEL 
  CHARACTER(4)  :: VARB, VARBX, VARBY
  INTEGER       :: IYR,IMO,IDA,IHR,IFH,IMN
  INTEGER       :: NXY,KL,KG,NEXP,KREC,KSUM,LEVEL
  REAL          :: PREC,VAR1

!---------------------------------------------------------------------
  INTERFACE

  SUBROUTINE PAKINP(RVAR,CVAR,NX,NY,NX1,NY1,LX,LY,PREC,NEXP,VAR1,KSUM)
  REAL,          INTENT(OUT)   :: rvar (:,:)     ! real data unpacked
  CHARACTER(1),  INTENT(IN)    :: cvar (:)       ! packed input of NX*NY
  INTEGER,       INTENT(IN)    :: nx,ny          ! size of input array  
  INTEGER,       INTENT(IN)    :: nx1,ny1        ! optional sub-grid left edge 
  INTEGER,       INTENT(IN)    :: lx,ly          ! length of sub-grid
  REAL,          INTENT(IN)    :: prec           ! precision of packed data 
  INTEGER,       INTENT(IN)    :: nexp           ! packing scaling exponent
  REAL,          INTENT(IN)    :: var1           ! value of array at (1,1)
  INTEGER,       INTENT(INOUT) :: ksum           ! rotating checksum 
  END SUBROUTINE pakinp

  SUBROUTINE U2GRID(UU,VV,NX,NY)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY      ! array dimensions
  REAL, INTENT(INOUT) :: UU(:,:)     ! u component
  REAL, INTENT(INOUT) :: VV(:,:)     ! v component
  END SUBROUTINE u2grid

  END INTERFACE
!---------------------------------------------------------------------

  IMN=0
  NXY=NXP*NYP
  ALLOCATE (CVAR(NXP*NYP))
  ALLOCATE (XVAR(NXP,NYP),YVAR(NXP,NYP))

  100 FORMAT(7I2,A4,I4,2E14.7)

  KREC=1
  DO WHILE (KREC.LT.NREC)
     READ(20,REC=KREC,ERR=900)LABEL
     READ(LABEL,'(14X,A4)')VARB

     IF(VARB(1:1).EQ.'U')THEN
!       first record of record pair gives u-component variable
        READ(20,REC=KREC)LABEL,CVAR
        READ(LABEL,100)IYR,IMO,IDA,IHR,IFH,KL,KG,VARBX,NEXP,PREC,VAR1
        KSUM=-1
        CALL PAKINP(XVAR,CVAR,NXP,NYP,1,1,NXP,NYP,PREC,NEXP,VAR1,KSUM)

!       second record in pair gives v-component variable
        KREC=KREC+1
        READ(20,REC=KREC)LABEL,CVAR
        READ(LABEL,100)IYR,IMO,IDA,IHR,IFH,KL,KG,VARBY,NEXP,PREC,VAR1
        KSUM=-1
        CALL PAKINP(YVAR,CVAR,NXP,NYP,1,1,NXP,NYP,PREC,NEXP,VAR1,KSUM)

        IF(VARBY(1:1).EQ.'V')THEN
!          both components available then rotate
           LEVEL=KL+1
           CALL U2GRID(XVAR,YVAR,NXP,NYP)

!          rewrite those records into the data file
           CALL PAKREC(20,XVAR,CVAR,NXP,NYP,NXY,VARBX,                     &
                       IYR,IMO,IDA,IHR,IMN,IFH,LEVEL,0)
           CALL PAKREC(20,YVAR,CVAR,NXP,NYP,NXY,VARBY,                     &
                       IYR,IMO,IDA,IHR,IMN,IFH,LEVEL,0)
!          WRITE(*,*)'Rotate: ',VARBX,VARBY,LEVEL
        ELSE
           WRITE(*,*)'Error ROTATE: record pair not U,V'
           STOP
        END IF
     END IF
     KREC=KREC+1

  END DO
  DEALLOCATE (CVAR,XVAR,YVAR)
  WRITE(*,*)'Wind components rotated to grid alignment' 
  RETURN

  900 WRITE(*,*)'Error ROTATE: reading input data'
  STOP

END SUBROUTINE rotate



!##############################################################################
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  U2GRID           CONVERT WINDS FROM TRUE TO GRID ORIENTATION
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!            CONVERTS THE COMPONENT WINDS FROM THE LAT/LON GRID RELATIVE
!            TO NORTH TO COMPONENTS RELATIVE TO THE ORIENTATION OF THE
!            OUTPUT GRID AT EACH NODE OF THE OUTPUT GRID
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 18 Feb 1997 (RRD) 
!                 15 Nov 2001 (RRD) - dynamic array allocation
!
! USAGE:  CALL U2GRID(UU,VV,NX,NY)
!   INPUT ARGUMENT LIST:  see below
!   OUTPUT ARGUMENT LIST: see below
!   INPUT FILES:          NONE
!   OUTPUT FILES:         NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE U2GRID(UU,VV,NX,NY)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NX, NY      ! array dimensions
  REAL, INTENT(INOUT) :: UU(:,:)     ! u component
  REAL, INTENT(INOUT) :: VV(:,:)     ! v component

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
