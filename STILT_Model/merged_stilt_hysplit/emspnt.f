!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  EMSPNT           EMiSsion puff/particle at a PoiNT
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   EMISSION INITIALIZATION STARTS A NEW PUFF OR PARTICLE IF THE
!   CURRENT MODEL TIME IS WITHIN THE TIME LIMITS SPECIFIED FOR START
!   A POLLUTANT RELEASE.  EMISSIONS CAN BE STARTED AGAIN AT QCYCLE
!   INTERVALS.  IF MULTIPLE RELEASE LOCATIONS ARE DEFINED THEN THE
!   EMISSIONS ARE UNIFORMLY DISTRIBUTED WITHIN A LAYER FROM THE LAST
!   STARTING HEIGHT TO THE CURRENT STARTING HEIGHT WHEN TWO RELEASE
!   ARE IN THE SAME LOCATION.  OTHERWISE IT IS A POINT SOURCE.
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 14 Feb 1997 (RRD)
!                 22 Apr 1999 (RRD) - added location specific emissions
!                 14 Apr 2000 (RRD) - generalized emission start stop test
!                 04 Sep 2000 (RRD) - fortran90 upgrade
!                 21 May 2001 (RRD) - multiple process emission cycle -DMPI
!                 30 Jul 2001 (RRD) - temporal emission variation file
!                 04 Oct 2001 (RRD) - simultaneous multiple meteorology
!                 22 Mar 2002 (RRD) - consistency with pm10 emission 
!                 09 Sep 2002 (RRD) - fortran coding standards
!                 15 Dec 2003 (RRD) - space/comma delimted format on emit
!                 02 Apr 2004 (RRD) - generic file unit numbers
!                 16 Jul 2004 (RRD) - emission on/off message
!                 22 Jul 2004 (RRD) - exceeding puff limit no exit
!                 10 Aug 2004 (RRD) - enhanced particle-puff conversions
!                 29 Nov 2004 (RRD) - four digit year in efile tm call
!                 11 May 2005 (RRD) - zero initial vertical variance
!                 24 Oct 2005 (RRD) - emission duration < time step
!                 06 Jan 2006 (RRD) - backward emission fix
!                 07 Mar 2006 (RRD) - removed file I/O option
!                 17 Oct 2007 (RRD) - revised emission start/stop tests 
!                 05 Nov 2007 (RRD) - refined short-time release correction
!                 17 Jan 2008 (RRD) - num_job test on minimum release rate
!                 23 Jan 2008 (RRD) - emit denial return code
!                 07 Aug 2008 (RRD) - refined INITD test
!
! USAGE:  CALL EMSPNT(SPOT,DIRT,NLOC,NUMTYP,KPM,INITD,DT,JET,NSORT,
!              MASS,XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,
!              QCYCLE,NUMPAR,MAXPAR,job_id,num_job,ichem,kret)
!
!   INPUT ARGUMENT LIST:    see below
!   OUTPUT ARGUMENT LIST:   see below
!   INPUT FILES:            none
!   OUTPUT FILES:           none
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE EMSPNT(SPOT,DIRT,NLOC,NUMTYP,KPM,INITD,DT,JET,NSORT,MASS,       &
                  XPOS,YPOS,ZPOS,SIGH,SIGU,SIGV,SIGW,HDWP,PAGE,PTYP,PGRD,  &
                  QCYCLE,NUMPAR,MAXPAR,job_id,num_job,ichem,kret,rseed,zsg,zsg_full,nlvl)

  USE funits
  use module_defgrid ! meteorology grid and file

  IMPLICIT NONE

  INCLUDE 'DEFCONC.INC' ! pollutant and concentration grid
  INCLUDE 'DEFSPOT.INC' ! multiple source information

!-------------------------------------------------------------------------------
! argument list variables
!-------------------------------------------------------------------------------

  TYPE(rset),INTENT(INOUT) :: spot(:) ! source location characteristics
  TYPE(pset),INTENT(INOUT) :: dirt(:) ! for each pollutant type 
  INTEGER,   INTENT(IN)    :: nloc    ! total number of source locations
  INTEGER,   INTENT(IN)    :: numtyp  ! number of pollutant types
  INTEGER,   INTENT(INOUT) :: kpm     ! number of puffs or particles
  INTEGER,   INTENT(IN)    :: initd   ! initial distribution type
  REAL,      INTENT(IN)    :: dt      ! time step (min)
  INTEGER,   INTENT(IN)    :: jet     ! current elapsed time (min)
  REAL,      INTENT(IN)    :: qcycle  ! optional emission cycle time in hours
  INTEGER,   INTENT(IN)    :: numpar  ! maximum number of particles permitted

  INTEGER,   INTENT(INOUT) :: nsort (:)   ! index of sorted elements
  REAL,      INTENT(INOUT) :: mass  (:,:) ! mass of pollutant (arbitrary units)
  REAL,      INTENT(INOUT) :: xpos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: ypos  (:)   ! horizontal position (grid units)
  REAL,      INTENT(INOUT) :: zpos  (:)   ! puff center height (sigma)
  REAL,      INTENT(INOUT) :: sigh (:)    ! horizontal puff sigma 
  REAL,      INTENT(INOUT) :: sigu (:)    ! turbulence u'2    
  REAL,      INTENT(INOUT) :: sigv (:)    ! turbulence v'2
  REAL,      INTENT(INOUT) :: sigw (:)    ! turbulence w'2 or vertical puff sigma
  INTEGER,   INTENT(INOUT) :: hdwp  (:)   ! Horizontal distribution pollutant
  INTEGER,   INTENT(INOUT) :: page  (:)   ! pollutant age since release (min)
  INTEGER,   INTENT(INOUT) :: ptyp  (:)   ! pollutant type index number
  INTEGER,   INTENT(INOUT) :: pgrd  (:)   ! meteorological grid of puff position
  INTEGER,   INTENT(IN)    :: maxpar      ! maximum particle number
  INTEGER,   INTENT(IN)    :: job_id      ! mpi implementation 
  INTEGER,   INTENT(IN)    :: num_job     ! mpi implementation
  INTEGER,   INTENT(IN)    :: ichem       ! chemistry options index
  INTEGER,   INTENT(OUT)   :: kret        ! emit denial return code
  INTEGER,   INTENT(IN)    :: rseed       ! random number seed
  INTEGER,   INTENT(IN)    :: nlvl        ! number of vertical levels
  real,      intent(in)    :: zsg(nlvl)   ! model vertical layer values
  real,      intent(in)    :: zsg_full(nlvl)   ! model vertical layer values

!-------------------------------------------------------------------------------
! internal variables
!-------------------------------------------------------------------------------

  REAL, PARAMETER  :: sigr   = 1.54      ! top-hat radius
  LOGICAL, allocatable          :: emit(:)   ! flag current emissions
  LOGICAL, allocatable          :: pmit(:), print_pmit(:)   ! previous emission flags for controlling emissions, diagnostics
  REAL, PARAMETER  :: PI     = 3.14159265358979

  INTEGER          :: ksb,ksd,initk,nemit
  REAL             :: qsumall,qtot,zpl,delz,qval,qstep
  INTEGER          :: np,np1,np2,kcycl,n,nrpts,kk,npar,ktp,kp,nphr,maxdim
  real             :: ran3
  real             :: xp1, xp2, yp1, yp2, zp1, zp2, dzfull
  integer          :: kklvl, klvl, kt, npar_kk

  REAL,ALLOCATABLE :: tfact(:)

  SAVE PMIT, print_pmit

!-------------------------------------------------------------------------------
! check to determine if emissions required
!-------------------------------------------------------------------------------

  IF(.NOT.ALLOCATED(tfact)) &
     ALLOCATE (tfact(numtyp))
  IF(.NOT.ALLOCATED(emit)) &
     ALLOCATE (emit(numtyp))
  IF(.NOT.ALLOCATED(pmit)) then
     ALLOCATE (pmit(numtyp))
     pmit(:) = .FALSE.
  end IF
  IF(.NOT.ALLOCATED(print_pmit)) then
     ALLOCATE (print_pmit(numtyp))
     print_pmit(:) = .FALSE.
  end IF
  
  KRET=0
  TFACT=1.0
  EMIT(:)=.FALSE.
  nemit=0
  DO KK=1,NUMTYP
     IF(INT(DT).GT.0)THEN
        KSB=JET-DIRT(KK)%START%MACC
        KSD=MAX(INT(DIRT(KK)%QHRS*60.0),INT(DT))
     ELSE
        KSB=DIRT(KK)%START%MACC-JET
        KSD=INT(MAX(ABS(DIRT(KK)%QHRS*60.0),ABS(DT)))
     END IF
     ! disable instantaneous emissions if they were emitted already
     ! (possible if abs(dt) in previous call was smaller than in current call,
     !  in which case ksb < ksd for both calls)
     if (pmit(kk) .and. ABS(DIRT(KK)%QHRS*60.0) .lt. abs(dt)) ksd=0
     IF(KSB.GE.0.AND.KSB.LT.KSD) THEN 
        nemit=nemit+1
        EMIT(kk)=.TRUE.
!       adjust for short emissions relative to the time step:
!       tfact is the fraction of the time step that is contained within the release window
!       Only do this for non-instantaneous emissions (ABS(DIRT(KK)%QHRS*60.0) .gt. abs(dt))
        IF(ABS(DIRT(KK)%QHRS*60.0) .gt. abs(dt) .and. KSD-KSB .LT. ABS(INT(DT))) &
             TFACT(KK)=(KSD-KSB)/ABS(DT)
     END IF
  END DO

  IF(any(EMIT.AND..NOT.PMIT)) WRITE(KF21,*)' NOTICE emspnt: emissions started, jet, kpm, nemit=',jet,kpm,nemit
  IF(any(.NOT.EMIT.AND.print_PMIT)) WRITE(KF21,*)' NOTICE emspnt: emissions terminated, jet, kpm=',jet,kpm
  do kk=1,numtyp
     ! reset previously emitted flag:
     print_PMIT(kk)=EMIT(kk)
     ! for non-instantaneous receptors, do not ever reset to false
     if (ABS(DIRT(KK)%QHRS*60.0) .lt. abs(dt)) then
        PMIT(kk)=PMIT(kk) .or. EMIT(kk)
     else
        PMIT(kk)=EMIT(kk)
     end if
  end do
  
!*********************************************
! CHG(11/04/02) don't emit particles when not wanted!
! dwen(20090317):change according to EMSPNT of STILT
      IF(all(.NOT. EMIT) .OR. NUMPAR .EQ. 0)RETURN
!  IF(.NOT.EMIT)RETURN
!*********************************************

!-------------------------------------------------------------------------------
! determine type of emission and number of locations
!-------------------------------------------------------------------------------

! multiple starting locations can either be point sources
! if each starting point is in a different location or vertical
! line sources if two starting points are in the same location.
! The line source is distributed between the given heights

  IF(NLOC.GT.1)THEN
     SPOT(1)%ZV=SPOT(1)%ZP
     DO N=2,NLOC
!****************************************
! JCL:(5/5/00)don't want weird line source
!dwen(20090317):comment out the following lines according to EMSPNT of STILT
!        IF(INT(SPOT(N)%XP*10000.0).EQ.                                     &
!           INT(SPOT(N-1)%XP*10000.0).AND.                                  &
!           INT(SPOT(N)%YP*10000.0).EQ.                                     &
!           INT(SPOT(N-1)%YP*10000.0))THEN
!
!!          when position the same move previous point release
!!          height into ZV (line source defined as ZV->ZP)
!!          then only emit at locations with ZV<>0
!           SPOT(N)%ZV=SPOT(N-1)%ZP
!           SPOT(N-1)%ZV=0.0
!        ELSE
!!          point source bottom equals release height
           SPOT(N)%ZV=SPOT(N)%ZP
!        END IF
!*************************************************
     END DO

!    count up the number of different x,y release locations with
     NRPTS=0
     DO N=1,NLOC
!       those set to zero are used to determine line source heights
        IF(SPOT(N)%ZV.NE.0.0)NRPTS=NRPTS+1
     END DO

  ELSE
     NRPTS=1
     SPOT(1)%ZV=SPOT(1)%ZP
  END IF

!-------------------------------------------------------------------------------
! set the number of units to emit
!-------------------------------------------------------------------------------

  INITK=INITD
  IF(INITK.GE.100) INITK=MOD(INITD/10,10)

  IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
!    gaussian or top-hat emissions
     NPAR=1
  ELSE
     QSUMall=0.0
     do n=1,nloc
!    total emissions - pollutant hours is summed because each
!    pollutant type is emitted as an independent particle
        if (spot(n)%ptyp .eq. 0 .and. ichem .ne. 1) then
           ! default: all pollutants released from all locations
           np1=1
           np2=numtyp
        else
           ! only one pollutant released from this location:
           np1 = spot(n)%ptyp
           if (spot(n)%ptyp .eq. 0) np1=n !ichem=1 option
           NP2=Np1
        end if
!!$        WRITE(KF21,*)' debug emspnt: n, np1, np2=',n, np1, np2
        DO KK=np1,np2
           QSUMall=QSUMall+MAX(ABS(DT/60.0),ABS(DIRT(KK)%QHRS))
        END DO
     end do
!    particle emission rate: number particles per hour per source
     NPHR=NINT(FLOAT(NUMPAR)/qsumall)
!    particle emissions per time step (requires at least one per job)
     NPAR=MAX(num_job, INT(FLOAT(NPHR)*ABS(DT/60.0)))
!!$     WRITE(KF21,*)' debug emspnt: numpar, qsumall, nphr, npar=',numpar, qsumall, nphr, npar
  END IF

!-------------------------------------------------------------------------------
! loop through the number of independent source locations
!-------------------------------------------------------------------------------

  MAXDIM = SIZE(mass,1)  ! number of pollutants on single particle
  NP=0                   ! internal particle counter
  kt=1                   ! time index for grid

nloop : DO N=1,NLOC

! check for source skip (vertical line source)
  IF(SPOT(N)%ZV.NE.0.0)THEN
     xp1=spot(n)%xp
     xp2=xp1
     yp1=spot(n)%yp
     yp2=yp1
     zp1=spot(n)%zp
     zp2=zp1
     if (spot(n)%dxyp .gt. 0) then
        xp1=max(2.,xp1-spot(n)%dxyp)
        xp2=min(grid(spot(n)%kg,kt)%nx-2.,xp2+spot(n)%dxyp)
        yp1=max(2.,yp1-spot(n)%dxyp)
        yp2=min(grid(spot(n)%kg,kt)%ny-2.,yp2+spot(n)%dxyp)
     endif
     if (spot(n)%dzp .gt. 0) then
        ! find vertical level:
        if (zp1 .le. 1 .and. zp1 .ge. zsg_full(1)) then
           klvl=1
           dzfull=1-zsg_full(1)
        else
           klvl=0
           do kklvl=2,nlvl
              if (zp1 .le. zsg_full(kklvl-1) .and. zp1 .ge. zsg_full(kklvl)) then
                 klvl=kklvl
                 exit
              end if
           end do
           if (klvl .eq. 0) stop 'emspnt: Unexpected SPOT%ZP > 1 or < zsg_full(nlvl)'
           dzfull=zsg_full(klvl-1)-zsg_full(klvl)
        end if
        zp1=max(zsg(nlvl),zp1-spot(n)%dzp*dzfull) !upper limit (numerically smaller)
        zp2=min(1.,zp2+spot(n)%dzp*dzfull) !lower limit (numerically larger)
     end if

! Each pollutant type will start its own trajectory unless this routine is 
! modified accordingly (MAXDIM >1). 

     if (spot(n)%ptyp .eq. 0 .and. ichem .ne. 1) then
        ! default: all pollutants released from all locations
        np1=1
        np2=numtyp
     else
        ! only one pollutant released from this location:
        np1 = spot(n)%ptyp
        if (spot(n)%ptyp .eq. 0) np1=n !ichem=1 option
        NP2=Np1
     end if

tloop : DO KK=NP1,NP2 

! all type emissions go into index 1 of particle mass array unless maxdim>1
  KTP=MIN(KK,MAXDIM)

! check if this pollutant requires a start
  IF(emit(kk))THEN       

     ! pro-rate emission rate by the fraction of the time step remaining in the release window
     ! avoids exceeding maximum allowable number of particles

     npar_kk=int(0.5+tfact(kk)*real(npar))


!    multiple emissions only for particles
ploop : DO KP=1,NPAR_kk

!       Use for simulations with multiple processors to only emit particles 
!       when process id even multiple of particle numbe. In a single process
!       environment num_job=1 and job_id=0.

        NP=NP+1
        IF(MOD(NP,num_job).NE.job_id)CYCLE ploop

        KPM=KPM+1
        IF(KPM.GT.NUMPAR)THEN
           KPM=NUMPAR
           WRITE(KF21,*)'Warning: emspnt - exceeding puff limit'
           KRET=1
           EXIT nloop
        END IF
        NSORT(KPM)=KPM

!       initial position from main program, plus optional random perturbation
        if (xp1 .eq. xp2) then
           XPOS(KPM)=SPOT(N)%XP
        else
           XPOS(KPM)=xp1 + ran3(rseed)*(xp2-xp1)
        end if
        if (yp1 .eq. yp2) then
           YPOS(KPM)=SPOT(N)%YP
        else
           YPOS(KPM)=yp1 + ran3(rseed)*(yp2-yp1)
        end if
        if (zp1 .eq. zp2) then
           ZPOS(KPM)=SPOT(N)%ZP
        else
           ZPOS(KPM)=zp1 + ran3(rseed)*(zp2-zp1)
        end if
!!$        write (KFJCLMSG,'(a,2i8,9g15.6)') 'emspnt: ptyp, kpm, xpos, ypos, zpos, xp1, xp2, yp1, yp2, zp1, zp2=',&
!!$             kk, kpm, xpos(kpm), ypos(kpm), zpos(kpm), xp1, xp2, yp1, yp2, zp1, zp2
        
!       horizontal variances all start at zero
        SIGU(KPM)=0.0
        SIGV(KPM)=0.0
        SIGW(KPM)=0.0

        IF(SPOT(N)%AREA.LE.0.0)THEN
!          points source defined for all species
           SIGH(KPM)=0.0
        ELSE
!          defined by source compute sigma for uniform radius
           SIGH(KPM)=SQRT(SPOT(N)%AREA/PI)/SIGR
        END IF

!       multiple locations get initial vertical distribution
!       when the first two starting points at the same x,y position
        DELZ=ABS(SPOT(N)%ZV-SPOT(N)%ZP)
        IF(DELZ.GT.0.0)THEN
           ZPL=MAX(SPOT(N)%ZV,SPOT(N)%ZP)
           IF(INITK.EQ.1.OR.INITK.EQ.2)THEN
!              puff variance set to layer depth
               SIGW(KPM)=ABS(DELZ/SIGR/2.0)
               ZPOS(KPM)=ZPL-DELZ/2.0
            ELSE
!              vertical variance set to zero
               SIGW(KPM)=0.0
!              particles get distributed in the layer
               ZPOS(KPM)=ZPL-FLOAT(KP)*DELZ/FLOAT(NPAR)
            END IF
         END IF

!        initial distribution (see main for definitions)
         HDWP(KPM)=INITD
!        initial age at zero
         PAGE(KPM)=0
!        pollutant type always defined from pollutant index
         PTYP(KPM)=KK
!        initial grid is the default startup grid from main
         PGRD(KPM)=SPOT(N)%KG

         IF(SPOT(N)%QTRM.EQ.0.0)THEN
!           emission defined by species
            QTOT=ABS(DIRT(KK)%QHRS)*DIRT(KK)%QRATE
         ELSE
!           emission defined by source location
            QTOT=ABS(DIRT(KK)%QHRS)*SPOT(N)%QTRM
         END IF

!        number of time steps in emission period
         QSTEP=MAX(1.0,60.0*ABS(DIRT(KK)%QHRS)/ABS(DT))
!        emission per time step
         QVAL=QTOT/QSTEP
!        divide amount over the number of units emitted
!         MASS(KTP,KPM)=TFACT(KK)*QVAL/NPAR (adjust npar_kk instead - see ploop above)
         MASS(KTP,KPM)=QVAL/NPAR

!     particle loop
      END DO ploop

!  start time test
   END IF

! pollutant type loop
  END DO tloop

! vertical line source skip
  END IF

! number of sources loop
  END DO nloop

!-------------------------------------------------------------------------------
! check for emission cycling
!-------------------------------------------------------------------------------

  DO KK=1,NUMTYP
!    test for end of emission cycle
     IF(INT(DT).GT.0)THEN
        KSB=JET+INT(DT)-DIRT(KK)%START%MACC
        KSD=MAX(INT(DIRT(KK)%QHRS*60.0),INT(DT))
     ELSE
        KSB=DIRT(KK)%START%MACC-JET-INT(DT)
        KSD=INT(MAX(ABS(DIRT(KK)%QHRS*60.0),ABS(DT)))
     END IF

     IF(KSB.GE.KSD)THEN        
!       optional restart of emissions at some later time
        KCYCL=NINT(SIGN(QCYCLE*60.0,DIRT(KK)%QHRS))
        DIRT(KK)%START%MACC=DIRT(KK)%START%MACC+KCYCL
     END IF
  END DO

END SUBROUTINE emspnt
