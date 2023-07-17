module unpack_subs

  private

  public:: unpack, idate_timdiff, id2pos

contains

  SUBROUTINE UNPACK(label,icx,mn,CPACK,RDATA,NX,NY,NEXP,PREC,VAR1,i1,i2,di,j1,j2,dj, &
       &            use_16,iprint)

    implicit none

    CHARACTER(52),INTENT(IN)  :: label
    CHARACTER(1),INTENT(IN)  :: CPACK(:)  
    INTEGER,     INTENT(IN)  :: NX,NY,NEXP,i1,i2,di,j1,j2,dj,icx,mn,iprint
    REAL,        INTENT(OUT) :: RDATA(NX,NY)   
    REAL,        INTENT(IN)  :: PREC,VAR1
    logical, intent (in) :: use_16

    real ::rmin, rmax, rmean

    real :: scale, vold, rmean_i, rmin_i, rmax_i
    integer :: indx, tstj, j, tsti, i, itmp, itmp1, itmp2
    integer :: k18

    ! only required when dealing with SUN F90 compiler
    ! replace ICHAR below with internally defined JCHAR function
    ! CHARACTER MYCHR*1
    ! JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

    if (.not. use_16) then
       SCALE=2.0**(7-NEXP)
    else
       scale=2.0**(15-NEXP)
    end if

    VOLD=VAR1
    INDX=0
    tstj=j1
    rmean = 0.
    rmean_i = 0.
    rmin_i = +1.e10
    rmax_i = -1.e10
    DO J=1,NY
       tsti=i1
       DO I=1,NX
          INDX=INDX+1
          ! the following 2 lines taken from pakinp.f90:
          if (.not. use_16) then
             ITMP=ICHAR(CPACK(INDX))
             IF(ITMP.LT.0)ITMP=ITMP+256
             RDATA(I,J)=(ITMP-127.)/SCALE+VOLD
          else
             ITMP1=ICHAR(CPACK(INDX))
             INDX=INDX+1
             ITMP2=ICHAR(CPACK(INDX))
             ITMP=ITMP2*256+ITMP1
             IF(ITMP.LT.0)ITMP=ITMP+65536
             RDATA(I,J)=((ITMP-32767.0)/SCALE)+VOLD
          end if

          VOLD=RDATA(I,J)
          ! the following line taken from pakinp.f90:
          IF(ABS(RDATA(I,J)) .LT. PREC) RDATA(I,J)=0.0
          rmean = rmean + rdata(i,j)
          if (i .le. nx-1 .and. j .le. ny-1) then
             rmean_i=rmean_i+rdata(i,j)
             rmin_i=min(rmin_i,rdata(i,j))
             rmax_i=max(rmax_i,rdata(i,j))
          endif
          IF(I.eq.tsti.AND.J.eq.tstj) then
             if (iprint .gt. 2) WRITE(*,'(3I5,E15.6)')J,I,ICHAR(CPACK(INDX)),RDATA(I,J)
             tsti = tsti + di
             if (tsti .gt. i2) tsti=-99
          endif
       END DO
       IF(J.eq.tstj) then
          tstj=tstj+dj
          if (tstj .gt. j2) tstj=-99
       endif
       VOLD=RDATA(1,J)
    END DO
    rmin = minval(rdata)
    rmax = maxval(rdata)
    rmean = rmean/(nx*ny)
    rmean_i = rmean_i/((nx-1.)*(ny-1.))
    k18 = 18
    if (label(11:11) .eq. 'X') k18=20
    if (iprint .gt. 2) WRITE(*,'(a,2i4,a,6g15.6)') &
         label(1:k18),icx,mn,': min,max,mean=',rmin,rmax,rmean, &
         rmin_i,rmax_i,rmean_i

    return

  END SUBROUTINE unpack

  subroutine idate_timdiff(ref_date,dates,diffs)
    implicit none
    integer, intent(in) :: ref_date(:), dates(:,:)
    integer, intent(out) :: diffs(:)

    integer :: datlen, ndates, idate
    integer :: iyr,imo,ida,ihr,imn, isec, riyr,rimo,rida,rihr,rimn,risec
    character (len=19) :: ndate, odate
    character (len=80) :: fmt='(i4,5(1x,i2))'

    datlen = size(ref_date)
    if (datlen .lt. 3 .or. datlen .ne. size(dates,1)) &
         stop 'invalid datlen (ref_date dim and/or dates dim)'
    ndates = min(size(dates,2),size(diffs))

    rihr=0
    rimn=0
    risec=0
    riyr=ref_date(1)
    if (riyr .lt. 100) then
       if (riyr .lt. 50) then
          riyr=riyr+2000
       else
          riyr=riyr+1900
       end if
    end if
    rimo=ref_date(2)
    rida=ref_date(3)
    if (datlen .ge. 4) rihr=ref_date(4)
    if (datlen .ge. 5) rimn=ref_date(5)
    if (datlen .ge. 6) risec=ref_date(6)

    ihr=0
    imn=0
    isec=0.

    write (odate,fmt) riyr,rimo,rida,rihr,rimn,risec
    do idate=1,ndates
       iyr=dates(1,idate)
       if (iyr .lt. 100) then
          if (iyr .lt. 50) then
             iyr=iyr+2000
          else
             iyr=iyr+1900
          end if
       end if
       imo=dates(2,idate)
       ida=dates(3,idate)
       if (datlen .ge. 4) ihr=dates(4,idate)
       if (datlen .ge. 5) imn=dates(5,idate)
       if (datlen .ge. 6) isec=dates(6,idate)
       write (ndate,fmt) iyr,imo,ida,ihr,imn,isec
       
       call geth_idts(ndate,odate,diffs(idate))
    end do
    return
  end subroutine idate_timdiff

  SUBROUTINE geth_idts (ndate, odate, idts)

    ! adapted from geth_idts in WPSV3.4.1/metgrid/src/module_date_pack.F

    IMPLICIT NONE

    !  From 2 input mdates ('YYYY-MM-DD HH:MM:SS.ffff'), 
    !  compute the time difference.

    !  on entry     -  ndate  -  the new hdate.
    !                  odate  -  the old hdate.

    !  on exit      -  idts    -  the change in time in seconds.

    CHARACTER (LEN=*) , INTENT(INOUT) :: ndate, odate
    INTEGER           , INTENT(OUT)   :: idts

    !  Local Variables

    !  yrnew    -  indicates the year associated with "ndate"
    !  yrold    -  indicates the year associated with "odate"
    !  monew    -  indicates the month associated with "ndate"
    !  moold    -  indicates the month associated with "odate"
    !  dynew    -  indicates the day associated with "ndate"
    !  dyold    -  indicates the day associated with "odate"
    !  hrnew    -  indicates the hour associated with "ndate"
    !  hrold    -  indicates the hour associated with "odate"
    !  minew    -  indicates the minute associated with "ndate"
    !  miold    -  indicates the minute associated with "odate"
    !  scnew    -  indicates the second associated with "ndate"
    !  scold    -  indicates the second associated with "odate"
    !  i        -  loop counter
    !  mday     -  a list assigning the number of days in each month

    CHARACTER (LEN=24) :: tdate
    INTEGER :: olen, nlen
    INTEGER :: yrnew, monew, dynew, hrnew, minew, scnew
    INTEGER :: yrold, moold, dyold, hrold, miold, scold
    INTEGER :: mday(12), i, newdys, olddys
    LOGICAL :: npass, opass
    INTEGER :: isign

    IF (odate.GT.ndate) THEN
       isign = -1
       tdate=ndate
       ndate=odate
       odate=tdate
    ELSE
       isign = 1
    END IF

    !  Assign the number of days in a months

    mday( 1) = 31
    mday( 2) = 28
    mday( 3) = 31
    mday( 4) = 30
    mday( 5) = 31
    mday( 6) = 30
    mday( 7) = 31
    mday( 8) = 31
    mday( 9) = 30
    mday(10) = 31
    mday(11) = 30
    mday(12) = 31

    !  Break down old hdate into parts

    hrold = 0
    miold = 0
    scold = 0
    olen = LEN(odate)

    READ(odate(1:4),  '(I4)') yrold
    READ(odate(6:7),  '(I2)') moold
    READ(odate(9:10), '(I2)') dyold
    IF (olen.GE.13) THEN
       READ(odate(12:13),'(I2)') hrold
       IF (olen.GE.16) THEN
          READ(odate(15:16),'(I2)') miold
          IF (olen.GE.19) THEN
             READ(odate(18:19),'(I2)') scold
          END IF
       END IF
    END IF

    !  Break down new hdate into parts

    hrnew = 0
    minew = 0
    scnew = 0
    nlen = LEN(ndate)

    READ(ndate(1:4),  '(I4)') yrnew
    READ(ndate(6:7),  '(I2)') monew
    READ(ndate(9:10), '(I2)') dynew
    IF (nlen.GE.13) THEN
       READ(ndate(12:13),'(I2)') hrnew
       IF (nlen.GE.16) THEN
          READ(ndate(15:16),'(I2)') minew
          IF (nlen.GE.19) THEN
             READ(ndate(18:19),'(I2)') scnew
          END IF
       END IF
    END IF

    !  Check that the dates make sense.

    npass = .true.
    opass = .true.

    !  Check that the month of NDATE makes sense.

    IF ((monew.GT.12).or.(monew.LT.1)) THEN
       PRINT*, 'GETH_IDTS:  Month of NDATE = ', monew
       npass = .false.
    END IF

    !  Check that the month of ODATE makes sense.

    IF ((moold.GT.12).or.(moold.LT.1)) THEN
       PRINT*, 'GETH_IDTS:  Month of ODATE = ', moold
       opass = .false.
    END IF

    !  Check that the day of NDATE makes sense.

    IF (monew.ne.2) THEN
       ! ...... For all months but February
       IF ((dynew.GT.mday(monew)).or.(dynew.LT.1)) THEN
          PRINT*, 'GETH_IDTS:  Day of NDATE = ', dynew
          npass = .false.
       END IF
    ELSE IF (monew.eq.2) THEN
       ! ...... For February
       IF ((dynew.GT.nfeb(yrnew)).OR.(dynew.LT.1)) THEN
          PRINT*, 'GETH_IDTS:  Day of NDATE = ', dynew
          npass = .false.
       END IF
    END IF

    !  Check that the day of ODATE makes sense.

    IF (moold.ne.2) THEN
       ! ...... For all months but February
       IF ((dyold.GT.mday(moold)).or.(dyold.LT.1)) THEN
          PRINT*, 'GETH_IDTS:  Day of ODATE = ', dyold
          opass = .false.
       END IF
    ELSE IF (moold.eq.2) THEN
       ! ....... For February
       IF ((dyold.GT.nfeb(yrold)).or.(dyold.LT.1)) THEN
          PRINT*, 'GETH_IDTS:  Day of ODATE = ', dyold
          opass = .false.
       END IF
    END IF

    !  Check that the hour of NDATE makes sense.

    IF ((hrnew.GT.23).or.(hrnew.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  Hour of NDATE = ', hrnew
       npass = .false.
    END IF

    !  Check that the hour of ODATE makes sense.

    IF ((hrold.GT.23).or.(hrold.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  Hour of ODATE = ', hrold
       opass = .false.
    END IF

    !  Check that the minute of NDATE makes sense.

    IF ((minew.GT.59).or.(minew.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  Minute of NDATE = ', minew
       npass = .false.
    END IF

    !  Check that the minute of ODATE makes sense.

    IF ((miold.GT.59).or.(miold.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  Minute of ODATE = ', miold
       opass = .false.
    END IF

    !  Check that the second of NDATE makes sense.

    IF ((scnew.GT.59).or.(scnew.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  SECOND of NDATE = ', scnew
       npass = .false.
    END IF

    !  Check that the second of ODATE makes sense.

    IF ((scold.GT.59).or.(scold.LT.0)) THEN
       PRINT*, 'GETH_IDTS:  Second of ODATE = ', scold
       opass = .false.
    END IF

    IF (.not. npass) THEN
       print*, 'Screwy NDATE: ',ndate(1:nlen)
       stop 'Bad date passed in to idate_timdiff'
    END IF

    IF (.not. opass) THEN
       print*, 'Screwy ODATE:',odate(1:olen)
       stop 'Bad date passed in to idate_timdiff'
    END IF

    !  Date Checks are completed.  Continue.

    !  Compute number of days from 1 January ODATE, 00:00:00 until ndate
    !  Compute number of hours from 1 January ODATE, 00:00:00 until ndate
    !  Compute number of minutes from 1 January ODATE, 00:00:00 until ndate

    newdys = 0
    DO i = yrold, yrnew - 1
       newdys = newdys + (365 + (nfeb(i)-28))
    END DO

    IF (monew .GT. 1) THEN
       mday(2) = nfeb(yrnew)
       DO i = 1, monew - 1
          newdys = newdys + mday(i)
       END DO
       mday(2) = 28
    END IF

    newdys = newdys + dynew-1

    !  Compute number of hours from 1 January ODATE, 00:00:00 until odate
    !  Compute number of minutes from 1 January ODATE, 00:00:00 until odate

    olddys = 0

    IF (moold .GT. 1) THEN
       mday(2) = nfeb(yrold)
       DO i = 1, moold - 1
          olddys = olddys + mday(i)
       END DO
       mday(2) = 28
    END IF

    olddys = olddys + dyold-1

    !  Determine the time difference in seconds

    idts = (newdys - olddys) * 86400
    idts = idts + (hrnew - hrold) * 3600
    idts = idts + (minew - miold) * 60
    idts = idts + (scnew - scold)

    IF (isign .eq. -1) THEN
       tdate=ndate
       ndate=odate
       odate=tdate
       idts = idts * isign
    END IF

  END SUBROUTINE geth_idts

  FUNCTION nfeb ( year ) RESULT (num_days)

    ! Compute the number of days in February for the given year

    IMPLICIT NONE

    INTEGER :: year
    INTEGER :: num_days

    num_days = 28 ! By default, February has 28 days ...
    IF (MOD(year,4).eq.0) THEN  
       num_days = 29  ! But every four years, it has 29 days ...
       IF (MOD(year,100).eq.0) THEN
          num_days = 28  ! Except every 100 years, when it has 28 days ...
          IF (MOD(year,400).eq.0) THEN
             num_days = 29  ! Except every 400 years, when it has 29 days.
          END IF
       END IF
    END IF

  END FUNCTION nfeb

  subroutine id2pos(receptors,nreceptors,rec_lon,rec_lat,rec_alt,rec_idate,ierr)

    ! Decode receptor string such as "2002x08x03x10x55x45.335Sx179.884Wx00030"
    ! into its components:
    !     yr  mon  day  hr  min      lat      lon alt
    !    2002  08   03  10   55  -45.335 -179.884  30
    ! The min and alt components (and associated delimiters) are optional, 
    !  and the precision of lat/lon is variable
    ! Assumes all receptors have the same format


    implicit none
    integer, intent(in) :: nreceptors
    character (len=*), intent(in) :: receptors(nreceptors)
    real, intent(out) :: rec_lon(nreceptors), rec_lat(nreceptors)
    integer, intent(out) :: rec_idate(5,nreceptors) ! yr,mo,da,hr,mn
    integer, intent(out) :: rec_alt(nreceptors) ! height AGL
    integer, intent(out) :: ierr ! >0: error; <0: warning; 0: ok

    character (len=1) :: delim='x',east='E',west='W',south='S',north='N' !
    integer, parameter :: min_delim=5,max_delim=8 ! Allow for missing min/alt, extra delimiter after agl
    integer :: missval=-999 ! Missing value for rec_alt in case of reading errors
    integer :: i, i1, i2, idelim(max_delim), i_southnorth,i_eastwest, ierr_io, irec
    integer :: ndelim, ndelim_date, nlen, i1_alt, i2_alt

    ! Find delimiters and east_west and south_north positions

    nlen=len_trim(receptors(1))
    ndelim=0
    ndelim_date=0
    idelim(:)=0
    i_southnorth=nlen+1
    i_eastwest=nlen+1
    do i=1,nlen
       if (receptors(1)(i:i) .eq. delim) then
          ndelim=ndelim+1
          if (ndelim .gt. max_delim) exit
          idelim(ndelim)=i
          if (i .lt. i_southnorth) ndelim_date=ndelim_date+1
       elseif (receptors(1)(i:i) .eq. south .or. receptors(1)(i:i) .eq. north) then
          i_southnorth=i
       elseif (receptors(1)(i:i) .eq. east .or. receptors(1)(i:i) .eq. west) then
          i_eastwest=i
       end if
    end do

    ! check output:
    ierr = 1 ! too few delimiters
    if (ndelim .lt. min_delim) return
    ierr = 2 ! no lat/lons
    if (i_southnorth .gt. nlen .or. i_eastwest .gt. nlen) return
    ierr = 3 ! wrong number of date delimiters: 4 if no minutes, 5 if minutes
    if (ndelim_date .lt. 4 .or. ndelim_date .gt. 5) return
    ierr=4 ! inconsistent delimiter and south_north positions
    if (i_southnorth .ne. idelim(ndelim_date+1)-1) return
    ierr = 0
    if (ndelim .gt. max_delim) then
       ierr = -1 ! too many delimiters, ignoring the extra ones
       ndelim=max_delim
    end if
    i1_alt=0
    i2_alt=-1
    if (ndelim .ge. ndelim_date+2) then
       if (i_eastwest .ne. idelim(ndelim_date+2)-1) then
          ierr=4 ! inconsistent delimiter and east_west positions
          return
       end if
       ! enough delimiters so there is a final alt field
       i=ndelim_date+2
       i1_alt=idelim(i)+1
       if (ndelim .gt. i) then
          i2_alt=idelim(i+1)-1  !final delimiters after alt field
       else
          i2_alt=nlen           !no final delimiters after alt field
       end if
    end if
    if (i1_alt .le. 0 .or. i2_alt .lt. i1_alt) ierr=-2 !cannot read alt
    
    do irec=1,nreceptors
       rec_idate(5,irec)=0 !default is 0 for minutes
       i1=1
       do i=1,ndelim_date
          i2=idelim(i)-1
          read(receptors(irec)(i1:i2),*,iostat=ierr_io) rec_idate(i,irec)
          i1=idelim(i)+1
          if (ierr_io .ne. 0) then
             ierr=5 !error reading date
             return
          endif
       end do
       i2=i_southnorth-1
       read(receptors(irec)(i1:i2),*,iostat=ierr_io) rec_lat(irec)
       if (ierr_io .ne. 0) then
          ierr=6 !error reading lat
          return
       endif
       if (receptors(irec)(i_southnorth:i_southnorth) .eq. south) &
            rec_lat(irec)=-rec_lat(irec)
       i=ndelim_date+1
       i1=idelim(i)+1
       i2=i_eastwest-1
       read(receptors(irec)(i1:i2),*,iostat=ierr_io) rec_lon(irec)
       if (ierr_io .ne. 0) then
          ierr=7 !error reading lon
          return
       endif
       if (receptors(irec)(i_eastwest:i_eastwest) .eq. west) &
            rec_lon(irec)=-rec_lon(irec)
       if (i1_alt .gt. 0 .and. i2_alt .ge. i1_alt) then
          read(receptors(irec)(i1:i2),*,iostat=ierr_io) rec_alt(irec)
          if (ierr_io .ne. 0) ierr=-2 !error reading alt
       end if
    end do
    
  end subroutine id2pos
  
end module unpack_subs
  
