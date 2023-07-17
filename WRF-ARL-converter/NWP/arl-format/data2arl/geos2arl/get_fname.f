! Use emacs f90-mode:   -*- mode:f90  -*-
module get_fname_module

  implicit none
  private

  public :: get_new_timestamp, get_fname, get_time

contains

  subroutine get_fname(geos_prefix, vfname, geos_sysver, valid_time, geos_suffix, &
       minutes_unit_factor, iloop, delta_per_loop, fname, ierr)

    implicit none

    character (len=*), intent(in) :: geos_prefix, geos_sysver, vfname, valid_time, geos_suffix
    integer, intent(in) :: minutes_unit_factor, iloop, delta_per_loop
    character (len=*), intent(out) :: fname
    integer, intent(out) :: ierr

    integer :: len_date, delta_time, tavg, offset_delta
    character (len=100) :: new_time

    offset_delta=0
    if (vfname(1:4) .eq. 'tavg') then
       read(vfname(5:5),'(i1)') tavg
       offset_delta=-60*minutes_unit_factor*tavg/2
    end if
    delta_time=offset_delta + (iloop-1)*delta_per_loop
    call get_new_timestamp(valid_time, new_time, delta_time)
    if (len_trim(new_time) .lt. 8) then
       ierr=1
       fname=' '
    else
       ierr=0
       if (trim(geos_sysver) .ne. 'none') then
          fname=trim(geos_prefix) // '.' // trim(vfname) // '.' // trim(geos_sysver) // '.' &
               // trim(new_time) // '.' // trim(geos_suffix)
       else
          fname=trim(geos_prefix) // '.' // trim(vfname) // '.' &
               // trim(new_time) // '.' // trim(geos_suffix)
       end if
    end if
    return
  end subroutine get_fname

  SUBROUTINE get_new_timestamp(old_time_stamp, new_time_stamp, delta_time)
    
    IMPLICIT NONE

  ! adapted from: geth_newdate (ndate, odate, idt) (WPSV3.4.1/metgrid/src/module_date_pack.F)
      
  !  From old time_stamp and delta-time, compute the new time_stamp
   
  !  on entry     -  old_time_stamp  -  the old time stamp
  !                  delta_time    -  the change in time (in days, hour, minutes, or seconds, depending on length of time_stamp)
   
  !  on exit      -  new_time_stamp  -  the new time stamp
  ! time_stamp format: yyyymmdd_hhmnss
  !                            [optional]
  !                           8              - days
  !                             11           - hours
  !                               13         - minutes
  !                                 15       - seconds
      
    INTEGER , INTENT(IN)           :: delta_time
    CHARACTER (LEN=*) , INTENT(OUT) :: new_time_stamp
    CHARACTER (LEN=*) , INTENT(IN)  :: old_time_stamp
      
       
      !  Local Variables
       
      !  yrold    -  indicates the year associated with "old_time_stamp"
      !  moold    -  indicates the month associated with "old_time_stamp"
      !  dyold    -  indicates the day associated with "old_time_stamp"
      !  hrold    -  indicates the hour associated with "old_time_stamp"
      !  miold    -  indicates the minute associated with "old_time_stamp"
      !  scold    -  indicates the second associated with "old_time_stamp"
       
      !  yrnew    -  indicates the year associated with "new_time_stamp"
      !  monew    -  indicates the month associated with "new_time_stamp"
      !  dynew    -  indicates the day associated with "new_time_stamp"
      !  hrnew    -  indicates the hour associated with "new_time_stamp"
      !  minew    -  indicates the minute associated with "new_time_stamp"
      !  scnew    -  indicates the second associated with "new_time_stamp"
       
      !  mday     -  a list assigning the number of days in each month
      
      !  i        -  loop counter
      !  nday     -  the integer number of days represented by "delta_time"
      !  nhour    -  the integer number of hours in "delta_time" after taking out
      !              all the whole days
      !  nmin     -  the integer number of minutes in "delta_time" after taking out
      !              all the whole days and whole hours.
      !  nsec     -  the integer number of minutes in "delta_time" after taking out
      !              all the whole days, whole hours, and whole minutes.
       
    INTEGER :: nlen, olen
    INTEGER :: yrnew, monew, dynew, hrnew, minew, scnew
    INTEGER :: yrold, moold, dyold, hrold, miold, scold
    integer, save :: mday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    INTEGER :: nday, nhour, nmin, nsec, i
    LOGICAL :: opass
    CHARACTER (LEN=1) :: sp

    !  Break down old hdate into parts

    olen = LEN_TRIM(old_time_stamp)
    IF (olen.GE.9) THEN
       sp = old_time_stamp(9:9)
    else
       sp = '_'
    END IF

    call get_time(old_time_stamp,yrold, moold, dyold, hrold, miold, scold, opass)

    !  Set the number of days in February for that year.

    mday(2) = nfeb(yrold)

    !  Compute the number of days, hours, minutes, and seconds in delta_time

    IF (olen.eq.15) THEN  !delta_time should be in seconds
       nday   = ABS(delta_time)/86400 ! Integer number of days in delta-time
       nhour  = MOD(ABS(delta_time),86400)/3600
       nmin   = MOD(ABS(delta_time),3600)/60
       nsec   = MOD(ABS(delta_time),60)
    ELSE IF (olen.eq.13) THEN !delta_time should be in minutes
       nday   = ABS(delta_time)/1440 ! Integer number of days in delta-time
       nhour  = MOD(ABS(delta_time),1440)/60
       nmin   = MOD(ABS(delta_time),60)
       nsec   = 0
    ELSE IF (olen.eq.11) THEN !delta_time should be in hours
       nday   = ABS(delta_time)/24 ! Integer number of days in delta-time
       nhour  = MOD(ABS(delta_time),24)
       nmin   = 0
       nsec   = 0
    ELSE IF (olen.eq.8) THEN !delta_time should be in days
       nday   = ABS(delta_time)/24 ! Integer number of days in delta-time
       nhour  = 0
       nmin   = 0
       nsec   = 0
    ELSE
       write (*,*) 'get_new_timestamp: Strange length for OLD_TIME_STAMP: ',olen
       opass=.FALSE.
    END IF

    new_time_stamp = ' '
    if (.not.opass .or. min(LEN(new_time_stamp),olen) .lt. 8) then
       write (*,*) 'get_new_timestamp: returning invalid new_time_stamp'
    else
       IF (delta_time.GE.0) THEN

          scnew = scold + nsec
          IF (scnew .GE. 60) THEN
             scnew = scnew - 60
             nmin  = nmin + 1
          END IF

          minew = miold + nmin
          IF (minew .GE. 60) THEN
             minew = minew - 60
             nhour  = nhour + 1
          END IF

          hrnew = hrold + nhour
          IF (hrnew .GE. 24) THEN
             hrnew = hrnew - 24
             nday  = nday + 1
          END IF

          dynew = dyold
          monew = moold
          yrnew = yrold
          DO i = 1, nday
             dynew = dynew + 1
             IF (dynew.GT.mday(monew)) THEN
                dynew = dynew - mday(monew)
                monew = monew + 1
                IF (monew .GT. 12) THEN
                   monew = 1
                   yrnew = yrnew + 1
                   ! If the year changes, recompute the number of days in February
                   mday(2) = nfeb(yrnew)
                END IF
             END IF
          END DO

       ELSE IF (delta_time.LT.0) THEN

          scnew = scold - nsec
          IF (scnew .LT. 00) THEN
             scnew = scnew + 60
             nmin  = nmin + 1
          END IF

          minew = miold - nmin
          IF (minew .LT. 00) THEN
             minew = minew + 60
             nhour  = nhour + 1
          END IF

          hrnew = hrold - nhour
          IF (hrnew .LT. 00) THEN
             hrnew = hrnew + 24
             nday  = nday + 1
          END IF

          dynew = dyold
          monew = moold
          yrnew = yrold
          DO i = 1, nday
             dynew = dynew - 1
             IF (dynew.eq.0) THEN
                monew = monew - 1
                IF (monew.eq.0) THEN
                   monew = 12
                   yrnew = yrnew - 1
                   ! If the year changes, recompute the number of days in February
                   mday(2) = nfeb(yrnew)
                END IF
                dynew = mday(monew)
             END IF
          END DO
       END IF

       !  Now construct the new mdate

       nlen = min(LEN(new_time_stamp),olen)

       IF (nlen.eq.15) THEN
          WRITE(new_time_stamp(1:15),'(I4,I2.2,I2.2,a1,I2.2,I2.2,I2.2)') yrnew, monew, dynew, sp, hrnew, minew, scnew
       ELSE IF (nlen.eq.13) THEN
          WRITE(new_time_stamp(1:13),'(I4,I2.2,I2.2,a1,I2.2,I2.2,I2.2)') yrnew, monew, dynew, sp, hrnew, minew
       ELSE IF (nlen.eq.11) THEN
          WRITE(new_time_stamp(1:11),'(I4,I2.2,I2.2,a1,I2.2,I2.2,I2.2)') yrnew, monew, dynew, sp, hrnew
       ELSE IF (nlen.eq.8) THEN
          WRITE(new_time_stamp(1:8),'(I4,I2.2,I2.2,a1,I2.2,I2.2,I2.2)') yrnew, monew, dynew
       else
          write (*,*) 'get_new_timestamp: returning invalid new_time_stamp because of invalid nlen=',nlen
       END IF
    end if
    
  end SUBROUTINE get_new_timestamp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  FUNCTION nfeb ( year ) RESULT (num_days)
   
      ! Compute the number of days in February for the given year
   
    IMPLICIT NONE
   
    INTEGER :: year
    INTEGER :: num_days

!!$#ifdef NO_LEAP_CALENDAR
!!$      num_days = 28 ! February always has 28 days for No Leap Calendar ...
!!$#else
    num_days = 28 ! By default, February has 28 days ...
    IF (MOD(year,4).eq.0) THEN  
       num_days = 29       ! But every four years, it has 29 days ...
       IF (MOD(year,100).eq.0) THEN
          num_days = 28    ! Except every 100 years, when it has 28 days ...
          IF (MOD(year,400).eq.0) THEN
             num_days = 29 ! Except every 400 years, when it has 29 days.
          END IF
       END IF
    END IF
!!$#endif
   
  END FUNCTION nfeb

  subroutine get_time(time_stamp,yr, mo, dy, hr, mi, sc, opass)

    CHARACTER (LEN=*) , INTENT(IN)  :: time_stamp
    INTEGER, intent(out) :: yr, mo, dy, hr, mi, sc
    LOGICAL, intent(out) :: opass
      
    integer, save :: mday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
       
      !  yr    -  indicates the year associated with "time_stamp"
      !  mo    -  indicates the month associated with "time_stamp"
      !  dy    -  indicates the day associated with "time_stamp"
      !  hr    -  indicates the hour associated with "time_stamp"
      !  mi    -  indicates the minute associated with "time_stamp"
      !  sc    -  indicates the second associated with "time_stamp"
       
    INTEGER :: olen, ierr

    !  Break down hdate into parts

    hr = 0
    mi = 0
    sc = 0
    olen = LEN_TRIM(time_stamp)
    opass = .TRUE.

    !  Use internal READ statements to convert the CHARACTER string
    !  date into INTEGER components.

    READ(time_stamp(1:8),  '(I4,I2,I2)', iostat=ierr) yr, mo, dy
    IF (olen.GE.11) THEN
       if (ierr .eq. 0) READ(time_stamp(10:11),'(I2)', iostat=ierr) hr
       IF (olen.GE.13) THEN
          if (ierr .eq. 0) READ(time_stamp(12:13),'(I2)', iostat=ierr) mi
          IF (olen.GE.15) THEN
             if (ierr .eq. 0) READ(time_stamp(14:15),'(I2)', iostat=ierr) sc
          END IF
       END IF
    END IF

    !  Set the number of days in February for that year.

    mday(2) = nfeb(yr)

    !  Check that TIME_STAMP makes sense.

    !  Check that the month of TIME_STAMP makes sense.

    IF ((mo.GT.12).or.(mo.LT.1)) THEN
       WRITE(*,*) 'get_time:  Month of TIME_STAMP = ', mo
       opass = .FALSE.
    END IF

    !  Check that the day of TIME_STAMP makes sense.

    if (mo .ge. 1 .and. mo .le. 12) then
       IF ((dy.GT.mday(mo)).or.(dy.LT.1)) THEN
          WRITE(*,*) 'get_time:  Day of TIME_STAMP = ', dy
          opass = .FALSE.
       END IF
    end if
    
    !  Check that the hour of TIME_STAMP makes sense.

    IF ((hr.GT.23).or.(hr.LT.0)) THEN
       WRITE(*,*) 'get_time:  Hour of TIME_STAMP = ', hr
       opass = .FALSE.
    END IF

    !  Check that the minute of TIME_STAMP makes sense.

    IF ((mi.GT.59).or.(mi.LT.0)) THEN
       WRITE(*,*) 'get_time:  Minute of TIME_STAMP = ', mi
       opass = .FALSE.
    END IF

    !  Check that the second of TIME_STAMP makes sense.

    IF ((sc.GT.59).or.(sc.LT.0)) THEN
       WRITE(*,*) 'get_time:  Second of TIME_STAMP = ', sc
       opass = .FALSE.
    END IF

    if (.not. opass) &
         write (*,*) 'Error reading date/time from string: ',trim(time_stamp)
    
    !  Date Checks are completed.  Continue.

  end subroutine get_time

end module get_fname_module

