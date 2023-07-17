!# CSU IDENTIFICATION : delta_time
!#     $Id: delta_time.f,v 1.2 2005/07/05 19:39:43 trn Exp $

!# PURPOSE : Perform date/time calculations

!# CHANGE LOG : 
!#      $Log: delta_time.f,v $
!#      Revision 1.2  2005/07/05 19:39:43  trn
!#      untabify
!#
!#      Revision 1.2  2005/07/05 19:39:43  trn
!#      untabify
!#
!#
!#      Revision 1.1  2005/07/05 19:26:03  trn
!#      Initial revision of WRF GRIB to ARL converter
!#      
!#      1.1 Based on: Revision 1.6  2004/02/20 16:46:09  jmh of native/delta_time.f
!#      1.2 same as:  Revision 1.7  2005/07/05 19:35:54  trn of native/delta_time.f
!#                  

!# CSU SPECIFICATION AND CONSTRAINTS:

!# REQUIREMENTS : 
!# compute time difference between two date/times, or
!# compute new date/time from old date/time plus a positive lead time

!# CONSTRAINTS : 
!# (1) date/times are coded as 14-character strings, as follows:
!#  YYYYMMDDhhmmss, where
!#    YYYY - 4-digit year; MM - 2 digit month; DD - 2-digit day of month
!#    hh - 2-digit hour (24-hr clock); 
!#    mm - 2-digit minutes; ss - 2-digit seconds
!# (2) time differences are given in units specified by the unit_code
!#     They are truncated (not rounded) from seconds to integer values
!#     of the specified units.

!# LANGUAGE : Fortran

!# CSU DESIGN :

!# INPUT/OUTPUT INTERFACE :

      MODULE delta_time_m   !#
!# GLOBAL AND SHARED DATA : None
      INTERFACE delta_time  !only needed for overloading        !#
         MODULE PROCEDURE delta_time_diff       !#t1,t2 -> delta t
         MODULE PROCEDURE delta_time_date       !#t1, delta t -> t2
      END INTERFACE     !#
      CONTAINS        !#
! ----------------------------------------------------------------------
      SUBROUTINE delta_time_diff (                                      &
     &     T1, T2, unit_code, time_diff, ierr)  !#

!#  * Routine gets passed two time strings and returns the
!#      difference in specified units.  This value is positive if the
!#      second time is later than the first.

      implicit none

      character*(*), intent(in) :: T1, T2       !#input date/time strings
      integer, intent(in) :: unit_code  !#unit_code for output:
      integer, intent(out) :: time_diff !#time difference
      integer, intent(out) :: ierr      !#error code:
!#   0 - no error; 
!#   1 - unsupported unit_code; 
!#   12-14 - invalid input date/time string(s)

!# DATA CONVERSION : None

!# ALGORITHMS : 

!# REFERENCES : None

!# LIMITATIONS : 


!# LOCAL DATA ELEMENTS : TBD

      integer low_year  !#
      integer :: year1, month1, day1, hour1, minute1, second1   !#
      integer :: year2, month2, day2, hour2, minute2, second2   !#
      integer unit_factor, time_dif1, time_dif2 !#

!# LOCAL DATA STRUCTURES : None

!# DATA FILES : None

!# LOGIC FLOW AND DETAILED ALGORITHM: 

      ierr = 0  !#no error
      if (unit_code .eq. -1) then       !#
         unit_factor = 1        !seconds        !#
      elseif (unit_code .eq. 0) then    !#
         unit_factor = 60       !minutes        !#
      elseif (unit_code .eq. 1) then    !#
         unit_factor = 60*60    !hours  !#
      elseif (unit_code .eq. 2) then    !#
         unit_factor = 24*60*60 !days   !#
      else
         print *,'Unsupported unit code=',unit_code
         ierr=1 !#Unsupported unit code
         return
      endif
!#~  
!#~  * Routine first calculates a low year, which is the minimum
!#~      year of T1 or T2

      call decode_datetime(T1,                                          &
     &     year1, month1, day1, hour1, minute1, second1, ierr)
      if (ierr .ne. 0) return
      call decode_datetime(T2,                                          &
     &     year2, month2, day2, hour2, minute2, second2, ierr)
      if (ierr .ne. 0) return

      low_year = min( year1, year2 )    !#

!#~  
!#~  * Then the dif_from_year routine calculates the dif in seconds
!#~      between T1 and T2 from the low year, Jan 1st, 00:00:00.

      call dif_from_year( year1, month1, day1, hour1, minute1, second1, &
     &     low_year, time_dif1 ) !#

      call dif_from_year( year2, month2, day2, hour2, minute2, second2, &
     &     low_year, time_dif2 ) !#

!#~      These two difs are subtracted and scaled for the result.

      time_diff = (time_dif2 - time_dif1) / unit_factor !#
      return
      end subroutine delta_time_diff    !#
! ----------------------------------------------------------------------
      SUBROUTINE delta_time_date (                                      &
     &     T1, unit_code, in_time_diff, out_T2, ierr)   !#

      implicit none

!#  * Routine gets passed one time string and a positive lead time
!#      (difference) in specified units.  Returns the date/time string
!#      of the later time.

      character*(*), intent(in) :: T1   !#initial time date/time string
      integer, intent(in) :: unit_code, in_time_diff !#lead-time: units, value
      character*(*), intent(out) :: out_T2      !#new date/time string
      integer, intent(out) :: ierr      !#error code
!#   0 - no error; 
!#   1 - unsupported unit_code; 
!#   4 - negative in_time_diff
!#   12-14 - invalid input date/time string(s)
!#   15 - error writing output date/time string

!#~  ***** Local vars
      integer unit_factor, idts, nday, nhour, nmin, nsec        !#
      integer year, month, day, hour, minute, second    !#
      integer yrnew, monew, dynew, hrnew, minew, scnew
      integer i !#
      integer mday( 12 )        !#

!#~  ***** Set up mday array, compute unit factor

      mday= (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      if (unit_code .eq. -1) then       !#
         unit_factor = 1        !seconds        !#
      elseif (unit_code .eq. 0) then    !#
         unit_factor = 60       !minutes        !#
      elseif (unit_code .eq. 1) then    !#
         unit_factor = 60*60    !hours  !#
      elseif (unit_code .eq. 2) then    !#
         unit_factor = 24*60*60 !days   !#
      else
         print *,'delta_time: Unsupported unit code=',unit_code
         ierr=1 !#Unsupported unit code
         return
      endif

      call decode_datetime(T1,                                          &
     &     year, month, day, hour, minute, second, ierr)
      if (ierr .ne. 0) return
      if (in_time_diff .lt. 0) then
         print *,'Negative time_diff not supported in delta_time'
         ierr=4
         return
      endif

!#   Convert the time difference into days, hours, minutes, seconds:
      idts = in_time_diff * unit_factor
      nday   = idts/86400  ! Integer number of days in delta-time
      nhour   = mod(idts,86400)/3600
      nmin   = mod(idts,3600)/60
      nsec   = mod(idts,60)

!#   Compute the new time:
      scnew = second + nsec
      if (scnew .ge. 60) then
         scnew = scnew - 60
         nmin  = nmin + 1
      end if
      minew = minute + nmin
      if (minew .ge. 60) then
         minew = minew - 60
         nhour  = nhour + 1
      end if
      hrnew = hour + nhour
      if (hrnew .ge. 24) then
         hrnew = hrnew - 24
         nday  = nday + 1
      end if

!#   Compute the new date:
      dynew = day
      monew = month
      yrnew = year

      do i = 1, nday
         dynew = dynew + 1

            mday(2) = 28
            if (mod(yrnew,4).eq.0) then
               mday(2) = 29
               if (mod(yrnew,100).eq.0) then
                  mday(2) = 28
                  if (mod(yrnew,400).eq.0) then
                     mday(2) = 29
                  endif
               endif
            endif

         if (dynew.gt.mday(monew)) then
            dynew = dynew - mday(monew)
            monew = monew + 1
            if (monew .gt. 12) then
               monew = 1
               yrnew = yrnew + 1
            endif
         endif
      enddo
!#
!#  Encode the output time string
!#

      write (out_T2, '(i4.4,5i2.2)', iostat=ierr)                       &
     &     yrnew, monew, dynew, hrnew, minew, scnew
      if (ierr .ne. 0) ierr = 15

      return
      end SUBROUTINE delta_time_date
! ----------------------------------------------------------------------
      subroutine dif_from_year( year, month, day, hour, minute, second, &
     &     ref_year, difference )       !#

!#~  * Routine returns seconds difference from ref_year, Jan 1st,
!#~      00:00:00, and the passed time.  Routine assumes ref_year
!#~      is earlier than the passed in time.
!#~  
!#~  * Routine adds years, days, months, seconds, etc., to the 
!#~      ref_year until it equals the passed in time, and returns
!#~      the number of seconds added.

      implicit none

!#~  ***** Argument vars
      integer, intent(in) :: year, month, day, hour, minute, second     !#
      integer, intent(in) :: ref_year   !#
      integer, intent(out) :: difference        !#

!#~  ***** Local vars
      integer i !#
      integer sec_per_day, sec_per_hour !#
      parameter( sec_per_hour = 60 * 60,                                &
     &           sec_per_day = 24 * sec_per_hour )      !#
      integer mday( 12 )        !#

!#~  ***** Set up mday array

      mday= (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

!#~  ***** Get dif from ref_year, 00/00, 00:00:00, to time T.

      difference = 0

!#~  ***** Add the years.  Leap years are 366 days, and they happen
!#~  *****   when the year is divisible by four, with the exception
!#~  *****   that the only century years that are leap years are
!#~  *****   those divisible by 400.

      do 10, i = ref_year, year-1
         if( ( mod( I, 4 ).eq.0 )  .and.                                &
     &       ( mod( I, 100 ).ne.0  .or.  mod( I, 400 ).eq.0 ) ) then
               difference = difference + 366 * sec_per_day
         else
               difference = difference + 365 * sec_per_day
         endif
 10   continue

!#~  ***** Add the month.  If a leap year and it is February,
!#~  *****  you must add an extra day.

      do 20, i = 1, month-1
         difference = difference + sec_per_day * mday( i )
         if( ( i.eq.2 )               .and.                              &
     &       ( mod( year, 4 ).eq.0 )  .and.                              &
     &       ( mod( year, 100 ).ne.0 .or. mod( year, 400 ).eq.0 ) )      &
     &        difference = difference + sec_per_day
 20   continue

!#~  ***** Now we can just add the seconds for the day, month, and
!#~  *****   second.  There is no day 0, so use day-1.  Hours go
!#~  *****   from 0 to 23.
!#~  ***** 

      difference = difference + ( (day-1) * sec_per_day )                &
     &                              + ( (hour) * sec_per_hour )          &
     &                              + ( (minute) * 60 )                  &
     &                              + second

      return
      end subroutine dif_from_year
! ----------------------------------------------------------------------
      subroutine decode_datetime (T,                                     &
     &     year, month, day, hour, minute, second, ierr) !#

!#  Decode date/time string, perform error checking

      character *(*), intent(in) :: T
      integer, intent(out) :: year, month, day, hour, minute, second    !#
      integer, intent(out) :: ierr
      
      ierr = 0

!#~  ***** Diagnostic check on the time string:
      if( len(t) .lt. 14 .or.                                              &
     &     index(T(1:14), ':').ne.0 .or. index(T(1:14), ' ').ne.0 .or.     &
     &    index(T(1:14), ',').ne.0 ) then
         ierr=12        !#date/time has wrong length or illegal characters
         return
      endif

!#~  ***** Do internal file reads for the variables

      read( T(1:14),'(I4,5i2)', iostat=ierr ) year, month, day, hour,       &
     &     minute, second
      if (ierr .ne. 0) then
         ierr=13 !#error reading input date/time string
         return
      endif

!#~  ***** Sanity check, no day 0, month 0:

      if( month.lt.1  .or.  month.gt.12  .or.                             &
     &    day.lt.1  .or.  day.gt.31  .or.                                 &
     &    hour .gt. 23  .or.  minute.gt.59  .or.  second.gt.59 )          &
     &     ierr = 14    !#invalid month, dat, hour, minute, second values

      return
      end subroutine decode_datetime

      end module delta_time_m
