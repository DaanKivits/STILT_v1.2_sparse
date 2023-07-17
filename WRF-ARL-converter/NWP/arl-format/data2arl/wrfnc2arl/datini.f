! Use emacs f90-mode:   -*- mode:f90  -*-
module datini_module    

! $Id: datini.f,v 1.3 2013/03/12 16:35:42 trn Exp $

contains

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
      character (len=52) :: LABEL
      integer :: k50
      character (len=256) :: fmt

      k50=50
      if (nlvl .ge. 100) k50=52

      inquire(iolength=lrectest) label(1:k50),cvar
      lrec = nxy+k50
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
      fmt = '(7I2,A4,I4,2E14.7)'
      if (k50 .ne. 50) fmt='(5i2,a1,i3,i2,A4,I4,2E14.7)'

!     index record
      if (k50 .eq. 50) then
         WRITE(LABEL(1:k50),'(7I2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,0,IG,'NULL',NEXP,PREC,VAR1
      else
         WRITE(LABEL(1:k50),'(5i2,a1,i3,i2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,'X',0,IG,'NULL',NEXP,PREC,VAR1
      end if
      
      MREC=MREC+1
      WRITE(KUNIT,REC=MREC)LABEL(1:k50),CVAR

      DO NL=1,NLVL
      DO NV=1,NVAR(NL)

         IL=NL-1
         if (k50 .eq. 50) then
            WRITE(LABEL(1:k50),'(7I2,A4,I4,2E14.7)') IY,IM,ID,IH,IC,IL,IG,'NULL',NEXP,PREC,VAR1
         else
            WRITE(LABEL(1:k50),'(5i2,a1,i3,i2,A4,I4,2E14.7)') IY,IM,ID,IH,IC,'X',IL,IG,'NULL',NEXP,PREC,VAR1
         end if
         MREC=MREC+1
         WRITE(KUNIT,REC=MREC)LABEL(1:k50),CVAR

      END DO
      END DO

      RETURN
      end SUBROUTINE DATINI
    end module datini_module
    
