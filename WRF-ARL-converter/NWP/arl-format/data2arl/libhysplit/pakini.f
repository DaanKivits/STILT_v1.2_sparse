!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  PAKINI           PAcKed data INItialization
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:96-06-01
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   PAKED DATA INITIALIZATION - WRITES NULL VARIABLES FIELDS TO ALL
!   THE RECORDS OF AN OUTPUT FILE FOR A PARTICULAR TIME PERIOD
!
! PROGRAM HISTORY LOG:
!   Last Revised: 14 Feb 1997 (RRD)
!                 02 Feb 2001 (RRD) - fortran90 upgrade
!                 18 Oct 2001 (RRD) - extended grid domains
!                 15 Nov 2001 (RRD) - replaced NULL with variable ID
!
! USAGE:  CALL PAKINI(KG,CVAR,NXY)
!
!   INPUT ARGUMENT LIST:   see below
!   OUTPUT ARGUMENT LIST:  see below
!   INPUT FILES:           none
!   OUTPUT FILES:          unit defined in common block
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE PAKINI(KG,CVAR,NXY)

  IMPLICIT NONE

!-------------------------------------------------------------------------------

  INTEGER,      INTENT(IN)  :: kg         ! current active grid number
  INTEGER,      INTENT(IN)  :: nxy        ! dimension of packed data
  CHARACTER(1), INTENT(OUT) :: cvar(nxy)  ! packed data array

  CHARACTER(52) :: label   
  INTEGER       :: k,kk,iy,im,ih
  INTEGER       :: mrec,nhrec,nexp,nlvl,nvar
  INTEGER       :: ng,ic,nv,nl,il,id
  REAL          :: prec,var1 
  integer       :: k50

!-------------------------------------------------------------------------------

  INCLUDE 'DEFPACK.INC'
  COMMON / PAKCOM / GV, NG

!-------------------------------------------------------------------------------

  IC=-1  ! forecast field defaults to missing (-1)

! packing values all default to zero
  NEXP=0
  PREC=0.0
  VAR1=0.0

  k50=50
  NLVL=GV(KG)%NLVL
  if (nlvl .ge. 100) k50=52

! record number of index record
  MREC=GV(KG)%MREC
  IF(MREC.LT.1)THEN
     WRITE(*,*)'*ERROR* pakini: index record <1 = ',MREC
     STOP
  ELSE
     WRITE(*,*)' NOTICE pakini: start initialization rec = ',MREC
  END IF

! initialize packed data array
  DO K=1,NXY
     CVAR(K)=' '
  END DO
  CVAR(NXY)=CHAR(13)

! set dates
  IY=GV(KG)%IY0
  IM=GV(KG)%IM0
  ID=GV(KG)%ID0
  IH=GV(KG)%IH0

! index record
  IF(GV(NG)%XGPT)THEN     
     if (k50 .eq. 50) then
        WRITE(LABEL(1:K50),'(6I2,A2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,0,GV(KG)%IGC,'INDX',NEXP,PREC,VAR1
     else
        WRITE(LABEL(1:K50),'(5I2,a1,i3,A2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,'X',0,GV(KG)%IGC,'INDX',NEXP,PREC,VAR1
     end if
  ELSE
     if (k50 .eq. 50) then
        WRITE(LABEL(1:K50),'(7I2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,0,GV(KG)%IG, 'INDX',NEXP,PREC,VAR1
     else
        WRITE(LABEL(1:K50),'(5I2,a1,i3,i2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,'X',0,GV(KG)%IG, 'INDX',NEXP,PREC,VAR1
     end if
  END IF
  MREC=MREC-1
  NHREC=GV(KG)%NHREC
  DO KK=1,NHREC
     MREC=MREC+1
     WRITE(GV(KG)%KUNIT,REC=MREC)LABEL(1:k50),CVAR
  END DO

! number of level loop
  DO NL=1,NLVL
!    level indicator (0=surface)
     IL=NL-1

!    number of variables per level
     NVAR=GV(KG)%NVAR(NL)
     DO NV=1,NVAR

        IF(GV(NG)%XGPT)THEN     
           if (k50 .eq. 50) then
              WRITE(LABEL(1:K50),'(6I2,A2,A4,I4,2E14.7)')IY,IM,ID,IH,IC,IL,GV(KG)%IGC,GV(KG)%VARB(NV,NL), &
                           NEXP,PREC,VAR1
           else
              WRITE(LABEL(1:K50),'(5I2,a1,i3,A2,A4,I4,2E14.7)') IY,IM,ID,IH,IC,'X',IL,GV(KG)%IGC,GV(KG)%VARB(NV,NL), &
                           NEXP,PREC,VAR1
           end if
        ELSE
           if (k50 .eq. 50) then
              WRITE(LABEL(1:K50),'(7I2,A4,I4,2E14.7)') IY,IM,ID,IH,IC,IL,GV(KG)%IG,GV(KG)%VARB(NV,NL), &
                           NEXP,PREC,VAR1
           else
              WRITE(LABEL(1:K50),'(5I2,a1,i3,i2,A4,I4,2E14.7)') IY,IM,ID,IH,IC,'X',IL,GV(KG)%IG,GV(KG)%VARB(NV,NL), &
                           NEXP,PREC,VAR1
           end if
        END IF

        MREC=MREC+1
        WRITE(GV(KG)%KUNIT,REC=MREC)LABEL(1:K50),CVAR

     END DO
  END DO
  WRITE(*,*)' NOTICE pakini:   end initialization rec = ',MREC

END SUBROUTINE pakini
