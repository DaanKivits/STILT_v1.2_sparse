!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:  DECODR           Decodes real input variables
!   PRGMMR:    ROLAND DRAXLER   ORG: R/ARL       DATE:02-07-12
!
! ABSTRACT:  THIS CODE WRITTEN AT THE AIR RESOURCES LABORATORY ...
!   DECODES THE VARIABLES ON AN INPUT LINE BY FIRST READING A 
!   CHARACTER STRING FOR THE ENTIRE LINE AND THEN INTERPRETING
!   THE CHARACTER STRING ONE VARIABLE AT A TIME UNTIL THE LAST 
!   INPUT VARIABLE IS PROCESSED
!
! PROGRAM HISTORY LOG:
!   LAST REVISED: 12 Jul 2002 (RRD) - initial version
!
! USAGE:  CALL DECODR(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5)
!
!   INPUT ARGUMENT LIST:     see below
!   OUTPUT ARGUMENT LIST:    see below
!   INPUT FILES:             none 
!   OUTPUT FILES:            none 
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM RS6000
!
!$$$

SUBROUTINE DECODR(IUNIT,VAR1,VAR2,VAR3,VAR4,VAR5,VAR6,VAR7,VAR8)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: IUNIT   ! unit number

  REAL,    OPTIONAL, INTENT(INOUT) :: VAR1
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR2
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR3
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR4
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR5
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR6
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR7
  REAL,    OPTIONAL, INTENT(INOUT) :: VAR8
  
  CHARACTER(130) :: LABEL         ! character string
  INTEGER       :: KLEN          ! string length
  INTEGER       :: KP            ! string position
  INTEGER       :: KV            ! variable number 
  INTEGER       :: KEND          ! relative end            
  INTEGER       :: NVAR          ! variable counter

  NVAR=0
! determine number of variables
  IF(PRESENT(VAR1))NVAR=NVAR+1
  IF(PRESENT(VAR2))NVAR=NVAR+1
  IF(PRESENT(VAR3))NVAR=NVAR+1
  IF(PRESENT(VAR4))NVAR=NVAR+1
  IF(PRESENT(VAR5))NVAR=NVAR+1
  IF(PRESENT(VAR6))NVAR=NVAR+1
  IF(PRESENT(VAR7))NVAR=NVAR+1
  IF(PRESENT(VAR8))NVAR=NVAR+1

! multiple input variables read on a single line
  READ(IUNIT,'(A)')LABEL

! check for optional label information following #
  KLEN=INDEX(LABEL,'#')-1
  IF(KLEN.LE.0)KLEN=130

! initialize string and variable pointers
  KP=1
  KV=0

! loop through entire string by variable
  DO WHILE (KP.LT.KLEN.AND.KV.LT.NVAR)

!    find next non-blank character
     DO WHILE (KP.LT.KLEN.AND.LABEL(KP:KP).EQ.' ')
        KP=KP+1
     END DO

!    find relative end of character field
     KEND=INDEX(LABEL(KP:),' ')

     IF(KEND.GT.0.AND.KP.LT.KLEN)THEN
        KV=KV+1

        IF(KV.EQ.1) READ(LABEL(KP:),*) VAR1
        IF(KV.EQ.2) READ(LABEL(KP:),*) VAR2
        IF(KV.EQ.3) READ(LABEL(KP:),*) VAR3
        IF(KV.EQ.4) READ(LABEL(KP:),*) VAR4
        IF(KV.EQ.5) READ(LABEL(KP:),*) VAR5
        IF(KV.EQ.6) READ(LABEL(KP:),*) VAR6
        IF(KV.EQ.7) READ(LABEL(KP:),*) VAR7
        IF(KV.EQ.8) READ(LABEL(KP:),*) VAR8

!       convert end back to absolute units
        KP=KP+KEND-1
     END IF

  END DO

END SUBROUTINE decodr 
