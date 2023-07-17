      PROGRAM FINDGRIB

!-----------------------------------------------------------------------
! finds the starting GRIB string at the beginning of each grib record
! and determines the record length between GRIB records from the byte
! count and the actual record length encoded in the binary data
!-----------------------------------------------------------------------

      CHARACTER fname*40, buff(256,2)*1

!-------------------------------------------------------------------------------
! When dealing with some F90 compilers, replace ICHAR below with JCHAR function
  CHARACTER(1)        :: mychr
  INTEGER             :: jchar
  JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)
!-------------------------------------------------------------------------------

      WRITE(*,*)'Enter file name: '
      READ(*,'(a)')fname
      OPEN(10,file=fname,form='unformatted',access='direct',recl=1)

      nrec=0

!     fill first half of data buffer
      krec=1
      call fcread(10,buff(1,1),krec,256,*900)

      kbyte0=1
  100 CONTINUE
      krec=krec+256
      call fcread(10,buff(1,2),krec,256,*900)

      DO K=1,256
         kbyte=krec-257+k

         IF(buff(k  ,1).EQ.'G'.AND.                                            &
            buff(k+1,1).EQ.'R'.AND.                                            &
            buff(k+2,1).EQ.'I'.AND.                                            &
            buff(k+3,1).EQ.'B')THEN

!           increment grib record counter
            nrec=nrec+1

            IF(nrec.gt.1)THEN
!              grib record length as difference between GRIB field
               klen1=kbyte-kbyte0
               WRITE(*,*)   'Starts:',kbyte0,                                  &
                         '   Length:',klen1,'   Coded length:',klen2
!              READ(*,*)
            END IF

!           next grib record length coded in binary data
            klen2= jchar(buff(k+6,1))+                                         &
                   jchar(buff(k+5,1))*256+                                     &
                   jchar(buff(k+4,1))*65536

!           diagnostic dump of characters
!           DO kk=0,7
!              WRITE(*,*)kk, buff(k+kk,1), jchar(buff(k+kk,1))
!           END DO
!           WRITE(*,*)

            kbyte0=kbyte
         END IF
      END DO

!     copy 2nd buffer to 1st buffer
      DO K=1,256
         buff(k,1)=buff(k,2)
      END DO
      GO TO 100

  900 CLOSE(10)
      WRITE(*,*)'Number of grib records in file: ',nrec
      END

!-----------------------------------------------------------------

      SUBROUTINE fcread(kunit,buff,krec,nbyte,*)
      CHARACTER buff(nbyte)*1
      DO k=1,nbyte
         READ(kunit,rec=krec+k-1,err=900)buff(k)
      END DO
      RETURN
  900 RETURN 1
      END
