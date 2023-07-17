!======================================================================
! convert ascii byte data to integer or real

SUBROUTINE nc2dat (handle,byte,klen,buff,nval,rval)

  IMPLICIT NONE
  
  REAL(4)       :: rval
  REAL(8)       :: dval
  INTEGER(2)    :: hval
  INTEGER(4)    :: k, handle, byte, klen, nval
  CHARACTER(1)  :: buff(klen)  
  CHARACTER(2)  :: half   
  CHARACTER(4)  :: full   
  CHARACTER(8)  :: double

  call fcptps(handle,byte,*900)
  call fcread(handle,buff,1,klen,*900)

  IF(klen.eq.2)THEN
     DO k=1,2
        half(k:k)=buff(k)
     END DO
     READ(half,'(a2)')hval
     nval=hval

  ELSEIF(klen.eq.4)THEN
     DO k=1,4
        full(k:k)=buff(k)
     END DO
     READ(full,'(a4)')nval
     READ(full,'(a4)')rval

  ELSEIF(klen.eq.8)THEN
     DO k=1,8
        double(k:k)=buff(k)
     END DO
     READ(double,'(a8)')dval
     rval=dval
  END IF

  RETURN

  900 nval=0
      rval=0.0

END SUBROUTINE nc2dat
