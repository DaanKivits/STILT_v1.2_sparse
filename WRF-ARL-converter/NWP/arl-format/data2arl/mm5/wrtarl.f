SUBROUTINE wrtarl (ix, jx, kx, ka, n0, n1, vchar0, vchar1, psfc,   &
           pblhgt, shflux, ust, swdown, precip, terrain, ua, va,   &
           wa, ta, pres, tke, qva, currentdate, time, gnum )

!-------------------------------------------------------------------------------
! Name:     Write ARL (wrtarl) format data
! Purpose:  Create output in ARL one-byte packed format.
! Notes:  
! Revised:  09 Aug 2000  Original version.  (RRD)
!           19 Dec 2001  Round time nearest hour (RRD)
!           09 Mar 2004  Added TKE field (RRD)
!-------------------------------------------------------------------------------

  USE file

  IMPLICIT NONE

  INTEGER,       INTENT(IN)    :: ix
  INTEGER,       INTENT(IN)    :: jx
  INTEGER,       INTENT(IN)    :: kx
  INTEGER,       INTENT(IN)    :: ka,n0,n1
  CHARACTER*4,   INTENT(IN)    :: vchar0     ( : )
  CHARACTER*4,   INTENT(IN)    :: vchar1     ( : )
  REAL,          INTENT(IN)    :: psfc       ( : , : )
  REAL,          INTENT(IN)    :: pblhgt     ( : , : )
  REAL,          INTENT(IN)    :: shflux     ( : , : )
  REAL,          INTENT(IN)    :: ust        ( : , : )
  REAL,          INTENT(IN)    :: swdown     ( : , : )
  REAL,          INTENT(IN)    :: precip     ( : , : )
  REAL,          INTENT(IN)    :: terrain    ( : , : )
  REAL,          INTENT(IN)    :: ua         ( : , : , : )
  REAL,          INTENT(IN)    :: va         ( : , : , : )
  REAL,          INTENT(IN)    :: wa         ( : , : , : )
  REAL,          INTENT(IN)    :: ta         ( : , : , : )
  REAL,          INTENT(IN)    :: pres       ( : , : , : )
  REAL,          INTENT(IN)    :: tke        ( : , : , : )
  REAL,          INTENT(IN)    :: qva        ( : , : , : )
  CHARACTER*24,  INTENT(IN)    :: currentdate
  REAL,          INTENT(IN)    :: time
  INTEGER,       INTENT(IN)    :: gnum

  INTEGER                      :: k,n
  LOGICAL                      :: ifirst     = .TRUE.
  REAL,          ALLOCATABLE   :: rvar       ( : , : )
  CHARACTER*1,   ALLOCATABLE   :: cvar       ( : )
  CHARACTER*1                  :: num
  CHARACTER*80                 :: fname
  INTEGER                      :: iy,im,id,ih,imin,isec,macc

  SAVE ifirst

!-------------------------------------------------------------------------------
! On first entry intialize packing subroutines and open output file
!-------------------------------------------------------------------------------

  IF( ifirst ) THEN
      write(num,'(I1)') (gnum/10)
      CALL PAKSET(iutarl,'CFG'//NUM//'_MM5',1,jx,ix,(kx+1))
      FNAME='DATA'//NUM//'.ARL'
      OPEN(iutarl,FILE=FNAME,RECL=(50+ix*jx),ACCESS='DIRECT',FORM='UNFORMATTED')
      ifirst = .FALSE.
  END IF

!-------------------------------------------------------------------------------
! Convert date string to current time                                
!-------------------------------------------------------------------------------

  READ(currentdate,FMT='(    i4.4)') iy
  READ(currentdate,FMT='( 5x,i2.2)') im
  READ(currentdate,FMT='( 8x,i2.2)') id  
  READ(currentdate,FMT='(11x,i2.2)') ih  
  READ(currentdate,FMT='(14x,i2.2)') imin 
  READ(currentdate,FMT='(17x,i2.2)') isec 
  iy = iy - INT(iy/100)*100

! round time because time of output not always on even hour/min   
  CALL TM2MIN(iy,im,id,ih,imin,macc)
  macc=15*NINT(macc/15.0) 
  CALL TM2DAY(macc,iy,im,id,ih,imin)

!-------------------------------------------------------------------------------
! Write out all the surface fields                                   
!-------------------------------------------------------------------------------

  ALLOCATE ( rvar(jx,ix) )
  ALLOCATE ( cvar(jx*ix) )

  DO n = 1, n0

     IF ( vchar0 (n) == 'PRSS' ) rvar = TRANSPOSE ( psfc    )
     IF ( vchar0 (n) == 'TPP1' ) rvar = TRANSPOSE ( precip  )
     IF ( vchar0 (n) == 'SHTF' ) rvar = TRANSPOSE ( shflux  )
     IF ( vchar0 (n) == 'USTR' ) rvar = TRANSPOSE ( ust     )
     IF ( vchar0 (n) == 'DSWF' ) rvar = TRANSPOSE ( swdown  )
     IF ( vchar0 (n) == 'SHGT' ) rvar = TRANSPOSE ( terrain )
     IF ( vchar0 (n) == 'MXHT' ) rvar = TRANSPOSE ( pblhgt  )

     CALL PAKREC(iutarl,rvar,cvar,jx,ix,(jx*ix),  &
                 vchar0(n),iy,im,id,ih,imin,INT(time/60.0),1,0)  

  END DO

!-------------------------------------------------------------------------------
! Write out all the upper air fields                                 
!-------------------------------------------------------------------------------

  DO k = kx, (kx-ka+1), -1
  DO n = 1, n1

     IF ( vchar1 (n) == 'UWND' ) rvar = TRANSPOSE ( ua  (:,:,k) )
     IF ( vchar1 (n) == 'VWND' ) rvar = TRANSPOSE ( va  (:,:,k) )
     IF ( vchar1 (n) == 'WWND' ) rvar = TRANSPOSE ( wa  (:,:,k) )
     IF ( vchar1 (n) == 'TEMP' ) rvar = TRANSPOSE ( ta  (:,:,k) )
     IF ( vchar1 (n) == 'TKEN' ) rvar = TRANSPOSE ( tke (:,:,k) )
     IF ( vchar1 (n) == 'SPHU' ) rvar = TRANSPOSE ( qva (:,:,k) )
     IF ( vchar1 (n) == 'PRES' ) rvar = TRANSPOSE ( pres(:,:,k) )

     CALL PAKREC(iutarl,rvar,cvar,jx,ix,(jx*ix),  &
                 vchar1(n),iy,im,id,ih,imin,INT(time/60.0),2+kx-k,0)  

  END DO
  END DO

  DEALLOCATE ( rvar )
  DEALLOCATE ( cvar )

!-------------------------------------------------------------------------------
! Close out time period by writing index record                      
!-------------------------------------------------------------------------------

  CALL PAKNDX(iutarl)


  RETURN

END SUBROUTINE wrtarl
