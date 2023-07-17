!===============================================================================
! analyze netCDF header for decoding information

SUBROUTINE ncscan (handle,varid,kdate,kdim,vdata,nbyte,begin,offset,scale,diag)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: handle
  CHARACTER(8), INTENT(IN)  :: varid
  INTEGER,      INTENT(OUT) :: kdate (3)
  INTEGER,      INTENT(OUT) :: kdim  (4)
  INTEGER,      INTENT(OUT) :: vdata 
  INTEGER,      INTENT(OUT) :: nbyte 
  INTEGER,      INTENT(OUT) :: begin 
  REAL,         INTENT(OUT) :: offset
  REAL,         INTENT(OUT) :: scale
  LOGICAL,      INTENT(IN)  :: diag

  LOGICAL        :: flag, vset
  REAL           :: rval
  INTEGER        :: nr, nv, na, ne, nd, nl 
  INTEGER        :: byte, nlen, klen, vdim, nval, rank, vlen
  INTEGER        :: nrec, ncor, ndim, natg, natr, nelm, ntyp, ncva, nvar, noff
  CHARACTER(1)   :: buff(50000)
  CHARACTER(512) :: label

  kdate= 0
  kdim = 1

  vdata= 0
  begin= 0
  nbyte= 0

  offset= 0.0 
  scale = 0.0

  byte=0
  flag=.false.
  vlen=index(varid,'&')-1

! file identification
  call nc2dat (handle,byte,4,buff,nval,rval)
  byte=byte+4

! number of records
  call nc2dat (handle,byte,4,buff,nrec,rval)
  if(diag) write(*,*)'Number of records: ',nrec
  byte=byte+4

  if(diag)then
     write(*,*)
     write(*,*)'-------------------------------------------> dim_array'
  end if

  call nc2dat (handle,byte,4,buff,ncor,rval)
  if(diag) write(*,*)'nc_dimension: ',ncor
  byte=byte+4

! number of dimensions
  call nc2dat (handle,byte,4,buff,ndim,rval)
  if(diag) write(*,*)'nelms: ',ndim
  byte=byte+4

  DO nd=1,ndim
!    string length
     call nc2dat (handle,byte,4,buff,nlen,rval)
     klen=4*(int((nlen-1)/4)+1)
     byte=byte+4
      
!    variable name
     call nc2dat (handle,byte,klen,buff,nval,rval)
     do nl=1,nlen 
        label(nl:nl)=buff(nl)
     end do
     if(diag) write(*,*)'  name: ',label(1:nlen) 
     byte=byte+klen
     
!    variable dimension
     call nc2dat (handle,byte,4,buff,vdim,rval)
     if(vdim.eq.0)vdim=nrec
     if(diag) write(*,*)'  dim_length: ',vdim 
     byte=byte+4

!---------------------------------------------------------
! grid dimensions
  IF(label(1:nlen).eq.'lon')  kdim(1)=vdim
  IF(label(1:nlen).eq.'lat')  kdim(2)=vdim
  IF(label(1:nlen).eq.'level')kdim(3)=vdim
  IF(label(1:nlen).eq.'time') kdim(4)=vdim
!---------------------------------------------------------
  END DO

  if(diag)then
     write(*,*)
     write(*,*)'-------------------------------------------> gatt_array'
  end if

  call nc2dat (handle,byte,4,buff,natg,rval)
  if(diag) write(*,*)'nc_attribute: ',natg
  byte=byte+4

! number of attributes
  call nc2dat (handle,byte,4,buff,natr,rval)
  if(diag) write(*,*)'nelems: ',natr
  byte=byte+4

  DO na=1,natr
!    string length
     call nc2dat (handle,byte,4,buff,nlen,rval)
     klen=4*(int((nlen-1)/4)+1)
     byte=byte+4
      
!    element name 
     call nc2dat (handle,byte,klen,buff,nval,rval)
     do nl=1,nlen 
        label(nl:nl)=buff(nl)
     end do
     if(diag) write(*,*)'  name: ',na,' - ',label(1:nlen) 
     byte=byte+klen

!---------------------------------------------------------
! starting date for data in this file
  IF(label(1:nlen).eq.'base_date')THEN
     flag=.true.
  ELSE
     flag=.false.
  END IF
!---------------------------------------------------------

!    element type      
     call nc2dat (handle,byte,4,buff,ntyp,rval)
     if(diag) write(*,*)'  nc_type: ',ntyp 
     byte=byte+4

!    number of elements in this attribute
     call nc2dat (handle,byte,4,buff,nelm,rval)
     if(diag) write(*,*)'  nelems: ',nelm
     byte=byte+4

     IF(ntyp.eq.1.or.ntyp.eq.2)THEN
!       character string  
        call nc2dat (handle,byte,nelm,buff,nval,rval)
        do ne=1,nelm 
           if(ichar(buff(ne)).eq.0)buff(ne)=' '
           if(ichar(buff(ne)).eq.10)buff(ne)=' '
           label(ne:ne)=buff(ne)
        end do
        if(diag) write(*,*)'  '//label(1:nelm)
        nlen=4*(int((nelm-1)/4)+1)
        byte=byte+nlen

     ELSEIF(ntyp.eq.3)THEN
!       short integer
        do ne=1,nelm 
           call nc2dat (handle,byte,2,buff,nval,rval)
           if(diag) write(*,*)'  Int2: ',nval
           byte=byte+2
!---------------------------------------------------------
! save all elements of the starting date (year month day)
  IF(flag) kdate(ne)=nval
!---------------------------------------------------------
        end do
        if(mod(nelm,2).ne.0)byte=byte+2

     ELSEIF(ntyp.eq.4)THEN
!       long integer
        do ne=1,nelm 
           call nc2dat (handle,byte,4,buff,nval,rval)
           if(diag) write(*,*)'  Int4: ',nval
           byte=byte+4
        end do

     ELSEIF(ntyp.eq.5)THEN
!       real number
        do ne=1,nelm 
           call nc2dat (handle,byte,4,buff,nval,rval)
           if(diag) write(*,*)'  Real: ',rval
           byte=byte+4
        end do

     ELSEIF(ntyp.eq.6)THEN
!       double precision
        do ne=1,nelm 
           call nc2dat (handle,byte,4,buff,nval,rval)
           if(diag) write(*,*)'  Dble: ',rval
           byte=byte+8
        end do
     END IF
  END DO

  if(diag)then
     write(*,*)
     write(*,*)'-------------------------------------------> var_array'
  end if

  call nc2dat (handle,byte,4,buff,ncva,rval)
  if(diag) write(*,*)'nc_variable: ',ncva
  byte=byte+4

! number of variables 
  call nc2dat (handle,byte,4,buff,nvar,rval)
  if(diag) write(*,*)'nelems: ',nvar
  byte=byte+4

  DO nv=1,nvar

     if(diag)then
        write(*,*)
        write(*,*)'-------------------------------------------> var ',nv 
     end if

!    string length
     call nc2dat (handle,byte,4,buff,nlen,rval)
     klen=4*( int( (nlen-1)/4 )+1  )
     byte=byte+4
      
!    variable name
     call nc2dat (handle,byte,klen,buff,nval,rval)
     do nl=1,nlen 
        label(nl:nl)=buff(nl)
     end do
     if(diag) write(*,*)'  name: ',nv,' - ',label(1:nlen) 
     byte=byte+klen

!---------------------------------------------------------
! requested data variable
  IF(label(1:nlen).eq.varid(1:vlen).or.varid(1:vlen).eq.'dummy')THEN
     flag=.true.
  ELSE
     flag=.false.
  END IF

! level information variable
  IF(label(1:nlen).eq.'level')THEN
     vset=.true.
  ELSE
     vset=.false.
  END IF
!---------------------------------------------------------

!    dimension rank for this element      
     call nc2dat (handle,byte,4,buff,rank,rval)
     if(diag) write(*,*)'  rank: ',rank 
     byte=byte+4

     do nr=1,rank
!       dimension id      
        call nc2dat (handle,byte,4,buff,nval,rval)
        if(diag) write(*,*)'  dimid -',nr,':',nval 
        byte=byte+4
     end do

     call nc2dat (handle,byte,4,buff,natg,rval)
     if(diag) write(*,*)'  vatt_attribute: ',natg
     byte=byte+4

!    number of attributes
     call nc2dat (handle,byte,4,buff,natr,rval)
     if(diag) write(*,*)'  nelems: ',natr
     byte=byte+4

     DO na=1,natr

        if(diag)then
           write(*,*)
           write(*,*)'-------------------------------------------> atr ',na 
        end if

!       string length
        call nc2dat (handle,byte,4,buff,nlen,rval)
        klen=4*(int((nlen-1)/4)+1)
        byte=byte+4
      
!       element name 
        call nc2dat (handle,byte,klen,buff,nval,rval)
        do nl=1,nlen 
           label(nl:nl)=buff(nl)
        end do
        if(diag) write(*,*)'    name: ',na,' - ',label(1:nlen) 
        byte=byte+klen

!       element type      
        call nc2dat (handle,byte,4,buff,ntyp,rval)
        if(diag) write(*,*)'    nc_type: ',ntyp 
        byte=byte+4

!       number of elements in this attribute
        call nc2dat (handle,byte,4,buff,nelm,rval)
        if(diag) write(*,*)'    nelems: ',nelm
        byte=byte+4

        IF(ntyp.eq.1.or.ntyp.eq.2)THEN
!          character string  
           call nc2dat (handle,byte,nelm,buff,nval,rval)
           do ne=1,nelm 
              if(ichar(buff(ne)).eq.0)buff(ne)=' '
              if(ichar(buff(ne)).eq.10)buff(ne)=' '
              label(ne:ne)=buff(ne)
           end do
           if(diag) write(*,*)'    ',label(1:nelm)
           nlen=4*(int((nelm-1)/4)+1)
           byte=byte+nlen

        ELSEIF(ntyp.eq.3)THEN
!          short integer
           do ne=1,nelm 
              call nc2dat (handle,byte,2,buff,nval,rval)
              if(diag) write(*,*)'    Int2: ',nval
              byte=byte+2
           end do
           if(mod(nelm,2).ne.0)byte=byte+2

        ELSEIF(ntyp.eq.4)THEN
!          long integer
           do ne=1,nelm 
              call nc2dat (handle,byte,4,buff,nval,rval)
              if(diag) write(*,*)'    Int4: ',nval
              byte=byte+4
           end do

        ELSEIF(ntyp.eq.5)THEN
!          real number
           do ne=1,nelm 
              call nc2dat (handle,byte,4,buff,nval,rval)
              if(diag) write(*,*)'    Real: ',rval
              byte=byte+4
           end do

        ELSEIF(ntyp.eq.6)THEN
!          double precision
           do ne=1,nelm 
              call nc2dat (handle,byte,4,buff,nval,rval)
              if(diag) write(*,*)'    Dble: ',rval
              byte=byte+8
           end do
        END IF

!---------------------------------------------------------
! packing information scale factor and offset
  IF(flag)THEN
     IF(label(1:nlen).eq.'add_offset')  offset=rval
     IF(label(1:nlen).eq.'scale_factor') scale=rval
  END IF
!---------------------------------------------------------
     END DO

!    element type      
     call nc2dat (handle,byte,4,buff,ntyp,rval)
     if(diag) write(*,*)'  nc_type: ',ntyp 
     byte=byte+4

!    number of elements in this attribute
     call nc2dat (handle,byte,4,buff,nelm,rval)
     if(diag) write(*,*)'  vsize: ',nelm
     byte=byte+4

!    offset in bytes                       
     call nc2dat (handle,byte,4,buff,noff,rval)
     if(diag) write(*,*)'  begin: ',noff

!---------------------------------------------------------
! byte position where data starts
  IF(flag)nbyte=nelm
  IF(flag)begin=noff
  IF(vset)vdata=noff
!---------------------------------------------------------
     byte=byte+4
  END DO

  if(diag)then
     write(*,*)
     write(*,*)'-------------------------------------------> data' 
     write(*,*)'At byte: ',byte
  end if

  900 CONTINUE

END SUBROUTINE ncscan
