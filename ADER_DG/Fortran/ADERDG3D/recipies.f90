
MODULE recipies_mod

   IMPLICIT NONE
   
   INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
   INTEGER, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
   INTEGER, PARAMETER :: NPAR_CUMSUM=16
   INTEGER, PARAMETER :: NPAR_CUMPROD=8
   INTEGER, PARAMETER :: NPAR_POLY=8
   INTEGER, PARAMETER :: NPAR_POLYTERM=8
   INTEGER, PARAMETER :: LGT = KIND(.true.)

   REAL   , PARAMETER :: PI=3.141592653589793238462643383279502884197
   REAL   , PARAMETER :: PIO2=1.57079632679489661923132169163975144209858
   REAL   , PARAMETER :: TWOPI=6.283185307179586476925286766559005768394
   REAL   , PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967
   REAL   , PARAMETER :: NEULER=0.5772156649015328606065120900824024310422
   TYPE sprs2_sp
      INTEGER                        :: n,len
      REAL   , DIMENSION(:), POINTER :: val
      INTEGER, DIMENSION(:), POINTER :: irow
      INTEGER, DIMENSION(:), POINTER :: jcol
   END TYPE sprs2_sp


   INTERFACE ei 
      MODULE PROCEDURE ei 
   END INTERFACE

   INTERFACE expint
      MODULE PROCEDURE expint 
   END INTERFACE

   INTERFACE indexx
      MODULE PROCEDURE indexx 
   END INTERFACE

   INTERFACE factln 
      MODULE PROCEDURE factln
   END INTERFACE

   INTERFACE golden
      MODULE PROCEDURE golden
   END INTERFACE

   INTERFACE bico 
      PURE ELEMENTAL REAL FUNCTION bico(n,k)
        INTEGER , INTENT(IN) :: k,n
      END FUNCTION 
   END INTERFACE
      
   INTERFACE pythag
       MODULE PROCEDURE pythag
   END INTERFACE

      INTERFACE array_copy
      MODULE PROCEDURE array_copy_r, array_copy_i
   END INTERFACE

   INTERFACE swap
      MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
         swap_cv,swap_cm, &
         masked_swap_rs,masked_swap_rv,masked_swap_rm
   END INTERFACE

   INTERFACE reallocate
      MODULE PROCEDURE reallocate_rv,reallocate_rm,&
         reallocate_iv,reallocate_im,reallocate_hv
   END INTERFACE

   INTERFACE imaxloc
      MODULE PROCEDURE imaxloc_r,imaxloc_i
   END INTERFACE
   INTERFACE assert
      MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
   END INTERFACE
   INTERFACE assert_eq
      MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
   END INTERFACE
   INTERFACE arth
      MODULE PROCEDURE arth_r, arth_i
   END INTERFACE
   INTERFACE geop
      MODULE PROCEDURE geop_r, geop_i, geop_c
   END INTERFACE
   INTERFACE cumsum
      MODULE PROCEDURE cumsum_r,cumsum_i
   END INTERFACE
   INTERFACE poly
      MODULE PROCEDURE poly_rr,poly_rrv,&
         poly_rc,poly_cc,poly_msk_rrv
   END INTERFACE
   INTERFACE poly_term
      MODULE PROCEDURE poly_term_rr,poly_term_cc
   END INTERFACE
   INTERFACE outerprod
      MODULE PROCEDURE outerprod_r
   END INTERFACE
   INTERFACE outerdiff
      MODULE PROCEDURE outerdiff_r,outerdiff_i
   END INTERFACE
   INTERFACE scatter_add
      MODULE PROCEDURE scatter_add_r
   END INTERFACE
   INTERFACE scatter_max
      MODULE PROCEDURE scatter_max_r
   END INTERFACE
   INTERFACE diagadd
      MODULE PROCEDURE diagadd_rv,diagadd_r
   END INTERFACE
   INTERFACE diagmult
      MODULE PROCEDURE diagmult_rv,diagmult_r
   END INTERFACE
   INTERFACE get_diag
      MODULE PROCEDURE get_diag_rv
   END INTERFACE
   INTERFACE put_diag
      MODULE PROCEDURE put_diag_rv, put_diag_r
   END INTERFACE

   INTERFACE gammln
      MODULE PROCEDURE gammln
   END INTERFACE

   INTERFACE gaulag
      MODULE PROCEDURE gaulag
   END INTERFACE

   INTERFACE gauleg
      MODULE PROCEDURE gauleg
   END INTERFACE

   INTERFACE gaujac
      MODULE PROCEDURE gaujac
   END INTERFACE

   INTERFACE polint
      MODULE PROCEDURE polint
   END INTERFACE

   INTERFACE ratint
      MODULE PROCEDURE ratint
   END INTERFACE
   
   INTERFACE svbksb
      MODULE PROCEDURE svbksb
   END INTERFACE

   INTERFACE svdcmp
      MODULE PROCEDURE svdcmp
   END INTERFACE

   INTERFACE splie2
      MODULE PROCEDURE splie2
   END INTERFACE

   INTERFACE splin2
      MODULE PROCEDURE splin2
   END INTERFACE

   INTERFACE spline
      MODULE PROCEDURE spline
   END INTERFACE

   INTERFACE splint
      MODULE PROCEDURE splint
   END INTERFACE
   
   INTERFACE hpsort
      MODULE PROCEDURE hpsorti,hpsortr,hpsortia 
   END INTERFACE 
   
   INTERFACE gaussj
      MODULE PROCEDURE gaussj
   END INTERFACE 
   
   INTERFACE rg
      MODULE PROCEDURE rg
   END INTERFACE 


   !---------------------------------------------------------------------------!
   PUBLIC  :: pythag
   PUBLIC  :: svbksb
   PUBLIC  :: svdcmp
   PUBLIC  :: splie2
   PUBLIC  :: splin2
   PUBLIC  :: spline
   PUBLIC  :: splint
   PUBLIC  :: polint
   PUBLIC  :: ratint
   PUBLIC  :: gammln
   PUBLIC  :: gauleg
   PUBLIC  :: gaulag 
   PUBLIC  :: gaujac
   PUBLIC  :: hpsort
   PUBLIC  :: gaussj  
   PUBLIC  :: rg
   !---------------------------------------------------------------------------!

CONTAINS
 
  !
  FUNCTION pythag(a,b)
  !--------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  REAL :: a,b
  REAL :: pythag
  REAL:: absa,absb
  !--------------------------------------------------------------------------
  INTENT(IN) :: a,b
  !--------------------------------------------------------------------------

   absa=abs(a)
   absb=abs(b)
   if (absa > absb) then
      pythag=absa*sqrt(1.0+(absb/absa)**2)
   else
      if (absb == 0.0) then
         pythag=0.0
      else
         pythag=absb*sqrt(1.0+(absa/absb)**2)
      end if
   end if
 END FUNCTION pythag

 SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
   REAL, DIMENSION(:), INTENT(IN) :: src
   REAL, DIMENSION(:), INTENT(OUT) :: dest
   INTEGER, INTENT(OUT) :: n_copied, n_not_copied
   n_copied=min(size(src),size(dest))
   n_not_copied=size(src)-n_copied
   dest(1:n_copied)=src(1:n_copied)
   END SUBROUTINE array_copy_r
!BL
   SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
   INTEGER, DIMENSION(:), INTENT(IN) :: src
   INTEGER, DIMENSION(:), INTENT(OUT) :: dest
   INTEGER, INTENT(OUT) :: n_copied, n_not_copied
   n_copied=min(size(src),size(dest))
   n_not_copied=size(src)-n_copied
   dest(1:n_copied)=src(1:n_copied)
   END SUBROUTINE array_copy_i
!BL
!BL
   SUBROUTINE swap_i(a,b)
   INTEGER, INTENT(INOUT) :: a,b
   INTEGER :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_i
!BL
   SUBROUTINE swap_r(a,b)
   REAL, INTENT(INOUT) :: a,b
   REAL :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_r
!BL
   SUBROUTINE swap_rv(a,b)
   REAL, DIMENSION(:), INTENT(INOUT) :: a,b
   REAL, DIMENSION(SIZE(a)) :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_rv
!BL
   SUBROUTINE swap_c(a,b)
   COMPLEX, INTENT(INOUT) :: a,b
   COMPLEX :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_c
!BL
   SUBROUTINE swap_cv(a,b)
   COMPLEX, DIMENSION(:), INTENT(INOUT) :: a,b
   COMPLEX, DIMENSION(SIZE(a)) :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_cv
!BL
   SUBROUTINE swap_cm(a,b)
   COMPLEX, DIMENSION(:,:), INTENT(INOUT) :: a,b
   COMPLEX, DIMENSION(size(a,1),size(a,2)) :: dum
   dum=a
   a=b
   b=dum
   END SUBROUTINE swap_cm
!BL
   SUBROUTINE masked_swap_rs(a,b,mask)
   REAL, INTENT(INOUT) :: a,b
   LOGICAL(LGT), INTENT(IN) :: mask
   REAL :: swp
   if (mask) then
      swp=a
      a=b
      b=swp
   end if
   END SUBROUTINE masked_swap_rs
!BL
   SUBROUTINE masked_swap_rv(a,b,mask)
   REAL, DIMENSION(:), INTENT(INOUT) :: a,b
   LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
   REAL, DIMENSION(size(a)) :: swp
   where (mask)
      swp=a
      a=b
      b=swp
   end where
   END SUBROUTINE masked_swap_rv
!BL
   SUBROUTINE masked_swap_rm(a,b,mask)
   REAL, DIMENSION(:,:), INTENT(INOUT) :: a,b
   LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
   REAL, DIMENSION(size(a,1),size(a,2)) :: swp
   where (mask)
      swp=a
      a=b
      b=swp
   end where
   END SUBROUTINE masked_swap_rm
!BL
!BL
   FUNCTION reallocate_rv(p,n)
   REAL, DIMENSION(:), POINTER :: p, reallocate_rv
   INTEGER, INTENT(IN) :: n
   INTEGER :: nold,ierr
   allocate(reallocate_rv(n),stat=ierr)
   if (ierr /= 0) call &
      nrerror('reallocate_rv: problem in attempt to allocate memory')
   if (.not. associated(p)) RETURN
   nold=size(p)
   reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
   deallocate(p)
   END FUNCTION reallocate_rv
!BL
   FUNCTION reallocate_iv(p,n)
   INTEGER, DIMENSION(:), POINTER :: p, reallocate_iv
   INTEGER, INTENT(IN) :: n
   INTEGER :: nold,ierr
   allocate(reallocate_iv(n),stat=ierr)
   if (ierr /= 0) call &
      nrerror('reallocate_iv: problem in attempt to allocate memory')
   if (.not. associated(p)) RETURN
   nold=size(p)
   reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
   deallocate(p)
   END FUNCTION reallocate_iv
!BL
   FUNCTION reallocate_hv(p,n)
   CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
   INTEGER, INTENT(IN) :: n
   INTEGER :: nold,ierr
   allocate(reallocate_hv(n),stat=ierr)
   if (ierr /= 0) call &
      nrerror('reallocate_hv: problem in attempt to allocate memory')
   if (.not. associated(p)) RETURN
   nold=size(p)
   reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
   deallocate(p)
   END FUNCTION reallocate_hv
!BL
   FUNCTION reallocate_rm(p,n,m)
   REAL, DIMENSION(:,:), POINTER :: p, reallocate_rm
   INTEGER, INTENT(IN) :: n,m
   INTEGER :: nold,mold,ierr
   allocate(reallocate_rm(n,m),stat=ierr)
   if (ierr /= 0) call &
      nrerror('reallocate_rm: problem in attempt to allocate memory')
   if (.not. associated(p)) RETURN
   nold=size(p,1)
   mold=size(p,2)
   reallocate_rm(1:min(nold,n),1:min(mold,m))=&
      p(1:min(nold,n),1:min(mold,m))
   deallocate(p)
   END FUNCTION reallocate_rm
!BL
   FUNCTION reallocate_im(p,n,m)
   INTEGER, DIMENSION(:,:), POINTER :: p, reallocate_im
   INTEGER, INTENT(IN) :: n,m
   INTEGER :: nold,mold,ierr
   allocate(reallocate_im(n,m),stat=ierr)
   if (ierr /= 0) call &
      nrerror('reallocate_im: problem in attempt to allocate memory')
   if (.not. associated(p)) RETURN
   nold=size(p,1)
   mold=size(p,2)
   reallocate_im(1:min(nold,n),1:min(mold,m))=&
      p(1:min(nold,n),1:min(mold,m))
   deallocate(p)
   END FUNCTION reallocate_im
!BL
   FUNCTION ifirstloc(mask)
   LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
   INTEGER :: ifirstloc
   INTEGER, DIMENSION(1) :: loc
   loc=maxloc(merge(1,0,mask))
   ifirstloc=loc(1)
   if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
   END FUNCTION ifirstloc
!BL
   FUNCTION imaxloc_r(arr)
   REAL, DIMENSION(:), INTENT(IN) :: arr
   INTEGER :: imaxloc_r
   INTEGER, DIMENSION(1) :: imax
   imax=maxloc(arr(:))
   imaxloc_r=imax(1)
   END FUNCTION imaxloc_r
!BL
   FUNCTION imaxloc_i(iarr)
   INTEGER, DIMENSION(:), INTENT(IN) :: iarr
   INTEGER, DIMENSION(1) :: imax
   INTEGER :: imaxloc_i
   imax=maxloc(iarr(:))
   imaxloc_i=imax(1)
   END FUNCTION imaxloc_i
!BL
   FUNCTION iminloc(arr)
   REAL, DIMENSION(:), INTENT(IN) :: arr
   INTEGER, DIMENSION(1) :: imin
   INTEGER :: iminloc
   imin=minloc(arr(:))
   iminloc=imin(1)
   END FUNCTION iminloc
!BL
   SUBROUTINE assert1(n1,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   LOGICAL, INTENT(IN) :: n1
   if (.not. n1) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
         string
      STOP 'program terminated by assert1'
   end if
   END SUBROUTINE assert1
!BL
   SUBROUTINE assert2(n1,n2,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   LOGICAL, INTENT(IN) :: n1,n2
   if (.not. (n1 .and. n2)) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
         string
      STOP 'program terminated by assert2'
   end if
   END SUBROUTINE assert2
!BL
   SUBROUTINE assert3(n1,n2,n3,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   LOGICAL, INTENT(IN) :: n1,n2,n3
   if (.not. (n1 .and. n2 .and. n3)) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
         string
      STOP 'program terminated by assert3'
   end if
   END SUBROUTINE assert3
!BL
   SUBROUTINE assert4(n1,n2,n3,n4,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   LOGICAL, INTENT(IN) :: n1,n2,n3,n4
   if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
         string
      STOP 'program terminated by assert4'
   end if
   END SUBROUTINE assert4
!BL
   SUBROUTINE assert_v(n,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   LOGICAL, DIMENSION(:), INTENT(IN) :: n
   if (.not. all(n)) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
         string
      STOP 'program terminated by assert_v'
   end if
   END SUBROUTINE assert_v
!BL
   FUNCTION assert_eq2(n1,n2,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, INTENT(IN) :: n1,n2
   INTEGER :: assert_eq2
   if (n1 == n2) then
      assert_eq2=n1
   else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
         string
      STOP 'program terminated by assert_eq2'
   end if
   END FUNCTION assert_eq2
!BL
   FUNCTION assert_eq3(n1,n2,n3,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, INTENT(IN) :: n1,n2,n3
   INTEGER :: assert_eq3
   if (n1 == n2 .and. n2 == n3) then
      assert_eq3=n1
   else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
         string
      STOP 'program terminated by assert_eq3'
   end if
   END FUNCTION assert_eq3
!BL
   FUNCTION assert_eq4(n1,n2,n3,n4,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, INTENT(IN) :: n1,n2,n3,n4
   INTEGER :: assert_eq4
   if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
      assert_eq4=n1
   else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
         string
      STOP 'program terminated by assert_eq4'
   end if
   END FUNCTION assert_eq4
!BL
   FUNCTION assert_eqn(nn,string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, DIMENSION(:), INTENT(IN) :: nn
   INTEGER :: assert_eqn
   if (all(nn(2:) == nn(1))) then
      assert_eqn=nn(1)
   else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
         string
      STOP 'program terminated by assert_eqn'
   end if
   END FUNCTION assert_eqn
!BL
   SUBROUTINE nrerror(string)
   CHARACTER(LEN=*), INTENT(IN) :: string
   write (*,*) 'nrerror: ',string
   STOP 'program terminated by nrerror'
   END SUBROUTINE nrerror
!BL
   FUNCTION arth_r(first,increment,n)
   REAL, INTENT(IN) :: first,increment
   INTEGER, INTENT(IN) :: n
   REAL, DIMENSION(n) :: arth_r
   INTEGER :: k,k2
   REAL :: temp
   if (n > 0) arth_r(1)=first
   if (n <= NPAR_ARTH) then
      do k=2,n
         arth_r(k)=arth_r(k-1)+increment
      end do
   else
      do k=2,NPAR2_ARTH
         arth_r(k)=arth_r(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
         if (k >= n) exit
         k2=k+k
         arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
         temp=temp+temp
         k=k2
      end do
   end if
   END FUNCTION arth_r
!BL
   FUNCTION arth_i(first,increment,n)
   INTEGER, INTENT(IN) :: first,increment,n
   INTEGER, DIMENSION(n) :: arth_i
   INTEGER :: k,k2,temp
   if (n > 0) arth_i(1)=first
   if (n <= NPAR_ARTH) then
      do k=2,n
         arth_i(k)=arth_i(k-1)+increment
      end do
   else
      do k=2,NPAR2_ARTH
         arth_i(k)=arth_i(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
         if (k >= n) exit
         k2=k+k
         arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
         temp=temp+temp
         k=k2
      end do
   end if
   END FUNCTION arth_i
!BL
!BL
   FUNCTION geop_r(first,factor,n)
   REAL, INTENT(IN) :: first,factor
   INTEGER, INTENT(IN) :: n
   REAL, DIMENSION(n) :: geop_r
   INTEGER :: k,k2
   REAL :: temp
   if (n > 0) geop_r(1)=first
   if (n <= NPAR_GEOP) then
      do k=2,n
         geop_r(k)=geop_r(k-1)*factor
      end do
   else
      do k=2,NPAR2_GEOP
         geop_r(k)=geop_r(k-1)*factor
      end do
      temp=factor**NPAR2_GEOP
      k=NPAR2_GEOP
      do
         if (k >= n) exit
         k2=k+k
         geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
         temp=temp*temp
         k=k2
      end do
   end if
   END FUNCTION geop_r
!BL
   FUNCTION geop_i(first,factor,n)
   INTEGER, INTENT(IN) :: first,factor,n
   INTEGER, DIMENSION(n) :: geop_i
   INTEGER :: k,k2,temp
   if (n > 0) geop_i(1)=first
   if (n <= NPAR_GEOP) then
      do k=2,n
         geop_i(k)=geop_i(k-1)*factor
      end do
   else
      do k=2,NPAR2_GEOP
         geop_i(k)=geop_i(k-1)*factor
      end do
      temp=factor**NPAR2_GEOP
      k=NPAR2_GEOP
      do
         if (k >= n) exit
         k2=k+k
         geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
         temp=temp*temp
         k=k2
      end do
   end if
   END FUNCTION geop_i
!BL
   FUNCTION geop_c(first,factor,n)
   COMPLEX, INTENT(IN) :: first,factor
   INTEGER, INTENT(IN) :: n
   COMPLEX, DIMENSION(n) :: geop_c
   INTEGER :: k,k2
   COMPLEX :: temp
   if (n > 0) geop_c(1)=first
   if (n <= NPAR_GEOP) then
      do k=2,n
         geop_c(k)=geop_c(k-1)*factor
      end do
   else
      do k=2,NPAR2_GEOP
         geop_c(k)=geop_c(k-1)*factor
      end do
      temp=factor**NPAR2_GEOP
      k=NPAR2_GEOP
      do
         if (k >= n) exit
         k2=k+k
         geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
         temp=temp*temp
         k=k2
      end do
   end if
   END FUNCTION geop_c
!BL
!BL
   RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
   REAL, DIMENSION(:), INTENT(IN) :: arr
   REAL, OPTIONAL, INTENT(IN) :: seed
   REAL, DIMENSION(size(arr)) :: ans
   INTEGER :: n,j
   REAL :: sd
   n=size(arr)
   if (n == 0) RETURN
   sd=0.0
   if (present(seed)) sd=seed
   ans(1)=arr(1)+sd
   if (n < NPAR_CUMSUM) then
      do j=2,n
         ans(j)=ans(j-1)+arr(j)
      end do
   else
      ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
      ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
   end if
   END FUNCTION cumsum_r
!BL
   RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
   INTEGER, DIMENSION(:), INTENT(IN) :: arr
   INTEGER, OPTIONAL, INTENT(IN) :: seed
   INTEGER, DIMENSION(size(arr)) :: ans
   INTEGER :: n,j,sd
   n=size(arr)
   if (n == 0) RETURN
   sd=0
   if (present(seed)) sd=seed
   ans(1)=arr(1)+sd
   if (n < NPAR_CUMSUM) then
      do j=2,n
         ans(j)=ans(j-1)+arr(j)
      end do
   else
      ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
      ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
   end if
   END FUNCTION cumsum_i
!BL
!BL
   RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
   REAL, DIMENSION(:), INTENT(IN) :: arr
   REAL, OPTIONAL, INTENT(IN) :: seed
   REAL, DIMENSION(size(arr)) :: ans
   INTEGER :: n,j
   REAL :: sd
   n=size(arr)
   if (n == 0) RETURN
   sd=1.0
   if (present(seed)) sd=seed
   ans(1)=arr(1)*sd
   if (n < NPAR_CUMPROD) then
      do j=2,n
         ans(j)=ans(j-1)*arr(j)
      end do
   else
      ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
      ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
   end if
   END FUNCTION cumprod
!BL
!BL
   FUNCTION poly_rr(x,coeffs)
   REAL, INTENT(IN) :: x
   REAL, DIMENSION(:), INTENT(IN) :: coeffs
   REAL :: poly_rr
   REAL :: pow
   REAL, DIMENSION(:), ALLOCATABLE :: vec
   INTEGER :: i,n,nn
   n=size(coeffs)
   if (n <= 0) then
      poly_rr=0.0
   else if (n < NPAR_POLY) then
      poly_rr=coeffs(n)
      do i=n-1,1,-1
         poly_rr=x*poly_rr+coeffs(i)
      end do
   else
      allocate(vec(n+1))
      pow=x
      vec(1:n)=coeffs
      do
         vec(n+1)=0.0
         nn=ishft(n+1,-1)
         vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
         if (nn == 1) exit
         pow=pow*pow
         n=nn
      end do
      poly_rr=vec(1)
      deallocate(vec)
   end if
   END FUNCTION poly_rr
!BL
   FUNCTION poly_rc(x,coeffs)
   COMPLEX, INTENT(IN) :: x
   REAL, DIMENSION(:), INTENT(IN) :: coeffs
   COMPLEX :: poly_rc
   COMPLEX :: pow
   COMPLEX, DIMENSION(:), ALLOCATABLE :: vec
   INTEGER :: i,n,nn
   n=size(coeffs)
   if (n <= 0) then
      poly_rc=0.0
   else if (n < NPAR_POLY) then
      poly_rc=coeffs(n)
      do i=n-1,1,-1
         poly_rc=x*poly_rc+coeffs(i)
      end do
   else
      allocate(vec(n+1))
      pow=x
      vec(1:n)=coeffs
      do
         vec(n+1)=0.0
         nn=ishft(n+1,-1)
         vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
         if (nn == 1) exit
         pow=pow*pow
         n=nn
      end do
      poly_rc=vec(1)
      deallocate(vec)
   end if
   END FUNCTION poly_rc
!BL
   FUNCTION poly_cc(x,coeffs)
   COMPLEX, INTENT(IN) :: x
   COMPLEX, DIMENSION(:), INTENT(IN) :: coeffs
   COMPLEX :: poly_cc
   COMPLEX :: pow
   COMPLEX, DIMENSION(:), ALLOCATABLE :: vec
   INTEGER :: i,n,nn
   n=size(coeffs)
   if (n <= 0) then
      poly_cc=0.0
   else if (n < NPAR_POLY) then
      poly_cc=coeffs(n)
      do i=n-1,1,-1
         poly_cc=x*poly_cc+coeffs(i)
      end do
   else
      allocate(vec(n+1))
      pow=x
      vec(1:n)=coeffs
      do
         vec(n+1)=0.0
         nn=ishft(n+1,-1)
         vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
         if (nn == 1) exit
         pow=pow*pow
         n=nn
      end do
      poly_cc=vec(1)
      deallocate(vec)
   end if
   END FUNCTION poly_cc
!BL
   FUNCTION poly_rrv(x,coeffs)
   REAL, DIMENSION(:), INTENT(IN) :: coeffs,x
   REAL, DIMENSION(size(x)) :: poly_rrv
   INTEGER :: i,n,m
   m=size(coeffs)
   n=size(x)
   if (m <= 0) then
      poly_rrv=0.0
   else if (m < n .or. m < NPAR_POLY) then
      poly_rrv=coeffs(m)
      do i=m-1,1,-1
         poly_rrv=x*poly_rrv+coeffs(i)
      end do
   else
      do i=1,n
         poly_rrv(i)=poly_rr(x(i),coeffs)
      end do
   end if
   END FUNCTION poly_rrv
!BL
   FUNCTION poly_msk_rrv(x,coeffs,mask)
   REAL, DIMENSION(:), INTENT(IN) :: coeffs,x
   LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
   REAL, DIMENSION(size(x)) :: poly_msk_rrv
   poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0)
   END FUNCTION poly_msk_rrv
!BL
   RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
   REAL, DIMENSION(:), INTENT(IN) :: a
   REAL, INTENT(IN) :: b
   REAL, DIMENSION(size(a)) :: u
   INTEGER :: n,j
   n=size(a)
   if (n <= 0) RETURN
   u(1)=a(1)
   if (n < NPAR_POLYTERM) then
      do j=2,n
         u(j)=a(j)+b*u(j-1)
      end do
   else
      u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
      u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
   end if
   END FUNCTION poly_term_rr
!BL
   RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
   COMPLEX, DIMENSION(:), INTENT(IN) :: a
   COMPLEX, INTENT(IN) :: b
   COMPLEX, DIMENSION(size(a)) :: u
   INTEGER :: n,j
   n=size(a)
   if (n <= 0) RETURN
   u(1)=a(1)
   if (n < NPAR_POLYTERM) then
      do j=2,n
         u(j)=a(j)+b*u(j-1)
      end do
   else
      u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
      u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
   end if
   END FUNCTION poly_term_cc
!BL
!BL
   FUNCTION zroots_unity(n,nn)
   INTEGER, INTENT(IN) :: n,nn
   COMPLEX, DIMENSION(nn) :: zroots_unity
   INTEGER :: k
   REAL :: theta
   zroots_unity(1)=1.0
   theta=TWOPI/n
   k=1
   do
      if (k >= nn) exit
      zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta))
      zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
         zroots_unity(2:min(k,nn-k))
      k=2*k
   end do
   END FUNCTION zroots_unity
!BL
   FUNCTION outerprod_r(a,b)
   REAL, DIMENSION(:), INTENT(IN) :: a,b
   REAL, DIMENSION(size(a),size(b)) :: outerprod_r
   outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerprod_r
!BL
   FUNCTION outerdiv(a,b)
   REAL, DIMENSION(:), INTENT(IN) :: a,b
   REAL, DIMENSION(size(a),size(b)) :: outerdiv
   outerdiv = spread(a,dim=2,ncopies=size(b)) / &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerdiv
!BL
   FUNCTION outersum(a,b)
   REAL, DIMENSION(:), INTENT(IN) :: a,b
   REAL, DIMENSION(size(a),size(b)) :: outersum
   outersum = spread(a,dim=2,ncopies=size(b)) + &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outersum
!BL
   FUNCTION outerdiff_r(a,b)
   REAL, DIMENSION(:), INTENT(IN) :: a,b
   REAL, DIMENSION(size(a),size(b)) :: outerdiff_r
   outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerdiff_r
!BL
   FUNCTION outerdiff_i(a,b)
   INTEGER, DIMENSION(:), INTENT(IN) :: a,b
   INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
   outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerdiff_i
!BL
   FUNCTION outerand(a,b)
   LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
   LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
   outerand = spread(a,dim=2,ncopies=size(b)) .and. &
      spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerand
!BL
   SUBROUTINE scatter_add_r(dest,source,dest_index)
   REAL, DIMENSION(:), INTENT(OUT) :: dest
   REAL, DIMENSION(:), INTENT(IN) :: source
   INTEGER, DIMENSION(:), INTENT(IN) :: dest_index
   INTEGER :: m,n,j,i
   n=assert_eq2(size(source),size(dest_index),'scatter_add_r')
   m=size(dest)
   do j=1,n
      i=dest_index(j)
      if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
   end do
   END SUBROUTINE scatter_add_r
   SUBROUTINE scatter_max_r(dest,source,dest_index)
   REAL, DIMENSION(:), INTENT(OUT) :: dest
   REAL, DIMENSION(:), INTENT(IN) :: source
   INTEGER, DIMENSION(:), INTENT(IN) :: dest_index
   INTEGER :: m,n,j,i
   n=assert_eq2(size(source),size(dest_index),'scatter_max_r')
   m=size(dest)
   do j=1,n
      i=dest_index(j)
      if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
   end do
   END SUBROUTINE scatter_max_r
!BL
   SUBROUTINE diagadd_rv(mat,diag)
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   REAL, DIMENSION(:), INTENT(IN) :: diag
   INTEGER :: j,n
   n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
   do j=1,n
      mat(j,j)=mat(j,j)+diag(j)
   end do
   END SUBROUTINE diagadd_rv
!BL
   SUBROUTINE diagadd_r(mat,diag)
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   REAL, INTENT(IN) :: diag
   INTEGER :: j,n
   n = min(size(mat,1),size(mat,2))
   do j=1,n
      mat(j,j)=mat(j,j)+diag
   end do
   END SUBROUTINE diagadd_r
!BL
   SUBROUTINE diagmult_rv(mat,diag)
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   REAL, DIMENSION(:), INTENT(IN) :: diag
   INTEGER :: j,n
   n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
   do j=1,n
      mat(j,j)=mat(j,j)*diag(j)
   end do
   END SUBROUTINE diagmult_rv
!BL
   SUBROUTINE diagmult_r(mat,diag)
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   REAL, INTENT(IN) :: diag
   INTEGER :: j,n
   n = min(size(mat,1),size(mat,2))
   do j=1,n
      mat(j,j)=mat(j,j)*diag
   end do
   END SUBROUTINE diagmult_r
!BL
   FUNCTION get_diag_rv(mat)
   REAL, DIMENSION(:,:), INTENT(IN) :: mat
   REAL, DIMENSION(size(mat,1)) :: get_diag_rv
   INTEGER :: j
   j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
   do j=1,size(mat,1)
      get_diag_rv(j)=mat(j,j)
   end do
   END FUNCTION get_diag_rv
!BL
   SUBROUTINE put_diag_rv(diagv,mat)
   REAL, DIMENSION(:), INTENT(IN) :: diagv
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   INTEGER :: j,n
   n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
   do j=1,n
      mat(j,j)=diagv(j)
   end do
   END SUBROUTINE put_diag_rv
!BL
   SUBROUTINE put_diag_r(scal,mat)
   REAL, INTENT(IN) :: scal
   REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
   INTEGER :: j,n
   n = min(size(mat,1),size(mat,2))
   do j=1,n
      mat(j,j)=scal
   end do
   END SUBROUTINE put_diag_r
!BL
   SUBROUTINE unit_matrix(mat)
   REAL, DIMENSION(:,:), INTENT(OUT) :: mat
   INTEGER :: i,n
   n=min(size(mat,1),size(mat,2))
   mat(:,:)=0.0
   do i=1,n
      mat(i,i)=1.0
   end do
   END SUBROUTINE unit_matrix
!BL
   FUNCTION upper_triangle(j,k,extra)
   INTEGER, INTENT(IN) :: j,k
   INTEGER, OPTIONAL, INTENT(IN) :: extra
   LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
   INTEGER :: n
   n=0
   if (present(extra)) n=extra
   upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
   END FUNCTION upper_triangle
!BL
   FUNCTION lower_triangle(j,k,extra)
   INTEGER, INTENT(IN) :: j,k
   INTEGER, OPTIONAL, INTENT(IN) :: extra
   LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
   INTEGER :: n
   n=0
   if (present(extra)) n=extra
   lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
   END FUNCTION lower_triangle
!BL
   FUNCTION vabs(v)
   REAL, DIMENSION(:), INTENT(IN) :: v
   REAL :: vabs
   vabs=sqrt(dot_product(v,v))
   END FUNCTION vabs

PURE ELEMENTAL FUNCTION gammln(xx)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  REAL             :: gammln,xx
  INTEGER          :: j
  REAL             :: ser,stp,tmp,x,y,cof(6)
  PARAMETER(stp=  2.5066282746310005  )
  PARAMETER(cof= (/           &
       76.18009172947146    , &
       -86.50532032941677   , &
       24.01409824083091    , &
       -1.231739572450155   , &
       .1208650973866179e-2 , &
       -.5395239384953e-5    /)       )
  !--------------------------------------------------------------------------
  INTENT(IN) :: xx
  !--------------------------------------------------------------------------

  x   = xx
  y   = x
  tmp = x+5.5d0
  tmp = (x+0.5d0)*LOG(tmp)-tmp
  ser = 1.000000000190015d0

  DO  j=1,6
     y  = y+1.d0
     ser= ser+cof(j)/y
  END DO

  gammln=tmp+LOG(stp*ser/x)

  RETURN

END FUNCTION gammln


PURE SUBROUTINE gauleg(x1,x2,x,w,n)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER     ::  n
  REAL        :: x1,x2,x(n),w(n)
  REAL        :: EPS
  INTEGER     :: i,j,m
  REAL        :: p1,p2,p3,pp,xl,xm,z,z1
  !--------------------------------------------------------------------------
  INTENT(IN)  :: x1,x2,n
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  PARAMETER (EPS=3.E-14)
  !--------------------------------------------------------------------------

  m  = (n+1)/2
  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)

  DO i=1,m
     z = COS(3.141592654*(i-.25)/(n+.5))
1 CONTINUE
     p1 = 1.
     p2 = 0.
     DO j = 1,n
        p3 = p2
        p2 = p1
        p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j
     END DO
     pp = n*(z*p1-p2)/(z*z-1.)
     z1 = z
     z  = z1-p1/pp
     IF(ABS(z-z1).GT.EPS)GOTO 1
     x(i)    = xm-xl*z
     x(n+1-i)= xm+xl*z
     w(i)    = 2.*xl/((1.-z*z)*pp*pp)
     w(n+1-i)= w(i)
  END DO
  RETURN
END SUBROUTINE gauleg

SUBROUTINE gaulag(x,w,alf)
    IMPLICIT NONE
    REAL, INTENT(IN) :: alf
    REAL, DIMENSION(:), INTENT(OUT) :: x,w
    REAL, PARAMETER :: EPS=3.0e-13 
    INTEGER :: its,j
    INTEGER :: n
    INTEGER, PARAMETER :: MAXIT=10
    REAL :: anu
    REAL, PARAMETER :: C1=9.084064e-01, C2=5.214976e-02, C3=2.579930e-03,C4=3.986126e-03
    REAL, DIMENSION(size(x)) :: rhs,r2,r3,theta
    REAL, DIMENSION(size(x)) :: p1,p2,p3,pp,z,z1
    LOGICAL, DIMENSION(size(x)) :: unfinished
    n=assert_eq(size(x),size(w),'gaulag')
    anu=4.0*n+2.0*alf+2.0
    !
    rhs=arth(4*n-1,-4,n)*PI/anu
    !
    r3=rhs**(1.0/3.0)
    r2=r3**2
    theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
    z=anu*cos(theta)**2
    unfinished=.true.
    do its=1,MAXIT
        where (unfinished)
            p1=1.0
            p2=0.0
        end where
        do j=1,n
            where (unfinished)
                p3=p2
                p2=p1
                p1=((2.0*j-1.0+alf-z)*p2-(j-1.0+alf)*p3)/j
            end where
        end do
        where (unfinished)
            pp=(n*p1-(n+alf)*p2)/z
            z1=z
            z=z1-p1/pp
            unfinished=(abs(z-z1) > EPS*z)
        end where
        if (.not. any(unfinished)) exit
    end do
    if (its == MAXIT+1) call nrerror('too many iterations in gaulag')
    x=z
    w=-exp(gammln(alf+n)-gammln(real(n)))/(pp*n*p2)
END SUBROUTINE gaulag


SUBROUTINE gaujac(x,w,n,alf,bet)
  !--------------------------------------------------------------------------
  INTEGER ::  n,MAXIT
  REAL    ::  alf,bet,w(n),x(n)
  REAL    :: EPS
  INTEGER :: i,its,j,pr
  REAL    :: alfbet,an,bn,r1,r2,r3
  REAL    :: a,b,c,p1,p2,p3,pp,temp,z,z1
  REAL    :: test
  !--------------------------------------------------------------------------
  INTENT(IN)  :: n, alf, bet
  INTENT(OUT) :: x,w
  !--------------------------------------------------------------------------
  test = 1. 
  MAXIT=50
  pr = PRECISION(test)
  eps = 0.1**REAL(pr)
  DO i=1,n
     IF(i.EQ.1)THEN
        an=alf/n
        bn=bet/n
        r1=(1.+alf)*(2.78/(4.+n*n)+.768*an/n)
        r2=1.+1.48*an+.96*bn+.452*an*an+.83*an*bn
        z=1.-r1/r2
     ELSE IF(i.EQ.2)THEN
        r1=(4.1+alf)/((1.+alf)*(1.+.156*alf))
        r2=1.+.06*(n-8.)*(1.+.12*alf)/n
        r3=1.+.012*bet*(1.+.25*ABS(alf))/n
        z=z-(1.-z)*r1*r2*r3
     ELSE IF(i.EQ.3)THEN
        r1=(1.67+.28*alf)/(1.+.37*alf)
        r2=1.+.22*(n-8.)/n
        r3=1.+8.*bet/((6.28+bet)*n*n)
        z=z-(x(1)-z)*r1*r2*r3
     ELSE IF(i.EQ.n-1)THEN
        r1=(1.+.235*bet)/(.766+.119*bet)
        r2=1./(1.+.639*(n-4.)/(1.+.71*(n-4.)))
        r3=1./(1.+20.*alf/((7.5+alf)*n*n))
        z=z+(z-x(n-3))*r1*r2*r3
     ELSE IF(i.EQ.n)THEN
        r1=(1.+.37*bet)/(1.67+.28*bet)
        r2=1./(1.+.22*(n-8.)/n)
        r3=1./(1.+8.*alf/((6.28+alf)*n*n))
        z=z+(z-x(n-2))*r1*r2*r3
     ELSE
        z=3.*x(i-1)-3.*x(i-2)+x(i-3)
     ENDIF
     alfbet=alf+bet
     DO its=1,MAXIT
        temp=2.0+alfbet
        p1=(alf-bet+temp*z)/2.0
        p2=1.0
        DO j=2,n
           p3=p2
           p2=p1
           temp=2*j+alfbet
           a=2*j*(j+alfbet)*(temp-2.0)
           b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z)
           c=2.0*(j-1+alf)*(j-1+bet)*temp
           p1=(b*p2-c*p3)/a
        END DO
        pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z))
        z1=z
        z=z1-p1/pp
        IF(ABS(z-z1).LE.EPS) GOTO 1
     END DO

     stop 'too many iterations in gaujac'
1    x(i)=z
     w(i)=EXP(gammln(alf+n)+gammln(bet+n)-gammln(n+1.)-gammln(n+alfbet+1.))*temp*2.**alfbet/(pp*p2)
  END DO
  RETURN
END SUBROUTINE gaujac

SUBROUTINE polint(xa,ya,x,y,dy)
!	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL, DIMENSION(:), INTENT(IN) :: xa,ya
	REAL, INTENT(IN) :: x
	REAL, INTENT(OUT) :: y,dy
	INTEGER :: m,n,ns
	REAL, DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))
	y=ya(ns)
	ns=ns-1
	do m=1,n-1
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)
		if (any(den(1:n-m) == 0.0)) &
			call nrerror('polint: calculation failure')
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m)
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
END SUBROUTINE polint

SUBROUTINE ratint(xa,ya,x,y,dy)
	IMPLICIT NONE
	REAL, DIMENSION(:), INTENT(IN) :: xa,ya
	REAL, INTENT(IN) :: x
	REAL, INTENT(OUT) :: y,dy
	INTEGER :: m,n,ns
	REAL, DIMENSION(size(xa)) :: c,d,dd,h,t
	REAL, PARAMETER :: TINY=1.0e-25
	n=assert_eq(size(xa),size(ya),'ratint')
	h=xa-x
	ns=iminloc(abs(h))
	y=ya(ns)
	if (x == xa(ns)) then
		dy=0.0
		RETURN
	end if
	c=ya
	d=ya+TINY
	ns=ns-1
	do m=1,n-1
		t(1:n-m)=(xa(1:n-m)-x)*d(1:n-m)/h(1+m:n)
		dd(1:n-m)=t(1:n-m)-c(2:n-m+1)
		if (any(dd(1:n-m) == 0.0)) &
			call nrerror('failure in ratint')
		dd(1:n-m)=(c(2:n-m+1)-d(1:n-m))/dd(1:n-m)
		d(1:n-m)=c(2:n-m+1)*dd(1:n-m)
		c(1:n-m)=t(1:n-m)*dd(1:n-m)
		if (2*ns < n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		end if
		y=y+dy
	end do
END SUBROUTINE ratint


  SUBROUTINE svbksb(u,w,v,b,x)
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  REAL, DIMENSION(:,:)    :: u,v
  REAL, DIMENSION(:)      :: w,b
  REAL,DIMENSION(:)       :: x
  INTEGER                 :: mdum,ndum
  REAL, DIMENSION(SIZE(x)):: tmp
  !--------------------------------------------------------------------------
  INTENT(IN)  :: u,v,w,b
  INTENT(OUT) :: x
  !--------------------------------------------------------------------------

  mdum=assert_eq(SIZE(u,1),SIZE(b),'svbksb: mdum')

  ndum=assert_eq((/SIZE(u,2),SIZE(v,1),SIZE(v,2),SIZE(w),SIZE(x)/),&
                 'svbksb: ndum')

  WHERE (ABS(w).GE.1e-10)
     tmp=MATMUL(b,u)/w
  ELSEWHERE
     tmp=0.0
  END WHERE

  x=MATMUL(v,tmp)

END SUBROUTINE svbksb


SUBROUTINE svdcmp(a,w,v)
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  REAL, DIMENSION(:,:)       :: a
  REAL, DIMENSION(:)         :: w
  REAL, DIMENSION(:,:)       :: v
  INTEGER                    :: i,its,j,k,l,m,n,nm
  REAL                       :: anorm,c,f,g,h,s,scale,x,y,z
  REAL, DIMENSION(size(a,1)) :: tempm
  REAL, DIMENSION(size(a,2)) :: rv1,tempn
  !--------------------------------------------------------------------------
  INTENT(INOUT) :: a
  INTENT(OUT)   :: w, v
  !--------------------------------------------------------------------------

  m=size(a,1)
  n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp')
  g=0.0
  scale=0.0
  do i=1,n
     l=i+1
     rv1(i)=scale*g
     g=0.0
     scale=0.0
     if (i <= m) then
        scale=sum(abs(a(i:m,i)))
        if (ABS(scale).GT.1e-12) then
           a(i:m,i)=a(i:m,i)/scale
           s=dot_product(a(i:m,i),a(i:m,i))
           f=a(i,i)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,i)=f-g
           tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
           a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
           a(i:m,i)=scale*a(i:m,i)
        end if
     end if
     w(i)=scale*g
     g=0.0
     scale=0.0
     if ((i <= m) .and. (i /= n)) then
        scale=sum(abs(a(i,l:n)))
        if (abs(scale).GT.1e-12) then
           a(i,l:n)=a(i,l:n)/scale
           s=dot_product(a(i,l:n),a(i,l:n))
           f=a(i,l)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,l)=f-g
           rv1(l:n)=a(i,l:n)/h
           tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
           a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
           a(i,l:n)=scale*a(i,l:n)
        end if
     end if
  end do
  anorm=maxval(abs(w)+abs(rv1))
  do i=n,1,-1
     if (i < n) then
        if (abs(g).GT.1e-12) then
           v(l:n,i)=(a(i,l:n)/a(i,l))/g
           tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
           v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
        end if
        v(i,l:n)=0.0
        v(l:n,i)=0.0
     end if
     v(i,i)=1.0
     g=rv1(i)
     l=i
  end do
  do i=min(m,n),1,-1
     l=i+1
     g=w(i)
     a(i,l:n)=0.0
     if (ABS(g).GT.1e-14) then
        g=1.0/g
        tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
        a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
        a(i:m,i)=a(i:m,i)*g
     else
        a(i:m,i)=0.0
     end if
     a(i,i)=a(i,i)+1.0
  end do
  do k=n,1,-1
     do its=1,30
        do l=k,1,-1
           nm=l-1
           if ((abs(rv1(l))+anorm) == anorm) exit
           if ((abs(w(nm))+anorm) == anorm) then
              c=0.0
              s=1.0
              do i=l,k
                 f=s*rv1(i)
                 rv1(i)=c*rv1(i)
                 if ((abs(f)+anorm) == anorm) exit
                 g=w(i)
                 h=pythag(f,g)
                 w(i)=h
                 h=1.0/h
                 c= (g*h)
                 s=-(f*h)
                 tempm(1:m)=a(1:m,nm)
                 a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                 a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
              end do
              exit
           end if
        end do
        z=w(k)
        if (l == k) then
           if (z < 0.0) then
              w(k)=-z
              v(1:n,k)=-v(1:n,k)
           end if
           exit
        end if
        if (its == 30) call nrerror('svdcmp: no convergence in svdcmp')
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
        g=pythag(f,1.0)
        f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c=1.0
        s=1.0
        do j=l,nm
           i=j+1
           g=rv1(i)
           y=w(i)
           h=s*g
           g=c*g
           z=pythag(f,h)
           rv1(j)=z
           c=f/z
           s=h/z
           f= (x*c)+(g*s)
           g=-(x*s)+(g*c)
           h=y*s
           y=y*c
           tempn(1:n)=v(1:n,j)
           v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
           v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
           z=pythag(f,h)
           w(j)=z
           if (abs(z).GT.1e-12) then
              z=1.0/z
              c=f*z
              s=h*z
           end if
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           tempm(1:m)=a(1:m,j)
           a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
           a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
        end do
        rv1(l)=0.0
        rv1(k)=f
        w(k)=x
     end do
  end do
END SUBROUTINE svdcmp


SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER                    :: j, k, m, n, NN
  REAL                       :: x1a(m), x2a(n), y2a(m,n), ya(m,n)
  PARAMETER (NN=1000000)
  REAL                       :: y2tmp(NN), ytmp(NN)

      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return

END SUBROUTINE splie2


SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER                    :: j, k, m, n, NN
  REAL                       :: x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
  PARAMETER (NN=1000000)
  REAL                       :: y2tmp(NN),ytmp(NN),yytmp(NN)

      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        CALL splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      CALL spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      CALL splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      
END SUBROUTINE splin2


SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------      
  INTEGER                    :: i, k, n, NMAX
  REAL                       :: yp1,ypn,x(n),y(n),y2(n)
  REAL                       :: p,qn,sig,un
  PARAMETER (NMAX=500000)
  REAL                       :: u(NMAX)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      
END SUBROUTINE spline


SUBROUTINE splint(xa,ya,y2a,n,x,y)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------    
  INTEGER                    :: n,k,khi,klo
  REAL                       :: a,b,h,x,y,xa(n),y2a(n),ya(n)

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
        print*,'bad xa input in splint'
        stop 
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      
END SUBROUTINE splint


SUBROUTINE hpsortr(n,ra) 
      INTEGER n
      REAL ra(n)
      INTEGER i,ir,j,l
      REAL rra
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
      goto 10
END SUBROUTINE 

SUBROUTINE hpsorti(n,ra) 
      INTEGER n
      INTEGER ra(n)
      INTEGER i,ir,j,l
      INTEGER rra
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
      goto 10
END SUBROUTINE 

SUBROUTINE hpsortia(n,m,ra) 
      INTEGER n,m
      INTEGER ra(m,n)
      INTEGER i,ir,j,l
      INTEGER rra(m) 
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(:,l)
        else
          rra=ra(:,ir)
          ra(:,ir)=ra(:,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(:,1)=rra
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(1,j).lt.ra(1,j+1))j=j+1
          endif
          if(rra(1).lt.ra(1,j))then
            ra(:,i)=ra(:,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(:,i)=rra 
      goto 10
END SUBROUTINE 

SUBROUTINE gaussj(a,b)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER, DIMENSION(size(a,1)) :: ipiv,indxr,indxc
    LOGICAL, DIMENSION(size(a,1)) :: lpiv
    REAL :: pivinv
    REAL, DIMENSION(size(a,1)) :: dumc
    INTEGER, TARGET :: irc(2)
    INTEGER :: i,l,n
    INTEGER, POINTER :: irow,icol
    n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
        lpiv = (ipiv == 0)
        irc=maxloc(abs(a),outerand(lpiv,lpiv))
        ipiv(icol)=ipiv(icol)+1
        if (ipiv(icol) > 1) then
            call nrerror('gaussj: singular matrix (1)')
        endif
        if (irow /= icol) then
            call swap(a(irow,:),a(icol,:))
            call swap(b(irow,:),b(icol,:))
        end if
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol) == 0.0) then
            call nrerror('gaussj: singular matrix (2)')
        endif
        pivinv=1.0/a(icol,icol)
        a(icol,icol)=1.0
        a(icol,:)=a(icol,:)*pivinv
        b(icol,:)=b(icol,:)*pivinv
        dumc=a(:,icol)
        a(:,icol)=0.0
        a(icol,icol)=pivinv
        if (icol .gt. 1) then
           a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
           b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
        endif
        if (icol .lt. size(a,1)) then
           a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
           b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
        endif
    end do
    do l=n,1,-1
        call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
END SUBROUTINE gaussj

SUBROUTINE rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
!
      integer n,nm,is1,is2,ierr,matz
      real a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
!        iv1  and  fv1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
!
   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)
   50 return
      
END SUBROUTINE rg      

SUBROUTINE balanc(nm,n,a,low,igh,scale)
!
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision a(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
!
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j),      j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      radix = 16.0d0
!
      b2 = radix * radix
      k = 1
      l = n
      go to 100
!     .......... in-line procedure for row and
!                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
!
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
!
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
!
   50 go to (80,130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
!
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.0d0) go to 120
  110    continue
!
         m = l
         iexc = 1
         go to 20
  120 continue
!
      go to 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1
!
  140 do 170 j = k, l
!
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.0d0) go to 170
  150    continue
!
         m = k
         iexc = 2
         go to 20
  170 continue
!     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
!
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(a(j,i))
            r = r + dabs(a(i,j))
  200    continue
!     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
!     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
!
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
!
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
!
  270 continue
!
      if (noconv) go to 190
!
  280 low = k
      igh = l
      return

END SUBROUTINE balanc

SUBROUTINE balbak(nm,n,low,igh,scale,m,z)
!
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),z(nm,m)
      double precision s
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
!
      do 110 i = low, igh
         s = scale(i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
  100    z(i,j) = z(i,j) * s
!
  110 continue
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
!
         do 130 j = 1, m
            s = z(i,j)
            z(i,j) = z(k,j)
            z(k,j) = s
  130    continue
!
  140 continue
!
  200 return
END SUBROUTINE balbak

SUBROUTINE cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
      double precision s,ars,ais,brs,bis
      s = dabs(br) + dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      
END SUBROUTINE cdiv 

SUBROUTINE elmhes(nm,n,low,igh,a,int)
!
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision a(nm,n)
      double precision x,y
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.0d0
         i = m
!
         do 100 j = m, igh
            if (dabs(a(j,mm1)) .le. dabs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
!
         int(m) = i
         if (i .eq. m) go to 130
!     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
!
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
!     .......... end interchange ..........
  130    if (x .eq. 0.0d0) go to 180
         mp1 = m + 1
!
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.0d0) go to 160
            y = y / x
            a(i,mm1) = y
!
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
!
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
!
  160    continue
!
  180 continue
!
  200 return

END SUBROUTINE elmhes
      
SUBROUTINE eltran(nm,n,low,igh,a,int,z)
!
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),z(nm,n)
      integer int(igh)
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper hessenberg form by  elmhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  elmhes.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z to identity matrix ..........
      do 80 j = 1, n
!
         do 60 i = 1, n
   60    z(i,j) = 0.0d0
!
         z(j,j) = 1.0d0
   80 continue
!
      kl = igh - low - 1
      if (kl .lt. 1) go to 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      do 140 mm = 1, kl
         mp = igh - mm
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    z(i,mp) = a(i,mp-1)
!
         i = int(mp)
         if (i .eq. mp) go to 140
!
         do 130 j = mp, igh
            z(mp,j) = z(i,j)
            z(i,j) = 0.0d0
  130    continue
!
         z(i,mp) = 1.0d0
  140 continue
!
  200 return

END SUBROUTINE eltran

SUBROUTINE hqr(nm,n,low,igh,h,wr,wi,ierr)
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.  information about
!          the transformations used in the reduction to hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated september 1989.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0d0
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + dabs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
!
      en = igh
      t = 0.0d0
      itn = 30*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = dabs(h(l-1,l-1)) + dabs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + dabs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = dabs(h(en,na)) + dabs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = dabs(p) + dabs(q) + dabs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = dabs(p)*(dabs(h(m-1,m-1)) + dabs(zz) + dabs(h(m+1,m+1)))
         tst2 = tst1 + dabs(h(m,m-1))*(dabs(q) + dabs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = dabs(p) + dabs(q) + dabs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.0d0
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt(dabs(q))
      x = x + t
      if (q .lt. 0.0d0) go to 320
!     .......... real pair ..........
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
END SUBROUTINE hqr


SUBROUTINE hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
!
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,igh,itn,its,low,mp2,enm2,ierr
      real h(nm,n),wr(n),wi(n),z(nm,n)
      real p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0d0
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.0d0
   50 continue
!
      en = igh
      t = 0.0d0
      itn = 50*n
!     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.0d0) s = norm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
!     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
!     .......... form exceptional shift ..........
      t = t + x
!
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75d0 * s
      y = x
      w = -0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = 0.0d0
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.0d0
  160 continue
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0d0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0d0) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
!     .......... row modification ..........
         do 200 j = k, n
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 210 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
!     .......... accumulate transformations ..........
         do 220 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
  220    continue
         go to 255
  225    continue
!     .......... row modification ..........
         do 230 j = k, n
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
         j = min0(en,k+3)
!     .......... column modification ..........
         do 240 i = 1, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
!     .......... accumulate transformations ..........
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
            z(i,k) = z(i,k) - p
            z(i,k+1) = z(i,k+1) - p * q
            z(i,k+2) = z(i,k+2) - p * r
  250    continue
  255    continue
!
  260 continue
!
      go to 70
!     .......... one root found ..........
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
      go to 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = sqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0d0) go to 320
!     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0d0) wr(en) = x - w / zz
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = sqrt(p*p+q*q)
      p = p / r
      q = q / r
!     .......... row modification ..........
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
!     .......... column modification ..........
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
!     .......... accumulate transformations ..........
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
!
      go to 330
!     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 if (norm .eq. 0.0d0) go to 1001
!     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q) 710, 600, 800
!     .......... real vector ..........
  600    m = en
         h(en,en) = 1.0d0
         if (na .eq. 0) go to 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = 0.0d0
!
            do 610 j = m, en
  610       r = r + h(i,j) * h(j,en)
!
            if (wi(i) .ge. 0.0d0) go to 630
            zz = w
            s = r
            go to 700
  630       m = i
            if (wi(i) .ne. 0.0d0) go to 640
            t = w
            if (t .ne. 0.0d0) go to 635
               tst1 = norm
               t = tst1
  632          t = 0.01d0 * t
               tst2 = norm + t
               if (tst2 .gt. tst1) go to 632
  635       h(i,en) = -r / t
            go to 680
!     .......... solve real equations ..........
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) go to 650
            h(i+1,en) = (-r - w * t) / x
            go to 680
  650       h(i+1,en) = (-s - y * t) / zz
!
!     .......... overflow control ..........
  680       t = abs(h(i,en))
            if (t .eq. 0.0d0) go to 700
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 700
            do 690 j = i, en
               h(j,en) = h(j,en)/t
  690       continue
!
  700    continue
!     .......... end real vector ..........
         go to 800
!     .......... complex vector ..........
  710    m = na
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
         if (abs(h(en,na)) .le. abs(h(na,en))) go to 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         go to 730
  720    call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
  730    h(en,na) = 0.0d0
         h(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
         do 795 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0d0
            sa = 0.0d0
!
            do 760 j = m, en
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. 0.0d0) go to 770
            zz = w
            r = ra
            s = sa
            go to 795
  770       m = i
            if (wi(i) .ne. 0.0d0) go to 780
            call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
            go to 790
!     .......... solve complex equations ..........
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * 2.0d0 * q
            if (vr .ne. 0.0d0 .or. vi .ne. 0.0d0) go to 784
               tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
               vr = tst1
  783          vr = 0.01d0 * vr
               tst2 = tst1 + vr
               if (tst2 .gt. tst1) go to 783
  784       call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))
            if (abs(x) .le. abs(zz) + abs(q)) go to 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            go to 790
  785       call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en))
!
!     .......... overflow control ..........
  790       t = max(abs(h(i,na)), abs(h(i,en)))
            if (t .eq. 0.0d0) go to 795
            tst1 = t
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 795
            do 792 j = i, en
               h(j,na) = h(j,na)/t
               h(j,en) = h(j,en)/t
  792       continue
!
  795    continue
!     .......... end complex vector ..........
  800 continue
!     .......... end back substitution.
!                vectors of isolated roots ..........
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 840
!
         do 820 j = i, n
  820    z(i,j) = h(i,j)
!
  840 continue
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
!
         do 880 i = low, igh
            zz = 0.0d0
!
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
!
            z(i,j) = zz
  880 continue
!
      go to 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
 
END SUBROUTINE hqr2

REAL FUNCTION expint(n,x)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  INTEGER  :: n
  REAL     :: x
  INTEGER  :: MAXIT
  REAL     :: EPS,BIG
  INTEGER  :: i,nm1
  REAL     :: a,b,c,d,del,fact,h
  !--------------------------------------------------------------------------
  INTENT(IN):: n,x
  !--------------------------------------------------------------------------
  PARAMETER(MAXIT=100)
  PARAMETER(EPS=epsilon(x))
  PARAMETER(BIG=huge(x)*EPS)
  !--------------------------------------------------------------------------

   call assert(n >= 0, x >= 0.0, (x > 0.0 .or. n > 1), &
      'expint args')
   if (n == 0) then
      expint=exp(-x)/x
      RETURN
   end if
   nm1=n-1
   if (x == 0.0) then
      expint=1.0/nm1
   else if (x > 1.0) then
      b=x+n
      c=BIG
      d=1.0/b
      h=d
      do i=1,MAXIT
         a=-i*(nm1+i)
         b=b+2.0
         d=1.0/(a*d+b)
         c=b+a/c
         del=c*d
         h=h*del
         if (abs(del-1.0) <= EPS) exit
      end do
      if (i > MAXIT) call nrerror('expint: continued fraction failed')
      expint=h*exp(-x)
   else
      if (nm1 /= 0) then
         expint=1.0/nm1
      else
         expint=-log(x)-NEULER
      end if
      fact=1.0
      do i=1,MAXIT
         fact=-fact*x/i
         if (i /= nm1) then
            del=-fact/(i-nm1)
         else
            del=fact*(-log(x)-NEULER+sum(1.0/arth(1,1,nm1)))
         end if
         expint=expint+del
         if (abs(del) < abs(expint)*EPS) exit
      end do
      if (i > MAXIT) call nrerror('expint: series failed')
   end if
   END FUNCTION expint

FUNCTION ei(x)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  REAL    :: x
  REAL    :: ei
  INTEGER :: MAXIT
  REAL    :: EPS,FPMIN
  INTEGER :: k
  REAL    :: fact,prev,sm,term
  !--------------------------------------------------------------------------
  INTENT(IN) :: x
  !--------------------------------------------------------------------------
  PARAMETER(MAXIT=100)
  PARAMETER(EPS=epsilon(x))
  PARAMETER(FPMIN=tiny(x)/EPS)
  !--------------------------------------------------------------------------

   call assert(x > 0.0, 'ei arg')
   if (x < FPMIN) then
      ei=log(x)+NEULER
   else if (x <= -log(EPS)) then
      sm=0.0
      fact=1.0
      do k=1,MAXIT
         fact=fact*x/k
         term=fact/k
         sm=sm+term
         if (term < EPS*sm) exit
      end do
      if (k > MAXIT) call nrerror('series failed in ei')
      ei=sm+log(x)+NEULER
   else
      sm=0.0
      term=1.0
      do k=1,MAXIT
         prev=term
         term=term*k/x
         if (term < EPS) exit
         if (term < prev) then
            sm=sm+term
         else
            sm=sm-prev
            exit
         end if
      end do
      if (k > MAXIT) call nrerror('asymptotic failed in ei')
      ei=exp(x)*(1.0+sm)/x
   end if
   END FUNCTION ei

FUNCTION plgndr_s(l,m,x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l,m
  REAL, INTENT(IN) :: x
  REAL :: plgndr_s
  INTEGER :: ll
  REAL :: pll,pmm,pmmp1,somx2
  call assert(m >= 0, m <= l, abs(x) <= 1.0, 'plgndr_s args')
  pmm=1.0
  if (m > 0) then
     somx2=sqrt((1.0-x)*(1.0+x))
     pmm=product(arth(1.0,2.0,m))*somx2**m
     if (mod(m,2) == 1) pmm=-pmm
  end if
  if (l == m) then
     plgndr_s=pmm
  else
     pmmp1=x*(2*m+1)*pmm
     if (l == m+1) then
        plgndr_s=pmmp1
     else
        do ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plgndr_s=pll
     end if
  end if
END FUNCTION plgndr_s


FUNCTION plgndr_v(l,m,x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l,m
  REAL, DIMENSION(:), INTENT(IN) :: x
  REAL, DIMENSION(size(x)) :: plgndr_v
  INTEGER :: ll
  REAL, DIMENSION(size(x)) :: pll,pmm,pmmp1,somx2
  call assert(m >= 0, m <= l, all(abs(x) <= 1.0), 'plgndr_v args')
  pmm=1.0
  if (m > 0) then
     somx2=sqrt((1.0-x)*(1.0+x))
     pmm=product(arth(1.0,2.0,m))*somx2**m
     if (mod(m,2) == 1) pmm=-pmm
  end if
  if (l == m) then
     plgndr_v=pmm
  else
     pmmp1=x*(2*m+1)*pmm
     if (l == m+1) then
        plgndr_v=pmmp1
     else
        do ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plgndr_v=pll
     end if
  end if
END FUNCTION plgndr_v

PURE ELEMENTAL REAL FUNCTION factln(n)
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------
  ! USES gammln
  !--------------------------------------------------------------------------
  INTEGER    :: n
  !--------------------------------------------------------------------------
  INTENT(IN) :: n
  !--------------------------------------------------------------------------
 
  factln=gammln(ABS(n)+1.)

  RETURN

END FUNCTION factln

SUBROUTINE indexx(n,arr,indx)

      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a

      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) THEN
           print *,'NSTACK too small in indexx'
           stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE indexx      

SUBROUTINE golden(ax,bx,cx,f,tol,xmin,fmin)
  REAL ax,bx,cx,tol,xmin,fmin,f,R,C
  EXTERNAL f
  PARAMETER (R=.61803399,C=1.-R)
  REAL f1,f2,x0,x1,x2,x3
  x0=ax
  x3=cx
  if(abs(cx-bx).gt.abs(bx-ax))then
     x1=bx
     x2=bx+C*(cx-bx)
  else
     x2=bx
     x1=bx-C*(bx-ax)
  endif
  f1=f(x1)
  f2=f(x2)
1 if(abs(x3-x0).gt.tol+tol*(abs(x1)+abs(x2)))then
     if(f2.lt.f1)then
        x0=x1
        x1=x2
        x2=R*x1+C*x3
        f1=f2
        f2=f(x2)
     else
        x3=x2
        x2=x1
        x1=R*x2+C*x0
        f2=f1
        f1=f(x1)
     endif
     goto 1
  endif
  if(f1.lt.f2)then
     fmin=f1
     xmin=x1
  else
     fmin=f2
     xmin=x2
  endif
  return
END SUBROUTINE golden

END MODULE recipies_mod 
