SUBROUTINE METRIC ( xc, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
#ifdef GRMHD
  USE typesDef, ONLY: ICType, aom, Mbh, P_eps
#else
  USE typesDef, ONLY : ICType  
#endif
#ifdef TWOPUNCTURES  
	USE TwoPunctures_mod, ONLY : TwoPunctures_Interpolate
#endif 
#ifdef RNSID 
    USE RNSID_mod, ONLY : RNSID_Interpolate
#endif
  IMPLICIT NONE 
  !
  REAL, DIMENSION(3), intent(IN) :: xc
  REAL                           :: lapse,gp,gm, phi 
  REAL,dimension(3)              :: shift
  REAL,dimension(3,3)            :: Kex, g_cov, g_contr 

  REAL :: x, y, z, r, z2, r2, aom2, det 
  REAL :: lx, ly, lz, HH, SS, detg, rho, r_hz, sig, theta, a, M    
  REAL :: st, st2, delta, rho2, sigma, zz, V0Z4(54),V070(70) 
  REAL :: transl(3), xGP(3) 
#ifndef GRMHD  
  REAL, PARAMETER :: aom = 0.9, Mbh = 1.0, P_eps = 1e-4    
#endif     
  !
  Kex = 0.0 
  phi = 1.0
  !
  SELECT CASE(ICType)
  CASE('CCZ4Puncture') 
    !   
    r = P_eps + SQRT( SUM(xc**2) ) 
    g_cov = 0.0
    g_cov(1,1) = 1.0 
    g_cov(2,2) = 1.0
    g_cov(3,3) = 1.0 
    g_cov = g_cov * (1.0+1./r)**4 
    g_contr = 0.0 
    g_contr(1,1) = 1./g_cov(1,1)  
    g_contr(2,2) = 1./g_cov(2,2)  
    g_contr(3,3) = 1./g_cov(3,3)  
    lapse = 1.0
    shift = 0.0 
    det = g_cov(1,1)*g_cov(2,2)*g_cov(3,3) 
    phi = det**(-1./6.) 
    ! 
  CASE('CCZ4TwoPunctures') 
    !
    transl = (/ -1.0, 0.0, 0.0 /) 
    xGP = xc - transl 
    !  
#ifdef TWOPUNCTURES  
        CALL TwoPunctures_Interpolate(xGP, V0Z4) 
#else
        PRINT *, ' TwoPunctures not available. Please compile with -DTWOPUNCTURES flag. '
#endif      
    !
    g_cov(1,1) = V0Z4(1) 
    g_cov(1,2) = V0Z4(2) 
    g_cov(2,1) = V0Z4(2) 
    g_cov(1,3) = V0Z4(3) 
    g_cov(3,1) = V0Z4(3) 
    g_cov(2,2) = V0Z4(4) 
    g_cov(2,3) = V0Z4(5) 
    g_cov(3,2) = V0Z4(5) 
    g_cov(3,3) = V0Z4(6) 
    !
    Kex(1,1) = V0Z4(7) 
    Kex(1,2) = V0Z4(8) 
    Kex(2,1) = V0Z4(8) 
    Kex(1,3) = V0Z4(9) 
    Kex(3,1) = V0Z4(9) 
    Kex(2,2) = V0Z4(10) 
    Kex(2,3) = V0Z4(11) 
    Kex(3,2) = V0Z4(11) 
    Kex(3,3) = V0Z4(12) 
    !
    CALL  MatrixInverse3x3(g_cov,g_contr,det)  
    lapse = V0Z4(17) 
    shift = V0Z4(18:20) 
    phi = det**(-1./6.) 
    ! 
  CASE('CCZ4GRMHDAccretion') 
    !
    ! Rotating black hole in Kerr-Schild Cartesian coordinates. See De Felice & Clarke Sect. 11.4
    x    = xc(1)
    y    = xc(2)
    z    = xc(3)
  
    z2   = z**2
    aom2 = aom**2
    r    = MAX( 1e-1, SQRT( (x**2 + y**2 + z**2 - aom2)/2.0 + SQRT(((x**2 + y**2 + z**2 - aom2)/2.0)**2 + z2*aom2)) )  
  
    r2   = r**2
  
    HH = Mbh*r2*r / (r2*r2 + aom2*z2)
    SS = 1.0 + 2.0*HH
    lx = (r*x + aom*y)/(r2 + aom2)  
    ly = (r*y - aom*x)/(r2 + aom2)
    lz = z/r  
  
    lapse   = 1.0/SQRT(SS)
    shift(1) = 2.0*HH/SS*lx         
    shift(2) = 2.0*HH/SS*ly         
    shift(3) = 2.0*HH/SS*lz        
     
    g_cov( 1, 1:3) = (/ 1.0 + 2.0*HH*lx**2, 2.0*HH*lx*ly,        2.0*HH*lx*lz       /)
    g_cov( 2, 1:3) = (/ 2.0*HH*lx*ly,       1.0 + 2.0*HH*ly**2,  2.0*HH*ly*lz       /)
    g_cov( 3, 1:3) = (/ 2.0*HH*lx*lz,       2.0*HH*ly*lz,        1.0 + 2.0*HH*lz**2 /)
  
    g_contr( 1, 1:3) = (/  1.0 + 2.0*HH*ly**2 + 2.0*HH*lz**2,   -2.0*HH*lx*ly,                       -2.0*HH*lx*lz                     /)
    g_contr( 2, 1:3) = (/ -2.0*HH*lx*ly,                         1.0 + 2.0*HH*lx**2 + 2.0*HH*lz**2,  -2.0*HH*ly*lz                     /)
    g_contr( 3, 1:3) = (/ -2.0*HH*lx*lz,                        -2.0*HH*ly*lz,                        1.0 + 2.0*HH*lx**2 + 2.0*HH*ly**2 /)
  
    g_contr = g_contr/SS
  
    gp = SQRT(SS)
    gm = 1.0/gp
    !
    phi = SS**(-1./6.) 
    !
  CASE('CCZ4RNSID')
    !
    transl = (/ 0.0, 0.0, 0.0 /) 
    xGP = xc - transl 
    !  
#ifdef RNSID   
        CALL RNSID_Interpolate(xGP, V070) 
#else
        PRINT *, ' RNSID not available. Please compile with -DRNSID flag. '
        STOP 
#endif      

    !
    g_cov(1,1) = V070(1) 
    g_cov(1,2) = V070(2) 
    g_cov(2,1) = V070(2) 
    g_cov(1,3) = V070(3) 
    g_cov(3,1) = V070(3) 
    g_cov(2,2) = V070(4) 
    g_cov(2,3) = V070(5) 
    g_cov(3,2) = V070(5) 
    g_cov(3,3) = V070(6) 
    !
    Kex(1,1) = V070(7) 
    Kex(1,2) = V070(8) 
    Kex(2,1) = V070(8) 
    Kex(1,3) = V070(9) 
    Kex(3,1) = V070(9) 
    Kex(2,2) = V070(10) 
    Kex(2,3) = V070(11) 
    Kex(3,2) = V070(11) 
    Kex(3,3) = V070(12) 
    !
    CALL  MatrixInverse3x3(g_cov,g_contr,det)  
    lapse = V070(17) 
    shift = V070(18:20) 
    phi = det**(-1./6.)       
    !  
  CASE('CCZ4Kerr2D')  
      !
      ! Rotating black hole in Kerr-Schild spherical coordinates. See Appendix B of Komissarov (2004) MNRAS, 350, 427
      !
      r  = xc(1)
      r2 = r*r
      !    
      st = SIN(xc(2))
      IF (st < 1.e-6) st = 1.
      st2 = st**2
  
      aom2 = aom**2
  
      delta = r2 - 2.0 * Mbh * r + aom2
      rho2  = r2 + aom2*(1.0 - st2)
      sigma = (r2 + aom2)**2 - aom2*delta*st2
      zz    = 2.0*r/rho2
  
      lapse    = 1.0 / sqrt(1.0 + zz)
      shift(1) = zz/(1.0 + zz)       
      shift(2) = 0.0         
      shift(3) = 0.0  
  
      g_cov( 1, 1:3) = (/ 1.0 + zz,             0.0,        -aom*st2*(1.0 + zz)   /)
      g_cov( 2, 1:3) = (/ 0.0,                  rho2,       0.0                   /)
      g_cov( 3, 1:3) = (/ -aom*st2*(1.0 + zz),  0.0,        (sigma/rho2)*st2      /)
      !
      det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
      phi = det**(-1./6.) 
      !
      g_contr( 1, 1:3) = (/  lapse**2 + aom2*st2/rho2,   0.0,                       aom/rho2         /)
      g_contr( 2, 1:3) = (/  0.0,                        1.0/rho2,                  0.0              /)
      g_contr( 3, 1:3) = (/  aom/rho2,                   0.0,                       1.0/(rho2*st2)   /)
  
      gp = rho2*st*sqrt(1.0 + zz)
      gm = 1.0/gp    
      !
  CASE('Z4Kerr2D')
      !
       r      = xc(1)   ! SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2)
       theta  = xc(2)   ! 
       a      = aom 
       M      = Mbh  
       r_hz   = Mbh + SQRT(Mbh**2 - aom**2) 
       !
       rho   = SQRT(r**2 + a**2*COS(theta)**2) 
       Delta = r**2 - 2*Mbh*r+a**2 
       Sig   = (r**2+a**2)**2 - a**2*Delta*SIN(theta)**2  
       zz    = 2.*M*r/rho**2
       !
       g_cov( 1, 1:3) = (/ 1.0+zz,                     0.0,     -a*SIN(theta)**2*(1.0+zz)   /) 
       g_cov( 2, 1:3) = (/ 0.0,                        rho**2,   0.0                        /) 
       g_cov( 3, 1:3) = (/ -a*SIN(theta)**2*(1.0+zz),  0.0,      sig*SIN(theta)**2/rho**2   /)       
       !
       CALL  MatrixInverse3x3(g_cov,g_contr,det)  
       !
       phi = det**(-1./6.)         
       !
       gp = SQRT(det)
       gm = 1.0/gp 
       !
       lapse = 1.0/SQRT(1+zz)                            ! alpha (lapse) 
       ! 
       shift(1) = zz/(1+zz)                              ! beta_1 (shift) 
       shift(2) = 0.0                                    ! beta_2 
       shift(3) = 0.0                                    ! beta_3        
       !
  CASE('GRMHDAccretion','GRMHDTorus','Z4Kerr3D','CCZ4Kerr3D')
      !
#ifndef Spherical
  ! Rotating black hole in Kerr-Schild Cartesian coordinates. See De Felice & Clarke Sect. 11.4
  x    = xc(1)
  y    = xc(2)
  z    = xc(3)
  
  z2   = z**2
  aom2 = aom**2
  r    = SQRT( (x**2 + y**2 + z**2 - aom2)/2.0 + SQRT(((x**2 + y**2 + z**2 - aom2)/2.0)**2 + z2*aom2))
  
  r2   = r**2
  
  HH = Mbh*r2*r / (r2*r2 + aom2*z2)
  SS = 1.0 + 2.0*HH
  lx = (r*x + aom*y)/(r2 + aom2)  
  ly = (r*y - aom*x)/(r2 + aom2)
  lz = z/r  
  
  lapse   = 1.0/SQRT(SS)
  shift(1) = 2.0*HH/SS*lx         
  shift(2) = 2.0*HH/SS*ly         
  shift(3) = 2.0*HH/SS*lz        
     
  g_cov( 1, 1:3) = (/ 1.0 + 2.0*HH*lx**2, 2.0*HH*lx*ly,        2.0*HH*lx*lz       /)
  g_cov( 2, 1:3) = (/ 2.0*HH*lx*ly,       1.0 + 2.0*HH*ly**2,  2.0*HH*ly*lz       /)
  g_cov( 3, 1:3) = (/ 2.0*HH*lx*lz,       2.0*HH*ly*lz,        1.0 + 2.0*HH*lz**2 /)
  
  det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
  phi = det**(-1./6.) 

  g_contr( 1, 1:3) = (/  1.0 + 2.0*HH*ly**2 + 2.0*HH*lz**2,   -2.0*HH*lx*ly,                       -2.0*HH*lx*lz                     /)
  g_contr( 2, 1:3) = (/ -2.0*HH*lx*ly,                         1.0 + 2.0*HH*lx**2 + 2.0*HH*lz**2,  -2.0*HH*ly*lz                     /)
  g_contr( 3, 1:3) = (/ -2.0*HH*lx*lz,                        -2.0*HH*ly*lz,                        1.0 + 2.0*HH*lx**2 + 2.0*HH*ly**2 /)
  
  g_contr = g_contr/SS
  
  gp = SQRT(SS)
  gm = 1.0/gp
  ! 
#else
  ! Rotating black hole in Kerr-Schild spherical coordinates. See Appendix B of Komissarov (2004) MNRAS, 350, 427
  r  = xc(1)
  r2 = r*r
  
  st = SIN(xc(2))
  IF (st < 1.e-6) st = 1.
  st2 = st**2
  
  aom2 = aom**2
  
  delta = r2 - 2.0 * Mbh * r + aom2
  rho2  = r2 + aom2*(1.0 - st2)
  sigma = (r2 + aom2)**2 - aom2*delta*st2
  zz    = 2.0*r/rho2
  
  lapse    = 1.0 / sqrt(1.0 + zz)
  shift(1) = zz/(1.0 + zz)       
  shift(2) = 0.0         
  shift(3) = 0.0  
  
  g_cov( 1, 1:3) = (/ 1.0 + zz,             0.0,        -aom*st2*(1.0 + zz)   /)
  g_cov( 2, 1:3) = (/ 0.0,                  rho2,       0.0                   /)
  g_cov( 3, 1:3) = (/ -aom*st2*(1.0 + zz),  0.0,        (sigma/rho2)*st2      /)
  
  g_contr( 1, 1:3) = (/  lapse**2 + aom2*st2/rho2,   0.0,                       aom/rho2         /)
  g_contr( 2, 1:3) = (/  0.0,                        1.0/rho2,                  0.0              /)
  g_contr( 3, 1:3) = (/  aom/rho2,                   0.0,                       1.0/(rho2*st2)   /)
  
  gp = rho2*st*sqrt(1.0 + zz)
  gm = 1.0/gp  
  !
#endif  
  CASE('GRMHDAlfven')
      !
      lapse    = 1.0  
      shift(1) = 0.0
      shift(2) = 0.0         
      shift(3) = 0.0  
      !
      g_cov = 0.
      g_contr = 0. 
      g_cov( 1, 1) = 1.0
      g_cov( 2, 2) = 1.0
      g_cov( 3, 3) = 1.0
      g_contr( 1, 1) = 1.0
      g_contr( 2, 2) = 1.0
      g_contr( 3, 3) = 1.0 
      !
      gp = 1.0
      gm = 1.0
      !
   END SELECT 
   ! 
END SUBROUTINE METRIC


SUBROUTINE DMETRIC (xc, dalpha, BB, DD, dphi)

  IMPLICIT NONE 
  !
  REAL, DIMENSION(3),     intent(IN)  :: xc
  REAL, DIMENSION(3),     intent(OUT) :: dalpha
  REAL, DIMENSION(3,3),   intent(OUT) :: BB
  REAL, DIMENSION(3,3,3), intent(OUT) :: DD
  REAL, DIMENSION(3),     intent(OUT) :: dphi
  INTEGER                             :: i 
  REAL                                :: xp(3), xm(3), gp(3,3), gm(3,3), ap, am, betap(3), betam(3), tempp, tempm, temp(3,3)  
  REAL                                :: xp1(3), xm1(3), gp1(3,3), gm1(3,3), ap1, am1, betap1(3), betam1(3) 
  REAL                                :: xp2(3), xm2(3), gp2(3,3), gm2(3,3), ap2, am2, betap2(3), betam2(3) 
  REAL                                :: alpha, beta(3), g_cov(3,3), phip, phim, phip1, phip2, phim1, phim2, Kex(3,3)  
  REAL, PARAMETER                     :: epsilon = 1e-7, eps4 = 1e-4   
  !
  ! Metric derivative computed with a simple second order central finite difference 
  !
  !DO i = 1, 3
  !  xp = xc
  !  xp(i) = xp(i)+epsilon 
  !  xm = xc
  !  xm(i) = xm(i)-epsilon 
  !  CALL METRIC ( xp, ap, tempp, tempm, betap, Kex, gp, temp, phip)
  !  CALL METRIC ( xm, am, tempp, tempm, betam, Kex, gm, temp, phim)
  !  dalpha(i) = (ap-am)/(2*epsilon) 
  !  BB(i,:)   = (betap(:)-betam(:))/(2*epsilon) 
  !  DD(i,:,:) = (gp(:,:)-gm(:,:))/(2*epsilon) 
  !  dphi(i)   = (phip - phim)/(2*epsilon) 
  !  CONTINUE 
  !ENDDO
  !
  ! Metric derivative computed with a fourth order central finite difference 
  !
  DO i = 1, 3
    xp1 = xc
    xp1(i) = xp1(i)+eps4 
    xm1 = xc
    xm1(i) = xm1(i)-eps4 
    xp2 = xc
    xp2(i) = xp2(i)+2*eps4 
    xm2 = xc
    xm2(i) = xm2(i)-2*eps4 
    CALL METRIC ( xp1, ap1, tempp, tempm, betap1, Kex, gp1, temp, phip1)
    CALL METRIC ( xm1, am1, tempp, tempm, betam1, Kex, gm1, temp, phim1)
    CALL METRIC ( xp2, ap2, tempp, tempm, betap2, Kex, gp2, temp, phip2)
    CALL METRIC ( xm2, am2, tempp, tempm, betam2, Kex, gm2, temp, phim2)
    dalpha(i) = ( 8.0*ap1      -8.0*am1       +am2       -ap2       )/(12.0*eps4) 
    BB(i,:)   = ( 8.0*betap1(:)-8.0*betam1(:) +betam2(:) -betap2(:) )/(12.0*eps4) 
    DD(i,:,:) = ( 8.0*gp1(:,:) -8.0*gm1(:,:)  +gm2(:,:)  -gp2(:,:)  )/(12.0*eps4) 
    dphi(i)   = ( 8.0*phip1    -8.0*phim1     +phim2     -phip2     )/(12.0*eps4) 
  ENDDO
  !
 END SUBROUTINE DMETRIC  
       
 SUBROUTINE MatrixInverse3x3(M,iM,det) 
    !---------------
    ! compute the determinant det of the NxN-matrix M
    !---------------
    IMPLICIT NONE
    ! input variables 
    REAL, INTENT(IN)   :: M(3,3)
    ! output variables
    REAL, INTENT(OUT)    :: iM(3,3)
    REAL, INTENT(OUT)    :: det
    ! output variables
    REAL    :: ComputeDet,Id(3,3)
    INTEGER :: i,j
    ! 
    det = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-M(2,1)*M(1,2)*M(3,3)+M(2,1)*M(1,3)*M(3,2)+M(3,1)*M(1,2)*M(2,3)-M(3,1)*M(1,3)*M(2,2) !ComputeDet(M,3)
    IF(det*det.LT.1e-20) THEN
        print *, 'FATAL ERROR: det = 0'
        STOP
    ENDIF
    !
    iM(1,1) =M(2,2)*M(3,3)-M(2,3)*M(3,2)
    iM(1,2) =M(1,3)*M(3,2)-M(1,2)*M(3,3)
    iM(1,3) =M(1,2)*M(2,3)-M(1,3)*M(2,2)
    iM(2,1) =M(2,3)*M(3,1)-M(2,1)*M(3,3)
    iM(2,2) =M(1,1)*M(3,3)-M(1,3)*M(3,1)
    iM(2,3) =M(1,3)*M(2,1)-M(1,1)*M(2,3)
    iM(3,1) =M(2,1)*M(3,2)-M(2,2)*M(3,1)
    iM(3,2) =M(1,2)*M(3,1)-M(1,1)*M(3,2)
    iM(3,3) =M(1,1)*M(2,2)-M(1,2)*M(2,1)
    iM = iM/det
    !
    !Id = MATMUL(M,iM)
    !DO i=1,3
    !    DO j=1,3
    !        IF(i.eq.j) THEN
    !            IF((Id(i,j)-1.)**2..GT.1e-18) THEN
    !                print *, 'FATAL ERROR: iM*M !=1'
    !                STOP
    !            ENDIF
    !        ELSE
    !            IF((Id(i,j)**2).GT.1e-18) THEN
    !                print *, 'FATAL ERROR: iM*M !=1'
    !                STOP
    !            ENDIF
    !        ENDIF
    !    ENDDO
    !ENDDO
    !
    CONTINUE
    !
    END SUBROUTINE MatrixInverse3x3
    
    
    
    
#ifdef TORUS  
    
MODULE Metric_mod 
    IMPLICIT NONE
    PRIVATE
    !
    INTERFACE METRIC_KSS  
        MODULE PROCEDURE METRIC_KSS
    END INTERFACE
    ! 
    INTERFACE Cart2SphMatrix_cov  
        MODULE PROCEDURE Cart2SphMatrix_cov
    END INTERFACE
    ! 
    INTERFACE Cart2Sph_cov  
        MODULE PROCEDURE Cart2Sph_cov
    END INTERFACE
    ! 
    PUBLIC :: METRIC_KSS,Cart2SphMatrix_cov,Cart2Sph_cov
    
    CONTAINS


SUBROUTINE METRIC_KSS( xc, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
#ifdef GRMHD
  USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
  IMPLICIT NONE 
  !
  REAL, DIMENSION(3), intent(IN) :: xc
  REAL                           :: lapse,gp,gm, phi 
  REAL,dimension(3)              :: shift
  REAL,dimension(3,3)            :: Kex, g_cov, g_contr 

  REAL :: x, y, z, r, z2, r2, aom2, det 
  REAL :: lx, ly, lz, HH, SS, detg, rho, r_hz, sig, theta, a, M    
  REAL :: st, st2, delta, rho2, sigma, zz, V0Z4(54),V070(70) 
  REAL :: transl(3), xGP(3) 
#ifndef GRMHD  
  REAL, PARAMETER :: aom = 0.9, Mbh = 1.0, P_eps = 1e-4    
#endif  
  !
  Kex = 0.0 
  phi = 1.0
  !
  ! Rotating black hole in Kerr-Schild spherical coordinates. See Appendix B of Komissarov (2004) MNRAS, 350, 427
  r  = xc(1)
  r2 = r*r
  
  st = SIN(xc(2))
  IF (st < 1.e-6) st = 1.
  st2 = st**2
  
  aom2 = aom**2
  
  delta = r2 - 2.0 * Mbh * r + aom2
  rho2  = r2 + aom2*(1.0 - st2)
  sigma = (r2 + aom2)**2 - aom2*delta*st2
  zz    = 2.0*r/rho2
  
  lapse    = 1.0 / sqrt(1.0 + zz)
  shift(1) = zz/(1.0 + zz)       
  shift(2) = 0.0         
  shift(3) = 0.0  
  
  g_cov( 1, 1:3) = (/ 1.0 + zz,             0.0,        -aom*st2*(1.0 + zz)   /)
  g_cov( 2, 1:3) = (/ 0.0,                  rho2,       0.0                   /)
  g_cov( 3, 1:3) = (/ -aom*st2*(1.0 + zz),  0.0,        (sigma/rho2)*st2      /)
  
  g_contr( 1, 1:3) = (/  lapse**2 + aom2*st2/rho2,   0.0,                       aom/rho2         /)
  g_contr( 2, 1:3) = (/  0.0,                        1.0/rho2,                  0.0              /)
  g_contr( 3, 1:3) = (/  aom/rho2,                   0.0,                       1.0/(rho2*st2)   /)
  
  gp = rho2*st*sqrt(1.0 + zz)
  gm = 1.0/gp  
  !
END SUBROUTINE METRIC_KSS


    
    
    SUBROUTINE Cart2Sph_cov(B_Cart,B_sph,x) 
#ifdef GRMHD
  USE typesDef, ONLY: aom 
#endif  
    IMPLICIT NONE
    !input variables
    REAL :: B_Cart(3),B_sph(3),x(3) 
    !local var.
    REAL :: A(3,3),theta,sintheta,costheta,cosphip,sinphip,phi,r,r2,ir2,raom,atanAR,datanAR 
    !!
    !! x=r*cos(theta)*cos(phi)
    !! y=r*cos(theta)*sin(phi)
    !! z=r*sin(theta) 
    !!
    !A(1,1:d) = (/DCOS(r(2))*DCOS(r(3)), - r(1)*DSIN(r(2))*DCOS(r(3)), - r(1)*DCOS(r(2))*DSIN(r(3))          /)
    !A(2,1:d) = (/DCOS(r(2))*DSIN(r(3)), - r(1)*DSIN(r(2))*DSIN(r(3)),   r(1)*DCOS(r(2))*DCOS(r(3))          /)
    !A(3,1:d) = (/DSIN(r(2))           ,   r(1)*DCOS(r(2))           ,   0.0                                 /)
    !
    IF(aom.eq.0.) THEN
        !
        ! x=r*sin(theta)*cos(phi)
        ! y=r*sin(theta)*sin(phi)
        ! z=r*cos(theta) 
        r = x(1)
        theta = x(2)
        phi = x(3)
        !
        sintheta = DSIN(theta)
        costheta = DCOS(theta)
        cosphip = DCOS(phi)
        sinphip = DSIN(phi)
        ! 
        ! wrong. we need the transpose of.
        !A(1,1:d) = (/sintheta*cosphip, + r*costheta*cosphip, - r*sintheta*sinphip          /)
        !A(2,1:d) = (/sintheta*sinphip, + r*costheta*sinphip,   r*sintheta*cosphip          /)
        !A(3,1:d) = (/costheta        , - r*sintheta        ,   0.0                                 /)
        !
        ! 
        A(1:3,1) = (/sintheta*cosphip, + r*costheta*cosphip, - r*sintheta*sinphip          /)
        A(1:3,2) = (/sintheta*sinphip, + r*costheta*sinphip,   r*sintheta*cosphip          /)
        A(1:3,3) = (/costheta        , - r*sintheta        ,   0.0                                 /)
        !
    ELSE
        ! ! Kerr-Schild Cartesian
        ! x=SQRT(r**2+aom**2)*sin(theta)*cos(phi-ATAN2(aom,r))
        ! y=SQRT(r**2+aom**2)*sin(theta)*SIN(phi-ATAN2(aom,r))
        ! z=r*COS(theta)
        !
        ! ! from Kerr-Schild spherical
        r = x(1)
        theta = x(2)
        phi = x(3)
        !
        raom = SQRT(r**2+aom**2)
        atanAR = ATAN2(aom,r)
        datanAR = 1.0/(1.+atanAR**2)
        r2 = r**2
        ir2 = 1./r2
        sintheta = DSIN(theta)
        costheta = DCOS(theta)
        cosphip = DCOS(phi-atanAR)
        sinphip = DSIN(phi-atanAR)
    
        ! wrong. we need the transpose of.
        !A(1,1:d) = (/    0.5/raom*2.0*r*sintheta*cosphip - raom*sintheta*sinphip*(-datanAR)*(-aom*ir2) , +raom*costheta*cosphip , -raom*sintheta*sinphip          /)
        !A(2,1:d) = (/    0.5/raom*2.0*r*sintheta*sinphip + raom*sintheta*cosphip*(-datanAR)*(-aom*ir2) , +raom*costheta*sinphip ,  raom*sintheta*cosphip           /)
        !A(3,1:d) = (/    costheta          ,  - r*sintheta         ,   0.0                                 /)
        !
        A(1:3,1) = (/    0.5/raom*2.0*r*sintheta*cosphip - raom*sintheta*sinphip*(-datanAR)*(-aom*ir2) , +raom*costheta*cosphip , -raom*sintheta*sinphip          /)
        A(1:3,2) = (/    0.5/raom*2.0*r*sintheta*sinphip + raom*sintheta*cosphip*(-datanAR)*(-aom*ir2) , +raom*costheta*sinphip ,  raom*sintheta*cosphip           /)
        A(1:3,3) = (/    costheta          ,  - r*sintheta         ,   0.0                                 /)
        !                                                       ,   raom*sintheta*cosphip                                                          ,   0.0                                 /)
    
    
    
    ENDIF
    !
    B_sph = MATMUL(A,B_Cart)
    !
    continue
    !
    END SUBROUTINE Cart2Sph_cov
    
    
    SUBROUTINE Cart2SphMatrix_cov(A,x) 
#ifdef GRMHD
  USE typesDef, ONLY: aom 
#endif  
    IMPLICIT NONE
    !input variables
    REAL :: A(3,3),x(3) 
    !local var.
    REAL :: theta,sintheta,costheta,cosphip,sinphip,phi,r,r2,ir2,raom,atanAR,datanAR 
    !! x=r*cos(theta)*cos(phi)
    !! x=r*cos(theta)*sin(phi)
    !! x=r*sin(theta) 
    !!
    !A(1,1:d) = (/DCOS(r(2))*DCOS(r(3)), - r(1)*DSIN(r(2))*DCOS(r(3)), - r(1)*DCOS(r(2))*DSIN(r(3))          /)
    !A(2,1:d) = (/DCOS(r(2))*DSIN(r(3)), - r(1)*DSIN(r(2))*DSIN(r(3)),   r(1)*DCOS(r(2))*DCOS(r(3))          /)
    !A(3,1:d) = (/DSIN(r(2))           ,   r(1)*DCOS(r(2))           ,   0.0                                 /)
    !
    !
    ! x=r*sin(theta)*cos(phi)
    ! x=r*sin(theta)*sin(phi)
    ! x=r*cos(theta) 
    !
    IF(aom.eq.0.) THEN
        !
        ! x=r*sin(theta)*cos(phi)
        ! y=r*sin(theta)*sin(phi)
        ! z=r*cos(theta) 
        r = x(1)
        theta = x(2)
        phi = x(3)
        !
        sintheta = DSIN(theta)
        costheta = DCOS(theta)
        cosphip = DCOS(phi)
        sinphip = DSIN(phi)
        !  
        ! wrong!! we need the transpose of. for COVARIANT! (this was correct for contravariant)
        !A(1,1:d) = (/sintheta*coship, + r*costheta*cosphip, - r*sintheta*sinphip          /)
        !A(2,1:d) = (/sintheta*sinphip, + r*costheta*sinphip,   r*sintheta*cosphip          /)
        !A(3,1:d) = (/costheta        , - r*sintheta        ,   0.0                                 /)
        !
        A(1:3,1) = (/sintheta*cosphip, + r*costheta*cosphip, - r*sintheta*sinphip          /)
        A(1:3,2) = (/sintheta*sinphip, + r*costheta*sinphip,   r*sintheta*cosphip          /)
        A(1:3,3) = (/costheta        , - r*sintheta        ,   0.0                                 /)
        !
    ELSE
        ! ! Kerr-Schild Cartesian
        ! x=SQRT(r**2+aom**2)*sin(theta)*cos(phi-ATAN2(aom,r))
        ! y=SQRT(r**2+aom**2)*sin(theta)*SIN(phi-ATAN2(aom,r))
        ! z=r*COS(theta)
        !
        ! ! from Kerr-Schild spherical
        r = x(1)
        theta = x(2)
        phi = x(3)
        !
        raom = SQRT(r**2+aom**2)
        atanAR = ATAN2(aom,r)
        datanAR = 1.0/(1.+atanAR**2)
        r2 = r**2
        ir2 = 1./r2
        sintheta = DSIN(theta)
        costheta = DCOS(theta)
        cosphip = DCOS(phi-atanAR)
        sinphip = DSIN(phi-atanAR)
    
        ! wrong. we need the transpose of.
        !A(1,1:d) = (/    0.5/raom*2.0*r*sintheta*cosphip - raom*sintheta*sinphip*(-datanAR)*(-aom*ir2) , +raom*costheta*cosphip , -raom*sintheta*sinphip          /)
        !A(2,1:d) = (/    0.5/raom*2.0*r*sintheta*sinphip + raom*sintheta*cosphip*(-datanAR)*(-aom*ir2) , +raom*costheta*sinphip ,  raom*sintheta*cosphip           /)
        !A(3,1:d) = (/    costheta          ,  - r*sintheta         ,   0.0                                 /)
        !
        A(1:3,1) = (/    0.5/raom*2.0*r*sintheta*cosphip - raom*sintheta*sinphip*(-datanAR)*(-aom*ir2) , +raom*costheta*cosphip , -raom*sintheta*sinphip          /)
        A(1:3,2) = (/    0.5/raom*2.0*r*sintheta*sinphip + raom*sintheta*cosphip*(-datanAR)*(-aom*ir2) , +raom*costheta*sinphip ,  raom*sintheta*cosphip           /)
        A(1:3,3) = (/    costheta          ,  - r*sintheta         ,   0.0                                 /)
        !
    
    ENDIF
    continue
    !
    END SUBROUTINE Cart2SphMatrix_cov

END MODULE Metric_mod

#endif