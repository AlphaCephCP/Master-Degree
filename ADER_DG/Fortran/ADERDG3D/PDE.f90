!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDEFlux(F,Q,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam)  
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
    REAL :: V(nVar) 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, phi, fa, k0, beta0(3)  
    REAL :: gamma1, rho, vx, vy, vz, bx, by, bz, ex, ey, ez, v2, b2, e2, lf, w, ww, uem, wwx, wwy, wwz 
    REAL :: gp, gm, g_cov(3,3), g_contr(3,3), vxB(3), vxB_contr(3), BV(3), BQ(3), vtr(3), vf(3), lapse, shift(3)    
    REAL :: Fij(3,3), vf_cov(3), Qv_contr(3), QB_contr(3), Bv_contr(3), psi 
    REAL :: QGRMHD(19), FGRMHD(19,d) 
    REAL :: A(3,3), GT(3,3), devG(3,3), AU(3), detA, Id(3,3), T, TT(3,3), falpha    
    !
    F = 0.0 
    !
#ifdef EULER     
    !
    ! 3D compressible Euler equations 
    !
    irho = 1.0/Q(1)
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
    ! 
    F(1,1) = Q(2) 
    F(2,1) = irho*Q(2)*Q(2) + p 
    F(3,1) = irho*Q(2)*Q(3)
    F(4,1) = irho*Q(2)*Q(4)
    F(5,1) = irho*Q(2)*(Q(5)+p)  
    !
    F(1,2) = Q(3) 
    F(2,2) = irho*Q(3)*Q(2)  
    F(3,2) = irho*Q(3)*Q(3) + p 
    F(4,2) = irho*Q(3)*Q(4)
    F(5,2) = irho*Q(3)*(Q(5)+p)  
    ! 
    F(1,3) = Q(4) 
    F(2,3) = irho*Q(4)*Q(2)  
    F(3,3) = irho*Q(4)*Q(3)  
    F(4,3) = irho*Q(4)*Q(4) + p
    F(5,3) = irho*Q(4)*(Q(5)+p)  
    !
#endif
    !
#ifdef ____ELASTICITY

    lam  = par(1)   
    mu   = par(2) 
    irho = 1./par(3) 
    !
    F(1,1) = - (lam+2*mu)*Q(7) 
    F(2,1) = - lam*Q(7) 
    F(3,1) = - lam*Q(7) 
    F(4,1) = - mu *Q(8) 
    F(5,1) = 0. 
    F(6,1) = - mu *Q(9) 
    F(7,1) = - irho *Q(1) 
    F(8,1) = - irho *Q(4) 
    F(9,1) = - irho *Q(6) 
    
    F(1,2) = - lam*Q(8) 
    F(2,2) = - (lam+2*mu)*Q(8) 
    F(3,2) = - lam*Q(8) 
    F(4,2) = - mu *Q(7) 
    F(5,2) = - mu *Q(9) 
    F(6,2) = 0. 
    F(7,2) = - irho *Q(4) 
    F(8,2) = - irho *Q(2) 
    F(9,2) = - irho *Q(5) 
    
    F(1,3) = - lam*Q(9) 
    F(2,3) = - lam*Q(9) 
    F(3,3) = - (lam+2*mu)*Q(9) 
    F(4,3) = 0. 
    F(5,3) = - mu *Q(8) 
    F(6,3) = - mu *Q(7) 
    F(7,3) = - irho *Q(4) 
    F(8,3) = - irho *Q(2) 
    F(9,3) = - irho *Q(5) 
       
    !F = 0.0
    !F(:,1) = Q(:)    
        
#endif 
    !          
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    !
   k1  = EQN%CCZ4k1  
   k2  = EQN%CCZ4k2  
   fff = EQN%CCZ4f 
   ggg = EQN%CCZ4g 
   e   = EQN%CCZ4e 
   c   = EQN%CCZ4c  
   ds  = EQN%CCZ4ds 
   mu  = EQN%CCZ4mu 
   xi  = EQN%CCZ4xi
   sk  = EQN%CCZ4sk

      alpha = Q(17)
      SELECT CASE(EQN%CCZ4LapseType) 
      CASE(0)  ! harmonic 
          fa = 1.0 
      CASE DEFAULT  ! 1 + log 
          fa = 2.0/alpha
      END SELECT 


    K0 = Q(54) 
    beta0 = 0.0 ! Q(55:57)    

     F = 0.0 
         
      F(13,1) = mu*(-Q(31)-Q(35))
      F(13,2) = mu*Q(30)
      F(13,3) = mu*Q(33)
      F(14,1) = mu*Q(28)
      F(14,2) = mu*(-Q(27)-Q(35))
      F(14,3) = mu*Q(34)
      F(15,1) = mu*Q(29)
      F(15,2) = mu*Q(32)
      F(15,3) = mu*(-Q(27)-Q(31))

      !FTens(24,1) = -Q(17)*fa*K0
      !FTens(25,2) = -Q(17)*fa*K0
      !FTens(26,3) = -Q(17)*fa*K0

      F(27,1) = -fff*Q(21)-beta0(1)
      F(28,2) = -fff*Q(21)-beta0(1)
      F(29,3) = -fff*Q(21)-beta0(1)
      F(30,1) = -fff*Q(22)-beta0(2)
      F(31,2) = -fff*Q(22)-beta0(2)
      F(32,3) = -fff*Q(22)-beta0(2)
      F(33,1) = -fff*Q(23)-beta0(3)
      F(34,2) = -fff*Q(23)-beta0(3)
      F(35,3) = -fff*Q(23)-beta0(3)
    !
#ifdef Z4GRMHD
    QGRMHD(1:9)   = Q(55:63)    ! MHD variables 
    QGRMHD(10)    = Q(17)       ! lapse 
    QGRMHD(11:13) = Q(18:20)    ! shift 
    QGRMHD(14:19) = Q(1:6)      ! metric 
    CALL PDEFluxGRMHD(FGRMHD,QGRMHD,par) 
    F(55:63,:) = FGRMHD(1:9,:) 
    !
    CONTINUE
    !
#endif 
    !
#endif 
    !
#ifdef GPR3D 
    !
    !  1   2    3    4    5    6   7   8   9   10  11  12  13  14    15    16     17  
    ! rho rhou rhov rhow rhoE A11 A12 A13 A21 A22 A23 A31 A32 A33  rho j1 rho j2 rho j3 
    ! rho u    v    w    p    A11 A12 A13 A21 A22 A23 A31 A32 A33      j1     j2     j3 
    CALL PDECons2Prim(V,Q,iErr)  
    ! Compute the pressure (hydrodynamic part) 
    p = V(5) 
    ! Compute the viscous stress tensor 
    A(1,:) = (/ V( 6), V( 7), V( 8) /) 
    A(2,:) = (/ V( 9), V(10), V(11) /)
    A(3,:) = (/ V(12), V(13), V(14) /)   
    AU = MATMUL( A, V(2:4) ) 
    detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
    ! 
    GT     = MATMUL( TRANSPOSE(A), A ) 
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id
    TT      = -EQN%rho0*detA*EQN%cs**2*MATMUL(GT,devG) 
    ! Compute the temperature from the ideal gas law 
    T = p/V(1)/EQN%cv/(EQN%gamma-1)  
    falpha  = EQN%alpha**2  
    ! 
    F = 0.0
    ! 
    F(1,1) = Q(2) 
    F(2,1) = Q(2)*V(2) + p - TT(1,1) 
    F(3,1) = Q(3)*V(2)     - TT(2,1) 
    F(4,1) = Q(4)*V(2)     - TT(3,1) 
    F(5,1) = V(2)*(Q(5)+p) - V(2)*TT(1,1) - V(3)*TT(2,1) - V(4)*TT(3,1) + falpha*V(15)*T   
    F(6,1)  = AU(1) 
    F(9,1)  = AU(2) 
    F(12,1) = AU(3) 
    F(15,1) = Q(15)*V(2) + T
    F(16,1) = Q(16)*V(2) 
    F(17,1) = Q(17)*V(2) 
    !
    F(1,2) = Q(3) 
    F(2,2) = Q(2)*V(3)      - TT(1,2) 
    F(3,2) = Q(3)*V(3) + p  - TT(2,2) 
    F(4,2) = Q(4)*V(3)      - TT(3,2) 
    F(5,2) = V(3)*(Q(5)+p) - V(2)*TT(1,2) - V(3)*TT(2,2) - V(4)*TT(3,2) + falpha*V(16)*T 
    F(7,2)  = AU(1) 
    F(10,2) = AU(2) 
    F(13,2) = AU(3) 
    F(15,2) = Q(15)*V(3) 
    F(16,2) = Q(16)*V(3) + T 
    F(17,2) = Q(17)*V(3)   
    !
    F(1,3) = Q(4) 
    F(2,3) = Q(2)*V(4)      - TT(1,3) 
    F(3,3) = Q(3)*V(4)      - TT(2,3) 
    F(4,3) = Q(4)*V(4) + p  - TT(3,3) 
    F(5,3) = V(4)*(Q(5)+p) - V(2)*TT(1,3) - V(3)*TT(2,3) - V(4)*TT(3,3) + falpha*V(17)*T   
    F(8,3)  = AU(1) 
    F(11,3) = AU(2) 
    F(14,3) = AU(3) 
    F(15,3) = Q(15)*V(4) 
    F(16,3) = Q(16)*V(4) 
    F(17,3) = Q(17)*V(4) + T 
    !
#endif 
    !
    ! ------------------------
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
    !
    ! We do not use any conservative fluxes any more in this PDE 
    !
#ifdef CCZ4EINSTEIN
    RETURN
#endif 
    !
#ifdef CCZ4GRMHD 
    alpha         = EXP(Q(17)) 
    phi           = EXP(Q(55)) 
    QGRMHD(1:9)   = Q(60:68)        ! MHD variables 
    QGRMHD(10)    = alpha           ! lapse 
    QGRMHD(11:13) = Q(18:20)        ! shift 
    QGRMHD(14:19) = Q(1:6)/phi**2   ! metric 
    CALL PDEFluxGRMHD(FGRMHD,QGRMHD,par) 
    F(60:68,:) = FGRMHD(1:9,:) 
    !
    CONTINUE
    !    
#endif 
    !
#ifdef BSSZ4GRMHD 
    alpha         = EXP(Q(17)) 
    phi           = EXP(Q(55)) 
    QGRMHD(1:9)   = Q(63:71)        ! MHD variables 
    QGRMHD(10)    = alpha           ! lapse 
    QGRMHD(11:13) = Q(18:20)        ! shift 
    QGRMHD(14:19) = Q(1:6)/phi**2   ! metric 
    CALL PDEFluxGRMHD(FGRMHD,QGRMHD,par) 
    F(63:71,:) = FGRMHD(1:9,:) 
    !
    CONTINUE
    !
#endif 

    !    
#endif 
    ! 
    ! ------------------------
    !
#ifdef SRMHD
    CALL PDECons2Prim(V,Q,iErr)
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    rho    = V(1)
    vx     = V(2)
    vy     = V(3)
    vz     = V(4)
    p      = V(5)
    bx     = V(6)
    by     = V(7)
    bz     = V(8)
    !
    ex     = - (vy*bz - vz*by)
    ey     = - (vz*bx - vx*bz)
    ez     = - (vx*by - vy*bx)
    !
    v2     = vx**2 + vy**2 + vz**2
    b2     = bx**2 + by**2 + bz**2
    e2     = ex**2 + ey**2 + ez**2
    lf     = 1.0/sqrt(1.0 - v2)
    w      = rho + gamma1*p
    ww     = w*lf**2
    uem    = 0.5*(b2 + e2)
    wwx    = ww*vx
    wwy    = ww*vy
    wwz    = ww*vz
    !
    F(1,1)   = vx*rho*lf
    F(2,1)   = wwx*vx - bx*bx - ex*ex + p + uem
    F(3,1)   = wwx*vy - bx*by - ex*ey
    F(4,1)   = wwx*vz - bx*bz - ex*ez 
    F(5,1)   = wwx + (ey*bz - ez*by) 
    F(6,1)   = V(9)
    F(7,1)   = -ez
    F(8,1)   = ey  
    F(9,1)   = EQN%DivCleaning_a**2*bx
    !
    F(1,2)   = vy*rho*lf
    F(2,2)   = wwy*vx - by*bx - ey*ex 
    F(3,2)   = wwy*vy - by*by - ey*ey + p + uem
    F(4,2)   = wwy*vz - by*bz - ey*ez 
    F(5,2)   = wwy + (ez*bx - ex*bz) 
    F(6,2)   = ez 
    F(7,2)   = V(9) 
    F(8,2)   = -ex   
    F(9,2)   = EQN%DivCleaning_a**2*by
    !
    F(1,3)   = vz*rho*lf
    F(2,3)   = wwz*vx - bz*bx - ez*ex 
    F(3,3)   = wwz*vy - bz*by - ez*ey 
    F(4,3)   = wwz*vz - bz*bz - ez*ez + p + uem
    F(5,3)   = wwz + (ex*by - ey*bx) 
    F(6,3)   = -ey  
    F(7,3)   = ex   
    F(8,3)   = V(9)   
    F(9,3)   = EQN%DivCleaning_a**2*bz    
#endif 
    !
#ifdef GRMHD
  !
  CALL PDEFluxGRMHD(F,Q,par) 
  RETURN 
  !  
#endif 
  !
END SUBROUTINE PDEFlux
!
!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDESource(S,Q,par,time) 
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), time   
    REAL, INTENT(OUT) :: S(nVar) 
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,iErr  
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3), eta 
    REAL :: g_cov(3,3), det, g_contr(3,3), Christoffel(3,3,3), dgup(3,3,3), faa, b0(3), dChristoffelSrc(3,3,3,3), RiemannSrc(3,3,3,3)
    REAL :: RicciSrc(3,3), gammaup(3), dtmetric(3,3), dtgup(3,3), dtChristoffelSrc(3,3,3), dtGammaUpSrc(3)
    REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25
    REAL :: QGRMHD(19), gradQGRMHD(19,d), BgradQGRMHD(19), V(nVar), T   
    REAL :: dtgamma(3,3), dtK(3,3), dtTheta, dtZ(3), dtalpha, dtGhat(3), dtbeta(3), dtbb(3), dtA(3), dtB(3,3), dtD(3,3,3) 
    REAL :: Kex(3,3), AA(3), BB(3,3), DD(3,3,3), Kmix(3,3), Kup(3,3), Z(3), Zup(3), theta, beta(3), nablaijalphasrc(3,3) 
    REAL :: traceK, Kupdown, nablaZSrc(3,3), RicciPlusNablaZSrc(3,3), RPlusTwoNablaZSrc, dKtempsrc(3), dtraceKSrc(3), QG(3), b(3)  
    REAL :: AM(3,3), G(3,3), devG(3,3), Id(3,3), deta, deta2, psim(3,3), temp2   
    !
    S = 0.0 
    !
#ifdef EULER     
    !
    ! 3D compressible Euler equations 
    !
    !S(1) = -1.0 ! -1.0e-1*Q(1)   ! density decay 
    !S(2)   =  -1.0*Q(2)   ! friction  
    !S(5)   =  DOT_PRODUCT(Q(2:4)/Q(1),S(2:4)) 
    !
#endif
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    !
    k1  = EQN%CCZ4k1  
    k2  = EQN%CCZ4k2  
    fff = EQN%CCZ4f 
    ggg = EQN%CCZ4g 
    e   = EQN%CCZ4e 
    c   = EQN%CCZ4c  
    ds  = EQN%CCZ4ds 
    mu  = EQN%CCZ4mu 
    eta = EQN%CCZ4eta 
    xi  = EQN%CCZ4xi 
    sk  = EQN%CCZ4sk 
    !
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
    g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
    g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
      
    alpha = Q(17)
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0  
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 

    K0    = Q(54) 
    beta0 = 0.0 ! Q(55:57) 
    b0    = 0.0 ! Q(58:60)

    Christoffel(1,1,1) = g_contr(1,1)*Q(36)+g_contr(1,2)*(2*Q(37)-Q(42))+g_contr(1,3)*(2*Q(38)-Q(48))
    Christoffel(1,1,2) = g_contr(2,1)*Q(36)+g_contr(2,2)*(2*Q(37)-Q(42))+g_contr(2,3)*(2*Q(38)-Q(48))
    Christoffel(1,1,3) = g_contr(3,1)*Q(36)+g_contr(3,2)*(2*Q(37)-Q(42))+g_contr(3,3)*(2*Q(38)-Q(48))
    Christoffel(1,2,1) = g_contr(1,1)*Q(42)+g_contr(1,2)*Q(39)+g_contr(1,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(1,2,2) = g_contr(2,1)*Q(42)+g_contr(2,2)*Q(39)+g_contr(2,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(1,2,3) = g_contr(3,1)*Q(42)+g_contr(3,2)*Q(39)+g_contr(3,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(1,3,1) = g_contr(1,1)*Q(48)+g_contr(1,2)*(Q(49)+Q(40)-Q(44))+g_contr(1,3)*Q(41)
    Christoffel(1,3,2) = g_contr(2,1)*Q(48)+g_contr(2,2)*(Q(49)+Q(40)-Q(44))+g_contr(2,3)*Q(41)
    Christoffel(1,3,3) = g_contr(3,1)*Q(48)+g_contr(3,2)*(Q(49)+Q(40)-Q(44))+g_contr(3,3)*Q(41)
    Christoffel(2,1,1) = g_contr(1,1)*Q(42)+g_contr(1,2)*Q(39)+g_contr(1,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(2,1,2) = g_contr(2,1)*Q(42)+g_contr(2,2)*Q(39)+g_contr(2,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(2,1,3) = g_contr(3,1)*Q(42)+g_contr(3,2)*Q(39)+g_contr(3,3)*(Q(44)+Q(40)-Q(49))
    Christoffel(2,2,1) = g_contr(1,1)*(2*Q(43)-Q(39))+g_contr(1,2)*Q(45)+g_contr(1,3)*(2*Q(46)-Q(51))
    Christoffel(2,2,2) = g_contr(2,1)*(2*Q(43)-Q(39))+g_contr(2,2)*Q(45)+g_contr(2,3)*(2*Q(46)-Q(51))
    Christoffel(2,2,3) = g_contr(3,1)*(2*Q(43)-Q(39))+g_contr(3,2)*Q(45)+g_contr(3,3)*(2*Q(46)-Q(51))
    Christoffel(2,3,1) = g_contr(1,1)*(Q(49)+Q(44)-Q(40))+g_contr(1,2)*Q(51)+g_contr(1,3)*Q(47)
    Christoffel(2,3,2) = g_contr(2,1)*(Q(49)+Q(44)-Q(40))+g_contr(2,2)*Q(51)+g_contr(2,3)*Q(47)
    Christoffel(2,3,3) = g_contr(3,1)*(Q(49)+Q(44)-Q(40))+g_contr(3,2)*Q(51)+g_contr(3,3)*Q(47)
    Christoffel(3,1,1) = g_contr(1,1)*Q(48)+g_contr(1,2)*(Q(49)+Q(40)-Q(44))+g_contr(1,3)*Q(41)
    Christoffel(3,1,2) = g_contr(2,1)*Q(48)+g_contr(2,2)*(Q(49)+Q(40)-Q(44))+g_contr(2,3)*Q(41)
    Christoffel(3,1,3) = g_contr(3,1)*Q(48)+g_contr(3,2)*(Q(49)+Q(40)-Q(44))+g_contr(3,3)*Q(41)
    Christoffel(3,2,1) = g_contr(1,1)*(Q(49)+Q(44)-Q(40))+g_contr(1,2)*Q(51)+g_contr(1,3)*Q(47)
    Christoffel(3,2,2) = g_contr(2,1)*(Q(49)+Q(44)-Q(40))+g_contr(2,2)*Q(51)+g_contr(2,3)*Q(47)
    Christoffel(3,2,3) = g_contr(3,1)*(Q(49)+Q(44)-Q(40))+g_contr(3,2)*Q(51)+g_contr(3,3)*Q(47)
    Christoffel(3,3,1) = g_contr(1,1)*(2*Q(50)-Q(41))+g_contr(1,2)*(2*Q(52)-Q(47))+g_contr(1,3)*Q(53)
    Christoffel(3,3,2) = g_contr(2,1)*(2*Q(50)-Q(41))+g_contr(2,2)*(2*Q(52)-Q(47))+g_contr(2,3)*Q(53)
    Christoffel(3,3,3) = g_contr(3,1)*(2*Q(50)-Q(41))+g_contr(3,2)*(2*Q(52)-Q(47))+g_contr(3,3)*Q(53)
      

    dgup(1,1,1) = -2*g_contr(1,1)**2*Q(36)-2*g_contr(1,1)*g_contr(2,1)*Q(37)-2*g_contr(1,1)*g_contr(3,1)*Q(38)-2*g_contr(1,2)*g_contr(1,1)*Q(37)-2*g_contr(1,2)*g_contr(2,1)*Q(39)-2*g_contr(1,2)*g_contr(3,1)*Q(40)-2*g_contr(1,3)*g_contr(1,1)*Q(38)-2*g_contr(1,3)*g_contr(2,1)*Q(40)-2*g_contr(1,3)*g_contr(3,1)*Q(41)
    dgup(1,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(36)-2*g_contr(1,1)*g_contr(2,2)*Q(37)-2*g_contr(1,1)*g_contr(3,2)*Q(38)-2*g_contr(1,2)**2*Q(37)-2*g_contr(1,2)*g_contr(2,2)*Q(39)-2*g_contr(1,2)*g_contr(3,2)*Q(40)-2*g_contr(1,3)*g_contr(1,2)*Q(38)-2*g_contr(1,3)*g_contr(2,2)*Q(40)-2*g_contr(1,3)*g_contr(3,2)*Q(41)
    dgup(1,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(36)-2*g_contr(1,1)*g_contr(2,3)*Q(37)-2*g_contr(1,1)*g_contr(3,3)*Q(38)-2*g_contr(1,2)*g_contr(1,3)*Q(37)-2*g_contr(1,2)*g_contr(2,3)*Q(39)-2*g_contr(1,2)*g_contr(3,3)*Q(40)-2*g_contr(1,3)**2*Q(38)-2*g_contr(1,3)*g_contr(2,3)*Q(40)-2*g_contr(1,3)*g_contr(3,3)*Q(41)
    dgup(1,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(36)-2*g_contr(2,1)**2*Q(37)-2*g_contr(2,1)*g_contr(3,1)*Q(38)-2*g_contr(1,1)*g_contr(2,2)*Q(37)-2*g_contr(2,2)*g_contr(2,1)*Q(39)-2*g_contr(2,2)*g_contr(3,1)*Q(40)-2*g_contr(2,3)*g_contr(1,1)*Q(38)-2*g_contr(2,3)*g_contr(2,1)*Q(40)-2*g_contr(2,3)*g_contr(3,1)*Q(41)
    dgup(1,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(36)-2*g_contr(2,1)*g_contr(2,2)*Q(37)-2*g_contr(2,1)*g_contr(3,2)*Q(38)-2*g_contr(2,2)*g_contr(1,2)*Q(37)-2*g_contr(2,2)**2*Q(39)-2*g_contr(2,2)*g_contr(3,2)*Q(40)-2*g_contr(2,3)*g_contr(1,2)*Q(38)-2*g_contr(2,3)*g_contr(2,2)*Q(40)-2*g_contr(2,3)*g_contr(3,2)*Q(41)
    dgup(1,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(36)-2*g_contr(2,1)*g_contr(2,3)*Q(37)-2*g_contr(2,1)*g_contr(3,3)*Q(38)-2*g_contr(2,2)*g_contr(1,3)*Q(37)-2*g_contr(2,2)*g_contr(2,3)*Q(39)-2*g_contr(2,2)*g_contr(3,3)*Q(40)-2*g_contr(2,3)*g_contr(1,3)*Q(38)-2*g_contr(2,3)**2*Q(40)-2*g_contr(2,3)*g_contr(3,3)*Q(41)
    dgup(1,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(36)-2*g_contr(3,1)*g_contr(2,1)*Q(37)-2*g_contr(3,1)**2*Q(38)-2*g_contr(3,2)*g_contr(1,1)*Q(37)-2*g_contr(3,2)*g_contr(2,1)*Q(39)-2*g_contr(3,2)*g_contr(3,1)*Q(40)-2*g_contr(1,1)*g_contr(3,3)*Q(38)-2*g_contr(3,3)*g_contr(2,1)*Q(40)-2*g_contr(3,3)*g_contr(3,1)*Q(41)
    dgup(1,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(36)-2*g_contr(3,1)*g_contr(2,2)*Q(37)-2*g_contr(3,1)*g_contr(3,2)*Q(38)-2*g_contr(3,2)*g_contr(1,2)*Q(37)-2*g_contr(3,2)*g_contr(2,2)*Q(39)-2*g_contr(3,2)**2*Q(40)-2*g_contr(3,3)*g_contr(1,2)*Q(38)-2*g_contr(2,2)*g_contr(3,3)*Q(40)-2*g_contr(3,3)*g_contr(3,2)*Q(41)
    dgup(1,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(36)-2*g_contr(3,1)*g_contr(2,3)*Q(37)-2*g_contr(3,1)*g_contr(3,3)*Q(38)-2*g_contr(3,2)*g_contr(1,3)*Q(37)-2*g_contr(3,2)*g_contr(2,3)*Q(39)-2*g_contr(3,2)*g_contr(3,3)*Q(40)-2*g_contr(3,3)*g_contr(1,3)*Q(38)-2*g_contr(3,3)*g_contr(2,3)*Q(40)-2*g_contr(3,3)**2*Q(41)
    dgup(2,1,1) = -2*g_contr(1,1)**2*Q(42)-2*g_contr(1,1)*g_contr(2,1)*Q(43)-2*g_contr(1,1)*g_contr(3,1)*Q(44)-2*g_contr(1,2)*g_contr(1,1)*Q(43)-2*g_contr(1,2)*g_contr(2,1)*Q(45)-2*g_contr(1,2)*g_contr(3,1)*Q(46)-2*g_contr(1,3)*g_contr(1,1)*Q(44)-2*g_contr(1,3)*g_contr(2,1)*Q(46)-2*g_contr(1,3)*g_contr(3,1)*Q(47)
    dgup(2,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(42)-2*g_contr(1,1)*g_contr(2,2)*Q(43)-2*g_contr(1,1)*g_contr(3,2)*Q(44)-2*g_contr(1,2)**2*Q(43)-2*g_contr(1,2)*g_contr(2,2)*Q(45)-2*g_contr(1,2)*g_contr(3,2)*Q(46)-2*g_contr(1,3)*g_contr(1,2)*Q(44)-2*g_contr(1,3)*g_contr(2,2)*Q(46)-2*g_contr(1,3)*g_contr(3,2)*Q(47)
    dgup(2,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(42)-2*g_contr(1,1)*g_contr(2,3)*Q(43)-2*g_contr(1,1)*g_contr(3,3)*Q(44)-2*g_contr(1,2)*g_contr(1,3)*Q(43)-2*g_contr(1,2)*g_contr(2,3)*Q(45)-2*g_contr(1,2)*g_contr(3,3)*Q(46)-2*g_contr(1,3)**2*Q(44)-2*g_contr(1,3)*g_contr(2,3)*Q(46)-2*g_contr(1,3)*g_contr(3,3)*Q(47)
    dgup(2,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(42)-2*g_contr(2,1)**2*Q(43)-2*g_contr(2,1)*g_contr(3,1)*Q(44)-2*g_contr(1,1)*g_contr(2,2)*Q(43)-2*g_contr(2,2)*g_contr(2,1)*Q(45)-2*g_contr(2,2)*g_contr(3,1)*Q(46)-2*g_contr(2,3)*g_contr(1,1)*Q(44)-2*g_contr(2,3)*g_contr(2,1)*Q(46)-2*g_contr(2,3)*g_contr(3,1)*Q(47)
    dgup(2,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(42)-2*g_contr(2,1)*g_contr(2,2)*Q(43)-2*g_contr(2,1)*g_contr(3,2)*Q(44)-2*g_contr(2,2)*g_contr(1,2)*Q(43)-2*g_contr(2,2)**2*Q(45)-2*g_contr(2,2)*g_contr(3,2)*Q(46)-2*g_contr(2,3)*g_contr(1,2)*Q(44)-2*g_contr(2,3)*g_contr(2,2)*Q(46)-2*g_contr(2,3)*g_contr(3,2)*Q(47)
    dgup(2,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(42)-2*g_contr(2,1)*g_contr(2,3)*Q(43)-2*g_contr(2,1)*g_contr(3,3)*Q(44)-2*g_contr(2,2)*g_contr(1,3)*Q(43)-2*g_contr(2,2)*g_contr(2,3)*Q(45)-2*g_contr(2,2)*g_contr(3,3)*Q(46)-2*g_contr(2,3)*g_contr(1,3)*Q(44)-2*g_contr(2,3)**2*Q(46)-2*g_contr(2,3)*g_contr(3,3)*Q(47)
    dgup(2,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(42)-2*g_contr(3,1)*g_contr(2,1)*Q(43)-2*g_contr(3,1)**2*Q(44)-2*g_contr(3,2)*g_contr(1,1)*Q(43)-2*g_contr(3,2)*g_contr(2,1)*Q(45)-2*g_contr(3,2)*g_contr(3,1)*Q(46)-2*g_contr(1,1)*g_contr(3,3)*Q(44)-2*g_contr(3,3)*g_contr(2,1)*Q(46)-2*g_contr(3,3)*g_contr(3,1)*Q(47)
    dgup(2,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(42)-2*g_contr(3,1)*g_contr(2,2)*Q(43)-2*g_contr(3,1)*g_contr(3,2)*Q(44)-2*g_contr(3,2)*g_contr(1,2)*Q(43)-2*g_contr(3,2)*g_contr(2,2)*Q(45)-2*g_contr(3,2)**2*Q(46)-2*g_contr(3,3)*g_contr(1,2)*Q(44)-2*g_contr(2,2)*g_contr(3,3)*Q(46)-2*g_contr(3,3)*g_contr(3,2)*Q(47)
    dgup(2,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(42)-2*g_contr(3,1)*g_contr(2,3)*Q(43)-2*g_contr(3,1)*g_contr(3,3)*Q(44)-2*g_contr(3,2)*g_contr(1,3)*Q(43)-2*g_contr(3,2)*g_contr(2,3)*Q(45)-2*g_contr(3,2)*g_contr(3,3)*Q(46)-2*g_contr(3,3)*g_contr(1,3)*Q(44)-2*g_contr(3,3)*g_contr(2,3)*Q(46)-2*g_contr(3,3)**2*Q(47)
    dgup(3,1,1) = -2*g_contr(1,1)**2*Q(48)-2*g_contr(1,1)*g_contr(2,1)*Q(49)-2*g_contr(1,1)*g_contr(3,1)*Q(50)-2*g_contr(1,2)*g_contr(1,1)*Q(49)-2*g_contr(1,2)*g_contr(2,1)*Q(51)-2*g_contr(1,2)*g_contr(3,1)*Q(52)-2*g_contr(1,3)*g_contr(1,1)*Q(50)-2*g_contr(1,3)*g_contr(2,1)*Q(52)-2*g_contr(1,3)*g_contr(3,1)*Q(53)
    dgup(3,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(48)-2*g_contr(1,1)*g_contr(2,2)*Q(49)-2*g_contr(1,1)*g_contr(3,2)*Q(50)-2*g_contr(1,2)**2*Q(49)-2*g_contr(1,2)*g_contr(2,2)*Q(51)-2*g_contr(1,2)*g_contr(3,2)*Q(52)-2*g_contr(1,3)*g_contr(1,2)*Q(50)-2*g_contr(1,3)*g_contr(2,2)*Q(52)-2*g_contr(1,3)*g_contr(3,2)*Q(53)
    dgup(3,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(48)-2*g_contr(1,1)*g_contr(2,3)*Q(49)-2*g_contr(1,1)*g_contr(3,3)*Q(50)-2*g_contr(1,2)*g_contr(1,3)*Q(49)-2*g_contr(1,2)*g_contr(2,3)*Q(51)-2*g_contr(1,2)*g_contr(3,3)*Q(52)-2*g_contr(1,3)**2*Q(50)-2*g_contr(1,3)*g_contr(2,3)*Q(52)-2*g_contr(1,3)*g_contr(3,3)*Q(53)
    dgup(3,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(48)-2*g_contr(2,1)**2*Q(49)-2*g_contr(2,1)*g_contr(3,1)*Q(50)-2*g_contr(1,1)*g_contr(2,2)*Q(49)-2*g_contr(2,2)*g_contr(2,1)*Q(51)-2*g_contr(2,2)*g_contr(3,1)*Q(52)-2*g_contr(2,3)*g_contr(1,1)*Q(50)-2*g_contr(2,3)*g_contr(2,1)*Q(52)-2*g_contr(2,3)*g_contr(3,1)*Q(53)
    dgup(3,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(48)-2*g_contr(2,1)*g_contr(2,2)*Q(49)-2*g_contr(2,1)*g_contr(3,2)*Q(50)-2*g_contr(2,2)*g_contr(1,2)*Q(49)-2*g_contr(2,2)**2*Q(51)-2*g_contr(2,2)*g_contr(3,2)*Q(52)-2*g_contr(2,3)*g_contr(1,2)*Q(50)-2*g_contr(2,3)*g_contr(2,2)*Q(52)-2*g_contr(2,3)*g_contr(3,2)*Q(53)
    dgup(3,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(48)-2*g_contr(2,1)*g_contr(2,3)*Q(49)-2*g_contr(2,1)*g_contr(3,3)*Q(50)-2*g_contr(2,2)*g_contr(1,3)*Q(49)-2*g_contr(2,2)*g_contr(2,3)*Q(51)-2*g_contr(2,2)*g_contr(3,3)*Q(52)-2*g_contr(2,3)*g_contr(1,3)*Q(50)-2*g_contr(2,3)**2*Q(52)-2*g_contr(2,3)*g_contr(3,3)*Q(53)
    dgup(3,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(48)-2*g_contr(3,1)*g_contr(2,1)*Q(49)-2*g_contr(3,1)**2*Q(50)-2*g_contr(3,2)*g_contr(1,1)*Q(49)-2*g_contr(3,2)*g_contr(2,1)*Q(51)-2*g_contr(3,2)*g_contr(3,1)*Q(52)-2*g_contr(1,1)*g_contr(3,3)*Q(50)-2*g_contr(3,3)*g_contr(2,1)*Q(52)-2*g_contr(3,3)*g_contr(3,1)*Q(53)
    dgup(3,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(48)-2*g_contr(3,1)*g_contr(2,2)*Q(49)-2*g_contr(3,1)*g_contr(3,2)*Q(50)-2*g_contr(3,2)*g_contr(1,2)*Q(49)-2*g_contr(3,2)*g_contr(2,2)*Q(51)-2*g_contr(3,2)**2*Q(52)-2*g_contr(3,3)*g_contr(1,2)*Q(50)-2*g_contr(2,2)*g_contr(3,3)*Q(52)-2*g_contr(3,3)*g_contr(3,2)*Q(53)
    dgup(3,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(48)-2*g_contr(3,1)*g_contr(2,3)*Q(49)-2*g_contr(3,1)*g_contr(3,3)*Q(50)-2*g_contr(3,2)*g_contr(1,3)*Q(49)-2*g_contr(3,2)*g_contr(2,3)*Q(51)-2*g_contr(3,2)*g_contr(3,3)*Q(52)-2*g_contr(3,3)*g_contr(1,3)*Q(50)-2*g_contr(3,3)*g_contr(2,3)*Q(52)-2*g_contr(3,3)**2*Q(53)

    dChristoffelSrc(1,1,1,1) = dgup(1,1,1)*Q(36)+dgup(1,1,2)*(2*Q(37)-Q(42))+dgup(1,1,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(1,1,1,2) = dgup(1,2,1)*Q(36)+dgup(1,2,2)*(2*Q(37)-Q(42))+dgup(1,2,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(1,1,1,3) = dgup(1,3,1)*Q(36)+dgup(1,3,2)*(2*Q(37)-Q(42))+dgup(1,3,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(1,1,2,1) = dgup(1,1,1)*Q(42)+dgup(1,1,2)*Q(39)+dgup(1,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,1,2,2) = dgup(1,2,1)*Q(42)+dgup(1,2,2)*Q(39)+dgup(1,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,1,2,3) = dgup(1,3,1)*Q(42)+dgup(1,3,2)*Q(39)+dgup(1,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,1,3,1) = dgup(1,1,1)*Q(48)+dgup(1,1,2)*(Q(49)+Q(40)-Q(44))+dgup(1,1,3)*Q(41)
    dChristoffelSrc(1,1,3,2) = dgup(1,2,1)*Q(48)+dgup(1,2,2)*(Q(49)+Q(40)-Q(44))+dgup(1,2,3)*Q(41)
    dChristoffelSrc(1,1,3,3) = dgup(1,3,1)*Q(48)+dgup(1,3,2)*(Q(49)+Q(40)-Q(44))+dgup(1,3,3)*Q(41)
    dChristoffelSrc(1,2,1,1) = dgup(1,1,1)*Q(42)+dgup(1,1,2)*Q(39)+dgup(1,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,2,1,2) = dgup(1,2,1)*Q(42)+dgup(1,2,2)*Q(39)+dgup(1,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,2,1,3) = dgup(1,3,1)*Q(42)+dgup(1,3,2)*Q(39)+dgup(1,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(1,2,2,1) = dgup(1,1,1)*(2*Q(43)-Q(39))+dgup(1,1,2)*Q(45)+dgup(1,1,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(1,2,2,2) = dgup(1,2,1)*(2*Q(43)-Q(39))+dgup(1,2,2)*Q(45)+dgup(1,2,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(1,2,2,3) = dgup(1,3,1)*(2*Q(43)-Q(39))+dgup(1,3,2)*Q(45)+dgup(1,3,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(1,2,3,1) = dgup(1,1,1)*(Q(49)+Q(44)-Q(40))+dgup(1,1,2)*Q(51)+dgup(1,1,3)*Q(47)
    dChristoffelSrc(1,2,3,2) = dgup(1,2,1)*(Q(49)+Q(44)-Q(40))+dgup(1,2,2)*Q(51)+dgup(1,2,3)*Q(47)
    dChristoffelSrc(1,2,3,3) = dgup(1,3,1)*(Q(49)+Q(44)-Q(40))+dgup(1,3,2)*Q(51)+dgup(1,3,3)*Q(47)
    dChristoffelSrc(1,3,1,1) = dgup(1,1,1)*Q(48)+dgup(1,1,2)*(Q(49)+Q(40)-Q(44))+dgup(1,1,3)*Q(41)
    dChristoffelSrc(1,3,1,2) = dgup(1,2,1)*Q(48)+dgup(1,2,2)*(Q(49)+Q(40)-Q(44))+dgup(1,2,3)*Q(41)
    dChristoffelSrc(1,3,1,3) = dgup(1,3,1)*Q(48)+dgup(1,3,2)*(Q(49)+Q(40)-Q(44))+dgup(1,3,3)*Q(41)
    dChristoffelSrc(1,3,2,1) = dgup(1,1,1)*(Q(49)+Q(44)-Q(40))+dgup(1,1,2)*Q(51)+dgup(1,1,3)*Q(47)
    dChristoffelSrc(1,3,2,2) = dgup(1,2,1)*(Q(49)+Q(44)-Q(40))+dgup(1,2,2)*Q(51)+dgup(1,2,3)*Q(47)
    dChristoffelSrc(1,3,2,3) = dgup(1,3,1)*(Q(49)+Q(44)-Q(40))+dgup(1,3,2)*Q(51)+dgup(1,3,3)*Q(47)
    dChristoffelSrc(1,3,3,1) = dgup(1,1,1)*(2*Q(50)-Q(41))+dgup(1,1,2)*(2*Q(52)-Q(47))+dgup(1,1,3)*Q(53)
    dChristoffelSrc(1,3,3,2) = dgup(1,2,1)*(2*Q(50)-Q(41))+dgup(1,2,2)*(2*Q(52)-Q(47))+dgup(1,2,3)*Q(53)
    dChristoffelSrc(1,3,3,3) = dgup(1,3,1)*(2*Q(50)-Q(41))+dgup(1,3,2)*(2*Q(52)-Q(47))+dgup(1,3,3)*Q(53)
    dChristoffelSrc(2,1,1,1) = dgup(2,1,1)*Q(36)+dgup(2,1,2)*(2*Q(37)-Q(42))+dgup(2,1,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(2,1,1,2) = dgup(2,2,1)*Q(36)+dgup(2,2,2)*(2*Q(37)-Q(42))+dgup(2,2,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(2,1,1,3) = dgup(2,3,1)*Q(36)+dgup(2,3,2)*(2*Q(37)-Q(42))+dgup(2,3,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(2,1,2,1) = dgup(2,1,1)*Q(42)+dgup(2,1,2)*Q(39)+dgup(2,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,1,2,2) = dgup(2,2,1)*Q(42)+dgup(2,2,2)*Q(39)+dgup(2,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,1,2,3) = dgup(2,3,1)*Q(42)+dgup(2,3,2)*Q(39)+dgup(2,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,1,3,1) = dgup(2,1,1)*Q(48)+dgup(2,1,2)*(Q(49)+Q(40)-Q(44))+dgup(2,1,3)*Q(41)
    dChristoffelSrc(2,1,3,2) = dgup(2,2,1)*Q(48)+dgup(2,2,2)*(Q(49)+Q(40)-Q(44))+dgup(2,2,3)*Q(41)
    dChristoffelSrc(2,1,3,3) = dgup(2,3,1)*Q(48)+dgup(2,3,2)*(Q(49)+Q(40)-Q(44))+dgup(2,3,3)*Q(41)
    dChristoffelSrc(2,2,1,1) = dgup(2,1,1)*Q(42)+dgup(2,1,2)*Q(39)+dgup(2,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,2,1,2) = dgup(2,2,1)*Q(42)+dgup(2,2,2)*Q(39)+dgup(2,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,2,1,3) = dgup(2,3,1)*Q(42)+dgup(2,3,2)*Q(39)+dgup(2,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(2,2,2,1) = dgup(2,1,1)*(2*Q(43)-Q(39))+dgup(2,1,2)*Q(45)+dgup(2,1,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(2,2,2,2) = dgup(2,2,1)*(2*Q(43)-Q(39))+dgup(2,2,2)*Q(45)+dgup(2,2,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(2,2,2,3) = dgup(2,3,1)*(2*Q(43)-Q(39))+dgup(2,3,2)*Q(45)+dgup(2,3,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(2,2,3,1) = dgup(2,1,1)*(Q(49)+Q(44)-Q(40))+dgup(2,1,2)*Q(51)+dgup(2,1,3)*Q(47)
    dChristoffelSrc(2,2,3,2) = dgup(2,2,1)*(Q(49)+Q(44)-Q(40))+dgup(2,2,2)*Q(51)+dgup(2,2,3)*Q(47)
    dChristoffelSrc(2,2,3,3) = dgup(2,3,1)*(Q(49)+Q(44)-Q(40))+dgup(2,3,2)*Q(51)+dgup(2,3,3)*Q(47)
    dChristoffelSrc(2,3,1,1) = dgup(2,1,1)*Q(48)+dgup(2,1,2)*(Q(49)+Q(40)-Q(44))+dgup(2,1,3)*Q(41)
    dChristoffelSrc(2,3,1,2) = dgup(2,2,1)*Q(48)+dgup(2,2,2)*(Q(49)+Q(40)-Q(44))+dgup(2,2,3)*Q(41)
    dChristoffelSrc(2,3,1,3) = dgup(2,3,1)*Q(48)+dgup(2,3,2)*(Q(49)+Q(40)-Q(44))+dgup(2,3,3)*Q(41)
    dChristoffelSrc(2,3,2,1) = dgup(2,1,1)*(Q(49)+Q(44)-Q(40))+dgup(2,1,2)*Q(51)+dgup(2,1,3)*Q(47)
    dChristoffelSrc(2,3,2,2) = dgup(2,2,1)*(Q(49)+Q(44)-Q(40))+dgup(2,2,2)*Q(51)+dgup(2,2,3)*Q(47)
    dChristoffelSrc(2,3,2,3) = dgup(2,3,1)*(Q(49)+Q(44)-Q(40))+dgup(2,3,2)*Q(51)+dgup(2,3,3)*Q(47)
    dChristoffelSrc(2,3,3,1) = dgup(2,1,1)*(2*Q(50)-Q(41))+dgup(2,1,2)*(2*Q(52)-Q(47))+dgup(2,1,3)*Q(53)
    dChristoffelSrc(2,3,3,2) = dgup(2,2,1)*(2*Q(50)-Q(41))+dgup(2,2,2)*(2*Q(52)-Q(47))+dgup(2,2,3)*Q(53)
    dChristoffelSrc(2,3,3,3) = dgup(2,3,1)*(2*Q(50)-Q(41))+dgup(2,3,2)*(2*Q(52)-Q(47))+dgup(2,3,3)*Q(53)
    dChristoffelSrc(3,1,1,1) = dgup(3,1,1)*Q(36)+dgup(3,1,2)*(2*Q(37)-Q(42))+dgup(3,1,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(3,1,1,2) = dgup(3,2,1)*Q(36)+dgup(3,2,2)*(2*Q(37)-Q(42))+dgup(3,2,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(3,1,1,3) = dgup(3,3,1)*Q(36)+dgup(3,3,2)*(2*Q(37)-Q(42))+dgup(3,3,3)*(2*Q(38)-Q(48))
    dChristoffelSrc(3,1,2,1) = dgup(3,1,1)*Q(42)+dgup(3,1,2)*Q(39)+dgup(3,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,1,2,2) = dgup(3,2,1)*Q(42)+dgup(3,2,2)*Q(39)+dgup(3,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,1,2,3) = dgup(3,3,1)*Q(42)+dgup(3,3,2)*Q(39)+dgup(3,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,1,3,1) = dgup(3,1,1)*Q(48)+dgup(3,1,2)*(Q(49)+Q(40)-Q(44))+dgup(3,1,3)*Q(41)
    dChristoffelSrc(3,1,3,2) = dgup(3,2,1)*Q(48)+dgup(3,2,2)*(Q(49)+Q(40)-Q(44))+dgup(3,2,3)*Q(41)
    dChristoffelSrc(3,1,3,3) = dgup(3,3,1)*Q(48)+dgup(3,3,2)*(Q(49)+Q(40)-Q(44))+dgup(3,3,3)*Q(41)
    dChristoffelSrc(3,2,1,1) = dgup(3,1,1)*Q(42)+dgup(3,1,2)*Q(39)+dgup(3,1,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,2,1,2) = dgup(3,2,1)*Q(42)+dgup(3,2,2)*Q(39)+dgup(3,2,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,2,1,3) = dgup(3,3,1)*Q(42)+dgup(3,3,2)*Q(39)+dgup(3,3,3)*(Q(44)+Q(40)-Q(49))
    dChristoffelSrc(3,2,2,1) = dgup(3,1,1)*(2*Q(43)-Q(39))+dgup(3,1,2)*Q(45)+dgup(3,1,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(3,2,2,2) = dgup(3,2,1)*(2*Q(43)-Q(39))+dgup(3,2,2)*Q(45)+dgup(3,2,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(3,2,2,3) = dgup(3,3,1)*(2*Q(43)-Q(39))+dgup(3,3,2)*Q(45)+dgup(3,3,3)*(2*Q(46)-Q(51))
    dChristoffelSrc(3,2,3,1) = dgup(3,1,1)*(Q(49)+Q(44)-Q(40))+dgup(3,1,2)*Q(51)+dgup(3,1,3)*Q(47)
    dChristoffelSrc(3,2,3,2) = dgup(3,2,1)*(Q(49)+Q(44)-Q(40))+dgup(3,2,2)*Q(51)+dgup(3,2,3)*Q(47)
    dChristoffelSrc(3,2,3,3) = dgup(3,3,1)*(Q(49)+Q(44)-Q(40))+dgup(3,3,2)*Q(51)+dgup(3,3,3)*Q(47)
    dChristoffelSrc(3,3,1,1) = dgup(3,1,1)*Q(48)+dgup(3,1,2)*(Q(49)+Q(40)-Q(44))+dgup(3,1,3)*Q(41)
    dChristoffelSrc(3,3,1,2) = dgup(3,2,1)*Q(48)+dgup(3,2,2)*(Q(49)+Q(40)-Q(44))+dgup(3,2,3)*Q(41)
    dChristoffelSrc(3,3,1,3) = dgup(3,3,1)*Q(48)+dgup(3,3,2)*(Q(49)+Q(40)-Q(44))+dgup(3,3,3)*Q(41)
    dChristoffelSrc(3,3,2,1) = dgup(3,1,1)*(Q(49)+Q(44)-Q(40))+dgup(3,1,2)*Q(51)+dgup(3,1,3)*Q(47)
    dChristoffelSrc(3,3,2,2) = dgup(3,2,1)*(Q(49)+Q(44)-Q(40))+dgup(3,2,2)*Q(51)+dgup(3,2,3)*Q(47)
    dChristoffelSrc(3,3,2,3) = dgup(3,3,1)*(Q(49)+Q(44)-Q(40))+dgup(3,3,2)*Q(51)+dgup(3,3,3)*Q(47)
    dChristoffelSrc(3,3,3,1) = dgup(3,1,1)*(2*Q(50)-Q(41))+dgup(3,1,2)*(2*Q(52)-Q(47))+dgup(3,1,3)*Q(53)
    dChristoffelSrc(3,3,3,2) = dgup(3,2,1)*(2*Q(50)-Q(41))+dgup(3,2,2)*(2*Q(52)-Q(47))+dgup(3,2,3)*Q(53)
    dChristoffelSrc(3,3,3,3) = dgup(3,3,1)*(2*Q(50)-Q(41))+dgup(3,3,2)*(2*Q(52)-Q(47))+dgup(3,3,3)*Q(53) 
      
    RiemannSrc(1,1,1,1) = 0
    RiemannSrc(1,1,1,2) = 0
    RiemannSrc(1,1,1,3) = 0
    RiemannSrc(1,1,2,1) = dChristoffelSrc(1,1,2,1)-dChristoffelSrc(2,1,1,1)+Christoffel(1,2,2)*Christoffel(2,1,1)-Christoffel(1,1,2)*Christoffel(2,2,1)+Christoffel(1,2,3)*Christoffel(3,1,1)-Christoffel(1,1,3)*Christoffel(3,2,1)
    RiemannSrc(1,1,2,2) = dChristoffelSrc(1,1,2,2)-dChristoffelSrc(2,1,1,2)+Christoffel(1,2,1)*Christoffel(1,1,2)-Christoffel(1,1,1)*Christoffel(1,2,2)+Christoffel(1,2,2)*Christoffel(2,1,2)-Christoffel(1,1,2)*Christoffel(2,2,2)+Christoffel(1,2,3)*Christoffel(3,1,2)-Christoffel(1,1,3)*Christoffel(3,2,2)
    RiemannSrc(1,1,2,3) = dChristoffelSrc(1,1,2,3)-dChristoffelSrc(2,1,1,3)+Christoffel(1,2,1)*Christoffel(1,1,3)-Christoffel(1,1,1)*Christoffel(1,2,3)+Christoffel(1,2,2)*Christoffel(2,1,3)-Christoffel(1,1,2)*Christoffel(2,2,3)+Christoffel(1,2,3)*Christoffel(3,1,3)-Christoffel(1,1,3)*Christoffel(3,2,3)
    RiemannSrc(1,1,3,1) = dChristoffelSrc(1,1,3,1)-dChristoffelSrc(3,1,1,1)+Christoffel(1,3,2)*Christoffel(2,1,1)-Christoffel(1,1,2)*Christoffel(2,3,1)+Christoffel(1,3,3)*Christoffel(3,1,1)-Christoffel(1,1,3)*Christoffel(3,3,1)
    RiemannSrc(1,1,3,2) = dChristoffelSrc(1,1,3,2)-dChristoffelSrc(3,1,1,2)+Christoffel(1,3,1)*Christoffel(1,1,2)-Christoffel(1,1,1)*Christoffel(1,3,2)+Christoffel(1,3,2)*Christoffel(2,1,2)-Christoffel(1,1,2)*Christoffel(2,3,2)+Christoffel(1,3,3)*Christoffel(3,1,2)-Christoffel(1,1,3)*Christoffel(3,3,2)
    RiemannSrc(1,1,3,3) = dChristoffelSrc(1,1,3,3)-dChristoffelSrc(3,1,1,3)+Christoffel(1,3,1)*Christoffel(1,1,3)-Christoffel(1,1,1)*Christoffel(1,3,3)+Christoffel(1,3,2)*Christoffel(2,1,3)-Christoffel(1,1,2)*Christoffel(2,3,3)+Christoffel(1,3,3)*Christoffel(3,1,3)-Christoffel(1,1,3)*Christoffel(3,3,3)
    RiemannSrc(1,2,1,1) = dChristoffelSrc(2,1,1,1)-dChristoffelSrc(1,1,2,1)+Christoffel(1,1,2)*Christoffel(2,2,1)-Christoffel(1,2,2)*Christoffel(2,1,1)+Christoffel(1,1,3)*Christoffel(3,2,1)-Christoffel(1,2,3)*Christoffel(3,1,1)
    RiemannSrc(1,2,1,2) = dChristoffelSrc(2,1,1,2)-dChristoffelSrc(1,1,2,2)+Christoffel(1,1,1)*Christoffel(1,2,2)-Christoffel(1,2,1)*Christoffel(1,1,2)+Christoffel(1,1,2)*Christoffel(2,2,2)-Christoffel(1,2,2)*Christoffel(2,1,2)+Christoffel(1,1,3)*Christoffel(3,2,2)-Christoffel(1,2,3)*Christoffel(3,1,2)
    RiemannSrc(1,2,1,3) = dChristoffelSrc(2,1,1,3)-dChristoffelSrc(1,1,2,3)+Christoffel(1,1,1)*Christoffel(1,2,3)-Christoffel(1,2,1)*Christoffel(1,1,3)+Christoffel(1,1,2)*Christoffel(2,2,3)-Christoffel(1,2,2)*Christoffel(2,1,3)+Christoffel(1,1,3)*Christoffel(3,2,3)-Christoffel(1,2,3)*Christoffel(3,1,3)
    RiemannSrc(1,2,2,1) = 0
    RiemannSrc(1,2,2,2) = 0
    RiemannSrc(1,2,2,3) = 0
    RiemannSrc(1,2,3,1) = dChristoffelSrc(2,1,3,1)-dChristoffelSrc(3,1,2,1)+Christoffel(1,3,2)*Christoffel(2,2,1)-Christoffel(1,2,2)*Christoffel(2,3,1)+Christoffel(1,3,3)*Christoffel(3,2,1)-Christoffel(1,2,3)*Christoffel(3,3,1)
    RiemannSrc(1,2,3,2) = dChristoffelSrc(2,1,3,2)-dChristoffelSrc(3,1,2,2)+Christoffel(1,3,1)*Christoffel(1,2,2)-Christoffel(1,2,1)*Christoffel(1,3,2)+Christoffel(1,3,2)*Christoffel(2,2,2)-Christoffel(1,2,2)*Christoffel(2,3,2)+Christoffel(1,3,3)*Christoffel(3,2,2)-Christoffel(1,2,3)*Christoffel(3,3,2)
    RiemannSrc(1,2,3,3) = dChristoffelSrc(2,1,3,3)-dChristoffelSrc(3,1,2,3)+Christoffel(1,3,1)*Christoffel(1,2,3)-Christoffel(1,2,1)*Christoffel(1,3,3)+Christoffel(1,3,2)*Christoffel(2,2,3)-Christoffel(1,2,2)*Christoffel(2,3,3)+Christoffel(1,3,3)*Christoffel(3,2,3)-Christoffel(1,2,3)*Christoffel(3,3,3)
    RiemannSrc(1,3,1,1) = dChristoffelSrc(3,1,1,1)-dChristoffelSrc(1,1,3,1)+Christoffel(1,1,2)*Christoffel(2,3,1)-Christoffel(1,3,2)*Christoffel(2,1,1)+Christoffel(1,1,3)*Christoffel(3,3,1)-Christoffel(1,3,3)*Christoffel(3,1,1)
    RiemannSrc(1,3,1,2) = dChristoffelSrc(3,1,1,2)-dChristoffelSrc(1,1,3,2)+Christoffel(1,1,1)*Christoffel(1,3,2)-Christoffel(1,3,1)*Christoffel(1,1,2)+Christoffel(1,1,2)*Christoffel(2,3,2)-Christoffel(1,3,2)*Christoffel(2,1,2)+Christoffel(1,1,3)*Christoffel(3,3,2)-Christoffel(1,3,3)*Christoffel(3,1,2)
    RiemannSrc(1,3,1,3) = dChristoffelSrc(3,1,1,3)-dChristoffelSrc(1,1,3,3)+Christoffel(1,1,1)*Christoffel(1,3,3)-Christoffel(1,3,1)*Christoffel(1,1,3)+Christoffel(1,1,2)*Christoffel(2,3,3)-Christoffel(1,3,2)*Christoffel(2,1,3)+Christoffel(1,1,3)*Christoffel(3,3,3)-Christoffel(1,3,3)*Christoffel(3,1,3)
    RiemannSrc(1,3,2,1) = dChristoffelSrc(3,1,2,1)-dChristoffelSrc(2,1,3,1)+Christoffel(1,2,2)*Christoffel(2,3,1)-Christoffel(1,3,2)*Christoffel(2,2,1)+Christoffel(1,2,3)*Christoffel(3,3,1)-Christoffel(1,3,3)*Christoffel(3,2,1)
    RiemannSrc(1,3,2,2) = dChristoffelSrc(3,1,2,2)-dChristoffelSrc(2,1,3,2)+Christoffel(1,2,1)*Christoffel(1,3,2)-Christoffel(1,3,1)*Christoffel(1,2,2)+Christoffel(1,2,2)*Christoffel(2,3,2)-Christoffel(1,3,2)*Christoffel(2,2,2)+Christoffel(1,2,3)*Christoffel(3,3,2)-Christoffel(1,3,3)*Christoffel(3,2,2)
    RiemannSrc(1,3,2,3) = dChristoffelSrc(3,1,2,3)-dChristoffelSrc(2,1,3,3)+Christoffel(1,2,1)*Christoffel(1,3,3)-Christoffel(1,3,1)*Christoffel(1,2,3)+Christoffel(1,2,2)*Christoffel(2,3,3)-Christoffel(1,3,2)*Christoffel(2,2,3)+Christoffel(1,2,3)*Christoffel(3,3,3)-Christoffel(1,3,3)*Christoffel(3,2,3)
    RiemannSrc(1,3,3,1) = 0
    RiemannSrc(1,3,3,2) = 0
    RiemannSrc(1,3,3,3) = 0
    RiemannSrc(2,1,1,1) = 0
    RiemannSrc(2,1,1,2) = 0
    RiemannSrc(2,1,1,3) = 0
    RiemannSrc(2,1,2,1) = dChristoffelSrc(1,2,2,1)-dChristoffelSrc(2,2,1,1)+Christoffel(2,2,1)*Christoffel(1,1,1)-Christoffel(2,1,1)*Christoffel(1,2,1)+Christoffel(2,2,2)*Christoffel(2,1,1)-Christoffel(2,1,2)*Christoffel(2,2,1)+Christoffel(2,2,3)*Christoffel(3,1,1)-Christoffel(2,1,3)*Christoffel(3,2,1)
    RiemannSrc(2,1,2,2) = dChristoffelSrc(1,2,2,2)-dChristoffelSrc(2,2,1,2)+Christoffel(1,1,2)*Christoffel(2,2,1)-Christoffel(1,2,2)*Christoffel(2,1,1)+Christoffel(2,2,3)*Christoffel(3,1,2)-Christoffel(2,1,3)*Christoffel(3,2,2)
    RiemannSrc(2,1,2,3) = dChristoffelSrc(1,2,2,3)-dChristoffelSrc(2,2,1,3)+Christoffel(2,2,1)*Christoffel(1,1,3)-Christoffel(2,1,1)*Christoffel(1,2,3)+Christoffel(2,2,2)*Christoffel(2,1,3)-Christoffel(2,1,2)*Christoffel(2,2,3)+Christoffel(2,2,3)*Christoffel(3,1,3)-Christoffel(2,1,3)*Christoffel(3,2,3)
    RiemannSrc(2,1,3,1) = dChristoffelSrc(1,2,3,1)-dChristoffelSrc(3,2,1,1)+Christoffel(2,3,1)*Christoffel(1,1,1)-Christoffel(2,1,1)*Christoffel(1,3,1)+Christoffel(2,3,2)*Christoffel(2,1,1)-Christoffel(2,1,2)*Christoffel(2,3,1)+Christoffel(2,3,3)*Christoffel(3,1,1)-Christoffel(2,1,3)*Christoffel(3,3,1)
    RiemannSrc(2,1,3,2) = dChristoffelSrc(1,2,3,2)-dChristoffelSrc(3,2,1,2)+Christoffel(1,1,2)*Christoffel(2,3,1)-Christoffel(1,3,2)*Christoffel(2,1,1)+Christoffel(2,3,3)*Christoffel(3,1,2)-Christoffel(2,1,3)*Christoffel(3,3,2)
    RiemannSrc(2,1,3,3) = dChristoffelSrc(1,2,3,3)-dChristoffelSrc(3,2,1,3)+Christoffel(2,3,1)*Christoffel(1,1,3)-Christoffel(2,1,1)*Christoffel(1,3,3)+Christoffel(2,3,2)*Christoffel(2,1,3)-Christoffel(2,1,2)*Christoffel(2,3,3)+Christoffel(2,3,3)*Christoffel(3,1,3)-Christoffel(2,1,3)*Christoffel(3,3,3)
    RiemannSrc(2,2,1,1) = dChristoffelSrc(2,2,1,1)-dChristoffelSrc(1,2,2,1)+Christoffel(2,1,1)*Christoffel(1,2,1)-Christoffel(2,2,1)*Christoffel(1,1,1)+Christoffel(2,1,2)*Christoffel(2,2,1)-Christoffel(2,2,2)*Christoffel(2,1,1)+Christoffel(2,1,3)*Christoffel(3,2,1)-Christoffel(2,2,3)*Christoffel(3,1,1)
    RiemannSrc(2,2,1,2) = dChristoffelSrc(2,2,1,2)-dChristoffelSrc(1,2,2,2)+Christoffel(1,2,2)*Christoffel(2,1,1)-Christoffel(1,1,2)*Christoffel(2,2,1)+Christoffel(2,1,3)*Christoffel(3,2,2)-Christoffel(2,2,3)*Christoffel(3,1,2)
    RiemannSrc(2,2,1,3) = dChristoffelSrc(2,2,1,3)-dChristoffelSrc(1,2,2,3)+Christoffel(2,1,1)*Christoffel(1,2,3)-Christoffel(2,2,1)*Christoffel(1,1,3)+Christoffel(2,1,2)*Christoffel(2,2,3)-Christoffel(2,2,2)*Christoffel(2,1,3)+Christoffel(2,1,3)*Christoffel(3,2,3)-Christoffel(2,2,3)*Christoffel(3,1,3)
    RiemannSrc(2,2,2,1) = 0
    RiemannSrc(2,2,2,2) = 0
    RiemannSrc(2,2,2,3) = 0
    RiemannSrc(2,2,3,1) = dChristoffelSrc(2,2,3,1)-dChristoffelSrc(3,2,2,1)+Christoffel(2,3,1)*Christoffel(1,2,1)-Christoffel(2,2,1)*Christoffel(1,3,1)+Christoffel(2,3,2)*Christoffel(2,2,1)-Christoffel(2,2,2)*Christoffel(2,3,1)+Christoffel(2,3,3)*Christoffel(3,2,1)-Christoffel(2,2,3)*Christoffel(3,3,1)
    RiemannSrc(2,2,3,2) = dChristoffelSrc(2,2,3,2)-dChristoffelSrc(3,2,2,2)+Christoffel(1,2,2)*Christoffel(2,3,1)-Christoffel(1,3,2)*Christoffel(2,2,1)+Christoffel(2,3,3)*Christoffel(3,2,2)-Christoffel(2,2,3)*Christoffel(3,3,2)
    RiemannSrc(2,2,3,3) = dChristoffelSrc(2,2,3,3)-dChristoffelSrc(3,2,2,3)+Christoffel(2,3,1)*Christoffel(1,2,3)-Christoffel(2,2,1)*Christoffel(1,3,3)+Christoffel(2,3,2)*Christoffel(2,2,3)-Christoffel(2,2,2)*Christoffel(2,3,3)+Christoffel(2,3,3)*Christoffel(3,2,3)-Christoffel(2,2,3)*Christoffel(3,3,3)
    RiemannSrc(2,3,1,1) = dChristoffelSrc(3,2,1,1)-dChristoffelSrc(1,2,3,1)+Christoffel(2,1,1)*Christoffel(1,3,1)-Christoffel(2,3,1)*Christoffel(1,1,1)+Christoffel(2,1,2)*Christoffel(2,3,1)-Christoffel(2,3,2)*Christoffel(2,1,1)+Christoffel(2,1,3)*Christoffel(3,3,1)-Christoffel(2,3,3)*Christoffel(3,1,1)
    RiemannSrc(2,3,1,2) = dChristoffelSrc(3,2,1,2)-dChristoffelSrc(1,2,3,2)+Christoffel(1,3,2)*Christoffel(2,1,1)-Christoffel(1,1,2)*Christoffel(2,3,1)+Christoffel(2,1,3)*Christoffel(3,3,2)-Christoffel(2,3,3)*Christoffel(3,1,2)
    RiemannSrc(2,3,1,3) = dChristoffelSrc(3,2,1,3)-dChristoffelSrc(1,2,3,3)+Christoffel(2,1,1)*Christoffel(1,3,3)-Christoffel(2,3,1)*Christoffel(1,1,3)+Christoffel(2,1,2)*Christoffel(2,3,3)-Christoffel(2,3,2)*Christoffel(2,1,3)+Christoffel(2,1,3)*Christoffel(3,3,3)-Christoffel(2,3,3)*Christoffel(3,1,3)
    RiemannSrc(2,3,2,1) = dChristoffelSrc(3,2,2,1)-dChristoffelSrc(2,2,3,1)+Christoffel(2,2,1)*Christoffel(1,3,1)-Christoffel(2,3,1)*Christoffel(1,2,1)+Christoffel(2,2,2)*Christoffel(2,3,1)-Christoffel(2,3,2)*Christoffel(2,2,1)+Christoffel(2,2,3)*Christoffel(3,3,1)-Christoffel(2,3,3)*Christoffel(3,2,1)
    RiemannSrc(2,3,2,2) = dChristoffelSrc(3,2,2,2)-dChristoffelSrc(2,2,3,2)+Christoffel(1,3,2)*Christoffel(2,2,1)-Christoffel(1,2,2)*Christoffel(2,3,1)+Christoffel(2,2,3)*Christoffel(3,3,2)-Christoffel(2,3,3)*Christoffel(3,2,2)
    RiemannSrc(2,3,2,3) = dChristoffelSrc(3,2,2,3)-dChristoffelSrc(2,2,3,3)+Christoffel(2,2,1)*Christoffel(1,3,3)-Christoffel(2,3,1)*Christoffel(1,2,3)+Christoffel(2,2,2)*Christoffel(2,3,3)-Christoffel(2,3,2)*Christoffel(2,2,3)+Christoffel(2,2,3)*Christoffel(3,3,3)-Christoffel(2,3,3)*Christoffel(3,2,3)
    RiemannSrc(2,3,3,1) = 0
    RiemannSrc(2,3,3,2) = 0
    RiemannSrc(2,3,3,3) = 0
    RiemannSrc(3,1,1,1) = 0
    RiemannSrc(3,1,1,2) = 0
    RiemannSrc(3,1,1,3) = 0
    RiemannSrc(3,1,2,1) = dChristoffelSrc(1,3,2,1)-dChristoffelSrc(2,3,1,1)+Christoffel(3,2,1)*Christoffel(1,1,1)-Christoffel(3,1,1)*Christoffel(1,2,1)+Christoffel(3,2,2)*Christoffel(2,1,1)-Christoffel(3,1,2)*Christoffel(2,2,1)+Christoffel(3,2,3)*Christoffel(3,1,1)-Christoffel(3,1,3)*Christoffel(3,2,1)
    RiemannSrc(3,1,2,2) = dChristoffelSrc(1,3,2,2)-dChristoffelSrc(2,3,1,2)+Christoffel(3,2,1)*Christoffel(1,1,2)-Christoffel(3,1,1)*Christoffel(1,2,2)+Christoffel(3,2,2)*Christoffel(2,1,2)-Christoffel(3,1,2)*Christoffel(2,2,2)+Christoffel(3,2,3)*Christoffel(3,1,2)-Christoffel(3,1,3)*Christoffel(3,2,2)
    RiemannSrc(3,1,2,3) = dChristoffelSrc(1,3,2,3)-dChristoffelSrc(2,3,1,3)+Christoffel(1,1,3)*Christoffel(3,2,1)-Christoffel(1,2,3)*Christoffel(3,1,1)+Christoffel(2,1,3)*Christoffel(3,2,2)-Christoffel(2,2,3)*Christoffel(3,1,2)
    RiemannSrc(3,1,3,1) = dChristoffelSrc(1,3,3,1)-dChristoffelSrc(3,3,1,1)+Christoffel(3,3,1)*Christoffel(1,1,1)-Christoffel(3,1,1)*Christoffel(1,3,1)+Christoffel(3,3,2)*Christoffel(2,1,1)-Christoffel(3,1,2)*Christoffel(2,3,1)+Christoffel(3,3,3)*Christoffel(3,1,1)-Christoffel(3,1,3)*Christoffel(3,3,1)
    RiemannSrc(3,1,3,2) = dChristoffelSrc(1,3,3,2)-dChristoffelSrc(3,3,1,2)+Christoffel(3,3,1)*Christoffel(1,1,2)-Christoffel(3,1,1)*Christoffel(1,3,2)+Christoffel(3,3,2)*Christoffel(2,1,2)-Christoffel(3,1,2)*Christoffel(2,3,2)+Christoffel(3,3,3)*Christoffel(3,1,2)-Christoffel(3,1,3)*Christoffel(3,3,2)
    RiemannSrc(3,1,3,3) = dChristoffelSrc(1,3,3,3)-dChristoffelSrc(3,3,1,3)+Christoffel(1,1,3)*Christoffel(3,3,1)-Christoffel(1,3,3)*Christoffel(3,1,1)+Christoffel(2,1,3)*Christoffel(3,3,2)-Christoffel(2,3,3)*Christoffel(3,1,2)
    RiemannSrc(3,2,1,1) = dChristoffelSrc(2,3,1,1)-dChristoffelSrc(1,3,2,1)+Christoffel(3,1,1)*Christoffel(1,2,1)-Christoffel(3,2,1)*Christoffel(1,1,1)+Christoffel(3,1,2)*Christoffel(2,2,1)-Christoffel(3,2,2)*Christoffel(2,1,1)+Christoffel(3,1,3)*Christoffel(3,2,1)-Christoffel(3,2,3)*Christoffel(3,1,1)
    RiemannSrc(3,2,1,2) = dChristoffelSrc(2,3,1,2)-dChristoffelSrc(1,3,2,2)+Christoffel(3,1,1)*Christoffel(1,2,2)-Christoffel(3,2,1)*Christoffel(1,1,2)+Christoffel(3,1,2)*Christoffel(2,2,2)-Christoffel(3,2,2)*Christoffel(2,1,2)+Christoffel(3,1,3)*Christoffel(3,2,2)-Christoffel(3,2,3)*Christoffel(3,1,2)
    RiemannSrc(3,2,1,3) = dChristoffelSrc(2,3,1,3)-dChristoffelSrc(1,3,2,3)+Christoffel(1,2,3)*Christoffel(3,1,1)-Christoffel(1,1,3)*Christoffel(3,2,1)+Christoffel(2,2,3)*Christoffel(3,1,2)-Christoffel(2,1,3)*Christoffel(3,2,2)
    RiemannSrc(3,2,2,1) = 0
    RiemannSrc(3,2,2,2) = 0
    RiemannSrc(3,2,2,3) = 0
    RiemannSrc(3,2,3,1) = dChristoffelSrc(2,3,3,1)-dChristoffelSrc(3,3,2,1)+Christoffel(3,3,1)*Christoffel(1,2,1)-Christoffel(3,2,1)*Christoffel(1,3,1)+Christoffel(3,3,2)*Christoffel(2,2,1)-Christoffel(3,2,2)*Christoffel(2,3,1)+Christoffel(3,3,3)*Christoffel(3,2,1)-Christoffel(3,2,3)*Christoffel(3,3,1)
    RiemannSrc(3,2,3,2) = dChristoffelSrc(2,3,3,2)-dChristoffelSrc(3,3,2,2)+Christoffel(3,3,1)*Christoffel(1,2,2)-Christoffel(3,2,1)*Christoffel(1,3,2)+Christoffel(3,3,2)*Christoffel(2,2,2)-Christoffel(3,2,2)*Christoffel(2,3,2)+Christoffel(3,3,3)*Christoffel(3,2,2)-Christoffel(3,2,3)*Christoffel(3,3,2)
    RiemannSrc(3,2,3,3) = dChristoffelSrc(2,3,3,3)-dChristoffelSrc(3,3,2,3)+Christoffel(1,2,3)*Christoffel(3,3,1)-Christoffel(1,3,3)*Christoffel(3,2,1)+Christoffel(2,2,3)*Christoffel(3,3,2)-Christoffel(2,3,3)*Christoffel(3,2,2)
    RiemannSrc(3,3,1,1) = dChristoffelSrc(3,3,1,1)-dChristoffelSrc(1,3,3,1)+Christoffel(3,1,1)*Christoffel(1,3,1)-Christoffel(3,3,1)*Christoffel(1,1,1)+Christoffel(3,1,2)*Christoffel(2,3,1)-Christoffel(3,3,2)*Christoffel(2,1,1)+Christoffel(3,1,3)*Christoffel(3,3,1)-Christoffel(3,3,3)*Christoffel(3,1,1)
    RiemannSrc(3,3,1,2) = dChristoffelSrc(3,3,1,2)-dChristoffelSrc(1,3,3,2)+Christoffel(3,1,1)*Christoffel(1,3,2)-Christoffel(3,3,1)*Christoffel(1,1,2)+Christoffel(3,1,2)*Christoffel(2,3,2)-Christoffel(3,3,2)*Christoffel(2,1,2)+Christoffel(3,1,3)*Christoffel(3,3,2)-Christoffel(3,3,3)*Christoffel(3,1,2)
    RiemannSrc(3,3,1,3) = dChristoffelSrc(3,3,1,3)-dChristoffelSrc(1,3,3,3)+Christoffel(1,3,3)*Christoffel(3,1,1)-Christoffel(1,1,3)*Christoffel(3,3,1)+Christoffel(2,3,3)*Christoffel(3,1,2)-Christoffel(2,1,3)*Christoffel(3,3,2)
    RiemannSrc(3,3,2,1) = dChristoffelSrc(3,3,2,1)-dChristoffelSrc(2,3,3,1)+Christoffel(3,2,1)*Christoffel(1,3,1)-Christoffel(3,3,1)*Christoffel(1,2,1)+Christoffel(3,2,2)*Christoffel(2,3,1)-Christoffel(3,3,2)*Christoffel(2,2,1)+Christoffel(3,2,3)*Christoffel(3,3,1)-Christoffel(3,3,3)*Christoffel(3,2,1)
    RiemannSrc(3,3,2,2) = dChristoffelSrc(3,3,2,2)-dChristoffelSrc(2,3,3,2)+Christoffel(3,2,1)*Christoffel(1,3,2)-Christoffel(3,3,1)*Christoffel(1,2,2)+Christoffel(3,2,2)*Christoffel(2,3,2)-Christoffel(3,3,2)*Christoffel(2,2,2)+Christoffel(3,2,3)*Christoffel(3,3,2)-Christoffel(3,3,3)*Christoffel(3,2,2)
    RiemannSrc(3,3,2,3) = dChristoffelSrc(3,3,2,3)-dChristoffelSrc(2,3,3,3)+Christoffel(1,3,3)*Christoffel(3,2,1)-Christoffel(1,2,3)*Christoffel(3,3,1)+Christoffel(2,3,3)*Christoffel(3,2,2)-Christoffel(2,2,3)*Christoffel(3,3,2)
    RiemannSrc(3,3,3,1) = 0
    RiemannSrc(3,3,3,2) = 0
    RiemannSrc(3,3,3,3) = 0 
      
    RicciSrc(1,1) = RiemannSrc(1,1,1,1)+RiemannSrc(1,2,1,2)+RiemannSrc(1,3,1,3)
    RicciSrc(1,2) = RiemannSrc(1,1,2,1)+RiemannSrc(1,2,2,2)+RiemannSrc(1,3,2,3)
    RicciSrc(1,3) = RiemannSrc(1,1,3,1)+RiemannSrc(1,2,3,2)+RiemannSrc(1,3,3,3)
    RicciSrc(2,1) = RiemannSrc(2,1,1,1)+RiemannSrc(2,2,1,2)+RiemannSrc(2,3,1,3)
    RicciSrc(2,2) = RiemannSrc(2,1,2,1)+RiemannSrc(2,2,2,2)+RiemannSrc(2,3,2,3)
    RicciSrc(2,3) = RiemannSrc(2,1,3,1)+RiemannSrc(2,2,3,2)+RiemannSrc(2,3,3,3)
    RicciSrc(3,1) = RiemannSrc(3,1,1,1)+RiemannSrc(3,2,1,2)+RiemannSrc(3,3,1,3)
    RicciSrc(3,2) = RiemannSrc(3,1,2,1)+RiemannSrc(3,2,2,2)+RiemannSrc(3,3,2,3)
    RicciSrc(3,3) = RiemannSrc(3,1,3,1)+RiemannSrc(3,2,3,2)+RiemannSrc(3,3,3,3)
      
    Gammaup(1) = g_contr(1,1)*Christoffel(1,1,1)+g_contr(1,2)*Christoffel(1,2,1)+g_contr(1,3)*Christoffel(1,3,1)+g_contr(2,1)*Christoffel(2,1,1)+g_contr(2,2)*Christoffel(2,2,1)+g_contr(2,3)*Christoffel(2,3,1)+g_contr(3,1)*Christoffel(3,1,1)+g_contr(3,2)*Christoffel(3,2,1)+g_contr(3,3)*Christoffel(3,3,1)
    Gammaup(2) = g_contr(1,1)*Christoffel(1,1,2)+g_contr(1,2)*Christoffel(1,2,2)+g_contr(1,3)*Christoffel(1,3,2)+g_contr(2,1)*Christoffel(2,1,2)+g_contr(2,2)*Christoffel(2,2,2)+g_contr(2,3)*Christoffel(2,3,2)+g_contr(3,1)*Christoffel(3,1,2)+g_contr(3,2)*Christoffel(3,2,2)+g_contr(3,3)*Christoffel(3,3,2)
    Gammaup(3) = g_contr(1,1)*Christoffel(1,1,3)+g_contr(1,2)*Christoffel(1,2,3)+g_contr(1,3)*Christoffel(1,3,3)+g_contr(2,1)*Christoffel(2,1,3)+g_contr(2,2)*Christoffel(2,2,3)+g_contr(2,3)*Christoffel(2,3,3)+g_contr(3,1)*Christoffel(3,1,3)+g_contr(3,2)*Christoffel(3,2,3)+g_contr(3,3)*Christoffel(3,3,3)

    dtmetric(1,1) = -2*Q(17)*Q(7)+2*g_cov(1,1)*Q(27)+2*Q(18)*Q(36)+2*g_cov(2,1)*Q(30)+2*Q(19)*Q(42)+2*g_cov(3,1)*Q(33)+2*Q(20)*Q(48)
    dtmetric(1,2) = -2*Q(17)*Q(8)+g_cov(1,2)*Q(27)+g_cov(1,1)*Q(28)+2*Q(18)*Q(37)+g_cov(2,2)*Q(30)+g_cov(2,1)*Q(31)+2*Q(19)*Q(43)+g_cov(3,2)*Q(33)+g_cov(3,1)*Q(34)+2*Q(20)*Q(49)
    dtmetric(1,3) = -2*Q(17)*Q(9)+g_cov(1,3)*Q(27)+g_cov(1,1)*Q(29)+2*Q(18)*Q(38)+g_cov(2,3)*Q(30)+g_cov(2,1)*Q(32)+2*Q(19)*Q(44)+g_cov(3,3)*Q(33)+g_cov(3,1)*Q(35)+2*Q(20)*Q(50)
    dtmetric(2,1) = -2*Q(17)*Q(8)+g_cov(1,2)*Q(27)+g_cov(1,1)*Q(28)+2*Q(18)*Q(37)+g_cov(2,2)*Q(30)+g_cov(2,1)*Q(31)+2*Q(19)*Q(43)+g_cov(3,2)*Q(33)+g_cov(3,1)*Q(34)+2*Q(20)*Q(49)
    dtmetric(2,2) = -2*Q(17)*Q(10)+2*g_cov(1,2)*Q(28)+2*Q(18)*Q(39)+2*g_cov(2,2)*Q(31)+2*Q(19)*Q(45)+2*g_cov(3,2)*Q(34)+2*Q(20)*Q(51)
    dtmetric(2,3) = -2*Q(17)*Q(11)+g_cov(1,3)*Q(28)+g_cov(1,2)*Q(29)+2*Q(18)*Q(40)+g_cov(2,3)*Q(31)+g_cov(2,2)*Q(32)+2*Q(19)*Q(46)+g_cov(3,3)*Q(34)+g_cov(3,2)*Q(35)+2*Q(20)*Q(52)
    dtmetric(3,1) = -2*Q(17)*Q(9)+g_cov(1,3)*Q(27)+g_cov(1,1)*Q(29)+2*Q(18)*Q(38)+g_cov(2,3)*Q(30)+g_cov(2,1)*Q(32)+2*Q(19)*Q(44)+g_cov(3,3)*Q(33)+g_cov(3,1)*Q(35)+2*Q(20)*Q(50)
    dtmetric(3,2) = -2*Q(17)*Q(11)+g_cov(1,3)*Q(28)+g_cov(1,2)*Q(29)+2*Q(18)*Q(40)+g_cov(2,3)*Q(31)+g_cov(2,2)*Q(32)+2*Q(19)*Q(46)+g_cov(3,3)*Q(34)+g_cov(3,2)*Q(35)+2*Q(20)*Q(52)
    dtmetric(3,3) = -2*Q(17)*Q(12)+2*g_cov(1,3)*Q(29)+2*Q(18)*Q(41)+2*g_cov(2,3)*Q(32)+2*Q(19)*Q(47)+2*g_cov(3,3)*Q(35)+2*Q(20)*Q(53)

    dtgup(1,1) = -g_contr(1,1)**2*dtmetric(1,1)-g_contr(1,1)*g_contr(1,2)*dtmetric(2,1)-g_contr(1,1)*g_contr(1,3)*dtmetric(3,1)-g_contr(1,2)*g_contr(1,1)*dtmetric(1,2)-g_contr(1,2)**2*dtmetric(2,2)-g_contr(1,2)*g_contr(1,3)*dtmetric(3,2)-g_contr(1,3)*g_contr(1,1)*dtmetric(1,3)-g_contr(1,3)*g_contr(1,2)*dtmetric(2,3)-g_contr(1,3)**2*dtmetric(3,3)
    dtgup(1,2) = -g_contr(1,1)*g_contr(2,1)*dtmetric(1,1)-g_contr(1,1)*g_contr(2,2)*dtmetric(2,1)-g_contr(1,1)*g_contr(2,3)*dtmetric(3,1)-g_contr(1,2)*g_contr(2,1)*dtmetric(1,2)-g_contr(1,2)*g_contr(2,2)*dtmetric(2,2)-g_contr(1,2)*g_contr(2,3)*dtmetric(3,2)-g_contr(1,3)*g_contr(2,1)*dtmetric(1,3)-g_contr(1,3)*g_contr(2,2)*dtmetric(2,3)-g_contr(1,3)*g_contr(2,3)*dtmetric(3,3)
    dtgup(1,3) = -g_contr(1,1)*g_contr(3,1)*dtmetric(1,1)-g_contr(1,1)*g_contr(3,2)*dtmetric(2,1)-g_contr(1,1)*g_contr(3,3)*dtmetric(3,1)-g_contr(1,2)*g_contr(3,1)*dtmetric(1,2)-g_contr(1,2)*g_contr(3,2)*dtmetric(2,2)-g_contr(1,2)*g_contr(3,3)*dtmetric(3,2)-g_contr(1,3)*g_contr(3,1)*dtmetric(1,3)-g_contr(1,3)*g_contr(3,2)*dtmetric(2,3)-g_contr(1,3)*g_contr(3,3)*dtmetric(3,3)
    dtgup(2,1) = -g_contr(1,1)*g_contr(2,1)*dtmetric(1,1)-g_contr(2,1)*g_contr(1,2)*dtmetric(2,1)-g_contr(2,1)*g_contr(1,3)*dtmetric(3,1)-g_contr(2,2)*g_contr(1,1)*dtmetric(1,2)-g_contr(1,2)*g_contr(2,2)*dtmetric(2,2)-g_contr(2,2)*g_contr(1,3)*dtmetric(3,2)-g_contr(2,3)*g_contr(1,1)*dtmetric(1,3)-g_contr(2,3)*g_contr(1,2)*dtmetric(2,3)-g_contr(1,3)*g_contr(2,3)*dtmetric(3,3)
    dtgup(2,2) = -g_contr(2,1)**2*dtmetric(1,1)-g_contr(2,1)*g_contr(2,2)*dtmetric(2,1)-g_contr(2,1)*g_contr(2,3)*dtmetric(3,1)-g_contr(2,2)*g_contr(2,1)*dtmetric(1,2)-g_contr(2,2)**2*dtmetric(2,2)-g_contr(2,2)*g_contr(2,3)*dtmetric(3,2)-g_contr(2,3)*g_contr(2,1)*dtmetric(1,3)-g_contr(2,3)*g_contr(2,2)*dtmetric(2,3)-g_contr(2,3)**2*dtmetric(3,3)
    dtgup(2,3) = -g_contr(2,1)*g_contr(3,1)*dtmetric(1,1)-g_contr(2,1)*g_contr(3,2)*dtmetric(2,1)-g_contr(2,1)*g_contr(3,3)*dtmetric(3,1)-g_contr(2,2)*g_contr(3,1)*dtmetric(1,2)-g_contr(2,2)*g_contr(3,2)*dtmetric(2,2)-g_contr(2,2)*g_contr(3,3)*dtmetric(3,2)-g_contr(2,3)*g_contr(3,1)*dtmetric(1,3)-g_contr(2,3)*g_contr(3,2)*dtmetric(2,3)-g_contr(2,3)*g_contr(3,3)*dtmetric(3,3)
    dtgup(3,1) = -g_contr(1,1)*g_contr(3,1)*dtmetric(1,1)-g_contr(3,1)*g_contr(1,2)*dtmetric(2,1)-g_contr(3,1)*g_contr(1,3)*dtmetric(3,1)-g_contr(3,2)*g_contr(1,1)*dtmetric(1,2)-g_contr(1,2)*g_contr(3,2)*dtmetric(2,2)-g_contr(3,2)*g_contr(1,3)*dtmetric(3,2)-g_contr(3,3)*g_contr(1,1)*dtmetric(1,3)-g_contr(3,3)*g_contr(1,2)*dtmetric(2,3)-g_contr(1,3)*g_contr(3,3)*dtmetric(3,3)
    dtgup(3,2) = -g_contr(2,1)*g_contr(3,1)*dtmetric(1,1)-g_contr(3,1)*g_contr(2,2)*dtmetric(2,1)-g_contr(3,1)*g_contr(2,3)*dtmetric(3,1)-g_contr(3,2)*g_contr(2,1)*dtmetric(1,2)-g_contr(2,2)*g_contr(3,2)*dtmetric(2,2)-g_contr(3,2)*g_contr(2,3)*dtmetric(3,2)-g_contr(3,3)*g_contr(2,1)*dtmetric(1,3)-g_contr(3,3)*g_contr(2,2)*dtmetric(2,3)-g_contr(2,3)*g_contr(3,3)*dtmetric(3,3)
    dtgup(3,3) = -g_contr(3,1)**2*dtmetric(1,1)-g_contr(3,1)*g_contr(3,2)*dtmetric(2,1)-g_contr(3,1)*g_contr(3,3)*dtmetric(3,1)-g_contr(3,2)*g_contr(3,1)*dtmetric(1,2)-g_contr(3,2)**2*dtmetric(2,2)-g_contr(3,2)*g_contr(3,3)*dtmetric(3,2)-g_contr(3,3)*g_contr(3,1)*dtmetric(1,3)-g_contr(3,3)*g_contr(3,2)*dtmetric(2,3)-g_contr(3,3)**2*dtmetric(3,3)

    s1 = -dtgup(1,2)*Q(42)+2*dtgup(1,2)*Q(37)-dtgup(1,3)*Q(48)+dtgup(1,1)*Q(36)+2*g_contr(1,1)*Q(38)*Q(33)-2*g_contr(1,2)*Q(17)*Q(24)*Q(8)+g_contr(1,2)*Q(17)*Q(25)*Q(7)+g_contr(1,3)*Q(17)*Q(26)*Q(7)-2*g_contr(1,3)*Q(17)*Q(24)*Q(9)-g_contr(1,1)*Q(17)*Q(24)*Q(7)+2*dtgup(1,3)*Q(38)+2*g_contr(1,3)*Q(41)*Q(33)+2*g_contr(1,3)*Q(38)*Q(35)+g_contr(1,3)*Q(29)*Q(36)+4*g_contr(1,3)*Q(27)*Q(38)-2*g_contr(1,3)*Q(48)*Q(27)-g_contr(1,3)*Q(32)*Q(42)-2*g_contr(1,3)*Q(49)*Q(30)
    dtChristoffelSRC(1,1,1) = s1-g_contr(1,3)*Q(35)*Q(48)+2*g_contr(1,2)*Q(40)*Q(33)+2*g_contr(1,2)*Q(38)*Q(34)+2*g_contr(1,2)*Q(33)*Q(49)+2*g_contr(1,2)*Q(39)*Q(30)+4*g_contr(1,2)*Q(27)*Q(37)+2*g_contr(1,2)*Q(37)*Q(31)+3*g_contr(1,1)*Q(27)*Q(36)+g_contr(1,1)*Q(33)*Q(48)+2*g_contr(1,3)*Q(30)*Q(44)-2*g_contr(1,2)*Q(44)*Q(33)+2*g_contr(1,3)*Q(40)*Q(30)+2*g_contr(1,3)*Q(37)*Q(32)+g_contr(1,2)*Q(28)*Q(36)-2*g_contr(1,2)*Q(42)*Q(27)-g_contr(1,2)*Q(31)*Q(42)-g_contr(1,2)*Q(34)*Q(48)+g_contr(1,1)*Q(30)*Q(42)+2*g_contr(1,1)*Q(37)*Q(30)
    s1 = -2*g_contr(2,3)*Q(49)*Q(30)+4*g_contr(2,3)*Q(27)*Q(38)+2*dtgup(2,2)*Q(37)-dtgup(2,3)*Q(48)+2*dtgup(2,3)*Q(38)+dtgup(2,1)*Q(36)-g_contr(2,3)*Q(35)*Q(48)-g_contr(2,1)*Q(17)*Q(24)*Q(7)-2*g_contr(2,2)*Q(17)*Q(24)*Q(8)+g_contr(2,2)*Q(17)*Q(25)*Q(7)+g_contr(2,3)*Q(17)*Q(26)*Q(7)-2*g_contr(2,3)*Q(17)*Q(24)*Q(9)-dtgup(2,2)*Q(42)+g_contr(2,1)*Q(30)*Q(42)+3*g_contr(2,1)*Q(27)*Q(36)+2*g_contr(2,1)*Q(37)*Q(30)+4*g_contr(2,2)*Q(27)*Q(37)+2*g_contr(2,2)*Q(39)*Q(30)
    dtChristoffelSRC(1,1,2) = s1+2*g_contr(2,2)*Q(33)*Q(49)+2*g_contr(2,2)*Q(38)*Q(34)+g_contr(2,1)*Q(33)*Q(48)+2*g_contr(2,1)*Q(38)*Q(33)+2*g_contr(2,2)*Q(37)*Q(31)+2*g_contr(2,2)*Q(40)*Q(33)+g_contr(2,2)*Q(28)*Q(36)-2*g_contr(2,2)*Q(42)*Q(27)-g_contr(2,2)*Q(31)*Q(42)-g_contr(2,2)*Q(34)*Q(48)-2*g_contr(2,2)*Q(44)*Q(33)+2*g_contr(2,3)*Q(30)*Q(44)+2*g_contr(2,3)*Q(37)*Q(32)+2*g_contr(2,3)*Q(40)*Q(30)+2*g_contr(2,3)*Q(38)*Q(35)+2*g_contr(2,3)*Q(41)*Q(33)+g_contr(2,3)*Q(29)*Q(36)-2*g_contr(2,3)*Q(48)*Q(27)-g_contr(2,3)*Q(32)*Q(42)
    s1 = 2*g_contr(3,2)*Q(38)*Q(34)+dtgup(3,1)*Q(36)-dtgup(3,3)*Q(48)-2*g_contr(3,2)*Q(17)*Q(24)*Q(8)+g_contr(3,2)*Q(17)*Q(25)*Q(7)+g_contr(3,3)*Q(17)*Q(26)*Q(7)-2*g_contr(3,3)*Q(17)*Q(24)*Q(9)+3*g_contr(3,1)*Q(27)*Q(36)-dtgup(3,2)*Q(42)+g_contr(3,1)*Q(30)*Q(42)+2*g_contr(3,1)*Q(37)*Q(30)+g_contr(3,1)*Q(33)*Q(48)+2*g_contr(3,1)*Q(38)*Q(33)+4*g_contr(3,2)*Q(27)*Q(37)+2*g_contr(3,2)*Q(37)*Q(31)+2*g_contr(3,2)*Q(39)*Q(30)+2*g_contr(3,2)*Q(33)*Q(49)+2*g_contr(3,2)*Q(40)*Q(33)
    dtChristoffelSRC(1,1,3) = s1+g_contr(3,2)*Q(28)*Q(36)-2*g_contr(3,2)*Q(42)*Q(27)-g_contr(3,2)*Q(31)*Q(42)-g_contr(3,2)*Q(34)*Q(48)-2*g_contr(3,2)*Q(44)*Q(33)+2*g_contr(3,3)*Q(30)*Q(44)+2*g_contr(3,3)*Q(37)*Q(32)+2*g_contr(3,3)*Q(40)*Q(30)+2*g_contr(3,3)*Q(38)*Q(35)+2*g_contr(3,3)*Q(41)*Q(33)+4*g_contr(3,3)*Q(27)*Q(38)+g_contr(3,3)*Q(29)*Q(36)-2*g_contr(3,3)*Q(48)*Q(27)+2*dtgup(3,3)*Q(38)-g_contr(3,3)*Q(32)*Q(42)-2*g_contr(3,3)*Q(49)*Q(30)-g_contr(3,3)*Q(35)*Q(48)+2*dtgup(3,2)*Q(37)-g_contr(3,1)*Q(17)*Q(24)*Q(7)
    s1 = g_contr(1,3)*Q(40)*Q(35)+2*g_contr(1,1)*Q(44)*Q(33)-g_contr(1,3)*Q(35)*Q(49)+g_contr(1,3)*Q(44)*Q(35)+g_contr(1,1)*Q(34)*Q(48)+2*g_contr(1,1)*Q(43)*Q(30)+2*g_contr(1,2)*Q(28)*Q(37)+g_contr(1,3)*Q(44)*Q(27)+2*g_contr(1,3)*Q(28)*Q(38)+g_contr(1,2)*Q(27)*Q(39)+g_contr(1,2)*Q(45)*Q(30)+2*g_contr(1,2)*Q(39)*Q(31)+g_contr(1,2)*Q(33)*Q(51)+2*g_contr(1,2)*Q(40)*Q(34)+g_contr(1,3)*Q(41)*Q(34)+g_contr(1,3)*Q(42)*Q(29)+g_contr(1,3)*Q(31)*Q(44)-g_contr(1,3)*Q(48)*Q(28)-g_contr(1,3)*Q(49)*Q(27)
    dtChristoffelSRC(1,2,1) = s1-g_contr(1,3)*Q(49)*Q(31)-g_contr(1,3)*Q(51)*Q(30)+g_contr(1,3)*Q(27)*Q(40)+g_contr(1,3)*Q(39)*Q(32)+g_contr(1,3)*Q(40)*Q(31)-dtgup(1,3)*Q(49)+dtgup(1,3)*Q(40)+dtgup(1,2)*Q(39)+dtgup(1,3)*Q(44)+2*g_contr(1,1)*Q(42)*Q(27)+g_contr(1,1)*Q(28)*Q(36)+g_contr(1,1)*Q(31)*Q(42)+2*g_contr(1,3)*Q(46)*Q(30)+g_contr(1,3)*Q(47)*Q(33)-g_contr(1,1)*Q(17)*Q(25)*Q(7)-g_contr(1,2)*Q(17)*Q(24)*Q(10)-g_contr(1,3)*Q(17)*Q(24)*Q(11)-g_contr(1,3)*Q(17)*Q(25)*Q(9)+g_contr(1,3)*Q(17)*Q(26)*Q(8)+dtgup(1,1)*Q(42)
    s1 = dtgup(2,1)*Q(42)-dtgup(2,3)*Q(49)+g_contr(2,2)*Q(45)*Q(30)-g_contr(2,3)*Q(17)*Q(24)*Q(11)-g_contr(2,3)*Q(17)*Q(25)*Q(9)+g_contr(2,3)*Q(17)*Q(26)*Q(8)+2*g_contr(2,1)*Q(44)*Q(33)+g_contr(2,1)*Q(28)*Q(36)+g_contr(2,3)*Q(47)*Q(33)+2*g_contr(2,2)*Q(39)*Q(31)+dtgup(2,2)*Q(39)+dtgup(2,3)*Q(44)+g_contr(2,3)*Q(41)*Q(34)+2*g_contr(2,1)*Q(43)*Q(30)+g_contr(2,3)*Q(31)*Q(44)+2*g_contr(2,3)*Q(46)*Q(30)+g_contr(2,3)*Q(44)*Q(35)-g_contr(2,3)*Q(49)*Q(27)-g_contr(2,3)*Q(49)*Q(31)
    dtChristoffelSRC(1,2,2) = s1+dtgup(2,3)*Q(40)-g_contr(2,3)*Q(51)*Q(30)-g_contr(2,3)*Q(35)*Q(49)+g_contr(2,3)*Q(39)*Q(32)+g_contr(2,3)*Q(40)*Q(31)+2*g_contr(2,1)*Q(42)*Q(27)+g_contr(2,1)*Q(31)*Q(42)+g_contr(2,1)*Q(34)*Q(48)+g_contr(2,2)*Q(27)*Q(39)+2*g_contr(2,2)*Q(28)*Q(37)+g_contr(2,2)*Q(33)*Q(51)+2*g_contr(2,3)*Q(28)*Q(38)+g_contr(2,3)*Q(40)*Q(35)+g_contr(2,3)*Q(42)*Q(29)+g_contr(2,3)*Q(44)*Q(27)-g_contr(2,1)*Q(17)*Q(25)*Q(7)-g_contr(2,2)*Q(17)*Q(24)*Q(10)-g_contr(2,3)*Q(48)*Q(28)+2*g_contr(2,2)*Q(40)*Q(34)+g_contr(2,3)*Q(27)*Q(40)
    s1 = dtgup(3,3)*Q(40)-g_contr(3,3)*Q(17)*Q(24)*Q(11)+g_contr(3,3)*Q(17)*Q(26)*Q(8)-g_contr(3,3)*Q(17)*Q(25)*Q(9)+dtgup(3,3)*Q(44)+2*g_contr(3,3)*Q(28)*Q(38)+g_contr(3,1)*Q(28)*Q(36)-dtgup(3,3)*Q(49)+dtgup(3,1)*Q(42)+dtgup(3,2)*Q(39)+2*g_contr(3,1)*Q(43)*Q(30)+g_contr(3,1)*Q(34)*Q(48)+2*g_contr(3,1)*Q(44)*Q(33)+2*g_contr(3,2)*Q(28)*Q(37)+g_contr(3,2)*Q(45)*Q(30)+2*g_contr(3,2)*Q(39)*Q(31)+g_contr(3,2)*Q(33)*Q(51)+2*g_contr(3,2)*Q(40)*Q(34)+g_contr(3,3)*Q(40)*Q(35)
    dtChristoffelSRC(1,2,3) = s1+g_contr(3,3)*Q(41)*Q(34)+g_contr(3,3)*Q(42)*Q(29)+g_contr(3,3)*Q(44)*Q(27)+g_contr(3,3)*Q(31)*Q(44)+2*g_contr(3,3)*Q(46)*Q(30)+g_contr(3,3)*Q(44)*Q(35)+g_contr(3,3)*Q(47)*Q(33)-g_contr(3,3)*Q(48)*Q(28)-g_contr(3,3)*Q(49)*Q(27)-g_contr(3,3)*Q(49)*Q(31)-g_contr(3,3)*Q(51)*Q(30)-g_contr(3,3)*Q(35)*Q(49)+g_contr(3,3)*Q(27)*Q(40)+g_contr(3,3)*Q(39)*Q(32)+g_contr(3,2)*Q(27)*Q(39)+g_contr(3,1)*Q(31)*Q(42)+2*g_contr(3,1)*Q(42)*Q(27)-g_contr(3,1)*Q(17)*Q(25)*Q(7)-g_contr(3,2)*Q(17)*Q(24)*Q(10)+g_contr(3,3)*Q(40)*Q(31)
    s1 = g_contr(1,2)*Q(27)*Q(40)-g_contr(1,1)*Q(17)*Q(26)*Q(7)-g_contr(1,2)*Q(17)*Q(24)*Q(11)+g_contr(1,2)*Q(17)*Q(25)*Q(9)-g_contr(1,2)*Q(17)*Q(26)*Q(8)-g_contr(1,3)*Q(17)*Q(24)*Q(12)+2*g_contr(1,1)*Q(48)*Q(27)-g_contr(1,2)*Q(47)*Q(33)+g_contr(1,1)*Q(29)*Q(36)+dtgup(1,2)*Q(40)+g_contr(1,1)*Q(32)*Q(42)-dtgup(1,2)*Q(44)+dtgup(1,1)*Q(48)+2*g_contr(1,1)*Q(49)*Q(30)+dtgup(1,2)*Q(49)+dtgup(1,3)*Q(41)+g_contr(1,1)*Q(35)*Q(48)+2*g_contr(1,1)*Q(50)*Q(33)+g_contr(1,2)*Q(40)*Q(35)
    dtChristoffelSRC(1,3,1) = s1+g_contr(1,2)*Q(41)*Q(34)-g_contr(1,2)*Q(42)*Q(29)-g_contr(1,2)*Q(44)*Q(27)-g_contr(1,2)*Q(31)*Q(44)-g_contr(1,2)*Q(44)*Q(35)+2*g_contr(1,2)*Q(29)*Q(37)+g_contr(1,2)*Q(48)*Q(28)+g_contr(1,2)*Q(49)*Q(27)+g_contr(1,2)*Q(49)*Q(31)+g_contr(1,2)*Q(51)*Q(30)+g_contr(1,2)*Q(35)*Q(49)+2*g_contr(1,2)*Q(52)*Q(33)+g_contr(1,2)*Q(39)*Q(32)+g_contr(1,2)*Q(40)*Q(31)+g_contr(1,3)*Q(27)*Q(41)+2*g_contr(1,3)*Q(29)*Q(38)+g_contr(1,3)*Q(30)*Q(47)+2*g_contr(1,3)*Q(40)*Q(32)+g_contr(1,3)*Q(53)*Q(33)+2*g_contr(1,3)*Q(41)*Q(35)
    s1 = -g_contr(2,2)*Q(17)*Q(24)*Q(11)-g_contr(2,2)*Q(44)*Q(35)-g_contr(2,2)*Q(47)*Q(33)+g_contr(2,2)*Q(48)*Q(28)+g_contr(2,2)*Q(49)*Q(27)+g_contr(2,2)*Q(49)*Q(31)+g_contr(2,2)*Q(51)*Q(30)-g_contr(2,2)*Q(31)*Q(44)-g_contr(2,1)*Q(17)*Q(26)*Q(7)+g_contr(2,1)*Q(29)*Q(36)+2*g_contr(2,2)*Q(29)*Q(37)-g_contr(2,2)*Q(44)*Q(27)+g_contr(2,3)*Q(27)*Q(41)+2*g_contr(2,2)*Q(52)*Q(33)+g_contr(2,2)*Q(27)*Q(40)+g_contr(2,2)*Q(39)*Q(32)+g_contr(2,2)*Q(40)*Q(31)+2*g_contr(2,3)*Q(29)*Q(38)+g_contr(2,3)*Q(30)*Q(47)
    dtChristoffelSRC(1,3,2) = s1+2*g_contr(2,3)*Q(40)*Q(32)+g_contr(2,3)*Q(53)*Q(33)+2*g_contr(2,3)*Q(41)*Q(35)+g_contr(2,2)*Q(35)*Q(49)-g_contr(2,2)*Q(17)*Q(26)*Q(8)+g_contr(2,2)*Q(17)*Q(25)*Q(9)+dtgup(2,2)*Q(40)+dtgup(2,3)*Q(41)+dtgup(2,2)*Q(49)+dtgup(2,1)*Q(48)-g_contr(2,3)*Q(17)*Q(24)*Q(12)+2*g_contr(2,1)*Q(48)*Q(27)+g_contr(2,1)*Q(32)*Q(42)+2*g_contr(2,1)*Q(49)*Q(30)+g_contr(2,1)*Q(35)*Q(48)+2*g_contr(2,1)*Q(50)*Q(33)+g_contr(2,2)*Q(40)*Q(35)+g_contr(2,2)*Q(41)*Q(34)-g_contr(2,2)*Q(42)*Q(29)-dtgup(2,2)*Q(44)
    s1 = g_contr(3,2)*Q(40)*Q(35)+g_contr(3,2)*Q(41)*Q(34)-g_contr(3,2)*Q(42)*Q(29)-g_contr(3,2)*Q(44)*Q(27)-g_contr(3,2)*Q(44)*Q(35)+2*g_contr(3,2)*Q(29)*Q(37)+g_contr(3,2)*Q(48)*Q(28)+g_contr(3,2)*Q(49)*Q(27)+g_contr(3,2)*Q(49)*Q(31)+g_contr(3,2)*Q(51)*Q(30)+g_contr(3,2)*Q(35)*Q(49)+2*g_contr(3,2)*Q(52)*Q(33)+g_contr(3,2)*Q(39)*Q(32)+g_contr(3,2)*Q(40)*Q(31)+g_contr(3,3)*Q(30)*Q(47)+2*g_contr(3,3)*Q(40)*Q(32)+g_contr(3,3)*Q(53)*Q(33)+g_contr(3,3)*Q(27)*Q(41)+dtgup(3,1)*Q(48)
    dtChristoffelSRC(1,3,3) = s1-g_contr(3,2)*Q(31)*Q(44)-g_contr(3,2)*Q(47)*Q(33)+2*g_contr(3,3)*Q(41)*Q(35)+g_contr(3,2)*Q(27)*Q(40)+dtgup(3,2)*Q(49)+dtgup(3,2)*Q(40)-dtgup(3,2)*Q(44)+2*g_contr(3,1)*Q(48)*Q(27)+2*g_contr(3,1)*Q(49)*Q(30)+dtgup(3,3)*Q(41)-g_contr(3,1)*Q(17)*Q(26)*Q(7)-g_contr(3,2)*Q(17)*Q(24)*Q(11)+g_contr(3,2)*Q(17)*Q(25)*Q(9)-g_contr(3,2)*Q(17)*Q(26)*Q(8)-g_contr(3,3)*Q(17)*Q(24)*Q(12)+g_contr(3,1)*Q(29)*Q(36)+g_contr(3,1)*Q(35)*Q(48)+2*g_contr(3,1)*Q(50)*Q(33)+2*g_contr(3,3)*Q(29)*Q(38)+g_contr(3,1)*Q(32)*Q(42)
    s1 = g_contr(1,3)*Q(40)*Q(35)+2*g_contr(1,1)*Q(44)*Q(33)-g_contr(1,3)*Q(35)*Q(49)+g_contr(1,3)*Q(44)*Q(35)+g_contr(1,1)*Q(34)*Q(48)+2*g_contr(1,1)*Q(43)*Q(30)+2*g_contr(1,2)*Q(28)*Q(37)+g_contr(1,3)*Q(44)*Q(27)+2*g_contr(1,3)*Q(28)*Q(38)+g_contr(1,2)*Q(27)*Q(39)+g_contr(1,2)*Q(45)*Q(30)+2*g_contr(1,2)*Q(39)*Q(31)+g_contr(1,2)*Q(33)*Q(51)+2*g_contr(1,2)*Q(40)*Q(34)+g_contr(1,3)*Q(41)*Q(34)+g_contr(1,3)*Q(42)*Q(29)+g_contr(1,3)*Q(31)*Q(44)-g_contr(1,3)*Q(48)*Q(28)-g_contr(1,3)*Q(49)*Q(27)
    dtChristoffelSRC(2,1,1) = s1-g_contr(1,3)*Q(49)*Q(31)-g_contr(1,3)*Q(51)*Q(30)+g_contr(1,3)*Q(27)*Q(40)+g_contr(1,3)*Q(39)*Q(32)+g_contr(1,3)*Q(40)*Q(31)-dtgup(1,3)*Q(49)+dtgup(1,3)*Q(40)+dtgup(1,2)*Q(39)+dtgup(1,3)*Q(44)+2*g_contr(1,1)*Q(42)*Q(27)+g_contr(1,1)*Q(28)*Q(36)+g_contr(1,1)*Q(31)*Q(42)+2*g_contr(1,3)*Q(46)*Q(30)+g_contr(1,3)*Q(47)*Q(33)-g_contr(1,1)*Q(17)*Q(25)*Q(7)-g_contr(1,2)*Q(17)*Q(24)*Q(10)-g_contr(1,3)*Q(17)*Q(24)*Q(11)-g_contr(1,3)*Q(17)*Q(25)*Q(9)+g_contr(1,3)*Q(17)*Q(26)*Q(8)+dtgup(1,1)*Q(42)
    s1 = dtgup(2,1)*Q(42)-dtgup(2,3)*Q(49)+g_contr(2,2)*Q(45)*Q(30)-g_contr(2,3)*Q(17)*Q(24)*Q(11)-g_contr(2,3)*Q(17)*Q(25)*Q(9)+g_contr(2,3)*Q(17)*Q(26)*Q(8)+2*g_contr(2,1)*Q(44)*Q(33)+g_contr(2,1)*Q(28)*Q(36)+g_contr(2,3)*Q(47)*Q(33)+2*g_contr(2,2)*Q(39)*Q(31)+dtgup(2,2)*Q(39)+dtgup(2,3)*Q(44)+g_contr(2,3)*Q(41)*Q(34)+2*g_contr(2,1)*Q(43)*Q(30)+g_contr(2,3)*Q(31)*Q(44)+2*g_contr(2,3)*Q(46)*Q(30)+g_contr(2,3)*Q(44)*Q(35)-g_contr(2,3)*Q(49)*Q(27)-g_contr(2,3)*Q(49)*Q(31)
    dtChristoffelSRC(2,1,2) = s1+dtgup(2,3)*Q(40)-g_contr(2,3)*Q(51)*Q(30)-g_contr(2,3)*Q(35)*Q(49)+g_contr(2,3)*Q(39)*Q(32)+g_contr(2,3)*Q(40)*Q(31)+2*g_contr(2,1)*Q(42)*Q(27)+g_contr(2,1)*Q(31)*Q(42)+g_contr(2,1)*Q(34)*Q(48)+g_contr(2,2)*Q(27)*Q(39)+2*g_contr(2,2)*Q(28)*Q(37)+g_contr(2,2)*Q(33)*Q(51)+2*g_contr(2,3)*Q(28)*Q(38)+g_contr(2,3)*Q(40)*Q(35)+g_contr(2,3)*Q(42)*Q(29)+g_contr(2,3)*Q(44)*Q(27)-g_contr(2,1)*Q(17)*Q(25)*Q(7)-g_contr(2,2)*Q(17)*Q(24)*Q(10)-g_contr(2,3)*Q(48)*Q(28)+2*g_contr(2,2)*Q(40)*Q(34)+g_contr(2,3)*Q(27)*Q(40)
    s1 = dtgup(3,3)*Q(40)-g_contr(3,3)*Q(17)*Q(24)*Q(11)+g_contr(3,3)*Q(17)*Q(26)*Q(8)-g_contr(3,3)*Q(17)*Q(25)*Q(9)+dtgup(3,3)*Q(44)+2*g_contr(3,3)*Q(28)*Q(38)+g_contr(3,1)*Q(28)*Q(36)-dtgup(3,3)*Q(49)+dtgup(3,1)*Q(42)+dtgup(3,2)*Q(39)+2*g_contr(3,1)*Q(43)*Q(30)+g_contr(3,1)*Q(34)*Q(48)+2*g_contr(3,1)*Q(44)*Q(33)+2*g_contr(3,2)*Q(28)*Q(37)+g_contr(3,2)*Q(45)*Q(30)+2*g_contr(3,2)*Q(39)*Q(31)+g_contr(3,2)*Q(33)*Q(51)+2*g_contr(3,2)*Q(40)*Q(34)+g_contr(3,3)*Q(40)*Q(35)
    dtChristoffelSRC(2,1,3) = s1+g_contr(3,3)*Q(41)*Q(34)+g_contr(3,3)*Q(42)*Q(29)+g_contr(3,3)*Q(44)*Q(27)+g_contr(3,3)*Q(31)*Q(44)+2*g_contr(3,3)*Q(46)*Q(30)+g_contr(3,3)*Q(44)*Q(35)+g_contr(3,3)*Q(47)*Q(33)-g_contr(3,3)*Q(48)*Q(28)-g_contr(3,3)*Q(49)*Q(27)-g_contr(3,3)*Q(49)*Q(31)-g_contr(3,3)*Q(51)*Q(30)-g_contr(3,3)*Q(35)*Q(49)+g_contr(3,3)*Q(27)*Q(40)+g_contr(3,3)*Q(39)*Q(32)+g_contr(3,2)*Q(27)*Q(39)+g_contr(3,1)*Q(31)*Q(42)+2*g_contr(3,1)*Q(42)*Q(27)-g_contr(3,1)*Q(17)*Q(25)*Q(7)-g_contr(3,2)*Q(17)*Q(24)*Q(10)+g_contr(3,3)*Q(40)*Q(31)
    s1 = 4*g_contr(1,1)*Q(31)*Q(43)+g_contr(1,1)*Q(45)*Q(30)-2*g_contr(1,1)*Q(17)*Q(25)*Q(8)+g_contr(1,1)*Q(17)*Q(24)*Q(10)-g_contr(1,2)*Q(17)*Q(25)*Q(10)+g_contr(1,3)*Q(17)*Q(26)*Q(10)-2*g_contr(1,3)*Q(17)*Q(25)*Q(11)+2*g_contr(1,3)*Q(47)*Q(34)-dtgup(1,1)*Q(39)+dtgup(1,2)*Q(45)+2*g_contr(1,1)*Q(44)*Q(34)+2*g_contr(1,1)*Q(46)*Q(33)-g_contr(1,1)*Q(27)*Q(39)-2*g_contr(1,1)*Q(39)*Q(31)-g_contr(1,1)*Q(33)*Q(51)-2*g_contr(1,1)*Q(40)*Q(34)+g_contr(1,2)*Q(28)*Q(39)+2*g_contr(1,2)*Q(43)*Q(28)
    dtChristoffelSRC(2,2,1) = s1+3*g_contr(1,2)*Q(31)*Q(45)+g_contr(1,2)*Q(34)*Q(51)+2*g_contr(1,2)*Q(46)*Q(34)+g_contr(1,3)*Q(32)*Q(45)-2*g_contr(1,3)*Q(51)*Q(31)-g_contr(1,3)*Q(35)*Q(51)-2*g_contr(1,3)*Q(49)*Q(28)+2*g_contr(1,3)*Q(28)*Q(40)+2*g_contr(1,3)*Q(43)*Q(29)+2*g_contr(1,3)*Q(44)*Q(28)+2*g_contr(1,1)*Q(42)*Q(28)+4*g_contr(1,3)*Q(31)*Q(46)-dtgup(1,3)*Q(51)+2*g_contr(1,3)*Q(46)*Q(35)-g_contr(1,3)*Q(29)*Q(39)+2*dtgup(1,1)*Q(43)+2*dtgup(1,3)*Q(46)+2*g_contr(1,1)*Q(43)*Q(27)+2*g_contr(1,1)*Q(34)*Q(49)
    s1 = 4*g_contr(2,1)*Q(31)*Q(43)+2*g_contr(2,1)*Q(43)*Q(27)+g_contr(2,1)*Q(45)*Q(30)+g_contr(2,2)*Q(34)*Q(51)+2*g_contr(2,2)*Q(46)*Q(34)+2*dtgup(2,1)*Q(43)+2*g_contr(2,3)*Q(47)*Q(34)-g_contr(2,3)*Q(29)*Q(39)+2*g_contr(2,1)*Q(44)*Q(34)+2*dtgup(2,3)*Q(46)+dtgup(2,2)*Q(45)-dtgup(2,1)*Q(39)+2*g_contr(2,1)*Q(42)*Q(28)+g_contr(2,3)*Q(32)*Q(45)-2*g_contr(2,3)*Q(51)*Q(31)-g_contr(2,3)*Q(35)*Q(51)-2*g_contr(2,3)*Q(49)*Q(28)+2*g_contr(2,3)*Q(28)*Q(40)
    dtChristoffelSRC(2,2,2) = s1+2*g_contr(2,3)*Q(43)*Q(29)+2*g_contr(2,3)*Q(44)*Q(28)+4*g_contr(2,3)*Q(31)*Q(46)+2*g_contr(2,3)*Q(46)*Q(35)+2*g_contr(2,1)*Q(34)*Q(49)+2*g_contr(2,1)*Q(46)*Q(33)-g_contr(2,1)*Q(27)*Q(39)-2*g_contr(2,1)*Q(39)*Q(31)-g_contr(2,1)*Q(33)*Q(51)-2*g_contr(2,1)*Q(40)*Q(34)+g_contr(2,2)*Q(28)*Q(39)+2*g_contr(2,2)*Q(43)*Q(28)+3*g_contr(2,2)*Q(31)*Q(45)-2*g_contr(2,1)*Q(17)*Q(25)*Q(8)+g_contr(2,1)*Q(17)*Q(24)*Q(10)-g_contr(2,2)*Q(17)*Q(25)*Q(10)+g_contr(2,3)*Q(17)*Q(26)*Q(10)-2*g_contr(2,3)*Q(17)*Q(25)*Q(11)-dtgup(2,3)*Q(51)
    s1 = 2*g_contr(3,3)*Q(47)*Q(34)+g_contr(3,2)*Q(28)*Q(39)+2*g_contr(3,1)*Q(42)*Q(28)+2*g_contr(3,1)*Q(43)*Q(27)+4*g_contr(3,1)*Q(31)*Q(43)+g_contr(3,1)*Q(45)*Q(30)+2*g_contr(3,1)*Q(34)*Q(49)+2*g_contr(3,1)*Q(44)*Q(34)+2*g_contr(3,1)*Q(46)*Q(33)-g_contr(3,1)*Q(27)*Q(39)-2*g_contr(3,1)*Q(39)*Q(31)-g_contr(3,1)*Q(33)*Q(51)-2*g_contr(3,1)*Q(40)*Q(34)+2*g_contr(3,2)*Q(43)*Q(28)+3*g_contr(3,2)*Q(31)*Q(45)+g_contr(3,2)*Q(34)*Q(51)+2*g_contr(3,2)*Q(46)*Q(34)-2*g_contr(3,3)*Q(51)*Q(31)
    dtChristoffelSRC(2,2,3) = s1-g_contr(3,3)*Q(35)*Q(51)-g_contr(3,3)*Q(29)*Q(39)-2*g_contr(3,3)*Q(49)*Q(28)+2*g_contr(3,3)*Q(28)*Q(40)+2*g_contr(3,3)*Q(43)*Q(29)+2*dtgup(3,3)*Q(46)+2*g_contr(3,3)*Q(44)*Q(28)+4*g_contr(3,3)*Q(31)*Q(46)+g_contr(3,3)*Q(17)*Q(26)*Q(10)-2*g_contr(3,3)*Q(17)*Q(25)*Q(11)-2*g_contr(3,1)*Q(17)*Q(25)*Q(8)+g_contr(3,1)*Q(17)*Q(24)*Q(10)-g_contr(3,2)*Q(17)*Q(25)*Q(10)+2*dtgup(3,1)*Q(43)-dtgup(3,1)*Q(39)+dtgup(3,2)*Q(45)-dtgup(3,3)*Q(51)+2*g_contr(3,3)*Q(46)*Q(35)+g_contr(3,3)*Q(32)*Q(45)
    s1 = dtgup(1,3)*Q(47)+dtgup(1,1)*Q(49)-dtgup(1,1)*Q(40)-g_contr(1,2)*Q(17)*Q(26)*Q(10)-g_contr(1,3)*Q(17)*Q(25)*Q(12)-g_contr(1,1)*Q(40)*Q(35)-g_contr(1,1)*Q(41)*Q(34)+2*g_contr(1,1)*Q(32)*Q(43)+g_contr(1,1)*Q(42)*Q(29)+g_contr(1,1)*Q(44)*Q(27)+g_contr(1,1)*Q(31)*Q(44)+g_contr(1,1)*Q(44)*Q(35)+g_contr(1,1)*Q(47)*Q(33)+g_contr(1,1)*Q(48)*Q(28)+g_contr(1,1)*Q(49)*Q(27)+g_contr(1,1)*Q(49)*Q(31)+g_contr(1,1)*Q(51)*Q(30)+g_contr(1,1)*Q(35)*Q(49)+2*g_contr(1,1)*Q(50)*Q(34)
    dtChristoffelSRC(2,3,1) = s1-g_contr(1,1)*Q(27)*Q(40)-g_contr(1,1)*Q(39)*Q(32)-g_contr(1,1)*Q(40)*Q(31)+g_contr(1,2)*Q(29)*Q(39)+2*g_contr(1,2)*Q(49)*Q(28)+g_contr(1,2)*Q(32)*Q(45)+2*g_contr(1,2)*Q(51)*Q(31)+g_contr(1,2)*Q(35)*Q(51)+g_contr(1,1)*Q(17)*Q(24)*Q(11)-g_contr(1,1)*Q(17)*Q(25)*Q(9)-g_contr(1,1)*Q(17)*Q(26)*Q(8)+dtgup(1,1)*Q(44)+2*g_contr(1,2)*Q(52)*Q(34)+g_contr(1,3)*Q(28)*Q(41)+2*g_contr(1,3)*Q(44)*Q(29)+g_contr(1,3)*Q(31)*Q(47)+2*g_contr(1,3)*Q(32)*Q(46)+g_contr(1,3)*Q(53)*Q(34)+2*g_contr(1,3)*Q(47)*Q(35)+dtgup(1,2)*Q(51)
    s1 = dtgup(2,1)*Q(44)-g_contr(2,1)*Q(41)*Q(34)+dtgup(2,2)*Q(51)-dtgup(2,1)*Q(40)+g_contr(2,1)*Q(47)*Q(33)+g_contr(2,1)*Q(48)*Q(28)+g_contr(2,1)*Q(49)*Q(27)+g_contr(2,3)*Q(28)*Q(41)+2*g_contr(2,1)*Q(32)*Q(43)+g_contr(2,1)*Q(42)*Q(29)+g_contr(2,1)*Q(31)*Q(44)+g_contr(2,1)*Q(44)*Q(35)+g_contr(2,1)*Q(49)*Q(31)+g_contr(2,1)*Q(35)*Q(49)-g_contr(2,1)*Q(27)*Q(40)-g_contr(2,1)*Q(39)*Q(32)-g_contr(2,1)*Q(40)*Q(31)+g_contr(2,2)*Q(29)*Q(39)+2*g_contr(2,2)*Q(49)*Q(28)
    dtChristoffelSRC(2,3,2) = s1+g_contr(2,2)*Q(32)*Q(45)+2*g_contr(2,2)*Q(51)*Q(31)+g_contr(2,2)*Q(35)*Q(51)+2*g_contr(2,2)*Q(52)*Q(34)+2*g_contr(2,3)*Q(44)*Q(29)+g_contr(2,3)*Q(31)*Q(47)+dtgup(2,1)*Q(49)+dtgup(2,3)*Q(47)+g_contr(2,1)*Q(17)*Q(24)*Q(11)-g_contr(2,1)*Q(17)*Q(25)*Q(9)-g_contr(2,1)*Q(17)*Q(26)*Q(8)-g_contr(2,2)*Q(17)*Q(26)*Q(10)-g_contr(2,3)*Q(17)*Q(25)*Q(12)+2*g_contr(2,3)*Q(32)*Q(46)+g_contr(2,3)*Q(53)*Q(34)+2*g_contr(2,3)*Q(47)*Q(35)+g_contr(2,1)*Q(51)*Q(30)-g_contr(2,1)*Q(40)*Q(35)+2*g_contr(2,1)*Q(50)*Q(34)+g_contr(2,1)*Q(44)*Q(27)
    s1 = dtgup(3,1)*Q(49)+g_contr(3,1)*Q(17)*Q(24)*Q(11)-g_contr(3,1)*Q(17)*Q(25)*Q(9)-g_contr(3,1)*Q(17)*Q(26)*Q(8)-g_contr(3,2)*Q(17)*Q(26)*Q(10)-g_contr(3,3)*Q(17)*Q(25)*Q(12)+g_contr(3,3)*Q(28)*Q(41)+dtgup(3,2)*Q(51)+dtgup(3,3)*Q(47)-dtgup(3,1)*Q(40)+g_contr(3,2)*Q(35)*Q(51)+2*g_contr(3,2)*Q(52)*Q(34)-g_contr(3,1)*Q(41)*Q(34)+2*g_contr(3,1)*Q(32)*Q(43)+g_contr(3,1)*Q(42)*Q(29)+g_contr(3,1)*Q(44)*Q(27)+g_contr(3,1)*Q(31)*Q(44)+g_contr(3,1)*Q(44)*Q(35)+g_contr(3,1)*Q(47)*Q(33)
    dtChristoffelSRC(2,3,3) = s1+g_contr(3,1)*Q(48)*Q(28)+g_contr(3,1)*Q(49)*Q(27)+g_contr(3,1)*Q(49)*Q(31)+g_contr(3,1)*Q(51)*Q(30)+g_contr(3,1)*Q(35)*Q(49)+2*g_contr(3,1)*Q(50)*Q(34)-g_contr(3,1)*Q(27)*Q(40)-g_contr(3,1)*Q(39)*Q(32)-g_contr(3,1)*Q(40)*Q(31)+g_contr(3,2)*Q(29)*Q(39)+2*g_contr(3,2)*Q(49)*Q(28)+g_contr(3,2)*Q(32)*Q(45)+2*g_contr(3,2)*Q(51)*Q(31)+2*g_contr(3,3)*Q(44)*Q(29)+g_contr(3,3)*Q(31)*Q(47)+2*g_contr(3,3)*Q(32)*Q(46)+g_contr(3,3)*Q(53)*Q(34)-g_contr(3,1)*Q(40)*Q(35)+2*g_contr(3,3)*Q(47)*Q(35)+dtgup(3,1)*Q(44)
    s1 = g_contr(1,2)*Q(27)*Q(40)-g_contr(1,1)*Q(17)*Q(26)*Q(7)-g_contr(1,2)*Q(17)*Q(24)*Q(11)+g_contr(1,2)*Q(17)*Q(25)*Q(9)-g_contr(1,2)*Q(17)*Q(26)*Q(8)-g_contr(1,3)*Q(17)*Q(24)*Q(12)+2*g_contr(1,1)*Q(48)*Q(27)-g_contr(1,2)*Q(47)*Q(33)+g_contr(1,1)*Q(29)*Q(36)+dtgup(1,2)*Q(40)+g_contr(1,1)*Q(32)*Q(42)-dtgup(1,2)*Q(44)+dtgup(1,1)*Q(48)+2*g_contr(1,1)*Q(49)*Q(30)+dtgup(1,2)*Q(49)+dtgup(1,3)*Q(41)+g_contr(1,1)*Q(35)*Q(48)+2*g_contr(1,1)*Q(50)*Q(33)+g_contr(1,2)*Q(40)*Q(35)
    dtChristoffelSRC(3,1,1) = s1+g_contr(1,2)*Q(41)*Q(34)-g_contr(1,2)*Q(42)*Q(29)-g_contr(1,2)*Q(44)*Q(27)-g_contr(1,2)*Q(31)*Q(44)-g_contr(1,2)*Q(44)*Q(35)+2*g_contr(1,2)*Q(29)*Q(37)+g_contr(1,2)*Q(48)*Q(28)+g_contr(1,2)*Q(49)*Q(27)+g_contr(1,2)*Q(49)*Q(31)+g_contr(1,2)*Q(51)*Q(30)+g_contr(1,2)*Q(35)*Q(49)+2*g_contr(1,2)*Q(52)*Q(33)+g_contr(1,2)*Q(39)*Q(32)+g_contr(1,2)*Q(40)*Q(31)+g_contr(1,3)*Q(27)*Q(41)+2*g_contr(1,3)*Q(29)*Q(38)+g_contr(1,3)*Q(30)*Q(47)+2*g_contr(1,3)*Q(40)*Q(32)+g_contr(1,3)*Q(53)*Q(33)+2*g_contr(1,3)*Q(41)*Q(35)
    s1 = -g_contr(2,2)*Q(17)*Q(24)*Q(11)-g_contr(2,2)*Q(44)*Q(35)-g_contr(2,2)*Q(47)*Q(33)+g_contr(2,2)*Q(48)*Q(28)+g_contr(2,2)*Q(49)*Q(27)+g_contr(2,2)*Q(49)*Q(31)+g_contr(2,2)*Q(51)*Q(30)-g_contr(2,2)*Q(31)*Q(44)-g_contr(2,1)*Q(17)*Q(26)*Q(7)+g_contr(2,1)*Q(29)*Q(36)+2*g_contr(2,2)*Q(29)*Q(37)-g_contr(2,2)*Q(44)*Q(27)+g_contr(2,3)*Q(27)*Q(41)+2*g_contr(2,2)*Q(52)*Q(33)+g_contr(2,2)*Q(27)*Q(40)+g_contr(2,2)*Q(39)*Q(32)+g_contr(2,2)*Q(40)*Q(31)+2*g_contr(2,3)*Q(29)*Q(38)+g_contr(2,3)*Q(30)*Q(47)
    dtChristoffelSRC(3,1,2) = s1+2*g_contr(2,3)*Q(40)*Q(32)+g_contr(2,3)*Q(53)*Q(33)+2*g_contr(2,3)*Q(41)*Q(35)+g_contr(2,2)*Q(35)*Q(49)-g_contr(2,2)*Q(17)*Q(26)*Q(8)+g_contr(2,2)*Q(17)*Q(25)*Q(9)+dtgup(2,2)*Q(40)+dtgup(2,3)*Q(41)+dtgup(2,2)*Q(49)+dtgup(2,1)*Q(48)-g_contr(2,3)*Q(17)*Q(24)*Q(12)+2*g_contr(2,1)*Q(48)*Q(27)+g_contr(2,1)*Q(32)*Q(42)+2*g_contr(2,1)*Q(49)*Q(30)+g_contr(2,1)*Q(35)*Q(48)+2*g_contr(2,1)*Q(50)*Q(33)+g_contr(2,2)*Q(40)*Q(35)+g_contr(2,2)*Q(41)*Q(34)-g_contr(2,2)*Q(42)*Q(29)-dtgup(2,2)*Q(44)
    s1 = g_contr(3,2)*Q(40)*Q(35)+g_contr(3,2)*Q(41)*Q(34)-g_contr(3,2)*Q(42)*Q(29)-g_contr(3,2)*Q(44)*Q(27)-g_contr(3,2)*Q(44)*Q(35)+2*g_contr(3,2)*Q(29)*Q(37)+g_contr(3,2)*Q(48)*Q(28)+g_contr(3,2)*Q(49)*Q(27)+g_contr(3,2)*Q(49)*Q(31)+g_contr(3,2)*Q(51)*Q(30)+g_contr(3,2)*Q(35)*Q(49)+2*g_contr(3,2)*Q(52)*Q(33)+g_contr(3,2)*Q(39)*Q(32)+g_contr(3,2)*Q(40)*Q(31)+g_contr(3,3)*Q(30)*Q(47)+2*g_contr(3,3)*Q(40)*Q(32)+g_contr(3,3)*Q(53)*Q(33)+g_contr(3,3)*Q(27)*Q(41)+dtgup(3,1)*Q(48)
    dtChristoffelSRC(3,1,3) = s1-g_contr(3,2)*Q(31)*Q(44)-g_contr(3,2)*Q(47)*Q(33)+2*g_contr(3,3)*Q(41)*Q(35)+g_contr(3,2)*Q(27)*Q(40)+dtgup(3,2)*Q(49)+dtgup(3,2)*Q(40)-dtgup(3,2)*Q(44)+2*g_contr(3,1)*Q(48)*Q(27)+2*g_contr(3,1)*Q(49)*Q(30)+dtgup(3,3)*Q(41)-g_contr(3,1)*Q(17)*Q(26)*Q(7)-g_contr(3,2)*Q(17)*Q(24)*Q(11)+g_contr(3,2)*Q(17)*Q(25)*Q(9)-g_contr(3,2)*Q(17)*Q(26)*Q(8)-g_contr(3,3)*Q(17)*Q(24)*Q(12)+g_contr(3,1)*Q(29)*Q(36)+g_contr(3,1)*Q(35)*Q(48)+2*g_contr(3,1)*Q(50)*Q(33)+2*g_contr(3,3)*Q(29)*Q(38)+g_contr(3,1)*Q(32)*Q(42)
    s1 = dtgup(1,3)*Q(47)+dtgup(1,1)*Q(49)-dtgup(1,1)*Q(40)-g_contr(1,2)*Q(17)*Q(26)*Q(10)-g_contr(1,3)*Q(17)*Q(25)*Q(12)-g_contr(1,1)*Q(40)*Q(35)-g_contr(1,1)*Q(41)*Q(34)+2*g_contr(1,1)*Q(32)*Q(43)+g_contr(1,1)*Q(42)*Q(29)+g_contr(1,1)*Q(44)*Q(27)+g_contr(1,1)*Q(31)*Q(44)+g_contr(1,1)*Q(44)*Q(35)+g_contr(1,1)*Q(47)*Q(33)+g_contr(1,1)*Q(48)*Q(28)+g_contr(1,1)*Q(49)*Q(27)+g_contr(1,1)*Q(49)*Q(31)+g_contr(1,1)*Q(51)*Q(30)+g_contr(1,1)*Q(35)*Q(49)+2*g_contr(1,1)*Q(50)*Q(34)
    dtChristoffelSRC(3,2,1) = s1-g_contr(1,1)*Q(27)*Q(40)-g_contr(1,1)*Q(39)*Q(32)-g_contr(1,1)*Q(40)*Q(31)+g_contr(1,2)*Q(29)*Q(39)+2*g_contr(1,2)*Q(49)*Q(28)+g_contr(1,2)*Q(32)*Q(45)+2*g_contr(1,2)*Q(51)*Q(31)+g_contr(1,2)*Q(35)*Q(51)+g_contr(1,1)*Q(17)*Q(24)*Q(11)-g_contr(1,1)*Q(17)*Q(25)*Q(9)-g_contr(1,1)*Q(17)*Q(26)*Q(8)+dtgup(1,1)*Q(44)+2*g_contr(1,2)*Q(52)*Q(34)+g_contr(1,3)*Q(28)*Q(41)+2*g_contr(1,3)*Q(44)*Q(29)+g_contr(1,3)*Q(31)*Q(47)+2*g_contr(1,3)*Q(32)*Q(46)+g_contr(1,3)*Q(53)*Q(34)+2*g_contr(1,3)*Q(47)*Q(35)+dtgup(1,2)*Q(51)
    s1 = dtgup(2,1)*Q(44)-g_contr(2,1)*Q(41)*Q(34)+dtgup(2,2)*Q(51)-dtgup(2,1)*Q(40)+g_contr(2,1)*Q(47)*Q(33)+g_contr(2,1)*Q(48)*Q(28)+g_contr(2,1)*Q(49)*Q(27)+g_contr(2,3)*Q(28)*Q(41)+2*g_contr(2,1)*Q(32)*Q(43)+g_contr(2,1)*Q(42)*Q(29)+g_contr(2,1)*Q(31)*Q(44)+g_contr(2,1)*Q(44)*Q(35)+g_contr(2,1)*Q(49)*Q(31)+g_contr(2,1)*Q(35)*Q(49)-g_contr(2,1)*Q(27)*Q(40)-g_contr(2,1)*Q(39)*Q(32)-g_contr(2,1)*Q(40)*Q(31)+g_contr(2,2)*Q(29)*Q(39)+2*g_contr(2,2)*Q(49)*Q(28)
    dtChristoffelSRC(3,2,2) = s1+g_contr(2,2)*Q(32)*Q(45)+2*g_contr(2,2)*Q(51)*Q(31)+g_contr(2,2)*Q(35)*Q(51)+2*g_contr(2,2)*Q(52)*Q(34)+2*g_contr(2,3)*Q(44)*Q(29)+g_contr(2,3)*Q(31)*Q(47)+dtgup(2,1)*Q(49)+dtgup(2,3)*Q(47)+g_contr(2,1)*Q(17)*Q(24)*Q(11)-g_contr(2,1)*Q(17)*Q(25)*Q(9)-g_contr(2,1)*Q(17)*Q(26)*Q(8)-g_contr(2,2)*Q(17)*Q(26)*Q(10)-g_contr(2,3)*Q(17)*Q(25)*Q(12)+2*g_contr(2,3)*Q(32)*Q(46)+g_contr(2,3)*Q(53)*Q(34)+2*g_contr(2,3)*Q(47)*Q(35)+g_contr(2,1)*Q(51)*Q(30)-g_contr(2,1)*Q(40)*Q(35)+2*g_contr(2,1)*Q(50)*Q(34)+g_contr(2,1)*Q(44)*Q(27)
    s1 = dtgup(3,1)*Q(49)+g_contr(3,1)*Q(17)*Q(24)*Q(11)-g_contr(3,1)*Q(17)*Q(25)*Q(9)-g_contr(3,1)*Q(17)*Q(26)*Q(8)-g_contr(3,2)*Q(17)*Q(26)*Q(10)-g_contr(3,3)*Q(17)*Q(25)*Q(12)+g_contr(3,3)*Q(28)*Q(41)+dtgup(3,2)*Q(51)+dtgup(3,3)*Q(47)-dtgup(3,1)*Q(40)+g_contr(3,2)*Q(35)*Q(51)+2*g_contr(3,2)*Q(52)*Q(34)-g_contr(3,1)*Q(41)*Q(34)+2*g_contr(3,1)*Q(32)*Q(43)+g_contr(3,1)*Q(42)*Q(29)+g_contr(3,1)*Q(44)*Q(27)+g_contr(3,1)*Q(31)*Q(44)+g_contr(3,1)*Q(44)*Q(35)+g_contr(3,1)*Q(47)*Q(33)
    dtChristoffelSRC(3,2,3) = s1+g_contr(3,1)*Q(48)*Q(28)+g_contr(3,1)*Q(49)*Q(27)+g_contr(3,1)*Q(49)*Q(31)+g_contr(3,1)*Q(51)*Q(30)+g_contr(3,1)*Q(35)*Q(49)+2*g_contr(3,1)*Q(50)*Q(34)-g_contr(3,1)*Q(27)*Q(40)-g_contr(3,1)*Q(39)*Q(32)-g_contr(3,1)*Q(40)*Q(31)+g_contr(3,2)*Q(29)*Q(39)+2*g_contr(3,2)*Q(49)*Q(28)+g_contr(3,2)*Q(32)*Q(45)+2*g_contr(3,2)*Q(51)*Q(31)+2*g_contr(3,3)*Q(44)*Q(29)+g_contr(3,3)*Q(31)*Q(47)+2*g_contr(3,3)*Q(32)*Q(46)+g_contr(3,3)*Q(53)*Q(34)-g_contr(3,1)*Q(40)*Q(35)+2*g_contr(3,3)*Q(47)*Q(35)+dtgup(3,1)*Q(44)
    s1 = -dtgup(1,1)*Q(41)-2*g_contr(1,2)*Q(17)*Q(26)*Q(11)-g_contr(1,3)*Q(17)*Q(26)*Q(12)-g_contr(1,1)*Q(27)*Q(41)+2*g_contr(1,2)*Q(51)*Q(32)+dtgup(1,3)*Q(53)+g_contr(1,1)*Q(17)*Q(24)*Q(12)-2*g_contr(1,1)*Q(17)*Q(26)*Q(9)+g_contr(1,2)*Q(17)*Q(25)*Q(12)-dtgup(1,2)*Q(47)+2*g_contr(1,1)*Q(48)*Q(29)-g_contr(1,2)*Q(31)*Q(47)+2*dtgup(1,1)*Q(50)-g_contr(1,1)*Q(30)*Q(47)-2*g_contr(1,1)*Q(40)*Q(32)-2*g_contr(1,1)*Q(41)*Q(35)+2*g_contr(1,1)*Q(50)*Q(27)+2*g_contr(1,1)*Q(32)*Q(44)
    dtChristoffelSRC(3,3,1) = s1+2*g_contr(1,1)*Q(49)*Q(32)+2*g_contr(1,1)*Q(52)*Q(30)+4*g_contr(1,1)*Q(35)*Q(50)+g_contr(1,1)*Q(53)*Q(33)-g_contr(1,2)*Q(28)*Q(41)-2*g_contr(1,2)*Q(44)*Q(29)+2*g_contr(1,2)*Q(29)*Q(40)+2*g_contr(1,2)*Q(49)*Q(29)+2*g_contr(1,2)*Q(50)*Q(28)+2*g_contr(1,2)*Q(52)*Q(31)-2*g_contr(1,2)*Q(47)*Q(35)+4*g_contr(1,2)*Q(35)*Q(52)+g_contr(1,2)*Q(53)*Q(34)+g_contr(1,3)*Q(29)*Q(41)+2*g_contr(1,3)*Q(50)*Q(29)+g_contr(1,3)*Q(32)*Q(47)+2*g_contr(1,3)*Q(52)*Q(32)+3*g_contr(1,3)*Q(35)*Q(53)+2*dtgup(1,2)*Q(52)
    s1 = -dtgup(2,2)*Q(47)+2*g_contr(2,3)*Q(52)*Q(32)+2*g_contr(2,2)*Q(29)*Q(40)+4*g_contr(2,1)*Q(35)*Q(50)+2*g_contr(2,2)*Q(49)*Q(29)-g_contr(2,1)*Q(27)*Q(41)+2*g_contr(2,1)*Q(32)*Q(44)-g_contr(2,1)*Q(30)*Q(47)-2*g_contr(2,1)*Q(40)*Q(32)-2*g_contr(2,1)*Q(41)*Q(35)+2*dtgup(2,2)*Q(52)+2*dtgup(2,1)*Q(50)+g_contr(2,1)*Q(17)*Q(24)*Q(12)+2*g_contr(2,1)*Q(48)*Q(29)+2*g_contr(2,1)*Q(52)*Q(30)+g_contr(2,1)*Q(53)*Q(33)-g_contr(2,2)*Q(31)*Q(47)-g_contr(2,2)*Q(28)*Q(41)
    dtChristoffelSRC(3,3,2) = s1-2*g_contr(2,2)*Q(44)*Q(29)+2*g_contr(2,2)*Q(50)*Q(28)+2*g_contr(2,2)*Q(51)*Q(32)+2*g_contr(2,2)*Q(52)*Q(31)-2*g_contr(2,2)*Q(47)*Q(35)+4*g_contr(2,2)*Q(35)*Q(52)+g_contr(2,2)*Q(53)*Q(34)+g_contr(2,3)*Q(29)*Q(41)+2*g_contr(2,3)*Q(50)*Q(29)+dtgup(2,3)*Q(53)-2*g_contr(2,1)*Q(17)*Q(26)*Q(9)+g_contr(2,2)*Q(17)*Q(25)*Q(12)-2*g_contr(2,2)*Q(17)*Q(26)*Q(11)-g_contr(2,3)*Q(17)*Q(26)*Q(12)+3*g_contr(2,3)*Q(35)*Q(53)+2*g_contr(2,1)*Q(50)*Q(27)+2*g_contr(2,1)*Q(49)*Q(32)+g_contr(2,3)*Q(32)*Q(47)-dtgup(2,1)*Q(41)
    s1 = 2*g_contr(3,2)*Q(52)*Q(31)-g_contr(3,1)*Q(27)*Q(41)-2*g_contr(3,1)*Q(40)*Q(32)-2*g_contr(3,1)*Q(41)*Q(35)+2*g_contr(3,1)*Q(48)*Q(29)+2*g_contr(3,1)*Q(50)*Q(27)+2*g_contr(3,1)*Q(32)*Q(44)-g_contr(3,2)*Q(31)*Q(47)-g_contr(3,2)*Q(28)*Q(41)-2*g_contr(3,2)*Q(44)*Q(29)+2*g_contr(3,2)*Q(29)*Q(40)+2*g_contr(3,2)*Q(49)*Q(29)+2*g_contr(3,2)*Q(50)*Q(28)+2*g_contr(3,2)*Q(51)*Q(32)-2*g_contr(3,2)*Q(47)*Q(35)+4*g_contr(3,2)*Q(35)*Q(52)+g_contr(3,3)*Q(29)*Q(41)+2*g_contr(3,3)*Q(50)*Q(29)
    dtChristoffelSRC(3,3,3) = s1+g_contr(3,3)*Q(32)*Q(47)+2*g_contr(3,3)*Q(52)*Q(32)+3*g_contr(3,3)*Q(35)*Q(53)-dtgup(3,1)*Q(41)-g_contr(3,1)*Q(30)*Q(47)+g_contr(3,2)*Q(53)*Q(34)+2*dtgup(3,2)*Q(52)+2*g_contr(3,1)*Q(49)*Q(32)+g_contr(3,2)*Q(17)*Q(25)*Q(12)-2*g_contr(3,2)*Q(17)*Q(26)*Q(11)-g_contr(3,3)*Q(17)*Q(26)*Q(12)+g_contr(3,1)*Q(17)*Q(24)*Q(12)-2*g_contr(3,1)*Q(17)*Q(26)*Q(9)+dtgup(3,3)*Q(53)+2*g_contr(3,1)*Q(52)*Q(30)+4*g_contr(3,1)*Q(35)*Q(50)+g_contr(3,1)*Q(53)*Q(33)-dtgup(3,2)*Q(47)+2*dtgup(3,1)*Q(50)

    s2 = 2*Q(17)*g_contr(1,2)*Christoffel(1,1,1)*Q(8)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,1)*Q(7)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,3)*Q(9)+2*Q(17)*g_contr(2,3)*Christoffel(2,1,1)*Q(9)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,1)*Q(7)+2*Q(17)*g_contr(1,3)*Christoffel(1,1,3)*Q(12)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,3)*Q(9)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,1)*Q(7)+2*Q(17)*g_contr(1,3)*Christoffel(1,1,2)*Q(11)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,1)*Q(7)+dtgup(3,3)*Christoffel(3,3,1)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,3)*Q(9)+g_contr(1,1)*dtChristoffelSRC(1,1,1)
    s1 = s2+dtgup(1,2)*Christoffel(1,2,1)+g_contr(3,3)*dtChristoffelSRC(3,3,1)+g_contr(1,3)*dtChristoffelSRC(1,3,1)+g_contr(3,1)*dtChristoffelSRC(3,1,1)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,1)*Q(7)+2*Q(17)*g_contr(2,2)*Christoffel(2,1,1)*Q(8)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,3)*Q(9)+2*Q(17)*g_contr(2,3)*Christoffel(2,1,3)*Q(12)+2*Q(17)*g_contr(3,2)*Christoffel(3,1,2)*Q(10)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,2)*Q(8)+2*Q(17)*g_contr(3,2)*Christoffel(3,1,3)*Q(11)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,3)*Q(9)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,1)*Q(7)+2*Q(17)*g_contr(3,3)*Christoffel(3,1,1)*Q(9)
    s2 = 2*Q(17)*g_contr(3,3)*Christoffel(3,1,3)*Q(12)+2*Q(17)*g_contr(3,3)*Christoffel(3,1,2)*Q(11)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,2)*Q(8)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,2)*Q(8)+2*Q(17)*g_contr(3,2)*Christoffel(3,1,1)*Q(8)+2*Q(17)*g_contr(1,2)*Christoffel(1,1,3)*Q(11)+dtgup(1,1)*Christoffel(1,1,1)+g_contr(2,2)*dtChristoffelSRC(2,2,1)+dtgup(2,1)*Christoffel(2,1,1)+dtgup(2,2)*Christoffel(2,2,1)+g_contr(2,1)*dtChristoffelSRC(2,1,1)+dtgup(1,3)*Christoffel(1,3,1)+s1+dtgup(3,2)*Christoffel(3,2,1)
    dtGammaupSrc(1) = s2+dtgup(3,1)*Christoffel(3,1,1)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,2)*Q(8)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,2)*Q(8)+2*Q(17)*g_contr(1,2)*Christoffel(1,1,2)*Q(10)+2*Q(17)*g_contr(2,2)*Christoffel(2,1,2)*Q(10)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,3)*Q(9)+2*Q(17)*g_contr(2,2)*Christoffel(2,1,3)*Q(11)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,2)*Q(8)+2*Q(17)*g_contr(2,3)*Christoffel(2,1,2)*Q(11)+2*Q(17)*g_contr(1,3)*Christoffel(1,1,1)*Q(9)+dtgup(2,3)*Christoffel(2,3,1)+g_contr(1,2)*dtChristoffelSRC(1,2,1)+g_contr(3,2)*dtChristoffelSRC(3,2,1)+g_contr(2,3)*dtChristoffelSRC(2,3,1)
    s2 = g_contr(1,3)*dtChristoffelSRC(1,3,2)+dtgup(1,3)*Christoffel(1,3,2)+2*Q(17)*g_contr(1,1)*Christoffel(1,2,3)*Q(9)+2*Q(17)*g_contr(1,1)*Christoffel(1,2,2)*Q(8)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,2)*Q(10)+2*Q(17)*g_contr(1,1)*Christoffel(1,2,1)*Q(7)-2*Q(17)*g_contr(2,1)*Christoffel(2,1,1)*Q(8)+2*Q(17)*g_contr(2,1)*Christoffel(2,2,1)*Q(7)+2*Q(17)*g_contr(1,3)*Christoffel(1,2,3)*Q(12)+2*Q(17)*g_contr(3,3)*Christoffel(3,2,2)*Q(11)+dtgup(2,2)*Christoffel(2,2,2)+2*Q(17)*g_contr(3,3)*Christoffel(3,2,3)*Q(12)+dtgup(3,2)*Christoffel(3,2,2)
    s1 = s2+dtgup(2,1)*Christoffel(2,1,2)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,3)*Q(11)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,2)*Q(10)+2*Q(17)*g_contr(1,3)*Christoffel(1,2,2)*Q(11)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,2)*Q(10)-2*Q(17)*g_contr(1,3)*Christoffel(1,3,1)*Q(8)+2*Q(17)*g_contr(1,3)*Christoffel(1,2,1)*Q(9)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,1)*Q(8)+2*Q(17)*g_contr(2,1)*Christoffel(2,2,2)*Q(8)-2*Q(17)*g_contr(2,1)*Christoffel(2,1,2)*Q(10)+2*Q(17)*g_contr(3,3)*Christoffel(3,2,1)*Q(9)+dtgup(1,1)*Christoffel(1,1,2)-2*Q(17)*g_contr(3,1)*Christoffel(3,1,2)*Q(10)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,3)*Q(11)
    s2 = 2*Q(17)*g_contr(2,3)*Christoffel(2,2,3)*Q(12)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,1)*Q(8)+2*Q(17)*g_contr(2,1)*Christoffel(2,2,3)*Q(9)-2*Q(17)*g_contr(2,1)*Christoffel(2,1,3)*Q(11)+2*Q(17)*g_contr(2,3)*Christoffel(2,2,1)*Q(9)+2*Q(17)*g_contr(2,3)*Christoffel(2,2,2)*Q(11)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,2)*Q(10)-2*Q(17)*g_contr(2,3)*Christoffel(2,3,1)*Q(8)-2*Q(17)*g_contr(3,3)*Christoffel(3,3,3)*Q(11)+2*Q(17)*g_contr(3,1)*Christoffel(3,2,2)*Q(8)+2*Q(17)*g_contr(3,1)*Christoffel(3,2,1)*Q(7)-2*Q(17)*g_contr(3,1)*Christoffel(3,1,1)*Q(8)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,3)*Q(11)+2*Q(17)*g_contr(3,1)*Christoffel(3,2,3)*Q(9)
    dtGammaupSrc(2) = s2-2*Q(17)*g_contr(3,1)*Christoffel(3,1,3)*Q(11)+g_contr(2,1)*dtChristoffelSRC(2,1,2)+g_contr(2,2)*dtChristoffelSRC(2,2,2)+g_contr(3,2)*dtChristoffelSRC(3,2,2)+g_contr(3,3)*dtChristoffelSRC(3,3,2)+g_contr(1,1)*dtChristoffelSRC(1,1,2)+dtgup(1,2)*Christoffel(1,2,2)+g_contr(3,1)*dtChristoffelSRC(3,1,2)+dtgup(3,1)*Christoffel(3,1,2)+dtgup(2,3)*Christoffel(2,3,2)+s1+dtgup(3,3)*Christoffel(3,3,2)+g_contr(1,2)*dtChristoffelSRC(1,2,2)+g_contr(2,3)*dtChristoffelSRC(2,3,2)
    s2 = dtgup(1,1)*Christoffel(1,1,3)+g_contr(3,2)*dtChristoffelSRC(3,2,3)+dtgup(2,1)*Christoffel(2,1,3)+dtgup(2,2)*Christoffel(2,2,3)+g_contr(1,3)*dtChristoffelSRC(1,3,3)+dtgup(1,3)*Christoffel(1,3,3)+dtgup(2,3)*Christoffel(2,3,3)+g_contr(3,1)*dtChristoffelSRC(3,1,3)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,1)*Q(9)+2*Q(17)*g_contr(3,1)*Christoffel(3,3,1)*Q(7)-2*Q(17)*g_contr(3,1)*Christoffel(3,1,1)*Q(9)+2*Q(17)*g_contr(3,2)*Christoffel(3,3,1)*Q(8)+2*Q(17)*g_contr(2,1)*Christoffel(2,3,3)*Q(9)
    s1 = s2+2*Q(17)*g_contr(3,2)*Christoffel(3,3,3)*Q(11)+2*Q(17)*g_contr(1,1)*Christoffel(1,3,3)*Q(9)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,3)*Q(12)+2*Q(17)*g_contr(3,1)*Christoffel(3,3,2)*Q(8)-2*Q(17)*g_contr(3,1)*Christoffel(3,1,3)*Q(12)+2*Q(17)*g_contr(3,1)*Christoffel(3,3,3)*Q(9)+2*Q(17)*g_contr(1,1)*Christoffel(1,3,2)*Q(8)+2*Q(17)*g_contr(1,1)*Christoffel(1,3,1)*Q(7)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,2)*Q(11)+dtgup(1,2)*Christoffel(1,2,3)+dtgup(3,2)*Christoffel(3,2,3)+g_contr(3,3)*dtChristoffelSRC(3,3,3)-2*Q(17)*g_contr(2,1)*Christoffel(2,1,2)*Q(11)+2*Q(17)*g_contr(1,2)*Christoffel(1,3,1)*Q(8)
    s2 = -2*Q(17)*g_contr(2,1)*Christoffel(2,1,1)*Q(9)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,1)*Q(9)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,1)*Q(9)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,3)*Q(12)-2*Q(17)*g_contr(2,2)*Christoffel(2,2,2)*Q(11)+2*Q(17)*g_contr(2,2)*Christoffel(2,3,2)*Q(10)-2*Q(17)*g_contr(3,2)*Christoffel(3,2,3)*Q(12)+2*Q(17)*g_contr(2,1)*Christoffel(2,3,1)*Q(7)-2*Q(17)*g_contr(2,1)*Christoffel(2,1,3)*Q(12)+2*Q(17)*g_contr(2,2)*Christoffel(2,3,3)*Q(11)-2*Q(17)*g_contr(1,1)*Christoffel(1,1,2)*Q(11)+dtgup(3,3)*Christoffel(3,3,3)+g_contr(1,2)*dtChristoffelSRC(1,2,3)+g_contr(2,2)*dtChristoffelSRC(2,2,3)
    dtGammaupSrc(3) = s2+g_contr(1,1)*dtChristoffelSRC(1,1,3)+2*Q(17)*g_contr(3,2)*Christoffel(3,3,2)*Q(10)+2*Q(17)*g_contr(2,1)*Christoffel(2,3,2)*Q(8)+2*Q(17)*g_contr(2,2)*Christoffel(2,3,1)*Q(8)+2*Q(17)*g_contr(1,2)*Christoffel(1,3,3)*Q(11)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,3)*Q(12)-2*Q(17)*g_contr(3,1)*Christoffel(3,1,2)*Q(11)+2*Q(17)*g_contr(1,2)*Christoffel(1,3,2)*Q(10)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,2)*Q(11)-2*Q(17)*g_contr(1,2)*Christoffel(1,2,1)*Q(9)+g_contr(2,1)*dtChristoffelSRC(2,1,3)+s1+dtgup(3,1)*Christoffel(3,1,3)+g_contr(2,3)*dtChristoffelSRC(2,3,3)      
    !  
    Kex(1,1) = Q(7) 
    Kex(1,2) = Q(8) 
    Kex(1,3) = Q(9) 
    Kex(2,1) = Q(8) 
    Kex(2,2) = Q(10) 
    Kex(2,3) = Q(11) 
    Kex(3,1) = Q(9) 
    Kex(3,2) = Q(11) 
    Kex(3,3) = Q(12) 
    !
    Kmix = MATMUL(g_contr, Kex)
    Kup  = MATMUL(g_contr, TRANSPOSE(Kmix))     
    !
    Z = (/ Q(13), Q(14), Q(15) /)
    Zup = MATMUL(g_contr, Z) 
    Theta = Q(16)
    !
    alpha = Q(17) 
    AA    = (/ Q(24), Q(25), Q(26) /) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    DO j = 1, 3 
     DO i = 1, 3 
        dtgamma(i,j) = SUM(beta(:)*2*DD(:,i,j))  
     ENDDO
    ENDDO     
    dtgamma = dtgamma - 2*alpha*Kex  
    DO j = 1, 3 
     DO i = 1, 3 
      DO k = 1, 3 
         dtgamma(i,j) = dtgamma(i,j) + g_cov(k,j)*BB(i,k)+g_cov(k,i)*BB(j,k) 
      ENDDO
     ENDDO
    ENDDO 
    !
    DO j = 1, 3 
     DO i = 1, 3 
      nablaijalphaSrc(i,j) = alpha*AA(j)*AA(i) 
      DO k = 1, 3 
         nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j) - Christoffel(i,j,k)*alpha*AA(k) 
      ENDDO
     ENDDO 
    ENDDO     
    !
    !RicciSrc = 0.0
    !DO n = 1, 3 
    ! DO m = 1, 3
    !  DO l = 1, 3
    !     RicciSrc(m,n) = RicciSrc(m,n) + RiemannSrc(m,l,n,l) 
    !  ENDDO
    ! ENDDO
    !ENDDO 
    !
    nablaZSrc = 0.0 
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3 
         nablaZSrc(i,j) = nablaZSrc(i,j) - Christoffel(i,j,k)*Z(k)  
      ENDDO
     ENDDO
    ENDDO 
    !
    RicciPlusNablaZSrc = RicciSrc + 1.0/ds**2*(nablaZSrc+TRANSPOSE(nablaZSrc)) 
    RPlusTwoNablaZSrc  = SUM( g_contr*RicciPlusNablaZSrc )   
    !
    traceK  = SUM(g_contr*Kex)  
    Kupdown = SUM(Kup*Kex) 
    dtK = -nablaijalphaSrc + alpha*RicciPlusNablaZSrc - 2*alpha*MATMUL(TRANSPOSE(Kmix),Kex)+alpha*(traceK-2*Theta)*Kex - alpha*k1*(1+k2)*Theta*g_cov 
    DO j = 1, 3 
     DO i = 1, 3 
      DO k = 1, 3 
         dtK(i,j) = dtK(i,j) + Kex(k,j)*BB(i,k) + Kex(k,i)*BB(j,k) 
      ENDDO
     ENDDO
    ENDDO 
    !
    dtTheta = 0.5*alpha*( e**2*RPlusTwoNablaZSrc + (traceK - 2*Theta)*traceK - Kupdown - 2.0*SUM(Zup*AA) - 2*k1*(2+k2)*Theta ) 
    !
    dKtempSrc = 0.0 
    DO j = 1, 3 
     DO i = 1, 3 
      DO l = 1, 3          
       DO m = 1, 3 
          dKtempSrc(i) = dKtempSrc(i) + g_contr(j,l)*(-Christoffel(j,l,m)*Kex(m,i)+Christoffel(j,i,m)*Kex(m,l))   
       ENDDO
      ENDDO
     ENDDO
    ENDDO 
    dtraceKSrc = 0.0 
    DO k = 1, 3 
     DO j = 1, 3
      DO i = 1, 3          
        dtraceKSrc(k) =  dtraceKSrc(k) + dgup(k,i,j)*Kex(i,j)
      ENDDO
     ENDDO
    ENDDO 
    !
    dtZ = alpha*( e**2*dKtempSrc -2*MATMUL(Kmix,Z) - Theta*AA - k1*Z ) + MATMUL(BB, Z) 
    !
    dtalpha = + SUM(beta*alpha*AA) - alpha**2*fa*(traceK-K0-2*Theta) 
    !
    QG = ggg*MATMUL( g_contr, AA - Z ) 
    DO j = 1, 3 
     DO i = 1, 3 
        DO n = 1, 3
         DO l = 1, 3 
           QG(i) = QG(i) + ggg*g_contr(i,j)*(-2*g_contr(n,l)*DD(l,n,j)) 
         ENDDO
        ENDDO
     ENDDO
    ENDDO 
    !            
    dtbeta  = sk*( MATMUL(TRANSPOSE(BB),beta) + fff*b - alpha**2*QG ) 
    !
    dtA = -alpha*fa*dtraceKSrc - alpha*AA*(fa+alpha*faa)*(traceK-K0-2*Theta) + MATMUL(BB, AA)      
    !
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3 
        dtD(k,i,j) = -alpha*AA(k)*Kex(i,j)
        DO l = 1, 3 
          dtD(k,i,j) = dtD(k,i,j) + BB(k,l)*DD(l,i,j) + DD(k,l,i)*BB(j,l) + DD(k,l,j)*BB(i,l) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO     
    !
    dtB = 0.0 
    !
    DO k = 1, 3 
     DO i = 1, 3 
       dtB(k,i) = dtB(k,i) - 2*(alpha**2*AA(k))*QG(i) 
       DO j = 1, 3 
         dtB(k,i) = dtB(k,i) - ggg*alpha**2*dgup(k,i,j)*(AA(j)-Z(j))    
         DO n = 1, 3
          DO l = 1, 3            
              dtB(k,i) = dtB(k,i) - ggg*alpha**2*g_contr(i,j)*( -2*dgup(k,n,l)*DD(l,n,j) ) - ggg*alpha**2*dgup(k,i,j)*( -2*g_contr(n,l)*DD(l,n,j) )  
          ENDDO
         ENDDO 
       ENDDO 
     ENDDO 
    ENDDO 
    dtB = dtB*sk 
    !
    S(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)     ! gamma_ij 
    S(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                             ! K_ij  
    S(13:15)  = dtZ(1:3)                                                                                     ! Z_i 
    S(16)     = dtTheta                                                                                      ! Theta 
    S(17)     = dtalpha                                                                                      ! alpha 
    S(18:20)  = dtbeta                                                                                       ! beta^i 
    S(21:23)  = 0.0                                                                                          ! b^i 
    S(24:26)  = dtA(1:3)                                                                                     ! A_k            
    S(27:35)  = 0.0                                                                                          ! B_k^i 
    S(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                 ! D_kij    
    S(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /) 
    S(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)      
    !
    CONTINUE
    !
    !INCLUDE 'Z4_SRC.inc' 
    !
    !CONTINUE 
    !
#ifdef Z4GRMHD 
    QGRMHD(1:9)   = Q(55:63)    ! hydro variables 
    QGRMHD(10)    = Q(17)       ! lapse 
    QGRMHD(11:13) = Q(18:20)    ! shift 
    QGRMHD(14:19) = Q(1:6)      ! metric 
    !
    gradQGRMHD(1:9,:)   = 0.0   ! hydro variables (not needed) 
    gradQGRMHD(10,:)    = Q(17)*(/ Q(24), Q(25), Q(26) /)       ! lapse 
    gradQGRMHD(11,:)    = (/ Q(27), Q(28), Q(29) /)             ! shift1 
    gradQGRMHD(12,:)    = (/ Q(30), Q(31), Q(32) /)             ! shift2 
    gradQGRMHD(13,:)    = (/ Q(33), Q(34), Q(35) /)             ! shift3 
    gradQGRMHD(14:19,1)  = 2*Q(36:41)                           ! metric 
    gradQGRMHD(14:19,2)  = 2*Q(42:47)                           ! metric 
    gradQGRMHD(14:19,3)  = 2*Q(48:53)                           ! metric 
    !
    CALL PDENCPGRMHD(BgradQGRMHD,QGRMHD,gradQGRMHD,par) 
    S(55:63) = -BgradQGRMHD(1:9)
#endif      
    !
#endif 
    !
#ifdef SRMHD
    ! S(9) = -EQN%DivCleaning_a*Q(9) 
#endif 
    !
#ifdef GPR3D 
    !    
    S = 0. 
    CALL PDECons2Prim(V,Q,iErr) 
    AM(1,:) = (/ V( 6), V( 7), V( 8) /) 
    AM(2,:) = (/ V( 9), V(10), V(11) /)
    AM(3,:) = (/ V(12), V(13), V(14) /)         
    G      = MATMUL( TRANSPOSE(AM), AM ) 
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
    detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2)     ! this is the determinant we have     
    detA2  = Q(1)/EQN%rho0                                ! this is the determinant we should have from the compatibility relation 
    psiM   = 3./(detA)*MATMUL(AM,devG)                    ! in the coefficient of the source term, we use the detA2 from the compatibility relation, to reduce nonlinearities in A 
    temp2  = detA2**(1./3.) 
    S(6:14) = -(/ psiM(1,1), psiM(1,2), psiM(1,3), & 
                  psiM(2,1), psiM(2,2), psiM(2,3), & 
                  psiM(3,1), psiM(3,2), psiM(3,3) /) /EQN%tau1*detA**(8./3.)   &      ! real relaxation term 
              -(detA-detA2)/EQN%tau1*Q(6:14)                                          ! artificial relaxation term to correct the errors of the ODE integrator 
    T = V(5)/V(1)/EQN%cv/(EQN%gamma-1)
    S(15:17) = -Q(15:17)/EQN%tau2 * T/V(1) ! the source term for the heat flux is very nice and simple    
    !
    ! S = 0.0 
    !
#endif 
    !
#ifdef GRMHD
    !
    !S = - 0.1*Q 
    !RETURN 
    !  
#endif
    !
END SUBROUTINE PDESource
!
! Nonconservative part of the PDE ( B(Q) * gradQ ) 
!    
SUBROUTINE PDENCP(BgradQ,Q,gradQ,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d), par(nParam)  
    REAL, INTENT(OUT) :: BgradQ(nVar) 
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,ip,iq,ii,jj,kk,ll,iErr,count    
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    REAL :: k1,k2,k3,fff,ggg,e,c,ds,xi,sk,sknl,bs,g_cov(3,3),g_contr(3,3),dgup(3,3,3)
    REAL :: det, alpha, fa, k0, dk0(3), beta0(3), b0(3), u(3), ialpha 
    REAL :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, eta, itau    
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar)  
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL :: v2,vf(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim 
    REAL :: QGRMHD(19), gradQGRMHD(19,d), BgradQGRMHD(19), ov(3) 
    REAL :: Christoffel(3,3,3), RiemannNCP(3,3,3,3), RiemannSrc(3,3,3,3), dChristoffelNCP(3,3,3,3), dChristoffel_tildeNCP(3,3,3,3), dChristoffelSrc(3,3,3,3), DD(3,3,3), dDD(3,3,3,3)  
    REAL :: AA(3), dAA(3,3), BB(3,3), dBB(3,3,3), beta(3), Kex(3,3), Kmix(3,3), Kup(3,3), Z(3), dZ(3,3), nablaZNCP(3,3), nablaZSrc(3,3), RplusNablaZNCP, RplusNablaZSrc  
    REAL :: Theta, dTheta(3), nablaijalphaNCP(3,3), nablaijalphaSrc(3,3), Ricci(3,3), RicciNCP(3,3), RicciSrc(3,3), dtraceK(3), dtraceKNCP(3), dKtempNCP(3), dZNCP(3,3), dZSrc(3,3)  
    REAL :: dtgamma(3,3), dtK(3,3), dK(3,3,3), dtTheta, dtZ(3), dtalpha, dtGhat(3), dtbeta(3), dtbb(3), dxbb(3,3), dtA(3), dtB(3,3), dtD(3,3,3)  
    REAL :: Aupdown, Aex(3,3), dAex(3,3,3), Amix(3,3), Aup(3,3), Ghat(3), Gtilde(3), dGhat(3,3), traceK, Kupdown, phi, dphi(3), PP(3), dPP(3,3), TwoNablaZNCP, TwoNablaZSrc, dKex(3,3,3)      
    REAL :: dGtildeSrc(3,3), dGtildeNCP(3,3), RiccitildeNCP(3,3), RicciphiNCP(3,3), RiccitildeSrc(3,3), RicciphiSrc(3,3), Mom(3), Ham, Pup(3), DDcontr(3)   
    REAL :: Christoffel_tilde(3,3,3), Christoffel_kind1(3,3,3), Zup(3), RicciPlusNablaZNCP(3,3), RicciPlusNablaZSrc(3,3), traceA, traceB, QG(3), b(3), faa, temp   
    REAL :: SecondOrderTermsNCP(3,3), SecondOrderTermsSrc(3,3), traceNCP, traceSrc, dtphi, dtTraceK, dtP(3), dtX(3), XX(3), dXX(3,3), nablaXNCP(3,3)     
    REAL :: RPlusTwoNablaZNCP, RNCP, RSrc, RPlusTwoNablaZSrc, nablanablaalpha, nablanablaalphaNCP, nablanablaalphaSrc, Riemann(3,3,3,3), dChristoffel(3,3,3,3), divAupNCP(3), divAupSrc(3) 
    !
    BgradQ = 0.0
    !
#ifdef EULER     
    !
    ! 3D compressible Euler equations 
    !
    BgradQ = 0. 
    !
#endif
    !
#ifdef ELASTICITY
    !
    lam  = Q(10)
    mu   = Q(11) 
    irho = 1./Q(12)
    IF(Q(13)<=1e-3) THEN
        BgradQ = 0.0 
        RETURN 
        !
    ELSE
        ialpha = 1./Q(13) 
        u      = Q(7:9)*ialpha 
    ENDIF 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)
    !
    BgradQ(1) = - (lam+2*mu)*Qx(7) + (lam+2*mu)*u(1)*Qx(13) 
    BgradQ(2) = - lam*Qx(7)        + lam*u(1)*Qx(13)
    BgradQ(3) = - lam*Qx(7)        + lam*u(1)*Qx(13) 
    BgradQ(4) = - mu *Qx(8)        + mu *u(2)*Qx(13) 
    BgradQ(5) =   0.0 
    BgradQ(6) = - mu *Qx(9)        + mu *u(3)*Qx(13) 
    BgradQ(7) = - irho * Qx(1) - 2*Q(1)*irho*Qx(13)  
    BgradQ(8) = - irho * Qx(4) - 2*Q(4)*irho*Qx(13)   
    BgradQ(9) = - irho * Qx(6) - 2*Q(6)*irho*Qx(13)       
    !
    BgradQ(1) = BgradQ(1) - lam*Qy(8)        + lam*u(2)*Qy(13) 
    BgradQ(2) = BgradQ(2) - (lam+2*mu)*Qy(8) + (lam+2*mu)*u(2)*Qy(13) 
    BgradQ(3) = BgradQ(3) - lam*Qy(8)        + lam*u(2)*Qy(13)  
    BgradQ(4) = BgradQ(4) - mu *Qy(7)        + mu *u(1)*Qy(13) 
    BgradQ(5) = BgradQ(5) - mu *Qy(9)        + mu *u(3)*Qy(13) 
    BgradQ(6) = BgradQ(6) - 0.0  
    BgradQ(7) = BgradQ(7) - irho * Qy(4) - 2*Q(4)*irho*Qy(13)  
    BgradQ(8) = BgradQ(8) - irho * Qy(2) - 2*Q(2)*irho*Qy(13)   
    BgradQ(9) = BgradQ(9) - irho * Qy(5) - 2*Q(5)*irho*Qy(13)      
    !
    BgradQ(1) = BgradQ(1) - lam*Qz(9)        + lam*u(3)*Qz(13) 
    BgradQ(2) = BgradQ(2) - lam*Qz(9)        + lam*u(3)*Qz(13) 
    BgradQ(3) = BgradQ(3) - (lam+2*mu)*Qz(9) + (lam+2*mu)*u(3)*Qz(13) 
    BgradQ(4) = BgradQ(4) - 0.0  
    BgradQ(5) = BgradQ(5) - mu *Qz(8)        + mu *u(2)*Qz(13)  
    BgradQ(6) = BgradQ(6) - mu *Qz(7)        + mu *u(1)*Qz(13) 
    BgradQ(7) = BgradQ(7) - irho * Qz(6) - 2*Q(6)*irho*Qz(13)  
    BgradQ(8) = BgradQ(8) - irho * Qz(5) - 2*Q(5)*irho*Qz(13)  
    BgradQ(9) = BgradQ(9) - irho * Qz(3) - 2*Q(3)*irho*Qz(13)     
    ! 
    ! advection equation 
    !BgradQ = 0. 
    !BgradQ(1,1) = Qx(1) 
    !BgradQ(1,2) = Qy(1) 
    !BgradQ(1,3) = Qz(1) 
    ! 
    ! acoustic wave equation 
!    BgradQ(1,1) = -Qx(2) 
!    BgradQ(1,2) = -Qy(3) 
!    BgradQ(1,3) = -Qz(4) 
!    BgradQ(2,1) = -Qx(1) 
!    BgradQ(3,2) = -Qy(1) 
!    BgradQ(4,3) = -Qz(1) 

    
#endif 

#ifdef ACOUSTIC 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)
    !
    ! acoustic wave equation 
    BgradQ(1) = -EQN%c0**2*( Qx(2) + Qy(3) + Qz(4) ) 
    BgradQ(2) = -Qx(1) 
    BgradQ(3) = -Qy(1) 
    BgradQ(4) = -Qz(1) 
    
    ! debug 
    !BgradQ = Qx 
    
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)

    k1 = EQN%CCZ4k1  
    k2 = EQN%CCZ4k2  
    fff = EQN%CCZ4f 
    ggg = EQN%CCZ4g 
    e  = EQN%CCZ4e 
    c  = EQN%CCZ4c  
    ds  = EQN%CCZ4ds 
    xi  = EQN%CCZ4xi 
    sk  = EQN%CCZ4sk
    !
    sknl = sk 
    IF(sk==0.0) THEN
       IF(MAXVAL(ABS(Q(18:20))) < 1e-6 ) THEN
           sknl = 0.0
       ELSE
           sknl = 1.0  
       ENDIF 
    ENDIF
    !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
    g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
    g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
      
    alpha = Q(17)
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
    END SELECT 

    K0    = Q(54)
    dK0   = sk*(/ Qx(54), Qy(54), Qz(54) /) 
    beta0 = 0.0 ! Q(55:57)
    b0    = 0.0 ! Q(58:60)

    !Christoffel(1,1,1) = g_contr(1,1)*Q(36)+g_contr(1,2)*(2*Q(37)-Q(42))+g_contr(1,3)*(2*Q(38)-Q(48))
    !Christoffel(1,1,2) = g_contr(2,1)*Q(36)+g_contr(2,2)*(2*Q(37)-Q(42))+g_contr(2,3)*(2*Q(38)-Q(48))
    !Christoffel(1,1,3) = g_contr(3,1)*Q(36)+g_contr(3,2)*(2*Q(37)-Q(42))+g_contr(3,3)*(2*Q(38)-Q(48))
    !Christoffel(1,2,1) = g_contr(1,1)*Q(42)+g_contr(1,2)*Q(39)+g_contr(1,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(1,2,2) = g_contr(2,1)*Q(42)+g_contr(2,2)*Q(39)+g_contr(2,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(1,2,3) = g_contr(3,1)*Q(42)+g_contr(3,2)*Q(39)+g_contr(3,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(1,3,1) = g_contr(1,1)*Q(48)+g_contr(1,2)*(Q(49)+Q(40)-Q(44))+g_contr(1,3)*Q(41)
    !Christoffel(1,3,2) = g_contr(2,1)*Q(48)+g_contr(2,2)*(Q(49)+Q(40)-Q(44))+g_contr(2,3)*Q(41)
    !Christoffel(1,3,3) = g_contr(3,1)*Q(48)+g_contr(3,2)*(Q(49)+Q(40)-Q(44))+g_contr(3,3)*Q(41)
    !Christoffel(2,1,1) = g_contr(1,1)*Q(42)+g_contr(1,2)*Q(39)+g_contr(1,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(2,1,2) = g_contr(2,1)*Q(42)+g_contr(2,2)*Q(39)+g_contr(2,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(2,1,3) = g_contr(3,1)*Q(42)+g_contr(3,2)*Q(39)+g_contr(3,3)*(Q(44)+Q(40)-Q(49))
    !Christoffel(2,2,1) = g_contr(1,1)*(2*Q(43)-Q(39))+g_contr(1,2)*Q(45)+g_contr(1,3)*(2*Q(46)-Q(51))
    !Christoffel(2,2,2) = g_contr(2,1)*(2*Q(43)-Q(39))+g_contr(2,2)*Q(45)+g_contr(2,3)*(2*Q(46)-Q(51))
    !Christoffel(2,2,3) = g_contr(3,1)*(2*Q(43)-Q(39))+g_contr(3,2)*Q(45)+g_contr(3,3)*(2*Q(46)-Q(51))
    !Christoffel(2,3,1) = g_contr(1,1)*(Q(49)+Q(44)-Q(40))+g_contr(1,2)*Q(51)+g_contr(1,3)*Q(47)
    !Christoffel(2,3,2) = g_contr(2,1)*(Q(49)+Q(44)-Q(40))+g_contr(2,2)*Q(51)+g_contr(2,3)*Q(47)
    !Christoffel(2,3,3) = g_contr(3,1)*(Q(49)+Q(44)-Q(40))+g_contr(3,2)*Q(51)+g_contr(3,3)*Q(47)
    !Christoffel(3,1,1) = g_contr(1,1)*Q(48)+g_contr(1,2)*(Q(49)+Q(40)-Q(44))+g_contr(1,3)*Q(41)
    !Christoffel(3,1,2) = g_contr(2,1)*Q(48)+g_contr(2,2)*(Q(49)+Q(40)-Q(44))+g_contr(2,3)*Q(41)
    !Christoffel(3,1,3) = g_contr(3,1)*Q(48)+g_contr(3,2)*(Q(49)+Q(40)-Q(44))+g_contr(3,3)*Q(41)
    !Christoffel(3,2,1) = g_contr(1,1)*(Q(49)+Q(44)-Q(40))+g_contr(1,2)*Q(51)+g_contr(1,3)*Q(47)
    !Christoffel(3,2,2) = g_contr(2,1)*(Q(49)+Q(44)-Q(40))+g_contr(2,2)*Q(51)+g_contr(2,3)*Q(47)
    !Christoffel(3,2,3) = g_contr(3,1)*(Q(49)+Q(44)-Q(40))+g_contr(3,2)*Q(51)+g_contr(3,3)*Q(47)
    !Christoffel(3,3,1) = g_contr(1,1)*(2*Q(50)-Q(41))+g_contr(1,2)*(2*Q(52)-Q(47))+g_contr(1,3)*Q(53)
    !Christoffel(3,3,2) = g_contr(2,1)*(2*Q(50)-Q(41))+g_contr(2,2)*(2*Q(52)-Q(47))+g_contr(2,3)*Q(53)
    !Christoffel(3,3,3) = g_contr(3,1)*(2*Q(50)-Q(41))+g_contr(3,2)*(2*Q(52)-Q(47))+g_contr(3,3)*Q(53)
    !  
    !
    !dgup(1,1,1) = -2*g_contr(1,1)**2*Q(36)-2*g_contr(1,1)*g_contr(2,1)*Q(37)-2*g_contr(1,1)*g_contr(3,1)*Q(38)-2*g_contr(1,2)*g_contr(1,1)*Q(37)-2*g_contr(1,2)*g_contr(2,1)*Q(39)-2*g_contr(1,2)*g_contr(3,1)*Q(40)-2*g_contr(1,3)*g_contr(1,1)*Q(38)-2*g_contr(1,3)*g_contr(2,1)*Q(40)-2*g_contr(1,3)*g_contr(3,1)*Q(41)
    !dgup(1,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(36)-2*g_contr(1,1)*g_contr(2,2)*Q(37)-2*g_contr(1,1)*g_contr(3,2)*Q(38)-2*g_contr(1,2)**2*Q(37)-2*g_contr(1,2)*g_contr(2,2)*Q(39)-2*g_contr(1,2)*g_contr(3,2)*Q(40)-2*g_contr(1,3)*g_contr(1,2)*Q(38)-2*g_contr(1,3)*g_contr(2,2)*Q(40)-2*g_contr(1,3)*g_contr(3,2)*Q(41)
    !dgup(1,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(36)-2*g_contr(1,1)*g_contr(2,3)*Q(37)-2*g_contr(1,1)*g_contr(3,3)*Q(38)-2*g_contr(1,2)*g_contr(1,3)*Q(37)-2*g_contr(1,2)*g_contr(2,3)*Q(39)-2*g_contr(1,2)*g_contr(3,3)*Q(40)-2*g_contr(1,3)**2*Q(38)-2*g_contr(1,3)*g_contr(2,3)*Q(40)-2*g_contr(1,3)*g_contr(3,3)*Q(41)
    !dgup(1,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(36)-2*g_contr(2,1)**2*Q(37)-2*g_contr(2,1)*g_contr(3,1)*Q(38)-2*g_contr(1,1)*g_contr(2,2)*Q(37)-2*g_contr(2,2)*g_contr(2,1)*Q(39)-2*g_contr(2,2)*g_contr(3,1)*Q(40)-2*g_contr(2,3)*g_contr(1,1)*Q(38)-2*g_contr(2,3)*g_contr(2,1)*Q(40)-2*g_contr(2,3)*g_contr(3,1)*Q(41)
    !dgup(1,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(36)-2*g_contr(2,1)*g_contr(2,2)*Q(37)-2*g_contr(2,1)*g_contr(3,2)*Q(38)-2*g_contr(2,2)*g_contr(1,2)*Q(37)-2*g_contr(2,2)**2*Q(39)-2*g_contr(2,2)*g_contr(3,2)*Q(40)-2*g_contr(2,3)*g_contr(1,2)*Q(38)-2*g_contr(2,3)*g_contr(2,2)*Q(40)-2*g_contr(2,3)*g_contr(3,2)*Q(41)
    !dgup(1,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(36)-2*g_contr(2,1)*g_contr(2,3)*Q(37)-2*g_contr(2,1)*g_contr(3,3)*Q(38)-2*g_contr(2,2)*g_contr(1,3)*Q(37)-2*g_contr(2,2)*g_contr(2,3)*Q(39)-2*g_contr(2,2)*g_contr(3,3)*Q(40)-2*g_contr(2,3)*g_contr(1,3)*Q(38)-2*g_contr(2,3)**2*Q(40)-2*g_contr(2,3)*g_contr(3,3)*Q(41)
    !dgup(1,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(36)-2*g_contr(3,1)*g_contr(2,1)*Q(37)-2*g_contr(3,1)**2*Q(38)-2*g_contr(3,2)*g_contr(1,1)*Q(37)-2*g_contr(3,2)*g_contr(2,1)*Q(39)-2*g_contr(3,2)*g_contr(3,1)*Q(40)-2*g_contr(1,1)*g_contr(3,3)*Q(38)-2*g_contr(3,3)*g_contr(2,1)*Q(40)-2*g_contr(3,3)*g_contr(3,1)*Q(41)
    !dgup(1,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(36)-2*g_contr(3,1)*g_contr(2,2)*Q(37)-2*g_contr(3,1)*g_contr(3,2)*Q(38)-2*g_contr(3,2)*g_contr(1,2)*Q(37)-2*g_contr(3,2)*g_contr(2,2)*Q(39)-2*g_contr(3,2)**2*Q(40)-2*g_contr(3,3)*g_contr(1,2)*Q(38)-2*g_contr(2,2)*g_contr(3,3)*Q(40)-2*g_contr(3,3)*g_contr(3,2)*Q(41)
    !dgup(1,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(36)-2*g_contr(3,1)*g_contr(2,3)*Q(37)-2*g_contr(3,1)*g_contr(3,3)*Q(38)-2*g_contr(3,2)*g_contr(1,3)*Q(37)-2*g_contr(3,2)*g_contr(2,3)*Q(39)-2*g_contr(3,2)*g_contr(3,3)*Q(40)-2*g_contr(3,3)*g_contr(1,3)*Q(38)-2*g_contr(3,3)*g_contr(2,3)*Q(40)-2*g_contr(3,3)**2*Q(41)
    !dgup(2,1,1) = -2*g_contr(1,1)**2*Q(42)-2*g_contr(1,1)*g_contr(2,1)*Q(43)-2*g_contr(1,1)*g_contr(3,1)*Q(44)-2*g_contr(1,2)*g_contr(1,1)*Q(43)-2*g_contr(1,2)*g_contr(2,1)*Q(45)-2*g_contr(1,2)*g_contr(3,1)*Q(46)-2*g_contr(1,3)*g_contr(1,1)*Q(44)-2*g_contr(1,3)*g_contr(2,1)*Q(46)-2*g_contr(1,3)*g_contr(3,1)*Q(47)
    !dgup(2,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(42)-2*g_contr(1,1)*g_contr(2,2)*Q(43)-2*g_contr(1,1)*g_contr(3,2)*Q(44)-2*g_contr(1,2)**2*Q(43)-2*g_contr(1,2)*g_contr(2,2)*Q(45)-2*g_contr(1,2)*g_contr(3,2)*Q(46)-2*g_contr(1,3)*g_contr(1,2)*Q(44)-2*g_contr(1,3)*g_contr(2,2)*Q(46)-2*g_contr(1,3)*g_contr(3,2)*Q(47)
    !dgup(2,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(42)-2*g_contr(1,1)*g_contr(2,3)*Q(43)-2*g_contr(1,1)*g_contr(3,3)*Q(44)-2*g_contr(1,2)*g_contr(1,3)*Q(43)-2*g_contr(1,2)*g_contr(2,3)*Q(45)-2*g_contr(1,2)*g_contr(3,3)*Q(46)-2*g_contr(1,3)**2*Q(44)-2*g_contr(1,3)*g_contr(2,3)*Q(46)-2*g_contr(1,3)*g_contr(3,3)*Q(47)
    !dgup(2,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(42)-2*g_contr(2,1)**2*Q(43)-2*g_contr(2,1)*g_contr(3,1)*Q(44)-2*g_contr(1,1)*g_contr(2,2)*Q(43)-2*g_contr(2,2)*g_contr(2,1)*Q(45)-2*g_contr(2,2)*g_contr(3,1)*Q(46)-2*g_contr(2,3)*g_contr(1,1)*Q(44)-2*g_contr(2,3)*g_contr(2,1)*Q(46)-2*g_contr(2,3)*g_contr(3,1)*Q(47)
    !dgup(2,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(42)-2*g_contr(2,1)*g_contr(2,2)*Q(43)-2*g_contr(2,1)*g_contr(3,2)*Q(44)-2*g_contr(2,2)*g_contr(1,2)*Q(43)-2*g_contr(2,2)**2*Q(45)-2*g_contr(2,2)*g_contr(3,2)*Q(46)-2*g_contr(2,3)*g_contr(1,2)*Q(44)-2*g_contr(2,3)*g_contr(2,2)*Q(46)-2*g_contr(2,3)*g_contr(3,2)*Q(47)
    !dgup(2,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(42)-2*g_contr(2,1)*g_contr(2,3)*Q(43)-2*g_contr(2,1)*g_contr(3,3)*Q(44)-2*g_contr(2,2)*g_contr(1,3)*Q(43)-2*g_contr(2,2)*g_contr(2,3)*Q(45)-2*g_contr(2,2)*g_contr(3,3)*Q(46)-2*g_contr(2,3)*g_contr(1,3)*Q(44)-2*g_contr(2,3)**2*Q(46)-2*g_contr(2,3)*g_contr(3,3)*Q(47)
    !dgup(2,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(42)-2*g_contr(3,1)*g_contr(2,1)*Q(43)-2*g_contr(3,1)**2*Q(44)-2*g_contr(3,2)*g_contr(1,1)*Q(43)-2*g_contr(3,2)*g_contr(2,1)*Q(45)-2*g_contr(3,2)*g_contr(3,1)*Q(46)-2*g_contr(1,1)*g_contr(3,3)*Q(44)-2*g_contr(3,3)*g_contr(2,1)*Q(46)-2*g_contr(3,3)*g_contr(3,1)*Q(47)
    !dgup(2,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(42)-2*g_contr(3,1)*g_contr(2,2)*Q(43)-2*g_contr(3,1)*g_contr(3,2)*Q(44)-2*g_contr(3,2)*g_contr(1,2)*Q(43)-2*g_contr(3,2)*g_contr(2,2)*Q(45)-2*g_contr(3,2)**2*Q(46)-2*g_contr(3,3)*g_contr(1,2)*Q(44)-2*g_contr(2,2)*g_contr(3,3)*Q(46)-2*g_contr(3,3)*g_contr(3,2)*Q(47)
    !dgup(2,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(42)-2*g_contr(3,1)*g_contr(2,3)*Q(43)-2*g_contr(3,1)*g_contr(3,3)*Q(44)-2*g_contr(3,2)*g_contr(1,3)*Q(43)-2*g_contr(3,2)*g_contr(2,3)*Q(45)-2*g_contr(3,2)*g_contr(3,3)*Q(46)-2*g_contr(3,3)*g_contr(1,3)*Q(44)-2*g_contr(3,3)*g_contr(2,3)*Q(46)-2*g_contr(3,3)**2*Q(47)
    !dgup(3,1,1) = -2*g_contr(1,1)**2*Q(48)-2*g_contr(1,1)*g_contr(2,1)*Q(49)-2*g_contr(1,1)*g_contr(3,1)*Q(50)-2*g_contr(1,2)*g_contr(1,1)*Q(49)-2*g_contr(1,2)*g_contr(2,1)*Q(51)-2*g_contr(1,2)*g_contr(3,1)*Q(52)-2*g_contr(1,3)*g_contr(1,1)*Q(50)-2*g_contr(1,3)*g_contr(2,1)*Q(52)-2*g_contr(1,3)*g_contr(3,1)*Q(53)
    !dgup(3,1,2) = -2*g_contr(1,1)*g_contr(1,2)*Q(48)-2*g_contr(1,1)*g_contr(2,2)*Q(49)-2*g_contr(1,1)*g_contr(3,2)*Q(50)-2*g_contr(1,2)**2*Q(49)-2*g_contr(1,2)*g_contr(2,2)*Q(51)-2*g_contr(1,2)*g_contr(3,2)*Q(52)-2*g_contr(1,3)*g_contr(1,2)*Q(50)-2*g_contr(1,3)*g_contr(2,2)*Q(52)-2*g_contr(1,3)*g_contr(3,2)*Q(53)
    !dgup(3,1,3) = -2*g_contr(1,1)*g_contr(1,3)*Q(48)-2*g_contr(1,1)*g_contr(2,3)*Q(49)-2*g_contr(1,1)*g_contr(3,3)*Q(50)-2*g_contr(1,2)*g_contr(1,3)*Q(49)-2*g_contr(1,2)*g_contr(2,3)*Q(51)-2*g_contr(1,2)*g_contr(3,3)*Q(52)-2*g_contr(1,3)**2*Q(50)-2*g_contr(1,3)*g_contr(2,3)*Q(52)-2*g_contr(1,3)*g_contr(3,3)*Q(53)
    !dgup(3,2,1) = -2*g_contr(2,1)*g_contr(1,1)*Q(48)-2*g_contr(2,1)**2*Q(49)-2*g_contr(2,1)*g_contr(3,1)*Q(50)-2*g_contr(1,1)*g_contr(2,2)*Q(49)-2*g_contr(2,2)*g_contr(2,1)*Q(51)-2*g_contr(2,2)*g_contr(3,1)*Q(52)-2*g_contr(2,3)*g_contr(1,1)*Q(50)-2*g_contr(2,3)*g_contr(2,1)*Q(52)-2*g_contr(2,3)*g_contr(3,1)*Q(53)
    !dgup(3,2,2) = -2*g_contr(2,1)*g_contr(1,2)*Q(48)-2*g_contr(2,1)*g_contr(2,2)*Q(49)-2*g_contr(2,1)*g_contr(3,2)*Q(50)-2*g_contr(2,2)*g_contr(1,2)*Q(49)-2*g_contr(2,2)**2*Q(51)-2*g_contr(2,2)*g_contr(3,2)*Q(52)-2*g_contr(2,3)*g_contr(1,2)*Q(50)-2*g_contr(2,3)*g_contr(2,2)*Q(52)-2*g_contr(2,3)*g_contr(3,2)*Q(53)
    !dgup(3,2,3) = -2*g_contr(2,1)*g_contr(1,3)*Q(48)-2*g_contr(2,1)*g_contr(2,3)*Q(49)-2*g_contr(2,1)*g_contr(3,3)*Q(50)-2*g_contr(2,2)*g_contr(1,3)*Q(49)-2*g_contr(2,2)*g_contr(2,3)*Q(51)-2*g_contr(2,2)*g_contr(3,3)*Q(52)-2*g_contr(2,3)*g_contr(1,3)*Q(50)-2*g_contr(2,3)**2*Q(52)-2*g_contr(2,3)*g_contr(3,3)*Q(53)
    !dgup(3,3,1) = -2*g_contr(3,1)*g_contr(1,1)*Q(48)-2*g_contr(3,1)*g_contr(2,1)*Q(49)-2*g_contr(3,1)**2*Q(50)-2*g_contr(3,2)*g_contr(1,1)*Q(49)-2*g_contr(3,2)*g_contr(2,1)*Q(51)-2*g_contr(3,2)*g_contr(3,1)*Q(52)-2*g_contr(1,1)*g_contr(3,3)*Q(50)-2*g_contr(3,3)*g_contr(2,1)*Q(52)-2*g_contr(3,3)*g_contr(3,1)*Q(53)
    !dgup(3,3,2) = -2*g_contr(3,1)*g_contr(1,2)*Q(48)-2*g_contr(3,1)*g_contr(2,2)*Q(49)-2*g_contr(3,1)*g_contr(3,2)*Q(50)-2*g_contr(3,2)*g_contr(1,2)*Q(49)-2*g_contr(3,2)*g_contr(2,2)*Q(51)-2*g_contr(3,2)**2*Q(52)-2*g_contr(3,3)*g_contr(1,2)*Q(50)-2*g_contr(2,2)*g_contr(3,3)*Q(52)-2*g_contr(3,3)*g_contr(3,2)*Q(53)
    !dgup(3,3,3) = -2*g_contr(3,1)*g_contr(1,3)*Q(48)-2*g_contr(3,1)*g_contr(2,3)*Q(49)-2*g_contr(3,1)*g_contr(3,3)*Q(50)-2*g_contr(3,2)*g_contr(1,3)*Q(49)-2*g_contr(3,2)*g_contr(2,3)*Q(51)-2*g_contr(3,2)*g_contr(3,3)*Q(52)-2*g_contr(3,3)*g_contr(1,3)*Q(50)-2*g_contr(3,3)*g_contr(2,3)*Q(52)-2*g_contr(3,3)**2*Q(53)
    !  
    Kex(1,1) = Q(7) 
    Kex(1,2) = Q(8) 
    Kex(1,3) = Q(9) 
    Kex(2,1) = Q(8) 
    Kex(2,2) = Q(10) 
    Kex(2,3) = Q(11) 
    Kex(3,1) = Q(9) 
    Kex(3,2) = Q(11) 
    Kex(3,3) = Q(12) 
    !
    dK(:,1,1) = gradQ(7,:) 
    dK(:,1,2) = gradQ(8,:) 
    dK(:,1,3) = gradQ(9,:) 
    dK(:,2,1) = gradQ(8,:) 
    dK(:,2,2) = gradQ(10,:) 
    dK(:,2,3) = gradQ(11,:) 
    dK(:,3,1) = gradQ(9,:) 
    dK(:,3,2) = gradQ(11,:) 
    dK(:,3,3) = gradQ(12,:) 
    !
    Kmix = MATMUL(g_contr, Kex)
    Kup  = MATMUL(g_contr, Kmix) 
    !
    Z = (/ Q(13), Q(14), Q(15) /)
    dZ(:,1) = gradQ(13,:)
    dZ(:,2) = gradQ(14,:)
    dZ(:,3) = gradQ(15,:)
    Theta = Q(16)
    dTheta = gradQ(16,:) 
    !
    alpha = Q(17) 
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
    !    
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dChristoffelNCP = 0.0
    !
    DO i = 1, 3 
     DO j = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3 
          DO l = 1, 3 
             dChristoffelNCP(k,i,j,m) = dChristoffelNCP(k,i,j,m) + g_contr(m,l)*(0.5*(dDD(k,i,j,l)+dDD(i,k,j,l))+0.5*(dDD(k,j,i,l)+dDD(j,k,i,l))-0.5*(dDD(k,l,i,j)+dDD(l,k,i,j))) 
           ENDDO
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
     ! 
     DO i = 1, 3  
      DO j = 1, 3  
         nablaijalphaNCP(i,j) = alpha*0.5*(dAA(i,j)+dAA(j,i)) 
      ENDDO
     ENDDO     
     !
     nablanablaalphaNCP = SUM( g_contr*nablaijalphaNCP ) 
     !
     DO i = 1, 3 
      DO j = 1, 3 
       DO m = 1, 3 
        DO k = 1, 3 
          RiemannNCP(i,k,j,m) = dChristoffelNCP(k,i,j,m)-dChristoffelNCP(j,i,k,m) 
        ENDDO 
       ENDDO 
      ENDDO 
     ENDDO 
     RicciNCP = 0.0 
     DO m = 1, 3  
      DO n = 1, 3  
       DO l = 1, 3  
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l) 
       ENDDO 
      ENDDO 
     ENDDO      
     !
     nablaZNCP      = dZ 
     RplusNablaZNCP = SUM( g_contr*( RicciNCP + nablaZNCP + TRANSPOSE(nablaZNCP) ) ) 
     DO k = 1, 3
        dtraceKNCP(k) = SUM(g_contr(:,:)*dK(k,:,:)) 
     ENDDO
     dKtempNCP = 0.0 
     DO i = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dKtempNCP(i) = dKtempNCP(i) + g_contr(j,l)*(dK(l,i,j)-dK(i,j,l))   
       ENDDO
      ENDDO
     ENDDO               
     !
     ! Main variables of the Z4 system 
     dtK = nablaijalphaNCP - alpha*( RicciNCP + nablaZNCP + TRANSPOSE(nablaZNCP) ) - beta(1)*dK(1,:,:) - beta(2)*dK(2,:,:) - beta(3)*dK(3,:,:)      ! extrinsic curvature 
     dtTheta = -0.5*alpha*RplusNablaZNCP - beta(1)*dTheta(1) - beta(2)*dTheta(2) - beta(3)*dTheta(3)                                                ! temporal Z 
     dtZ     = -e**2*alpha*dKtempNCP -alpha*dTheta - beta(1)*dZ(1,:) - beta(2)*dZ(2,:) - beta(3)*dZ(3,:)                                            ! spatial Z 
     ! Auxiliary variables 
     dtA = +alpha*fa*( dtraceKNCP(:) -dK0(:) - 2*dTheta(:) ) - beta(1)*dAA(1,:) - beta(2)*dAA(2,:) - beta(3)*dAA(3,:) 
     dtB = 0.0 
     dtD = alpha*dK 
     DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            dtD(k,i,j) = dtD(k,i,j) - 0.5*(g_cov(m,i)*0.5*sknl*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*sknl*(dBB(k,i,m)+dBB(i,k,m)) )  
        ENDDO
       ENDDO
      ENDDO
     ENDDO 
     !
    DO i = 1, 3 
     DO k = 1, 3     
      DO j = 1, 3 
         dtB(k,i) = dtB(k,i) - ggg*alpha**2*g_contr(i,j)*( 0.5*(dAA(k,j)+dAA(j,k)) - dZ(k,j) )  
         dtB(k,i) = dtB(k,i) -   c*alpha**2*g_contr(i,j)*( 0.5*(dAA(k,j)-dAA(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            dtB(k,i) = dtB(k,i) -ggg*alpha**2*g_contr(i,j)*( -2*g_contr(n,l)*0.5*(dDD(k,l,n,j)+dDD(l,k,n,j)) )
            dtB(k,i) = dtB(k,i)   +c*alpha**2*g_contr(i,j)*( -2*g_contr(n,l)*0.5*(dDD(k,l,n,j)-dDD(l,k,n,j)) ) 
          ENDDO 
         ENDDO 
         ! 
       ENDDO 
      ENDDO 
     ENDDO 
     ! 
     dtB = -sk*dtB
     !
     dtD = dtD - beta(1)*dDD(1,:,:,:) - beta(2)*dDD(2,:,:,:) - beta(3)*dDD(3,:,:,:)
     !
     BgradQ(1:6)    = 0.0    ! gamma_ij 
     BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /) 
     BgradQ(13:15)  = dtZ(1:3) 
     BgradQ(16)     = dtTheta 
     BgradQ(17)     = 0.0        ! alpha 
     BgradQ(18:20)  = 0.0        ! beta^i 
     BgradQ(21:23)  = 0.0        ! b^i 
     BgradQ(24:26)  = dtA(1:3) 
     BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /) 
     BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /) 
     BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /) 
     BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)      
     ! 
     !AQx = 0.0
     !BQy = 0.0
     !CQz = 0.0       
     !INCLUDE 'Z4_NCP.inc' 
     !! 
     !! to get a steady solution 
     !!  
     !AQx(24) = AQx(24) - alpha*fa*dK0(1) 
     !BQy(25) = BQy(25) - alpha*fa*dK0(2) 
     !CQz(26) = CQz(26) - alpha*fa*dK0(3)       
     !!
     !BgradQ = AQx + BQy + CQz 
     !!
     !RETURN 
     !
     CONTINUE
     !
#ifdef Z4GRMHD
    ! This coupling leads to a non-hyperbolic system (Jordan blocks!) due to the gradients of alpha, beta and gamma in the MHD part
    ! which have no counterpart with gradients of MHD variables in the Z4 system. 
!    QGRMHD(1:9)   = Q(55:63)    ! hydro variables 
!    QGRMHD(10)    = Q(17)       ! lapse 
!    QGRMHD(11:13) = Q(18:20)    ! shift 
!    QGRMHD(14:19) = Q(1:6)      ! metric 
!    !
!    gradQGRMHD(1:9,:)   = gradQ(55:63,:)    ! hydro variables 
!    gradQGRMHD(10,:)    = gradQ(17,:)       ! lapse 
!    gradQGRMHD(11:13,:) = gradQ(18:20,:)    ! shift 
!    gradQGRMHD(14:19,:) = gradQ(1:6,:)      ! metric 
!    !
!    CALL PDENCPGRMHD(BgradQGRMHD,QGRMHD,gradQGRMHD,par) 
!    BgradQ(55:63) = BgradQGRMHD(1:9) 
    ! 
#endif 
      !  
#endif
    ! ---------------------------------
    !
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)

    k1   = EQN%CCZ4k1  
    k2   = EQN%CCZ4k2  
    k3   = EQN%CCZ4k3   
    fff  = EQN%CCZ4f 
    ggg  = EQN%CCZ4g 
    e    = EQN%CCZ4e 
    itau = EQN%CCZ4itau  
    eta  = EQN%CCZ4eta  
    c    = EQN%CCZ4c 
    mu   = EQN%CCZ4mu 
    ds   = EQN%CCZ4ds 
    bs   = EQN%CCZ4bs 
    xi   = EQN%CCZ4xi 
    sk   = EQN%CCZ4sk
    !
    ! These are the tilde quantities, so be careful !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    ! This determinant should be close to unity, since we use the conformal decomposition 
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    !g_cov = det**(-1./3.) * g_cov 
    !det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
    
    g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
    g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
    g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
    g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
    g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
    g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
    g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
    g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
    g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
      
    alpha = EXP(MAX(-20.,MIN(20.,Q(17))))  
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 
    ! 
    K0    = Q(59)
    dK0   = 0.0 ! sk*gradQ(59,:) 
    !  
    Aex(1,1) = Q(7) 
    Aex(1,2) = Q(8) 
    Aex(1,3) = Q(9) 
    Aex(2,1) = Q(8) 
    Aex(2,2) = Q(10) 
    Aex(2,3) = Q(11) 
    Aex(3,1) = Q(9) 
    Aex(3,2) = Q(11) 
    Aex(3,3) = Q(12) 
    !
    traceA = SUM(g_contr*Aex) 
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = MATMUL(g_contr, Aex)
    Aup  = MATMUL(g_contr, TRANSPOSE(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    b = Q(21:23) 
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   = EXP(MAX(-20.,MIN(20.,Q(55))))   
    dphi  = gradQ(55,:) 
    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
    !
    dBB = dBB*sk    
    ! 
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi**2 + 1./3.*traceK*g_cov/phi**2 
    Kmix = MATMUL( phi**2*g_contr, Kex  ) 
    Kup  = MATMUL( phi**2*g_contr, Kmix ) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       !Christoffel_kind1(i,j,k) =  DD(i,j,k)+DD(j,i,k)-DD(k,i,j)     ! this definition does not work ! 
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
          Gtilde(i)                = Gtilde(i)+2*g_contr(i,j)*g_contr(k,l)*DD(l,j,k) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    Z   = 0.5*MATMUL( g_cov, Ghat - Gtilde ) 
    Zup = MATMUL(phi**2*g_contr, Z) 
    !
    ! 
    !
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          !dChristoffel(k,i,ip,m) = 0 
          dChristoffelNCP(k,i,ip,m) = 0 
          dChristoffel_tildeNCP(k,i,ip,m) = 0 
          DO l = 1, 3 
            !dChristoffel(k,i,ip,m) = dChristoffel(k,i,ip,m) + g_contr(m,l)*(0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)))         & 
            !                                                - g_contr(m,l)*(g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)))  & 
            !                                                 +dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            ! 
            dChristoffelNCP(k,i,ip,m) = dChristoffelNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )         & 
                                                                  - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)) ) 
            ! 
            dChristoffel_tildeNCP(k,i,ip,m) = dChristoffel_tildeNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )           
            ! 
            !dChristoffelSrc(k,i,ip,m) = dChristoffelSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    !Riemann = 0.0 
    !RiemannSrc = 0.0 
    RiemannNCP = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          !Riemann(i,k,ip,m)    = dChristoffel(k,i,ip,m)-dChristoffel(ip,i,k,m)
          RiemannNCP(i,k,ip,m) = dChristoffelNCP(k,i,ip,m)-dChristoffelNCP(ip,i,k,m)
          !RiemannSrc(i,k,ip,m) = dChristoffelSrc(k,i,ip,m)-dChristoffelSrc(ip,i,k,m) 
          DO j = 1, 3
           !Riemann(i,k,ip,m)    = Riemann(i,k,ip,m)    + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
           !RiemannSrc(i,k,ip,m) = RiemannSrc(i,k,ip,m) + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    !Ricci = 0.0 
    RicciNCP = 0.0 
    !RicciSrc = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         !Ricci(m,n) = Ricci(m,n) + Riemann(m,l,n,l)  
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l)  
         !RicciSrc(m,n) = RicciSrc(m,n) + RiemannSrc(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO    
    !
    RNCP = phi**2*SUM(g_contr*RicciNCP) 
    !
    dGtildeNCP = 0.0
    !DO k = 1, 3 
    ! DO j = 1, 3 
    !  DO l = 1, 3 
    !   DO m = 1, 3
    !    DO n = 1, 3 
    !     dGtildeNCP(j,k) = dGtildeNCP(j,k) + 2.0*g_contr(k,l)*g_contr(m,n)*0.5*(dDD(j,n,l,m)+dDD(n,j,l,m)) 
    !    ENDDO
    !   ENDDO 
    !  ENDDO 
    ! ENDDO 
    !ENDDO 
    !
    ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible...      
    !
    DO i = 1, 3 
     DO k = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dGtildeNCP(k,i) = dGtildeNCP(k,i) + g_contr(j,l)*dChristoffel_tildeNCP(k,j,l,i) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    !
    !!RiccitildeNCP = 0.0 
    !!RicciphiNCP   = 0.0 
    !!RiccitildeSrc = 0.0 
    !!RicciphiSrc   = 0.0 
    !!DO j = 1, 3 
    !! DO i = 1, 3 
    !!    RiccitildeNCP(i,j) = 0 
    !!    DO l = 1, 3
    !!     DO m = 1, 3
    !!        RiccitildeNCP(i,j) = RiccitildeNCP(i,j)-g_contr(l,m)*0.5*(dDD(l,m,i,j)+dDD(m,l,i,j))
    !!     ENDDO
    !!    ENDDO 
    !!    DO k = 1, 3 
    !!        RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGtildeNCP(j,k)+g_cov(k,j)*dGtildeNCP(i,k)) 
    !!        !RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGhat(j,k)+g_cov(k,j)*dGhat(i,k)) 
    !!    ENDDO
    !! ENDDO
    !!ENDDO            
    !!! 
    !!DO j = 1, 3 
    !! DO i = 1, 3 
    !!   RicciphiNCP(i,j) = 0.5*(dPP(i,j)+dPP(j,i))
    !!   DO k = 1, 3 
    !!    DO l = 1, 3 
    !!       RicciphiNCP(i,j) = RicciphiNCP(i,j)+g_cov(i,j)*g_contr(l,k)*0.5*(dPP(k,l)+dPP(l,k))
    !!    ENDDO
    !!   ENDDO 
    !! ENDDO
    !!ENDDO    
    !
    dZNCP = 0.0
    !dZSrc = 0.0 
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3    
        dZNCP(k,i) = dZNCP(k,i) + 0.5*g_cov(i,j)*(dGhat(k,j)-dGtildeNCP(k,j))  
        !dZSrc(k,i) = dZSrc(k,i) + DD(k,i,j)*(Ghat(j)-Gtilde(j)) + 0.5*g_cov(i,j)*( -dGtildeSrc(k,j) )
       ENDDO 
      ENDDO 
    ENDDO     
    !
    DO j = 1, 3 
     DO i = 1, 3 
      nablaZNCP(i,j) = dZNCP(i,j)
      !nablaZSrc(i,j) = dZSrc(i,j)
      !DO k = 1, 3 
      !  nablaZSrc(i,j) = nablaZSrc(i,j) - Christoffel(i,j,k)*Z(k) 
      !ENDDO 
     ENDDO
    ENDDO    
    !
    !RicciPlusNablaZNCP = RiccitildeNCP + RicciphiNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    !RicciPlusNablaZSrc = RiccitildeSrc + RicciphiSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RicciPlusNablaZNCP = RicciNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    !RicciPlusNablaZSrc = RicciSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RPlusTwoNablaZNCP = phi**2*SUM(g_contr*RicciPlusNablaZNCP) 
    !RPlusTwoNablaZSrc = phi**2*SUM(g_contr*RicciPlusNablaZSrc) 
    !
    !Kupdown = SUM(Kex*Kup) 
    !temp = RPlusTwoNablaZNCP + RPlusTwoNablaZSrc - KupDown + traceK**2 
    !temp = SUM(phi**2*g_contr*Ricci) - KupDown + traceK**2 
    !
    nablaijalphaNCP = 0.0
    !nablaijalphaSrc = 0.0
    DO j = 1, 3 
     DO i = 1, 3 
       nablaijalphaNCP(i,j) = alpha*0.5*( dAA(i,j)+dAA(j,i) ) 
       !nablaijalphaSrc(i,j) = alpha*AA(j)*AA(i) 
       !DO k = 1, 3 
       !  nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j) - Christoffel(i,j,k)*alpha*AA(k)  
       !ENDDO
     ENDDO
    ENDDO 
    nablanablaalphaNCP = phi**2*SUM( g_contr*nablaijalphaNCP ) 
    !nablanablaalphaSrc = phi**2*SUM( g_contr*nablaijalphaSrc ) 
!    nablanablaalpha = 0.0
!    DO j = 1, 3 
!     DO i = 1, 3
!         nablanablaalpha = nablanablaalpha + phi**2*( g_contr(i,j)*AA(j)*AA(i)*alpha + alpha*g_contr(i,j)*dAA(i,j) - 2*g_contr(i,j)*SUM(g_contr(:,:)*TRANSPOSE(DD(:,j,:))*AA(i)*alpha)- g_contr(i,j)*PP(i)*AA(j)*alpha )  
!     ENDDO
!    ENDDO 
    !
    SecondOrderTermsNCP = -nablaijalphaNCP + alpha*RicciPlusNablaZNCP 
    !SecondOrderTermsSrc = -nablaijalphaSrc + alpha*RicciPlusNablaZSrc 
    traceNCP = SUM( g_contr*SecondOrderTermsNCP ) 
    SecondOrderTermsNCP = SecondOrderTermsNCP - 1./3.*g_cov*traceNCP 
    !traceSrc = SUM( g_contr*SecondOrderTermsSrc ) 
    !SecondOrderTermsSrc = SecondOrderTermsSrc - 1./3.*g_cov*traceSrc 
    !
    ! Now assemble all this terrible stuff... 
    !
    dtgamma = 0.0 
    !
    ! Main variables of the CCZ4 system 
    dtK = phi**2*SecondOrderTermsNCP + beta(1)*dAex(1,:,:) + beta(2)*dAex(2,:,:) + beta(3)*dAex(3,:,:)      ! extrinsic curvature
    !dtK = dtK + phi**2*SecondOrderTermsSrc + alpha*Aex*(traceK-2*Theta) - 2*alpha*MATMUL(Aex,Amix) 
!    DO j = 1, 3 
!     DO i = 1, 3 
!      DO k = 1, 3 
!         dtK(i,j) = dtK(i,j) + Aex(k,i)*BB(j,k) + Aex(k,j)*BB(i,k) - 2./3.*Aex(i,j)*BB(k,k) 
!      ENDDO
!     ENDDO
!    ENDDO 
    !
    dtTraceK = -nablanablaalphaNCP + alpha*RPlusTwoNablaZNCP + SUM(beta(:)*dtraceK(:)) 
    !
    traceB = BB(1,1) + BB(2,2) + BB(3,3) 
    dtphi   = 0.0 
    dtalpha = 0.0 

    Aupdown = SUM(Aex*Aup) 
    dtTheta = 0.5*alpha*e**2*( RplusTwoNablaZNCP ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)        ! *** original cleaning *** 
    !dtTheta = 0.5*e**2*( RplusTwoNablaZNCP ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)             ! *** use turbo cleaning here *** original 0.5*alpha*e**2*RplusTwoNablaZNCP 
    !
    divAupNCP = 0.0
    DO i = 1, 3
        DO j = 1, 3
         DO l = 1, 3
          DO k = 1, 3    
            divAupNCP(i) = divAupNCP(i) + g_contr(i,l)*g_contr(j,k)*dAex(j,l,k) 
          ENDDO
         ENDDO
        ENDDO        
    ENDDO     
    DO i = 1, 3 
        Mom(i) = - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) + divAupNCP(i)  
    ENDDO 
    !
    !!dKex = 0.0
    !!DO j = 1, 3
    !! DO i = 1, 3
    !!  DO k = 1, 3
    !!      dKex(k,i,j) = 1.0/phi**2*( dAex(k,i,j) + 1./3.*dtraceK(k)*g_cov(i,j) ) 
    !!  ENDDO
    !! ENDDO
    !!ENDDO 
    !!! 
    !!Mom(:) = 0.0
    !!DO ii = 1, 3
    !!    DO jj = 1, 3
    !!        DO ll = 1, 3
    !!            Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*(dKex(ll,ii,jj) - dKex(ii,jj,ll))
    !!            !DO mm = 1, 3
    !!            !    Mom(ii) = Mom(ii) + g_contr(jj,ll)*( - Christoffel(jj,ll,mm)*Kex(mm,ii) + Christoffel(jj,ii,mm)*Kex(mm,ll)) 
    !!            !ENDDO
    !!        ENDDO     
    !!    ENDDO
    !!ENDDO     
    !!Mom = MATMUL( g_contr, Mom )     
    !
    DO i = 1, 3
        dtGhat(i) = - 4./3.*alpha*SUM(g_contr(i,:)*dtraceK(:))     &  
        !dtGhat(i)  =  2.0*alpha*( -divAupNCP(i) + ds**2*Mom(i) )  &          
                    + 2.0*alpha*SUM( g_contr(:,i)*( dTheta(:)  ) ) &                    
                    + beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i) 
        DO l = 1, 3
         DO k = 1, 3
             dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) 
         ENDDO
        ENDDO         
    ENDDO 
    DO k = 1, 3 
        ov(k) = 2*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) )           ! here we can use the constraint that trace A tilde = 0. 
    ENDDO
    dtGhat = dtGhat + sk*MATMUL(g_contr,ov)                         ! Ghat is an "up" vector, so we need to multiply with g_contr 
    !
    dtbb = xi*dtGhat + bs*( beta(1)*gradQ(21:23,1) + beta(2)*gradQ(21:23,2) + beta(3)*gradQ(21:23,3) - beta(1)*gradQ(14:16,1) - beta(2)*gradQ(14:16,2) - beta(3)*gradQ(14:16,3) ) 
    dtbb = sk*dtbb  
    !
    dtbeta  = 0.0    
    !
    ! Auxiliary variables 
    dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:)  
    DO k = 1, 3 
        dtA(k) = dtA(k) - sk*alpha*fa*( SUM(g_contr(:,:)*dAex(k,:,:)) )   ! here we can use the constraint that trace A tilde = 0. 
    ENDDO
    !
    ! We have removed the conservative fluxes for CCZ4, so put all the stuff into the NCP and FusedNCP 
    dtB(:,1) = fff*gradQ(21,:)  
    dtB(:,2) = fff*gradQ(22,:)  
    dtB(:,3) = fff*gradQ(23,:)  
    ! #ordB1#     
    ! for the ordering constraints, we have to check whether they should be multiplied by alpha**2, or not... 
    DO i = 1, 3 
     DO k = 1, 3 
      DO j = 1, 3 
         !dtB(k,i) = -mu*g_contr(i,j)*dZNCP(k,j)    
         dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( (dPP(k,j)-dPP(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            dtB(k,i) = dtB(k,i) - mu*alpha**2*g_contr(i,j)*g_contr(n,l)*( dDD(k,l,j,n)-dDD(l,k,j,n) )   
          ENDDO 
         ENDDO
      ENDDO 
      !
      ENDDO 
    ENDDO 
    dtB = dtB + bs*( beta(1)*dBB(1,:,:) + beta(2)*dBB(2,:,:) + beta(3)*dBB(3,:,:) ) 
    ! New stuff 1, which makes appear a Lie derivative and a second order ordering constraint 
    !DO i = 1, 3
    ! DO k = 1, 3 
    !   dtB(k,i) = dtB(k,i) + bs*( beta(1)*dBB(1,k,i) + beta(2)*dBB(2,k,i) + beta(3)*dBB(3,k,i) - beta(1)*dBB(k,1,i) - beta(2)*dBB(k,2,i) - beta(3)*dBB(k,3,i) ) 
    ! ENDDO
    !ENDDO 
    dtB = dtB*sk 
    !
    dtD = -alpha*dAex  
    DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            dtD(k,i,j) = dtD(k,i,j) + ( 0.5*(g_cov(m,i)*0.5*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./3.*g_cov(i,j)*0.5*(dBB(k,m,m)+dBB(m,k,m)) ) 
            DO n = 1, 3
                dtD(k,i,j) = dtD(k,i,j) + 1./3*alpha*g_cov(i,j)*g_contr(n,m)*dAex(k,n,m)     ! explicitly remove the trace of tilde A again 
            ENDDO 
        ENDDO        
       ENDDO
      ENDDO
    ENDDO 
    dtD = dtD + beta(1)*dDD(1,:,:,:) + beta(2)*dDD(2,:,:,:) + beta(3)*dDD(3,:,:,:)
    !
    dtP = + beta(1)*dPP(1,:) + beta(2)*dPP(2,:) + beta(3)*dPP(3,:)    
    DO k = 1, 3 
     dtP(k) = dtP(k) + 1./3.*alpha*dtraceK(k) + sk*1./3.*alpha*SUM(g_contr(:,:)*dAex(k,:,:))  ! use the fact that trace A tilde = 0 
     DO i = 1, 3 
          dtP(k) = dtP(k) - 1./3.*0.5*(dBB(k,i,i)+dBB(i,k,i))  
     ENDDO
    ENDDO 
    !
    BgradQ(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)          ! \tilde \gamma_ij 
    BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                                  ! \tilde A_ij 
    BgradQ(13)     = dtTheta                                                                                           ! Theta       
    BgradQ(14:16)  = dtGhat(1:3)                                                                                       ! \hat \Gamma^i           
    BgradQ(17)     = dtalpha                                                                                           ! log alpha 
    BgradQ(18:20)  = dtbeta                                                                                            ! beta^i 
    BgradQ(21:23)  = dtbb                                                                                              ! b^i 
    BgradQ(24:26)  = dtA(1:3)                                                                                          ! A_k       
    BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /)    ! B_k^i 
    BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                      ! D_kij 
    BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /)                      ! D_kij 
    BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)                      ! D_kij 
    BgradQ(54)     = dtTraceK                                                                                          ! traceK 
    BgradQ(55)     = dtphi                                                                                             ! log phi 
    BgradQ(56:58)  = dtP                                                                                               ! P_k 
    !
    BgradQ = -BgradQ ! change sign, since we work on the left hand side in PDENCP 
    !   
    CONTINUE
    !
#endif
    ! 
    ! ---------------------------------
    !
#if defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)

    k1   = EQN%CCZ4k1  
    k2   = EQN%CCZ4k2  
    k3   = EQN%CCZ4k3   
    fff  = EQN%CCZ4f 
    ggg  = EQN%CCZ4g 
    eta  = EQN%CCZ4eta 
    itau = EQN%CCZ4itau  
    e    = EQN%CCZ4e 
    c    = EQN%CCZ4c 
    mu   = EQN%CCZ4mu 
    ds   = EQN%CCZ4ds 
    bs   = EQN%CCZ4bs 
    xi   = EQN%CCZ4xi 
    sk   = EQN%CCZ4sk
    !
    ! These are the tilde quantities, so be careful !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    ! This determinant should be unity, since we use the conformal decomposition 
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    g_cov = det**(-1./3.) * g_cov 
    det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
    
    g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
    g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
    g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
    g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
    g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
    g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
    g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
    g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
    g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
      
    alpha = EXP(MAX(-20.,MIN(20.,Q(17)))) 
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 
    ! 
    K0    = Q(62)
    dK0   = sk*gradQ(62,:) 
    !  
    Aex(1,1) = Q(7) 
    Aex(1,2) = Q(8) 
    Aex(1,3) = Q(9) 
    Aex(2,1) = Q(8) 
    Aex(2,2) = Q(10) 
    Aex(2,3) = Q(11) 
    Aex(3,1) = Q(9) 
    Aex(3,2) = Q(11) 
    Aex(3,3) = Q(12) 
    !
    traceA = SUM(g_contr*Aex) 
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = MATMUL(g_contr, Aex)
    Aup  = MATMUL(g_contr, TRANSPOSE(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    XX = (/ Q(59), Q(60), Q(61) /) 
    dXX(:,1) = gradQ(59,:) 
    dXX(:,2) = gradQ(60,:) 
    dXX(:,3) = gradQ(61,:) 
    !
    b = Q(21:23) 
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   =  EXP(MAX(-20.,MIN(20.,Q(55))))  
    dphi  = gradQ(55,:) 
    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
    !
    dBB = dBB*sk    
    ! 
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi**2 + 1./3.*traceK*g_cov/phi**2 
    Kmix = MATMUL( phi**2*g_contr, Kex  ) 
    Kup  = MATMUL( phi**2*g_contr, Kmix ) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       !Christoffel_kind1(i,j,k) =  DD(i,j,k)+DD(j,i,k)-DD(k,i,j)     ! this definition does not work ! 
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
          !Gtilde(i)                = Gtilde(i)+2*g_contr(i,j)*g_contr(k,l)*DD(l,j,k) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    !
    DO i = 1, 3
     DO j = 1, 3
      DO l = 1, 3
          Gtilde(i) = Gtilde(i) + g_contr(j,l)*Christoffel_tilde(j,l,i) 
      ENDDO
     ENDDO     
    ENDDO    
    !
    Z   = 0.5*MATMUL( g_cov, Ghat - Gtilde ) 
    Zup = MATMUL(phi**2*g_contr, Z) 
    !
    !dChristoffel    = 0.0 
    dChristoffelNCP = 0.0
    dChristoffelSrc = 0.0 
    dChristoffel_tildeNCP = 0.0
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          DO l = 1, 3 
            !dChristoffel(k,i,ip,m) = dChristoffel(k,i,ip,m) + g_contr(m,l)*(0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)))         & 
            !                                                - g_contr(m,l)*(g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)))  & 
            !                                                 +dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            ! 
            dChristoffelNCP(k,i,ip,m) = dChristoffelNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )         & 
                                                                  - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)) ) 
            !
            dChristoffel_tildeNCP(k,i,ip,m) = dChristoffel_tildeNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) ) 
            ! 
            !dChristoffelOrd(k,i,ip,m) = dChristoffelOrd(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)-dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)-dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)-dDD(l,k,i,ip)) )         & 
            !                                                      - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)-dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)-dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)-dPP(l,k)) )             
            !
            dChristoffelSrc(k,i,ip,m) = dChristoffelSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            !
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    !Riemann = 0.0 
    RiemannSrc = 0.0 
    RiemannNCP = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          !Riemann(i,k,ip,m)    = dChristoffel(k,i,ip,m)-dChristoffel(ip,i,k,m)
          RiemannNCP(i,k,ip,m) = dChristoffelNCP(k,i,ip,m)-dChristoffelNCP(ip,i,k,m)
          RiemannSrc(i,k,ip,m) = dChristoffelSrc(k,i,ip,m)-dChristoffelSrc(ip,i,k,m) 
          DO j = 1, 3
           !Riemann(i,k,ip,m)    = Riemann(i,k,ip,m)    + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
           RiemannSrc(i,k,ip,m) = RiemannSrc(i,k,ip,m) + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    !Ricci = 0.0 
    RicciNCP = 0.0 
    RicciSrc = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         !Ricci(m,n) = Ricci(m,n) + Riemann(m,l,n,l)  
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l)  
         RicciSrc(m,n) = RicciSrc(m,n) + RiemannSrc(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO    
    !
    RNCP = phi**2*SUM(g_contr*RicciNCP) 
    RSrc = phi**2*SUM(g_contr*RicciSrc) 
    !
    dGtildeNCP = 0.0
    dGtildeSrc = 0.0
    !DO k = 1, 3 
    ! DO j = 1, 3 
    !  DO l = 1, 3 
    !   DO m = 1, 3
    !    DO n = 1, 3 
    !     dGtildeNCP(j,k) = dGtildeNCP(j,k) + 2.0*g_contr(k,l)*g_contr(m,n)*0.5*(dDD(j,n,l,m)+dDD(n,j,l,m)) 
    !     dGtildeSrc(j,k) = dGtildeSrc(j,k) + 2.0*( g_contr(m,n)*DD(n,l,m)*dgup(j,k,l) + g_contr(k,l)*DD(n,l,m)*dgup(j,m,n) )         
    !    ENDDO
    !   ENDDO 
    !  ENDDO 
    ! ENDDO 
    !ENDDO 
    !    
    ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible... 
    !
    DO i = 1, 3 
     DO k = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dGtildeNCP(k,i) = dGtildeNCP(k,i)                                        + g_contr(j,l)*dChristoffel_tildeNCP(k,j,l,i) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    ! If you want the original computation of the Ricci tensor according to the CCZ4 paper of Alic et al. 2012, use the following version. 
    ! By default, however, we compute the Ricci tensor ab definitionem from the Riemann tensor and the Christoffel symbols. 
    !
    DO j = 1, 3 
     DO i = 1, 3 
        RiccitildeNCP(i,j) = 0 
        DO l = 1, 3
         DO m = 1, 3
            RiccitildeNCP(i,j) = RiccitildeNCP(i,j)-g_contr(l,m)*0.5*(dDD(l,m,i,j)+dDD(m,l,i,j))
         ENDDO
        ENDDO 
        DO k = 1, 3 
            RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGtildeNCP(j,k)+g_cov(k,j)*dGtildeNCP(i,k)) 
            !RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGhat(j,k)+g_cov(k,j)*dGhat(i,k)) 
        ENDDO
     ENDDO
    ENDDO    
        
    DO j = 1, 3 
     DO i = 1, 3         
        RiccitildeSrc(i,j) = 0
        DO k = 1, 3 
          RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + 0.5*(g_cov(k,i)*dGtildeSrc(j,k)+g_cov(k,j)*dGtildeSrc(i,k)) 
          RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + Ghat(k)*0.5*(Christoffel_kind1(i,j,k)+Christoffel_kind1(j,i,k)) 
          DO l = 1, 3 
           DO m = 1, 3 
            RiccitildeSrc(i,j) = RiccitildeSrc(i,j)+g_contr(l,m)*(Christoffel_tilde(l,i,k)*Christoffel_kind1(j,k,m)+Christoffel_tilde(l,j,k)*Christoffel_kind1(i,k,m)+Christoffel_tilde(i,m,k)*Christoffel_kind1(k,j,l))
           ENDDO
          ENDDO
        ENDDO
     ENDDO
    ENDDO
    
    
    DO j = 1, 3 
     DO i = 1, 3 
       RicciphiNCP(i,j) = 0.5*(dPP(i,j)+dPP(j,i))
       DO k = 1, 3 
        DO l = 1, 3 
           RicciphiNCP(i,j) = RicciphiNCP(i,j)+g_cov(i,j)*g_contr(l,k)*0.5*(dPP(k,l)+dPP(l,k))
        ENDDO
       ENDDO 
     ENDDO
    ENDDO
    
    Pup = MATMUL(g_contr,PP) 
    DO m = 1, 3
        DDcontr(m) = SUM(g_contr*DD(:,m,:)) 
    ENDDO 
    
    DO i = 1, 3 
     DO j = 1, 3 
        RicciphiSrc(i,j) = PP(i)*PP(j) 
        DO k = 1, 3 
         RicciphiSrc(i,j) = RicciphiSrc(i,j) - Christoffel_tilde(i,j,k)*PP(k) - 2*g_cov(i,j)*DDcontr(k)*Pup(k)       
         DO l = 1, 3 
            RicciphiSrc(i,j) = RicciphiSrc(i,j) - g_cov(i,j)*g_contr(l,k)*PP(l)*PP(k) !-g_contr(k,l)*PP(k)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))   
         ENDDO
        ENDDO
     ENDDO
    ENDDO    
    !
    dZNCP = 0.0 
    dZSrc = 0.0 
    DO i = 1, 3
     DO k = 1, 3    
      DO j = 1, 3 
        dZNCP(k,i) = dZNCP(k,i) + 0.5*g_cov(i,j)*( dGhat(k,j)-dGtildeNCP(k,j) )  
        dZSrc(k,i) = dZSrc(k,i) + DD(k,i,j)*(Ghat(j)-Gtilde(j)) + 0.5*g_cov(i,j)*( -dGtildeSrc(k,j) )
       ENDDO 
      ENDDO 
    ENDDO     
    !
    nablaZNCP = dZNCP 
    nablaXNCP = dXX  
    nablaZSrc = 0.0 
    DO j = 1, 3 
     DO i = 1, 3 
      nablaZSrc(i,j) = dZSrc(i,j)
      DO k = 1, 3 
        nablaZSrc(i,j) = nablaZSrc(i,j) - Christoffel(i,j,k)*Z(k) 
      ENDDO 
     ENDDO
    ENDDO    
    !
    RicciPlusNablaZNCP = RiccitildeNCP + RicciphiNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) ! + ( nablaXNCP + TRANSPOSE(nablaXNCP) ) 
    !
    !RicciPlusNablaZNCP = RicciNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    !RicciPlusNablaZSrc = RicciSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RPlusTwoNablaZNCP = phi**2*SUM(g_contr*RicciPlusNablaZNCP) 
    !
    !Kupdown = SUM(Kex*Kup) 
    !temp = RPlusTwoNablaZNCP + RPlusTwoNablaZSrc - KupDown + traceK**2 
    !temp = SUM(phi**2*g_contr*Ricci) - KupDown + traceK**2 
    !
    nablaijalphaNCP = 0.0
    nablaijalphaSrc = 0.0
    DO j = 1, 3 
     DO i = 1, 3 
       nablaijalphaNCP(i,j) = alpha*0.5*( dAA(i,j)+dAA(j,i) ) 
       nablaijalphaSrc(i,j) = alpha*AA(j)*AA(i) 
       DO k = 1, 3 
         nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j) - Christoffel(i,j,k)*alpha*AA(k)  
       ENDDO
     ENDDO
    ENDDO 
    nablanablaalphaNCP = phi**2*SUM( g_contr*nablaijalphaNCP ) 
    nablanablaalphaSrc = phi**2*SUM( g_contr*nablaijalphaSrc ) 
!    nablanablaalpha = 0.0
!    DO j = 1, 3 
!     DO i = 1, 3
!         nablanablaalpha = nablanablaalpha + phi**2*( g_contr(i,j)*AA(j)*AA(i)*alpha + alpha*g_contr(i,j)*dAA(i,j) - 2*g_contr(i,j)*SUM(g_contr(:,:)*TRANSPOSE(DD(:,j,:))*AA(i)*alpha)- g_contr(i,j)*PP(i)*AA(j)*alpha )  
!     ENDDO
!    ENDDO 
    !
    SecondOrderTermsNCP = -nablaijalphaNCP + alpha*RicciPlusNablaZNCP 
    SecondOrderTermsSrc = -nablaijalphaSrc + alpha*RicciPlusNablaZSrc 
    traceNCP = SUM( g_contr*SecondOrderTermsNCP ) 
    SecondOrderTermsNCP = SecondOrderTermsNCP - 1./3.*g_cov*traceNCP 
    traceSrc = SUM( g_contr*SecondOrderTermsSrc ) 
    SecondOrderTermsSrc = SecondOrderTermsSrc - 1./3.*g_cov*traceSrc 
    !
    ! Now assemble all this terrible stuff... 
    !
    dtgamma = 0.0 
    !
    ! Main variables of the CCZ4 system 
    dtK = phi**2*SecondOrderTermsNCP + beta(1)*dAex(1,:,:) + beta(2)*dAex(2,:,:) + beta(3)*dAex(3,:,:)      ! extrinsic curvature
    !
    dtTraceK = - nablanablaalphaNCP + alpha*( RPlusTwoNablaZNCP ) + SUM(beta(:)*dtraceK(:)) 
    !
    dtphi   = 0.0 
    dtalpha = 0.0 

    Aupdown = SUM(Aex*Aup) 
    ! *** original 
    dtTheta = 0.5*alpha*e**2*(RplusTwoNablaZNCP ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)               ! temporal Z 
    !
    dKex = 0.0
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
          dKex(k,i,j) = 1.0/phi**2*( dAex(k,i,j) + 1./3.*dtraceK(k)*g_cov(i,j) ) 
      ENDDO
     ENDDO
    ENDDO 
    !
    Mom(:) = 0.0
    DO ii = 1, 3
        DO jj = 1, 3
            DO ll = 1, 3
                Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*(dKex(ll,ii,jj) - dKex(ii,jj,ll))
            ENDDO     
        ENDDO
    ENDDO   
    !
    dtX = alpha*ds**2*Mom + beta(1)*dXX(1,:) + beta(2)*dXX(2,:) + beta(3)*dXX(3,:) 
    !
    DO i = 1, 3
        dtGhat(i) = 2*alpha*( - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) )    &
                    + 2*alpha*SUM( g_contr(:,i)*( dTheta(:)   ) )       & 
                    + beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i)  
        DO l = 1, 3
         DO k = 1, 3
             dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) 
         ENDDO
        ENDDO         
    ENDDO 
    DO k = 1, 3 
        ov(k) = + 2*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) )    ! here we can use the constraint that trace A tilde = 0.         
    ENDDO        
    dtGhat = dtGhat + MATMUL(g_contr,ov)                                                  ! add the ordering constraint "up" (raised by g_contr) 
    !
    dtbb = xi*dtGhat    !  <= be careful, this damping term -eta*b may be dangerous for the gamma driver, since it may kill the waves that you want !  
    !
    ! Add the following terms if you want shift convection in the PDE for b^i 
    dtbb = dtbb + bs*( - beta(1)*gradQ(14:16,1) - beta(2)*gradQ(14:16,2) - beta(3)*gradQ(14:16,3)  + beta(1)*gradQ(21:23,1) + beta(2)*gradQ(21:23,2) + beta(3)*gradQ(21:23,3)  ) !     
    dtbb = sk*dtbb 
    !
    dtbeta  = 0.0 
    ! Add the following term if you want to have shift convection in the PDE for beta^i 
    !
    ! Auxiliary variables 
    dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:) 
    DO k = 1, 3 
        dtA(k) = dtA(k) - alpha*fa*( SUM(g_contr(:,:)*dAex(k,:,:)) )   ! here we can use the constraint that trace A tilde = 0. 
    ENDDO    
    ! 
    ! #xordB1# 
    ! for the ordering constraints, we have to check whether they should be multiplied by alpha**2, or not... 
    dtB = 0.0 
    DO i = 1, 3 
     DO k = 1, 3     
      DO j = 1, 3 
         !dtB(k,i) = +mu*alpha*g_contr(i,j)*( dZNCP(k,j) + dZSrc(k,j) ) + mu*dgup(k,i,j)*Z(j)      
         dtB(k,i) = dtB(k,i) + mu*alpha*g_contr(i,j)*( (dPP(k,j)-dPP(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            !dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( -2*g_contr(n,l)*0.5*(dDD(k,l,n,j)-dDD(l,k,n,j)) ) 
            dtB(k,i) = dtB(k,i) - mu*alpha*g_contr(i,j)*g_contr(n,l)*( dDD(k,l,j,n)-dDD(l,k,j,n) )   
          ENDDO 
         ENDDO 
      ENDDO 
      !
      ENDDO 
    ENDDO 
    dtB = dtB + bs*( beta(1)*dBB(1,:,:) + beta(2)*dBB(2,:,:) + beta(3)*dBB(3,:,:) ) 
    dtB = dtB*sk 
    !
    dtD = -alpha*dAex  
    DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            dtD(k,i,j) = dtD(k,i,j) + ( 0.5*(g_cov(m,i)*0.5*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./3.*g_cov(i,j)*0.5*(dBB(k,m,m)+dBB(m,k,m)) ) 
            DO n = 1, 3
                dtD(k,i,j) = dtD(k,i,j) + 1./3*alpha*g_cov(i,j)*g_contr(n,m)*dAex(k,n,m)   ! explicitly remove the trace of tilde A again 
            ENDDO 
        ENDDO        
       ENDDO
      ENDDO
    ENDDO     
    !
    dtD = dtD + beta(1)*dDD(1,:,:,:) + beta(2)*dDD(2,:,:,:) + beta(3)*dDD(3,:,:,:)
    !
    dtP = + beta(1)*dPP(1,:) + beta(2)*dPP(2,:) + beta(3)*dPP(3,:)    
    DO k = 1, 3 
     dtP(k) = dtP(k) + 1./3.*alpha*dtraceK(k) + 1./3.*alpha*( SUM(g_contr(:,:)*dAex(k,:,:))  )  
     DO i = 1, 3 
          dtP(k) = dtP(k) - 1./3.*0.5*(dBB(k,i,i)+dBB(i,k,i)) 
     ENDDO
    ENDDO 
    !
    BgradQ(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)          ! \tilde \gamma_ij 
    BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                                  ! \tilde A_ij 
    BgradQ(13)     = dtTheta                                                                                           ! Theta       
    BgradQ(14:16)  = dtGhat(1:3)                                                                                       ! \hat \Gamma^i           
    BgradQ(17)     = dtalpha                                                                                           ! log alpha 
    BgradQ(18:20)  = dtbeta                                                                                            ! beta^i 
    BgradQ(21:23)  = dtbb                                                                                              ! b^i 
    BgradQ(24:26)  = dtA(1:3)                                                                                          ! A_k       
    BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /)    ! B_k^i 
    BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                      ! D_kij 
    BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /)                      ! D_kij 
    BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)                      ! D_kij 
    BgradQ(54)     = dtTraceK                                                                                          ! traceK 
    BgradQ(55)     = dtphi                                                                                             ! log phi 
    BgradQ(56:58)  = dtP                                                                                               ! P_k 
    BgradQ(59:61)  = dtX                                                                                               ! X_k 
    !
    BgradQ = -BgradQ    ! here, we have to change sign, since we work on the left hand side 
    !
    RETURN
    !
#endif 
    !
    ! *********************************
    !
#ifdef GRMHD
    !
    CALL PDENCPGRMHD(BgradQ,Q,gradQ,par) 
    RETURN
    !
#endif 
    !            
#ifdef GPR3D

    AQx = 0. 
    BQy = 0. 
    CQz = 0. 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)    
    !
    u = Q(2:4)/Q(1)  
    ! 
    AQx(6) =  -u(2)*Qx(7) - u(3)*Qx(8) 
    BQy(6) =  +u(2)*Qy(6) 
    CQz(6) =  +u(3)*Qz(6)  
    ! 
    AQx(7) =  +u(1)*Qx(7) 
    BQy(7) =  -u(1)*Qy(6) - u(3)*Qy(8) 
    CQz(7) =  +u(3)*Qz(7) 
    !
    AQx(8) =  +u(1)*Qx(8) 
    BQy(8) =  +u(2)*Qy(8) 
    CQz(8) =  -u(1)*Qz(6) - u(2)*Qz(7) 
    ! 
    AQx(9) =  -u(2)*Qx(10) - u(3)*Qx(11) 
    BQy(9) =  +u(2)*Qy(9) 
    CQz(9) =  +u(3)*Qz(9) 
    !
    AQx(10) = +u(1)*Qx(10)     
    BQy(10) = -u(1)*Qy(9) - u(3)*Qy(11) 
    CQz(10) = +u(3)*Qz(10) 
    !
    AQx(11) = +u(1)*Qx(11) 
    BQy(11) = +u(2)*Qy(11) 
    CQz(11) = -u(1)*Qz(9) - u(2)*Qz(10) 
    !
    AQx(12) = -u(2)*Qx(13) - u(3)*Qx(14)  
    BQy(12) = +u(2)*Qy(12)
    CQz(12) = +u(3)*Qz(12) 
    !
    AQx(13) = +u(1)*Qx(13)  
    BQy(13) = -u(1)*Qy(12) - u(3)*Qy(14) 
    CQz(13) = +u(3)*Qz(13) 
    !
    AQx(14) = +u(1)*Qx(14) 
    BQy(14) = +u(2)*Qy(14) 
    CQz(14) = -u(1)*Qz(12) - u(2)*Qz(13) 
    !
    BgradQ = AQx + BQy + CQz 
    !
#endif 
    !
END SUBROUTINE PDENCP     

SUBROUTINE PDEMatrixB(Bn,Q,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: Bn(nVar,nVar) 
    ! Local variables 
    INTEGER :: i, ml(2)  
    REAL    :: p, irho, lam, mu, ialpha, mval   
    REAL    :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)  
    REAL    :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Bnfile(nVar,nVar), BnDiff(nVar,nVar)   
    REAL    :: k1, k2, fff, ggg, e, ds, cs, xi, sk, sknl, alpha, fa, k0, beta0(3), b0(3)   
    REAL    :: g_cov(3,3), g_contr(3,3), det, Christoffel(3,3,3), dgup(3,3,3), uv(3), bb2   
    REAL    :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL    :: v2,vf(3),uem,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim    
    REAL    :: ev(nVar,d), ncp(nVar), src(nVar), vp(nVar)   
    REAL    :: BnGRMHD(19,19), QGRMHD(19) 
    INTEGER :: j,k,l,m,iErr, count    
    !
    Bn = 0.0
    ! 
    B1 = 0. 
    B2 = 0. 
    B3 = 0. 
    !
#ifdef EULER     
    !
    ! 3D compressible Euler equations 
    !
    B1 = 0. 
    B2 = 0. 
    B3 = 0. 
    !
#endif
    !
#ifdef ELASTICITY
    ! 
    lam  = Q(10) 
    mu   = Q(11) 
    irho = 1./Q(12)
    !
    IF(Q(13)<=1e-3) THEN        
        Bn = 0.0
        RETURN 
    ELSE
        ialpha = 1./Q(13)
        uv     = Q(7:9)*ialpha 
    ENDIF 
    !
    B1(1,7)  = - (lam+2*mu) 
    B1(1,13) = + (lam+2*mu)*uv(1)  
    B1(2,7)  = - lam 
    B1(2,13) = + lam*uv(1) 
    B1(3,7)  = - lam 
    B1(3,13) = + lam*uv(1)  
    B1(4,8)  = - mu
    B1(4,13) = + mu*uv(2)
    B1(6,9)  = - mu
    B1(6,13) = + mu*uv(3)
    B1(7,1)  = - irho
    B1(7,13) = - 2*Q(1)*irho  
    B1(8,4)  = - irho
    B1(8,13) = - 2*Q(4)*irho    
    B1(9,6)  = - irho 
    B1(9,13) = - 2*Q(6)*irho        
    !
    B2(1,8)  = - lam 
    B2(1,13) = + lam*uv(2) 
    B2(2,8)  = - (lam+2*mu) 
    B2(2,13) = + (lam+2*mu)*uv(2) 
    B2(3,8)  = - lam 
    B2(3,13) = + lam*uv(2) 
    B2(4,7)  = - mu  
    B2(4,13) = + mu*uv(1)  
    B2(5,9)  = - mu 
    B2(5,13) = + mu*uv(3) 
    B2(7,4)  = - irho 
    B2(7,13) = - 2*Q(4)*irho 
    B2(8,2)  = - irho 
    B2(8,13) = - 2*Q(2)*irho 
    B2(9,5)  = - irho 
    B2(9,13) = - 2*Q(5)*irho       
    !
    B3(1,9)  = - lam 
    B3(1,13) = + lam*uv(3) 
    B3(2,9)  = - lam 
    B3(2,13) = + lam*uv(3) 
    B3(3,9)  = - (lam+2*mu) 
    B3(3,13) = + (lam+2*mu)*uv(3)
    B3(5,8)  = - mu  
    B3(5,13) = + mu*uv(2)  
    B3(6,7)  = - mu 
    B3(6,13) = + mu*uv(1) 
    B3(7,6)  = - irho 
    B3(7,13) = - 2*Q(6)*irho 
    B3(8,5)  = - irho 
    B3(8,13) = - 2*Q(5)*irho 
    B3(9,3)  = - irho 
    B3(9,13) = - 2*Q(3)*irho     
    
!    B1 = 0. 
!    B2 = 0. 
!    B3 = 0. 
!    B1(1,1) = 1.0 
!    B2(1,1) = 1.0 
!    B3(1,1) = 1.0 

!    B1(1,2) = -1.0 
!    B2(1,3) = -1.0 
!    B3(1,4) = -1.0  
!    B1(2,1) = -1.0  
!    B2(3,1) = -1.0  
!    B3(4,1) = -1.0  

    !            
    Bn = B1*nv(1) + B2*nv(2) + B3*nv(3) 
    !        
#endif 
    !
#ifdef ACOUSTIC
    ! 
    B1 = 0. 
    B2 = 0. 
    B3 = 0. 
    ! 
    B1(1,2) = -1.0*EQN%c0**2 
    B2(1,3) = -1.0*EQN%c0**2 
    B3(1,4) = -1.0*EQN%c0**2  
    B1(2,1) = -1.0  
    B2(3,1) = -1.0  
    B3(4,1) = -1.0  
    !            
    Bn = B1*nv(1) + B2*nv(2) + B3*nv(3) 
    !        
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) || defined(CCZ4GRMHD) || defined(CCZ4EINSTEIN)  || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)  
    !
    ! we compute the matrix by applying the matrix-vector-product to all unit vectors 
    Bn = 0.0
    ev = 0.0 
    !CALL PDENCP(src,Q,ev,par) 
    !CALL PDEFusedSrcNCP(src,Q,ev,par,0.0) 
    CONTINUE
    DO i = 1, nVar
        ev = 0.0
        ev(i,:) = nv 
        CALL PDENCP(Bn(:,i),Q,ev,par) 
        !CALL PDEFusedSrcNCP(Bn(:,i),Q,ev,par,0.0) 
        !Bn(:,i) = -(Bn(:,i) - src(:))   
    ENDDO    
    !
    !OPEN(UNIT=334,FILE='Bn.dat')
    !WRITE(333,*) Bn 
    !CLOSE(333) 
    !
    !OPEN(UNIT=334,FILE='Bn.dat')
    !READ(333,*) BnFile 
    !CLOSE(333) 
    !!
    !BnDiff = Bn - BnFile 
    !ml = MAXLOC(ABS(BnDiff))
    !mval = MAXVAL(ABS(BnDiff)) 
    !
    !CONTINUE 
    !
    RETURN 
    !    
#endif 
    !
#ifdef GPR3D
    
    B1 = 0. 
    B2 = 0. 
    B3 = 0. 
    !
    VP = 0.0
    VP(2:4) = Q(2:4)/Q(1) 
    ! 
    B1(6,7) = -VP(3)
    B1(6,8) = -VP(4)
    B2(6,6) = +VP(3)
    B3(6,6) = +VP(4) 
    ! 
    B1(7,7) = +VP(2)
    B2(7,6) = -VP(2)
    B2(7,8) = -VP(4)
    B3(7,7) = +VP(4)
    !
    B1(8,8) = +VP(2)
    B2(8,8) = +VP(3)
    B3(8,6) = -VP(2)
    B3(8,7) = -VP(3)
    ! 
    B1(9,10) = -VP(3)
    B1(9,11) = -VP(4)
    B2(9,9)  = +VP(3)
    B3(9,9)  = +VP(4) 
    !
    B1(10,10) = +VP(2)
    B2(10, 9) = -VP(2)
    B2(10,11) = -VP(4)
    B3(10,10) = +VP(4)
    !
    B1(11,11) = +VP(2)
    B2(11,11) = +VP(3)
    B3(11, 9) = -VP(2)
    B3(11,10) = -VP(3)
    !
    B1(12,13) = -VP(3)
    B1(12,14) = -VP(4)
    B2(12,12) = +VP(3)
    B3(12,12) = +VP(4)
    !
    B1(13,13) = +VP(2)
    B2(13,12) = -VP(2)
    B2(13,14) = -VP(4)
    B3(13,13) = +VP(4)
    !
    B1(14,14) = +VP(2)
    B2(14,14) = +VP(3)
    B3(14,12) = -VP(2)
    B3(14,13) = -VP(3)
    
    Bn = B1*nv(1) + B2*nv(2) + B3*nv(3) 

#endif 
    ! 
END SUBROUTINE PDEMatrixB 
    
    
SUBROUTINE PDEEigenvalues(Lambda,Q,par,nv) 
    USE typesDef, ONLY : nVar, nParam, d, nDim, EQN 
    USE recipies_mod, ONLY : RG 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d) 
    REAL, INTENT(OUT) :: Lambda(nVar) 
    ! Local variables 
    INTEGER :: i, j, iErr, itemp(nVar) 
    REAL    :: p, u, c, cp, cs, rho0, irho, lam, mu, alpha   
    REAL    :: dfdQ(nVar,nVar), ImLambda(nVar), rtemp(nVar), vp(nVar) 
    REAL    :: AA(nVar,nVar), BB(nVar,nVar), nvv(d), fpp(nVar,d), fmm(nVar,d)    
    REAL    :: Qp(nVar), Qm(nVar), eps, t1, R(nVar,nVar), iR(nVar,nVar)
    REAL    :: AMaple(nVar,nVar), Adiff(nVar,nVar), V(nVar)  
    REAL :: ComputeDet, psi, BV(3),BV_contr(3)
    REAL :: b2_4, lf, lf2M, VdotB 
    REAL :: g_contr(3,3), g_cov(3,3)
    REAl :: shift(3), vf(3), vf_cov(3)
    REAL :: b2,cs2,den,gg
    REAL :: lapse, gp, gm 
    REAL :: rho, gamma1,sft,v2,vn,w ! ,vx,vy,vz
    !
#ifdef EULER     
    !
    u = ( Q(2)*nv(1) + Q(3)*nv(2) + Q(4)*nv(3) )/Q(1)       ! normal velocity 
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    c = SQRT(EQN%gamma*p/Q(1))                              ! sound speed
    !
    Lambda = (/ u-c, u, u, u, u+c /)                        ! The eigenvalues of the Euler equations 
    !
#endif 
    !
#ifdef ELASTICITY
    ! get the local Lamé constants from the local material parameters 
    lam  = Q(10)   
    mu   = Q(11) 
    irho = 1./Q(12)
    !
    cp = SQRT((lam+2*mu)*irho) 
    cs = SQRT(mu*irho) 
    Lambda = 0.0 
    ! 
    Lambda(1) = -cp
    Lambda(2) = -cs 
    Lambda(3) = -cs 
    Lambda(4) = 0. 
    Lambda(5) = 0. 
    Lambda(6) = 0. 
    Lambda(7) = +cs
    Lambda(8) = +cs
    Lambda(9) = +cp

#endif 
    !

#ifdef ACOUSTIC
    Lambda = 0.0 
    Lambda(1) = -EQN%c0 
    Lambda(2) = +EQN%c0 
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) || defined(GRMHD) || defined(CCZ4EINSTEIN)  || defined(CCZ4GRMHD)  || defined(BSSZ4EINSTEIN)  || defined(BSSZ4GRMHD)  
    ! 
#ifdef GRMHD
    Lambda = 0.0
    Lambda(1) = -1.0    
    Lambda(2) = +1.0    
    RETURN 
    !
    ! Eigenvalues as given by Olindo 
    CALL PDECons2Prim(V,Q,iErr)
    rho    = V(1)
    vf_cov = V(2:4)
    p      = V(5)
    !
    !BV(1:3) = V(6:8)   ! wrong!
    BV_contr(1:3) = V(6:8)
    psi = V(9)
    lapse = V(10)
    shift = V(11:13)
    !
    !gammaij = V(14:19) 
    g_cov(1,1) = V(14)
    g_cov(1,2) = V(15)
    g_cov(1,3) = V(16)
    g_cov(2,2) = V(17)
    g_cov(2,3) = V(18)
    g_cov(3,3) = V(19)
    g_cov(2,1) = V(15)
    g_cov(3,1) = V(16)
    g_cov(3,2) = V(18)
    !
    CALL MatrixInverse3x3(g_cov,g_contr,gp)
    gp = SQRT(gp)
    gm = 1./gp
    !  evaluate contr. cov. variables
    vf      = MATMUL(g_contr,vf_cov)
    BV      = MATMUL(g_cov,BV_contr)
    !  evaluate useful quantities
    v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
    lf     = 1.0/sqrt(1.0 - v2)
    lf2m   = 1.0 - v2
    b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
    VdotB     = vf(1)*BV(1) + vf(2)*BV(2) + vf(3)*BV(3)
    b2_4 = b2*lf2m + VdotB  ! this is b^2
    gamma1 = EQN%gamma/(EQN%gamma-1.0) 
    w      = rho + gamma1*p + b2_4 ! this is rho*h + b^2
    cs2    = (EQN%gamma*p + b2_4)/w
    !
    vn     = vf(1)*nv(1) + vf(2)*nv(2) + vf(3)*nv(3)
    sft    = shift(1)*nv(1) + shift(2)*nv(2) + shift(3)*nv(3) 
    gg     = g_contr(1,1)*ABS(nv(1)) + g_contr(2,2)*ABS(nv(2)) + g_contr(3,3)*ABS(nv(3))
    den    = 1.0/(1.0 - v2*cs2)
    IF(SUM(nv**2).EQ.0.) THEN  
        u = SQRT( v2) 
        WRITE(*,*)'Impossible error!'
        STOP
    ELSE
        u = vn 
    ENDIF
    Lambda(1)   = ( u*(1.0-cs2) - SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den
    Lambda(2:4) = u
    Lambda(5)   = ( u*(1.0-cs2) + SQRT( cs2*lf2m*( (1.0-v2*cs2)*gg - u**2*(1.0-cs2) )) )*den 
    Lambda(1:5)   = lapse*Lambda(1:5) - sft
    !
    Lambda(6:) = 0.
    !
    ! ADD MAGNETIC FIELD!!!!
    ! see ECHO1 and CAFE (from Alejandro)
    !
    Lambda(9) =  EQN%DivCleaning_a   !1.  !EQN%ch
    !
#else
#if defined(CCZ4EINSTEIN)  || defined(CCZ4GRMHD) 
    alpha = MAX( 1.0, EXP(Q(17)) )*MAX( 1.0, EXP(Q(55)) )/MIN( SQRT(Q(1)), SQRT(Q(4)), SQRT(Q(6)) ) 
#else
    alpha = 1.0 
#endif
    Lambda = 0.0 
    Lambda(1) = -alpha*MAX(SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) - DOT_PRODUCT(Q(18:20),nv(:))   ! MAX( SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) + SQRT(SUM(Q(18:20)**2)) 
    Lambda(2) = +alpha*MAX(SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) - DOT_PRODUCT(Q(18:20),nv(:))   ! MAX( SQRT(2.0), EQN%CCZ4e, EQN%CCZ4ds ) + SQRT(SUM(Q(18:20)**2)) 
#endif 
    !
    IF(0==1) THEN
        dfdQ  = 0.0 
        !nvv = (/ 0.0, 1.0, 0.0 /) 
        !DO j = 1, nVar 
        !   Qp = Q 
        !   Qm = Q 
        !   eps = MAX( 1e-6, 1e-6*ABS(Q(j)) ) 
        !   Qp(j) = Qp(j) + eps
        !   Qm(j) = Qm(j) - eps 
        !   CALL PDEFlux(fpp,Qp,par) 
        !   CALL PDEFlux(fmm,Qm,par)
        !   dfdQ(1:nVar-nParam,j) = (fpp(1:nVar,1)-fmm(1:nVar,1))/(2*eps)*nvv(1) + (fpp(1:nVar,2)-fmm(1:nVar,2))/(2*eps)*nvv(2) + (fpp(1:nVar,3)-fmm(1:nVar,3))/(2*eps)*nvv(3)  
        !ENDDO 
        CALL PDEMatrixB(BB,Q,nv,par) 
        AA = dfdQ + BB 
        !!
        !!  Compute the eigenstructure numerically 
        !!where(abs(AA)<1e-13) AA = 0.0  
        CALL RG(nVar,nVar,AA,Lambda,ImLambda,1,R,itemp,rtemp,ierr)   
        !RETURN 
        !!where(abs(R)<1e-13) R = 0.0 
        CALL MatrixInverse(nVar,R,iR) 
        t1 = MAXVAL(ABS(iR))       
        !
        nvv = (/ 1., 0., 0. /) 
        CALL PDEMatrixA(AA,Q,nvv,par)
        !OPEN(UNIT=333,FILE='A1.dat',RECL=5000) 
        !DO i = 1, nVar
        !    WRITE(333,*) AA(i,:)
        !ENDDO
        !CLOSE(333) 
        CALL PDEEigenvectors(R,Lambda,iR,Q,nvv,par) 
        t1 = MAXVAL(ABS(iR))         
        !
        !AMaple = MATMUL( R, MATMUL( Lambda, iR ) )         
        !Adiff = AMaple - AA 
        !t1 = MAXVAL( ABS(Adiff) ) 
        !
        return
        !! 
        !AMaple = 0.0 
        !OPEN( UNIT = 333, FILE='asys.txt', RECL=3000 ) 
        !READ(333,*) 
        !DO i = 1, 58 
        ! DO j = 1, 58 
        !     READ(333,*) AMaple(i,j) 
        ! ENDDO
        !ENDDO 
        !CLOSE(333) 
        !Adiff = AMaple - AA 
        !t1 = MAXVAL(ABS(Adiff)) 
        !
        nvv = (/ 0., 1., 0. /) 
        CALL PDEMatrixA(AA,Q,nvv,par) 
        OPEN(UNIT=333,FILE='A2.dat',RECL=5000) 
        DO i = 1, nVar
            WRITE(333,*) AA(i,:)
        ENDDO
        CLOSE(333) 
        CALL PDEEigenvectors(R,Lambda,iR,Q,nvv,par)
        ! 
        t1 = MAXVAL(ABS(iR))         
        !
        nvv = (/ 0., 0., 1. /) 
        CALL PDEMatrixA(AA,Q,nvv,par) 
        OPEN(UNIT=333,FILE='A3.dat',RECL=5000) 
        DO i = 1, nVar
            WRITE(333,*) AA(i,:)
        ENDDO
        CLOSE(333) 
        CALL PDEEigenvectors(R,Lambda,iR,Q,nvv,par)
        !
        t1 = MAXVAL(ABS(iR))         
        !
        CONTINUE 
        !
    ENDIF
    !
#endif 
    !
#ifdef SRMHD
    Lambda = 1.0 
#endif 
    !
!!#ifdef GRMHD
!!    Lambda = 1.0 
!!#endif 
    !
#ifdef GPR3D
    
    CALL PDECons2Prim(Vp,Q,iErr)  
    
    p = Vp(5) 
    c = SQRT( EQN%gamma*(p+EQN%p0)/Vp(1) + 4./3.*EQN%cs**2 + EQN%alpha**2 )  
    IF(SUM(nv**2).EQ.0.) THEN  
        u = SQRT(Vp(2)**2 + Vp(3)**2 + Vp(4)**2 ) 
    ELSE
        u = Vp(2)*nv(1) + Vp(3)*nv(2) + Vp(4)*nv(3)  
    ENDIF   
    ! 
    Lambda = 0. 
    Lambda(1)  = u-c 
    Lambda(2)  = u 
    Lambda(3)  = u+c     
    !
#endif 
    !
END SUBROUTINE PDEEigenvalues

SUBROUTINE PDEMatrixA(An,Q,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: An(nVar,nVar) 
    ! Local variables 
    INTEGER :: j 
    REAL    :: dfdQ(nVar,nVar), Qp(nVar), Qm(nVar), Fpp(nVar,d), Fmm(nVar,d), Bn(nVar,nVar) 
    REAL    :: eps 
    !     
    dfdQ  = 0.0 
    DO j = 1, nVar 
        Qp = Q 
        Qm = Q 
        eps = MAX( 1e-6, 1e-6*ABS(Q(j)) ) 
        Qp(j) = Qp(j) + eps
        Qm(j) = Qm(j) - eps 
        CALL PDEFlux(Fpp,Qp,par) 
        CALL PDEFlux(Fmm,Qm,par)
        dfdQ(1:nVar-nParam,j) = (fpp(1:nVar,1)-fmm(1:nVar,1))/(2*eps)*nv(1) + (fpp(1:nVar,2)-fmm(1:nVar,2))/(2*eps)*nv(2) + (fpp(1:nVar,3)-fmm(1:nVar,3))/(2*eps)*nv(3)  
    ENDDO 
    CALL PDEMatrixB(Bn,Q,nv,par) 
    An = dfdQ + Bn 
    !
    !WHERE(abs(An)<1e-12)
    !    An = 0.0
    !ENDWHERE     
    !
END SUBROUTINE PDEMatrixA 

SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    USE recipies_mod, ONLY : RG 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: R(nVar,nVar), L(nVar), iR(nVar,nVar) 
    ! Local variables 
    INTEGER :: i,j,n,m,iErr,itemp(nVar), itemp47(47), dynvar(35), dynvar47(47), ml(2) 
    REAL    :: An(nVar,nVar), Qp(nVar), Qm(nVar), Fpp(nVar,d), Fmm(nVar,d), Bn(nVar,nVar)         
    REAL    :: t1, ImLambda(nVar), rtemp(nVar), rtemp47(47), ImLambda47(47), Lambda47(47), ATest(35,35) 
    REAL    :: s1, s2, s3, e, fff, mu, alpha, phi, g_cov(3,3), g_contr(3,3), det, R35(35,35), iR35(35,35), L35(35,35), A35(35,35), Adiff(35,35)
    REAL    :: temp35(35,35), temp47(47,47), R47(47,47), L47(47,47), iR47(47,47), A47(47,47), Adiff47(47,47), Atest47(47,47), Id47(47,47), beta(3)  
    REAL    :: mv,mir, LM(nVar,nVar), Afull(nVar,nVar), dABig(nVar,nVar), Lv47(47), Rman47(47,47)  
    REAL, EXTERNAL :: frac 
    CHARACTER(LEN=10)  :: VarName(nVar)
    CHARACTER(LEN=10)  :: Name35(35), Name47(47)  
    !         
#ifdef CCZ4EINSTEIN 
    ! 
    IF(EQN%CCZ4sk==0) THEN    
        ! 
        CALL PDEMatrixA(An,Q,nv,par) 
        !
        dynvar(1:35) = (/ 7,8,9,10,11,12,54,13,14,15,16,24,25,26,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,56,57,58 /) 
        !
        CALL PDEVarName(VarName)
        DO i = 1, 35
            n = dynvar(i) 
            Name35(i) = VarName(n) 
        ENDDO    
        !
        e     = EQN%CCZ4e 
        phi   = EXP(MAX(-20.,MIN(20.,Q(55))))  
        alpha = EXP(MAX(-20.,MIN(20.,Q(17))))     
        !
        DO i = 1, 35 
            DO j = 1, 35 
            n=dynvar(i) 
            m=dynvar(j) 
            ATest(i,j)  = An(n,m) 
            ENDDO
        ENDDO     
        !
        g_cov(1,1) = Q(1)
        g_cov(1,2) = Q(2)
        g_cov(1,3) = Q(3)
        g_cov(2,1) = Q(2)
        g_cov(2,2) = Q(4)
        g_cov(2,3) = Q(5)
        g_cov(3,1) = Q(3)
        g_cov(3,2) = Q(5)
        g_cov(3,3) = Q(6)
        det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
        g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
        g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
        g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
        g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
        g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
        g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
        g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
        g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
        g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
        !
        L35 = 0.0 
        L35(22,22) = +SQRT(g_contr(1,1))*phi*alpha*e 
        L35(23,23) = -SQRT(g_contr(1,1))*phi*alpha*e 
        DO i = 24, 29     
            L35(i,i) = +SQRT(g_contr(1,1))*phi*alpha 
        ENDDO
        DO i = 30, 35
            L35(i,i) = -SQRT(g_contr(1,1))*phi*alpha
        ENDDO 
        !
            R35(1,1) = 0      
            R35(1,2) = 0
            R35(1,3) = 0
            R35(1,4) = 0
            R35(1,5) = 0
            R35(1,6) = 0
            R35(1,7) = g_cov(1,1)/g_cov(3,3)
            R35(1,8) = 0
            R35(1,9) = 0
            R35(1,10) = 0
            R35(1,11) = 0
            R35(1,12) = 0
            R35(1,13) = 0
            R35(1,14) = 0
            R35(1,15) = 0
            R35(1,16) = 0
            R35(1,17) = 0
            R35(1,18) = 0
            R35(1,19) = 0
            R35(1,20) = 0
            R35(1,21) = 0
            R35(1,22) = -1/sqrt(g_contr(1,1))*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))*e*phi
            R35(1,23) = 1/sqrt(g_contr(1,1))*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))*e*phi
            R35(1,24) = -2*phi*g_contr(2,3)/sqrt(g_contr(1,1))
            R35(1,25) = -phi*g_contr(2,2)/sqrt(g_contr(1,1))
            R35(1,26) = -2*phi*g_contr(1,3)/sqrt(g_contr(1,1))
            R35(1,27) = -2*phi*g_contr(1,2)/sqrt(g_contr(1,1))
            R35(1,28) = 0
            R35(1,29) = -phi*g_contr(3,3)/sqrt(g_contr(1,1))
            R35(1,30) = 2*phi*g_contr(2,3)/sqrt(g_contr(1,1))
            R35(1,31) = phi*g_contr(2,2)/sqrt(g_contr(1,1))
            R35(1,32) = 2*phi*g_contr(1,3)/sqrt(g_contr(1,1))
            R35(1,33) = 2*phi*g_contr(1,2)/sqrt(g_contr(1,1))
            R35(1,34) = 0
            R35(1,35) = phi*g_contr(3,3)/sqrt(g_contr(1,1))
            R35(2,1) = 0
            R35(2,2) = 0
            R35(2,3) = 0
            R35(2,4) = 0
            R35(2,5) = 0
            R35(2,6) = 0
            R35(2,7) = g_cov(1,2)/g_cov(3,3)
            R35(2,8) = 0
            R35(2,9) = 0
            R35(2,10) = 0
            R35(2,11) = 0
            R35(2,12) = 0
            R35(2,13) = 0
            R35(2,14) = 0
            R35(2,15) = 0
            R35(2,16) = 0
            R35(2,17) = 0
            R35(2,18) = 0
            R35(2,19) = 0
            R35(2,20) = 0
            R35(2,21) = 0
            R35(2,22) = sqrt(g_contr(1,1))*g_cov(1,2)*e*phi
            R35(2,23) = -sqrt(g_contr(1,1))*g_cov(1,2)*e*phi
            R35(2,24) = 0
            R35(2,25) = 0
            R35(2,26) = 0
            R35(2,27) = phi*sqrt(g_contr(1,1))
            R35(2,28) = 0
            R35(2,29) = 0
            R35(2,30) = 0
            R35(2,31) = 0
            R35(2,32) = 0
            R35(2,33) = -phi*sqrt(g_contr(1,1))
            R35(2,34) = 0
            R35(2,35) = 0
            R35(3,1) = 0
            R35(3,2) = 0
            R35(3,3) = 0
            R35(3,4) = 0
            R35(3,5) = 0
            R35(3,6) = 0
            R35(3,7) = g_cov(1,3)/g_cov(3,3)
            R35(3,8) = 0
            R35(3,9) = 0
            R35(3,10) = 0
            R35(3,11) = 0
            R35(3,12) = 0
            R35(3,13) = 0
            R35(3,14) = 0
            R35(3,15) = 0
            R35(3,16) = 0
            R35(3,17) = 0
            R35(3,18) = 0
            R35(3,19) = 0
            R35(3,20) = 0
            R35(3,21) = 0
            R35(3,22) = sqrt(g_contr(1,1))*g_cov(1,3)*e*phi
            R35(3,23) = -sqrt(g_contr(1,1))*g_cov(1,3)*e*phi
            R35(3,24) = 0
            R35(3,25) = 0
            R35(3,26) = phi*sqrt(g_contr(1,1))
            R35(3,27) = 0
            R35(3,28) = 0
            R35(3,29) = 0
            R35(3,30) = 0
            R35(3,31) = 0
            R35(3,32) = -phi*sqrt(g_contr(1,1))
            R35(3,33) = 0
            R35(3,34) = 0
            R35(3,35) = 0
            R35(4,1) = 0
            R35(4,2) = 0
            R35(4,3) = 0
            R35(4,4) = 0
            R35(4,5) = 0
            R35(4,6) = 0
            R35(4,7) = g_cov(2,2)/g_cov(3,3)
            R35(4,8) = 0
            R35(4,9) = 0
            R35(4,10) = 0
            R35(4,11) = 0
            R35(4,12) = 0
            R35(4,13) = 0
            R35(4,14) = 0
            R35(4,15) = 0
            R35(4,16) = 0
            R35(4,17) = 0
            R35(4,18) = 0
            R35(4,19) = 0
            R35(4,20) = 0
            R35(4,21) = 0
            R35(4,22) = sqrt(g_contr(1,1))*g_cov(2,2)*e*phi
            R35(4,23) = -sqrt(g_contr(1,1))*g_cov(2,2)*e*phi
            R35(4,24) = 0
            R35(4,25) = phi*sqrt(g_contr(1,1))
            R35(4,26) = 0
            R35(4,27) = 0
            R35(4,28) = 0
            R35(4,29) = 0
            R35(4,30) = 0
            R35(4,31) = -phi*sqrt(g_contr(1,1))
            R35(4,32) = 0
            R35(4,33) = 0
            R35(4,34) = 0
            R35(4,35) = 0
            R35(5,1) = 0
            R35(5,2) = 0
            R35(5,3) = 0
            R35(5,4) = 0
            R35(5,5) = 0
            R35(5,6) = 0
            R35(5,7) = g_cov(2,3)/g_cov(3,3)
            R35(5,8) = 0
            R35(5,9) = 0
            R35(5,10) = 0
            R35(5,11) = 0
            R35(5,12) = 0
            R35(5,13) = 0
            R35(5,14) = 0
            R35(5,15) = 0
            R35(5,16) = 0
            R35(5,17) = 0
            R35(5,18) = 0
            R35(5,19) = 0
            R35(5,20) = 0
            R35(5,21) = 0
            R35(5,22) = sqrt(g_contr(1,1))*g_cov(2,3)*e*phi
            R35(5,23) = -sqrt(g_contr(1,1))*g_cov(2,3)*e*phi
            R35(5,24) = phi*sqrt(g_contr(1,1))
            R35(5,25) = 0
            R35(5,26) = 0
            R35(5,27) = 0
            R35(5,28) = 0
            R35(5,29) = 0
            R35(5,30) = -phi*sqrt(g_contr(1,1))
            R35(5,31) = 0
            R35(5,32) = 0
            R35(5,33) = 0
            R35(5,34) = 0
            R35(5,35) = 0
            R35(6,1) = 0
            R35(6,2) = 0
            R35(6,3) = 0
            R35(6,4) = 0
            R35(6,5) = 0
            R35(6,6) = 0
            R35(6,7) = 1
            R35(6,8) = 0
            R35(6,9) = 0
            R35(6,10) = 0
            R35(6,11) = 0
            R35(6,12) = 0
            R35(6,13) = 0
            R35(6,14) = 0
            R35(6,15) = 0
            R35(6,16) = 0
            R35(6,17) = 0
            R35(6,18) = 0
            R35(6,19) = 0
            R35(6,20) = 0
            R35(6,21) = 0
            R35(6,22) = sqrt(g_contr(1,1))*g_cov(3,3)*e*phi
            R35(6,23) = -sqrt(g_contr(1,1))*g_cov(3,3)*e*phi
            R35(6,24) = 0
            R35(6,25) = 0
            R35(6,26) = 0
            R35(6,27) = 0
            R35(6,28) = 0
            R35(6,29) = phi*sqrt(g_contr(1,1))
            R35(6,30) = 0
            R35(6,31) = 0
            R35(6,32) = 0
            R35(6,33) = 0
            R35(6,34) = 0
            R35(6,35) = -phi*sqrt(g_contr(1,1))
            R35(7,1) = 0
            R35(7,2) = 0
            R35(7,3) = 0
            R35(7,4) = 0
            R35(7,5) = 0
            R35(7,6) = 0
            R35(7,7) = 0
            R35(7,8) = 0
            R35(7,9) = 0
            R35(7,10) = 0
            R35(7,11) = 0
            R35(7,12) = 0
            R35(7,13) = 0
            R35(7,14) = 0
            R35(7,15) = 0
            R35(7,16) = 0
            R35(7,17) = 0
            R35(7,18) = 0
            R35(7,19) = 0
            R35(7,20) = 0
            R35(7,21) = 0
            R35(7,22) = -3.D0/2.D0*sqrt(g_contr(1,1))*(e**2-1)*e*phi
            R35(7,23) = 3.D0/2.D0*sqrt(g_contr(1,1))*(e**2-1)*e*phi
            R35(7,24) = 0
            R35(7,25) = 0
            R35(7,26) = 0
            R35(7,27) = 0
            R35(7,28) = 0
            R35(7,29) = 0
            R35(7,30) = 0
            R35(7,31) = 0
            R35(7,32) = 0
            R35(7,33) = 0
            R35(7,34) = 0
            R35(7,35) = 0
            R35(8,1) = 2*g_contr(1,1)*g_contr(1,2)
            R35(8,2) = 2*g_contr(1,3)**2
            R35(8,3) = 0
            R35(8,4) = 0
            R35(8,5) = 2*g_contr(1,2)*g_contr(1,3)
            R35(8,6) = g_contr(1,1)*g_contr(1,3)
            R35(8,7) = 0
            R35(8,8) = 0
            R35(8,9) = 2*g_contr(1,1)*g_contr(1,3)
            R35(8,10) = -g_contr(1,3)*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))
            R35(8,11) = 2*g_contr(1,2)**2
            R35(8,12) = 2*g_contr(1,2)*g_contr(1,3)
            R35(8,13) = 0
            R35(8,14) = -g_contr(1,2)*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))
            R35(8,15) = 0
            R35(8,16) = 0
            R35(8,17) = 0
            R35(8,18) = g_contr(1,1)*g_contr(1,2)
            R35(8,19) = -g_contr(1,1)*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))
            R35(8,20) = 0
            R35(8,21) = g_contr(1,1)**2
            R35(8,22) = (3*e**2-7)*g_contr(1,1)
            R35(8,23) = (3*e**2-7)*g_contr(1,1)
            R35(8,24) = 0
            R35(8,25) = 0
            R35(8,26) = 0
            R35(8,27) = 0
            R35(8,28) = -4*g_contr(1,1)
            R35(8,29) = 0
            R35(8,30) = 0
            R35(8,31) = 0
            R35(8,32) = 0
            R35(8,33) = 0
            R35(8,34) = -4*g_contr(1,1)
            R35(8,35) = 0
            R35(9,1) = 2*g_contr(1,1)*g_contr(2,2)
            R35(9,2) = 2*g_contr(1,3)*g_contr(2,3)
            R35(9,3) = 0
            R35(9,4) = 0
            R35(9,5) = 2*g_contr(1,3)*g_contr(2,2)
            R35(9,6) = g_contr(1,2)*g_contr(1,3)
            R35(9,7) = 0
            R35(9,8) = g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
            R35(9,9) = 2*g_contr(1,1)*g_contr(2,3)
            R35(9,10) = -g_contr(2,3)*(g_cov(1,2)*g_contr(1,2)+g_cov(2,2)*g_contr(2,2)+g_cov(2,3)*g_contr(2,3))+g_cov(1,1)*g_contr(1,3)*g_contr(1,2)
            R35(9,11) = 2*g_contr(1,2)*g_contr(2,2)
            R35(9,12) = 2*g_contr(1,2)*g_contr(2,3)
            R35(9,13) = 0
            R35(9,14) = -g_contr(2,2)*(g_cov(1,2)*g_contr(1,2)+g_cov(2,2)*g_contr(2,2)+g_cov(2,3)*g_contr(2,3))+g_cov(1,1)*g_contr(1,2)**2
            R35(9,15) = 0
            R35(9,16) = 0
            R35(9,17) = -g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
            R35(9,18) = g_contr(1,2)**2
            R35(9,19) = -g_contr(1,2)*(g_cov(1,2)*g_contr(1,2)+g_cov(2,2)*g_contr(2,2)+g_cov(2,3)*g_contr(2,3))+g_cov(1,1)*g_contr(1,2)*g_contr(1,1)
            R35(9,20) = 0
            R35(9,21) = g_contr(1,1)*g_contr(1,2)
            R35(9,22) = (3*e**2-7)*g_contr(1,2)
            R35(9,23) = (3*e**2-7)*g_contr(1,2)
            R35(9,24) = 0
            R35(9,25) = 0
            R35(9,26) = 0
            R35(9,27) = 0
            R35(9,28) = -4*g_contr(1,2)
            R35(9,29) = 0
            R35(9,30) = 0
            R35(9,31) = 0
            R35(9,32) = 0
            R35(9,33) = 0
            R35(9,34) = -4*g_contr(1,2)
            R35(9,35) = 0
            R35(10,1) = 2*g_contr(1,1)*g_contr(3,2)
            R35(10,2) = 2*g_contr(1,3)*g_contr(3,3)
            R35(10,3) = 0
            R35(10,4) = 0
            R35(10,5) = 2*g_contr(1,3)*g_contr(3,2)
            R35(10,6) = g_contr(1,3)**2
            R35(10,7) = 0
            R35(10,8) = -g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
            R35(10,9) = 2*g_contr(1,1)*g_contr(3,3)
            R35(10,10) = -g_contr(3,3)*(g_cov(1,3)*g_contr(1,3)+g_cov(2,3)*g_contr(2,3)+g_cov(3,3)*g_contr(3,3))+g_cov(1,1)*g_contr(1,3)**2
            R35(10,11) = 2*g_contr(1,2)*g_contr(3,2)
            R35(10,12) = 2*g_contr(1,2)*g_contr(3,3)
            R35(10,13) = 0
            R35(10,14) = -g_contr(2,3)*(g_cov(1,3)*g_contr(1,3)+g_cov(2,3)*g_contr(2,3)+g_cov(3,3)*g_contr(3,3))+g_cov(1,1)*g_contr(1,3)*g_contr(1,2)
            R35(10,15) = 0
            R35(10,16) = 0
            R35(10,17) = g_cov(2,2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
            R35(10,18) = g_contr(1,2)*g_contr(1,3)
            R35(10,19) = -g_contr(1,3)*(g_cov(1,3)*g_contr(1,3)+g_cov(2,3)*g_contr(2,3)+g_cov(3,3)*g_contr(3,3))+g_cov(1,1)*g_contr(1,3)*g_contr(1,1)
            R35(10,20) = 0
            R35(10,21) = g_contr(1,1)*g_contr(1,3)
            R35(10,22) = (3*e**2-7)*g_contr(1,3)
            R35(10,23) = (3*e**2-7)*g_contr(1,3)
            R35(10,24) = 0
            R35(10,25) = 0
            R35(10,26) = 0
            R35(10,27) = 0
            R35(10,28) = -4*g_contr(1,3)
            R35(10,29) = 0
            R35(10,30) = 0
            R35(10,31) = 0
            R35(10,32) = 0
            R35(10,33) = 0
            R35(10,34) = -4*g_contr(1,3)
            R35(10,35) = 0
            R35(11,1) = 0
            R35(11,2) = 0
            R35(11,3) = 0
            R35(11,4) = 0
            R35(11,5) = 0
            R35(11,6) = 0
            R35(11,7) = 0
            R35(11,8) = -g_contr(1,2)/g_contr(1,1)
            R35(11,9) = 0
            R35(11,10) = 0
            R35(11,11) = 0
            R35(11,12) = 0
            R35(11,13) = 0
            R35(11,14) = 0
            R35(11,15) = 0
            R35(11,16) = 0
            R35(11,17) = -g_contr(1,3)/g_contr(1,1)
            R35(11,18) = 0
            R35(11,19) = 0
            R35(11,20) = 0
            R35(11,21) = 0
            R35(11,22) = -3
            R35(11,23) = -3
            R35(11,24) = 0
            R35(11,25) = 0
            R35(11,26) = 0
            R35(11,27) = 0
            R35(11,28) = -3
            R35(11,29) = 0
            R35(11,30) = 0
            R35(11,31) = 0
            R35(11,32) = 0
            R35(11,33) = 0
            R35(11,34) = -3
            R35(11,35) = 0
            R35(12,1) = 0
            R35(12,2) = 0
            R35(12,3) = 0
            R35(12,4) = 0
            R35(12,5) = 0
            R35(12,6) = 0
            R35(12,7) = 0
            R35(12,8) = 1
            R35(12,9) = 0
            R35(12,10) = 0
            R35(12,11) = 0
            R35(12,12) = 0
            R35(12,13) = 0
            R35(12,14) = 0
            R35(12,15) = 0
            R35(12,16) = 0
            R35(12,17) = 0
            R35(12,18) = 0
            R35(12,19) = 0
            R35(12,20) = 0
            R35(12,21) = 0
            R35(12,22) = 0
            R35(12,23) = 0
            R35(12,24) = 0
            R35(12,25) = 0
            R35(12,26) = 0
            R35(12,27) = 0
            R35(12,28) = 0
            R35(12,29) = 0
            R35(12,30) = 0
            R35(12,31) = 0
            R35(12,32) = 0
            R35(12,33) = 0
            R35(12,34) = 0
            R35(12,35) = 0
            R35(13,1) = 0
            R35(13,2) = 0
            R35(13,3) = 0
            R35(13,4) = 0
            R35(13,5) = 0
            R35(13,6) = 0
            R35(13,7) = 0
            R35(13,8) = 0
            R35(13,9) = 0
            R35(13,10) = 0
            R35(13,11) = 0
            R35(13,12) = 0
            R35(13,13) = 0
            R35(13,14) = 0
            R35(13,15) = 0
            R35(13,16) = 0
            R35(13,17) = 1
            R35(13,18) = 0
            R35(13,19) = 0
            R35(13,20) = 0
            R35(13,21) = 0
            R35(13,22) = 0
            R35(13,23) = 0
            R35(13,24) = 0
            R35(13,25) = 0
            R35(13,26) = 0
            R35(13,27) = 0
            R35(13,28) = 0
            R35(13,29) = 0
            R35(13,30) = 0
            R35(13,31) = 0
            R35(13,32) = 0
            R35(13,33) = 0
            R35(13,34) = 0
            R35(13,35) = 0
            R35(14,1) = 0
            R35(14,2) = 0
            R35(14,3) = 0
            R35(14,4) = 0
            R35(14,5) = 0
            R35(14,6) = 0
            R35(14,7) = 0
            R35(14,8) = 0
            R35(14,9) = 0
            R35(14,10) = 0
            R35(14,11) = 0
            R35(14,12) = 0
            R35(14,13) = 0
            R35(14,14) = 0
            R35(14,15) = 0
            R35(14,16) = 0
            R35(14,17) = 0
            R35(14,18) = 0
            R35(14,19) = 0
            R35(14,20) = 0
            R35(14,21) = 1
            R35(14,22) = -1/g_contr(1,1)*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))
            R35(14,23) = -1/g_contr(1,1)*(2*g_cov(1,1)*g_contr(1,1)+3*g_cov(1,2)*g_contr(1,2)+3*g_cov(1,3)*g_contr(1,3))
            R35(14,24) = -2*g_contr(2,3)/g_contr(1,1)
            R35(14,25) = -g_contr(2,2)/g_contr(1,1)
            R35(14,26) = -2*g_contr(1,3)/g_contr(1,1)
            R35(14,27) = -2*g_contr(1,2)/g_contr(1,1)
            R35(14,28) = 0
            R35(14,29) = -g_contr(3,3)/g_contr(1,1)
            R35(14,30) = -2*g_contr(2,3)/g_contr(1,1)
            R35(14,31) = -g_contr(2,2)/g_contr(1,1)
            R35(14,32) = -2*g_contr(1,3)/g_contr(1,1)
            R35(14,33) = -2*g_contr(1,2)/g_contr(1,1)
            R35(14,34) = 0
            R35(14,35) = -g_contr(3,3)/g_contr(1,1)
            R35(15,1) = 1
            R35(15,2) = 0
            R35(15,3) = 0
            R35(15,4) = 0
            R35(15,5) = 0
            R35(15,6) = 0
            R35(15,7) = 0
            R35(15,8) = 0
            R35(15,9) = 0
            R35(15,10) = 0
            R35(15,11) = 0
            R35(15,12) = 0
            R35(15,13) = 0
            R35(15,14) = 0
            R35(15,15) = 0
            R35(15,16) = 0
            R35(15,17) = 0
            R35(15,18) = 0
            R35(15,19) = 0
            R35(15,20) = 0
            R35(15,21) = 0
            R35(15,22) = g_cov(1,2)
            R35(15,23) = g_cov(1,2)
            R35(15,24) = 0
            R35(15,25) = 0
            R35(15,26) = 0
            R35(15,27) = 1
            R35(15,28) = 0
            R35(15,29) = 0
            R35(15,30) = 0
            R35(15,31) = 0
            R35(15,32) = 0
            R35(15,33) = 1
            R35(15,34) = 0
            R35(15,35) = 0
            R35(16,1) = 0
            R35(16,2) = 0
            R35(16,3) = 0
            R35(16,4) = 0
            R35(16,5) = 0
            R35(16,6) = 0
            R35(16,7) = 0
            R35(16,8) = 0
            R35(16,9) = 1
            R35(16,10) = 0
            R35(16,11) = 0
            R35(16,12) = 0
            R35(16,13) = 0
            R35(16,14) = 0
            R35(16,15) = 0
            R35(16,16) = 0
            R35(16,17) = 0
            R35(16,18) = 0
            R35(16,19) = 0
            R35(16,20) = 0
            R35(16,21) = 0
            R35(16,22) = g_cov(1,3)
            R35(16,23) = g_cov(1,3)
            R35(16,24) = 0
            R35(16,25) = 0
            R35(16,26) = 1
            R35(16,27) = 0
            R35(16,28) = 0
            R35(16,29) = 0
            R35(16,30) = 0
            R35(16,31) = 0
            R35(16,32) = 1
            R35(16,33) = 0
            R35(16,34) = 0
            R35(16,35) = 0
            R35(17,1) = 0
            R35(17,2) = 0
            R35(17,3) = -g_contr(1,2)/g_contr(1,1)
            R35(17,4) = 0
            R35(17,5) = 0
            R35(17,6) = 0
            R35(17,7) = 0
            R35(17,8) = 0
            R35(17,9) = 0
            R35(17,10) = g_cov(2,2)*g_contr(1,3)/g_contr(1,1)
            R35(17,11) = 0
            R35(17,12) = 0
            R35(17,13) = -g_contr(1,3)/g_contr(1,1)
            R35(17,14) = g_cov(2,2)*g_contr(1,2)/g_contr(1,1)
            R35(17,15) = 0
            R35(17,16) = 0
            R35(17,17) = 0
            R35(17,18) = 0
            R35(17,19) = g_cov(2,2)
            R35(17,20) = 0
            R35(17,21) = 0
            R35(17,22) = g_cov(2,2)
            R35(17,23) = g_cov(2,2)
            R35(17,24) = 0
            R35(17,25) = 1
            R35(17,26) = 0
            R35(17,27) = 0
            R35(17,28) = 0
            R35(17,29) = 0
            R35(17,30) = 0
            R35(17,31) = 1
            R35(17,32) = 0
            R35(17,33) = 0
            R35(17,34) = 0
            R35(17,35) = 0
            R35(18,1) = 0
            R35(18,2) = 0
            R35(18,3) = 0
            R35(18,4) = -g_contr(1,3)/g_contr(1,1)
            R35(18,5) = 0
            R35(18,6) = 0
            R35(18,7) = 0
            R35(18,8) = 0
            R35(18,9) = 0
            R35(18,10) = g_cov(2,3)*g_contr(1,3)/g_contr(1,1)
            R35(18,11) = 0
            R35(18,12) = 0
            R35(18,13) = 0
            R35(18,14) = g_cov(2,3)*g_contr(1,2)/g_contr(1,1)
            R35(18,15) = 0
            R35(18,16) = -g_contr(1,2)/g_contr(1,1)
            R35(18,17) = 0
            R35(18,18) = 0
            R35(18,19) = g_cov(2,3)
            R35(18,20) = 0
            R35(18,21) = 0
            R35(18,22) = g_cov(2,3)
            R35(18,23) = g_cov(2,3)
            R35(18,24) = 1
            R35(18,25) = 0
            R35(18,26) = 0
            R35(18,27) = 0
            R35(18,28) = 0
            R35(18,29) = 0
            R35(18,30) = 1
            R35(18,31) = 0
            R35(18,32) = 0
            R35(18,33) = 0
            R35(18,34) = 0
            R35(18,35) = 0
            R35(19,1) = 0
            R35(19,2) = 0
            R35(19,3) = 0
            R35(19,4) = 0
            R35(19,5) = 0
            R35(19,6) = 0
            R35(19,7) = 0
            R35(19,8) = 0
            R35(19,9) = 0
            R35(19,10) = g_cov(3,3)*g_contr(1,3)/g_contr(1,1)
            R35(19,11) = 0
            R35(19,12) = 0
            R35(19,13) = 0
            R35(19,14) = g_cov(3,3)*g_contr(1,2)/g_contr(1,1)
            R35(19,15) = -g_contr(1,2)/g_contr(1,1)
            R35(19,16) = 0
            R35(19,17) = 0
            R35(19,18) = 0
            R35(19,19) = g_cov(3,3)
            R35(19,20) = -g_contr(1,3)/g_contr(1,1)
            R35(19,21) = 0
            R35(19,22) = g_cov(3,3)
            R35(19,23) = g_cov(3,3)
            R35(19,24) = 0
            R35(19,25) = 0
            R35(19,26) = 0
            R35(19,27) = 0
            R35(19,28) = 0
            R35(19,29) = 1
            R35(19,30) = 0
            R35(19,31) = 0
            R35(19,32) = 0
            R35(19,33) = 0
            R35(19,34) = 0
            R35(19,35) = 1
            R35(20,1) = 0
            R35(20,2) = 0
            R35(20,3) = 0
            R35(20,4) = 0
            R35(20,5) = 0
            R35(20,6) = 0
            R35(20,7) = 0
            R35(20,8) = 0
            R35(20,9) = 0
            R35(20,10) = 0
            R35(20,11) = 0
            R35(20,12) = 0
            R35(20,13) = 0
            R35(20,14) = 0
            R35(20,15) = 0
            R35(20,16) = 0
            R35(20,17) = 0
            R35(20,18) = 1
            R35(20,19) = 0
            R35(20,20) = 0
            R35(20,21) = 0
            R35(20,22) = 0
            R35(20,23) = 0
            R35(20,24) = 0
            R35(20,25) = 0
            R35(20,26) = 0
            R35(20,27) = 0
            R35(20,28) = 0
            R35(20,29) = 0
            R35(20,30) = 0
            R35(20,31) = 0
            R35(20,32) = 0
            R35(20,33) = 0
            R35(20,34) = 0
            R35(20,35) = 0
            R35(21,1) = 0
            R35(21,2) = 0
            R35(21,3) = 0
            R35(21,4) = 0
            R35(21,5) = 0
            R35(21,6) = 0
            R35(21,7) = 0
            R35(21,8) = 0
            R35(21,9) = 0
            R35(21,10) = 0
            R35(21,11) = 1
            R35(21,12) = 0
            R35(21,13) = 0
            R35(21,14) = 0
            R35(21,15) = 0
            R35(21,16) = 0
            R35(21,17) = 0
            R35(21,18) = 0
            R35(21,19) = 0
            R35(21,20) = 0
            R35(21,21) = 0
            R35(21,22) = 0
            R35(21,23) = 0
            R35(21,24) = 0
            R35(21,25) = 0
            R35(21,26) = 0
            R35(21,27) = 0
            R35(21,28) = 0
            R35(21,29) = 0
            R35(21,30) = 0
            R35(21,31) = 0
            R35(21,32) = 0
            R35(21,33) = 0
            R35(21,34) = 0
            R35(21,35) = 0
            R35(22,1) = 0
            R35(22,2) = 0
            R35(22,3) = 0
            R35(22,4) = 0
            R35(22,5) = 0
            R35(22,6) = 0
            R35(22,7) = 0
            R35(22,8) = 0
            R35(22,9) = 0
            R35(22,10) = 0
            R35(22,11) = 0
            R35(22,12) = 1
            R35(22,13) = 0
            R35(22,14) = 0
            R35(22,15) = 0
            R35(22,16) = 0
            R35(22,17) = 0
            R35(22,18) = 0
            R35(22,19) = 0
            R35(22,20) = 0
            R35(22,21) = 0
            R35(22,22) = 0
            R35(22,23) = 0
            R35(22,24) = 0
            R35(22,25) = 0
            R35(22,26) = 0
            R35(22,27) = 0
            R35(22,28) = 0
            R35(22,29) = 0
            R35(22,30) = 0
            R35(22,31) = 0
            R35(22,32) = 0
            R35(22,33) = 0
            R35(22,34) = 0
            R35(22,35) = 0
            R35(23,1) = 0
            R35(23,2) = 0
            R35(23,3) = 1
            R35(23,4) = 0
            R35(23,5) = 0
            R35(23,6) = 0
            R35(23,7) = 0
            R35(23,8) = 0
            R35(23,9) = 0
            R35(23,10) = 0
            R35(23,11) = 0
            R35(23,12) = 0
            R35(23,13) = 0
            R35(23,14) = 0
            R35(23,15) = 0
            R35(23,16) = 0
            R35(23,17) = 0
            R35(23,18) = 0
            R35(23,19) = 0
            R35(23,20) = 0
            R35(23,21) = 0
            R35(23,22) = 0
            R35(23,23) = 0
            R35(23,24) = 0
            R35(23,25) = 0
            R35(23,26) = 0
            R35(23,27) = 0
            R35(23,28) = 0
            R35(23,29) = 0
            R35(23,30) = 0
            R35(23,31) = 0
            R35(23,32) = 0
            R35(23,33) = 0
            R35(23,34) = 0
            R35(23,35) = 0
            R35(24,1) = 0
            R35(24,2) = 0
            R35(24,3) = 0
            R35(24,4) = 0
            R35(24,5) = 0
            R35(24,6) = 0
            R35(24,7) = 0
            R35(24,8) = 0
            R35(24,9) = 0
            R35(24,10) = 0
            R35(24,11) = 0
            R35(24,12) = 0
            R35(24,13) = 0
            R35(24,14) = 0
            R35(24,15) = 0
            R35(24,16) = 1
            R35(24,17) = 0
            R35(24,18) = 0
            R35(24,19) = 0
            R35(24,20) = 0
            R35(24,21) = 0
            R35(24,22) = 0
            R35(24,23) = 0
            R35(24,24) = 0
            R35(24,25) = 0
            R35(24,26) = 0
            R35(24,27) = 0
            R35(24,28) = 0
            R35(24,29) = 0
            R35(24,30) = 0
            R35(24,31) = 0
            R35(24,32) = 0
            R35(24,33) = 0
            R35(24,34) = 0
            R35(24,35) = 0
            R35(25,1) = 0
            R35(25,2) = 0
            R35(25,3) = 0
            R35(25,4) = 0
            R35(25,5) = 0
            R35(25,6) = 0
            R35(25,7) = 0
            R35(25,8) = 0
            R35(25,9) = 0
            R35(25,10) = 0
            R35(25,11) = 0
            R35(25,12) = 0
            R35(25,13) = 0
            R35(25,14) = 0
            R35(25,15) = 1
            R35(25,16) = 0
            R35(25,17) = 0
            R35(25,18) = 0
            R35(25,19) = 0
            R35(25,20) = 0
            R35(25,21) = 0
            R35(25,22) = 0
            R35(25,23) = 0
            R35(25,24) = 0
            R35(25,25) = 0
            R35(25,26) = 0
            R35(25,27) = 0
            R35(25,28) = 0
            R35(25,29) = 0
            R35(25,30) = 0
            R35(25,31) = 0
            R35(25,32) = 0
            R35(25,33) = 0
            R35(25,34) = 0
            R35(25,35) = 0
            R35(26,1) = 0
            R35(26,2) = 0
            R35(26,3) = 0
            R35(26,4) = 0
            R35(26,5) = 0
            R35(26,6) = 1
            R35(26,7) = 0
            R35(26,8) = 0
            R35(26,9) = 0
            R35(26,10) = 0
            R35(26,11) = 0
            R35(26,12) = 0
            R35(26,13) = 0
            R35(26,14) = 0
            R35(26,15) = 0
            R35(26,16) = 0
            R35(26,17) = 0
            R35(26,18) = 0
            R35(26,19) = 0
            R35(26,20) = 0
            R35(26,21) = 0
            R35(26,22) = 0
            R35(26,23) = 0
            R35(26,24) = 0
            R35(26,25) = 0
            R35(26,26) = 0
            R35(26,27) = 0
            R35(26,28) = 0
            R35(26,29) = 0
            R35(26,30) = 0
            R35(26,31) = 0
            R35(26,32) = 0
            R35(26,33) = 0
            R35(26,34) = 0
            R35(26,35) = 0
            R35(27,1) = 0
            R35(27,2) = 0
            R35(27,3) = 0
            R35(27,4) = 0
            R35(27,5) = 1
            R35(27,6) = 0
            R35(27,7) = 0
            R35(27,8) = 0
            R35(27,9) = 0
            R35(27,10) = 0
            R35(27,11) = 0
            R35(27,12) = 0
            R35(27,13) = 0
            R35(27,14) = 0
            R35(27,15) = 0
            R35(27,16) = 0
            R35(27,17) = 0
            R35(27,18) = 0
            R35(27,19) = 0
            R35(27,20) = 0
            R35(27,21) = 0
            R35(27,22) = 0
            R35(27,23) = 0
            R35(27,24) = 0
            R35(27,25) = 0
            R35(27,26) = 0
            R35(27,27) = 0
            R35(27,28) = 0
            R35(27,29) = 0
            R35(27,30) = 0
            R35(27,31) = 0
            R35(27,32) = 0
            R35(27,33) = 0
            R35(27,34) = 0
            R35(27,35) = 0
            R35(28,1) = 0
            R35(28,2) = 1
            R35(28,3) = 0
            R35(28,4) = 0
            R35(28,5) = 0
            R35(28,6) = 0
            R35(28,7) = 0
            R35(28,8) = 0
            R35(28,9) = 0
            R35(28,10) = 0
            R35(28,11) = 0
            R35(28,12) = 0
            R35(28,13) = 0
            R35(28,14) = 0
            R35(28,15) = 0
            R35(28,16) = 0
            R35(28,17) = 0
            R35(28,18) = 0
            R35(28,19) = 0
            R35(28,20) = 0
            R35(28,21) = 0
            R35(28,22) = 0
            R35(28,23) = 0
            R35(28,24) = 0
            R35(28,25) = 0
            R35(28,26) = 0
            R35(28,27) = 0
            R35(28,28) = 0
            R35(28,29) = 0
            R35(28,30) = 0
            R35(28,31) = 0
            R35(28,32) = 0
            R35(28,33) = 0
            R35(28,34) = 0
            R35(28,35) = 0
            R35(29,1) = 0
            R35(29,2) = 0
            R35(29,3) = 0
            R35(29,4) = 0
            R35(29,5) = 0
            R35(29,6) = 0
            R35(29,7) = 0
            R35(29,8) = 0
            R35(29,9) = 0
            R35(29,10) = 0
            R35(29,11) = 0
            R35(29,12) = 0
            R35(29,13) = 1
            R35(29,14) = 0
            R35(29,15) = 0
            R35(29,16) = 0
            R35(29,17) = 0
            R35(29,18) = 0
            R35(29,19) = 0
            R35(29,20) = 0
            R35(29,21) = 0
            R35(29,22) = 0
            R35(29,23) = 0
            R35(29,24) = 0
            R35(29,25) = 0
            R35(29,26) = 0
            R35(29,27) = 0
            R35(29,28) = 0
            R35(29,29) = 0
            R35(29,30) = 0
            R35(29,31) = 0
            R35(29,32) = 0
            R35(29,33) = 0
            R35(29,34) = 0
            R35(29,35) = 0
            R35(30,1) = 0
            R35(30,2) = 0
            R35(30,3) = 0
            R35(30,4) = 1
            R35(30,5) = 0
            R35(30,6) = 0
            R35(30,7) = 0
            R35(30,8) = 0
            R35(30,9) = 0
            R35(30,10) = 0
            R35(30,11) = 0
            R35(30,12) = 0
            R35(30,13) = 0
            R35(30,14) = 0
            R35(30,15) = 0
            R35(30,16) = 0
            R35(30,17) = 0
            R35(30,18) = 0
            R35(30,19) = 0
            R35(30,20) = 0
            R35(30,21) = 0
            R35(30,22) = 0
            R35(30,23) = 0
            R35(30,24) = 0
            R35(30,25) = 0
            R35(30,26) = 0
            R35(30,27) = 0
            R35(30,28) = 0
            R35(30,29) = 0
            R35(30,30) = 0
            R35(30,31) = 0
            R35(30,32) = 0
            R35(30,33) = 0
            R35(30,34) = 0
            R35(30,35) = 0
            R35(31,1) = 0
            R35(31,2) = 0
            R35(31,3) = 0
            R35(31,4) = 0
            R35(31,5) = 0
            R35(31,6) = 0
            R35(31,7) = 0
            R35(31,8) = 0
            R35(31,9) = 0
            R35(31,10) = 0
            R35(31,11) = 0
            R35(31,12) = 0
            R35(31,13) = 0
            R35(31,14) = 0
            R35(31,15) = 0
            R35(31,16) = 0
            R35(31,17) = 0
            R35(31,18) = 0
            R35(31,19) = 0
            R35(31,20) = 1
            R35(31,21) = 0
            R35(31,22) = 0
            R35(31,23) = 0
            R35(31,24) = 0
            R35(31,25) = 0
            R35(31,26) = 0
            R35(31,27) = 0
            R35(31,28) = 0
            R35(31,29) = 0
            R35(31,30) = 0
            R35(31,31) = 0
            R35(31,32) = 0
            R35(31,33) = 0
            R35(31,34) = 0
            R35(31,35) = 0
            R35(32,1) = 0
            R35(32,2) = 0
            R35(32,3) = 0
            R35(32,4) = 0
            R35(32,5) = 0
            R35(32,6) = 0
            R35(32,7) = 0
            R35(32,8) = 0
            R35(32,9) = 0
            R35(32,10) = 0
            R35(32,11) = 0
            R35(32,12) = 0
            R35(32,13) = 0
            R35(32,14) = 0
            R35(32,15) = 0
            R35(32,16) = 0
            R35(32,17) = 0
            R35(32,18) = 0
            R35(32,19) = 0
            R35(32,20) = 0
            R35(32,21) = 0
            R35(32,22) = -3*sqrt(g_contr(1,1))*e*phi
            R35(32,23) = 3*sqrt(g_contr(1,1))*e*phi
            R35(32,24) = 0
            R35(32,25) = 0
            R35(32,26) = 0
            R35(32,27) = 0
            R35(32,28) = -3*phi*sqrt(g_contr(1,1))
            R35(32,29) = 0
            R35(32,30) = 0
            R35(32,31) = 0
            R35(32,32) = 0
            R35(32,33) = 0
            R35(32,34) = 3*phi*sqrt(g_contr(1,1))
            R35(32,35) = 0
            R35(33,1) = 0
            R35(33,2) = 0
            R35(33,3) = 0
            R35(33,4) = 0
            R35(33,5) = 0
            R35(33,6) = 0
            R35(33,7) = 0
            R35(33,8) = 0
            R35(33,9) = 0
            R35(33,10) = 0
            R35(33,11) = 0
            R35(33,12) = 0
            R35(33,13) = 0
            R35(33,14) = 0
            R35(33,15) = 0
            R35(33,16) = 0
            R35(33,17) = 0
            R35(33,18) = 0
            R35(33,19) = 1
            R35(33,20) = 0
            R35(33,21) = 0
            R35(33,22) = 1
            R35(33,23) = 1
            R35(33,24) = 0
            R35(33,25) = 0
            R35(33,26) = 0
            R35(33,27) = 0
            R35(33,28) = 1
            R35(33,29) = 0
            R35(33,30) = 0
            R35(33,31) = 0
            R35(33,32) = 0
            R35(33,33) = 0
            R35(33,34) = 1
            R35(33,35) = 0
            R35(34,1) = 0
            R35(34,2) = 0
            R35(34,3) = 0
            R35(34,4) = 0
            R35(34,5) = 0
            R35(34,6) = 0
            R35(34,7) = 0
            R35(34,8) = 0
            R35(34,9) = 0
            R35(34,10) = 0
            R35(34,11) = 0
            R35(34,12) = 0
            R35(34,13) = 0
            R35(34,14) = 1
            R35(34,15) = 0
            R35(34,16) = 0
            R35(34,17) = 0
            R35(34,18) = 0
            R35(34,19) = 0
            R35(34,20) = 0
            R35(34,21) = 0
            R35(34,22) = 0
            R35(34,23) = 0
            R35(34,24) = 0
            R35(34,25) = 0
            R35(34,26) = 0
            R35(34,27) = 0
            R35(34,28) = 0
            R35(34,29) = 0
            R35(34,30) = 0
            R35(34,31) = 0
            R35(34,32) = 0
            R35(34,33) = 0
            R35(34,34) = 0
            R35(34,35) = 0
            R35(35,1) = 0
            R35(35,2) = 0
            R35(35,3) = 0
            R35(35,4) = 0
            R35(35,5) = 0
            R35(35,6) = 0
            R35(35,7) = 0
            R35(35,8) = 0
            R35(35,9) = 0
            R35(35,10) = 1
            R35(35,11) = 0
            R35(35,12) = 0
            R35(35,13) = 0
            R35(35,14) = 0
            R35(35,15) = 0
            R35(35,16) = 0
            R35(35,17) = 0
            R35(35,18) = 0
            R35(35,19) = 0
            R35(35,20) = 0
            R35(35,21) = 0
            R35(35,22) = 0
            R35(35,23) = 0
            R35(35,24) = 0
            R35(35,25) = 0
            R35(35,26) = 0
            R35(35,27) = 0
            R35(35,28) = 0
            R35(35,29) = 0
            R35(35,30) = 0
            R35(35,31) = 0
            R35(35,32) = 0
            R35(35,33) = 0
            R35(35,34) = 0
            R35(35,35) = 0
            
            temp35 = R35
            temp35(7,:) = R35(32,:)
            DO i = 8, 32  
              temp35(i,:) = R35(i-1,:) 
            ENDDO
            R35 = temp35 
        !
        !CALL MatrixInverse(35,R35,iR35)  
        !
            iR35(1,1) = 0      
            iR35(1,2) = 0
            iR35(1,3) = 0
            iR35(1,4) = 0
            iR35(1,5) = 0
            iR35(1,6) = 0
            iR35(1,7) = 0
            iR35(1,8) = g_cov(1,2)/g_contr(1,1)/2
            iR35(1,9) = g_cov(2,2)/g_contr(1,1)/2
            iR35(1,10) = g_cov(2,3)/g_contr(1,1)/2
            iR35(1,11) = g_cov(1,2)/3
            iR35(1,12) = -1/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)/g_contr(1,1)/3
            iR35(1,13) = g_cov(1,2)*g_contr(1,3)/g_contr(1,1)/3
            iR35(1,14) = 0
            iR35(1,15) = 0
            iR35(1,16) = 0
            iR35(1,17) = 0
            iR35(1,18) = 0
            iR35(1,19) = 0
            iR35(1,20) = 0
            iR35(1,21) = -g_contr(1,2)/g_contr(1,1)
            iR35(1,22) = 0
            iR35(1,23) = 0
            iR35(1,24) = 0
            iR35(1,25) = 0
            iR35(1,26) = 0
            iR35(1,27) = -g_contr(1,3)/g_contr(1,1)
            iR35(1,28) = 0
            iR35(1,29) = 0
            iR35(1,30) = 0
            iR35(1,31) = 0
            iR35(1,32) = 0
            iR35(1,33) = g_cov(1,2)
            iR35(1,34) = 1/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)/g_contr(1,1)
            iR35(1,35) = g_cov(1,2)*g_contr(1,3)/g_contr(1,1)
            iR35(2,1) = 0
            iR35(2,2) = 0
            iR35(2,3) = 0
            iR35(2,4) = 0
            iR35(2,5) = 0
            iR35(2,6) = 0
            iR35(2,7) = 0
            iR35(2,8) = 0
            iR35(2,9) = 0
            iR35(2,10) = 0
            iR35(2,11) = 0
            iR35(2,12) = 0
            iR35(2,13) = 0
            iR35(2,14) = 0
            iR35(2,15) = 0
            iR35(2,16) = 0
            iR35(2,17) = 0
            iR35(2,18) = 0
            iR35(2,19) = 0
            iR35(2,20) = 0
            iR35(2,21) = 0
            iR35(2,22) = 0
            iR35(2,23) = 0
            iR35(2,24) = 0
            iR35(2,25) = 0
            iR35(2,26) = 0
            iR35(2,27) = 0
            iR35(2,28) = 1
            iR35(2,29) = 0
            iR35(2,30) = 0
            iR35(2,31) = 0
            iR35(2,32) = 0
            iR35(2,33) = 0
            iR35(2,34) = 0
            iR35(2,35) = 0
            iR35(3,1) = 0
            iR35(3,2) = 0
            iR35(3,3) = 0
            iR35(3,4) = 0
            iR35(3,5) = 0
            iR35(3,6) = 0
            iR35(3,7) = 0
            iR35(3,8) = 0
            iR35(3,9) = 0
            iR35(3,10) = 0
            iR35(3,11) = 0
            iR35(3,12) = 0
            iR35(3,13) = 0
            iR35(3,14) = 0
            iR35(3,15) = 0
            iR35(3,16) = 0
            iR35(3,17) = 0
            iR35(3,18) = 0
            iR35(3,19) = 0
            iR35(3,20) = 0
            iR35(3,21) = 0
            iR35(3,22) = 0
            iR35(3,23) = 1
            iR35(3,24) = 0
            iR35(3,25) = 0
            iR35(3,26) = 0
            iR35(3,27) = 0
            iR35(3,28) = 0
            iR35(3,29) = 0
            iR35(3,30) = 0
            iR35(3,31) = 0
            iR35(3,32) = 0
            iR35(3,33) = 0
            iR35(3,34) = 0
            iR35(3,35) = 0
            iR35(4,1) = 0
            iR35(4,2) = 0
            iR35(4,3) = 0
            iR35(4,4) = 0
            iR35(4,5) = 0
            iR35(4,6) = 0
            iR35(4,7) = 0
            iR35(4,8) = 0
            iR35(4,9) = 0
            iR35(4,10) = 0
            iR35(4,11) = 0
            iR35(4,12) = 0
            iR35(4,13) = 0
            iR35(4,14) = 0
            iR35(4,15) = 0
            iR35(4,16) = 0
            iR35(4,17) = 0
            iR35(4,18) = 0
            iR35(4,19) = 0
            iR35(4,20) = 0
            iR35(4,21) = 0
            iR35(4,22) = 0
            iR35(4,23) = 0
            iR35(4,24) = 0
            iR35(4,25) = 0
            iR35(4,26) = 0
            iR35(4,27) = 0
            iR35(4,28) = 0
            iR35(4,29) = 0
            iR35(4,30) = 1
            iR35(4,31) = 0
            iR35(4,32) = 0
            iR35(4,33) = 0
            iR35(4,34) = 0
            iR35(4,35) = 0
            iR35(5,1) = 0
            iR35(5,2) = 0
            iR35(5,3) = 0
            iR35(5,4) = 0
            iR35(5,5) = 0
            iR35(5,6) = 0
            iR35(5,7) = 0
            iR35(5,8) = 0
            iR35(5,9) = 0
            iR35(5,10) = 0
            iR35(5,11) = 0
            iR35(5,12) = 0
            iR35(5,13) = 0
            iR35(5,14) = 0
            iR35(5,15) = 0
            iR35(5,16) = 0
            iR35(5,17) = 0
            iR35(5,18) = 0
            iR35(5,19) = 0
            iR35(5,20) = 0
            iR35(5,21) = 0
            iR35(5,22) = 0
            iR35(5,23) = 0
            iR35(5,24) = 0
            iR35(5,25) = 0
            iR35(5,26) = 0
            iR35(5,27) = 1
            iR35(5,28) = 0
            iR35(5,29) = 0
            iR35(5,30) = 0
            iR35(5,31) = 0
            iR35(5,32) = 0
            iR35(5,33) = 0
            iR35(5,34) = 0
            iR35(5,35) = 0
            iR35(6,1) = 0
            iR35(6,2) = 0
            iR35(6,3) = 0
            iR35(6,4) = 0
            iR35(6,5) = 0
            iR35(6,6) = 0
            iR35(6,7) = 0
            iR35(6,8) = 0
            iR35(6,9) = 0
            iR35(6,10) = 0
            iR35(6,11) = 0
            iR35(6,12) = 0
            iR35(6,13) = 0
            iR35(6,14) = 0
            iR35(6,15) = 0
            iR35(6,16) = 0
            iR35(6,17) = 0
            iR35(6,18) = 0
            iR35(6,19) = 0
            iR35(6,20) = 0
            iR35(6,21) = 0
            iR35(6,22) = 0
            iR35(6,23) = 0
            iR35(6,24) = 0
            iR35(6,25) = 0
            iR35(6,26) = 1
            iR35(6,27) = 0
            iR35(6,28) = 0
            iR35(6,29) = 0
            iR35(6,30) = 0
            iR35(6,31) = 0
            iR35(6,32) = 0
            iR35(6,33) = 0
            iR35(6,34) = 0
            iR35(6,35) = 0
            iR35(7,1) = g_cov(3,3)*g_contr(1,1)/3
            iR35(7,2) = 2.D0/3.D0*g_cov(3,3)*g_contr(1,2)
            iR35(7,3) = 2.D0/3.D0*g_cov(3,3)*g_contr(1,3)
            iR35(7,4) = g_cov(3,3)*g_contr(2,2)/3
            iR35(7,5) = 2.D0/3.D0*g_cov(3,3)*g_contr(2,3)
            iR35(7,6) = g_cov(3,3)*g_contr(3,3)/3
            iR35(7,7) = 0
            iR35(7,8) = 0
            iR35(7,9) = 0
            iR35(7,10) = 0
            iR35(7,11) = 0
            iR35(7,12) = 0
            iR35(7,13) = 0
            iR35(7,14) = 0
            iR35(7,15) = 0
            iR35(7,16) = 0
            iR35(7,17) = 0
            iR35(7,18) = 0
            iR35(7,19) = 0
            iR35(7,20) = 0
            iR35(7,21) = 0
            iR35(7,22) = 0
            iR35(7,23) = 0
            iR35(7,24) = 0
            iR35(7,25) = 0
            iR35(7,26) = 0
            iR35(7,27) = 0
            iR35(7,28) = 0
            iR35(7,29) = 0
            iR35(7,30) = 0
            iR35(7,31) = 0
            iR35(7,32) = 0
            iR35(7,33) = 0
            iR35(7,34) = 0
            iR35(7,35) = 0
            iR35(8,1) = 0
            iR35(8,2) = 0
            iR35(8,3) = 0
            iR35(8,4) = 0
            iR35(8,5) = 0
            iR35(8,6) = 0
            iR35(8,7) = 0
            iR35(8,8) = 0
            iR35(8,9) = 0
            iR35(8,10) = 0
            iR35(8,11) = 0
            iR35(8,12) = 1
            iR35(8,13) = 0
            iR35(8,14) = 0
            iR35(8,15) = 0
            iR35(8,16) = 0
            iR35(8,17) = 0
            iR35(8,18) = 0
            iR35(8,19) = 0
            iR35(8,20) = 0
            iR35(8,21) = 0
            iR35(8,22) = 0
            iR35(8,23) = 0
            iR35(8,24) = 0
            iR35(8,25) = 0
            iR35(8,26) = 0
            iR35(8,27) = 0
            iR35(8,28) = 0
            iR35(8,29) = 0
            iR35(8,30) = 0
            iR35(8,31) = 0
            iR35(8,32) = 0
            iR35(8,33) = 0
            iR35(8,34) = 0
            iR35(8,35) = 0
            iR35(9,1) = 0
            iR35(9,2) = 0
            iR35(9,3) = 0
            iR35(9,4) = 0
            iR35(9,5) = 0
            iR35(9,6) = 0
            iR35(9,7) = 0
            iR35(9,8) = g_cov(1,3)/g_contr(1,1)/2
            iR35(9,9) = g_cov(2,3)/g_contr(1,1)/2
            iR35(9,10) = g_cov(3,3)/g_contr(1,1)/2
            iR35(9,11) = g_cov(1,3)/3
            iR35(9,12) = g_cov(1,3)*g_contr(1,2)/g_contr(1,1)/3
            iR35(9,13) = -1/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,3)/g_contr(1,1)/3
            iR35(9,14) = 0
            iR35(9,15) = 0
            iR35(9,16) = 0
            iR35(9,17) = 0
            iR35(9,18) = 0
            iR35(9,19) = 0
            iR35(9,20) = 0
            iR35(9,21) = 0
            iR35(9,22) = -g_contr(1,2)/g_contr(1,1)
            iR35(9,23) = 0
            iR35(9,24) = 0
            iR35(9,25) = 0
            iR35(9,26) = 0
            iR35(9,27) = 0
            iR35(9,28) = -g_contr(1,3)/g_contr(1,1)
            iR35(9,29) = 0
            iR35(9,30) = 0
            iR35(9,31) = 0
            iR35(9,32) = 0
            iR35(9,33) = g_cov(1,3)
            iR35(9,34) = g_cov(1,3)*g_contr(1,2)/g_contr(1,1)
            iR35(9,35) = 1/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,3)/g_contr(1,1)
            iR35(10,1) = 0
            iR35(10,2) = 0
            iR35(10,3) = 0
            iR35(10,4) = 0
            iR35(10,5) = 0
            iR35(10,6) = 0
            iR35(10,7) = 0
            iR35(10,8) = 0
            iR35(10,9) = 0
            iR35(10,10) = 0
            iR35(10,11) = 0
            iR35(10,12) = 0
            iR35(10,13) = 0
            iR35(10,14) = 0
            iR35(10,15) = 0
            iR35(10,16) = 0
            iR35(10,17) = 0
            iR35(10,18) = 0
            iR35(10,19) = 0
            iR35(10,20) = 0
            iR35(10,21) = 0
            iR35(10,22) = 0
            iR35(10,23) = 0
            iR35(10,24) = 0
            iR35(10,25) = 0
            iR35(10,26) = 0
            iR35(10,27) = 0
            iR35(10,28) = 0
            iR35(10,29) = 0
            iR35(10,30) = 0
            iR35(10,31) = 0
            iR35(10,32) = 0
            iR35(10,33) = 0
            iR35(10,34) = 0
            iR35(10,35) = 1
            iR35(11,1) = 0
            iR35(11,2) = 0
            iR35(11,3) = 0
            iR35(11,4) = 0
            iR35(11,5) = 0
            iR35(11,6) = 0
            iR35(11,7) = 0
            iR35(11,8) = 0
            iR35(11,9) = 0
            iR35(11,10) = 0
            iR35(11,11) = 0
            iR35(11,12) = 0
            iR35(11,13) = 0
            iR35(11,14) = 0
            iR35(11,15) = 0
            iR35(11,16) = 0
            iR35(11,17) = 0
            iR35(11,18) = 0
            iR35(11,19) = 0
            iR35(11,20) = 0
            iR35(11,21) = 1
            iR35(11,22) = 0
            iR35(11,23) = 0
            iR35(11,24) = 0
            iR35(11,25) = 0
            iR35(11,26) = 0
            iR35(11,27) = 0
            iR35(11,28) = 0
            iR35(11,29) = 0
            iR35(11,30) = 0
            iR35(11,31) = 0
            iR35(11,32) = 0
            iR35(11,33) = 0
            iR35(11,34) = 0
            iR35(11,35) = 0
            iR35(12,1) = 0
            iR35(12,2) = 0
            iR35(12,3) = 0
            iR35(12,4) = 0
            iR35(12,5) = 0
            iR35(12,6) = 0
            iR35(12,7) = 0
            iR35(12,8) = 0
            iR35(12,9) = 0
            iR35(12,10) = 0
            iR35(12,11) = 0
            iR35(12,12) = 0
            iR35(12,13) = 0
            iR35(12,14) = 0
            iR35(12,15) = 0
            iR35(12,16) = 0
            iR35(12,17) = 0
            iR35(12,18) = 0
            iR35(12,19) = 0
            iR35(12,20) = 0
            iR35(12,21) = 0
            iR35(12,22) = 1
            iR35(12,23) = 0
            iR35(12,24) = 0
            iR35(12,25) = 0
            iR35(12,26) = 0
            iR35(12,27) = 0
            iR35(12,28) = 0
            iR35(12,29) = 0
            iR35(12,30) = 0
            iR35(12,31) = 0
            iR35(12,32) = 0
            iR35(12,33) = 0
            iR35(12,34) = 0
            iR35(12,35) = 0
            iR35(13,1) = 0
            iR35(13,2) = 0
            iR35(13,3) = 0
            iR35(13,4) = 0
            iR35(13,5) = 0
            iR35(13,6) = 0
            iR35(13,7) = 0
            iR35(13,8) = 0
            iR35(13,9) = 0
            iR35(13,10) = 0
            iR35(13,11) = 0
            iR35(13,12) = 0
            iR35(13,13) = 0
            iR35(13,14) = 0
            iR35(13,15) = 0
            iR35(13,16) = 0
            iR35(13,17) = 0
            iR35(13,18) = 0
            iR35(13,19) = 0
            iR35(13,20) = 0
            iR35(13,21) = 0
            iR35(13,22) = 0
            iR35(13,23) = 0
            iR35(13,24) = 0
            iR35(13,25) = 0
            iR35(13,26) = 0
            iR35(13,27) = 0
            iR35(13,28) = 0
            iR35(13,29) = 1
            iR35(13,30) = 0
            iR35(13,31) = 0
            iR35(13,32) = 0
            iR35(13,33) = 0
            iR35(13,34) = 0
            iR35(13,35) = 0
            iR35(14,1) = 0
            iR35(14,2) = 0
            iR35(14,3) = 0
            iR35(14,4) = 0
            iR35(14,5) = 0
            iR35(14,6) = 0
            iR35(14,7) = 0
            iR35(14,8) = 0
            iR35(14,9) = 0
            iR35(14,10) = 0
            iR35(14,11) = 0
            iR35(14,12) = 0
            iR35(14,13) = 0
            iR35(14,14) = 0
            iR35(14,15) = 0
            iR35(14,16) = 0
            iR35(14,17) = 0
            iR35(14,18) = 0
            iR35(14,19) = 0
            iR35(14,20) = 0
            iR35(14,21) = 0
            iR35(14,22) = 0
            iR35(14,23) = 0
            iR35(14,24) = 0
            iR35(14,25) = 0
            iR35(14,26) = 0
            iR35(14,27) = 0
            iR35(14,28) = 0
            iR35(14,29) = 0
            iR35(14,30) = 0
            iR35(14,31) = 0
            iR35(14,32) = 0
            iR35(14,33) = 0
            iR35(14,34) = 1
            iR35(14,35) = 0
            iR35(15,1) = 0
            iR35(15,2) = 0
            iR35(15,3) = 0
            iR35(15,4) = 0
            iR35(15,5) = 0
            iR35(15,6) = 0
            iR35(15,7) = 0
            iR35(15,8) = 0
            iR35(15,9) = 0
            iR35(15,10) = 0
            iR35(15,11) = 0
            iR35(15,12) = 0
            iR35(15,13) = 0
            iR35(15,14) = 0
            iR35(15,15) = 0
            iR35(15,16) = 0
            iR35(15,17) = 0
            iR35(15,18) = 0
            iR35(15,19) = 0
            iR35(15,20) = 0
            iR35(15,21) = 0
            iR35(15,22) = 0
            iR35(15,23) = 0
            iR35(15,24) = 0
            iR35(15,25) = 1
            iR35(15,26) = 0
            iR35(15,27) = 0
            iR35(15,28) = 0
            iR35(15,29) = 0
            iR35(15,30) = 0
            iR35(15,31) = 0
            iR35(15,32) = 0
            iR35(15,33) = 0
            iR35(15,34) = 0
            iR35(15,35) = 0
            iR35(16,1) = 0
            iR35(16,2) = 0
            iR35(16,3) = 0
            iR35(16,4) = 0
            iR35(16,5) = 0
            iR35(16,6) = 0
            iR35(16,7) = 0
            iR35(16,8) = 0
            iR35(16,9) = 0
            iR35(16,10) = 0
            iR35(16,11) = 0
            iR35(16,12) = 0
            iR35(16,13) = 0
            iR35(16,14) = 0
            iR35(16,15) = 0
            iR35(16,16) = 0
            iR35(16,17) = 0
            iR35(16,18) = 0
            iR35(16,19) = 0
            iR35(16,20) = 0
            iR35(16,21) = 0
            iR35(16,22) = 0
            iR35(16,23) = 0
            iR35(16,24) = 1
            iR35(16,25) = 0
            iR35(16,26) = 0
            iR35(16,27) = 0
            iR35(16,28) = 0
            iR35(16,29) = 0
            iR35(16,30) = 0
            iR35(16,31) = 0
            iR35(16,32) = 0
            iR35(16,33) = 0
            iR35(16,34) = 0
            iR35(16,35) = 0
            iR35(17,1) = 0
            iR35(17,2) = 0
            iR35(17,3) = 0
            iR35(17,4) = 0
            iR35(17,5) = 0
            iR35(17,6) = 0
            iR35(17,7) = 0
            iR35(17,8) = 0
            iR35(17,9) = 0
            iR35(17,10) = 0
            iR35(17,11) = 0
            iR35(17,12) = 0
            iR35(17,13) = 1
            iR35(17,14) = 0
            iR35(17,15) = 0
            iR35(17,16) = 0
            iR35(17,17) = 0
            iR35(17,18) = 0
            iR35(17,19) = 0
            iR35(17,20) = 0
            iR35(17,21) = 0
            iR35(17,22) = 0
            iR35(17,23) = 0
            iR35(17,24) = 0
            iR35(17,25) = 0
            iR35(17,26) = 0
            iR35(17,27) = 0
            iR35(17,28) = 0
            iR35(17,29) = 0
            iR35(17,30) = 0
            iR35(17,31) = 0
            iR35(17,32) = 0
            iR35(17,33) = 0
            iR35(17,34) = 0
            iR35(17,35) = 0
            iR35(18,1) = 0
            iR35(18,2) = 0
            iR35(18,3) = 0
            iR35(18,4) = 0
            iR35(18,5) = 0
            iR35(18,6) = 0
            iR35(18,7) = 0
            iR35(18,8) = 0
            iR35(18,9) = 0
            iR35(18,10) = 0
            iR35(18,11) = 0
            iR35(18,12) = 0
            iR35(18,13) = 0
            iR35(18,14) = 0
            iR35(18,15) = 0
            iR35(18,16) = 0
            iR35(18,17) = 0
            iR35(18,18) = 0
            iR35(18,19) = 0
            iR35(18,20) = 1
            iR35(18,21) = 0
            iR35(18,22) = 0
            iR35(18,23) = 0
            iR35(18,24) = 0
            iR35(18,25) = 0
            iR35(18,26) = 0
            iR35(18,27) = 0
            iR35(18,28) = 0
            iR35(18,29) = 0
            iR35(18,30) = 0
            iR35(18,31) = 0
            iR35(18,32) = 0
            iR35(18,33) = 0
            iR35(18,34) = 0
            iR35(18,35) = 0
            iR35(19,1) = 0
            iR35(19,2) = 0
            iR35(19,3) = 0
            iR35(19,4) = 0
            iR35(19,5) = 0
            iR35(19,6) = 0
            iR35(19,7) = 0
            iR35(19,8) = 0
            iR35(19,9) = 0
            iR35(19,10) = 0
            iR35(19,11) = 1.D0/3.D0
            iR35(19,12) = g_contr(1,2)/g_contr(1,1)/3
            iR35(19,13) = g_contr(1,3)/g_contr(1,1)/3
            iR35(19,14) = 0
            iR35(19,15) = 0
            iR35(19,16) = 0
            iR35(19,17) = 0
            iR35(19,18) = 0
            iR35(19,19) = 0
            iR35(19,20) = 0
            iR35(19,21) = 0
            iR35(19,22) = 0
            iR35(19,23) = 0
            iR35(19,24) = 0
            iR35(19,25) = 0
            iR35(19,26) = 0
            iR35(19,27) = 0
            iR35(19,28) = 0
            iR35(19,29) = 0
            iR35(19,30) = 0
            iR35(19,31) = 0
            iR35(19,32) = 0
            iR35(19,33) = 1
            iR35(19,34) = 0
            iR35(19,35) = 0
            iR35(20,1) = 0
            iR35(20,2) = 0
            iR35(20,3) = 0
            iR35(20,4) = 0
            iR35(20,5) = 0
            iR35(20,6) = 0
            iR35(20,7) = 0
            iR35(20,8) = 0
            iR35(20,9) = 0
            iR35(20,10) = 0
            iR35(20,11) = 0
            iR35(20,12) = 0
            iR35(20,13) = 0
            iR35(20,14) = 0
            iR35(20,15) = 0
            iR35(20,16) = 0
            iR35(20,17) = 0
            iR35(20,18) = 0
            iR35(20,19) = 0
            iR35(20,20) = 0
            iR35(20,21) = 0
            iR35(20,22) = 0
            iR35(20,23) = 0
            iR35(20,24) = 0
            iR35(20,25) = 0
            iR35(20,26) = 0
            iR35(20,27) = 0
            iR35(20,28) = 0
            iR35(20,29) = 0
            iR35(20,30) = 0
            iR35(20,31) = 1
            iR35(20,32) = 0
            iR35(20,33) = 0
            iR35(20,34) = 0
            iR35(20,35) = 0
            iR35(21,1) = 0
            iR35(21,2) = 0
            iR35(21,3) = 0
            iR35(21,4) = 0
            iR35(21,5) = 0
            iR35(21,6) = 0
            iR35(21,7) = 0
            iR35(21,8) = -(g_cov(1,2)*g_contr(1,2)+g_cov(1,3)*g_contr(1,3))/g_contr(1,1)**2
            iR35(21,9) = g_cov(1,2)/g_contr(1,1)
            iR35(21,10) = g_cov(1,3)/g_contr(1,1)
            iR35(21,11) = -1/g_contr(1,1)+g_cov(1,1)/3
            iR35(21,12) = g_cov(1,1)*g_contr(1,2)/g_contr(1,1)/3
            iR35(21,13) = g_cov(1,1)*g_contr(1,3)/g_contr(1,1)/3
            iR35(21,14) = 1
            iR35(21,15) = 2*g_contr(1,2)/g_contr(1,1)
            iR35(21,16) = 2*g_contr(1,3)/g_contr(1,1)
            iR35(21,17) = g_contr(2,2)/g_contr(1,1)
            iR35(21,18) = 2*g_contr(2,3)/g_contr(1,1)
            iR35(21,19) = g_contr(3,3)/g_contr(1,1)
            iR35(21,20) = 0
            iR35(21,21) = 2*g_contr(1,2)**2/g_contr(1,1)**2
            iR35(21,22) = 2*g_contr(1,2)*g_contr(1,3)/g_contr(1,1)**2
            iR35(21,23) = g_contr(1,2)*g_contr(2,2)/g_contr(1,1)**2
            iR35(21,24) = 2*g_contr(1,2)*g_contr(2,3)/g_contr(1,1)**2
            iR35(21,25) = g_contr(1,2)*g_contr(3,3)/g_contr(1,1)**2
            iR35(21,26) = 0
            iR35(21,27) = 2*g_contr(1,2)*g_contr(1,3)/g_contr(1,1)**2
            iR35(21,28) = 2*g_contr(1,3)**2/g_contr(1,1)**2
            iR35(21,29) = g_contr(1,3)*g_contr(2,2)/g_contr(1,1)**2
            iR35(21,30) = 2*g_contr(1,3)*g_contr(2,3)/g_contr(1,1)**2
            iR35(21,31) = g_contr(1,3)*g_contr(3,3)/g_contr(1,1)**2
            iR35(21,32) = 0
            iR35(21,33) = (-2-g_cov(1,2)*g_contr(1,2)-g_cov(1,3)*g_contr(1,3))/g_contr(1,1)
            iR35(21,34) = -4*g_contr(1,2)/g_contr(1,1)**2+g_cov(1,1)*g_contr(1,2)/g_contr(1,1)
            iR35(21,35) = -4*g_contr(1,3)/g_contr(1,1)**2+g_cov(1,1)*g_contr(1,3)/g_contr(1,1)
            iR35(22,1) = 0
            iR35(22,2) = 0
            iR35(22,3) = 0
            iR35(22,4) = 0
            iR35(22,5) = 0
            iR35(22,6) = 0
            iR35(22,7) = -1/sqrt(g_contr(1,1))/(e**2-1)/e/phi/3
            iR35(22,8) = 1/g_contr(1,1)/(e**2-1)/6
            iR35(22,9) = 0
            iR35(22,10) = 0
            iR35(22,11) = 0
            iR35(22,12) = 0
            iR35(22,13) = 0
            iR35(22,14) = -g_contr(1,1)/(e**2-1)/6
            iR35(22,15) = -g_contr(1,2)/(e**2-1)/3
            iR35(22,16) = -g_contr(1,3)/(e**2-1)/3
            iR35(22,17) = -g_contr(2,2)/(e**2-1)/6
            iR35(22,18) = -g_contr(2,3)/(e**2-1)/3
            iR35(22,19) = -g_contr(3,3)/(e**2-1)/6
            iR35(22,20) = -g_contr(1,2)/(e**2-1)/6
            iR35(22,21) = -g_contr(1,2)**2/g_contr(1,1)/(e**2-1)/3
            iR35(22,22) = -g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(22,23) = -g_contr(1,2)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(22,24) = -g_contr(1,2)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(22,25) = -g_contr(1,2)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(22,26) = -g_contr(1,3)/(e**2-1)/6
            iR35(22,27) = -g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(22,28) = -g_contr(1,3)**2/g_contr(1,1)/(e**2-1)/3
            iR35(22,29) = -g_contr(1,3)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(22,30) = -g_contr(1,3)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(22,31) = -g_contr(1,3)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(22,32) = 0
            iR35(22,33) = 2.D0/3.D0/(e**2-1)
            iR35(22,34) = 2.D0/3.D0*g_contr(1,2)/g_contr(1,1)/(e**2-1)
            iR35(22,35) = 2.D0/3.D0*g_contr(1,3)/g_contr(1,1)/(e**2-1)
            iR35(23,1) = 0
            iR35(23,2) = 0
            iR35(23,3) = 0
            iR35(23,4) = 0
            iR35(23,5) = 0
            iR35(23,6) = 0
            iR35(23,7) = 1/sqrt(g_contr(1,1))/(e**2-1)/e/phi/3
            iR35(23,8) = 1/g_contr(1,1)/(e**2-1)/6
            iR35(23,9) = 0
            iR35(23,10) = 0
            iR35(23,11) = 0
            iR35(23,12) = 0
            iR35(23,13) = 0
            iR35(23,14) = -g_contr(1,1)/(e**2-1)/6
            iR35(23,15) = -g_contr(1,2)/(e**2-1)/3
            iR35(23,16) = -g_contr(1,3)/(e**2-1)/3
            iR35(23,17) = -g_contr(2,2)/(e**2-1)/6
            iR35(23,18) = -g_contr(2,3)/(e**2-1)/3
            iR35(23,19) = -g_contr(3,3)/(e**2-1)/6
            iR35(23,20) = -g_contr(1,2)/(e**2-1)/6
            iR35(23,21) = -g_contr(1,2)**2/g_contr(1,1)/(e**2-1)/3
            iR35(23,22) = -g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(23,23) = -g_contr(1,2)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(23,24) = -g_contr(1,2)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(23,25) = -g_contr(1,2)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(23,26) = -g_contr(1,3)/(e**2-1)/6
            iR35(23,27) = -g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(23,28) = -g_contr(1,3)**2/g_contr(1,1)/(e**2-1)/3
            iR35(23,29) = -g_contr(1,3)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(23,30) = -g_contr(1,3)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(23,31) = -g_contr(1,3)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(23,32) = 0
            iR35(23,33) = 2.D0/3.D0/(e**2-1)
            iR35(23,34) = 2.D0/3.D0*g_contr(1,2)/g_contr(1,1)/(e**2-1)
            iR35(23,35) = 2.D0/3.D0*g_contr(1,3)/g_contr(1,1)/(e**2-1)
            iR35(24,1) = -g_cov(2,3)*sqrt(g_contr(1,1))/phi/6
            iR35(24,2) = -g_cov(2,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(24,3) = -g_cov(2,3)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(24,4) = -g_cov(2,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(24,5) = -(-3+2*g_cov(2,3)*g_contr(2,3))/phi/sqrt(g_contr(1,1))/6
            iR35(24,6) = -g_cov(2,3)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(24,7) = g_cov(2,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(24,8) = -g_cov(2,3)/(e**2-1)/g_contr(1,1)/6
            iR35(24,9) = 0
            iR35(24,10) = 0
            iR35(24,11) = -g_cov(2,3)/6
            iR35(24,12) = -g_cov(2,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(24,13) = -g_cov(2,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(24,14) = g_cov(2,3)/(e**2-1)*g_contr(1,1)/6
            iR35(24,15) = g_cov(2,3)/(e**2-1)*g_contr(1,2)/3
            iR35(24,16) = g_cov(2,3)/(e**2-1)*g_contr(1,3)/3
            iR35(24,17) = g_cov(2,3)/(e**2-1)*g_contr(2,2)/6
            iR35(24,18) = 1.D0/2.D0+g_cov(2,3)/(e**2-1)*g_contr(2,3)/3
            iR35(24,19) = g_cov(2,3)/(e**2-1)*g_contr(3,3)/6
            iR35(24,20) = g_cov(2,3)/(e**2-1)*g_contr(1,2)/6
            iR35(24,21) = g_cov(2,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(24,22) = g_cov(2,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(24,23) = g_cov(2,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(24,24) = g_contr(1,2)/g_contr(1,1)/2+g_cov(2,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(24,25) = g_cov(2,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(24,26) = g_cov(2,3)/(e**2-1)*g_contr(1,3)/6
            iR35(24,27) = g_cov(2,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(24,28) = g_cov(2,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(24,29) = g_cov(2,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(24,30) = g_contr(1,3)/g_contr(1,1)/2+g_cov(2,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(24,31) = g_cov(2,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(24,32) = 0
            iR35(24,33) = -(1+3*e**2)*g_cov(2,3)/(e**2-1)/6
            iR35(24,34) = -(1+3*e**2)*g_cov(2,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(24,35) = -(1+3*e**2)*g_cov(2,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(25,1) = -g_cov(2,2)*sqrt(g_contr(1,1))/phi/6
            iR35(25,2) = -g_cov(2,2)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(25,3) = -g_cov(2,2)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(25,4) = -(-3+g_cov(2,2)*g_contr(2,2))/phi/sqrt(g_contr(1,1))/6
            iR35(25,5) = -g_cov(2,2)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(25,6) = -g_cov(2,2)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(25,7) = g_cov(2,2)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(25,8) = -g_cov(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(25,9) = 0
            iR35(25,10) = 0
            iR35(25,11) = -g_cov(2,2)/6
            iR35(25,12) = -g_cov(2,2)*g_contr(1,2)/g_contr(1,1)/6
            iR35(25,13) = -g_cov(2,2)*g_contr(1,3)/g_contr(1,1)/6
            iR35(25,14) = g_cov(2,2)/(e**2-1)*g_contr(1,1)/6
            iR35(25,15) = g_cov(2,2)/(e**2-1)*g_contr(1,2)/3
            iR35(25,16) = g_cov(2,2)/(e**2-1)*g_contr(1,3)/3
            iR35(25,17) = 1.D0/2.D0+g_cov(2,2)/(e**2-1)*g_contr(2,2)/6
            iR35(25,18) = g_cov(2,2)/(e**2-1)*g_contr(2,3)/3
            iR35(25,19) = g_cov(2,2)/(e**2-1)*g_contr(3,3)/6
            iR35(25,20) = g_cov(2,2)/(e**2-1)*g_contr(1,2)/6
            iR35(25,21) = g_cov(2,2)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(25,22) = g_cov(2,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(25,23) = g_contr(1,2)/g_contr(1,1)/2+g_cov(2,2)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(25,24) = g_cov(2,2)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(25,25) = g_cov(2,2)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(25,26) = g_cov(2,2)/(e**2-1)*g_contr(1,3)/6
            iR35(25,27) = g_cov(2,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(25,28) = g_cov(2,2)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(25,29) = g_contr(1,3)/g_contr(1,1)/2+g_cov(2,2)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(25,30) = g_cov(2,2)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(25,31) = g_cov(2,2)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(25,32) = 0
            iR35(25,33) = -(1+3*e**2)*g_cov(2,2)/(e**2-1)/6
            iR35(25,34) = -(1+3*e**2)*g_cov(2,2)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(25,35) = -(1+3*e**2)*g_cov(2,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(26,1) = -g_cov(1,3)*sqrt(g_contr(1,1))/phi/6
            iR35(26,2) = -g_cov(1,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(26,3) = -(-3+2*g_cov(1,3)*g_contr(1,3))/phi/sqrt(g_contr(1,1))/6
            iR35(26,4) = -g_cov(1,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(26,5) = -g_cov(1,3)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(26,6) = -g_cov(1,3)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(26,7) = g_cov(1,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(26,8) = -g_cov(1,3)*(3*e**2-1)/(e**2-1)/g_contr(1,1)/12
            iR35(26,9) = -g_cov(2,3)/g_contr(1,1)/4
            iR35(26,10) = -g_cov(3,3)/g_contr(1,1)/4
            iR35(26,11) = -g_cov(1,3)/6
            iR35(26,12) = -g_cov(1,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(26,13) = 1/g_contr(1,1)/4-g_cov(1,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(26,14) = g_cov(1,3)/(e**2-1)*g_contr(1,1)/6
            iR35(26,15) = g_cov(1,3)/(e**2-1)*g_contr(1,2)/3
            iR35(26,16) = 1.D0/2.D0+g_cov(1,3)/(e**2-1)*g_contr(1,3)/3
            iR35(26,17) = g_cov(1,3)/(e**2-1)*g_contr(2,2)/6
            iR35(26,18) = g_cov(1,3)/(e**2-1)*g_contr(2,3)/3
            iR35(26,19) = g_cov(1,3)/(e**2-1)*g_contr(3,3)/6
            iR35(26,20) = g_cov(1,3)/(e**2-1)*g_contr(1,2)/6
            iR35(26,21) = g_cov(1,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(26,22) = g_contr(1,2)/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(26,23) = g_cov(1,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(26,24) = g_cov(1,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(26,25) = g_cov(1,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(26,26) = g_cov(1,3)/(e**2-1)*g_contr(1,3)/6
            iR35(26,27) = g_cov(1,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(26,28) = g_contr(1,3)/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(26,29) = g_cov(1,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(26,30) = g_cov(1,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(26,31) = g_cov(1,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(26,32) = 0
            iR35(26,33) = -(1+3*e**2)*g_cov(1,3)/(e**2-1)/6
            iR35(26,34) = -(1+3*e**2)*g_cov(1,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(26,35) = -(1+3*e**2)*g_cov(1,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6-1/g_contr(1,1)/4
            iR35(27,1) = -g_cov(1,2)*sqrt(g_contr(1,1))/phi/6
            iR35(27,2) = -(-3+2*g_cov(1,2)*g_contr(1,2))/phi/sqrt(g_contr(1,1))/6
            iR35(27,3) = -g_cov(1,2)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(27,4) = -g_cov(1,2)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(27,5) = -g_cov(1,2)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(27,6) = -g_cov(1,2)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(27,7) = g_cov(1,2)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(27,8) = -g_cov(1,2)*(3*e**2-1)/(e**2-1)/g_contr(1,1)/12
            iR35(27,9) = -g_cov(2,2)/g_contr(1,1)/4
            iR35(27,10) = -g_cov(2,3)/g_contr(1,1)/4
            iR35(27,11) = -g_cov(1,2)/6
            iR35(27,12) = 1/g_contr(1,1)/4-g_cov(1,2)*g_contr(1,2)/g_contr(1,1)/6
            iR35(27,13) = -g_cov(1,2)*g_contr(1,3)/g_contr(1,1)/6
            iR35(27,14) = g_cov(1,2)/(e**2-1)*g_contr(1,1)/6
            iR35(27,15) = 1.D0/2.D0+g_cov(1,2)/(e**2-1)*g_contr(1,2)/3
            iR35(27,16) = g_cov(1,2)/(e**2-1)*g_contr(1,3)/3
            iR35(27,17) = g_cov(1,2)/(e**2-1)*g_contr(2,2)/6
            iR35(27,18) = g_cov(1,2)/(e**2-1)*g_contr(2,3)/3
            iR35(27,19) = g_cov(1,2)/(e**2-1)*g_contr(3,3)/6
            iR35(27,20) = g_cov(1,2)/(e**2-1)*g_contr(1,2)/6
            iR35(27,21) = g_contr(1,2)/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(27,22) = g_cov(1,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(27,23) = g_cov(1,2)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(27,24) = g_cov(1,2)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(27,25) = g_cov(1,2)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(27,26) = g_cov(1,2)/(e**2-1)*g_contr(1,3)/6
            iR35(27,27) = g_contr(1,3)/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(27,28) = g_cov(1,2)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(27,29) = g_cov(1,2)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(27,30) = g_cov(1,2)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(27,31) = g_cov(1,2)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(27,32) = 0
            iR35(27,33) = -(1+3*e**2)*g_cov(1,2)/(e**2-1)/6
            iR35(27,34) = -(1+3*e**2)*g_cov(1,2)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6-1/g_contr(1,1)/4
            iR35(27,35) = -(1+3*e**2)*g_cov(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(28,1) = 0
            iR35(28,2) = 0
            iR35(28,3) = 0
            iR35(28,4) = 0
            iR35(28,5) = 0
            iR35(28,6) = 0
            iR35(28,7) = 1/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(28,8) = -1/g_contr(1,1)/(e**2-1)/6
            iR35(28,9) = 0
            iR35(28,10) = 0
            iR35(28,11) = -1.D0/6.D0
            iR35(28,12) = -g_contr(1,2)/g_contr(1,1)/6
            iR35(28,13) = -g_contr(1,3)/g_contr(1,1)/6
            iR35(28,14) = g_contr(1,1)/(e**2-1)/6
            iR35(28,15) = g_contr(1,2)/(e**2-1)/3
            iR35(28,16) = g_contr(1,3)/(e**2-1)/3
            iR35(28,17) = g_contr(2,2)/(e**2-1)/6
            iR35(28,18) = g_contr(2,3)/(e**2-1)/3
            iR35(28,19) = g_contr(3,3)/(e**2-1)/6
            iR35(28,20) = g_contr(1,2)/(e**2-1)/6
            iR35(28,21) = g_contr(1,2)**2/g_contr(1,1)/(e**2-1)/3
            iR35(28,22) = g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(28,23) = g_contr(1,2)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(28,24) = g_contr(1,2)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(28,25) = g_contr(1,2)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(28,26) = g_contr(1,3)/(e**2-1)/6
            iR35(28,27) = g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(28,28) = g_contr(1,3)**2/g_contr(1,1)/(e**2-1)/3
            iR35(28,29) = g_contr(1,3)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(28,30) = g_contr(1,3)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(28,31) = g_contr(1,3)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(28,32) = -1/phi/sqrt(g_contr(1,1))/6
            iR35(28,33) = -2.D0/3.D0/(e**2-1)
            iR35(28,34) = -2.D0/3.D0*g_contr(1,2)/g_contr(1,1)/(e**2-1)
            iR35(28,35) = -2.D0/3.D0*g_contr(1,3)/g_contr(1,1)/(e**2-1)
            iR35(29,1) = -g_cov(3,3)*sqrt(g_contr(1,1))/phi/6
            iR35(29,2) = -g_cov(3,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(29,3) = -g_cov(3,3)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(29,4) = -g_cov(3,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(29,5) = -g_cov(3,3)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(29,6) = -(-3+g_cov(3,3)*g_contr(3,3))/phi/sqrt(g_contr(1,1))/6
            iR35(29,7) = g_cov(3,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(29,8) = -g_cov(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(29,9) = 0
            iR35(29,10) = 0
            iR35(29,11) = -g_cov(3,3)/6
            iR35(29,12) = -g_cov(3,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(29,13) = -g_cov(3,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(29,14) = g_cov(3,3)/(e**2-1)*g_contr(1,1)/6
            iR35(29,15) = g_cov(3,3)/(e**2-1)*g_contr(1,2)/3
            iR35(29,16) = g_cov(3,3)/(e**2-1)*g_contr(1,3)/3
            iR35(29,17) = g_cov(3,3)/(e**2-1)*g_contr(2,2)/6
            iR35(29,18) = g_cov(3,3)/(e**2-1)*g_contr(2,3)/3
            iR35(29,19) = 1.D0/2.D0+g_cov(3,3)/(e**2-1)*g_contr(3,3)/6
            iR35(29,20) = g_cov(3,3)/(e**2-1)*g_contr(1,2)/6
            iR35(29,21) = g_cov(3,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(29,22) = g_cov(3,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(29,23) = g_cov(3,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(29,24) = g_cov(3,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(29,25) = g_contr(1,2)/g_contr(1,1)/2+g_cov(3,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(29,26) = g_cov(3,3)/(e**2-1)*g_contr(1,3)/6
            iR35(29,27) = g_cov(3,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(29,28) = g_cov(3,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(29,29) = g_cov(3,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(29,30) = g_cov(3,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(29,31) = g_contr(1,3)/g_contr(1,1)/2+g_cov(3,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(29,32) = 0
            iR35(29,33) = -(1+3*e**2)*g_cov(3,3)/(e**2-1)/6
            iR35(29,34) = -(1+3*e**2)*g_cov(3,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(29,35) = -(1+3*e**2)*g_cov(3,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(30,1) = g_cov(2,3)*sqrt(g_contr(1,1))/phi/6
            iR35(30,2) = g_cov(2,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(30,3) = g_cov(2,3)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(30,4) = g_cov(2,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(30,5) = (-3+2*g_cov(2,3)*g_contr(2,3))/phi/sqrt(g_contr(1,1))/6
            iR35(30,6) = g_cov(2,3)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(30,7) = -g_cov(2,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(30,8) = -g_cov(2,3)/(e**2-1)/g_contr(1,1)/6
            iR35(30,9) = 0
            iR35(30,10) = 0
            iR35(30,11) = -g_cov(2,3)/6
            iR35(30,12) = -g_cov(2,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(30,13) = -g_cov(2,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(30,14) = g_cov(2,3)/(e**2-1)*g_contr(1,1)/6
            iR35(30,15) = g_cov(2,3)/(e**2-1)*g_contr(1,2)/3
            iR35(30,16) = g_cov(2,3)/(e**2-1)*g_contr(1,3)/3
            iR35(30,17) = g_cov(2,3)/(e**2-1)*g_contr(2,2)/6
            iR35(30,18) = 1.D0/2.D0+g_cov(2,3)/(e**2-1)*g_contr(2,3)/3
            iR35(30,19) = g_cov(2,3)/(e**2-1)*g_contr(3,3)/6
            iR35(30,20) = g_cov(2,3)/(e**2-1)*g_contr(1,2)/6
            iR35(30,21) = g_cov(2,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(30,22) = g_cov(2,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(30,23) = g_cov(2,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(30,24) = g_contr(1,2)/g_contr(1,1)/2+g_cov(2,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(30,25) = g_cov(2,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(30,26) = g_cov(2,3)/(e**2-1)*g_contr(1,3)/6
            iR35(30,27) = g_cov(2,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(30,28) = g_cov(2,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(30,29) = g_cov(2,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(30,30) = g_contr(1,3)/g_contr(1,1)/2+g_cov(2,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(30,31) = g_cov(2,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(30,32) = 0
            iR35(30,33) = -(1+3*e**2)*g_cov(2,3)/(e**2-1)/6
            iR35(30,34) = -(1+3*e**2)*g_cov(2,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(30,35) = -(1+3*e**2)*g_cov(2,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(31,1) = g_cov(2,2)*sqrt(g_contr(1,1))/phi/6
            iR35(31,2) = g_cov(2,2)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(31,3) = g_cov(2,2)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(31,4) = (-3+g_cov(2,2)*g_contr(2,2))/phi/sqrt(g_contr(1,1))/6
            iR35(31,5) = g_cov(2,2)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(31,6) = g_cov(2,2)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(31,7) = -g_cov(2,2)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(31,8) = -g_cov(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(31,9) = 0
            iR35(31,10) = 0
            iR35(31,11) = -g_cov(2,2)/6
            iR35(31,12) = -g_cov(2,2)*g_contr(1,2)/g_contr(1,1)/6
            iR35(31,13) = -g_cov(2,2)*g_contr(1,3)/g_contr(1,1)/6
            iR35(31,14) = g_cov(2,2)/(e**2-1)*g_contr(1,1)/6
            iR35(31,15) = g_cov(2,2)/(e**2-1)*g_contr(1,2)/3
            iR35(31,16) = g_cov(2,2)/(e**2-1)*g_contr(1,3)/3
            iR35(31,17) = 1.D0/2.D0+g_cov(2,2)/(e**2-1)*g_contr(2,2)/6
            iR35(31,18) = g_cov(2,2)/(e**2-1)*g_contr(2,3)/3
            iR35(31,19) = g_cov(2,2)/(e**2-1)*g_contr(3,3)/6
            iR35(31,20) = g_cov(2,2)/(e**2-1)*g_contr(1,2)/6
            iR35(31,21) = g_cov(2,2)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(31,22) = g_cov(2,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(31,23) = g_contr(1,2)/g_contr(1,1)/2+g_cov(2,2)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(31,24) = g_cov(2,2)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(31,25) = g_cov(2,2)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(31,26) = g_cov(2,2)/(e**2-1)*g_contr(1,3)/6
            iR35(31,27) = g_cov(2,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(31,28) = g_cov(2,2)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(31,29) = g_contr(1,3)/g_contr(1,1)/2+g_cov(2,2)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(31,30) = g_cov(2,2)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(31,31) = g_cov(2,2)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(31,32) = 0
            iR35(31,33) = -(1+3*e**2)*g_cov(2,2)/(e**2-1)/6
            iR35(31,34) = -(1+3*e**2)*g_cov(2,2)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(31,35) = -(1+3*e**2)*g_cov(2,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(32,1) = g_cov(1,3)*sqrt(g_contr(1,1))/phi/6
            iR35(32,2) = g_cov(1,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(32,3) = (-3+2*g_cov(1,3)*g_contr(1,3))/phi/sqrt(g_contr(1,1))/6
            iR35(32,4) = g_cov(1,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(32,5) = g_cov(1,3)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(32,6) = g_cov(1,3)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(32,7) = -g_cov(1,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(32,8) = -g_cov(1,3)*(3*e**2-1)/(e**2-1)/g_contr(1,1)/12
            iR35(32,9) = -g_cov(2,3)/g_contr(1,1)/4
            iR35(32,10) = -g_cov(3,3)/g_contr(1,1)/4
            iR35(32,11) = -g_cov(1,3)/6
            iR35(32,12) = -g_cov(1,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(32,13) = 1/g_contr(1,1)/4-g_cov(1,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(32,14) = g_cov(1,3)/(e**2-1)*g_contr(1,1)/6
            iR35(32,15) = g_cov(1,3)/(e**2-1)*g_contr(1,2)/3
            iR35(32,16) = 1.D0/2.D0+g_cov(1,3)/(e**2-1)*g_contr(1,3)/3
            iR35(32,17) = g_cov(1,3)/(e**2-1)*g_contr(2,2)/6
            iR35(32,18) = g_cov(1,3)/(e**2-1)*g_contr(2,3)/3
            iR35(32,19) = g_cov(1,3)/(e**2-1)*g_contr(3,3)/6
            iR35(32,20) = g_cov(1,3)/(e**2-1)*g_contr(1,2)/6
            iR35(32,21) = g_cov(1,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(32,22) = g_contr(1,2)/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(32,23) = g_cov(1,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(32,24) = g_cov(1,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(32,25) = g_cov(1,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(32,26) = g_cov(1,3)/(e**2-1)*g_contr(1,3)/6
            iR35(32,27) = g_cov(1,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(32,28) = g_contr(1,3)/g_contr(1,1)/2+g_cov(1,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(32,29) = g_cov(1,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(32,30) = g_cov(1,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(32,31) = g_cov(1,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(32,32) = 0
            iR35(32,33) = -(1+3*e**2)*g_cov(1,3)/(e**2-1)/6
            iR35(32,34) = -(1+3*e**2)*g_cov(1,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(32,35) = -(1+3*e**2)*g_cov(1,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6-1/g_contr(1,1)/4
            iR35(33,1) = g_cov(1,2)*sqrt(g_contr(1,1))/phi/6
            iR35(33,2) = (-3+2*g_cov(1,2)*g_contr(1,2))/phi/sqrt(g_contr(1,1))/6
            iR35(33,3) = g_cov(1,2)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(33,4) = g_cov(1,2)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(33,5) = g_cov(1,2)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(33,6) = g_cov(1,2)*g_contr(3,3)/phi/sqrt(g_contr(1,1))/6
            iR35(33,7) = -g_cov(1,2)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(33,8) = -g_cov(1,2)*(3*e**2-1)/(e**2-1)/g_contr(1,1)/12
            iR35(33,9) = -g_cov(2,2)/g_contr(1,1)/4
            iR35(33,10) = -g_cov(2,3)/g_contr(1,1)/4
            iR35(33,11) = -g_cov(1,2)/6
            iR35(33,12) = 1/g_contr(1,1)/4-g_cov(1,2)*g_contr(1,2)/g_contr(1,1)/6
            iR35(33,13) = -g_cov(1,2)*g_contr(1,3)/g_contr(1,1)/6
            iR35(33,14) = g_cov(1,2)/(e**2-1)*g_contr(1,1)/6
            iR35(33,15) = 1.D0/2.D0+g_cov(1,2)/(e**2-1)*g_contr(1,2)/3
            iR35(33,16) = g_cov(1,2)/(e**2-1)*g_contr(1,3)/3
            iR35(33,17) = g_cov(1,2)/(e**2-1)*g_contr(2,2)/6
            iR35(33,18) = g_cov(1,2)/(e**2-1)*g_contr(2,3)/3
            iR35(33,19) = g_cov(1,2)/(e**2-1)*g_contr(3,3)/6
            iR35(33,20) = g_cov(1,2)/(e**2-1)*g_contr(1,2)/6
            iR35(33,21) = g_contr(1,2)/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(33,22) = g_cov(1,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(33,23) = g_cov(1,2)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(33,24) = g_cov(1,2)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(33,25) = g_cov(1,2)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(33,26) = g_cov(1,2)/(e**2-1)*g_contr(1,3)/6
            iR35(33,27) = g_contr(1,3)/g_contr(1,1)/2+g_cov(1,2)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(33,28) = g_cov(1,2)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(33,29) = g_cov(1,2)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(33,30) = g_cov(1,2)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(33,31) = g_cov(1,2)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(33,32) = 0
            iR35(33,33) = -(1+3*e**2)*g_cov(1,2)/(e**2-1)/6
            iR35(33,34) = -(1+3*e**2)*g_cov(1,2)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6-1/g_contr(1,1)/4
            iR35(33,35) = -(1+3*e**2)*g_cov(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            iR35(34,1) = 0
            iR35(34,2) = 0
            iR35(34,3) = 0
            iR35(34,4) = 0
            iR35(34,5) = 0
            iR35(34,6) = 0
            iR35(34,7) = -1/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(34,8) = -1/g_contr(1,1)/(e**2-1)/6
            iR35(34,9) = 0
            iR35(34,10) = 0
            iR35(34,11) = -1.D0/6.D0
            iR35(34,12) = -g_contr(1,2)/g_contr(1,1)/6
            iR35(34,13) = -g_contr(1,3)/g_contr(1,1)/6
            iR35(34,14) = g_contr(1,1)/(e**2-1)/6
            iR35(34,15) = g_contr(1,2)/(e**2-1)/3
            iR35(34,16) = g_contr(1,3)/(e**2-1)/3
            iR35(34,17) = g_contr(2,2)/(e**2-1)/6
            iR35(34,18) = g_contr(2,3)/(e**2-1)/3
            iR35(34,19) = g_contr(3,3)/(e**2-1)/6
            iR35(34,20) = g_contr(1,2)/(e**2-1)/6
            iR35(34,21) = g_contr(1,2)**2/g_contr(1,1)/(e**2-1)/3
            iR35(34,22) = g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(34,23) = g_contr(1,2)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(34,24) = g_contr(1,2)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(34,25) = g_contr(1,2)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(34,26) = g_contr(1,3)/(e**2-1)/6
            iR35(34,27) = g_contr(1,2)*g_contr(1,3)/g_contr(1,1)/(e**2-1)/3
            iR35(34,28) = g_contr(1,3)**2/g_contr(1,1)/(e**2-1)/3
            iR35(34,29) = g_contr(1,3)*g_contr(2,2)/g_contr(1,1)/(e**2-1)/6
            iR35(34,30) = g_contr(1,3)*g_contr(2,3)/g_contr(1,1)/(e**2-1)/3
            iR35(34,31) = g_contr(1,3)*g_contr(3,3)/g_contr(1,1)/(e**2-1)/6
            iR35(34,32) = 1/phi/sqrt(g_contr(1,1))/6
            iR35(34,33) = -2.D0/3.D0/(e**2-1)
            iR35(34,34) = -2.D0/3.D0*g_contr(1,2)/g_contr(1,1)/(e**2-1)
            iR35(34,35) = -2.D0/3.D0*g_contr(1,3)/g_contr(1,1)/(e**2-1)
            iR35(35,1) = g_cov(3,3)*sqrt(g_contr(1,1))/phi/6
            iR35(35,2) = g_cov(3,3)*g_contr(1,2)/phi/sqrt(g_contr(1,1))/3
            iR35(35,3) = g_cov(3,3)*g_contr(1,3)/phi/sqrt(g_contr(1,1))/3
            iR35(35,4) = g_cov(3,3)*g_contr(2,2)/phi/sqrt(g_contr(1,1))/6
            iR35(35,5) = g_cov(3,3)*g_contr(2,3)/phi/sqrt(g_contr(1,1))/3
            iR35(35,6) = (-3+g_cov(3,3)*g_contr(3,3))/phi/sqrt(g_contr(1,1))/6
            iR35(35,7) = -g_cov(3,3)/phi/(e**2-1)/sqrt(g_contr(1,1))/3
            iR35(35,8) = -g_cov(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(35,9) = 0
            iR35(35,10) = 0
            iR35(35,11) = -g_cov(3,3)/6
            iR35(35,12) = -g_cov(3,3)*g_contr(1,2)/g_contr(1,1)/6
            iR35(35,13) = -g_cov(3,3)*g_contr(1,3)/g_contr(1,1)/6
            iR35(35,14) = g_cov(3,3)/(e**2-1)*g_contr(1,1)/6
            iR35(35,15) = g_cov(3,3)/(e**2-1)*g_contr(1,2)/3
            iR35(35,16) = g_cov(3,3)/(e**2-1)*g_contr(1,3)/3
            iR35(35,17) = g_cov(3,3)/(e**2-1)*g_contr(2,2)/6
            iR35(35,18) = g_cov(3,3)/(e**2-1)*g_contr(2,3)/3
            iR35(35,19) = 1.D0/2.D0+g_cov(3,3)/(e**2-1)*g_contr(3,3)/6
            iR35(35,20) = g_cov(3,3)/(e**2-1)*g_contr(1,2)/6
            iR35(35,21) = g_cov(3,3)*g_contr(1,2)**2/(e**2-1)/g_contr(1,1)/3
            iR35(35,22) = g_cov(3,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(35,23) = g_cov(3,3)*g_contr(1,2)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(35,24) = g_cov(3,3)*g_contr(1,2)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(35,25) = g_contr(1,2)/g_contr(1,1)/2+g_cov(3,3)*g_contr(1,2)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(35,26) = g_cov(3,3)/(e**2-1)*g_contr(1,3)/6
            iR35(35,27) = g_cov(3,3)*g_contr(1,2)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/3
            iR35(35,28) = g_cov(3,3)*g_contr(1,3)**2/(e**2-1)/g_contr(1,1)/3
            iR35(35,29) = g_cov(3,3)*g_contr(1,3)*g_contr(2,2)/(e**2-1)/g_contr(1,1)/6
            iR35(35,30) = g_cov(3,3)*g_contr(1,3)*g_contr(2,3)/(e**2-1)/g_contr(1,1)/3
            iR35(35,31) = g_contr(1,3)/g_contr(1,1)/2+g_cov(3,3)*g_contr(1,3)*g_contr(3,3)/(e**2-1)/g_contr(1,1)/6
            iR35(35,32) = 0
            iR35(35,33) = -(1+3*e**2)*g_cov(3,3)/(e**2-1)/6
            iR35(35,34) = -(1+3*e**2)*g_cov(3,3)*g_contr(1,2)/(e**2-1)/g_contr(1,1)/6
            iR35(35,35) = -(1+3*e**2)*g_cov(3,3)*g_contr(1,3)/(e**2-1)/g_contr(1,1)/6
            !
            temp35 = iR35
            temp35(:,7) = iR35(:,32)
            DO i = 8, 32  
              temp35(:,i) = iR35(:,i-1) 
            ENDDO
            iR35 = temp35 
                        !
            A35 = MATMUL( R35, iR35)
            CONTINUE
        !
        A35 = MATMUL( R35, MATMUL(L35, iR35) )        
        !
        Adiff = ABS( A35 - ATest ) 
        mv = MAXVAL(Adiff)
        ml = MAXLOC(Adiff)
        mir = MAXVAL(ABS(iR35)) 
        !
        LM = 0.0
        R  = 0.0 
        iR = 0.0 
        DO i = 1, nVar
            R(i,i)  = 1.0
            iR(i,i) = 1.0
        ENDDO 
        DO i = 1, 35
            n = dynvar(i)
            LM(n,n) = L35(i,i) 
            DO j = 1, 35
                m = dynvar(j) 
                R(n,m)  = R35(i,j) 
                iR(n,m) = iR35(i,j)
            ENDDO
        ENDDO 
        Afull = MATMUL( R, MATMUL(LM, iR) )    
        dABig = ABS( Afull - An ) 
        mv = MAXVAL(dABig)
        ml = MAXLOC(dABig)
        mir = MAXVAL(ABS(iR))     
        !
        CONTINUE
    ELSE
        !    
        !WHERE(ABS(An)<1e-12) 
        !    An = 0.0
        !ENDWHERE    
        !CALL RG(nVar,nVar,An,L,ImLambda,1,R,itemp,rtemp,ierr)   
        !
        t1 = frac(1.,2.) 
        !
        CALL PDEMatrixA(An,Q,nv,par) 
        !
        dynvar47(:) = (/ 7, 8, 9, 10, 11, 12, 54, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58 /) 
        !
        CALL PDEVarName(VarName)
        DO i = 1, 47
            n = dynvar47(i) 
            Name47(i) = VarName(n) 
        ENDDO    
        !
        e     = EQN%CCZ4e 
        fff   = EQN%CCZ4f 
        mu    = EQN%CCZ4mu 
        phi   = EXP(MAX(-20.,MIN(20.,Q(55))))  
        alpha = EXP(MAX(-20.,MIN(20.,Q(17))))     
        !
        DO i = 1, 47
            DO j = 1, 47 
                n=dynvar47(i) 
                m=dynvar47(j) 
                ATest47(i,j)  = An(n,m) 
            ENDDO
        ENDDO     
        !
        g_cov(1,1) = Q(1)
        g_cov(1,2) = Q(2)
        g_cov(1,3) = Q(3)
        g_cov(2,1) = Q(2)
        g_cov(2,2) = Q(4)
        g_cov(2,3) = Q(5)
        g_cov(3,1) = Q(3)
        g_cov(3,2) = Q(5)
        g_cov(3,3) = Q(6)
        det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
        g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
        g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
        g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
        g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
        g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
        g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
        g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
        g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
        g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det         
        !
        L47 = 0.0 
      L47(1,1) = 0
      L47(2,2) = 0
      L47(3,3) = 0
      L47(4,4) = 0
      L47(5,5) = 0
      L47(6,6) = 0
      L47(7,7) = 0
      L47(8,8) = 0
      L47(9,9) = 0
      L47(10,10) = 0
      L47(11,11) = 0
      L47(12,12) = 0
      L47(13,13) = 0
      L47(14,14) = 0
      L47(15,15) = 0
      L47(16,16) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(17,17) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(18,18) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(19,19) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(20,20) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(21,21) = sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(22,22) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(23,23) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(24,24) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(25,25) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(26,26) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(27,27) = -sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*phi*alpha
      L47(28,28) = sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*phi
      L47(29,29) = -sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*phi
      L47(30,30) = 2.D0/3.D0*sqrt(3.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(31,31) = -2.D0/3.D0*sqrt(3.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(32,32) = sqrt(2.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(33,33) = sqrt(2.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(34,34) = -sqrt(2.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(35,35) = -sqrt(2.D0)*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(36,36) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(37,37) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(38,38) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(39,39) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      L47(40,40) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(41,41) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(42,42) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(43,43) = sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(44,44) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(45,45) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(46,46) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2
      L47(47,47) = -sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*alpha/2             
      !
      beta = (/ Q(18), Q(19), Q(20) /)  
      DO i = 1, 47
          L47(i,i) = L47(i,i) - beta(1) 
      ENDDO      
        !
      R47(1,1) = 0
      R47(1,2) = 0
      R47(1,3) = 0
      R47(1,4) = 0
      R47(1,5) = 0
      R47(1,6) = 0
      R47(1,7) = (13*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-13*g_cov(1,1)*g_cov(2,3)**2-30)/alpha/30
      R47(1,8) = 0
      R47(1,9) = 0
      R47(1,10) = 0
      R47(1,11) = 0
      R47(1,12) = 0
      R47(1,13) = 0
      R47(1,14) = 0
      R47(1,15) = 0
      R47(1,16) = fff*phi*(g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/3
      R47(1,17) = -2*phi*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,18) = -2*phi*(-g_cov(1,1)*g_cov(2,3)+g_cov(1,2)*g_cov(1,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,19) = phi*(-g_cov(1,1)*g_cov(3,3)+g_cov(1,3)**2)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,20) = -phi*(g_cov(1,1)*g_cov(2,2)-g_cov(1,2)**2)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,21) = 2*phi*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,22) = 2*phi*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,23) = 2*phi*(-g_cov(1,1)*g_cov(2,3)+g_cov(1,2)*g_cov(1,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,24) = phi*(g_cov(1,1)*g_cov(2,2)-g_cov(1,2)**2)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,25) = -2*phi*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,26) = -phi*(-g_cov(1,1)*g_cov(3,3)+g_cov(1,3)**2)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,27) = -fff*phi*(g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/3
      R47(1,28) = -sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(2*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*fff-3*g_cov(2,2)*g_cov(3,3)*phi**2*alpha*g_cov(1,1)-2*g_cov(1,1)*g_cov(2,3)**2*fff+3*g_cov(1,1)*g_cov(2,3)**2*phi**2*alpha-6*fff+9*phi**2*alpha)/alpha/phi/3
      R47(1,29) = sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(2*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*fff-3*g_cov(2,2)*g_cov(3,3)*phi**2*alpha*g_cov(1,1)-2*g_cov(1,1)*g_cov(2,3)**2*fff+3*g_cov(1,1)*g_cov(2,3)**2*phi**2*alpha-6*fff+9*phi**2*alpha)/alpha/phi/3
      R47(1,30) = 0
      R47(1,31) = 0
      R47(1,32) = 3*mu*phi**2*sqrt(2.D0)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,33) = -3*mu*phi**2*sqrt(2.D0)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,34) = 3*mu*phi**2*sqrt(2.D0)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,35) = -3*mu*phi**2*sqrt(2.D0)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(1,36) = 0
      R47(1,37) = 0
      R47(1,38) = 0
      R47(1,39) = 0
      R47(1,40) = 2.D0/3.D0*(16*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*fff-3*g_cov(2,2)*g_cov(3,3)**2*mu*alpha**2*g_cov(1,1)+9*mu*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-48*fff*g_cov(1,2)**2*g_cov(3,3)**2-18*mu*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-16*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*fff+96*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+3*g_cov(3,3)*mu*alpha**2*g_cov(1,1)*g_cov(2,3)**2-48*g_cov(1,3)**2*g_cov(2,3)**2*fff+9*mu*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(1,41) = -4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*g_cov(1,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(1,42) = -4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*g_cov(1,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(1,43) = 4*mu*phi**2*(g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-2*g_cov(1,2)**2*g_cov(3,3)+2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-1)/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(1,44) = -4*mu*phi**2*(g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-2*g_cov(1,2)**2*g_cov(3,3)+2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-1)/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(1,45) = -2.D0/3.D0*(16*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*fff-3*g_cov(2,2)*g_cov(3,3)**2*mu*alpha**2*g_cov(1,1)+9*mu*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-48*fff*g_cov(1,2)**2*g_cov(3,3)**2-18*mu*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-16*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*fff+96*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+3*g_cov(3,3)*mu*alpha**2*g_cov(1,1)*g_cov(2,3)**2-48*g_cov(1,3)**2*g_cov(2,3)**2*fff+9*mu*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(1,46) = 4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*g_cov(1,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(1,47) = 4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*g_cov(1,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,1) = 0
      R47(2,2) = 0
      R47(2,3) = 0
      R47(2,4) = 0
      R47(2,5) = 0
      R47(2,6) = 0
      R47(2,7) = 13.D0/30.D0*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,2)/alpha
      R47(2,8) = 0
      R47(2,9) = 0
      R47(2,10) = 0
      R47(2,11) = 0
      R47(2,12) = 0
      R47(2,13) = 0
      R47(2,14) = 0
      R47(2,15) = 0
      R47(2,16) = phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,2)/3
      R47(2,17) = phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(2,18) = 0
      R47(2,19) = 0
      R47(2,20) = 0
      R47(2,21) = 0
      R47(2,22) = -phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(2,23) = 0
      R47(2,24) = 0
      R47(2,25) = 0
      R47(2,26) = 0
      R47(2,27) = -phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,2)/3
      R47(2,28) = -(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*g_cov(1,2)*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(2,29) = (-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*g_cov(1,2)*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(2,30) = 0
      R47(2,31) = 0
      R47(2,32) = -3.D0/2.D0*mu*phi**2*sqrt(2.D0)/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,33) = 0
      R47(2,34) = 0
      R47(2,35) = 3.D0/2.D0*mu*phi**2*sqrt(2.D0)/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,36) = 0
      R47(2,37) = 0
      R47(2,38) = 0
      R47(2,39) = 0
      R47(2,40) = -2.D0/3.D0*phi**2*(-6*g_cov(2,2)*alpha**2*g_cov(1,2)*mu*g_cov(3,3)**2+32*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)*fff+9*g_cov(1,3)*g_cov(2,2)*alpha**2*mu*g_cov(2,3)*g_cov(3,3)-48*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,3)*g_cov(2,3)+6*alpha**2*g_cov(1,2)*mu*g_cov(3,3)*g_cov(2,3)**2-32*g_cov(3,3)*g_cov(1,2)*fff*g_cov(2,3)**2-9*g_cov(1,3)*alpha**2*mu*g_cov(2,3)**3+48*g_cov(1,3)*g_cov(2,3)**3*fff)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(2,41) = 2*(g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-2*g_cov(1,3)*g_cov(2,3)**2+g_cov(3,3)*g_cov(1,2)*g_cov(2,3))*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,42) = 2*(g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-2*g_cov(1,3)*g_cov(2,3)**2+g_cov(3,3)*g_cov(1,2)*g_cov(2,3))*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,43) = 4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*g_cov(2,2)*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,44) = -4*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*g_cov(2,2)*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,45) = 2.D0/3.D0*phi**2*(-6*g_cov(2,2)*alpha**2*g_cov(1,2)*mu*g_cov(3,3)**2+32*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)*fff+9*g_cov(1,3)*g_cov(2,2)*alpha**2*mu*g_cov(2,3)*g_cov(3,3)-48*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,3)*g_cov(2,3)+6*alpha**2*g_cov(1,2)*mu*g_cov(3,3)*g_cov(2,3)**2-32*g_cov(3,3)*g_cov(1,2)*fff*g_cov(2,3)**2-9*g_cov(1,3)*alpha**2*mu*g_cov(2,3)**3+48*g_cov(1,3)*g_cov(2,3)**3*fff)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(2,46) = -2*(g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-2*g_cov(1,3)*g_cov(2,3)**2+g_cov(3,3)*g_cov(1,2)*g_cov(2,3))*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(2,47) = -2*(g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-2*g_cov(1,3)*g_cov(2,3)**2+g_cov(3,3)*g_cov(1,2)*g_cov(2,3))*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,1) = 0
      R47(3,2) = 0
      R47(3,3) = 0
      R47(3,4) = 0
      R47(3,5) = 0
      R47(3,6) = 0
      R47(3,7) = 13.D0/30.D0*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,3)/alpha
      R47(3,8) = 0
      R47(3,9) = 0
      R47(3,10) = 0
      R47(3,11) = 0
      R47(3,12) = 0
      R47(3,13) = 0
      R47(3,14) = 0
      R47(3,15) = 0
      R47(3,16) = phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,3)/3
      R47(3,17) = 0
      R47(3,18) = 0
      R47(3,19) = 0
      R47(3,20) = 0
      R47(3,21) = phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(3,22) = 0
      R47(3,23) = 0
      R47(3,24) = 0
      R47(3,25) = -phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(3,26) = 0
      R47(3,27) = -phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(1,3)/3
      R47(3,28) = -sqrt(2.D0)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*g_cov(1,3)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(3,29) = sqrt(2.D0)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*g_cov(1,3)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(3,30) = 0
      R47(3,31) = 0
      R47(3,32) = 0
      R47(3,33) = -3.D0/2.D0*mu*phi**2*sqrt(2.D0)/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,34) = 3.D0/2.D0*mu*phi**2*sqrt(2.D0)/(mu-2*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,35) = 0
      R47(3,36) = 0
      R47(3,37) = 0
      R47(3,38) = 0
      R47(3,39) = 0
      R47(3,40) = -2.D0/3.D0*phi**2*g_cov(1,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(3,41) = -2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,42) = -2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,43) = 4*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,44) = -4*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,45) = 2.D0/3.D0*phi**2*g_cov(1,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(3,46) = 2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(3,47) = 2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,1) = 0
      R47(4,2) = 0
      R47(4,3) = 0
      R47(4,4) = 0
      R47(4,5) = 0
      R47(4,6) = 0
      R47(4,7) = 13.D0/30.D0*g_cov(2,2)*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/alpha
      R47(4,8) = 0
      R47(4,9) = 0
      R47(4,10) = 0
      R47(4,11) = 0
      R47(4,12) = 0
      R47(4,13) = 0
      R47(4,14) = 0
      R47(4,15) = 0
      R47(4,16) = phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,2)/3
      R47(4,17) = 0
      R47(4,18) = 0
      R47(4,19) = phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(4,20) = 0
      R47(4,21) = 0
      R47(4,22) = 0
      R47(4,23) = 0
      R47(4,24) = 0
      R47(4,25) = 0
      R47(4,26) = -phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(4,27) = -phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,2)/3
      R47(4,28) = -g_cov(2,2)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(4,29) = g_cov(2,2)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(4,30) = 0
      R47(4,31) = 0
      R47(4,32) = 0
      R47(4,33) = 0
      R47(4,34) = 0
      R47(4,35) = 0
      R47(4,36) = 0
      R47(4,37) = 0
      R47(4,38) = 0
      R47(4,39) = 0
      R47(4,40) = 2.D0/3.D0*phi**2*(6*mu*g_cov(2,2)**2*g_cov(3,3)**2*alpha**2-32*fff*g_cov(3,3)**2*g_cov(2,2)**2-15*mu*g_cov(2,2)*g_cov(3,3)*g_cov(2,3)**2*alpha**2+80*fff*g_cov(2,2)*g_cov(2,3)**2*g_cov(3,3)+9*mu*alpha**2*g_cov(2,3)**4-48*fff*g_cov(2,3)**4)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/g_cov(3,3)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(4,41) = 4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*g_cov(2,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,42) = 4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*g_cov(2,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,43) = -4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,2)*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,44) = 4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,2)*mu*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,45) = -2.D0/3.D0*phi**2*(6*mu*g_cov(2,2)**2*g_cov(3,3)**2*alpha**2-32*fff*g_cov(3,3)**2*g_cov(2,2)**2-15*mu*g_cov(2,2)*g_cov(3,3)*g_cov(2,3)**2*alpha**2+80*fff*g_cov(2,2)*g_cov(2,3)**2*g_cov(3,3)+9*mu*alpha**2*g_cov(2,3)**4-48*fff*g_cov(2,3)**4)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/g_cov(3,3)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(4,46) = -4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*g_cov(2,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(4,47) = -4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*g_cov(2,3)*phi**2/g_cov(3,3)/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(5,1) = 0
      R47(5,2) = 0
      R47(5,3) = 0
      R47(5,4) = 0
      R47(5,5) = 0
      R47(5,6) = 0
      R47(5,7) = 13.D0/30.D0*g_cov(2,3)*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/alpha
      R47(5,8) = 0
      R47(5,9) = 0
      R47(5,10) = 0
      R47(5,11) = 0
      R47(5,12) = 0
      R47(5,13) = 0
      R47(5,14) = 0
      R47(5,15) = 0
      R47(5,16) = phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,3)/3
      R47(5,17) = 0
      R47(5,18) = phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(5,19) = 0
      R47(5,20) = 0
      R47(5,21) = 0
      R47(5,22) = 0
      R47(5,23) = -phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(5,24) = 0
      R47(5,25) = 0
      R47(5,26) = 0
      R47(5,27) = -phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(2,3)/3
      R47(5,28) = -g_cov(2,3)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(5,29) = g_cov(2,3)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))*sqrt(2.D0)/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/alpha/phi/3
      R47(5,30) = 0
      R47(5,31) = 0
      R47(5,32) = 0
      R47(5,33) = 0
      R47(5,34) = 0
      R47(5,35) = 0
      R47(5,36) = 0
      R47(5,37) = 0
      R47(5,38) = 0
      R47(5,39) = 0
      R47(5,40) = -2.D0/3.D0*phi**2*g_cov(2,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(5,41) = 2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(5,42) = 2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(5,43) = 0
      R47(5,44) = 0
      R47(5,45) = 2.D0/3.D0*phi**2*g_cov(2,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(5,46) = -2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(5,47) = -2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(6,1) = 0
      R47(6,2) = 0
      R47(6,3) = 0
      R47(6,4) = 0
      R47(6,5) = 0
      R47(6,6) = 0
      R47(6,7) = 13.D0/30.D0*g_cov(3,3)*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/alpha
      R47(6,8) = 0
      R47(6,9) = 0
      R47(6,10) = 0
      R47(6,11) = 0
      R47(6,12) = 0
      R47(6,13) = 0
      R47(6,14) = 0
      R47(6,15) = 0
      R47(6,16) = phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(3,3)/3
      R47(6,17) = 0
      R47(6,18) = 0
      R47(6,19) = 0
      R47(6,20) = phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(6,21) = 0
      R47(6,22) = 0
      R47(6,23) = 0
      R47(6,24) = -phi*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(6,25) = 0
      R47(6,26) = 0
      R47(6,27) = -phi*fff*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*g_cov(3,3)/3
      R47(6,28) = -sqrt(2.D0)*g_cov(3,3)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha/3
      R47(6,29) = sqrt(2.D0)*g_cov(3,3)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha/3
      R47(6,30) = 0
      R47(6,31) = 0
      R47(6,32) = 0
      R47(6,33) = 0
      R47(6,34) = 0
      R47(6,35) = 0
      R47(6,36) = 0
      R47(6,37) = 0
      R47(6,38) = 0
      R47(6,39) = 0
      R47(6,40) = -2.D0/3.D0*phi**2*g_cov(3,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(6,41) = 0
      R47(6,42) = 0
      R47(6,43) = 4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(6,44) = -4*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu*phi**2/(mu-4*phi**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(6,45) = 2.D0/3.D0*phi**2*g_cov(3,3)*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(6,46) = 0
      R47(6,47) = 0
      R47(7,1) = 0
      R47(7,2) = 0
      R47(7,3) = 0
      R47(7,4) = 0
      R47(7,5) = 0
      R47(7,6) = 0
      R47(7,7) = -13.D0/10.D0*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/alpha
      R47(7,8) = 0
      R47(7,9) = 0
      R47(7,10) = 0
      R47(7,11) = 0
      R47(7,12) = 0
      R47(7,13) = 0
      R47(7,14) = 0
      R47(7,15) = 0
      R47(7,16) = -(4*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-4*fff*g_cov(2,3)**2+3*alpha**2*phi**2*g_cov(2,3)**2)*phi/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(7,17) = 0
      R47(7,18) = 0
      R47(7,19) = 0
      R47(7,20) = 0
      R47(7,21) = 0
      R47(7,22) = 0
      R47(7,23) = 0
      R47(7,24) = 0
      R47(7,25) = 0
      R47(7,26) = 0
      R47(7,27) = (4*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-4*fff*g_cov(2,3)**2+3*alpha**2*phi**2*g_cov(2,3)**2)*phi/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(7,28) = sqrt(2.D0)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha
      R47(7,29) = -sqrt(2.D0)*(-3*g_cov(2,2)*phi**2*alpha*g_cov(3,3)+3*phi**2*alpha*g_cov(2,3)**2-2*fff*g_cov(2,3)**2+2*g_cov(2,2)*fff*g_cov(3,3))/sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha
      R47(7,30) = 0
      R47(7,31) = 0
      R47(7,32) = 0
      R47(7,33) = 0
      R47(7,34) = 0
      R47(7,35) = 0
      R47(7,36) = 0
      R47(7,37) = 0
      R47(7,38) = 0
      R47(7,39) = 0
      R47(7,40) = 2*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(7,41) = 0
      R47(7,42) = 0
      R47(7,43) = 0
      R47(7,44) = 0
      R47(7,45) = -2*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(7,46) = 0
      R47(7,47) = 0
      R47(8,1) = 0
      R47(8,2) = 0
      R47(8,3) = 0
      R47(8,4) = 0
      R47(8,5) = 0
      R47(8,6) = 0
      R47(8,7) = -(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/alpha/2
      R47(8,8) = 0
      R47(8,9) = 0
      R47(8,10) = 0
      R47(8,11) = 0
      R47(8,12) = 0
      R47(8,13) = 0
      R47(8,14) = 0
      R47(8,15) = 0
      R47(8,16) = -(4*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-4*fff*g_cov(2,3)**2+3*alpha**2*phi**2*g_cov(2,3)**2)*phi/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(8,17) = 0
      R47(8,18) = 0
      R47(8,19) = 0
      R47(8,20) = 0
      R47(8,21) = 0
      R47(8,22) = 0
      R47(8,23) = 0
      R47(8,24) = 0
      R47(8,25) = 0
      R47(8,26) = 0
      R47(8,27) = (4*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-4*fff*g_cov(2,3)**2+3*alpha**2*phi**2*g_cov(2,3)**2)*phi/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(8,28) = 0
      R47(8,29) = 0
      R47(8,30) = 0
      R47(8,31) = 0
      R47(8,32) = 0
      R47(8,33) = 0
      R47(8,34) = 0
      R47(8,35) = 0
      R47(8,36) = 0
      R47(8,37) = 0
      R47(8,38) = 0
      R47(8,39) = 0
      R47(8,40) = (3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(8,41) = 0
      R47(8,42) = 0
      R47(8,43) = 0
      R47(8,44) = 0
      R47(8,45) = -(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-16*g_cov(2,2)*fff*g_cov(3,3)-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2)*phi**2/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(8,46) = 0
      R47(8,47) = 0
      R47(9,1) = -(-g_cov(1,2)**3*g_cov(3,3)**3+3*g_cov(1,2)**2*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,3)-3*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)**2*g_cov(2,3)**2+g_cov(1,3)**3*g_cov(2,3)**3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,2) = g_cov(3,3)**2*g_cov(1,1)*g_cov(2,2)**2-2*g_cov(2,2)*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)-3*g_cov(2,2)*g_cov(3,3)+g_cov(2,3)**4*g_cov(1,1)+3*g_cov(2,3)**2
      R47(9,3) = 2*g_cov(1,3)*g_cov(2,2)-g_cov(1,3)*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2+g_cov(1,3)*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)+g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)-2*g_cov(1,2)*g_cov(2,3)
      R47(9,4) = 2*(g_cov(1,3)**2*g_cov(2,3)**3*g_cov(1,2)+g_cov(1,2)**2*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)-g_cov(1,1)*g_cov(1,3)*g_cov(2,3)**4-g_cov(1,3)*g_cov(2,3)**2-3*g_cov(1,2)**2*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)**2-2*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*g_cov(2,2)+2*g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3+g_cov(1,2)**3*g_cov(3,3)**2*g_cov(2,3)+2*g_cov(3,3)*g_cov(1,2)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,5) = g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,3)**3*g_cov(1,1)-2*g_cov(1,3)*g_cov(2,3)-g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)+g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)+2*g_cov(1,2)*g_cov(3,3)
      R47(9,6) = 0
      R47(9,7) = 0
      R47(9,8) = -2*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)+2*g_cov(3,3)*g_cov(2,3)**2*g_cov(1,2)+2*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(2,3)-2*g_cov(1,3)*g_cov(2,3)**3
      R47(9,9) = -2*(-g_cov(1,1)*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**3+g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*g_cov(1,3)*g_cov(2,2)**2+g_cov(1,3)*g_cov(2,3)**3*g_cov(1,2)**2-g_cov(1,3)*g_cov(2,3)*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)-g_cov(1,3)*g_cov(2,2)*g_cov(2,3)-g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)**2+g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,2)*g_cov(2,3)**2+g_cov(1,2)**3*g_cov(2,2)*g_cov(3,3)**2+g_cov(1,2)*g_cov(3,3)*g_cov(2,2)-g_cov(1,2)**3*g_cov(3,3)*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,10) = 2*g_cov(2,2)*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)-2*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)-2*g_cov(2,3)**3*g_cov(1,2)+2*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)
      R47(9,11) = 0
      R47(9,12) = (g_cov(1,1)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)**3-g_cov(1,1)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)**2-g_cov(1,2)**2*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)**2-g_cov(1,3)*g_cov(2,2)**2+g_cov(1,3)*g_cov(2,2)*g_cov(1,2)**2*g_cov(2,3)**2-g_cov(1,1)*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)*g_cov(2,2)**2+g_cov(1,2)*g_cov(1,1)*g_cov(2,2)*g_cov(2,3)**3+g_cov(1,2)**3*g_cov(2,2)*g_cov(3,3)*g_cov(2,3)+g_cov(1,2)*g_cov(2,3)*g_cov(2,2)-g_cov(1,2)**3*g_cov(2,3)**3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,13) = g_cov(3,3)**2*g_cov(2,2)**2-2*g_cov(2,2)*g_cov(2,3)**2*g_cov(3,3)+g_cov(2,3)**4
      R47(9,14) = (g_cov(1,3)**2*g_cov(2,3)**3*g_cov(1,2)+g_cov(1,2)**2*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)-g_cov(1,1)*g_cov(1,3)*g_cov(2,3)**4-g_cov(1,3)*g_cov(2,3)**2-3*g_cov(1,2)**2*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)**2-2*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*g_cov(2,2)+2*g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3+g_cov(1,2)**3*g_cov(3,3)**2*g_cov(2,3)+2*g_cov(3,3)*g_cov(1,2)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,15) = -(-g_cov(1,1)*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**3+g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*g_cov(1,3)*g_cov(2,2)**2+g_cov(1,3)*g_cov(2,3)**3*g_cov(1,2)**2-g_cov(1,3)*g_cov(2,3)*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)-g_cov(1,3)*g_cov(2,2)*g_cov(2,3)-g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)**2+g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,2)*g_cov(2,3)**2+g_cov(1,2)**3*g_cov(2,2)*g_cov(3,3)**2+g_cov(1,2)*g_cov(3,3)*g_cov(2,2)-g_cov(1,2)**3*g_cov(3,3)*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,16) = alpha**2*phi**2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,17) = 0
      R47(9,18) = 0
      R47(9,19) = 0
      R47(9,20) = 0
      R47(9,21) = 0
      R47(9,22) = 0
      R47(9,23) = 0
      R47(9,24) = 0
      R47(9,25) = 0
      R47(9,26) = 0
      R47(9,27) = alpha**2*phi**2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(9,28) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(9,29) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(9,30) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(9,31) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(9,32) = 0
      R47(9,33) = 0
      R47(9,34) = 0
      R47(9,35) = 0
      R47(9,36) = 0
      R47(9,37) = 0
      R47(9,38) = 0
      R47(9,39) = 0
      R47(9,40) = -alpha**2*(7*g_cov(2,2)*g_cov(3,3)*mu-7*g_cov(2,3)**2*mu+32*g_cov(2,3)**2*phi**2-32*g_cov(2,2)*g_cov(3,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(9,41) = 0
      R47(9,42) = 0
      R47(9,43) = 0
      R47(9,44) = 0
      R47(9,45) = -alpha**2*(7*g_cov(2,2)*g_cov(3,3)*mu-7*g_cov(2,3)**2*mu+32*g_cov(2,3)**2*phi**2-32*g_cov(2,2)*g_cov(3,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(9,46) = 0
      R47(9,47) = 0
      R47(10,1) = -(g_cov(1,3)**4*g_cov(2,3)**4-4*g_cov(1,2)*g_cov(3,3)*g_cov(2,3)**3*g_cov(1,3)**3+6*g_cov(1,2)**2*g_cov(3,3)**2*g_cov(1,3)**2*g_cov(2,3)**2+2*g_cov(1,3)**2*g_cov(2,3)**2*g_cov(3,3)-4*g_cov(1,2)**3*g_cov(3,3)**3*g_cov(1,3)*g_cov(2,3)-4*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*g_cov(3,3)**2+g_cov(1,2)**4*g_cov(3,3)**4+2*g_cov(1,2)**2*g_cov(3,3)**3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,2) = g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,3)**3*g_cov(1,1)-g_cov(1,3)*g_cov(2,3)-g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)+g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)+g_cov(1,2)*g_cov(3,3)
      R47(10,3) = (g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,2)*g_cov(1,3)*g_cov(1,1)-2*g_cov(2,2)*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(3,3)+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**4-g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)+2*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)**2*g_cov(2,3)**3+g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)-g_cov(1,1)**2*g_cov(2,3)**5-g_cov(2,3)**3*g_cov(1,1)+g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(10,4) = 2*(g_cov(1,3)**3*g_cov(1,2)*g_cov(2,3)**4-g_cov(1,3)**2*g_cov(2,3)**5*g_cov(1,1)-4*g_cov(1,2)**2*g_cov(3,3)*g_cov(2,3)**3*g_cov(1,3)**2-g_cov(1,2)*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)**2-g_cov(1,2)**3*g_cov(3,3)**3*g_cov(1,3)*g_cov(2,2)-3*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**2+5*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)**4+6*g_cov(1,2)**3*g_cov(3,3)**2*g_cov(2,3)**2*g_cov(1,3)+2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(3,3)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**3+3*g_cov(1,1)*g_cov(3,3)**3*g_cov(1,2)**2*g_cov(2,3)*g_cov(2,2)+g_cov(2,3)*g_cov(1,1)*g_cov(3,3)**2*g_cov(2,2)-g_cov(1,1)**2*g_cov(3,3)*g_cov(2,3)**5-4*g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)**2*g_cov(2,3)**3-2*g_cov(2,3)**3*g_cov(1,1)*g_cov(3,3)-2*g_cov(1,2)**4*g_cov(3,3)**3*g_cov(2,3)-2*g_cov(1,2)**2*g_cov(3,3)**2*g_cov(2,3)-g_cov(2,3)*g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,5) = -(g_cov(2,3)**4*g_cov(1,1)*g_cov(1,3)**2+2*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-4*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3*g_cov(1,2)*g_cov(1,3)-g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)**2-g_cov(2,2)*g_cov(1,2)**2*g_cov(1,1)*g_cov(3,3)**3+g_cov(3,3)*g_cov(2,3)**4*g_cov(1,1)**2+2*g_cov(3,3)**2*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2+g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)+g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(10,6) = -g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(10,7) = 0
      R47(10,8) = 2*g_cov(1,3)**2*g_cov(2,3)**2-4*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+2*g_cov(1,2)**2*g_cov(3,3)**2+2*g_cov(3,3)
      R47(10,9) = -2*(g_cov(2,3)**4*g_cov(1,2)**2*g_cov(1,3)**2+2*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)**3*g_cov(1,3)*g_cov(2,3)-2*g_cov(2,3)**5*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)-2*g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-g_cov(1,3)*g_cov(2,3)**3*g_cov(1,2)+g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3*g_cov(1,2)*g_cov(1,3)-4*g_cov(2,3)**3*g_cov(1,2)**3*g_cov(3,3)*g_cov(1,3)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)**2*g_cov(2,3)**2+g_cov(2,2)**2*g_cov(3,3)**3*g_cov(1,1)*g_cov(1,2)**2+g_cov(3,3)**2*g_cov(1,1)*g_cov(2,2)**2-2*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)**2-2*g_cov(2,2)*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)-g_cov(2,2)*g_cov(3,3)-g_cov(2,2)*g_cov(3,3)**3*g_cov(1,2)**4-2*g_cov(3,3)*g_cov(2,3)**4*g_cov(1,1)**2*g_cov(2,2)-3*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2+g_cov(1,1)**2*g_cov(2,3)**6+2*g_cov(2,3)**4*g_cov(1,1)*g_cov(1,2)**2*g_cov(3,3)+g_cov(2,3)**4*g_cov(1,1)+2*g_cov(2,3)**2*g_cov(1,2)**4*g_cov(3,3)**2+3*g_cov(2,3)**2*g_cov(1,2)**2*g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,10) = -2*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)+2*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+2*g_cov(2,3)**3*g_cov(1,1)-2*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)
      R47(10,11) = g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(10,12) = -(2*g_cov(1,2)**3*g_cov(3,3)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**4+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**3*g_cov(3,3)**2-2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**2*g_cov(2,3)**2*g_cov(3,3)-g_cov(1,2)**3*g_cov(1,3)*g_cov(3,3)**2*g_cov(2,2)**2+g_cov(1,2)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)-g_cov(1,2)*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)-g_cov(1,3)*g_cov(2,3)**4*g_cov(1,2)**3-g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)*g_cov(2,2)**3+2*g_cov(1,1)**2*g_cov(3,3)*g_cov(2,2)**2*g_cov(2,3)**3+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)*g_cov(3,3)**2*g_cov(2,2)**2+g_cov(1,2)**2*g_cov(3,3)*g_cov(2,2)*g_cov(2,3)+g_cov(2,3)*g_cov(2,2)-2*g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**3*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**5+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**5-g_cov(2,3)**3*g_cov(1,2)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,13) = -g_cov(3,3)**2*g_cov(1,2)*g_cov(2,2)+g_cov(3,3)*g_cov(2,3)**2*g_cov(1,2)+g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(2,3)-g_cov(1,3)*g_cov(2,3)**3
      R47(10,14) = (g_cov(1,3)**3*g_cov(1,2)*g_cov(2,3)**4-g_cov(1,3)**2*g_cov(2,3)**3-g_cov(1,3)**2*g_cov(2,3)**5*g_cov(1,1)-4*g_cov(1,2)**2*g_cov(3,3)*g_cov(2,3)**3*g_cov(1,3)**2-2*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)**2-g_cov(1,2)**3*g_cov(3,3)**3*g_cov(1,3)*g_cov(2,2)-3*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**2+5*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)**4+6*g_cov(1,2)**3*g_cov(3,3)**2*g_cov(2,3)**2*g_cov(1,3)+5*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(3,3)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**3+3*g_cov(1,1)*g_cov(3,3)**3*g_cov(1,2)**2*g_cov(2,3)*g_cov(2,2)+2*g_cov(2,3)*g_cov(1,1)*g_cov(3,3)**2*g_cov(2,2)-g_cov(1,1)**2*g_cov(3,3)*g_cov(2,3)**5-4*g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)**2*g_cov(2,3)**3-3*g_cov(2,3)**3*g_cov(1,1)*g_cov(3,3)-2*g_cov(1,2)**4*g_cov(3,3)**3*g_cov(2,3)-3*g_cov(1,2)**2*g_cov(3,3)**2*g_cov(2,3)-2*g_cov(2,3)*g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,15) = -(g_cov(2,3)**4*g_cov(1,2)**2*g_cov(1,3)**2-2*g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3*g_cov(1,2)*g_cov(1,3)+2*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)**3*g_cov(1,3)*g_cov(2,3)-2*g_cov(2,3)**5*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)-4*g_cov(2,3)**3*g_cov(1,2)**3*g_cov(3,3)*g_cov(1,3)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)**2*g_cov(2,3)**2+g_cov(2,2)**2*g_cov(3,3)**3*g_cov(1,1)*g_cov(1,2)**2-2*g_cov(3,3)*g_cov(2,3)**4*g_cov(1,1)**2*g_cov(2,2)-3*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2-g_cov(2,2)*g_cov(3,3)**3*g_cov(1,2)**4-g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)**2+g_cov(1,1)**2*g_cov(2,3)**6+2*g_cov(2,3)**4*g_cov(1,1)*g_cov(1,2)**2*g_cov(3,3)+2*g_cov(2,3)**2*g_cov(1,2)**4*g_cov(3,3)**2+2*g_cov(2,3)**2*g_cov(1,2)**2*g_cov(3,3)-g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(10,16) = alpha**2*phi**2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(10,17) = 0
      R47(10,18) = 0
      R47(10,19) = 0
      R47(10,20) = 0
      R47(10,21) = 0
      R47(10,22) = 0
      R47(10,23) = 0
      R47(10,24) = 0
      R47(10,25) = 0
      R47(10,26) = 0
      R47(10,27) = alpha**2*phi**2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(10,28) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(10,29) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(10,30) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(10,31) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(10,32) = -g_cov(3,3)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(10,33) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(10,34) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(10,35) = -g_cov(3,3)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(10,36) = -2*g_cov(2,3)
      R47(10,37) = 2*g_cov(3,3)
      R47(10,38) = -2*g_cov(2,3)
      R47(10,39) = 2*g_cov(3,3)
      R47(10,40) = -alpha**2*(7*g_cov(1,3)*g_cov(2,3)*mu-7*g_cov(1,2)*g_cov(3,3)*mu+32*g_cov(1,2)*g_cov(3,3)*phi**2-32*g_cov(1,3)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(10,41) = 0
      R47(10,42) = 0
      R47(10,43) = 0
      R47(10,44) = 0
      R47(10,45) = -alpha**2*(7*g_cov(1,3)*g_cov(2,3)*mu-7*g_cov(1,2)*g_cov(3,3)*mu+32*g_cov(1,2)*g_cov(3,3)*phi**2-32*g_cov(1,3)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(10,46) = 0
      R47(10,47) = 0
      R47(11,1) = (g_cov(1,3)**3*g_cov(1,2)*g_cov(2,3)**4-g_cov(1,3)**2*g_cov(2,3)**5*g_cov(1,1)+g_cov(1,3)**2*g_cov(2,3)**3-4*g_cov(1,2)**2*g_cov(3,3)*g_cov(2,3)**3*g_cov(1,3)**2-g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(3,3)+5*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)**4-3*g_cov(1,1)*g_cov(1,2)*g_cov(3,3)**2*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**2-g_cov(1,2)**3*g_cov(3,3)**3*g_cov(1,3)*g_cov(2,2)+6*g_cov(1,2)**3*g_cov(3,3)**2*g_cov(2,3)**2*g_cov(1,3)+3*g_cov(1,1)*g_cov(3,3)**3*g_cov(1,2)**2*g_cov(2,3)*g_cov(2,2)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**3-g_cov(1,1)**2*g_cov(3,3)*g_cov(2,3)**5-4*g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)**2*g_cov(2,3)**3-g_cov(2,3)**3*g_cov(1,1)*g_cov(3,3)-2*g_cov(1,2)**4*g_cov(3,3)**3*g_cov(2,3)-g_cov(1,2)**2*g_cov(3,3)**2*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,2) = g_cov(1,3)*g_cov(2,2)-g_cov(1,3)*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2+g_cov(1,3)*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)+g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)-g_cov(1,2)*g_cov(2,3)
      R47(11,3) = (g_cov(2,2)**3*g_cov(1,1)**2*g_cov(3,3)**2-2*g_cov(2,2)**2*g_cov(1,1)**2*g_cov(3,3)*g_cov(2,3)**2-g_cov(2,2)**2*g_cov(1,1)*g_cov(3,3)**2*g_cov(1,2)**2-g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2-g_cov(2,2)+g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)+2*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2+g_cov(2,2)*g_cov(2,3)**4*g_cov(1,1)**2-g_cov(1,1)*g_cov(1,2)**2*g_cov(2,3)**4)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(11,4) = -2*(g_cov(2,3)**4*g_cov(1,2)**2*g_cov(1,3)**2-2*g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3*g_cov(1,2)*g_cov(1,3)+2*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)**3*g_cov(1,3)*g_cov(2,3)-g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)-2*g_cov(2,3)**5*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)-4*g_cov(2,3)**3*g_cov(1,2)**3*g_cov(3,3)*g_cov(1,3)+g_cov(1,3)*g_cov(2,3)**3*g_cov(1,2)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)**2*g_cov(2,3)**2+g_cov(2,2)**2*g_cov(3,3)**3*g_cov(1,1)*g_cov(1,2)**2-2*g_cov(3,3)*g_cov(2,3)**4*g_cov(1,1)**2*g_cov(2,2)-3*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2-g_cov(2,2)*g_cov(3,3)**3*g_cov(1,2)**4+g_cov(1,1)**2*g_cov(2,3)**6+2*g_cov(2,3)**4*g_cov(1,1)*g_cov(1,2)**2*g_cov(3,3)+2*g_cov(2,3)**2*g_cov(1,2)**4*g_cov(3,3)**2+g_cov(2,3)**2*g_cov(1,2)**2*g_cov(3,3)-g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,5) = (g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,2)*g_cov(1,3)*g_cov(1,1)-2*g_cov(2,2)*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(3,3)+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)**4-g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)+2*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)**2*g_cov(2,3)**3+g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)-g_cov(1,1)**2*g_cov(2,3)**5-g_cov(2,3)**3*g_cov(1,1)+g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(11,6) = g_cov(2,2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(11,7) = 0
      R47(11,8) = -2*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)+2*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+2*g_cov(2,3)**3*g_cov(1,1)-2*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)
      R47(11,9) = -2*(2*g_cov(1,2)**3*g_cov(3,3)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**4+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**3*g_cov(3,3)**2-2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**2*g_cov(2,3)**2*g_cov(3,3)-g_cov(1,2)**3*g_cov(1,3)*g_cov(3,3)**2*g_cov(2,2)**2-g_cov(1,3)*g_cov(2,3)**4*g_cov(1,2)**3-g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)*g_cov(2,2)**3+2*g_cov(1,1)**2*g_cov(3,3)*g_cov(2,2)**2*g_cov(2,3)**3+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)*g_cov(3,3)**2*g_cov(2,2)**2+g_cov(2,3)*g_cov(2,2)-2*g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**3*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**5+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**5)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,10) = 2*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2-2*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)-2*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)+2*g_cov(1,2)**2*g_cov(2,3)**2
      R47(11,11) = -g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(11,12) = -(g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)**4-2*g_cov(1,1)**2*g_cov(3,3)*g_cov(2,2)**3*g_cov(2,3)**2-2*g_cov(1,1)*g_cov(1,2)**2*g_cov(2,2)**3*g_cov(3,3)**2-g_cov(2,2)**2+4*g_cov(1,1)*g_cov(2,2)**2*g_cov(1,2)**2*g_cov(2,3)**2*g_cov(3,3)+g_cov(1,2)**4*g_cov(3,3)**2*g_cov(2,2)**2+g_cov(1,1)**2*g_cov(2,2)**2*g_cov(2,3)**4-2*g_cov(1,1)*g_cov(2,2)*g_cov(1,2)**2*g_cov(2,3)**4-2*g_cov(1,2)**4*g_cov(2,2)*g_cov(2,3)**2*g_cov(3,3)+g_cov(1,2)**4*g_cov(2,3)**4)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,13) = g_cov(2,2)*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)-g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)-g_cov(2,3)**3*g_cov(1,2)+g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)
      R47(11,14) = -(g_cov(2,3)**4*g_cov(1,2)**2*g_cov(1,3)**2-2*g_cov(2,2)**2*g_cov(3,3)**2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**3*g_cov(1,2)*g_cov(1,3)+2*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,2)**3*g_cov(1,3)*g_cov(2,3)-2*g_cov(2,3)**5*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)-4*g_cov(2,3)**3*g_cov(1,2)**3*g_cov(3,3)*g_cov(1,3)+g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,2)**2*g_cov(2,3)**2+g_cov(2,2)**2*g_cov(3,3)**3*g_cov(1,1)*g_cov(1,2)**2-2*g_cov(3,3)*g_cov(2,3)**4*g_cov(1,1)**2*g_cov(2,2)-3*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*g_cov(2,3)**2*g_cov(1,2)**2-g_cov(2,2)*g_cov(3,3)**3*g_cov(1,2)**4-g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)**2+g_cov(1,1)**2*g_cov(2,3)**6+2*g_cov(2,3)**4*g_cov(1,1)*g_cov(1,2)**2*g_cov(3,3)+2*g_cov(2,3)**2*g_cov(1,2)**4*g_cov(3,3)**2+2*g_cov(2,3)**2*g_cov(1,2)**2*g_cov(3,3)-g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,15) = -(2*g_cov(1,2)**3*g_cov(3,3)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)*g_cov(2,3)**4+g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**3*g_cov(3,3)**2-2*g_cov(1,1)*g_cov(1,2)*g_cov(1,3)*g_cov(2,2)**2*g_cov(2,3)**2*g_cov(3,3)-g_cov(1,2)**3*g_cov(1,3)*g_cov(3,3)**2*g_cov(2,2)**2-g_cov(1,2)*g_cov(2,3)**2*g_cov(1,3)*g_cov(2,2)+g_cov(1,2)*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)-g_cov(1,3)*g_cov(2,3)**4*g_cov(1,2)**3-g_cov(3,3)**2*g_cov(1,1)**2*g_cov(2,3)*g_cov(2,2)**3+2*g_cov(1,1)**2*g_cov(3,3)*g_cov(2,2)**2*g_cov(2,3)**3+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)*g_cov(3,3)**2*g_cov(2,2)**2-g_cov(1,2)**2*g_cov(3,3)*g_cov(2,2)*g_cov(2,3)+g_cov(2,3)*g_cov(2,2)-2*g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**3*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)**2*g_cov(2,2)*g_cov(2,3)**5+g_cov(1,2)**2*g_cov(1,1)*g_cov(2,3)**5+g_cov(2,3)**3*g_cov(1,2)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(11,16) = -alpha**2*phi**2*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(11,17) = 0
      R47(11,18) = 0
      R47(11,19) = 0
      R47(11,20) = 0
      R47(11,21) = 0
      R47(11,22) = 0
      R47(11,23) = 0
      R47(11,24) = 0
      R47(11,25) = 0
      R47(11,26) = 0
      R47(11,27) = -alpha**2*phi**2*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(11,28) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(11,29) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(11,30) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(11,31) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(11,32) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(11,33) = -g_cov(2,2)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(11,34) = -g_cov(2,2)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(11,35) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(11,36) = 2*g_cov(2,2)
      R47(11,37) = -2*g_cov(2,3)
      R47(11,38) = 2*g_cov(2,2)
      R47(11,39) = -2*g_cov(2,3)
      R47(11,40) = alpha**2*(7*g_cov(1,3)*g_cov(2,2)*mu-7*g_cov(1,2)*g_cov(2,3)*mu-32*g_cov(1,3)*g_cov(2,2)*phi**2+32*g_cov(1,2)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(11,41) = 0
      R47(11,42) = 0
      R47(11,43) = 0
      R47(11,44) = 0
      R47(11,45) = alpha**2*(7*g_cov(1,3)*g_cov(2,2)*mu-7*g_cov(1,2)*g_cov(2,3)*mu-32*g_cov(1,3)*g_cov(2,2)*phi**2+32*g_cov(1,2)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(11,46) = 0
      R47(11,47) = 0
      R47(12,1) = 0
      R47(12,2) = 0
      R47(12,3) = 0
      R47(12,4) = 0
      R47(12,5) = 0
      R47(12,6) = 0
      R47(12,7) = 0
      R47(12,8) = 0
      R47(12,9) = 0
      R47(12,10) = 0
      R47(12,11) = 0
      R47(12,12) = 0
      R47(12,13) = 0
      R47(12,14) = 0
      R47(12,15) = 0
      R47(12,16) = alpha**2*phi**2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(12,17) = 0
      R47(12,18) = 0
      R47(12,19) = 0
      R47(12,20) = 0
      R47(12,21) = 0
      R47(12,22) = 0
      R47(12,23) = 0
      R47(12,24) = 0
      R47(12,25) = 0
      R47(12,26) = 0
      R47(12,27) = alpha**2*phi**2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(12,28) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(12,29) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(12,30) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(12,31) = -4*g_cov(2,2)*g_cov(3,3)+4*g_cov(2,3)**2
      R47(12,32) = 0
      R47(12,33) = 0
      R47(12,34) = 0
      R47(12,35) = 0
      R47(12,36) = 0
      R47(12,37) = 0
      R47(12,38) = 0
      R47(12,39) = 0
      R47(12,40) = -alpha**2*(7*g_cov(2,2)*g_cov(3,3)*mu-7*g_cov(2,3)**2*mu+32*g_cov(2,3)**2*phi**2-32*g_cov(2,2)*g_cov(3,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(12,41) = 0
      R47(12,42) = 0
      R47(12,43) = 0
      R47(12,44) = 0
      R47(12,45) = -alpha**2*(7*g_cov(2,2)*g_cov(3,3)*mu-7*g_cov(2,3)**2*mu+32*g_cov(2,3)**2*phi**2-32*g_cov(2,2)*g_cov(3,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(12,46) = 0
      R47(12,47) = 0
      R47(13,1) = g_cov(3,3)**2*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,2) = 0
      R47(13,3) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/fff
      R47(13,4) = -2*g_cov(2,3)*g_cov(3,3)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,5) = -g_cov(3,3)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/fff
      R47(13,6) = 0
      R47(13,7) = 0
      R47(13,8) = 0
      R47(13,9) = (g_cov(2,2)*g_cov(3,3)+g_cov(2,3)**2)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,10) = 0
      R47(13,11) = 0
      R47(13,12) = -g_cov(2,3)*g_cov(2,2)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,13) = 0
      R47(13,14) = -g_cov(2,3)*g_cov(3,3)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,15) = g_cov(2,3)**2*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(13,16) = alpha**2*phi**2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(13,17) = 0
      R47(13,18) = 0
      R47(13,19) = 0
      R47(13,20) = 0
      R47(13,21) = 0
      R47(13,22) = 0
      R47(13,23) = 0
      R47(13,24) = 0
      R47(13,25) = 0
      R47(13,26) = 0
      R47(13,27) = alpha**2*phi**2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(13,28) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(13,29) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(13,30) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(13,31) = 4*g_cov(1,2)*g_cov(3,3)-4*g_cov(1,3)*g_cov(2,3)
      R47(13,32) = -g_cov(3,3)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(13,33) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(13,34) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(13,35) = -g_cov(3,3)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(13,36) = -2*g_cov(2,3)
      R47(13,37) = 2*g_cov(3,3)
      R47(13,38) = -2*g_cov(2,3)
      R47(13,39) = 2*g_cov(3,3)
      R47(13,40) = -alpha**2*(7*g_cov(1,3)*g_cov(2,3)*mu-7*g_cov(1,2)*g_cov(3,3)*mu+32*g_cov(1,2)*g_cov(3,3)*phi**2-32*g_cov(1,3)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(13,41) = 0
      R47(13,42) = 0
      R47(13,43) = 0
      R47(13,44) = 0
      R47(13,45) = -alpha**2*(7*g_cov(1,3)*g_cov(2,3)*mu-7*g_cov(1,2)*g_cov(3,3)*mu+32*g_cov(1,2)*g_cov(3,3)*phi**2-32*g_cov(1,3)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(13,46) = 0
      R47(13,47) = 0
      R47(14,1) = -g_cov(2,3)*g_cov(3,3)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,2) = 0
      R47(14,3) = -g_cov(2,2)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/fff
      R47(14,4) = (g_cov(2,2)*g_cov(3,3)+g_cov(2,3)**2)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,5) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/fff
      R47(14,6) = 0
      R47(14,7) = 0
      R47(14,8) = 0
      R47(14,9) = -2*g_cov(2,3)*g_cov(2,2)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,10) = 0
      R47(14,11) = 0
      R47(14,12) = g_cov(2,2)**2*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,13) = 0
      R47(14,14) = g_cov(2,3)**2*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,15) = -g_cov(2,3)*g_cov(2,2)*alpha**2*mu/fff/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(14,16) = -alpha**2*phi**2*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(14,17) = 0
      R47(14,18) = 0
      R47(14,19) = 0
      R47(14,20) = 0
      R47(14,21) = 0
      R47(14,22) = 0
      R47(14,23) = 0
      R47(14,24) = 0
      R47(14,25) = 0
      R47(14,26) = 0
      R47(14,27) = -alpha**2*phi**2*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(14,28) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(14,29) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(14,30) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(14,31) = -4*g_cov(1,2)*g_cov(2,3)+4*g_cov(1,3)*g_cov(2,2)
      R47(14,32) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(14,33) = -g_cov(2,2)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(14,34) = -g_cov(2,2)*alpha**2*mu/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(14,35) = g_cov(2,3)*mu*alpha**2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(14,36) = 2*g_cov(2,2)
      R47(14,37) = -2*g_cov(2,3)
      R47(14,38) = 2*g_cov(2,2)
      R47(14,39) = -2*g_cov(2,3)
      R47(14,40) = alpha**2*(7*g_cov(1,3)*g_cov(2,2)*mu-7*g_cov(1,2)*g_cov(2,3)*mu-32*g_cov(1,3)*g_cov(2,2)*phi**2+32*g_cov(1,2)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(14,41) = 0
      R47(14,42) = 0
      R47(14,43) = 0
      R47(14,44) = 0
      R47(14,45) = alpha**2*(7*g_cov(1,3)*g_cov(2,2)*mu-7*g_cov(1,2)*g_cov(2,3)*mu-32*g_cov(1,3)*g_cov(2,2)*phi**2+32*g_cov(1,2)*g_cov(2,3)*phi**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(14,46) = 0
      R47(14,47) = 0
      R47(15,1) = 0
      R47(15,2) = 0
      R47(15,3) = 0
      R47(15,4) = 0
      R47(15,5) = 0
      R47(15,6) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(15,7) = 0
      R47(15,8) = 0
      R47(15,9) = 0
      R47(15,10) = 0
      R47(15,11) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(15,12) = 0
      R47(15,13) = 0
      R47(15,14) = 0
      R47(15,15) = 0
      R47(15,16) = 0
      R47(15,17) = 0
      R47(15,18) = 0
      R47(15,19) = 0
      R47(15,20) = 0
      R47(15,21) = 0
      R47(15,22) = 0
      R47(15,23) = 0
      R47(15,24) = 0
      R47(15,25) = 0
      R47(15,26) = 0
      R47(15,27) = 0
      R47(15,28) = 2*(-3*phi**2*alpha+2*fff)/alpha**2/phi**2
      R47(15,29) = 2*(-3*phi**2*alpha+2*fff)/alpha**2/phi**2
      R47(15,30) = 0
      R47(15,31) = 0
      R47(15,32) = 0
      R47(15,33) = 0
      R47(15,34) = 0
      R47(15,35) = 0
      R47(15,36) = 0
      R47(15,37) = 0
      R47(15,38) = 0
      R47(15,39) = 0
      R47(15,40) = 0
      R47(15,41) = 0
      R47(15,42) = 0
      R47(15,43) = 0
      R47(15,44) = 0
      R47(15,45) = 0
      R47(15,46) = 0
      R47(15,47) = 0
      R47(16,1) = 0
      R47(16,2) = 0
      R47(16,3) = 0
      R47(16,4) = 0
      R47(16,5) = 0
      R47(16,6) = 0
      R47(16,7) = 0
      R47(16,8) = 0
      R47(16,9) = 0
      R47(16,10) = 0
      R47(16,11) = 1
      R47(16,12) = 0
      R47(16,13) = 0
      R47(16,14) = 0
      R47(16,15) = 0
      R47(16,16) = 0
      R47(16,17) = 0
      R47(16,18) = 0
      R47(16,19) = 0
      R47(16,20) = 0
      R47(16,21) = 0
      R47(16,22) = 0
      R47(16,23) = 0
      R47(16,24) = 0
      R47(16,25) = 0
      R47(16,26) = 0
      R47(16,27) = 0
      R47(16,28) = 0
      R47(16,29) = 0
      R47(16,30) = 0
      R47(16,31) = 0
      R47(16,32) = 0
      R47(16,33) = 0
      R47(16,34) = 0
      R47(16,35) = 0
      R47(16,36) = 0
      R47(16,37) = 0
      R47(16,38) = 0
      R47(16,39) = 0
      R47(16,40) = 0
      R47(16,41) = 0
      R47(16,42) = 0
      R47(16,43) = 0
      R47(16,44) = 0
      R47(16,45) = 0
      R47(16,46) = 0
      R47(16,47) = 0
      R47(17,1) = 0
      R47(17,2) = 0
      R47(17,3) = 0
      R47(17,4) = 0
      R47(17,5) = 0
      R47(17,6) = 1
      R47(17,7) = 0
      R47(17,8) = 0
      R47(17,9) = 0
      R47(17,10) = 0
      R47(17,11) = 0
      R47(17,12) = 0
      R47(17,13) = 0
      R47(17,14) = 0
      R47(17,15) = 0
      R47(17,16) = 0
      R47(17,17) = 0
      R47(17,18) = 0
      R47(17,19) = 0
      R47(17,20) = 0
      R47(17,21) = 0
      R47(17,22) = 0
      R47(17,23) = 0
      R47(17,24) = 0
      R47(17,25) = 0
      R47(17,26) = 0
      R47(17,27) = 0
      R47(17,28) = 0
      R47(17,29) = 0
      R47(17,30) = 0
      R47(17,31) = 0
      R47(17,32) = 0
      R47(17,33) = 0
      R47(17,34) = 0
      R47(17,35) = 0
      R47(17,36) = 0
      R47(17,37) = 0
      R47(17,38) = 0
      R47(17,39) = 0
      R47(17,40) = 0
      R47(17,41) = 0
      R47(17,42) = 0
      R47(17,43) = 0
      R47(17,44) = 0
      R47(17,45) = 0
      R47(17,46) = 0
      R47(17,47) = 0
      R47(18,1) = 0
      R47(18,2) = 0
      R47(18,3) = 0
      R47(18,4) = 0
      R47(18,5) = 0
      R47(18,6) = 0
      R47(18,7) = -(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(18,8) = 0
      R47(18,9) = 0
      R47(18,10) = 0
      R47(18,11) = 0
      R47(18,12) = 0
      R47(18,13) = 0
      R47(18,14) = 0
      R47(18,15) = 0
      R47(18,16) = -phi*alpha*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff
      R47(18,17) = 0
      R47(18,18) = 0
      R47(18,19) = 0
      R47(18,20) = 0
      R47(18,21) = 0
      R47(18,22) = 0
      R47(18,23) = 0
      R47(18,24) = 0
      R47(18,25) = 0
      R47(18,26) = 0
      R47(18,27) = phi*alpha*sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff
      R47(18,28) = 2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha
      R47(18,29) = -2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))/phi/alpha
      R47(18,30) = 2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(18,31) = -2*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(18,32) = -3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(18,33) = 3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(18,34) = -3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(18,35) = 3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(18,36) = 0
      R47(18,37) = 0
      R47(18,38) = 0
      R47(18,39) = 0
      R47(18,40) = mu*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-3*g_cov(2,3)**2*alpha**2*mu-12*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2+12*alpha**2*phi**2*g_cov(2,3)**2-2*g_cov(2,2)*fff*g_cov(3,3)+2*fff*g_cov(2,3)**2)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(18,41) = 0
      R47(18,42) = 0
      R47(18,43) = 0
      R47(18,44) = 0
      R47(18,45) = -mu*(3*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-3*g_cov(2,3)**2*alpha**2*mu-12*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2+12*alpha**2*phi**2*g_cov(2,3)**2-2*g_cov(2,2)*fff*g_cov(3,3)+2*fff*g_cov(2,3)**2)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(18,46) = 0
      R47(18,47) = 0
      R47(19,1) = 0
      R47(19,2) = 0
      R47(19,3) = 0
      R47(19,4) = 0
      R47(19,5) = 0
      R47(19,6) = 0
      R47(19,7) = 0
      R47(19,8) = 0
      R47(19,9) = 0
      R47(19,10) = 0
      R47(19,11) = 0
      R47(19,12) = 0
      R47(19,13) = 0
      R47(19,14) = 0
      R47(19,15) = 0
      R47(19,16) = 0
      R47(19,17) = 0
      R47(19,18) = 0
      R47(19,19) = 0
      R47(19,20) = 0
      R47(19,21) = 0
      R47(19,22) = 0
      R47(19,23) = 0
      R47(19,24) = 0
      R47(19,25) = 0
      R47(19,26) = 0
      R47(19,27) = 0
      R47(19,28) = 0
      R47(19,29) = 0
      R47(19,30) = 0
      R47(19,31) = 0
      R47(19,32) = 3*mu*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(19,33) = 0
      R47(19,34) = 0
      R47(19,35) = -3*mu*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(19,36) = 0
      R47(19,37) = 0
      R47(19,38) = 0
      R47(19,39) = 0
      R47(19,40) = 0
      R47(19,41) = 0
      R47(19,42) = 0
      R47(19,43) = 0
      R47(19,44) = 0
      R47(19,45) = 0
      R47(19,46) = 0
      R47(19,47) = 0
      R47(20,1) = 0
      R47(20,2) = 0
      R47(20,3) = 0
      R47(20,4) = 0
      R47(20,5) = 0
      R47(20,6) = 0
      R47(20,7) = 0
      R47(20,8) = 0
      R47(20,9) = 0
      R47(20,10) = 0
      R47(20,11) = 0
      R47(20,12) = 0
      R47(20,13) = 0
      R47(20,14) = 0
      R47(20,15) = 0
      R47(20,16) = 0
      R47(20,17) = 0
      R47(20,18) = 0
      R47(20,19) = 0
      R47(20,20) = 0
      R47(20,21) = 0
      R47(20,22) = 0
      R47(20,23) = 0
      R47(20,24) = 0
      R47(20,25) = 0
      R47(20,26) = 0
      R47(20,27) = 0
      R47(20,28) = 0
      R47(20,29) = 0
      R47(20,30) = 0
      R47(20,31) = 0
      R47(20,32) = 0
      R47(20,33) = 3*mu*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(20,34) = -3*mu*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(20,35) = 0
      R47(20,36) = 0
      R47(20,37) = 0
      R47(20,38) = 0
      R47(20,39) = 0
      R47(20,40) = 0
      R47(20,41) = 0
      R47(20,42) = 0
      R47(20,43) = 0
      R47(20,44) = 0
      R47(20,45) = 0
      R47(20,46) = 0
      R47(20,47) = 0
      R47(21,1) = 0
      R47(21,2) = 0
      R47(21,3) = 0
      R47(21,4) = 0
      R47(21,5) = 0
      R47(21,6) = 0
      R47(21,7) = -(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(21,8) = 0
      R47(21,9) = 0
      R47(21,10) = 0
      R47(21,11) = 0
      R47(21,12) = 0
      R47(21,13) = 0
      R47(21,14) = 0
      R47(21,15) = 0
      R47(21,16) = -phi*alpha*fff*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,17) = 0
      R47(21,18) = 0
      R47(21,19) = 0
      R47(21,20) = 0
      R47(21,21) = 0
      R47(21,22) = 0
      R47(21,23) = 0
      R47(21,24) = 0
      R47(21,25) = 0
      R47(21,26) = 0
      R47(21,27) = phi*alpha*fff*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,28) = 2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/phi/alpha/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,29) = -2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/phi/alpha/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,30) = 2*fff*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(21,31) = -2*fff*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(21,32) = -alpha*mu*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*(-6*g_cov(1,3)**2*g_cov(2,3)**2*fff+3*mu*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2+12*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)-6*mu*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-6*fff*g_cov(1,2)**2*g_cov(3,3)**2+3*mu*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-fff*g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(21,33) = -alpha*mu*sqrt(2.D0)*(6*fff*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)-3*mu*g_cov(2,3)**2*alpha**2*g_cov(1,2)*g_cov(1,3)+3*mu*g_cov(2,2)*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)-6*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+6*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*fff-3*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*alpha**2*mu+3*g_cov(2,3)**3*g_cov(1,1)*alpha**2*mu-6*g_cov(2,3)**3*g_cov(1,1)*fff-5*fff*g_cov(2,3)+3*g_cov(2,3)*mu*alpha**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(21,34) = alpha*mu*sqrt(2.D0)*(6*fff*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)-3*mu*g_cov(2,3)**2*alpha**2*g_cov(1,2)*g_cov(1,3)+3*mu*g_cov(2,2)*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)-6*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+6*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*fff-3*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*alpha**2*mu+3*g_cov(2,3)**3*g_cov(1,1)*alpha**2*mu-6*g_cov(2,3)**3*g_cov(1,1)*fff-5*fff*g_cov(2,3)+3*g_cov(2,3)*mu*alpha**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(21,35) = alpha*mu*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*(-6*g_cov(1,3)**2*g_cov(2,3)**2*fff+3*mu*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2+12*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)-6*mu*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)-6*fff*g_cov(1,2)**2*g_cov(3,3)**2+3*mu*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-fff*g_cov(3,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(21,36) = 2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,37) = -2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,38) = -2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,39) = 2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(21,40) = 2*alpha*(-3*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+3*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2+12*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2+9*mu*g_cov(1,2)*g_cov(3,3)*fff-9*mu*g_cov(1,3)*g_cov(2,3)*fff-32*phi**2*fff*g_cov(1,2)*g_cov(3,3)+32*phi**2*fff*g_cov(1,3)*g_cov(2,3))/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(21,41) = -2*mu*g_cov(1,3)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(21,42) = 0
      R47(21,43) = 2*mu*g_cov(1,2)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(21,44) = -2*mu*g_cov(1,2)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(21,45) = -2*alpha*(-3*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+3*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2+12*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2+9*mu*g_cov(1,2)*g_cov(3,3)*fff-9*mu*g_cov(1,3)*g_cov(2,3)*fff-32*phi**2*fff*g_cov(1,2)*g_cov(3,3)+32*phi**2*fff*g_cov(1,3)*g_cov(2,3))/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(21,46) = 2*mu*g_cov(1,3)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(21,47) = 0
      R47(22,1) = 0
      R47(22,2) = 0
      R47(22,3) = 0
      R47(22,4) = 0
      R47(22,5) = 0
      R47(22,6) = 0
      R47(22,7) = 0
      R47(22,8) = 0
      R47(22,9) = 0
      R47(22,10) = 0
      R47(22,11) = 0
      R47(22,12) = 0
      R47(22,13) = 0
      R47(22,14) = 0
      R47(22,15) = 0
      R47(22,16) = 0
      R47(22,17) = 0
      R47(22,18) = 0
      R47(22,19) = 0
      R47(22,20) = 0
      R47(22,21) = 0
      R47(22,22) = 0
      R47(22,23) = 0
      R47(22,24) = 0
      R47(22,25) = 0
      R47(22,26) = 0
      R47(22,27) = 0
      R47(22,28) = 0
      R47(22,29) = 0
      R47(22,30) = 0
      R47(22,31) = 0
      R47(22,32) = 3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,33) = 0
      R47(22,34) = 0
      R47(22,35) = -3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,36) = 0
      R47(22,37) = 0
      R47(22,38) = 0
      R47(22,39) = 0
      R47(22,40) = -(3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(22,41) = -2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,42) = 0
      R47(22,43) = 2*g_cov(2,2)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,44) = -2*g_cov(2,2)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,45) = (3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(22,46) = 2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(22,47) = 0
      R47(23,1) = 0
      R47(23,2) = 0
      R47(23,3) = 0
      R47(23,4) = 0
      R47(23,5) = 0
      R47(23,6) = 0
      R47(23,7) = 0
      R47(23,8) = 0
      R47(23,9) = 0
      R47(23,10) = 0
      R47(23,11) = 0
      R47(23,12) = 0
      R47(23,13) = 0
      R47(23,14) = 0
      R47(23,15) = 0
      R47(23,16) = 0
      R47(23,17) = 0
      R47(23,18) = 0
      R47(23,19) = 0
      R47(23,20) = 0
      R47(23,21) = 0
      R47(23,22) = 0
      R47(23,23) = 0
      R47(23,24) = 0
      R47(23,25) = 0
      R47(23,26) = 0
      R47(23,27) = 0
      R47(23,28) = 0
      R47(23,29) = 0
      R47(23,30) = 0
      R47(23,31) = 0
      R47(23,32) = 0
      R47(23,33) = 3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,34) = -3*mu*alpha*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,35) = 0
      R47(23,36) = 0
      R47(23,37) = 0
      R47(23,38) = 0
      R47(23,39) = 0
      R47(23,40) = 0
      R47(23,41) = -2*g_cov(3,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,42) = 0
      R47(23,43) = 2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,44) = -2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,45) = 0
      R47(23,46) = 2*g_cov(3,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(23,47) = 0
      R47(24,1) = 0
      R47(24,2) = 0
      R47(24,3) = 0
      R47(24,4) = 0
      R47(24,5) = 0
      R47(24,6) = 0
      R47(24,7) = -g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2) 
      R47(24,8) = 0
      R47(24,9) = 0
      R47(24,10) = 0
      R47(24,11) = 0
      R47(24,12) = 0
      R47(24,13) = 0
      R47(24,14) = 0
      R47(24,15) = 0
      R47(24,16) = phi*alpha*fff*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,17) = 0
      R47(24,18) = 0
      R47(24,19) = 0
      R47(24,20) = 0
      R47(24,21) = 0
      R47(24,22) = 0
      R47(24,23) = 0
      R47(24,24) = 0
      R47(24,25) = 0
      R47(24,26) = 0
      R47(24,27) = -phi*alpha*fff*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/sqrt(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,28) = -2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/phi/alpha/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,29) = 2*fff*sqrt(2.D0)*sqrt(alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2))*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/phi/alpha/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,30) = -2*fff*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(24,31) = 2*fff*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(3.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)
      R47(24,32) = -alpha*mu*sqrt(2.D0)*(6*fff*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)-3*mu*g_cov(2,3)**2*alpha**2*g_cov(1,2)*g_cov(1,3)+3*mu*g_cov(2,2)*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)-6*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+6*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*fff-3*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*alpha**2*mu+3*g_cov(2,3)**3*g_cov(1,1)*alpha**2*mu-6*g_cov(2,3)**3*g_cov(1,1)*fff-5*fff*g_cov(2,3)+3*g_cov(2,3)*mu*alpha**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(24,33) = -alpha*mu*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*(-6*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2*fff+3*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2*alpha**2*mu-3*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)*alpha**2*mu+6*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)*fff-3*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)*alpha**2*mu+6*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)*fff+5*g_cov(2,2)*fff-3*g_cov(2,2)*alpha**2*mu+3*mu*alpha**2*g_cov(1,2)**2*g_cov(2,3)**2-6*g_cov(1,2)**2*g_cov(2,3)**2*fff)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(24,34) = alpha*mu*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)*(-6*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2*fff+3*g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2*alpha**2*mu-3*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)*alpha**2*mu+6*g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)*fff-3*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)*alpha**2*mu+6*g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)*fff+5*g_cov(2,2)*fff-3*g_cov(2,2)*alpha**2*mu+3*mu*alpha**2*g_cov(1,2)**2*g_cov(2,3)**2-6*g_cov(1,2)**2*g_cov(2,3)**2*fff)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(24,35) = alpha*mu*sqrt(2.D0)*(6*fff*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)-3*mu*g_cov(2,3)**2*alpha**2*g_cov(1,2)*g_cov(1,3)+3*mu*g_cov(2,2)*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)-6*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)+6*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*fff-3*g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)*alpha**2*mu+3*g_cov(2,3)**3*g_cov(1,1)*alpha**2*mu-6*g_cov(2,3)**3*g_cov(1,1)*fff-5*fff*g_cov(2,3)+3*g_cov(2,3)*mu*alpha**2)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/(alpha**2*mu-2*fff)
      R47(24,36) = -2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,37) = 2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,38) = 2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,39) = -2*sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*fff)*g_cov(2,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(24,40) = -alpha*(3*mu**2*alpha**2*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)-6*mu**2*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(2,3)+3*mu**2*alpha**2*g_cov(2,3)**2*g_cov(1,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)**2*phi**2+24*mu*alpha**2*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)*phi**2-12*mu*alpha**2*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*phi**2-16*mu*g_cov(1,3)*g_cov(2,3)**2*fff-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*fff+18*mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)*fff+64*g_cov(1,3)*g_cov(2,3)**2*phi**2*fff-64*g_cov(3,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3))/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(24,41) = 2*mu*g_cov(2,3)*g_cov(1,3)*alpha/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(24,42) = 2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*alpha/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(24,43) = -2*mu*alpha*g_cov(1,3)*g_cov(2,2)/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(24,44) = 2*mu*alpha*g_cov(1,3)*g_cov(2,2)/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(24,45) = alpha*(3*mu**2*alpha**2*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)-6*mu**2*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(2,3)+3*mu**2*alpha**2*g_cov(2,3)**2*g_cov(1,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)**2*phi**2+24*mu*alpha**2*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)*phi**2-12*mu*alpha**2*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*phi**2-16*mu*g_cov(1,3)*g_cov(2,3)**2*fff-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*fff+18*mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)*fff+64*g_cov(1,3)*g_cov(2,3)**2*phi**2*fff-64*g_cov(3,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3))/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(24,46) = -2*mu*g_cov(2,3)*g_cov(1,3)*alpha/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(24,47) = -2*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu*alpha/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,1) = 0
      R47(25,2) = 0
      R47(25,3) = 0
      R47(25,4) = 0
      R47(25,5) = 0
      R47(25,6) = 0
      R47(25,7) = 0
      R47(25,8) = 0
      R47(25,9) = 0
      R47(25,10) = 0
      R47(25,11) = 0
      R47(25,12) = 0
      R47(25,13) = 0
      R47(25,14) = 0
      R47(25,15) = 0
      R47(25,16) = 0
      R47(25,17) = 0
      R47(25,18) = 0
      R47(25,19) = 0
      R47(25,20) = 0
      R47(25,21) = 0
      R47(25,22) = 0
      R47(25,23) = 0
      R47(25,24) = 0
      R47(25,25) = 0
      R47(25,26) = 0
      R47(25,27) = 0
      R47(25,28) = 0
      R47(25,29) = 0
      R47(25,30) = 0
      R47(25,31) = 0
      R47(25,32) = -3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,33) = 0
      R47(25,34) = 0
      R47(25,35) = 3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,36) = 0
      R47(25,37) = 0
      R47(25,38) = 0
      R47(25,39) = 0
      R47(25,40) = g_cov(2,3)*(3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(25,41) = 2*alpha*g_cov(2,3)**2*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,42) = -2*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,43) = -2*alpha*g_cov(2,2)*g_cov(2,3)*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,44) = 2*alpha*g_cov(2,2)*g_cov(2,3)*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,45) = -g_cov(2,3)*(3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)*alpha/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)
      R47(25,46) = -2*alpha*g_cov(2,3)**2*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(25,47) = 2*alpha*(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu/g_cov(3,3)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,1) = 0
      R47(26,2) = 0
      R47(26,3) = 0
      R47(26,4) = 0
      R47(26,5) = 0
      R47(26,6) = 0
      R47(26,7) = 0
      R47(26,8) = 0
      R47(26,9) = 0
      R47(26,10) = 0
      R47(26,11) = 0
      R47(26,12) = 0
      R47(26,13) = 0
      R47(26,14) = 0
      R47(26,15) = 0
      R47(26,16) = 0
      R47(26,17) = 0
      R47(26,18) = 0
      R47(26,19) = 0
      R47(26,20) = 0
      R47(26,21) = 0
      R47(26,22) = 0
      R47(26,23) = 0
      R47(26,24) = 0
      R47(26,25) = 0
      R47(26,26) = 0
      R47(26,27) = 0
      R47(26,28) = 0
      R47(26,29) = 0
      R47(26,30) = 0
      R47(26,31) = 0
      R47(26,32) = 0
      R47(26,33) = -3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,34) = 3*mu*alpha*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))*sqrt(2.D0)/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,35) = 0
      R47(26,36) = 0
      R47(26,37) = 0
      R47(26,38) = 0
      R47(26,39) = 0
      R47(26,40) = 0
      R47(26,41) = 2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,42) = 0
      R47(26,43) = -2*g_cov(2,2)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,44) = 2*g_cov(2,2)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,45) = 0
      R47(26,46) = -2*g_cov(2,3)*alpha*mu/sqrt((g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*mu)
      R47(26,47) = 0
      R47(27,1) = 0
      R47(27,2) = 0
      R47(27,3) = 0
      R47(27,4) = 0
      R47(27,5) = 0
      R47(27,6) = 0
      R47(27,7) = 0
      R47(27,8) = 0
      R47(27,9) = 0
      R47(27,10) = 0
      R47(27,11) = 0
      R47(27,12) = 0
      R47(27,13) = 1
      R47(27,14) = 0
      R47(27,15) = 0
      R47(27,16) = 0
      R47(27,17) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(27,18) = -2*(-g_cov(1,1)*g_cov(2,3)+g_cov(1,2)*g_cov(1,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,19) = (-g_cov(1,1)*g_cov(3,3)+g_cov(1,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,20) = -(g_cov(1,1)*g_cov(2,2)-g_cov(1,2)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,21) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(27,22) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(27,23) = -2*(-g_cov(1,1)*g_cov(2,3)+g_cov(1,2)*g_cov(1,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,24) = -(g_cov(1,1)*g_cov(2,2)-g_cov(1,2)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,25) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(27,26) = (-g_cov(1,1)*g_cov(3,3)+g_cov(1,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,27) = 0
      R47(27,28) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,29) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,30) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,31) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,32) = -(-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*fff*g_cov(1,1)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)*fff-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)*alpha**2+g_cov(1,3)*g_cov(2,3)*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)-4*fff*g_cov(1,3)*g_cov(2,3)**3*phi**2*g_cov(1,1)+2*mu*alpha**2*g_cov(1,3)*g_cov(1,1)*g_cov(2,3)**3*phi**2-g_cov(1,3)*g_cov(2,3)**3*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(1,3)*g_cov(2,3)**3*g_cov(1,1)*fff+6*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2+10*mu*g_cov(1,3)*g_cov(2,3)*fff-8*phi**2*fff*g_cov(1,3)*g_cov(2,3)-6*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-4*g_cov(2,2)*g_cov(3,3)**2*phi**2*g_cov(1,1)*g_cov(1,2)*fff+2*mu*g_cov(2,2)*g_cov(3,3)**2*phi**2*g_cov(1,1)*g_cov(1,2)*alpha**2-g_cov(1,2)*g_cov(3,3)**2*g_cov(2,2)*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)**2*fff*g_cov(1,1)*g_cov(1,2)-2*mu*g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**2*fff+4*fff*g_cov(1,2)*g_cov(3,3)*phi**2*g_cov(1,1)*g_cov(2,3)**2-2*g_cov(1,2)*g_cov(3,3)*mu*alpha**2*phi**2*g_cov(1,1)*g_cov(2,3)**2+g_cov(1,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+8*phi**2*fff*g_cov(1,2)*g_cov(3,3)+6*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)-6*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-10*mu*g_cov(1,2)*g_cov(3,3)*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(27,33) = (6*g_cov(2,2)*mu*alpha**2*g_cov(1,3)*phi**2+10*mu*g_cov(1,3)*g_cov(2,2)*fff-8*g_cov(1,3)*g_cov(2,2)*phi**2*fff-6*g_cov(2,2)*mu**2*alpha**2*g_cov(1,3)-2*mu*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*fff*g_cov(1,1)+4*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*fff*phi**2*g_cov(1,1)-2*mu*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*alpha**2+g_cov(1,3)*g_cov(2,2)**2*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)-4*g_cov(2,2)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)**2*fff+2*mu*g_cov(2,2)*g_cov(1,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2-g_cov(1,3)*g_cov(2,2)*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+2*mu*g_cov(2,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,1)*fff-4*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,2)*g_cov(2,3)*phi**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*alpha**2-g_cov(1,2)*g_cov(2,3)*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,1)*g_cov(1,2)*g_cov(2,3)-2*mu*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)*fff+4*phi**2*fff*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)-2*mu*alpha**2*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)*phi**2+g_cov(1,2)*g_cov(2,3)**3*mu**2*alpha**2*g_cov(1,1)+8*phi**2*g_cov(1,2)*g_cov(2,3)*fff+6*g_cov(1,2)*g_cov(2,3)*mu**2*alpha**2-6*g_cov(1,2)*g_cov(2,3)*mu*phi**2*alpha**2-10*g_cov(1,2)*g_cov(2,3)*mu*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(27,34) = (6*g_cov(2,2)*mu*alpha**2*g_cov(1,3)*phi**2+10*mu*g_cov(1,3)*g_cov(2,2)*fff-8*g_cov(1,3)*g_cov(2,2)*phi**2*fff-6*g_cov(2,2)*mu**2*alpha**2*g_cov(1,3)-2*mu*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*fff*g_cov(1,1)+4*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*fff*phi**2*g_cov(1,1)-2*mu*g_cov(2,2)**2*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*alpha**2+g_cov(1,3)*g_cov(2,2)**2*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)-4*g_cov(2,2)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)**2*fff+2*mu*g_cov(2,2)*g_cov(1,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2-g_cov(1,3)*g_cov(2,2)*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+2*mu*g_cov(2,2)*g_cov(1,3)*g_cov(2,3)**2*g_cov(1,1)*fff-4*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,2)*g_cov(2,3)*phi**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(2,3)*g_cov(1,1)*alpha**2-g_cov(1,2)*g_cov(2,3)*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)*fff*g_cov(1,1)*g_cov(1,2)*g_cov(2,3)-2*mu*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)*fff+4*phi**2*fff*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)-2*mu*alpha**2*g_cov(1,2)*g_cov(2,3)**3*g_cov(1,1)*phi**2+g_cov(1,2)*g_cov(2,3)**3*mu**2*alpha**2*g_cov(1,1)+8*phi**2*g_cov(1,2)*g_cov(2,3)*fff+6*g_cov(1,2)*g_cov(2,3)*mu**2*alpha**2-6*g_cov(1,2)*g_cov(2,3)*mu*phi**2*alpha**2-10*g_cov(1,2)*g_cov(2,3)*mu*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(27,35) = -(-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*fff*g_cov(1,1)*g_cov(2,3)+4*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)*fff-2*mu*g_cov(2,2)*g_cov(3,3)*g_cov(1,3)*phi**2*g_cov(1,1)*g_cov(2,3)*alpha**2+g_cov(1,3)*g_cov(2,3)*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)-4*fff*g_cov(1,3)*g_cov(2,3)**3*phi**2*g_cov(1,1)+2*mu*alpha**2*g_cov(1,3)*g_cov(1,1)*g_cov(2,3)**3*phi**2-g_cov(1,3)*g_cov(2,3)**3*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(1,3)*g_cov(2,3)**3*g_cov(1,1)*fff+6*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2+10*mu*g_cov(1,3)*g_cov(2,3)*fff-8*phi**2*fff*g_cov(1,3)*g_cov(2,3)-6*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-4*g_cov(2,2)*g_cov(3,3)**2*phi**2*g_cov(1,1)*g_cov(1,2)*fff+2*mu*g_cov(2,2)*g_cov(3,3)**2*phi**2*g_cov(1,1)*g_cov(1,2)*alpha**2-g_cov(1,2)*g_cov(3,3)**2*g_cov(2,2)*mu**2*alpha**2*g_cov(1,1)+2*mu*g_cov(2,2)*g_cov(3,3)**2*fff*g_cov(1,1)*g_cov(1,2)-2*mu*g_cov(1,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)**2*fff+4*fff*g_cov(1,2)*g_cov(3,3)*phi**2*g_cov(1,1)*g_cov(2,3)**2-2*g_cov(1,2)*g_cov(3,3)*mu*alpha**2*phi**2*g_cov(1,1)*g_cov(2,3)**2+g_cov(1,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+8*phi**2*fff*g_cov(1,2)*g_cov(3,3)+6*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)-6*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-10*mu*g_cov(1,2)*g_cov(3,3)*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(27,36) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(27,37) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(27,38) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(27,39) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(27,40) = (-32*mu*g_cov(1,3)**2*g_cov(2,3)**2*fff+6*mu**2*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2-12*mu*alpha**2*phi**2*g_cov(2,3)**2*g_cov(1,3)**2+64*phi**2*fff*g_cov(1,3)**2*g_cov(2,3)**2-128*g_cov(3,3)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*phi**2*fff-12*mu**2*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+64*mu*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+24*mu*alpha**2*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+mu**2*g_cov(2,2)*g_cov(3,3)**2*alpha**2*g_cov(1,1)+4*mu*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*fff-8*mu*g_cov(2,2)*g_cov(3,3)**2*alpha**2*phi**2*g_cov(1,1)-mu**2*g_cov(3,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2+8*mu*g_cov(3,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2-4*mu*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*fff-32*mu*fff*g_cov(1,2)**2*g_cov(3,3)**2+6*mu**2*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-12*mu*alpha**2*g_cov(1,2)**2*g_cov(3,3)**2*phi**2+64*g_cov(3,3)**2*g_cov(1,2)**2*phi**2*fff+24*g_cov(3,3)*mu*alpha**2*phi**2+4*g_cov(3,3)*mu*fff-6*g_cov(3,3)*mu**2*alpha**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/mu/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,41) = -4*g_cov(1,3)*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,42) = -4*g_cov(1,3)*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,43) = 4*(2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*mu-4*phi**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*mu-2*phi**2*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2*mu+2*phi**2*g_cov(1,1)*g_cov(2,3)**2-2*mu*g_cov(1,2)**2*g_cov(3,3)+4*phi**2*g_cov(1,2)**2*g_cov(3,3)-mu+2*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,44) = 4*(2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*mu-4*phi**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*mu-2*phi**2*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2*mu+2*phi**2*g_cov(1,1)*g_cov(2,3)**2-2*mu*g_cov(1,2)**2*g_cov(3,3)+4*phi**2*g_cov(1,2)**2*g_cov(3,3)-mu+2*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,45) = (-32*mu*g_cov(1,3)**2*g_cov(2,3)**2*fff+6*mu**2*alpha**2*g_cov(1,3)**2*g_cov(2,3)**2-12*mu*alpha**2*phi**2*g_cov(2,3)**2*g_cov(1,3)**2+64*phi**2*fff*g_cov(1,3)**2*g_cov(2,3)**2-128*g_cov(3,3)*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*phi**2*fff-12*mu**2*g_cov(3,3)*alpha**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+64*mu*fff*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+24*mu*alpha**2*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)+mu**2*g_cov(2,2)*g_cov(3,3)**2*alpha**2*g_cov(1,1)+4*mu*g_cov(2,2)*g_cov(3,3)**2*g_cov(1,1)*fff-8*mu*g_cov(2,2)*g_cov(3,3)**2*alpha**2*phi**2*g_cov(1,1)-mu**2*g_cov(3,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2+8*mu*g_cov(3,3)*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2-4*mu*g_cov(2,3)**2*g_cov(1,1)*g_cov(3,3)*fff-32*mu*fff*g_cov(1,2)**2*g_cov(3,3)**2+6*mu**2*g_cov(3,3)**2*alpha**2*g_cov(1,2)**2-12*mu*alpha**2*g_cov(1,2)**2*g_cov(3,3)**2*phi**2+64*g_cov(3,3)**2*g_cov(1,2)**2*phi**2*fff+24*g_cov(3,3)*mu*alpha**2*phi**2+4*g_cov(3,3)*mu*fff-6*g_cov(3,3)*mu**2*alpha**2)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/mu/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,46) = -4*g_cov(1,3)*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(27,47) = -4*g_cov(1,3)*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(28,1) = 0
      R47(28,2) = 0
      R47(28,3) = 0
      R47(28,4) = 0
      R47(28,5) = 0
      R47(28,6) = 0
      R47(28,7) = 0
      R47(28,8) = 1
      R47(28,9) = 0
      R47(28,10) = 0
      R47(28,11) = 0
      R47(28,12) = 0
      R47(28,13) = 0
      R47(28,14) = 0
      R47(28,15) = 0
      R47(28,16) = 0
      R47(28,17) = 1
      R47(28,18) = 0
      R47(28,19) = 0
      R47(28,20) = 0
      R47(28,21) = 0
      R47(28,22) = 1
      R47(28,23) = 0
      R47(28,24) = 0
      R47(28,25) = 0
      R47(28,26) = 0
      R47(28,27) = 0
      R47(28,28) = g_cov(1,2)
      R47(28,29) = g_cov(1,2)
      R47(28,30) = g_cov(1,2)
      R47(28,31) = g_cov(1,2)
      R47(28,32) = -(2*mu**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2*phi**2-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*fff+8*g_cov(1,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3)+4*mu*g_cov(1,2)**2*g_cov(3,3)*fff-2*mu**2*g_cov(1,2)**2*g_cov(3,3)*alpha**2-8*phi**2*fff*g_cov(1,2)**2*g_cov(3,3)+4*g_cov(3,3)*mu*alpha**2*g_cov(1,2)**2*phi**2-4*phi**2*fff-4*mu*fff+3*mu**2*alpha**2)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(28,33) = g_cov(1,2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(28,34) = g_cov(1,2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(28,35) = -(2*mu**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2*phi**2-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*fff+8*g_cov(1,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3)+4*mu*g_cov(1,2)**2*g_cov(3,3)*fff-2*mu**2*g_cov(1,2)**2*g_cov(3,3)*alpha**2-8*phi**2*fff*g_cov(1,2)**2*g_cov(3,3)+4*g_cov(3,3)*mu*alpha**2*g_cov(1,2)**2*phi**2-4*phi**2*fff-4*mu*fff+3*mu**2*alpha**2)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(28,36) = 0
      R47(28,37) = 1
      R47(28,38) = 0
      R47(28,39) = 1
      R47(28,40) = -(-11*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+40*mu*g_cov(1,2)*g_cov(3,3)*fff+28*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-64*phi**2*fff*g_cov(1,2)*g_cov(3,3)+9*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)+64*phi**2*fff*g_cov(1,3)*g_cov(2,3)-48*mu*g_cov(1,3)*g_cov(2,3)*fff-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2)/g_cov(3,3)/mu/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/2
      R47(28,41) = (2*mu*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-3*g_cov(1,3)*g_cov(2,3)**2*mu+mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)-4*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*phi**2+4*phi**2*g_cov(1,3)*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/g_cov(3,3)/(mu-4*phi**2)
      R47(28,42) = (2*mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)+mu*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-3*g_cov(1,3)*g_cov(2,3)**2*mu+4*phi**2*g_cov(1,3)*g_cov(2,3)**2-4*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/g_cov(3,3)/(mu-4*phi**2)
      R47(28,43) = g_cov(2,2)*(3*g_cov(1,3)*g_cov(2,3)*mu-3*g_cov(1,2)*g_cov(3,3)*mu+4*g_cov(1,2)*g_cov(3,3)*phi**2-4*g_cov(1,3)*g_cov(2,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(28,44) = g_cov(2,2)*(3*g_cov(1,3)*g_cov(2,3)*mu-3*g_cov(1,2)*g_cov(3,3)*mu+4*g_cov(1,2)*g_cov(3,3)*phi**2-4*g_cov(1,3)*g_cov(2,3)*phi**2)/(mu-4*phi**2)/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(28,45) = -(-11*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+40*mu*g_cov(1,2)*g_cov(3,3)*fff+28*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-64*phi**2*fff*g_cov(1,2)*g_cov(3,3)+9*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)+64*phi**2*fff*g_cov(1,3)*g_cov(2,3)-48*mu*g_cov(1,3)*g_cov(2,3)*fff-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2)/g_cov(3,3)/mu/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/2
      R47(28,46) = (2*mu*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-3*g_cov(1,3)*g_cov(2,3)**2*mu+mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)-4*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)*phi**2+4*phi**2*g_cov(1,3)*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/g_cov(3,3)/(mu-4*phi**2)
      R47(28,47) = (2*mu*g_cov(3,3)*g_cov(1,2)*g_cov(2,3)+mu*g_cov(1,3)*g_cov(2,2)*g_cov(3,3)-3*g_cov(1,3)*g_cov(2,3)**2*mu+4*phi**2*g_cov(1,3)*g_cov(2,3)**2-4*g_cov(3,3)*phi**2*g_cov(1,2)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/g_cov(3,3)/(mu-4*phi**2)
      R47(29,1) = 0
      R47(29,2) = 0
      R47(29,3) = 0
      R47(29,4) = 0
      R47(29,5) = 0
      R47(29,6) = 0
      R47(29,7) = 0
      R47(29,8) = 0
      R47(29,9) = 0
      R47(29,10) = 1
      R47(29,11) = 0
      R47(29,12) = 0
      R47(29,13) = 0
      R47(29,14) = 0
      R47(29,15) = 0
      R47(29,16) = 0
      R47(29,17) = 0
      R47(29,18) = 0
      R47(29,19) = 0
      R47(29,20) = 0
      R47(29,21) = 1
      R47(29,22) = 0
      R47(29,23) = 0
      R47(29,24) = 0
      R47(29,25) = 1
      R47(29,26) = 0
      R47(29,27) = 0
      R47(29,28) = g_cov(1,3)
      R47(29,29) = g_cov(1,3)
      R47(29,30) = g_cov(1,3)
      R47(29,31) = g_cov(1,3)
      R47(29,32) = -g_cov(1,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,33) = (-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*fff+2*mu**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2+8*g_cov(1,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3)-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2*phi**2-4*g_cov(2,2)*g_cov(3,3)*mu*g_cov(1,1)*fff+2*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)+8*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*phi**2*fff-4*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2-2*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+4*mu*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2+4*g_cov(1,1)*g_cov(2,3)**2*mu*fff-8*g_cov(1,1)*g_cov(2,3)**2*phi**2*fff-2*mu**2*g_cov(1,2)**2*g_cov(3,3)*alpha**2+4*g_cov(3,3)*mu*alpha**2*g_cov(1,2)**2*phi**2+4*mu*g_cov(1,2)**2*g_cov(3,3)*fff-8*phi**2*fff*g_cov(1,2)**2*g_cov(3,3)+4*mu*alpha**2*phi**2-4*phi**2*fff-5*mu**2*alpha**2+8*mu*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(29,34) = (-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*fff+2*mu**2*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2+8*g_cov(1,3)*phi**2*fff*g_cov(1,2)*g_cov(2,3)-4*mu*g_cov(1,2)*g_cov(1,3)*g_cov(2,3)*alpha**2*phi**2-4*g_cov(2,2)*g_cov(3,3)*mu*g_cov(1,1)*fff+2*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2*g_cov(1,1)+8*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*phi**2*fff-4*g_cov(1,1)*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2-2*mu**2*alpha**2*g_cov(1,1)*g_cov(2,3)**2+4*mu*alpha**2*g_cov(1,1)*g_cov(2,3)**2*phi**2+4*g_cov(1,1)*g_cov(2,3)**2*mu*fff-8*g_cov(1,1)*g_cov(2,3)**2*phi**2*fff-2*mu**2*g_cov(1,2)**2*g_cov(3,3)*alpha**2+4*g_cov(3,3)*mu*alpha**2*g_cov(1,2)**2*phi**2+4*mu*g_cov(1,2)**2*g_cov(3,3)*fff-8*phi**2*fff*g_cov(1,2)**2*g_cov(3,3)+4*mu*alpha**2*phi**2-4*phi**2*fff-5*mu**2*alpha**2+8*mu*fff)/(alpha**2*mu-2*fff)/(mu-2*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)/2
      R47(29,35) = -g_cov(1,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,36) = 1
      R47(29,37) = 0
      R47(29,38) = 1
      R47(29,39) = 0
      R47(29,40) = g_cov(1,3)
      R47(29,41) = -(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,42) = -2*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,43) = (3*g_cov(1,3)*g_cov(2,2)*mu-3*g_cov(1,2)*g_cov(2,3)*mu-4*g_cov(1,3)*g_cov(2,2)*phi**2+4*g_cov(1,2)*g_cov(2,3)*phi**2)/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,44) = (3*g_cov(1,3)*g_cov(2,2)*mu-3*g_cov(1,2)*g_cov(2,3)*mu-4*g_cov(1,3)*g_cov(2,2)*phi**2+4*g_cov(1,2)*g_cov(2,3)*phi**2)/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,45) = g_cov(1,3)
      R47(29,46) = -(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))*mu/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(29,47) = -2*(-2*g_cov(1,3)*g_cov(2,3)*phi**2+g_cov(1,3)*g_cov(2,3)*mu-g_cov(1,2)*g_cov(3,3)*mu+2*g_cov(1,2)*g_cov(3,3)*phi**2)/(mu-4*phi**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,1) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(30,2) = g_cov(2,2)
      R47(30,3) = -g_cov(2,2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,4) = 0
      R47(30,5) = g_cov(2,2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,6) = 0
      R47(30,7) = 0
      R47(30,8) = 0
      R47(30,9) = 0
      R47(30,10) = 0
      R47(30,11) = 0
      R47(30,12) = 0
      R47(30,13) = 0
      R47(30,14) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(30,15) = 0
      R47(30,16) = 0
      R47(30,17) = 0
      R47(30,18) = 0
      R47(30,19) = 1
      R47(30,20) = 0
      R47(30,21) = 0
      R47(30,22) = 0
      R47(30,23) = 0
      R47(30,24) = 0
      R47(30,25) = 0
      R47(30,26) = 1
      R47(30,27) = 0
      R47(30,28) = g_cov(2,2)
      R47(30,29) = g_cov(2,2)
      R47(30,30) = g_cov(2,2)
      R47(30,31) = g_cov(2,2)
      R47(30,32) = -g_cov(2,2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,33) = g_cov(2,2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,34) = g_cov(2,2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,35) = -g_cov(2,2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(30,36) = 0
      R47(30,37) = 0
      R47(30,38) = 0
      R47(30,39) = 0
      R47(30,40) = (4*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2-8*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-12*g_cov(2,2)*fff*g_cov(3,3))/g_cov(3,3)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(30,41) = 2*mu*g_cov(2,3)/(mu-4*phi**2)/g_cov(3,3)
      R47(30,42) = 2*mu*g_cov(2,3)/(mu-4*phi**2)/g_cov(3,3)
      R47(30,43) = -2*g_cov(2,2)*mu/(mu-4*phi**2)/g_cov(3,3)
      R47(30,44) = -2*g_cov(2,2)*mu/(mu-4*phi**2)/g_cov(3,3)
      R47(30,45) = (4*g_cov(2,2)*g_cov(3,3)*alpha**2*mu-3*g_cov(2,3)**2*alpha**2*mu+16*fff*g_cov(2,3)**2-8*g_cov(2,2)*alpha**2*g_cov(3,3)*phi**2-12*g_cov(2,2)*fff*g_cov(3,3))/g_cov(3,3)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(30,46) = 2*mu*g_cov(2,3)/(mu-4*phi**2)/g_cov(3,3)
      R47(30,47) = 2*mu*g_cov(2,3)/(mu-4*phi**2)/g_cov(3,3)
      R47(31,1) = 0
      R47(31,2) = g_cov(2,3)
      R47(31,3) = -g_cov(2,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,4) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(31,5) = g_cov(2,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,6) = 0
      R47(31,7) = 0
      R47(31,8) = 0
      R47(31,9) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(31,10) = 0
      R47(31,11) = 0
      R47(31,12) = 0
      R47(31,13) = 0
      R47(31,14) = 0
      R47(31,15) = 0
      R47(31,16) = 0
      R47(31,17) = 0
      R47(31,18) = 1
      R47(31,19) = 0
      R47(31,20) = 0
      R47(31,21) = 0
      R47(31,22) = 0
      R47(31,23) = 1
      R47(31,24) = 0
      R47(31,25) = 0
      R47(31,26) = 0
      R47(31,27) = 0
      R47(31,28) = g_cov(2,3)
      R47(31,29) = g_cov(2,3)
      R47(31,30) = g_cov(2,3)
      R47(31,31) = g_cov(2,3)
      R47(31,32) = -g_cov(2,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,33) = g_cov(2,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,34) = g_cov(2,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,35) = -g_cov(2,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(31,36) = 0
      R47(31,37) = 0
      R47(31,38) = 0
      R47(31,39) = 0
      R47(31,40) = g_cov(2,3)
      R47(31,41) = mu/(mu-4*phi**2)
      R47(31,42) = mu/(mu-4*phi**2)
      R47(31,43) = 0
      R47(31,44) = 0
      R47(31,45) = g_cov(2,3)
      R47(31,46) = mu/(mu-4*phi**2)
      R47(31,47) = mu/(mu-4*phi**2)
      R47(32,1) = 0
      R47(32,2) = g_cov(3,3)
      R47(32,3) = -g_cov(3,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,4) = 0
      R47(32,5) = g_cov(3,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,6) = 0
      R47(32,7) = 0
      R47(32,8) = 0
      R47(32,9) = 0
      R47(32,10) = 0
      R47(32,11) = 0
      R47(32,12) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(32,13) = 0
      R47(32,14) = 0
      R47(32,15) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(32,16) = 0
      R47(32,17) = 0
      R47(32,18) = 0
      R47(32,19) = 0
      R47(32,20) = 1
      R47(32,21) = 0
      R47(32,22) = 0
      R47(32,23) = 0
      R47(32,24) = 1
      R47(32,25) = 0
      R47(32,26) = 0
      R47(32,27) = 0
      R47(32,28) = g_cov(3,3)
      R47(32,29) = g_cov(3,3)
      R47(32,30) = g_cov(3,3)
      R47(32,31) = g_cov(3,3)
      R47(32,32) = -g_cov(3,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,33) = g_cov(3,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,34) = g_cov(3,3)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,35) = -g_cov(3,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(32,36) = 0
      R47(32,37) = 0
      R47(32,38) = 0
      R47(32,39) = 0
      R47(32,40) = g_cov(3,3)
      R47(32,41) = 0
      R47(32,42) = 0
      R47(32,43) = 2*mu/(mu-4*phi**2)
      R47(32,44) = 2*mu/(mu-4*phi**2)
      R47(32,45) = g_cov(3,3)
      R47(32,46) = 0
      R47(32,47) = 0
      R47(33,1) = (g_cov(1,2)**2*g_cov(3,3)**2-2*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+g_cov(1,3)**2*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(33,2) = 0
      R47(33,3) = 0
      R47(33,4) = 2*(-g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)+g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)-g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)+g_cov(2,3)**3*g_cov(1,1)+g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(33,5) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,6) = 0
      R47(33,7) = 0
      R47(33,8) = 0
      R47(33,9) = 0
      R47(33,10) = 0
      R47(33,11) = 0
      R47(33,12) = 0
      R47(33,13) = 0
      R47(33,14) = 0
      R47(33,15) = (g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2-g_cov(2,2)-g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)-g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)+g_cov(1,2)**2*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(33,16) = 0
      R47(33,17) = 0
      R47(33,18) = 0
      R47(33,19) = 0
      R47(33,20) = 0
      R47(33,21) = 0
      R47(33,22) = 0
      R47(33,23) = 0
      R47(33,24) = 0
      R47(33,25) = 0
      R47(33,26) = 0
      R47(33,27) = 0
      R47(33,28) = 0
      R47(33,29) = 0
      R47(33,30) = 0
      R47(33,31) = 0
      R47(33,32) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,33) = 0
      R47(33,34) = 0
      R47(33,35) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,36) = 0
      R47(33,37) = 0
      R47(33,38) = 0
      R47(33,39) = 0
      R47(33,40) = -(-3*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+16*mu*g_cov(1,2)*g_cov(3,3)*fff+12*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-64*phi**2*fff*g_cov(1,2)*g_cov(3,3)-16*mu*g_cov(1,3)*g_cov(2,3)*fff+64*phi**2*fff*g_cov(1,3)*g_cov(2,3)+3*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2)/g_cov(3,3)/mu/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(33,41) = -2*g_cov(2,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,42) = 2*g_cov(1,3)/g_cov(3,3)
      R47(33,43) = 2*g_cov(2,2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,44) = 2*g_cov(2,2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,45) = -(-3*mu**2*alpha**2*g_cov(1,2)*g_cov(3,3)+16*mu*g_cov(1,2)*g_cov(3,3)*fff+12*mu*alpha**2*g_cov(1,2)*g_cov(3,3)*phi**2-64*phi**2*fff*g_cov(1,2)*g_cov(3,3)-16*mu*g_cov(1,3)*g_cov(2,3)*fff+64*phi**2*fff*g_cov(1,3)*g_cov(2,3)+3*mu**2*alpha**2*g_cov(1,3)*g_cov(2,3)-12*mu*alpha**2*g_cov(1,3)*g_cov(2,3)*phi**2)/g_cov(3,3)/mu/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)
      R47(33,46) = -2*g_cov(2,3)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))/g_cov(3,3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(33,47) = 2*g_cov(1,3)/g_cov(3,3)
      R47(34,1) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(34,2) = 0
      R47(34,3) = 0
      R47(34,4) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(34,5) = 0
      R47(34,6) = 0
      R47(34,7) = 0
      R47(34,8) = 0
      R47(34,9) = 0
      R47(34,10) = 0
      R47(34,11) = 0
      R47(34,12) = 0
      R47(34,13) = 0
      R47(34,14) = 0
      R47(34,15) = 0
      R47(34,16) = 0
      R47(34,17) = 0
      R47(34,18) = 0
      R47(34,19) = 0
      R47(34,20) = 0
      R47(34,21) = 0
      R47(34,22) = 0
      R47(34,23) = 0
      R47(34,24) = 0
      R47(34,25) = 0
      R47(34,26) = 0
      R47(34,27) = 0
      R47(34,28) = 0
      R47(34,29) = 0
      R47(34,30) = 0
      R47(34,31) = 0
      R47(34,32) = g_cov(1,2)
      R47(34,33) = 0
      R47(34,34) = 0
      R47(34,35) = g_cov(1,2)
      R47(34,36) = 0
      R47(34,37) = 0
      R47(34,38) = 0
      R47(34,39) = 0
      R47(34,40) = (3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)/mu/2
      R47(34,41) = g_cov(2,3)/g_cov(3,3)
      R47(34,42) = g_cov(2,3)/g_cov(3,3)
      R47(34,43) = -g_cov(2,2)/g_cov(3,3)
      R47(34,44) = -g_cov(2,2)/g_cov(3,3)
      R47(34,45) = (3*g_cov(2,2)*g_cov(3,3)*mu**2*alpha**2-3*mu**2*alpha**2*g_cov(2,3)**2-12*g_cov(2,2)*g_cov(3,3)*mu*alpha**2*phi**2+12*mu*alpha**2*g_cov(2,3)**2*phi**2-16*g_cov(2,2)*g_cov(3,3)*mu*fff+16*g_cov(2,3)**2*mu*fff+64*g_cov(2,2)*g_cov(3,3)*phi**2*fff-64*g_cov(2,3)**2*phi**2*fff)/(-8*alpha**2*phi**2+4*fff+alpha**2*mu)/g_cov(3,3)/mu/2
      R47(34,46) = g_cov(2,3)/g_cov(3,3)
      R47(34,47) = g_cov(2,3)/g_cov(3,3)
      R47(35,1) = 0
      R47(35,2) = 0
      R47(35,3) = 0
      R47(35,4) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(35,5) = 0
      R47(35,6) = 0
      R47(35,7) = 0
      R47(35,8) = 0
      R47(35,9) = 0
      R47(35,10) = 0
      R47(35,11) = 0
      R47(35,12) = 0
      R47(35,13) = 0
      R47(35,14) = 0
      R47(35,15) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(35,16) = 0
      R47(35,17) = 0
      R47(35,18) = 0
      R47(35,19) = 0
      R47(35,20) = 0
      R47(35,21) = 0
      R47(35,22) = 0
      R47(35,23) = 0
      R47(35,24) = 0
      R47(35,25) = 0
      R47(35,26) = 0
      R47(35,27) = 0
      R47(35,28) = 0
      R47(35,29) = 0
      R47(35,30) = 0
      R47(35,31) = 0
      R47(35,32) = g_cov(1,3)
      R47(35,33) = 0
      R47(35,34) = 0
      R47(35,35) = g_cov(1,3)
      R47(35,36) = 0
      R47(35,37) = 0
      R47(35,38) = 0
      R47(35,39) = 0
      R47(35,40) = 0
      R47(35,41) = 0
      R47(35,42) = 1
      R47(35,43) = 0
      R47(35,44) = 0
      R47(35,45) = 0
      R47(35,46) = 0
      R47(35,47) = 1
      R47(36,1) = 1
      R47(36,2) = 0
      R47(36,3) = 0
      R47(36,4) = 0
      R47(36,5) = 0
      R47(36,6) = 0
      R47(36,7) = 0
      R47(36,8) = 0
      R47(36,9) = 0
      R47(36,10) = 0
      R47(36,11) = 0
      R47(36,12) = 0
      R47(36,13) = 0
      R47(36,14) = 0
      R47(36,15) = 0
      R47(36,16) = 0
      R47(36,17) = 0
      R47(36,18) = 0
      R47(36,19) = 0
      R47(36,20) = 0
      R47(36,21) = 0
      R47(36,22) = 0
      R47(36,23) = 0
      R47(36,24) = 0
      R47(36,25) = 0
      R47(36,26) = 0
      R47(36,27) = 0
      R47(36,28) = 0
      R47(36,29) = 0
      R47(36,30) = 0
      R47(36,31) = 0
      R47(36,32) = g_cov(2,2)
      R47(36,33) = 0
      R47(36,34) = 0
      R47(36,35) = g_cov(2,2)
      R47(36,36) = 0
      R47(36,37) = 0
      R47(36,38) = 0
      R47(36,39) = 0
      R47(36,40) = 0
      R47(36,41) = 0
      R47(36,42) = 0
      R47(36,43) = 0
      R47(36,44) = 0
      R47(36,45) = 0
      R47(36,46) = 0
      R47(36,47) = 0
      R47(37,1) = 0
      R47(37,2) = 0
      R47(37,3) = 0
      R47(37,4) = 1
      R47(37,5) = 0
      R47(37,6) = 0
      R47(37,7) = 0
      R47(37,8) = 0
      R47(37,9) = 0
      R47(37,10) = 0
      R47(37,11) = 0
      R47(37,12) = 0
      R47(37,13) = 0
      R47(37,14) = 0
      R47(37,15) = 0
      R47(37,16) = 0
      R47(37,17) = 0
      R47(37,18) = 0
      R47(37,19) = 0
      R47(37,20) = 0
      R47(37,21) = 0
      R47(37,22) = 0
      R47(37,23) = 0
      R47(37,24) = 0
      R47(37,25) = 0
      R47(37,26) = 0
      R47(37,27) = 0
      R47(37,28) = 0
      R47(37,29) = 0
      R47(37,30) = 0
      R47(37,31) = 0
      R47(37,32) = g_cov(2,3)
      R47(37,33) = 0
      R47(37,34) = 0
      R47(37,35) = g_cov(2,3)
      R47(37,36) = 0
      R47(37,37) = 0
      R47(37,38) = 0
      R47(37,39) = 0
      R47(37,40) = 0
      R47(37,41) = 0
      R47(37,42) = 0
      R47(37,43) = 0
      R47(37,44) = 0
      R47(37,45) = 0
      R47(37,46) = 0
      R47(37,47) = 0
      R47(38,1) = 0
      R47(38,2) = 0
      R47(38,3) = 0
      R47(38,4) = 0
      R47(38,5) = 0
      R47(38,6) = 0
      R47(38,7) = 0
      R47(38,8) = 0
      R47(38,9) = 0
      R47(38,10) = 0
      R47(38,11) = 0
      R47(38,12) = 0
      R47(38,13) = 0
      R47(38,14) = 0
      R47(38,15) = 1
      R47(38,16) = 0
      R47(38,17) = 0
      R47(38,18) = 0
      R47(38,19) = 0
      R47(38,20) = 0
      R47(38,21) = 0
      R47(38,22) = 0
      R47(38,23) = 0
      R47(38,24) = 0
      R47(38,25) = 0
      R47(38,26) = 0
      R47(38,27) = 0
      R47(38,28) = 0
      R47(38,29) = 0
      R47(38,30) = 0
      R47(38,31) = 0
      R47(38,32) = g_cov(3,3)
      R47(38,33) = 0
      R47(38,34) = 0
      R47(38,35) = g_cov(3,3)
      R47(38,36) = 0
      R47(38,37) = 0
      R47(38,38) = 0
      R47(38,39) = 0
      R47(38,40) = 0
      R47(38,41) = 0
      R47(38,42) = 0
      R47(38,43) = 0
      R47(38,44) = 0
      R47(38,45) = 0
      R47(38,46) = 0
      R47(38,47) = 0
      R47(39,1) = 0
      R47(39,2) = 0
      R47(39,3) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(39,4) = 0
      R47(39,5) = 0
      R47(39,6) = 0
      R47(39,7) = 0
      R47(39,8) = 0
      R47(39,9) = 2*(-g_cov(1,3)*g_cov(2,3)**2*g_cov(1,2)+g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,2)-g_cov(2,2)*g_cov(3,3)*g_cov(1,1)*g_cov(2,3)+g_cov(2,3)**3*g_cov(1,1)+g_cov(2,3))/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(39,10) = 0
      R47(39,11) = 0
      R47(39,12) = (g_cov(1,1)*g_cov(3,3)*g_cov(2,2)**2-g_cov(2,2)-g_cov(1,1)*g_cov(2,3)**2*g_cov(2,2)-g_cov(2,2)*g_cov(1,2)**2*g_cov(3,3)+g_cov(1,2)**2*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(39,13) = 0
      R47(39,14) = (g_cov(1,2)**2*g_cov(3,3)**2-2*g_cov(1,2)*g_cov(3,3)*g_cov(1,3)*g_cov(2,3)+g_cov(1,3)**2*g_cov(2,3)**2)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)**2
      R47(39,15) = 0
      R47(39,16) = 0
      R47(39,17) = 0
      R47(39,18) = 0
      R47(39,19) = 0
      R47(39,20) = 0
      R47(39,21) = 0
      R47(39,22) = 0
      R47(39,23) = 0
      R47(39,24) = 0
      R47(39,25) = 0
      R47(39,26) = 0
      R47(39,27) = 0
      R47(39,28) = 0
      R47(39,29) = 0
      R47(39,30) = 0
      R47(39,31) = 0
      R47(39,32) = 0
      R47(39,33) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(39,34) = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)**2-3)/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)
      R47(39,35) = 0
      R47(39,36) = 0
      R47(39,37) = 0
      R47(39,38) = 0
      R47(39,39) = 0
      R47(39,40) = 0
      R47(39,41) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(39,42) = 0
      R47(39,43) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(39,44) = 2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(39,45) = 0
      R47(39,46) = -2/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(39,47) = 0
      R47(40,1) = 0
      R47(40,2) = 0
      R47(40,3) = 0
      R47(40,4) = 0
      R47(40,5) = 0
      R47(40,6) = 0
      R47(40,7) = 0
      R47(40,8) = 0
      R47(40,9) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(40,10) = 0
      R47(40,11) = 0
      R47(40,12) = 0
      R47(40,13) = 0
      R47(40,14) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(40,15) = 0
      R47(40,16) = 0
      R47(40,17) = 0
      R47(40,18) = 0
      R47(40,19) = 0
      R47(40,20) = 0
      R47(40,21) = 0
      R47(40,22) = 0
      R47(40,23) = 0
      R47(40,24) = 0
      R47(40,25) = 0
      R47(40,26) = 0
      R47(40,27) = 0
      R47(40,28) = 0
      R47(40,29) = 0
      R47(40,30) = 0
      R47(40,31) = 0
      R47(40,32) = 0
      R47(40,33) = g_cov(1,2)
      R47(40,34) = g_cov(1,2)
      R47(40,35) = 0
      R47(40,36) = 0
      R47(40,37) = 0
      R47(40,38) = 0
      R47(40,39) = 0
      R47(40,40) = 0
      R47(40,41) = 1
      R47(40,42) = 0
      R47(40,43) = 0
      R47(40,44) = 0
      R47(40,45) = 0
      R47(40,46) = 1
      R47(40,47) = 0
      R47(41,1) = 0
      R47(41,2) = 0
      R47(41,3) = 0
      R47(41,4) = 0
      R47(41,5) = 0
      R47(41,6) = 0
      R47(41,7) = 0
      R47(41,8) = 0
      R47(41,9) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(41,10) = 0
      R47(41,11) = 0
      R47(41,12) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(41,13) = 0
      R47(41,14) = 0
      R47(41,15) = 0
      R47(41,16) = 0
      R47(41,17) = 0
      R47(41,18) = 0
      R47(41,19) = 0
      R47(41,20) = 0
      R47(41,21) = 0
      R47(41,22) = 0
      R47(41,23) = 0
      R47(41,24) = 0
      R47(41,25) = 0
      R47(41,26) = 0
      R47(41,27) = 0
      R47(41,28) = 0
      R47(41,29) = 0
      R47(41,30) = 0
      R47(41,31) = 0
      R47(41,32) = 0
      R47(41,33) = g_cov(1,3)
      R47(41,34) = g_cov(1,3)
      R47(41,35) = 0
      R47(41,36) = 0
      R47(41,37) = 0
      R47(41,38) = 0
      R47(41,39) = 0
      R47(41,40) = 0
      R47(41,41) = 0
      R47(41,42) = 0
      R47(41,43) = 1
      R47(41,44) = 1
      R47(41,45) = 0
      R47(41,46) = 0
      R47(41,47) = 0
      R47(42,1) = 0
      R47(42,2) = 0
      R47(42,3) = 0
      R47(42,4) = 0
      R47(42,5) = 0
      R47(42,6) = 0
      R47(42,7) = 0
      R47(42,8) = 0
      R47(42,9) = 0
      R47(42,10) = 0
      R47(42,11) = 0
      R47(42,12) = 0
      R47(42,13) = 0
      R47(42,14) = 1
      R47(42,15) = 0
      R47(42,16) = 0
      R47(42,17) = 0
      R47(42,18) = 0
      R47(42,19) = 0
      R47(42,20) = 0
      R47(42,21) = 0
      R47(42,22) = 0
      R47(42,23) = 0
      R47(42,24) = 0
      R47(42,25) = 0
      R47(42,26) = 0
      R47(42,27) = 0
      R47(42,28) = 0
      R47(42,29) = 0
      R47(42,30) = 0
      R47(42,31) = 0
      R47(42,32) = 0
      R47(42,33) = g_cov(2,2)
      R47(42,34) = g_cov(2,2)
      R47(42,35) = 0
      R47(42,36) = 0
      R47(42,37) = 0
      R47(42,38) = 0
      R47(42,39) = 0
      R47(42,40) = 0
      R47(42,41) = 0
      R47(42,42) = 0
      R47(42,43) = 0
      R47(42,44) = 0
      R47(42,45) = 0
      R47(42,46) = 0
      R47(42,47) = 0
      R47(43,1) = 0
      R47(43,2) = 0
      R47(43,3) = 0
      R47(43,4) = 0
      R47(43,5) = 0
      R47(43,6) = 0
      R47(43,7) = 0
      R47(43,8) = 0
      R47(43,9) = 1
      R47(43,10) = 0
      R47(43,11) = 0
      R47(43,12) = 0
      R47(43,13) = 0
      R47(43,14) = 0
      R47(43,15) = 0
      R47(43,16) = 0
      R47(43,17) = 0
      R47(43,18) = 0
      R47(43,19) = 0
      R47(43,20) = 0
      R47(43,21) = 0
      R47(43,22) = 0
      R47(43,23) = 0
      R47(43,24) = 0
      R47(43,25) = 0
      R47(43,26) = 0
      R47(43,27) = 0
      R47(43,28) = 0
      R47(43,29) = 0
      R47(43,30) = 0
      R47(43,31) = 0
      R47(43,32) = 0
      R47(43,33) = g_cov(2,3)
      R47(43,34) = g_cov(2,3)
      R47(43,35) = 0
      R47(43,36) = 0
      R47(43,37) = 0
      R47(43,38) = 0
      R47(43,39) = 0
      R47(43,40) = 0
      R47(43,41) = 0
      R47(43,42) = 0
      R47(43,43) = 0
      R47(43,44) = 0
      R47(43,45) = 0
      R47(43,46) = 0
      R47(43,47) = 0
      R47(44,1) = 0
      R47(44,2) = 0
      R47(44,3) = 0
      R47(44,4) = 0
      R47(44,5) = 0
      R47(44,6) = 0
      R47(44,7) = 0
      R47(44,8) = 0
      R47(44,9) = 0
      R47(44,10) = 0
      R47(44,11) = 0
      R47(44,12) = 1
      R47(44,13) = 0
      R47(44,14) = 0
      R47(44,15) = 0
      R47(44,16) = 0
      R47(44,17) = 0
      R47(44,18) = 0
      R47(44,19) = 0
      R47(44,20) = 0
      R47(44,21) = 0
      R47(44,22) = 0
      R47(44,23) = 0
      R47(44,24) = 0
      R47(44,25) = 0
      R47(44,26) = 0
      R47(44,27) = 0
      R47(44,28) = 0
      R47(44,29) = 0
      R47(44,30) = 0
      R47(44,31) = 0
      R47(44,32) = 0
      R47(44,33) = g_cov(3,3)
      R47(44,34) = g_cov(3,3)
      R47(44,35) = 0
      R47(44,36) = 0
      R47(44,37) = 0
      R47(44,38) = 0
      R47(44,39) = 0
      R47(44,40) = 0
      R47(44,41) = 0
      R47(44,42) = 0
      R47(44,43) = 0
      R47(44,44) = 0
      R47(44,45) = 0
      R47(44,46) = 0
      R47(44,47) = 0
      R47(45,1) = 0
      R47(45,2) = 1
      R47(45,3) = 0
      R47(45,4) = 0
      R47(45,5) = 0
      R47(45,6) = 0
      R47(45,7) = 0
      R47(45,8) = 0
      R47(45,9) = 0
      R47(45,10) = 0
      R47(45,11) = 0
      R47(45,12) = 0
      R47(45,13) = 0
      R47(45,14) = 0
      R47(45,15) = 0
      R47(45,16) = (-alpha**2*phi**2+fff) 
      R47(45,17) = 0
      R47(45,18) = 0
      R47(45,19) = 0
      R47(45,20) = 0
      R47(45,21) = 0
      R47(45,22) = 0
      R47(45,23) = 0
      R47(45,24) = 0
      R47(45,25) = 0
      R47(45,26) = 0
      R47(45,27) = (-alpha**2*phi**2+fff) 
      R47(45,28) = 1
      R47(45,29) = 1
      R47(45,30) = 1
      R47(45,31) = 1
      R47(45,32) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(45,33) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(45,34) = 1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))
      R47(45,35) = -1/(g_cov(2,2)*g_cov(3,3)-g_cov(2,3)**2)*(-g_cov(1,2)*g_cov(3,3)+g_cov(1,3)*g_cov(2,3))
      R47(45,36) = 0
      R47(45,37) = 0
      R47(45,38) = 0
      R47(45,39) = 0
      R47(45,40) = 1
      R47(45,41) = 0
      R47(45,42) = 0
      R47(45,43) = 0
      R47(45,44) = 0
      R47(45,45) = 1
      R47(45,46) = 0
      R47(45,47) = 0
      R47(46,1) = 0
      R47(46,2) = 0
      R47(46,3) = 0
      R47(46,4) = 0
      R47(46,5) = 1
      R47(46,6) = 0
      R47(46,7) = 0
      R47(46,8) = 0
      R47(46,9) = 0
      R47(46,10) = 0
      R47(46,11) = 0
      R47(46,12) = 0
      R47(46,13) = 0
      R47(46,14) = 0
      R47(46,15) = 0
      R47(46,16) = 0
      R47(46,17) = 0
      R47(46,18) = 0
      R47(46,19) = 0
      R47(46,20) = 0
      R47(46,21) = 0
      R47(46,22) = 0
      R47(46,23) = 0
      R47(46,24) = 0
      R47(46,25) = 0
      R47(46,26) = 0
      R47(46,27) = 0
      R47(46,28) = 0
      R47(46,29) = 0
      R47(46,30) = 0
      R47(46,31) = 0
      R47(46,32) = 1
      R47(46,33) = 0
      R47(46,34) = 0
      R47(46,35) = 1
      R47(46,36) = 0
      R47(46,37) = 0
      R47(46,38) = 0
      R47(46,39) = 0
      R47(46,40) = 0
      R47(46,41) = 0
      R47(46,42) = 0
      R47(46,43) = 0
      R47(46,44) = 0
      R47(46,45) = 0
      R47(46,46) = 0
      R47(46,47) = 0
      R47(47,1) = 0
      R47(47,2) = 0
      R47(47,3) = 1
      R47(47,4) = 0
      R47(47,5) = 0
      R47(47,6) = 0
      R47(47,7) = 0
      R47(47,8) = 0
      R47(47,9) = 0
      R47(47,10) = 0
      R47(47,11) = 0
      R47(47,12) = 0
      R47(47,13) = 0
      R47(47,14) = 0
      R47(47,15) = 0
      R47(47,16) = 0
      R47(47,17) = 0
      R47(47,18) = 0
      R47(47,19) = 0
      R47(47,20) = 0
      R47(47,21) = 0
      R47(47,22) = 0
      R47(47,23) = 0
      R47(47,24) = 0
      R47(47,25) = 0
      R47(47,26) = 0
      R47(47,27) = 0
      R47(47,28) = 0
      R47(47,29) = 0
      R47(47,30) = 0
      R47(47,31) = 0
      R47(47,32) = 0
      R47(47,33) = 1
      R47(47,34) = 1
      R47(47,35) = 0
      R47(47,36) = 0
      R47(47,37) = 0
      R47(47,38) = 0
      R47(47,39) = 0
      R47(47,40) = 0
      R47(47,41) = 0
      R47(47,42) = 0
      R47(47,43) = 0
      R47(47,44) = 0
      R47(47,45) = 0
      R47(47,46) = 0
      R47(47,47) = 0       
        !
        CALL RG(47,47,Atest47,Lambda47,ImLambda47,1,R47,itemp47,rtemp47,ierr)   
        ! 
        t1 = MAXVAL(ABS(R47)) 
        !
        CALL MatrixInverse(47,R47,iR47) 
        t1 = MAXVAL(ABS(iR47)) 
        Id47 = 0.0 
        DO i = 1, 47 
            Id47(i,i) = 1.0
        ENDDO        
        !
        adiff47 = MATMUL( R47, iR47 ) - Id47 
        t1 = MAXVAL(ABS(adiff47))  
        continue
        !        
        adiff47 = MATMUL( Atest47, R47 )  - MATMUL( R47, L47 ) 
        t1 = MAXVAL(ABS(adiff47))  
        continue
        ! 
        A47 = MATMUL( R47, MATMUL(L47, iR47) ) 
        adiff47 = Atest47 - A47         
        !
        t1 = MAXVAL(ABS(adiff47))  
        continue
        DO i = 1, 47
            t1 = MAXVAL(ABS(adiff47(i,:))) 
            continue
        ENDDO        
        continue
        ! 
    ENDIF
    !
    !IF( t1 > 1e6 ) THEN
    !    PRINT *, ' Problem, iR = ', iR
    !    STOP 
    !ENDIF    
    !
#endif     
    !
END SUBROUTINE PDEEigenvectors      

SUBROUTINE PDEAbsA(absA,Q,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    USE recipies_mod, ONLY : RG 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: absA(nVar,nVar) 
    ! Local variables 
    INTEGER :: i 
    REAL    :: R(nVar,nVar), L(nVar), iR(nVar,nVar), absL(nVar,nVar) 
    REAL    :: t1, smax  
    !     
    CALL PDEEigenvectors(R,L,iR,Q,nv,par) 
    !
    absL = 0.0 
    DO i = 1, nVar
        absL(i,i) = ABS(L(i)) 
    ENDDO    
    !
    absA = MATMUL( R, MATMUL( absL, iR ) ) 
    !
    t1 = MAXVAL(ABS(absA)) 
    !
    IF( t1 > 1e4 ) THEN
        !PRINT *, ' Problem ' 
        !CONTINUE
        smax = MAXVAL( ABS(L) ) 
        absA = 0.0        
        DO i = 1, nVar
            absA(i,i) = smax
        ENDDO        
    ENDIF    
    !
    !WHERE(ABS(absA)<1e-12)
    !    absA = 0.0
    !ENDWHERE
    !
END SUBROUTINE PDEAbsA       
    
SUBROUTINE PDECons2Prim(V,Q,iErr)
    USE typesDef, ONLY : nVar, nParam, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)     :: Q(nVar)     ! vector of conserved quantities 
    REAL, INTENT(OUT)    :: V(nVar)     ! primitive variables 
    INTEGER, INTENT(OUT) :: iErr        ! error flag 
    ! Local variables 
    REAL                 :: p 
    REAL                 :: gamma1, gam, sb, dr, eps, sb2, sx, sy, sz, e, bx, by, bz, s2, b2 
    REAL                 :: x1, x2, x3, v2
    REAL                 :: w, rho, vx, vy, vz, den, vb 
    LOGICAL              :: FAILED
    REAL, PARAMETER      :: tol = 1e-8, third=1.0/3.0, p_floor = 1.0e-5, rho_floor = 1.0e-4    
    REAL                 :: RTSAFE_C2P_RMHD1 
    REAL                 :: lapse, shift(3), psi, gammaij(6), g_cov(3,3), g_contr(3,3), gp, gm, dd 
    REAL                 :: Qloc(nVar), bv(3), sm_cov(3), sm(3), bv_contr(3), vf(3), vf_cov(3)   
    REAL                 :: A(3,3), G(3,3), devG(3,3), Id(3,3), tempp(3,3), evv, ehh  
    REAL                 :: QGRMHD(19), VGRMHD(19)  
    !    
    iErr = 0     
    !
#ifdef EULER         
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    V(1) = Q(1)             ! fluid density 
    V(2:4) = Q(2:4)/Q(1)    ! fluid velocity 
    V(5)   = p              ! fluid pressure 
#endif 
    ! 
#ifdef ELASTICITY
    V = Q 
#endif 
    !
#ifdef ACOUSTIC 
    V = Q 
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    V = Q 
    !
#ifdef Z4GRMHD
    QGRMHD(1:9)   = Q(55:63)    ! hydro variables 
    QGRMHD(10)    = Q(17)       ! lapse 
    QGRMHD(11:13) = Q(18:20)    ! shift 
    QGRMHD(14:19) = Q(1:6)      ! metric 
    CALL PDECons2PrimGRMHD(VGRMHD,QGRMHD,iErr) 
    V(55:63) = VGRMHD(1:9) 
#endif     
#endif
    !
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    V = Q
    V(17) = EXP(MAX(-20.,MIN(20.,Q(17))))  
    V(55) = EXP(MAX(-20.,MIN(20.,Q(55))))  
    !
#ifdef CCZ4GRMHD
    QGRMHD(1:9)   = Q(60:68)    ! hydro variables 
    QGRMHD(10)    = EXP(Q(17))  ! lapse 
    QGRMHD(11:13) = Q(18:20)    ! shift 
    QGRMHD(14:19) = Q(1:6)/( EXP(Q(55)) )**2 ! metric 
    CALL PDECons2PrimGRMHD(VGRMHD,QGRMHD,iErr) 
    V(60:68) = VGRMHD(1:9) 
#endif     
#endif
    !
#if defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
    V = Q
    V(17) = EXP(Q(17))
    V(55) = EXP(Q(55)) 
    !
#ifdef BSSZ4GRMHD
    QGRMHD(1:9)   = Q(63:71)    ! hydro variables 
    QGRMHD(10)    = EXP(Q(17))  ! lapse 
    QGRMHD(11:13) = Q(18:20)    ! shift 
    QGRMHD(14:19) = Q(1:6)/( EXP(Q(55)) )**2 ! metric 
    CALL PDECons2PrimGRMHD(VGRMHD,QGRMHD,iErr) 
    V(63:71) = VGRMHD(1:9) 
#endif     
#endif

    !
#ifdef SRMHD
  gamma1 = EQN%gamma/(EQN%gamma - 1.0)
  gam    = 1.0/gamma1
  dr   = Q(1)
  sx   = Q(2)
  sy   = Q(3)
  sz   = Q(4)
  e    = Q(5) 
  bx   = Q(6)
  by   = Q(7)
  bz   = Q(8)
  !
  s2   = sx*sx + sy*sy + sz*sz
  b2   = bx*bx + by*by + bz*bz
  sb   = sx*bx + sy*by + sz*bz
  sb2  = sb**2
  eps  = 1.e-10
  !
  x1   = 0.
  x2   = 1.-eps
  v2   = RTSAFE_C2P_RMHD1(x1,x2,tol,gam,dr,e,s2,b2,sb2,w,FAILED)
  !
  IF (FAILED) THEN
     iErr = -1
     p    = p_floor
     rho  = rho_floor
     vx   = 0.0
     vy   = 0.0
     vz   = 0.0
     bx   = bx
     by   = by
     bz   = bz
  ELSE
     den  = 1.0/(w+b2)
     vb   = sb/w
     !
     rho  = dr*sqrt(1.-v2)
     vx   = (sx + vb*bx)*den
     vy   = (sy + vb*by)*den
     vz   = (sz + vb*bz)*den
     p    = max(1.e-15, gam*(w*(1.-v2)-rho))
  ENDIF

  V(1:9) = (/ rho, vx, vy, vz, p, bx, by, bz, Q(9) /)
#endif 
    !
#ifdef GRMHD
  !
  CALL PDECons2PrimGRMHD(V,Q,iErr) 
  RETURN 
  !  
#endif
    !
#ifdef GPR3D
    
    ! Peshkov-Romenski model with heat conduction 
    V(1)   = Q(1)        ! rho 
    V(2:4) = Q(2:4)/Q(1) ! u, v 
    V(6:14) = Q(6:14)      ! A11, A12, A13, A21, A22, A23, A31, A32, A33 
    V(15:17) = Q(15:17)/Q(1) ! j1, j2, j3  
    e      = Q(5)/Q(1) - 0.5*(V(2)**2+V(3)**2+V(4)**2) -0.5*(V(15)**2+V(16)**2+V(17)**2)*EQN%alpha**2      ! e = eh + ev 
    A(1,:) = (/ V( 6), V( 7),  V( 8) /) 
    A(2,:) = (/ V( 9), V(10),  V(11) /)
    A(3,:) = (/ V(12), V(13),  V(14) /)         
    G      = MATMUL( TRANSPOSE(A), A ) 
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
    tempp  = MATMUL( TRANSPOSE(devG), devG ) 
    evv    = EQN%cs**2/4.*(tempp(1,1)+tempp(2,2)+tempp(3,3))         
    ehh    = e-evv  
    V(5)   = ehh*(EQN%gamma-1.0)*V(1) - EQN%gamma*EQN%p0 

#endif 
    !
END SUBROUTINE PDECons2Prim 
        
SUBROUTINE PDEPrim2Cons(Q,V)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)      :: V(nVar)     ! primitive variables 
    REAL, INTENT(OUT)     :: Q(nVar)     ! vector of conserved quantities 
    ! Local variables 
    REAL                  :: rho, vx, vy, vz, p, bx, by, bz, ex, ey, ez, v2, b2, e2    
    REAL                  :: lf, gamma1, w, ww, uem 
    REAL                  :: gp, gm, g_cov(3,3), g_contr(3,3), bv(3), vf_cov(3), psi, lapse, shift(3)
    REAL                  :: vf(3), bv_contr(3), qv_contr(3), qb_contr(3), vxb(3), vb_cov(3), b2_cov, vxb_contr(3) 
    REAL                  :: QGRMHD(19), VGRMHD(19)
    REAL                  :: A(3,3), G(3,3), Id(3,3), devG(3,3), eh, evv, T, falpha, Temp(3,3)  
    !   
#ifdef EULER 
    !
    Q(1)   = V(1)           ! fluid density 
    Q(2:4) = V(1)*V(2:4)    ! momentum 
    Q(5)   = V(5)/(EQN%gamma-1) + 0.5*V(1)*SUM(V(2:4)**2)   ! total energy = internal energy + kinetic energy 
    !
#endif 
    !
#ifdef ELASTICITY 
    Q = V 
#endif 
    !
#ifdef ACOUSTIC 
    Q = V 
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    Q = V 
    !
#ifdef Z4GRMHD
    VGRMHD(1:9)   = V(55:63)    ! hydro variables 
    VGRMHD(10)    = V(17)       ! lapse 
    VGRMHD(11:13) = V(18:20)    ! shift 
    VGRMHD(14:19) = V(1:6)      ! metric 
    CALL PDEPrim2ConsGRMHD(QGRMHD,VGRMHD) 
    Q(55:63) = QGRMHD(1:9) 
#endif 
    ! 
#endif
    !
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    Q = V   
    Q(17) = LOG(V(17))
    Q(55) = LOG(V(55)) 
    !
#ifdef CCZ4GRMHD
    VGRMHD(1:9)   = V(60:68)    ! hydro variables 
    VGRMHD(10)    = V(17)       ! lapse 
    VGRMHD(11:13) = V(18:20)    ! shift 
    VGRMHD(14:19) = V(1:6)/( V(55) )**2   ! metric 
    CALL PDEPrim2ConsGRMHD(QGRMHD,VGRMHD) 
    Q(60:68) = QGRMHD(1:9) 
#endif 
#endif

#if defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
    Q = V   
    Q(17) = LOG(V(17))
    Q(55) = LOG(V(55)) 
    !
#ifdef BSSZ4GRMHD
    VGRMHD(1:9)   = V(63:71)    ! hydro variables 
    VGRMHD(10)    = EXP(V(17))  ! lapse 
    VGRMHD(11:13) = V(18:20)    ! shift 
    VGRMHD(14:19) = V(1:6)/( EXP(V(55)) )**2   ! metric 
    CALL PDEPrim2ConsGRMHD(QGRMHD,VGRMHD) 
    Q(63:71) = QGRMHD(1:9) 
#endif 
#endif

    !
#ifdef SRMHD
    rho    = V(1)
    vx     = V(2)
    vy     = V(3)
    vz     = V(4)
    p      = V(5)
    bx     = V(6)
    by     = V(7)
    bz     = V(8)
    !
    ex     = - (vy*bz - vz*by)
    ey     = - (vz*bx - vx*bz)
    ez     = - (vx*by - vy*bx)
    !
    v2     = vx**2 + vy**2 + vz**2
    b2     = bx**2 + by**2 + bz**2
    e2     = ex**2 + ey**2 + ez**2
    !
    IF (v2 > 1.0) THEN
        WRITE(*,*)'Superluminal velocity in PDEPrim2Cons!!'
        STOP
    ENDIF
    lf     = 1.0 / sqrt(1.0 - v2)
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    w      = rho + gamma1*p
    ww     = w*lf**2
    uem    = 0.5*(b2+e2)
    !
    Q(1)   = rho*lf
    Q(2)   = ww*vx + (ey*bz - ez*by)
    Q(3)   = ww*vy + (ez*bx - ex*bz)
    Q(4)   = ww*vz + (ex*by - ey*bx)
    Q(5)   = ww - p + uem 
    !
    Q(6)   = bx
    Q(7)   = by
    Q(8)   = bz
    Q(9)   = V(9)  
#endif 
    !
#ifdef GRMHD
  !
  CALL PDEPrim2ConsGRMHD(Q,V) 
  RETURN 
  !
#endif
    !
#ifdef GPR3D
     
    ! Peshkov-Romenski model in 3D with heat conduction 
    Q(1)   = V(1)        ! rho 
    Q(2:4) = V(1)*V(2:4) ! rho*u, rho*v, rho*w  
    Q(6:14) = V(6:14)    ! A11, A12, A13, A21, A22, A23, A31, A32, A33 
    eh     = (V(5)+EQN%gamma*EQN%p0)/(EQN%gamma-1)/V(1) 
    A(1,:) = (/ V(6),  V(7),  V(8)  /) 
    A(2,:) = (/ V(9),  V(10), V(11) /)
    A(3,:) = (/ V(12), V(13), V(14) /)         
    G      = MATMUL( TRANSPOSE(A), A ) 
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id  
    temp   = MATMUL( TRANSPOSE(devG), devG ) 
    evv     = EQN%cs**2/4.*(temp(1,1)+temp(2,2)+temp(3,3)) 
    ! Compute the temperature from the ideal gas law 
    p = V(5) 
    T = p/V(1)/EQN%cv/(EQN%gamma-1)  
    falpha  = EQN%alpha**2  
    
    Q(5)     = V(1)*(eh+evv) + 0.5*V(1)*(V(2)**2+V(3)**2+V(4)**2) + 0.5*V(1)*falpha*(V(15)**2 + V(16)**2 + V(17)**2) ! total energy rhoE     
    Q(15:17) = V(15:17)*V(1) ! rho j 

#endif 
    !
END SUBROUTINE PDEPrim2Cons      
            
SUBROUTINE PDEVarName(Name) 
    USE typesDef, ONLY : nVar, d, EQN 
    IMPLICIT NONE       
    ! Argument list 
    CHARACTER(LEN=10), INTENT(OUT) :: Name(nVar)
    !
#ifdef EULER 
    !
    Name(1) = 'rho' 
    Name(2) = 'u' 
    Name(3) = 'v'
    Name(4) = 'w' 
    Name(5) = 'p' 
    !
#endif 
    !
#ifdef ELASTICITY 
    Name(1)  = 'sxx' 
    Name(2)  = 'syy' 
    Name(3)  = 'szz'
    Name(4)  = 'sxy' 
    Name(5)  = 'syz' 
    Name(6)  = 'sxz' 
    Name(7)  = 'u' 
    Name(8)  = 'v' 
    Name(9)  = 'w' 
    Name(10) = 'lambda'
    Name(11) = 'mu'
    Name(12) = 'rho' 
    Name(13) = 'alpha'
    Name(14) = 'xi' 
    !
!    Name(1) = 'p' 
!    Name(2) = 'u' 
!    Name(3) = 'v'
!    Name(4) = 'w' 

#endif 
    !
#ifdef ACOUSTIC 
    Name(1) = 'p' 
    Name(2) = 'u' 
    Name(3) = 'v'
    Name(4) = 'w' 
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
  Name(:)  = 'Unknown'
  Name(1)  = 'g11' 
  Name(2)  = 'g12' 
  Name(3)  = 'g13'
  Name(4)  = 'g22'
  Name(5)  = 'g23' 
  Name(6)  = 'g33' 
  Name(7)  = 'K11' 
  Name(8)  = 'K12' 
  Name(9)  = 'K13'
  Name(10) = 'K22'
  Name(11) = 'K23' 
  Name(12) = 'K33'   
  Name(13) = 'Z1'
  Name(14) = 'Z2'
  Name(15) = 'Z3'
  Name(16) = 'Theta'
  Name(17) = 'lapse'
  Name(18) = 'shift1'
  Name(19) = 'shift2'
  Name(20) = 'shift3'
  Name(21) = 'b1'
  Name(22) = 'b2'
  Name(23) = 'b3'
  Name(24) = 'A1'
  Name(25) = 'A2'
  Name(26) = 'A3'
  Name(27) = 'B11'
  Name(28) = 'B21'
  Name(29) = 'B31'
  Name(30) = 'B12'
  Name(31) = 'B22'
  Name(32) = 'B32'
  Name(33) = 'B13'
  Name(34) = 'B23'
  Name(35) = 'B33' 
  Name(36) = 'D111'
  Name(37) = 'D112'
  Name(38) = 'D113'
  Name(39) = 'D122'
  Name(40) = 'D123'
  Name(41) = 'D133'
  Name(42) = 'D211'
  Name(43) = 'D212'
  Name(44) = 'D213'
  Name(45) = 'D222'
  Name(46) = 'D223'
  Name(47) = 'D233'
  Name(48) = 'D311'
  Name(49) = 'D312'
  Name(50) = 'D313'
  Name(51) = 'D322'
  Name(52) = 'D323'
  Name(53) = 'D333'
  Name(54) = 'K0'  
#ifdef Z4GRMHD
  Name(55) = 'rho' 
  Name(56) = 'u' 
  Name(57) = 'v'
  Name(58) = 'w'
  Name(59) = 'p' 
  Name(60) = 'bx' 
  Name(61) = 'by' 
  Name(62) = 'bz'
  Name(63) = 'phi'    
#endif 
    !
#endif 
    !
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
  Name(:)  = 'Unknown'
  Name(1)  = 'g11' 
  Name(2)  = 'g12' 
  Name(3)  = 'g13'
  Name(4)  = 'g22'
  Name(5)  = 'g23' 
  Name(6)  = 'g33' 
  Name(7)  = 'A11' 
  Name(8)  = 'A12' 
  Name(9)  = 'A13'
  Name(10) = 'A22'
  Name(11) = 'A23' 
  Name(12) = 'A33'   
  Name(13) = 'Theta' 
  Name(14) = 'G1'
  Name(15) = 'G2'
  Name(16) = 'G3'
  Name(17) = 'lapse'
  Name(18) = 'shift1'
  Name(19) = 'shift2'
  Name(20) = 'shift3'
  Name(21) = 'b1'
  Name(22) = 'b2'
  Name(23) = 'b3'
  Name(24) = 'A1'
  Name(25) = 'A2'
  Name(26) = 'A3'
  Name(27) = 'B11'
  Name(28) = 'B21'
  Name(29) = 'B31'
  Name(30) = 'B12'
  Name(31) = 'B22'
  Name(32) = 'B32'
  Name(33) = 'B13'
  Name(34) = 'B23'
  Name(35) = 'B33' 
  Name(36) = 'D111'
  Name(37) = 'D112'
  Name(38) = 'D113'
  Name(39) = 'D122'
  Name(40) = 'D123'
  Name(41) = 'D133'
  Name(42) = 'D211'
  Name(43) = 'D212'
  Name(44) = 'D213'
  Name(45) = 'D222'
  Name(46) = 'D223'
  Name(47) = 'D233'
  Name(48) = 'D311'
  Name(49) = 'D312'
  Name(50) = 'D313'
  Name(51) = 'D322'
  Name(52) = 'D323'
  Name(53) = 'D333'
  Name(54) = 'K'  
  Name(55) = 'phi'  
  Name(56) = 'P1'  
  Name(57) = 'P2'  
  Name(58) = 'P3'    
  Name(59) = 'K0'  
#ifdef CCZ4GRMHD 
  Name(60) = 'rho' 
  Name(61) = 'u' 
  Name(62) = 'v'
  Name(63) = 'w'
  Name(64) = 'p' 
  Name(65) = 'bx' 
  Name(66) = 'by' 
  Name(67) = 'bz'
  Name(68) = 'psi' 
#endif 
    !
#endif 
  !
#if defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
  Name(:)  = 'Unknown'
  Name(1)  = 'g11' 
  Name(2)  = 'g12' 
  Name(3)  = 'g13'
  Name(4)  = 'g22'
  Name(5)  = 'g23' 
  Name(6)  = 'g33' 
  Name(7)  = 'A11' 
  Name(8)  = 'A12' 
  Name(9)  = 'A13'
  Name(10) = 'A22'
  Name(11) = 'A23' 
  Name(12) = 'A33'   
  Name(13) = 'Theta' 
  Name(14) = 'G1'
  Name(15) = 'G2'
  Name(16) = 'G3'
  Name(17) = 'lapse'
  Name(18) = 'shift1'
  Name(19) = 'shift2'
  Name(20) = 'shift3'
  Name(21) = 'b1'
  Name(22) = 'b2'
  Name(23) = 'b3'
  Name(24) = 'A1'
  Name(25) = 'A2'
  Name(26) = 'A3'
  Name(27) = 'B11'
  Name(28) = 'B21'
  Name(29) = 'B31'
  Name(30) = 'B12'
  Name(31) = 'B22'
  Name(32) = 'B32'
  Name(33) = 'B13'
  Name(34) = 'B23'
  Name(35) = 'B33' 
  Name(36) = 'D111'
  Name(37) = 'D112'
  Name(38) = 'D113'
  Name(39) = 'D122'
  Name(40) = 'D123'
  Name(41) = 'D133'
  Name(42) = 'D211'
  Name(43) = 'D212'
  Name(44) = 'D213'
  Name(45) = 'D222'
  Name(46) = 'D223'
  Name(47) = 'D233'
  Name(48) = 'D311'
  Name(49) = 'D312'
  Name(50) = 'D313'
  Name(51) = 'D322'
  Name(52) = 'D323'
  Name(53) = 'D333'
  Name(54) = 'K'  
  Name(55) = 'phi'  
  Name(56) = 'P1'  
  Name(57) = 'P2'  
  Name(58) = 'P3'    
  Name(59) = 'X1'
  Name(60) = 'X2'
  Name(61) = 'X3'
  Name(62) = 'K0'  
#ifdef BSSZ4GRMHD 
  Name(63) = 'rho' 
  Name(64) = 'u' 
  Name(65) = 'v'
  Name(66) = 'w'
  Name(67) = 'p' 
  Name(68) = 'bx' 
  Name(69) = 'by' 
  Name(70) = 'bz'
  Name(71) = 'psi' 
#endif 
    !
#endif 
  !
#ifdef SRMHD
  Name(1) = 'rho' 
  Name(2) = 'u' 
  Name(3) = 'v'
  Name(4) = 'w'
  Name(5) = 'p' 
  Name(6) = 'bx' 
  Name(7) = 'by' 
  Name(8) = 'bz'
  Name(9) = 'phi'    
#endif 
    !
#ifdef GRMHD
  Name(1) = 'rho' 
  Name(2) = 'u' 
  Name(3) = 'v'
  Name(4) = 'w'
  Name(5) = 'p' 
  Name(6) = 'Bx' 
  Name(7) = 'By' 
  Name(8) = 'Bz' 
  Name(9) = 'psi' 
  Name(10) = 'alpha' 
  Name(11) = 'beta^1' 
  Name(12) = 'beta^2' 
  Name(13) = 'beta^3' 
  Name(14) = 'gamma_11' 
  Name(15) = 'gamma_12' 
  Name(16) = 'gamma_13' 
  Name(17) = 'gamma_22' 
  Name(18) = 'gamma_23' 
  Name(19) = 'gamma_33' 
#endif 
    !
#ifdef GPR3D 
    Name(1)  = 'rho'
    Name(2)  = 'u'
    Name(3)  = 'v'
    Name(4)  = 'w'
    Name(5)  = 'p'
    Name(6)  = 'A11'
    Name(7)  = 'A12'
    Name(8)  = 'A13'
    Name(9)  = 'A21'
    Name(10) = 'A22'
    Name(11) = 'A23'
    Name(12) = 'A31'
    Name(13) = 'A32'
    Name(14) = 'A33'
    Name(15) = 'j1'
    Name(16) = 'j2'
    Name(17) = 'j3'
#endif 
    !
END SUBROUTINE PDEVarName
    
SUBROUTINE PDEParName(Name) 
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE       
    ! Argument list 
    CHARACTER(LEN=10), INTENT(INOUT) :: Name(nParam)
    !
#ifdef EULER 
    !
#endif 
    !
#ifdef ELASTICITY 
    !Name(1) = 'lambda' 
    !Name(2) = 'mu' 
    !Name(3) = 'rho'
#endif 
    !
END SUBROUTINE PDEParName

!
! This subroutine checks if the state vector fulfills certain physical positivity criteria (e.g. for density and pressure)
! It can also hand back a "cured" state vector, where some floor values are used. However, this is to be considered as a 
! measure of last resort, and belongs to the "dirty trick" toolbox, which should never be opened, unless you are really desperate!    
!
SUBROUTINE PDEAssurePositivity(iErr,Q,Qout)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)     :: Q(nVar) 
    REAL, INTENT(OUT)    :: Qout(nVar) 
    INTEGER, INTENT(OUT) :: iErr
    ! Local variables 
    REAL :: p, irho, lam, mu 
    !
    iErr = 0 
    Qout = Q ! we do not apply dirty tricks here at this point 
    !
#ifdef EULER     
    !
    ! 3D compressible Euler equations 
    !
    p = (EQN%gamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
    ! 
    IF( Q(1).LE.1e-12) THEN
        iErr = -1 
    ENDIF 
    IF(p.LE.1e-12) THEN
        iErr = -2 
    ENDIF    
    !
#endif
    !
#ifdef ELASTICITY
    IF(Q(13)<-1e-12) THEN
        iErr = -1
    ENDIF
    IF(Q(13)>1.0+1e-12) THEN
        iErr = -2
    ENDIF
    IF(Q(14)<-1e-12) THEN
        iErr = -3
    ENDIF
    IF(Q(14)>1.0+1e-12) THEN
        iErr = -4
    ENDIF    
#endif     
    !           
END SUBROUTINE PDEAssurePositivity    
    
SUBROUTINE SourceTimeFunction(sigma,t,waveform)
    USE typesDef 
    IMPLICIT NONE
    ! Argument declaration
    INTEGER, INTENT(IN) :: waveform  
    REAL, INTENT(IN)  :: t
    REAL, INTENT(OUT) :: sigma(nVar)  
    ! Local variables 
    REAL :: a1, a2, tD 
    REAL :: Pi 
    sigma = 0. 
#ifdef ELASTICITY 
    Pi = ACOS(-1.0) 
    SELECT CASE(waveform)
    CASE(1)
        !sigma(7) = t*t !
        sigma(8) = SIN(2*Pi*t/0.5) 
    CASE(2)
        sigma(8) = EXP(-0.5*t**2/0.1**2)   
    CASE(3) 
        !a1 = -1.0       
        !a2 = -1000.0    
        !tD = 0.0      
        !sigma(8) = a1*(0.5+a2*(t-tD)**2)*EXP(a2*(t-tD)**2)        
        a1 = -2000.0/2200.  
        a2 = -(Pi*14.5)**2   
        tD = 0.08  
        sigma(8) = a1*(0.5+a2*(t-tD)**2)*EXP(a2*(t-tD)**2)        
    END SELECT 
#endif 
!        a1 = -1 !-2000.0/2200.0 ! -1.0 
!        a2 = -50. !-(Pi*14.5)**2 !-50.0  
!        tD = 0.08  ! 0.1 
!        sigma(1) = a1*(0.5+a2*(t-tD)**2)*EXP(a2*(t-tD)**2)        
END SUBROUTINE SourceTimeFunction     
    
    
SUBROUTINE PDEFusedSrcNCP(Src_BgradQ,Q,gradQ,par,time)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d), par(nParam), time 
    REAL, INTENT(OUT) :: Src_BgradQ(nVar) 
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,ip,iq,ii,jj,kk,ll,mm,iErr,count    
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar), BgradQ(nVar), src(nVar)  
    REAL :: k1,k2,k3,fff,ggg,e,c,ds,xi,sk,sknl,bs,g_cov(3,3),g_contr(3,3),dgup(3,3,3)
    REAL :: det, alpha, fa, k0, dk0(3), beta0(3), b0(3), u(3), ialpha 
    REAL :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar) 
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL :: v2,vf(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim 
    REAL :: VGRMHD(19), QGRMHD(19), gradQGRMHD(19,d), BgradQGRMHD(19), eta, itau, Ham, Mom(3), dKex(3,3,3), ov(3)     
    REAL :: Christoffel(3,3,3), RiemannNCP(3,3,3,3), RiemannSrc(3,3,3,3), dChristoffelNCP(3,3,3,3), dChristoffelSrc(3,3,3,3), DD(3,3,3), dDD(3,3,3,3), dChristoffelOrd(3,3,3,3)   
    REAL :: dChristoffel_tildeNCP(3,3,3,3), dChristoffel_tildeSrc(3,3,3,3) 
    REAL :: R, AA(3), dAA(3,3), BB(3,3), dBB(3,3,3), beta(3), Kex(3,3), Kmix(3,3), Kup(3,3), Z(3), dZ(3,3), nablaZNCP(3,3), nablaZSrc(3,3), RplusNablaZNCP, RplusNablaZSrc  
    REAL :: Theta, dTheta(3), nablaijalphaNCP(3,3), nablaijalphaSrc(3,3), Ricci(3,3), RicciNCP(3,3), RicciSrc(3,3), dtraceK(3), dtraceKNCP(3), dKtempNCP(3), dZNCP(3,3), dZSrc(3,3)  
    REAL :: dtgamma(3,3), dtK(3,3), dK(3,3,3), dtTheta, dtZ(3), dtalpha, dtGhat(3), dtbeta(3), dtbb(3), dtA(3), dtB(3,3), dtD(3,3,3)  
    REAL :: Aupdown, Aex(3,3), dAex(3,3,3), Amix(3,3), Aup(3,3), Ghat(3), Gtilde(3), dGhat(3,3), traceK, Kupdown, phi, dphi(3), PP(3), dPP(3,3), Pup(3), DDcontr(3) 
    REAL :: dGtildeSrc(3,3), dGtildeNCP(3,3), RiccitildeNCP(3,3), RicciphiNCP(3,3), RiccitildeSrc(3,3), RicciphiSrc(3,3)  
    REAL :: Christoffel_tilde(3,3,3), Christoffel_kind1(3,3,3), Zup(3), RicciPlusNablaZNCP(3,3), RicciPlusNablaZSrc(3,3), traceA, traceB, QG(3), b(3), faa, temp   
    REAL :: SecondOrderTermsNCP(3,3), SecondOrderTermsSrc(3,3), traceNCP, traceSrc, dtphi, dtTraceK, dtP(3)    
    REAL :: RNCP, RSrc, RPlusTwoNablaZNCP, RPlusTwoNablaZSrc, nablanablaalpha, nablanablaalphaNCP, nablanablaalphaSrc, Riemann(3,3,3,3), dChristoffel(3,3,3,3) 
    REAL :: TwoNablaZNCP, TwoNablaZSrc,  divAupNCP(3), divAupSrc(3), XX(3), dXX(3,3), nablaXNCP(3,3), nablaXSrc(3,3), dtX(3)  
    ! Matter contributions
    REAL :: sm(3), Sij(3,3), Sij_contr(3,3), sm_contr(3), S, tau, dens, bv_cov(3), sqrdet 
    REAL :: srctraceK, srcTheta
    REAL :: SrcK(3,3), SrcGhat(3)  

    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    !
    BgradQ = 0.0 
    !
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)

    k1   = EQN%CCZ4k1  
    k2   = EQN%CCZ4k2  
    k3   = EQN%CCZ4k3   
    fff  = EQN%CCZ4f 
    ggg  = EQN%CCZ4g 
    eta  = EQN%CCZ4eta 
    itau = EQN%CCZ4itau  
    e    = EQN%CCZ4e 
    c    = EQN%CCZ4c 
    mu   = EQN%CCZ4mu 
    ds   = EQN%CCZ4ds 
    bs   = EQN%CCZ4bs 
    xi   = EQN%CCZ4xi 
    sk   = EQN%CCZ4sk
    !
    ! These are the tilde quantities, so be careful !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    ! This determinant should be unity, since we use the conformal decomposition 
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    !g_cov = det**(-1./3.) * g_cov 
    !det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
    
    g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
    g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
    g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
    g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
    g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
    g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
    g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
    g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
    g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
      
    alpha = EXP(MAX(-20.,MIN(20.,Q(17)))) 
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 
    ! 
    K0    = Q(59)
    dK0   = sk*gradQ(59,:) 
    !  
    Aex(1,1) = Q(7) 
    Aex(1,2) = Q(8) 
    Aex(1,3) = Q(9) 
    Aex(2,1) = Q(8) 
    Aex(2,2) = Q(10) 
    Aex(2,3) = Q(11) 
    Aex(3,1) = Q(9) 
    Aex(3,2) = Q(11) 
    Aex(3,3) = Q(12) 
    !
    traceA = SUM(g_contr*Aex) 
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = MATMUL(g_contr, Aex)
    Aup  = MATMUL(g_contr, TRANSPOSE(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    b = Q(21:23) 
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   =  EXP(MAX(-20.,MIN(20.,Q(55))))  
    dphi  = gradQ(55,:) 
    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
    !
    dBB = dBB*sk    
    ! 
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi**2 + 1./3.*traceK*g_cov/phi**2 
    Kmix = MATMUL( phi**2*g_contr, Kex  ) 
    Kup  = MATMUL( phi**2*g_contr, TRANSPOSE(Kmix) ) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       !Christoffel_kind1(i,j,k) =  DD(i,j,k)+DD(j,i,k)-DD(k,i,j)     ! this definition does not work ! 
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
          !Gtilde(i)                = Gtilde(i)+2*g_contr(i,j)*g_contr(k,l)*DD(l,j,k) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    !
    DO i = 1, 3
     DO j = 1, 3
      DO l = 1, 3
          Gtilde(i) = Gtilde(i) + g_contr(j,l)*Christoffel_tilde(j,l,i) 
      ENDDO
     ENDDO     
    ENDDO    
    !
    Z   = 0.5*MATMUL( g_cov, Ghat - Gtilde ) 
    Zup = MATMUL(phi**2*g_contr, Z) 
    !
    !dChristoffel    = 0.0 
    dChristoffelNCP = 0.0
    dChristoffelSrc = 0.0 
    dChristoffel_tildeNCP = 0.0
    dChristoffel_tildeSrc = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          DO l = 1, 3 
            !dChristoffel(k,i,ip,m) = dChristoffel(k,i,ip,m) + g_contr(m,l)*(0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)))         & 
            !                                                - g_contr(m,l)*(g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)))  & 
            !                                                 +dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            ! 
            dChristoffelNCP(k,i,ip,m) = dChristoffelNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )         & 
                                                                  - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)) ) 
            !
            dChristoffel_tildeNCP(k,i,ip,m) = dChristoffel_tildeNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) ) 
            ! 
            !dChristoffelOrd(k,i,ip,m) = dChristoffelOrd(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)-dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)-dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)-dDD(l,k,i,ip)) )         & 
            !                                                      - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)-dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)-dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)-dPP(l,k)) )             
            !
            dChristoffelSrc(k,i,ip,m) = dChristoffelSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            !
            dChristoffel_tildeSrc(k,i,ip,m) = dChristoffel_tildeSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) 
            !             
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    !Riemann = 0.0 
    RiemannSrc = 0.0 
    RiemannNCP = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          !Riemann(i,k,ip,m)    = dChristoffel(k,i,ip,m)-dChristoffel(ip,i,k,m)
          RiemannNCP(i,k,ip,m) = dChristoffelNCP(k,i,ip,m)-dChristoffelNCP(ip,i,k,m)
          RiemannSrc(i,k,ip,m) = dChristoffelSrc(k,i,ip,m)-dChristoffelSrc(ip,i,k,m) 
          DO j = 1, 3
           !Riemann(i,k,ip,m)    = Riemann(i,k,ip,m)    + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
           RiemannSrc(i,k,ip,m) = RiemannSrc(i,k,ip,m) + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    !Ricci = 0.0 
    RicciNCP = 0.0 
    RicciSrc = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         !Ricci(m,n) = Ricci(m,n) + Riemann(m,l,n,l)  
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l)  
         RicciSrc(m,n) = RicciSrc(m,n) + RiemannSrc(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO    
    !
    !R    = phi**2*SUM(g_contr*Ricci) 
    RNCP = phi**2*SUM(g_contr*RicciNCP) 
    RSrc = phi**2*SUM(g_contr*RicciSrc) 
    !
    dGtildeNCP = 0.0
    dGtildeSrc = 0.0
    !DO k = 1, 3 
    ! DO j = 1, 3 
    !  DO l = 1, 3 
    !   DO m = 1, 3
    !    DO n = 1, 3 
    !     dGtildeNCP(j,k) = dGtildeNCP(j,k) + 2.0*g_contr(k,l)*g_contr(m,n)*0.5*(dDD(j,n,l,m)+dDD(n,j,l,m)) 
    !     dGtildeSrc(j,k) = dGtildeSrc(j,k) + 2.0*( g_contr(m,n)*DD(n,l,m)*dgup(j,k,l) + g_contr(k,l)*DD(n,l,m)*dgup(j,m,n) )         
    !    ENDDO
    !   ENDDO 
    !  ENDDO 
    ! ENDDO 
    !ENDDO 
    !    
    ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible... 
    ! 
    DO i = 1, 3 
     DO k = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dGtildeSrc(k,i) = dGtildeSrc(k,i) + dgup(k,j,l)*Christoffel_tilde(j,l,i) + g_contr(j,l)*dChristoffel_tildeSrc(k,j,l,i) 
           dGtildeNCP(k,i) = dGtildeNCP(k,i)                                        + g_contr(j,l)*dChristoffel_tildeNCP(k,j,l,i) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    ! If you want the original computation of the Ricci tensor according to the CCZ4 paper of Alic et al. 2012, use the following version. 
    ! By default, however, we compute the Ricci tensor ab definitionem from the Riemann tensor and the Christoffel symbols. 
    !
    !!DO j = 1, 3 
    !! DO i = 1, 3 
    !!    RiccitildeNCP(i,j) = 0 
    !!    DO l = 1, 3
    !!     DO m = 1, 3
    !!        RiccitildeNCP(i,j) = RiccitildeNCP(i,j)-g_contr(l,m)*0.5*(dDD(l,m,i,j)+dDD(m,l,i,j))
    !!     ENDDO
    !!    ENDDO 
    !!    DO k = 1, 3 
    !!        RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGtildeNCP(j,k)+g_cov(k,j)*dGtildeNCP(i,k)) 
    !!        !RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGhat(j,k)+g_cov(k,j)*dGhat(i,k)) 
    !!    ENDDO
    !! ENDDO
    !!ENDDO    
    !!    
    !!DO j = 1, 3 
    !! DO i = 1, 3         
    !!    RiccitildeSrc(i,j) = 0
    !!    DO k = 1, 3 
    !!      RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + 0.5*(g_cov(k,i)*dGtildeSrc(j,k)+g_cov(k,j)*dGtildeSrc(i,k)) 
    !!      RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + Ghat(k)*0.5*(Christoffel_kind1(i,j,k)+Christoffel_kind1(j,i,k)) 
    !!      DO l = 1, 3 
    !!       DO m = 1, 3 
    !!        RiccitildeSrc(i,j) = RiccitildeSrc(i,j)+g_contr(l,m)*(Christoffel_tilde(l,i,k)*Christoffel_kind1(j,k,m)+Christoffel_tilde(l,j,k)*Christoffel_kind1(i,k,m)+Christoffel_tilde(i,m,k)*Christoffel_kind1(k,j,l))
    !!       ENDDO
    !!      ENDDO
    !!    ENDDO
    !! ENDDO
    !!ENDDO
    !!
    !!
    !!DO j = 1, 3 
    !! DO i = 1, 3 
    !!   RicciphiNCP(i,j) = 0.5*(dPP(i,j)+dPP(j,i))
    !!   DO k = 1, 3 
    !!    DO l = 1, 3 
    !!       RicciphiNCP(i,j) = RicciphiNCP(i,j)+g_cov(i,j)*g_contr(l,k)*0.5*(dPP(k,l)+dPP(l,k))
    !!    ENDDO
    !!   ENDDO 
    !! ENDDO
    !!ENDDO
    !!
    !!Pup = MATMUL(g_contr,PP) 
    !!DO m = 1, 3
    !!    DDcontr(m) = SUM(g_contr*DD(:,m,:)) 
    !!ENDDO 
    !!
    !!DO i = 1, 3 
    !! DO j = 1, 3 
    !!    RicciphiSrc(i,j) = PP(i)*PP(j) 
    !!    DO k = 1, 3 
    !!     RicciphiSrc(i,j) = RicciphiSrc(i,j) - Christoffel_tilde(i,j,k)*PP(k) - 2*g_cov(i,j)*DDcontr(k)*Pup(k)       
    !!     DO l = 1, 3 
    !!        RicciphiSrc(i,j) = RicciphiSrc(i,j) - g_cov(i,j)*g_contr(l,k)*PP(l)*PP(k) !-g_contr(k,l)*PP(k)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))   
    !!     ENDDO
    !!    ENDDO
    !! ENDDO
    !!ENDDO    
    !
    dZNCP = 0.0 
    dZSrc = 0.0 
    DO i = 1, 3
     DO k = 1, 3    
      DO j = 1, 3 
        dZNCP(k,i) = dZNCP(k,i) + 0.5*g_cov(i,j)*( dGhat(k,j)-dGtildeNCP(k,j) )  
        dZSrc(k,i) = dZSrc(k,i) + DD(k,i,j)*(Ghat(j)-Gtilde(j)) + 0.5*g_cov(i,j)*( -dGtildeSrc(k,j) )
       ENDDO 
      ENDDO 
    ENDDO     
    !
    nablaZNCP = dZNCP 
    nablaZSrc = 0.0 
    DO j = 1, 3 
     DO i = 1, 3 
      nablaZSrc(i,j) = dZSrc(i,j)
      DO k = 1, 3 
        nablaZSrc(i,j) = nablaZSrc(i,j) - Christoffel(i,j,k)*Z(k) 
      ENDDO 
     ENDDO
    ENDDO    
    !
    !RicciPlusNablaZNCP = RiccitildeNCP + RicciphiNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    !RicciPlusNablaZSrc = RiccitildeSrc + RicciphiSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RicciPlusNablaZNCP = RicciNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    RicciPlusNablaZSrc = RicciSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RPlusTwoNablaZNCP = phi**2*SUM(g_contr*RicciPlusNablaZNCP) 
    RPlusTwoNablaZSrc = phi**2*SUM(g_contr*RicciPlusNablaZSrc) 
    !
    !Kupdown = SUM(Kex*Kup) 
    !temp = RPlusTwoNablaZNCP + RPlusTwoNablaZSrc - KupDown + traceK**2 
    !temp = SUM(phi**2*g_contr*Ricci) - KupDown + traceK**2 
    !
    nablaijalphaNCP = 0.0
    nablaijalphaSrc = 0.0
    DO j = 1, 3 
     DO i = 1, 3 
       nablaijalphaNCP(i,j) = alpha*0.5*( dAA(i,j)+dAA(j,i) ) 
       nablaijalphaSrc(i,j) = alpha*AA(j)*AA(i) 
       DO k = 1, 3 
         nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j) - Christoffel(i,j,k)*alpha*AA(k)  
       ENDDO
     ENDDO
    ENDDO 
    nablanablaalphaNCP = phi**2*SUM( g_contr*nablaijalphaNCP ) 
    nablanablaalphaSrc = phi**2*SUM( g_contr*nablaijalphaSrc ) 
!    nablanablaalpha = 0.0
!    DO j = 1, 3 
!     DO i = 1, 3
!         nablanablaalpha = nablanablaalpha + phi**2*( g_contr(i,j)*AA(j)*AA(i)*alpha + alpha*g_contr(i,j)*dAA(i,j) - 2*g_contr(i,j)*SUM(g_contr(:,:)*TRANSPOSE(DD(:,j,:))*AA(i)*alpha)- g_contr(i,j)*PP(i)*AA(j)*alpha )  
!     ENDDO
!    ENDDO 
    !
    SecondOrderTermsNCP = -nablaijalphaNCP + alpha*RicciPlusNablaZNCP 
    SecondOrderTermsSrc = -nablaijalphaSrc + alpha*RicciPlusNablaZSrc 
    traceNCP = SUM( g_contr*SecondOrderTermsNCP ) 
    SecondOrderTermsNCP = SecondOrderTermsNCP - 1./3.*g_cov*traceNCP 
    traceSrc = SUM( g_contr*SecondOrderTermsSrc ) 
    SecondOrderTermsSrc = SecondOrderTermsSrc - 1./3.*g_cov*traceSrc 
    !
    ! Now assemble all this terrible stuff... 
    !
    dtgamma = - 2*alpha*Aex - itau*(det-1.0)*g_cov 
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3
          dtgamma(i,j) = dtgamma(i,j) + g_cov(k,i)*BB(j,k) + g_cov(k,j)*BB(i,k) - 2./3.*g_cov(i,j)*BB(k,k) + beta(k)*2*DD(k,i,j) 
      ENDDO
     ENDDO
    ENDDO 
    !
    ! Main variables of the CCZ4 system 
    dtK = phi**2*SecondOrderTermsNCP + beta(1)*dAex(1,:,:) + beta(2)*dAex(2,:,:) + beta(3)*dAex(3,:,:)      ! extrinsic curvature
    dtK = dtK + phi**2*SecondOrderTermsSrc + alpha*Aex*(traceK-2*Theta) - 2*alpha*MATMUL(Aex,Amix) - itau*g_cov*traceA 
    DO j = 1, 3 
     DO i = 1, 3 
      DO k = 1, 3 
         dtK(i,j) = dtK(i,j) + Aex(k,i)*BB(j,k) + Aex(k,j)*BB(i,k) - 2./3.*Aex(i,j)*BB(k,k) 
      ENDDO
     ENDDO
    ENDDO 
    !
    dtTraceK = - nablanablaalphaNCP - nablanablaalphaSrc + alpha*( RPlusTwoNablaZNCP + RPlusTwoNablaZSrc + traceK**2 - 2*Theta*traceK ) - 3*alpha*k1*(1+k2)*Theta + SUM(beta(:)*dtraceK(:)) 
    !
    traceB = BB(1,1) + BB(2,2) + BB(3,3) 
    dtphi   = beta(1)*PP(1) + beta(2)*PP(2) + beta(3)*PP(3) + 1./3.*alpha*traceK - 1./3.*traceB 
    dtalpha = -alpha*fa*(traceK-K0-c*2*Theta) + beta(1)*AA(1) + beta(2)*AA(2) + beta(3)*AA(3) 

    Aupdown = SUM(Aex*Aup) 
    ! *** original 
    dtTheta =  0.5*alpha*e**2*(RplusTwoNablaZNCP + RplusTwoNablaZSrc) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)    &            ! temporal Z 
             + 0.5*alpha*e**2*( - Aupdown + 2./3.*traceK**2 ) - alpha*Theta*traceK - SUM(Zup*alpha*AA) - alpha*k1*(2+k2)*Theta  
    ! *** use turbo cleaning here, i.e. without multiplying with alpha *** 
    !dtTheta = 0.5*e**2*( RplusTwoNablaZNCP + RplusTwoNablaZSrc ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)    &                ! temporal Z 
    !         + 0.5*e**2*( - Aupdown + 2./3.*traceK**2 ) - c*alpha*Theta*traceK - SUM(Zup*alpha*AA) - k2*Theta 
    !
    divAupNCP = 0.0
    divAupSrc = 0.0 
    DO i = 1, 3
        DO j = 1, 3
         DO l = 1, 3
          DO k = 1, 3    
            divAupNCP(i) = divAupNCP(i) + g_contr(i,l)*g_contr(j,k)*dAex(j,l,k) 
            divAupSrc(i) = divAupSrc(i) + ( dgup(j,i,l)*g_contr(j,k) + g_contr(i,l)*dgup(j,j,k) )*Aex(l,k) 
          ENDDO
         ENDDO
        ENDDO        
    ENDDO 
    ! 
    DO i = 1, 3 
        Mom(i) = SUM(Christoffel_tilde(:,:,i)*Aup(:,:)) - 3.0*SUM(Aup(i,:)*PP(:)) - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) + divAupNCP(i) + divAupSrc(i)  
    ENDDO 
    !
    !!dKex = 0.0
    !!DO j = 1, 3
    !! DO i = 1, 3
    !!  DO k = 1, 3
    !!      dKex(k,i,j) = 1.0/phi**2*( dAex(k,i,j) - 2.0*Aex(i,j)*PP(k) + 1./3.*dtraceK(k)*g_cov(i,j) + 2./3.*traceK*DD(k,i,j) - 2./3.*traceK*g_cov(i,j)*PP(k) ) 
    !!  ENDDO
    !! ENDDO
    !!ENDDO 
    !!!
    !!Mom(:) = 0.0
    !!DO ii = 1, 3
    !!    DO jj = 1, 3
    !!        DO ll = 1, 3
    !!            Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*(dKex(ll,ii,jj) - dKex(ii,jj,ll))
    !!            DO mm = 1, 3
    !!                Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*( - Christoffel(jj,ll,mm)*Kex(mm,ii) + Christoffel(jj,ii,mm)*Kex(mm,ll)) 
    !!            ENDDO
    !!        ENDDO     
    !!    ENDDO
    !!ENDDO   
    !!Mom = phi**2*MATMUL( g_contr, Mom ) 
    !
    Kupdown = SUM(Kex*Kup) 
    Ham = RPlusTwoNablaZNCP + RPlusTwoNablaZSrc - KupDown + traceK**2     
    !
    DO i = 1, 3
        dtGhat(i) = 2*alpha*( SUM(Christoffel_tilde(:,:,i)*Aup(:,:)) - 3.0*SUM(Aup(i,:)*PP(:)) - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) )    &
        !dtGhat(i)  =  2*alpha*( -divAupNCP(i) - divAupSrc(i) + ds**2*Mom(i) )                           & 
                    + 2*alpha*SUM( g_contr(:,i)*( dTheta(:) - Theta*AA(:) - 2./3.*traceK*Z(:)  ) )  & 
                    - 2*SUM( Aup(i,:)*alpha*AA(:) ) + 2./3.*Gtilde(i)*traceB - 2*alpha*k1*SUM(g_contr(i,:)*Z(:)) - SUM(Gtilde(:)*BB(:,i))   &
                   + beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i) 
        DO l = 1, 3
         DO k = 1, 3
             dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) + & 
                                     2*k3*( 2./3.*g_contr(i,l)*Z(l)*BB(k,k) - g_contr(l,k)*Z(l)*BB(k,i) ) 
         ENDDO
        ENDDO         
    ENDDO 
    DO k = 1, 3 
        ov(k) = + 2*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )   ! here we can use the constraint that trace A tilde = 0.         
    ENDDO        
    dtGhat = dtGhat + sk*MATMUL(g_contr,ov)                                               ! the above ordering constraint is "down", so we must raise the index via g_contr. 
    !
    dtbb = xi*dtGhat - eta*b                                                              !  <= be careful, this damping term -eta*b may be dangerous for the gamma driver, since it may kill the waves that you want !  
    ! Add the following terms if you want shift convection in the PDE for b^i 
    dtbb = dtbb + bs*( - beta(1)*gradQ(14:16,1) - beta(2)*gradQ(14:16,2) - beta(3)*gradQ(14:16,3)  + beta(1)*gradQ(21:23,1) + beta(2)*gradQ(21:23,2) + beta(3)*gradQ(21:23,3)  )   !         
    dtbb = sk*dtbb 
    !
    QG = ggg*MATMUL( g_contr, AA - Z ) 
    DO j = 1, 3 
     DO i = 1, 3 
        DO n = 1, 3
         DO l = 1, 3 
           QG(i) = QG(i) + ggg*g_contr(i,j)*(-2*g_contr(n,l)*DD(l,n,j)) 
         ENDDO
        ENDDO
     ENDDO
    ENDDO 
    !            
    dtbeta  = + fff*b  
    ! Add the following term if you want to have shift convection in the PDE for beta^i 
    ! Do not add it if you want a real Lie derivative for beta. In this case, the advection term cancels out. 
    dtbeta = dtbeta + bs*( beta(1)*BB(1,:) + beta(2)*BB(2,:) + beta(3)*BB(3,:) )      
    dtbeta = sk*dtbeta 
    !
    ! Auxiliary variables 
    dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:) - alpha*AA*(fa+alpha*faa)*(traceK-K0-c*2*Theta) + MATMUL(BB, AA) 
    DO k = 1, 3 
        dtA(k) = dtA(k) - sk*alpha*fa*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )   ! here we can use the constraint that trace A tilde = 0. 
    ENDDO    
    ! 
    ! In CCZ4 we have completely removed all the conservative fluxes. 
    dtB(:,1) = fff*gradQ(21,:)  
    dtB(:,2) = fff*gradQ(22,:)  
    dtB(:,3) = fff*gradQ(23,:)  
    ! #ordB2# 
    ! for the ordering constraints, we have to check whether they should be multiplied by alpha**2, or not... 
    DO i = 1, 3 
     DO k = 1, 3     
      DO j = 1, 3 
         !dtB(k,i) = +mu*alpha*g_contr(i,j)*( dZNCP(k,j) + dZSrc(k,j) ) + mu*dgup(k,i,j)*Z(j)      
         dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( (dPP(k,j)-dPP(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            !dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( -2*g_contr(n,l)*0.5*(dDD(k,l,n,j)-dDD(l,k,n,j)) ) 
            dtB(k,i) = dtB(k,i) - mu*alpha**2*g_contr(i,j)*g_contr(n,l)*( dDD(k,l,j,n)-dDD(l,k,j,n) )   
          ENDDO 
         ENDDO 
      ENDDO 
      !
      ENDDO 
    ENDDO 
    !
    dtB = dtB + bs*( beta(1)*dBB(1,:,:) + beta(2)*dBB(2,:,:) + beta(3)*dBB(3,:,:) + MATMUL(BB,BB) ) 
    ! 
    ! New stuff 2, which makes appear a Lie derivative and a second order ordering constraint 
    !DO i = 1, 3
    ! DO k = 1, 3 
    !   dtB(k,i) = dtB(k,i) + ( beta(1)*dBB(1,k,i) + beta(2)*dBB(2,k,i) + beta(3)*dBB(3,k,i) - beta(1)*dBB(k,1,i) - beta(2)*dBB(k,2,i) - beta(3)*dBB(k,3,i) )    !  + MATMUL(BB,BB) 
    ! ENDDO
    !ENDDO 
    dtB = dtB*sk 
    !
    dtD = -alpha*dAex  
    DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            dtD(k,i,j) = dtD(k,i,j) + ( 0.5*(g_cov(m,i)*0.5*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./3.*g_cov(i,j)*0.5*(dBB(k,m,m)+dBB(m,k,m)) ) 
            DO n = 1, 3
                dtD(k,i,j) = dtD(k,i,j) + 1./3*alpha*g_cov(i,j)*g_contr(n,m)*dAex(k,n,m) + 1./3.*alpha*g_cov(i,j)*dgup(k,n,m)*Aex(n,m)   ! explicitly remove the trace of tilde A again 
            ENDDO 
        ENDDO        
       ENDDO
      ENDDO
    ENDDO 
    ! 
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3 
        dtD(k,i,j) = dtD(k,i,j) - alpha*AA(k)*Aex(i,j) 
        DO l = 1, 3 
          dtD(k,i,j) = dtD(k,i,j) + BB(k,l)*DD(l,i,j) + DD(k,l,i)*BB(j,l) + DD(k,l,j)*BB(i,l) - 2./3.*DD(k,i,j)*BB(l,l) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO         
    !
    dtD = dtD + beta(1)*dDD(1,:,:,:) + beta(2)*dDD(2,:,:,:) + beta(3)*dDD(3,:,:,:)
    !
    dtP = MATMUL(BB,PP) + beta(1)*dPP(1,:) + beta(2)*dPP(2,:) + beta(3)*dPP(3,:)    
    DO k = 1, 3 
     dtP(k) = dtP(k) + 1./3.*alpha*dtraceK(k) + 1./3.*alpha*AA(k)*traceK + sk*1./3.*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )  
     DO i = 1, 3 
          dtP(k) = dtP(k) - 1./3.*0.5*(dBB(k,i,i)+dBB(i,k,i)) 
     ENDDO
    ENDDO 
    !
#ifdef CCZ4GRMHD 
    ! 
    ! Add algebraic matter source terms to CCZ4 
    sqrdet   = phi**(-6.0/2.0) ! square root of determinant
    dens     = QGRMHD(1)   / sqrdet
    sm       = QGRMHD(2:4) / sqrdet
    sm_contr = MATMUL(phi**2*g_contr,sm) 
    tau      = QGRMHD(5)   / sqrdet
    Bv_cov   = QGRMHD(6:8) / sqrdet
    Bv_contr = MATMUL(phi**2*g_contr,Bv_cov)
    E        = tau + dens            ! U = E
    ! prepare the PDECons2PrimGRMHD.
    QGRMHD(1:9)   = Q(60:68)      ! hydro variables
    QGRMHD(10)    = alpha         ! lapse
    QGRMHD(11:13) = Q(18:20)      ! shift
    QGRMHD(14:19) = Q(1:6)/phi**2 ! metric (without tilde).
    CALL PDECons2PrimGRMHD(VGRMHD,QGRMHD,iErr)
    rho    = VGRMHD(1)
    vf_cov = VGRMHD(2:4)
    p      = VGRMHD(5)
    ! compute the lorentz factor
    vf       = MATMUL(phi**2*g_contr,vf_cov)
    v2       = SUM(vf*vf_cov)
    lf       = 1.0/sqrt(1.0 - v2)
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    w      = rho + gamma1*p   ! rho*hentalpy
    ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
    b2     = SUM(Bv_cov*Bv_contr)**2/lf**2 + SUM(Bv_contr*Bv_cov)**2 ! magnetic pressure
    !
    if(rho>0.9*1e-2) then
        continue
    endif 
    !
    DO j = 1, 3
      DO i = 1, 3
        Sij(i,j) = ww*vf_cov(i)*vf_cov(j) + (p + b2/2)*g_cov(i,j)/(phi**2) &
		 - Bv_cov(i)*Bv_cov(j)/lf**2 - SUM(Bv_contr*vf_cov)*vf_cov(i)*Bv_cov(j)
      ENDDO
    ENDDO
    !
    Sij_contr = MATMUL(phi**2*g_contr,MATMUL(phi**2*g_contr,Sij)) 
    S = SUM(phi**2*g_contr*Sij) ! Trace of Sij

    ! Compute Source contributions
    SrcK        = -(phi**2)*8*pi*alpha*( Sij - 1./3.*g_cov/phi**2 * S ) ! Tracefree
    dtK         = dtK + SrcK

    SrcTraceK   = 4*pi*alpha*(S - 3*E)
    dtTraceK    = dtTraceK + SrcTraceK
    
    SrcGhat     = -16*pi*alpha*MATMUL(g_contr, sm)
    dtGhat      = dtGhat + SrcGhat
    
    SrcTheta    = -8*pi*alpha*E
    dtTheta     = dtTheta + SrcTheta
    ! 
#endif 
    !
    BgradQ(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)          ! \tilde \gamma_ij 
    BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                                  ! \tilde A_ij 
    BgradQ(13)     = dtTheta                                                                                           ! Theta       
    BgradQ(14:16)  = dtGhat(1:3)                                                                                       ! \hat \Gamma^i           
    BgradQ(17)     = dtalpha                                                                                           ! log alpha 
    BgradQ(18:20)  = dtbeta                                                                                            ! beta^i 
    BgradQ(21:23)  = dtbb                                                                                              ! b^i 
    BgradQ(24:26)  = dtA(1:3)                                                                                          ! A_k       
    BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /)    ! B_k^i 
    BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                      ! D_kij 
    BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /)                      ! D_kij 
    BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)                      ! D_kij 
    BgradQ(54)     = dtTraceK                                                                                          ! traceK 
    BgradQ(55)     = dtphi                                                                                             ! log phi 
    BgradQ(56:58)  = dtP                                                                                               ! P_k 
    !
    Src_BgradQ = BgradQ    ! here, we do not have to change sign, since we work on the right hand side in the fused subroutine 
    !
#ifdef CCZ4GRMHD 
    !
    ! Add the metric terms to the GRMHD equations 
    !
    alpha = EXP(Q(17)) 
    phi   = EXP(Q(55)) 
    QGRMHD(1:9)   = Q(60:68)        ! hydro variables 
    QGRMHD(10)    = alpha           ! lapse 
    QGRMHD(11:13) = Q(18:20)        ! shift 
    QGRMHD(14:19) = Q(1:6)/phi**2   ! metric 
    !
    gradQGRMHD(1:9,:)   = 0.0   ! hydro variables (not needed) 
    gradQGRMHD(10,:)    = alpha*(/ Q(24), Q(25), Q(26) /)           ! lapse 
    gradQGRMHD(11,:)    = (/ Q(27), Q(28), Q(29) /)                 ! shift1 
    gradQGRMHD(12,:)    = (/ Q(30), Q(31), Q(32) /)                 ! shift2 
    gradQGRMHD(13,:)    = (/ Q(33), Q(34), Q(35) /)                 ! shift3 
    gradQGRMHD(14:19,1)  = 2.0/phi**2*( Q(36:41) - PP(1)*Q(1:6) )   ! metric 
    gradQGRMHD(14:19,2)  = 2.0/phi**2*( Q(42:47) - PP(2)*Q(1:6) )   ! metric 
    gradQGRMHD(14:19,3)  = 2.0/phi**2*( Q(48:53) - PP(3)*Q(1:6) )   ! metric 
    !
    CALL PDENCPGRMHD(BgradQGRMHD,QGRMHD,gradQGRMHD,par) 
    src_bgradQ(60:68) = -BgradQGRMHD(1:9)   ! Change sign, since the NCPGRMHD is on the left, but the fused stuff is on the right. 
    !
    src_bgradQ(64) = sqrdet*( alpha*SUM(Kex*Sij_contr) - SUM(sm_contr*aa*alpha) ) ! Energy with dynamical extr. curvature     
    !
    ! CCZ4end 
    !
    CONTINUE
    ! 
#endif      
    !
    RETURN
    !
#endif 
    !
    !
#if defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)

    k1   = EQN%CCZ4k1  
    k2   = EQN%CCZ4k2  
    k3   = EQN%CCZ4k3   
    fff  = EQN%CCZ4f 
    ggg  = EQN%CCZ4g 
    eta  = EQN%CCZ4eta 
    itau = EQN%CCZ4itau  
    e    = EQN%CCZ4e 
    c    = EQN%CCZ4c 
    mu   = EQN%CCZ4mu 
    ds   = EQN%CCZ4ds 
    bs   = EQN%CCZ4bs 
    xi   = EQN%CCZ4xi 
    sk   = EQN%CCZ4sk
    !
    ! These are the tilde quantities, so be careful !    
    g_cov(1,1) = Q(1)
    g_cov(1,2) = Q(2)
    g_cov(1,3) = Q(3)
    g_cov(2,1) = Q(2)
    g_cov(2,2) = Q(4)
    g_cov(2,3) = Q(5)
    g_cov(3,1) = Q(3)
    g_cov(3,2) = Q(5)
    g_cov(3,3) = Q(6)
    ! This determinant should be unity, since we use the conformal decomposition 
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    g_cov = det**(-1./3.) * g_cov 
    det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
    
    g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
    g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
    g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
    g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
    g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
    g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
    g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
    g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
    g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
      
    alpha = EXP(MAX(-20.,MIN(20.,Q(17)))) 
    SELECT CASE(EQN%CCZ4LapseType) 
    CASE(0)  ! harmonic 
        fa = 1.0 
        faa = 0.0 
    CASE DEFAULT  ! 1 + log 
        fa = 2.0/alpha
        faa = -2.0/alpha**2   
    END SELECT 
    ! 
    K0    = Q(62)
    dK0   = sk*gradQ(62,:) 
    !  
    Aex(1,1) = Q(7) 
    Aex(1,2) = Q(8) 
    Aex(1,3) = Q(9) 
    Aex(2,1) = Q(8) 
    Aex(2,2) = Q(10) 
    Aex(2,3) = Q(11) 
    Aex(3,1) = Q(9) 
    Aex(3,2) = Q(11) 
    Aex(3,3) = Q(12) 
    !
    traceA = SUM(g_contr*Aex) 
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = MATMUL(g_contr, Aex)
    Aup  = MATMUL(g_contr, TRANSPOSE(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    XX = (/ Q(59), Q(60), Q(61) /) 
    dXX(:,1) = gradQ(59,:) 
    dXX(:,2) = gradQ(60,:) 
    dXX(:,3) = gradQ(61,:) 
    !
    b = Q(21:23) 
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   =  EXP(MAX(-20.,MIN(20.,Q(55))))  
    dphi  = gradQ(55,:) 
    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
    !
    dBB = dBB*sk    
    ! 
    DD(1,1,1)=Q(36) 
    DD(1,1,2)=Q(37) 
    DD(1,1,3)=Q(38) 
    DD(1,2,1)=Q(37) 
    DD(1,2,2)=Q(39) 
    DD(1,2,3)=Q(40)
    DD(1,3,1)=Q(38) 
    DD(1,3,2)=Q(40) 
    DD(1,3,3)=Q(41)
    ! 
    DD(2,1,1)=Q(42) 
    DD(2,1,2)=Q(43) 
    DD(2,1,3)=Q(44) 
    DD(2,2,1)=Q(43) 
    DD(2,2,2)=Q(45) 
    DD(2,2,3)=Q(46)
    DD(2,3,1)=Q(44) 
    DD(2,3,2)=Q(46) 
    DD(2,3,3)=Q(47)
    !
    DD(3,1,1)=Q(48) 
    DD(3,1,2)=Q(49) 
    DD(3,1,3)=Q(50) 
    DD(3,2,1)=Q(49) 
    DD(3,2,2)=Q(51) 
    DD(3,2,3)=Q(52)
    DD(3,3,1)=Q(50) 
    DD(3,3,2)=Q(52) 
    DD(3,3,3)=Q(53)
    !
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi**2 + 1./3.*traceK*g_cov/phi**2 
    Kmix = MATMUL( phi**2*g_contr, Kex  ) 
    Kup  = MATMUL( phi**2*g_contr, Kmix ) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       !Christoffel_kind1(i,j,k) =  DD(i,j,k)+DD(j,i,k)-DD(k,i,j)     ! this definition does not work ! 
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
          !Gtilde(i)                = Gtilde(i)+2*g_contr(i,j)*g_contr(k,l)*DD(l,j,k) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    !
    DO i = 1, 3
     DO j = 1, 3
      DO l = 1, 3
          Gtilde(i) = Gtilde(i) + g_contr(j,l)*Christoffel_tilde(j,l,i) 
      ENDDO
     ENDDO     
    ENDDO    
    !
    Z   = 0.5*MATMUL( g_cov, Ghat - Gtilde ) 
    Zup = MATMUL(phi**2*g_contr, Z) 
    !
    !dChristoffel    = 0.0 
    dChristoffelNCP = 0.0
    dChristoffelSrc = 0.0 
    dChristoffel_tildeNCP = 0.0
    dChristoffel_tildeSrc = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          DO l = 1, 3 
            !dChristoffel(k,i,ip,m) = dChristoffel(k,i,ip,m) + g_contr(m,l)*(0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)))         & 
            !                                                - g_contr(m,l)*(g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)))  & 
            !                                                 +dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            ! 
            dChristoffelNCP(k,i,ip,m) = dChristoffelNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) )         & 
                                                                  - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)) ) 
            !
            dChristoffel_tildeNCP(k,i,ip,m) = dChristoffel_tildeNCP(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)) ) 
            ! 
            !dChristoffelOrd(k,i,ip,m) = dChristoffelOrd(k,i,ip,m) + g_contr(m,l)*( 0.5*(dDD(k,i,ip,l)-dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)-dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)-dDD(l,k,i,ip)) )         & 
            !                                                      - g_contr(m,l)*( g_cov(ip,l)*0.5*(dPP(k,i)-dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)-dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)-dPP(l,k)) )             
            !
            dChristoffelSrc(k,i,ip,m) = dChristoffelSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
            !
            dChristoffel_tildeSrc(k,i,ip,m) = dChristoffel_tildeSrc(k,i,ip,m) + dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) 
            !             
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    !Riemann = 0.0 
    RiemannSrc = 0.0 
    RiemannNCP = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          !Riemann(i,k,ip,m)    = dChristoffel(k,i,ip,m)-dChristoffel(ip,i,k,m)
          RiemannNCP(i,k,ip,m) = dChristoffelNCP(k,i,ip,m)-dChristoffelNCP(ip,i,k,m)
          RiemannSrc(i,k,ip,m) = dChristoffelSrc(k,i,ip,m)-dChristoffelSrc(ip,i,k,m) 
          DO j = 1, 3
           !Riemann(i,k,ip,m)    = Riemann(i,k,ip,m)    + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
           RiemannSrc(i,k,ip,m) = RiemannSrc(i,k,ip,m) + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    !Ricci = 0.0 
    RicciNCP = 0.0 
    RicciSrc = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         !Ricci(m,n) = Ricci(m,n) + Riemann(m,l,n,l)  
         RicciNCP(m,n) = RicciNCP(m,n) + RiemannNCP(m,l,n,l)  
         RicciSrc(m,n) = RicciSrc(m,n) + RiemannSrc(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO    
    !
    RNCP = phi**2*SUM(g_contr*RicciNCP) 
    RSrc = phi**2*SUM(g_contr*RicciSrc) 
    !
    dGtildeNCP = 0.0
    dGtildeSrc = 0.0
    !DO k = 1, 3 
    ! DO j = 1, 3 
    !  DO l = 1, 3 
    !   DO m = 1, 3
    !    DO n = 1, 3 
    !     dGtildeNCP(j,k) = dGtildeNCP(j,k) + 2.0*g_contr(k,l)*g_contr(m,n)*0.5*(dDD(j,n,l,m)+dDD(n,j,l,m)) 
    !     dGtildeSrc(j,k) = dGtildeSrc(j,k) + 2.0*( g_contr(m,n)*DD(n,l,m)*dgup(j,k,l) + g_contr(k,l)*DD(n,l,m)*dgup(j,m,n) )         
    !    ENDDO
    !   ENDDO 
    !  ENDDO 
    ! ENDDO 
    !ENDDO 
    !    
    ! Here we directly compute the derivative of Gtilde from its original definition as contracted Christoffel symbol,
    ! without assuming unit determinant of the conformal metric. Back to the roots, and as few assumptions as possible... 
    !
    DO i = 1, 3 
     DO k = 1, 3
      DO j = 1, 3
       DO l = 1, 3
           dGtildeSrc(k,i) = dGtildeSrc(k,i) + dgup(k,j,l)*Christoffel_tilde(j,l,i) + g_contr(j,l)*dChristoffel_tildeSrc(k,j,l,i) 
           dGtildeNCP(k,i) = dGtildeNCP(k,i)                                        + g_contr(j,l)*dChristoffel_tildeNCP(k,j,l,i) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    ! If you want the original computation of the Ricci tensor according to the CCZ4 paper of Alic et al. 2012, use the following version. 
    ! By default, however, we compute the Ricci tensor ab definitionem from the Riemann tensor and the Christoffel symbols. 
    !
    DO j = 1, 3 
     DO i = 1, 3 
        RiccitildeNCP(i,j) = 0 
        DO l = 1, 3
         DO m = 1, 3
            RiccitildeNCP(i,j) = RiccitildeNCP(i,j)-g_contr(l,m)*0.5*(dDD(l,m,i,j)+dDD(m,l,i,j))
         ENDDO
        ENDDO 
        DO k = 1, 3 
            RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGtildeNCP(j,k)+g_cov(k,j)*dGtildeNCP(i,k)) 
            !RiccitildeNCP(i,j) = RiccitildeNCP(i,j) + 0.5*(g_cov(k,i)*dGhat(j,k)+g_cov(k,j)*dGhat(i,k)) 
        ENDDO
     ENDDO
    ENDDO    
        
    DO j = 1, 3 
     DO i = 1, 3         
        RiccitildeSrc(i,j) = 0
        DO k = 1, 3 
          RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + 0.5*(g_cov(k,i)*dGtildeSrc(j,k)+g_cov(k,j)*dGtildeSrc(i,k)) 
          RiccitildeSrc(i,j) = RiccitildeSrc(i,j) + Ghat(k)*0.5*(Christoffel_kind1(i,j,k)+Christoffel_kind1(j,i,k)) 
          DO l = 1, 3 
           DO m = 1, 3 
            RiccitildeSrc(i,j) = RiccitildeSrc(i,j)+g_contr(l,m)*(Christoffel_tilde(l,i,k)*Christoffel_kind1(j,k,m)+Christoffel_tilde(l,j,k)*Christoffel_kind1(i,k,m)+Christoffel_tilde(i,m,k)*Christoffel_kind1(k,j,l))
           ENDDO
          ENDDO
        ENDDO
     ENDDO
    ENDDO
    
    
    DO j = 1, 3 
     DO i = 1, 3 
       RicciphiNCP(i,j) = 0.5*(dPP(i,j)+dPP(j,i))
       DO k = 1, 3 
        DO l = 1, 3 
           RicciphiNCP(i,j) = RicciphiNCP(i,j)+g_cov(i,j)*g_contr(l,k)*0.5*(dPP(k,l)+dPP(l,k))
        ENDDO
       ENDDO 
     ENDDO
    ENDDO
    
    Pup = MATMUL(g_contr,PP) 
    DO m = 1, 3
        DDcontr(m) = SUM(g_contr*DD(:,m,:)) 
    ENDDO 
    
    DO i = 1, 3 
     DO j = 1, 3 
        RicciphiSrc(i,j) = PP(i)*PP(j) 
        DO k = 1, 3 
         RicciphiSrc(i,j) = RicciphiSrc(i,j) - Christoffel_tilde(i,j,k)*PP(k) - 2*g_cov(i,j)*DDcontr(k)*Pup(k)       
         DO l = 1, 3 
            RicciphiSrc(i,j) = RicciphiSrc(i,j) - g_cov(i,j)*g_contr(l,k)*PP(l)*PP(k) !-g_contr(k,l)*PP(k)*(DD(i,j,l)+DD(j,i,l)-DD(l,i,j))   
         ENDDO
        ENDDO
     ENDDO
    ENDDO    
    !
    dZNCP = 0.0 
    dZSrc = 0.0 
    DO i = 1, 3
     DO k = 1, 3    
      DO j = 1, 3 
        dZNCP(k,i) = dZNCP(k,i) + 0.5*g_cov(i,j)*( dGhat(k,j)-dGtildeNCP(k,j) )  
        dZSrc(k,i) = dZSrc(k,i) + DD(k,i,j)*(Ghat(j)-Gtilde(j)) + 0.5*g_cov(i,j)*( -dGtildeSrc(k,j) )
       ENDDO 
      ENDDO 
    ENDDO     
    !
    nablaZNCP = dZNCP 
    nablaXNCP = dXX  
    nablaZSrc = 0.0 
    nablaXSrc = 0.0 
    DO j = 1, 3 
     DO i = 1, 3 
      nablaZSrc(i,j) = dZSrc(i,j)
      DO k = 1, 3 
        nablaZSrc(i,j) = nablaZSrc(i,j) - Christoffel(i,j,k)*Z(k) 
        nablaXSrc(i,j) = nablaXSrc(i,j) - Christoffel(i,j,k)*XX(k) 
      ENDDO 
     ENDDO
    ENDDO    
    !
    RicciPlusNablaZNCP = RiccitildeNCP + RicciphiNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) + ( nablaXNCP + TRANSPOSE(nablaXNCP) ) 
    RicciPlusNablaZSrc = RiccitildeSrc + RicciphiSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) + ( nablaXSrc + TRANSPOSE(nablaXSrc) ) 
    !
    !RicciPlusNablaZNCP = RicciNCP + ( nablaZNCP + TRANSPOSE(nablaZNCP) ) 
    !RicciPlusNablaZSrc = RicciSrc + ( nablaZSrc + TRANSPOSE(nablaZSrc) ) 
    !
    RPlusTwoNablaZNCP = phi**2*SUM(g_contr*RicciPlusNablaZNCP) 
    RPlusTwoNablaZSrc = phi**2*SUM(g_contr*RicciPlusNablaZSrc) 
    !
    !Kupdown = SUM(Kex*Kup) 
    !temp = RPlusTwoNablaZNCP + RPlusTwoNablaZSrc - KupDown + traceK**2 
    !temp = SUM(phi**2*g_contr*Ricci) - KupDown + traceK**2 
    !
    nablaijalphaNCP = 0.0
    nablaijalphaSrc = 0.0
    DO j = 1, 3 
     DO i = 1, 3 
       nablaijalphaNCP(i,j) = alpha*0.5*( dAA(i,j)+dAA(j,i) ) 
       nablaijalphaSrc(i,j) = alpha*AA(j)*AA(i) 
       DO k = 1, 3 
         nablaijalphaSrc(i,j) = nablaijalphaSrc(i,j) - Christoffel(i,j,k)*alpha*AA(k)  
       ENDDO
     ENDDO
    ENDDO 
    nablanablaalphaNCP = phi**2*SUM( g_contr*nablaijalphaNCP ) 
    nablanablaalphaSrc = phi**2*SUM( g_contr*nablaijalphaSrc ) 
!    nablanablaalpha = 0.0
!    DO j = 1, 3 
!     DO i = 1, 3
!         nablanablaalpha = nablanablaalpha + phi**2*( g_contr(i,j)*AA(j)*AA(i)*alpha + alpha*g_contr(i,j)*dAA(i,j) - 2*g_contr(i,j)*SUM(g_contr(:,:)*TRANSPOSE(DD(:,j,:))*AA(i)*alpha)- g_contr(i,j)*PP(i)*AA(j)*alpha )  
!     ENDDO
!    ENDDO 
    !
    SecondOrderTermsNCP = -nablaijalphaNCP + alpha*RicciPlusNablaZNCP 
    SecondOrderTermsSrc = -nablaijalphaSrc + alpha*RicciPlusNablaZSrc 
    traceNCP = SUM( g_contr*SecondOrderTermsNCP ) 
    SecondOrderTermsNCP = SecondOrderTermsNCP - 1./3.*g_cov*traceNCP 
    traceSrc = SUM( g_contr*SecondOrderTermsSrc ) 
    SecondOrderTermsSrc = SecondOrderTermsSrc - 1./3.*g_cov*traceSrc 
    !
    ! Now assemble all this terrible stuff... 
    !
    dtgamma = - 2*alpha*Aex - itau*(det-1.0)*g_cov 
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3
          dtgamma(i,j) = dtgamma(i,j) + g_cov(k,i)*BB(j,k) + g_cov(k,j)*BB(i,k) - 2./3.*g_cov(i,j)*BB(k,k) + beta(k)*2*DD(k,i,j) 
      ENDDO
     ENDDO
    ENDDO 
    !
    ! Main variables of the CCZ4 system 
    dtK = phi**2*SecondOrderTermsNCP + beta(1)*dAex(1,:,:) + beta(2)*dAex(2,:,:) + beta(3)*dAex(3,:,:)      ! extrinsic curvature
    dtK = dtK + phi**2*SecondOrderTermsSrc + alpha*Aex*(traceK-c*2*Theta) - 2*alpha*MATMUL(Aex,Amix) - itau*g_cov*traceA 
    DO j = 1, 3 
     DO i = 1, 3 
      DO k = 1, 3 
         dtK(i,j) = dtK(i,j) + Aex(k,i)*BB(j,k) + Aex(k,j)*BB(i,k) - 2./3.*Aex(i,j)*BB(k,k) 
      ENDDO
     ENDDO
    ENDDO 
    !
    dtTraceK = - nablanablaalphaNCP - nablanablaalphaSrc + alpha*( RPlusTwoNablaZNCP + RPlusTwoNablaZSrc + traceK**2 - c*2*Theta*traceK ) - 3*alpha*k1*Theta + SUM(beta(:)*dtraceK(:)) 
    !
    traceB = BB(1,1) + BB(2,2) + BB(3,3) 
    dtphi   = beta(1)*PP(1) + beta(2)*PP(2) + beta(3)*PP(3) + 1./3.*alpha*traceK - 1./3.*traceB 
    dtalpha = -alpha*fa*(traceK-K0-c*2*Theta) + beta(1)*AA(1) + beta(2)*AA(2) + beta(3)*AA(3) 

    Aupdown = SUM(Aex*Aup) 
    ! *** original 
    dtTheta = 0.5*alpha*e**2*(RplusTwoNablaZNCP + RplusTwoNablaZSrc) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)    &            ! temporal Z 
             + 0.5*alpha*e**2*( - Aupdown + 2./3.*traceK**2 ) - c*alpha*Theta*traceK - SUM(Zup*alpha*AA) - alpha*k2*Theta  
    ! *** use turbo cleaning here, i.e. without multiplying with alpha *** 
    !dtTheta = 0.5*e**2*( RplusTwoNablaZNCP + RplusTwoNablaZSrc ) + beta(1)*dTheta(1) + beta(2)*dTheta(2) + beta(3)*dTheta(3)    &                ! temporal Z 
    !         + 0.5*e**2*( - Aupdown + 2./3.*traceK**2 ) - c*alpha*Theta*traceK - SUM(Zup*alpha*AA) - k2*Theta 
    !
    !!divAupNCP = 0.0
    !!divAupSrc = 0.0 
    !!DO i = 1, 3
    !!    DO j = 1, 3
    !!     DO l = 1, 3
    !!      DO k = 1, 3    
    !!        divAupNCP(i) = divAupNCP(i) + g_contr(i,l)*g_contr(j,k)*dAex(j,l,k) 
    !!        divAupSrc(i) = divAupSrc(i) + ( dgup(j,i,l)*g_contr(j,k) + g_contr(i,l)*dgup(j,j,k) )*Aex(l,k) 
    !!      ENDDO
    !!     ENDDO
    !!    ENDDO        
    !!ENDDO 
    !!! 
    !!DO i = 1, 3 
    !!    Mom(i) = SUM(Christoffel_tilde(:,:,i)*Aup(:,:)) - 3.0*SUM(Aup(i,:)*PP(:)) - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) + divAupNCP(i) + divAupSrc(i)  
    !!ENDDO 
    !
    dKex = 0.0
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
          dKex(k,i,j) = 1.0/phi**2*( dAex(k,i,j) - 2.0*Aex(i,j)*PP(k) + 1./3.*dtraceK(k)*g_cov(i,j) + 2./3.*traceK*DD(k,i,j) - 2./3.*traceK*g_cov(i,j)*PP(k) ) 
      ENDDO
     ENDDO
    ENDDO 
    !
    Mom(:) = 0.0
    DO ii = 1, 3
        DO jj = 1, 3
            DO ll = 1, 3
                Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*(dKex(ll,ii,jj) - dKex(ii,jj,ll))
                DO mm = 1, 3
                    Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*( - Christoffel(jj,ll,mm)*Kex(mm,ii) + Christoffel(jj,ii,mm)*Kex(mm,ll)) 
                ENDDO
            ENDDO     
        ENDDO
    ENDDO   
    !
    dtX = alpha*ds**2*Mom + beta(1)*dXX(1,:) + beta(2)*dXX(2,:) + beta(3)*dXX(3,:) + MATMUL(BB,XX) 
    !
    DO i = 1, 3
        dtGhat(i) = 2*alpha*( SUM(Christoffel_tilde(:,:,i)*Aup(:,:)) - 3.0*SUM(Aup(i,:)*PP(:)) - 2./3.*SUM(g_contr(i,:)*dtraceK(:)) )    &
        !dtGhat(i)  =  2*alpha*( -divAupNCP(i) - divAupSrc(i) + ds**2*Mom(i) )                           & 
                    + 2*alpha*SUM( g_contr(:,i)*( dTheta(:) - Theta*AA(:) - 2./3.*traceK*Z(:)  ) )   & 
                    - 2*SUM( Aup(i,:)*alpha*AA(:) ) + 2./3.*Gtilde(i)*traceB - 2*alpha*k1*SUM(g_contr(i,:)*Z(:)) - SUM(Gtilde(:)*BB(:,i))   &
                   + beta(1)*dGhat(1,i) + beta(2)*dGhat(2,i) + beta(3)*dGhat(3,i) 
        DO l = 1, 3
         DO k = 1, 3
             dtGhat(i) = dtGhat(i) + g_contr(k,l)*0.5*(dBB(k,l,i)+dBB(l,k,i)) + 1./3*g_contr(i,k)*0.5*(dBB(k,l,l)+dBB(l,k,l)) + & 
                                     2*k3*( 2./3.*g_contr(i,l)*Z(l)*BB(k,k) - g_contr(l,k)*Z(l)*BB(k,i) ) 
         ENDDO
        ENDDO         
    ENDDO 
    DO k = 1, 3 
        ov(k) = + 2*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )   ! here we can use the constraint that trace A tilde = 0.         
    ENDDO        
    dtGhat = dtGhat + MATMUL(g_contr,ov)                                                  ! add the ordering constraint "up" (raised by g_contr) 
    !
    dtbb = xi*dtGhat - eta*b  !  <= be careful, this damping term -eta*b may be dangerous for the gamma driver, since it may kill the waves that you want !  
    !
    ! Add the following terms if you want shift convection in the PDE for b^i 
    dtbb = dtbb + bs*( - beta(1)*gradQ(14:16,1) - beta(2)*gradQ(14:16,2) - beta(3)*gradQ(14:16,3)  + beta(1)*gradQ(21:23,1) + beta(2)*gradQ(21:23,2) + beta(3)*gradQ(21:23,3)  ) !     
    dtbb = sk*dtbb 
    !
    QG = ggg*MATMUL( g_contr, AA - Z ) 
    DO j = 1, 3 
     DO i = 1, 3 
        DO n = 1, 3
         DO l = 1, 3 
           QG(i) = QG(i) + ggg*g_contr(i,j)*(-2*g_contr(n,l)*DD(l,n,j)) 
         ENDDO
        ENDDO
     ENDDO
    ENDDO 
    !            
    dtbeta  = + fff*b  
    ! Add the following term if you want to have shift convection in the PDE for beta^i 
    dtbeta = dtbeta + bs*( beta(1)*BB(1,:) + beta(2)*BB(2,:) + beta(3)*BB(3,:) )      
    dtbeta = sk*dtbeta 
    !
    ! Auxiliary variables 
    dtA = -alpha*fa*( dtraceK(:) -dK0(:) - c*2*dTheta(:) ) + beta(1)*dAA(1,:) + beta(2)*dAA(2,:) + beta(3)*dAA(3,:) - alpha*AA*(fa+alpha*faa)*(traceK-K0-c*2*Theta) + MATMUL(BB, AA) 
    DO k = 1, 3 
        dtA(k) = dtA(k) - alpha*fa*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )   ! here we can use the constraint that trace A tilde = 0. 
    ENDDO    
    ! 
    ! #xordB2# 
    ! for the ordering constraints, we have to check whether they should be multiplied by alpha**2, or not... 
    dtB = 0.0 
    DO i = 1, 3 
     DO k = 1, 3     
      DO j = 1, 3 
         !dtB(k,i) = +mu*alpha*g_contr(i,j)*( dZNCP(k,j) + dZSrc(k,j) ) + mu*dgup(k,i,j)*Z(j)      
         dtB(k,i) = dtB(k,i) + mu*alpha*g_contr(i,j)*( (dPP(k,j)-dPP(j,k)) )  
         DO n = 1, 3 
          DO l = 1, 3 
            !dtB(k,i) = dtB(k,i) + mu*alpha**2*g_contr(i,j)*( -2*g_contr(n,l)*0.5*(dDD(k,l,n,j)-dDD(l,k,n,j)) ) 
            dtB(k,i) = dtB(k,i) - mu*alpha*g_contr(i,j)*g_contr(n,l)*( dDD(k,l,j,n)-dDD(l,k,j,n) )   
          ENDDO 
         ENDDO 
      ENDDO 
      !
      ENDDO 
    ENDDO 
    dtB = dtB + bs*( beta(1)*dBB(1,:,:) + beta(2)*dBB(2,:,:) + beta(3)*dBB(3,:,:) ) 
    dtB = dtB*sk 
    !
    dtD = -alpha*dAex  
    DO i = 1, 3
      DO j = 1, 3
       DO k = 1, 3 
        DO m = 1, 3
            dtD(k,i,j) = dtD(k,i,j) + ( 0.5*(g_cov(m,i)*0.5*(dBB(k,j,m)+dBB(j,k,m))+g_cov(m,j)*0.5*(dBB(k,i,m)+dBB(i,k,m)) ) - 1./3.*g_cov(i,j)*0.5*(dBB(k,m,m)+dBB(m,k,m)) ) 
            DO n = 1, 3
                dtD(k,i,j) = dtD(k,i,j) + 1./3*alpha*g_cov(i,j)*g_contr(n,m)*dAex(k,n,m) + 1./3.*alpha*g_cov(i,j)*dgup(k,n,m)*Aex(n,m)   ! explicitly remove the trace of tilde A again 
            ENDDO 
        ENDDO        
       ENDDO
      ENDDO
    ENDDO 
    ! 
    DO j = 1, 3 
     DO i = 1, 3
      DO k = 1, 3 
        dtD(k,i,j) = dtD(k,i,j) - alpha*AA(k)*Aex(i,j) 
        DO l = 1, 3 
          dtD(k,i,j) = dtD(k,i,j) + BB(k,l)*DD(l,i,j) + DD(k,l,i)*BB(j,l) + DD(k,l,j)*BB(i,l) - 2./3.*DD(k,i,j)*BB(l,l) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO         
    !
    dtD = dtD + beta(1)*dDD(1,:,:,:) + beta(2)*dDD(2,:,:,:) + beta(3)*dDD(3,:,:,:)
    !
    dtP = MATMUL(BB,PP) + beta(1)*dPP(1,:) + beta(2)*dPP(2,:) + beta(3)*dPP(3,:)    
    DO k = 1, 3 
     dtP(k) = dtP(k) + 1./3.*alpha*dtraceK(k) + 1./3.*alpha*AA(k)*traceK + 1./3.*alpha*( SUM(g_contr(:,:)*dAex(k,:,:)) + SUM(dgup(k,:,:)*Aex(:,:)) )  
     DO i = 1, 3 
          dtP(k) = dtP(k) - 1./3.*0.5*(dBB(k,i,i)+dBB(i,k,i)) 
     ENDDO
    ENDDO 
    !
    BgradQ(1:6)    = (/ dtgamma(1,1), dtgamma(1,2), dtgamma(1,3), dtgamma(2,2), dtgamma(2,3), dtgamma(3,3) /)          ! \tilde \gamma_ij 
    BgradQ(7:12)   = (/ dtK(1,1), dtK(1,2), dtK(1,3), dtK(2,2), dtK(2,3), dtK(3,3) /)                                  ! \tilde A_ij 
    BgradQ(13)     = dtTheta                                                                                           ! Theta       
    BgradQ(14:16)  = dtGhat(1:3)                                                                                       ! \hat \Gamma^i           
    BgradQ(17)     = dtalpha                                                                                           ! log alpha 
    BgradQ(18:20)  = dtbeta                                                                                            ! beta^i 
    BgradQ(21:23)  = dtbb                                                                                              ! b^i 
    BgradQ(24:26)  = dtA(1:3)                                                                                          ! A_k       
    BgradQ(27:35)  = (/ dtB(1,1), dtB(2,1), dtB(3,1), dtB(1,2), dtB(2,2), dtB(3,2), dtB(1,3), dtB(2,3), dtB(3,3) /)    ! B_k^i 
    BgradQ(36:41)  = (/ dtD(1,1,1), dtD(1,1,2), dtD(1,1,3), dtD(1,2,2), dtD(1,2,3), dtD(1,3,3) /)                      ! D_kij 
    BgradQ(42:47)  = (/ dtD(2,1,1), dtD(2,1,2), dtD(2,1,3), dtD(2,2,2), dtD(2,2,3), dtD(2,3,3) /)                      ! D_kij 
    BgradQ(48:53)  = (/ dtD(3,1,1), dtD(3,1,2), dtD(3,1,3), dtD(3,2,2), dtD(3,2,3), dtD(3,3,3) /)                      ! D_kij 
    BgradQ(54)     = dtTraceK                                                                                          ! traceK 
    BgradQ(55)     = dtphi                                                                                             ! log phi 
    BgradQ(56:58)  = dtP                                                                                               ! P_k 
    BgradQ(59:61)  = dtX                                                                                               ! X_k 
    !
    Src_BgradQ = BgradQ    ! here, we do not have to change sign, since we work on the right hand side in the fused subroutine 
    !
    RETURN
    !
#endif 
    ! 
    !
    ! Default setting. Call PDENCP and PDESource separately. 
    !
    CALL PDENCP(BgradQ,Q,gradQ,par) 
    CALL PDESource(src,Q,par,time)  
    src_bgradQ = src - BgradQ 
    ! 
    CONTINUE     
    !            
END SUBROUTINE PDEFusedSrcNCP 
    
REAL FUNCTION frac(a,b)
    REAL :: a, b
    frac = a/b 
END FUNCTION frac 