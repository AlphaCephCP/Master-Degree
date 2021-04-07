!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDEFluxPrim(F,V,Q,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: V(nVar),  Q(nVar), par(nParam)  
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
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
#ifdef GRMHD
    !
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    rho    = V(1)
    vf_cov = V(2:4)
    p      = V(5)
    !
    BV_contr(1:3) = V(6:8)  ! B is contravariant
    QB_contr(1:3) = Q(6:8)  ! B is contravariant
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
    ! 
    !CALL METRIC(x, lapse, gp, gm, shift, g_cov, g_contr)
    !
    vf     = MATMUL(g_contr,vf_cov)
    Qv_contr = MATMUL(g_contr,Q(2:4))
    BQ = MATMUL(g_cov,QB_contr)
    BV = MATMUL(g_cov,BV_contr)
    vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
    vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
    vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
    !vxB_contr = MATMUL(g_contr,vxB(1:3))
    vxB_contr(1) = vf(2)*BV_contr(3) - vf(3)*BV_contr(2)
    vxB_contr(2) = vf(3)*BV_contr(1) - vf(1)*BV_contr(3)
    vxB_contr(3) = vf(1)*BV_contr(2) - vf(2)*BV_contr(1) 
    !
    !v2     = vx**2 + vy**2 + vz**2
    !b2     = bx**2 + by**2 + bz**2
    !e2     = ex**2 + ey**2 + ez**2 
    v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
    e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
    b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
    !
    uem    = 0.5*(b2 + e2) 
    !
    lf     = 1.0/sqrt(1.0 - v2)
    w      = rho + gamma1*p   ! rho*hentalpy
    ww     = w*lf**2
    wwx    = ww*vf(1)
    wwy    = ww*vf(2)
    wwz    = ww*vf(3) 
    !
    ! transport velocity
    Vtr(1:3) = lapse*vf(1:3)-shift(1:3)
    !
    !    Fij(1,1:3) =  (/ f1, g1, h1 /) for Q(6)  
    !    Fij(2,1:3) =  (/ f2, g2, h2 /) for Q(7)  ... without divergence cleaning
    !    Fij(3,1:3) =  (/ f3, g3, h3 /) for Q(8)  
    !
    DO m=1,3
        DO i=1,3
            Fij(i,m) = -Vtr(i)*QB_contr(m)+Vtr(m)*QB_contr(i)  ! Fij(i,i) = 0 !!!! ! NB: this is contravariant !!!!
        ENDDO
    ENDDO
    !
    F(1,1)   = vf(1)*Q(1) !rho*lf   !Q(1) ! rho*lf 
    F(2,1)   = wwx*vf_cov(1) - vxB_contr(1)*vxB(1) - BV_contr(1)*BV(1) + p + uem
    F(3,1)   = wwx*vf_cov(2) - vxB_contr(1)*vxB(2) - BV_contr(1)*BV(2) 
    F(4,1)   = wwx*vf_cov(3) - vxB_contr(1)*vxB(3) - BV_contr(1)*BV(3) 
    F(5,1)   = Qv_contr(1)-F(1,1) !wwx - f(1)         ! ADD MAGNETIC COMPONENT
    ! ADD MAGNETIC FIELD and DIV. CLEANING
    F(6,1)   = Fij(1,1) + V(9) !V(9)
    F(7,1)   = Fij(2,1) !-ez
    F(8,1)   = Fij(3,1) !ey  
    F(9,1)   = EQN%DivCleaning_a**2*QB_contr(1)
    !  lapse&shift&metric fluxes 
    F(10:19,1) = 0.
    !
    !
    F(1,2)   = vf(2)*Q(1) !rho*lf   ! rho*lf
    F(2,2)   = wwy*vf_cov(1) - vxB_contr(2)*vxB(1) - BV_contr(2)*BV(1) 
    F(3,2)   = wwy*vf_cov(2) - vxB_contr(2)*vxB(2) - BV_contr(2)*BV(2) + p + uem
    F(4,2)   = wwy*vf_cov(3) - vxB_contr(2)*vxB(3) - BV_contr(2)*BV(3) 
    F(5,2)   = Qv_contr(2)-F(1,2) !wwy - g(1)   ! ADD MAGNETIC COMPONENT
    ! ADD MAGNETIC FIELD and DIV. CLEANING
    F(6,2)   = Fij(1,2) !ez 
    F(7,2)   = Fij(2,2) + V(9) !V(9) 
    F(8,2)   = Fij(3,2) !-ex   
    F(9,2)   = EQN%DivCleaning_a**2*QB_contr(2)
    !  lapse&shift&metric fluxes 
    F(10:19,2) = 0.
    !
    !
    !
    F(1,3)   = vf(3)*Q(1) !rho*lf   ! rho*lf   !
    F(2,3)   = wwz*vf_cov(1) - vxB_contr(3)*vxB(1) - BV_contr(3)*BV(1) 
    F(3,3)   = wwz*vf_cov(2) - vxB_contr(3)*vxB(2) - BV_contr(3)*BV(2)   
    F(4,3)   = wwz*vf_cov(3) - vxB_contr(3)*vxB(3) - BV_contr(3)*BV(3) + p + uem
    F(5,3)   = Qv_contr(3)-F(1,3) ! wwz - h(1)   !ADD MAGNETIC COMPONENT
    ! ADD MAGNETIC FIELD and DIV. CLEANING
    F(6,3)   = Fij(1,3)  !-ey  
    F(7,3)   = Fij(2,3) !ex   
    F(8,3)   = Fij(3,3) + V(9) !V(9)   
    F(9,3)   = EQN%DivCleaning_a**2*QB_contr(3)
    !  lapse&shift&metric fluxes 
    F(10:19,3) = 0.
    ! 
    ! - - - - - - - - - REWRITE THE FOLLOWING 
    F(2:4,1)   = F(2:4,1)*gp
    F(2:4,2)   = F(2:4,2)*gp
    F(2:4,3)   = F(2:4,3)*gp
    ! Remember that Q(:) below contains already the factor gp, which is ok!
    F(1:5,1)   = lapse*F(1:5,1) - shift(1)*Q(1:5)
    F(1:5,2)   = lapse*F(1:5,2) - shift(2)*Q(1:5)
    F(1:5,3)   = lapse*F(1:5,3) - shift(3)*Q(1:5)
    !
    CONTINUE      
  !       
  !  
#endif 
  !
END SUBROUTINE PDEFluxPrim 
!
!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDESourcePrim(S,Q,par,time) 
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
END SUBROUTINE PDESourcePrim 
!
! Nonconservative part of the PDE ( B(Q) * gradQ ) 
!    
SUBROUTINE PDENCPPrim(BgradQ,Vc,gradQ,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Vc(nVar), gradQ(nVar,d), par(nParam)  
    REAL, INTENT(OUT) :: BgradQ(nVar) 
    ! Local variables 
    INTEGER :: i,j,k,l,m,n,ip,iq,ii,jj,kk,ll,iErr,count    
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar), Q(nVar)  
    REAL :: k1,k2,k3,fff,ggg,e,c,ds,xi,sk,sknl,bs,g_cov(3,3),g_contr(3,3),dgup(3,3,3)
    REAL :: det, alpha, fa, k0, dk0(3), beta0(3), b0(3), u(3), ialpha 
    REAL :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, eta, itau    
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar)  
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL :: v2,vf(3),uem,b2,e2,gp,gm,lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim 
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
#ifdef GRMHD     
    !
    CALL PDEPrim2Cons(Q,Vc) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)
    !
    lapse = Q(10)
    !gS(1:3) = Q(2:4) 
    shift = Q(11:13)
    !
    gammaij = Q(14:19) 
    g_cov(1,1) = Q(14)
    g_cov(1,2) = Q(15)
    g_cov(1,3) = Q(16)
    g_cov(2,2) = Q(17)
    g_cov(2,3) = Q(18)
    g_cov(3,3) = Q(19)
    g_cov(2,1) = Q(15)
    g_cov(3,1) = Q(16)
    g_cov(3,2) = Q(18)
    !
    delta = 0.
    DO i=1,3
        delta(i,i) = 1.0 
    ENDDO 
    !
    CALL MatrixInverse3x3(g_cov,g_contr,gp)
    gp = SQRT(gp)
    gm = 1./gp
    ! 
    gamma1 = EQN%gamma/(EQN%gamma-1.0)
    rho    = Vc(1)
    vf_cov = Vc(2:4)
    p      = Vc(5)
    !
    !BV(1:3) = Vc(6:8)
    BV_contr(1:3) = Vc(6:8) ! contravariant!
    QB_contr = Q(6:8)
    BV = MATMUL(g_cov,BV_contr)
    psi = Vc(9) 
    ! 
    vf       = MATMUL(g_contr,vf_cov)
    Qv_contr = MATMUL(g_contr,Q(2:4))
    !
    vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
    vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
    vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
    !
    !vxB_contr = MATMUL(g_contr,vxB(1:3))  THIS IS WRONG! vxB is a pseudo vector.
    vxB_contr(1) = vf(2)*BV_contr(3) - vf(3)*BV_contr(2)
    vxB_contr(2) = vf(3)*BV_contr(1) - vf(1)*BV_contr(3)
    vxB_contr(3) = vf(1)*BV_contr(2) - vf(2)*BV_contr(1) 
    !
    !v2     = vx**2 + vy**2 + vz**2
    !b2     = bx**2 + by**2 + bz**2
    !e2     = ex**2 + ey**2 + ez**2 
    v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
    e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
    b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
    !
    uem    = 0.5*(b2 + e2) 
    !
    !
    S_contr = MATMUL(g_contr,Q(2:4))
    !gv_contr = MATMUL(g_contr,gv)
    lf     = 1.0/sqrt(1.0 - v2)
    w      = rho + gamma1*p   ! rho*hentalpy
    ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
    !
    !DO j=1,3
    AQx = 0.
    BQy = 0.
    CQz = 0.
    ! 
      count=0
      DO i=1,3 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                AQx(1+j) = AQx(1+j) - Q(1+i)*Qx(10+i)  ! Q(11:13)  shift(i) or shift_contr(i)
                AQx(5) = AQx(5) - gp*W_ij*Qx(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                BQy(1+j) = BQy(1+j) - Q(1+i)*Qy(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                BQy(5) = BQy(5) - gp*W_ij*Qy(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
                !
                CQz(1+j) = CQz(1+j) - Q(1+i)*Qz(10+i)   ! Q(11:13)  shift(i) or shift_contr(i)
                CQz(5) = CQz(5) - gp*W_ij*Qz(10+i)     ! Q(11:13)  shift(i) or shift_contr(i)
                !
          DO m=1,3
            IF(m.GE.i) THEN  
                count=count+1
                !
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)-BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                AQx(1+j) = AQx(1+j) - 0.5*gp*lapse*Wim*Qx(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                AQx(5) = AQx(5) - 0.5*gp*Wim*shift(j)*Qx(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                BQy(1+j) = BQy(1+j) - 0.5*gp*lapse*Wim*Qy(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                BQy(5) = BQy(5) - 0.5*gp*Wim*shift(j)*Qy(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                CQz(1+j) = CQz(1+j) - 0.5*gp*lapse*Wim*Qz(13+count)  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                CQz(5) = CQz(5) - 0.5*gp*Wim*shift(j)*Qz(13+count)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
            ENDIF
          ENDDO
      ENDDO
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    AQx(1+j) = AQx(1+j) + (Q(5)+Q(1))*Qx(10)    ! Q(10) or lapse
    AQx(5) = AQx(5) + S_contr(j)*Qx(10)         !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    BQy(1+j) = BQy(1+j) + (Q(5)+Q(1))*Qy(10)    ! Q(10) or lapse
    BQy(5) = BQy(5) + S_contr(j)*Qy(10)         !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    CQz(1+j) = CQz(1+j) + (Q(5)+Q(1))*Qz(10)    ! Q(10) or lapse
    CQz(5) = CQz(5) + S_contr(j)*Qz(10)         !  Q(10) or lapse
    !
    BgradQ = AQx + BQy + CQz 
    CONTINUE    
    !     
#endif    
    !
END SUBROUTINE PDENCPPrim      

SUBROUTINE PDEFusedSrcNCPPrim(Src_BgradQ,V,gradQ,par,time)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: V(nVar), gradQ(nVar,d), par(nParam), time 
    REAL, INTENT(OUT) :: Src_BgradQ(nVar) 
    ! Local variables 
    REAL :: BgradQ(nVar), src(nVar) 
    ! 
    CALL PDENCP(BgradQ,V,gradQ,par) 
    CALL PDESource(src,V,par,time)  
    src_bgradQ = src - BgradQ 
    !
    CONTINUE     
    !            
END SUBROUTINE PDEFusedSrcNCPPrim 
   
    
SUBROUTINE PDEMatrixBPrim(Bn,V,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: V(nVar), par(nParam), nv(d)   
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
END SUBROUTINE PDEMatrixBPrim 
    
    
SUBROUTINE PDEEigenvaluesPrim(Lambda,V,par,nv) 
    USE typesDef, ONLY : nVar, nParam, d, nDim, EQN 
    USE recipies_mod, ONLY : RG 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: V(nVar), par(nParam), nv(d) 
    REAL, INTENT(OUT) :: Lambda(nVar) 
    ! Local variables 
    INTEGER :: i, j, iErr, itemp(nVar) 
    REAL    :: p, u, c, cp, cs, rho0, irho, lam, mu, alpha   
    REAL    :: dfdQ(nVar,nVar), ImLambda(nVar), rtemp(nVar), vp(nVar) 
    REAL    :: AA(nVar,nVar), BB(nVar,nVar), nvv(d), fpp(nVar,d), fmm(nVar,d)    
    REAL    :: Qp(nVar), Qm(nVar), eps, t1, R(nVar,nVar), iR(nVar,nVar)
    REAL    :: AMaple(nVar,nVar), Adiff(nVar,nVar), Q(nVar)  
    REAL :: ComputeDet, psi, BV(3),BV_contr(3)
    REAL :: b2_4, lf, lf2M, VdotB 
    REAL :: g_contr(3,3), g_cov(3,3)
    REAl :: shift(3), vf(3), vf_cov(3)
    REAL :: b2,cs2,den,gg
    REAL :: lapse, gp, gm 
    REAL :: rho, gamma1,sft,v2,vn,w ! ,vx,vy,vz
    !
    Lambda = 0.0 
    !
    !
END SUBROUTINE PDEEigenvaluesPrim 

SUBROUTINE PDEMatrixCPrim(Cn,V,nv,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: V(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: Cn(nVar,nVar) 
    ! Local variables 
    INTEGER :: j 
    REAL    :: dfdQ(nVar,nVar), Qp(nVar), Qm(nVar), Fpp(nVar,d), Fmm(nVar,d), Bn(nVar,nVar) 
    REAL    :: eps 
    !  
    Cn = 0.0 
    !
END SUBROUTINE PDEMatrixCPrim 

SUBROUTINE PDEdVdtPrim(dVdt,dQdt,V)
    USE typesDef, ONLY : nVar, EQN 
    !    
    REAL :: dVdt(nVar), dQdt(nVar), dVdQ(nVar,nVar), dQdV(nVar,nVar), V(nVar), test(nVar,nVar)  
    REAL :: iLOF, kappa 
    !
    !dVdQ = 0.0
    !dQdV = 0.0 
    !DO i = 1, nVar
    !    dVdQ(i,i) = 1.0 
    !    dQdV(i,i) = 1.0 
    !ENDDO        
    !
    dVdt = dQdt 
    !
#ifdef GRMHD
    !
    kappa = EQN%gamma 
    iLOF = SQRT( 1.0 - V(2)**2 - V(3)**2 - V(4)**2 ) 
    ! 
    dQdV(1,1) = 1/iLOF      
    dQdV(1,2) = V(1)/iLOF**3*V(2)
    dQdV(1,3) = V(1)/iLOF**3*V(3)
    dQdV(1,4) = V(1)/iLOF**3*V(4)
    dQdV(1,5) = 0
    dQdV(1,6) = 0
    dQdV(1,7) = 0
    dQdV(1,8) = 0
    dQdV(1,9) = 0
    dQdV(2,1) = V(2)/iLOF**2
    dQdV(2,2) = -(-2*V(1)*kappa+2*V(1)*kappa*V(3)**2+2*V(1)*kappa*V(4)**2+V(1)*kappa*iLOF**2+2*V(1)-2*V(1)*V(3)**2-2*V(1)*V(4)**2-V(1)*iLOF**2-2*kappa*V(5)+2*kappa*V(5)*V(3)**2+2*kappa*V(5)*V(4)**2+kappa*V(5)*iLOF**2-V(8)**2*iLOF**4*kappa+V(8)**2*iLOF**4-V(7)**2*iLOF**4*kappa+V(7)**2*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(2,3) = (2*V(2)*V(3)*V(1)*kappa-2*V(3)*V(1)*V(2)+2*V(2)*V(3)*kappa*V(5)-V(6)*V(7)*iLOF**4*kappa+V(6)*V(7)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(2,4) = (2*V(2)*V(4)*V(1)*kappa-2*V(4)*V(1)*V(2)+2*V(2)*V(4)*kappa*V(5)-V(6)*V(8)*iLOF**4*kappa+V(6)*V(8)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(2,5) = kappa/(kappa-1)*V(2)/iLOF**2
    dQdV(2,6) = -V(4)*V(8)-V(3)*V(7)
    dQdV(2,7) = 2*V(2)*V(7)-V(3)*V(6)
    dQdV(2,8) = 2*V(2)*V(8)-V(4)*V(6)
    dQdV(2,9) = 0
    dQdV(3,1) = V(3)/iLOF**2
    dQdV(3,2) = (2*V(2)*V(3)*V(1)*kappa-2*V(3)*V(1)*V(2)+2*V(2)*V(3)*kappa*V(5)-V(6)*V(7)*iLOF**4*kappa+V(6)*V(7)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(3,3) = (2*V(1)*kappa*V(3)**2-2*V(1)*V(3)**2+2*kappa*V(5)*V(3)**2+V(1)*kappa*iLOF**2-V(1)*iLOF**2+kappa*V(5)*iLOF**2+V(6)**2*iLOF**4*kappa-V(6)**2*iLOF**4+V(8)**2*iLOF**4*kappa-V(8)**2*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(3,4) = (2*V(4)*V(3)*V(1)*kappa-2*V(4)*V(1)*V(3)+2*V(4)*V(3)*kappa*V(5)-V(7)*V(8)*iLOF**4*kappa+V(7)*V(8)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(3,5) = kappa/(kappa-1)*V(3)/iLOF**2
    dQdV(3,6) = 2*V(3)*V(6)-V(2)*V(7)
    dQdV(3,7) = -V(2)*V(6)-V(4)*V(8)
    dQdV(3,8) = 2*V(3)*V(8)-V(4)*V(7)
    dQdV(3,9) = 0
    dQdV(4,1) = V(4)/iLOF**2
    dQdV(4,2) = (2*V(2)*V(4)*V(1)*kappa-2*V(4)*V(1)*V(2)+2*V(2)*V(4)*kappa*V(5)-V(6)*V(8)*iLOF**4*kappa+V(6)*V(8)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(4,3) = (2*V(4)*V(3)*V(1)*kappa-2*V(4)*V(1)*V(3)+2*V(4)*V(3)*kappa*V(5)-V(7)*V(8)*iLOF**4*kappa+V(7)*V(8)*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(4,4) = (2*V(1)*kappa*V(4)**2-2*V(1)*V(4)**2+2*kappa*V(5)*V(4)**2+V(1)*kappa*iLOF**2-V(1)*iLOF**2+kappa*V(5)*iLOF**2+V(7)**2*iLOF**4*kappa-V(7)**2*iLOF**4+V(6)**2*iLOF**4*kappa-V(6)**2*iLOF**4)/(kappa-1)/iLOF**4
    dQdV(4,5) = kappa/(kappa-1)*V(4)/iLOF**2
    dQdV(4,6) = 2*V(4)*V(6)-V(2)*V(8)
    dQdV(4,7) = 2*V(4)*V(7)-V(3)*V(8)
    dQdV(4,8) = -V(3)*V(7)-V(2)*V(6)
    dQdV(4,9) = 0
    dQdV(5,1) = 1/iLOF**2
    dQdV(5,2) = (2*V(2)*V(1)*kappa-2*V(1)*V(2)+2*V(2)*kappa*V(5)-V(8)*iLOF**4*V(4)*V(6)*kappa+V(8)*iLOF**4*V(4)*V(6)+V(8)**2*iLOF**4*V(2)*kappa-V(8)**2*iLOF**4*V(2)+V(7)**2*iLOF**4*V(2)*kappa-V(7)**2*iLOF**4*V(2)-V(7)*iLOF**4*V(3)*V(6)*kappa+V(7)*iLOF**4*V(3)*V(6))/(kappa-1)/iLOF**4
    dQdV(5,3) = (2*V(3)*V(1)*kappa-2*V(1)*V(3)+2*V(3)*kappa*V(5)+V(8)**2*iLOF**4*V(3)*kappa-V(8)**2*iLOF**4*V(3)-V(8)*iLOF**4*V(4)*V(7)*kappa+V(8)*iLOF**4*V(4)*V(7)-V(6)*iLOF**4*V(2)*V(7)*kappa+V(6)*iLOF**4*V(2)*V(7)+V(6)**2*iLOF**4*V(3)*kappa-V(6)**2*iLOF**4*V(3))/(kappa-1)/iLOF**4
    dQdV(5,4) = (2*V(4)*V(1)*kappa-2*V(1)*V(4)+2*V(4)*kappa*V(5)-V(7)*iLOF**4*V(3)*V(8)*kappa+V(7)*iLOF**4*V(3)*V(8)+V(7)**2*iLOF**4*V(4)*kappa-V(7)**2*iLOF**4*V(4)+V(6)**2*iLOF**4*V(4)*kappa-V(6)**2*iLOF**4*V(4)-V(6)*iLOF**4*V(2)*V(8)*kappa+V(6)*iLOF**4*V(2)*V(8))/(kappa-1)/iLOF**4
    dQdV(5,5) = -(-kappa+kappa*iLOF**2-iLOF**2)/(kappa-1)/iLOF**2
    dQdV(5,6) = V(6)+V(4)**2*V(6)-V(4)*V(2)*V(8)-V(3)*V(2)*V(7)+V(3)**2*V(6)
    dQdV(5,7) = V(7)-V(4)*V(3)*V(8)+V(4)**2*V(7)+V(2)**2*V(7)-V(2)*V(3)*V(6)
    dQdV(5,8) = V(8)+V(3)**2*V(8)-V(3)*V(4)*V(7)-V(2)*V(4)*V(6)+V(2)**2*V(8)
    dQdV(5,9) = 0
    dQdV(6,1) = 0
    dQdV(6,2) = 0
    dQdV(6,3) = 0
    dQdV(6,4) = 0
    dQdV(6,5) = 0
    dQdV(6,6) = 1
    dQdV(6,7) = 0
    dQdV(6,8) = 0
    dQdV(6,9) = 0
    dQdV(7,1) = 0
    dQdV(7,2) = 0
    dQdV(7,3) = 0
    dQdV(7,4) = 0
    dQdV(7,5) = 0
    dQdV(7,6) = 0
    dQdV(7,7) = 1
    dQdV(7,8) = 0
    dQdV(7,9) = 0
    dQdV(8,1) = 0
    dQdV(8,2) = 0
    dQdV(8,3) = 0
    dQdV(8,4) = 0
    dQdV(8,5) = 0
    dQdV(8,6) = 0
    dQdV(8,7) = 0
    dQdV(8,8) = 1
    dQdV(8,9) = 0
    dQdV(9,1) = 0
    dQdV(9,2) = 0
    dQdV(9,3) = 0
    dQdV(9,4) = 0
    dQdV(9,5) = 0
    dQdV(9,6) = 0
    dQdV(9,7) = 0
    dQdV(9,8) = 0
    dQdV(9,9) = 1    
    !
    t1 = dqdv(3,3)*dqdv(4,4)      
    t2 = dqdv(2,5)*dqdv(5,2)
    t4 = dqdv(5,2)*dqdv(2,3)
    t5 = dqdv(3,4)*dqdv(4,5)
    t7 = dqdv(5,2)*dqdv(2,4)
    t8 = dqdv(4,3)*dqdv(3,5)
    t10 = dqdv(2,2)*dqdv(4,4)
    t11 = dqdv(5,3)*dqdv(3,5)
    t13 = dqdv(2,2)*dqdv(4,3)
    t14 = dqdv(3,4)*dqdv(5,5)
    t16 = dqdv(2,2)*dqdv(3,3)
    t17 = dqdv(4,4)*dqdv(5,5)
    t19 = dqdv(5,4)*dqdv(4,5)
    t21 = dqdv(3,3)*dqdv(5,2)
    t22 = dqdv(2,4)*dqdv(4,5)
    t24 = dqdv(4,4)*dqdv(3,2)
    t25 = dqdv(2,3)*dqdv(5,5)
    t27 = dqdv(2,2)*dqdv(5,3)
    t29 = dqdv(2,2)*dqdv(5,4)
    t31 = dqdv(4,3)*dqdv(3,2)
    t32 = dqdv(2,4)*dqdv(5,5)
    t34 = -t2*t1-t5*t4-t8*t7-t11*t10-t14*t13+t17*t16-t19*t16+t22*t21-t25*t24+t5*t27+t8*t29+t32*t31
    t35 = dqdv(4,4)*dqdv(5,3)
    t36 = dqdv(2,5)*dqdv(3,2)
    t38 = dqdv(4,3)*dqdv(3,4)
    t40 = dqdv(5,4)*dqdv(3,2)
    t41 = dqdv(2,3)*dqdv(4,5)
    t43 = dqdv(4,4)*dqdv(5,2)
    t44 = dqdv(2,3)*dqdv(3,5)
    t46 = dqdv(5,4)*dqdv(4,2)
    t48 = dqdv(4,2)*dqdv(2,3)
    t50 = dqdv(5,3)*dqdv(3,2)
    t52 = dqdv(5,3)*dqdv(3,4)
    t53 = dqdv(2,5)*dqdv(4,2)
    t55 = dqdv(5,4)*dqdv(4,3)
    t57 = dqdv(3,3)*dqdv(4,2)
    t59 = dqdv(3,3)*dqdv(5,4)
    t61 = dqdv(4,2)*dqdv(2,4)
    t63 = t36*t35+t2*t38+t41*t40+t44*t43-t44*t46+t14*t48-t22*t50-t53*t52-t36*t55-t32*t57+t53*t59+t11*t61
    t65 = dqdv(1,1)*dqdv(2,2)
    t68 = dqdv(2,1)*dqdv(4,5)
    t69 = dqdv(1,2)*t68
    t71 = dqdv(1,4)*dqdv(2,1)
    t72 = dqdv(5,5)*t71
    t74 = dqdv(4,1)*dqdv(1,4)
    t75 = dqdv(2,5)*dqdv(5,3)
    t76 = dqdv(3,2)*t75
    t78 = dqdv(2,3)*dqdv(3,2)
    t79 = dqdv(5,5)*t78
    t81 = dqdv(3,5)*t4
    t87 = dqdv(4,4)*dqdv(3,1)
    t88 = dqdv(2,3)*dqdv(1,2)
    t91 = dqdv(1,4)*dqdv(3,1)
    t92 = dqdv(4,5)*t91
    t94 = dqdv(1,1)*dqdv(5,4)
    t95 = dqdv(3,5)*t48
    t97 = dqdv(1,1)*dqdv(4,4)
    t99 = -dqdv(3,5)*t35*t65+t69*t59+t72*t57-t76*t74+t79*t74-t81*t74-dqdv(5,5)*t38*t65+dqdv(4,5)*t52*t65+dqdv(5,5)*t88*t87+t92*t4-t95*t94+t81*t97
    t100 = dqdv(5,4)*dqdv(3,1)
    t101 = dqdv(2,5)*dqdv(1,3)
    t104 = dqdv(4,2)*dqdv(2,1)
    t105 = dqdv(5,3)*dqdv(1,4)
    t106 = dqdv(3,5)*t105
    t108 = dqdv(5,5)*t91
    t110 = dqdv(2,5)*dqdv(4,3)
    t111 = dqdv(3,2)*t110
    t113 = dqdv(4,5)*t71
    t115 = dqdv(1,3)*dqdv(3,1)
    t116 = dqdv(4,5)*t115
    t118 = dqdv(2,1)*dqdv(5,5)
    t119 = dqdv(1,2)*t118
    t121 = dqdv(5,3)*dqdv(3,1)
    t122 = dqdv(2,4)*dqdv(1,2)
    t123 = dqdv(4,5)*t122
    t125 = dqdv(4,1)*dqdv(2,5)
    t126 = dqdv(1,2)*t125
    t128 = dqdv(5,2)*dqdv(2,1)
    t129 = dqdv(4,3)*dqdv(1,4)
    t130 = dqdv(3,5)*t129
    t132 = dqdv(5,1)*dqdv(1,3)
    t133 = dqdv(2,5)*dqdv(3,4)
    t134 = dqdv(4,2)*t133
    t136 = dqdv(3,5)*t132
    t138 = -dqdv(4,2)*t101*t100-t106*t104+t108*t13-t111*t94-t113*t21-t116*t7-t119*t1+t123*t121-t126*t59+t130*t128+t134*t132-t136*t61
    t140 = dqdv(1,3)*dqdv(4,1)
    t141 = dqdv(3,5)*t140
    t144 = dqdv(4,3)*dqdv(3,1)
    t145 = dqdv(5,5)*t122
    t147 = dqdv(2,4)*dqdv(3,2)
    t148 = dqdv(5,5)*t147
    t150 = dqdv(5,4)*dqdv(4,1)
    t151 = dqdv(3,5)*t88
    t158 = dqdv(3,1)*dqdv(2,5)
    t159 = dqdv(1,2)*t158
    t161 = dqdv(5,1)*dqdv(1,4)
    t163 = dqdv(5,1)*dqdv(2,5)
    t164 = dqdv(1,2)*t163
    t166 = t141*t7-t108*t48-t145*t144-t148*t140+t151*t150-dqdv(5,5)*t74*t16+t126*t52+t136*t10+t119*t38-t159*t35+t95*t161-t164*t38
    t170 = dqdv(1,1)*dqdv(3,3)
    t173 = dqdv(2,5)*dqdv(1,4)
    t174 = dqdv(5,2)*t173
    t176 = dqdv(4,5)*t147
    t178 = dqdv(1,3)*dqdv(2,1)
    t179 = dqdv(3,5)*t178
    t181 = dqdv(2,2)*dqdv(4,1)
    t183 = dqdv(5,1)*dqdv(1,2)
    t184 = dqdv(4,3)*dqdv(2,4)
    t185 = dqdv(3,5)*t184
    t187 = dqdv(4,2)*t173
    t189 = dqdv(4,4)*dqdv(5,1)
    t191 = dqdv(3,4)*dqdv(1,3)
    t192 = dqdv(5,5)*t191
    t195 = t111*t161+dqdv(3,5)*t55*t65+dqdv(4,5)*t7*t170-t174*t144+t176*t132-t179*t43+t106*t181+t185*t183+t187*t121-t151*t189-t192*t104-t141*t29
    t200 = dqdv(3,2)*t101
    t207 = dqdv(1,1)*dqdv(4,2)
    t208 = dqdv(3,4)*dqdv(2,3)
    t209 = dqdv(5,5)*t208
    t211 = dqdv(1,1)*dqdv(5,3)
    t213 = dqdv(4,5)*t78
    t219 = -dqdv(4,5)*t88*t100+t200*t150-dqdv(4,5)*t59*t65+dqdv(5,5)*t1*t65+t159*t55+t209*t207-t176*t211+t213*t94+t113*t50-t134*t211+t76*t97-t69*t52
    t220 = dqdv(1,1)*dqdv(4,3)
    t221 = dqdv(5,2)*t133
    t225 = dqdv(2,1)*dqdv(3,5)
    t226 = dqdv(1,2)*t225
    t230 = dqdv(1,1)*dqdv(5,2)
    t231 = dqdv(4,5)*t208
    t233 = dqdv(5,5)*t115
    t235 = dqdv(2,2)*dqdv(5,1)
    t238 = dqdv(2,5)*dqdv(4,4)
    t243 = dqdv(4,1)*dqdv(1,2)
    t244 = dqdv(5,3)*dqdv(2,4)
    t245 = dqdv(3,5)*t244
    t247 = t221*t220-dqdv(4,5)*t178*t40-t226*t55+dqdv(4,5)*t161*t16-t231*t230-t233*t10-t130*t235+t233*t61-dqdv(5,2)*t238*t170-t200*t189+t226*t35-t245*t243
    t251 = dqdv(3,3)*dqdv(5,1)
    t255 = dqdv(4,5)*t191
    t261 = dqdv(3,3)*dqdv(4,1)
    t264 = -t72*t31-t209*t243-t187*t251-t123*t251-t221*t140-t255*t235+t192*t181+t116*t29-t213*t161+t255*t128+t174*t261+t179*t46
    t277 = dqdv(2,5)*dqdv(5,4)
    t282 = -t185*t230+t148*t220+t164*t1+t245*t207-t79*t97+dqdv(5,5)*t178*t24-t92*t27+dqdv(5,2)*t101*t87-dqdv(5,5)*t61*t170+dqdv(4,2)*t277*t170+t231*t183+t145*t261
    t286 = 1/(t99+t138+t166+t195+t219+t247+t264+t282)
    t288 = dqdv(3,3)*dqdv(5,5)
    t289 = dqdv(1,2)*dqdv(4,4)
    t291 = dqdv(4,5)*dqdv(1,2)
    t293 = dqdv(3,3)*dqdv(1,4)
    t294 = dqdv(4,5)*dqdv(5,2)
    t296 = dqdv(1,4)*dqdv(4,2)
    t298 = dqdv(1,3)*dqdv(3,2)
    t300 = dqdv(4,4)*dqdv(1,3)
    t301 = dqdv(3,5)*dqdv(5,2)
    t303 = dqdv(3,5)*dqdv(1,2)
    t306 = dqdv(5,5)*dqdv(1,2)
    t311 = dqdv(1,3)*dqdv(3,5)
    t314 = dqdv(5,5)*dqdv(1,3)
    t315 = dqdv(3,4)*dqdv(4,2)
    t317 = dqdv(1,4)*dqdv(4,5)
    t319 = dqdv(5,5)*dqdv(1,4)
    t322 = t289*t288-t291*t59+t294*t293-t296*t288-t298*t17+t301*t300-t303*t35+t296*t11-t306*t38+t291*t52-t301*t129+t298*t19-t46*t311+t303*t55+t315*t314-t50*t317+t31*t319-t294*t191
    t324 = dqdv(2,2)*dqdv(5,5)
    t327 = dqdv(4,5)*dqdv(1,3)
    t329 = dqdv(2,2)*dqdv(1,4)
    t330 = dqdv(4,5)*dqdv(5,3)
    t334 = dqdv(5,2)*dqdv(1,3)
    t336 = dqdv(4,2)*dqdv(5,3)
    t345 = dqdv(1,2)*dqdv(2,5)
    t348 = dqdv(4,2)*dqdv(1,3)
    t350 = -t300*t324+t129*t324+t327*t29-t330*t329-t75*t289+t88*t17+t334*t238+t336*t173-t184*t306-t48*t319+t330*t122+t314*t61+t4*t317-t88*t19-t327*t7+t55*t345-t129*t2-t348*t277
    t353 = dqdv(3,5)*dqdv(1,4)
    t355 = dqdv(2,2)*dqdv(1,3)
    t356 = dqdv(3,5)*dqdv(5,4)
    t360 = dqdv(3,3)*dqdv(2,5)
    t361 = dqdv(5,2)*dqdv(1,4)
    t363 = dqdv(3,3)*dqdv(1,2)
    t368 = dqdv(3,2)*dqdv(1,4)
    t377 = t293*t324-t353*t27+t356*t355-t191*t324-t122*t288-t361*t360+t277*t363-t40*t101+t191*t2+t208*t306+t368*t75-t52*t345+t122*t11+t147*t314-t7*t311-t319*t78+t353*t4-t356*t88
    t379 = dqdv(4,5)*dqdv(3,3)
    t398 = t379*t329+t311*t10-t5*t355-t8*t329+t345*t1-t22*t363-t53*t293-t311*t61+t5*t88+t44*t296+t8*t122-t44*t289-t36*t300-t345*t38+t22*t298+t36*t129+t53*t191-t41*t368
    t400 = dqdv(3,6)*t327
    t402 = dqdv(3,6)*t88
    t404 = dqdv(4,6)*t356
    t406 = dqdv(2,6)*dqdv(3,4)
    t407 = dqdv(4,2)*t406
    t409 = dqdv(5,6)*t22
    t413 = dqdv(4,6)*t122
    t415 = dqdv(2,6)*dqdv(4,5)
    t416 = dqdv(5,2)*t415
    t418 = dqdv(5,6)*t36
    t420 = dqdv(2,6)*dqdv(1,3)
    t421 = dqdv(3,2)*t420
    t423 = dqdv(3,6)*t314
    t425 = dqdv(2,6)*t291
    t427 = dqdv(2,6)*t303
    t429 = dqdv(4,6)*t40
    t431 = dqdv(4,6)*t353
    t434 = dqdv(3,6)*t48
    t436 = dqdv(3,6)*t129
    t438 = -t400*t7+t402*t17+t404*t88-t407*t314-t409*t363+dqdv(4,6)*t361*t360+t413*t288-t416*t293+t418*t129+t421*t17-t423*t10-t425*t52-t427*t55+t429*t101-t431*t4+t409*t298-t434*t319-t436*t2
    t439 = dqdv(5,6)*t8
    t442 = dqdv(5,6)*t5
    t445 = dqdv(4,6)*t191
    t447 = dqdv(3,6)*t330
    t449 = dqdv(5,6)*t61
    t452 = dqdv(2,6)*dqdv(3,5)
    t453 = dqdv(5,2)*t452
    t455 = dqdv(2,6)*dqdv(1,4)
    t456 = dqdv(4,2)*t455
    t458 = dqdv(4,6)*t7
    t461 = dqdv(3,6)*t336
    t466 = dqdv(3,6)*t75
    t468 = dqdv(4,6)*t147
    t470 = dqdv(5,6)*t317
    t472 = -t439*t329+t431*t27-t442*t355+t400*t29+t445*t324-t447*t329-t449*t311-t402*t19+t453*t129+t456*t288+t458*t311-t445*t2+t461*t173+t442*t88+dqdv(5,6)*t345*t1-t466*t289-t468*t314+t470*t16
    t475 = dqdv(3,6)*t4
    t477 = dqdv(2,6)*t306
    t481 = dqdv(5,6)*t44
    t483 = dqdv(4,6)*t277
    t488 = dqdv(5,6)*t38
    t493 = dqdv(5,6)*t78
    t496 = dqdv(4,6)*t208
    t498 = dqdv(3,6)*t184
    t500 = dqdv(5,6)*t53
    t503 = t447*t122+t475*t317-t477*t1-dqdv(3,6)*t348*t277+t481*t296-t483*t363+t425*t59-dqdv(4,6)*t368*t75-t488*t345-t418*t300+t436*t324+t427*t35-t493*t317+t423*t61-t496*t306-t498*t306+t500*t191+t439*t122
    t505 = dqdv(4,3)*dqdv(2,6)*dqdv(3,2)
    t508 = dqdv(5,3)*dqdv(2,6)*dqdv(3,2)
    t510 = dqdv(4,6)*t319
    t517 = dqdv(5,4)*dqdv(2,6)*dqdv(4,2)
    t520 = dqdv(3,6)*t55
    t525 = dqdv(4,6)*t52
    t529 = dqdv(5,6)*t311
    t533 = -t505*t319+t508*t317-t510*t16+t477*t38-t413*t11-t456*t11+t416*t191+t517*t311-t421*t19+t520*t345+dqdv(3,6)*t334*t238-t500*t293+t525*t345-t453*t300-t404*t355+t529*t10+t510*t78-t481*t289
    t538 = dqdv(5,3)*dqdv(2,7)*dqdv(3,2)
    t542 = dqdv(3,7)*t327
    t544 = dqdv(5,7)*t8
    t546 = dqdv(4,7)*t353
    t548 = dqdv(3,7)*t184
    t551 = dqdv(4,3)*dqdv(2,7)*dqdv(3,2)
    t553 = dqdv(5,7)*t5
    t555 = dqdv(3,7)*t330
    t557 = dqdv(3,7)*t55
    t559 = dqdv(4,7)*t356
    t561 = dqdv(4,7)*t52
    t564 = dqdv(5,4)*dqdv(2,7)*dqdv(4,2)
    t566 = dqdv(2,7)*dqdv(4,5)
    t567 = dqdv(5,2)*t566
    t569 = dqdv(3,7)*t48
    t571 = dqdv(5,7)*t44
    t575 = dqdv(4,7)*t319
    t577 = t538*t317+dqdv(3,7)*t334*t238+t542*t29-t544*t329+t546*t27-t548*t306-t551*t319-t553*t355-t555*t329+t557*t345-t559*t355+t561*t345+t564*t311+t567*t191-t569*t319-t571*t289-dqdv(3,7)*t348*t277-t575*t16
    t578 = dqdv(3,7)*t314
    t580 = dqdv(5,7)*t311
    t582 = dqdv(2,7)*dqdv(1,3)
    t583 = dqdv(3,2)*t582
    t585 = dqdv(3,7)*t129
    t588 = dqdv(4,7)*t40
    t593 = dqdv(5,7)*t36
    t595 = dqdv(2,7)*dqdv(1,4)
    t596 = dqdv(4,2)*t595
    t601 = dqdv(3,7)*t4
    t603 = dqdv(4,7)*t208
    t605 = dqdv(4,7)*t122
    t607 = dqdv(2,7)*t303
    t610 = -t578*t10+t580*t10-t583*t19+t585*t324+t544*t122+t588*t101-t546*t4+dqdv(4,7)*t361*t360+t593*t129+t596*t288-t593*t300+t571*t296+t583*t17+t601*t317-t603*t306+t605*t288+t607*t35+t578*t61
    t613 = dqdv(2,7)*t306
    t615 = dqdv(4,7)*t191
    t617 = dqdv(4,7)*t277
    t619 = dqdv(5,7)*t22
    t622 = dqdv(2,7)*dqdv(3,5)
    t623 = dqdv(5,2)*t622
    t626 = dqdv(5,7)*t78
    t629 = dqdv(2,7)*t291
    t631 = dqdv(3,7)*t88
    t633 = dqdv(3,7)*t336
    t638 = dqdv(5,7)*t53
    t641 = dqdv(4,7)*t7
    t643 = t559*t88+t613*t38-t615*t2-t617*t363-t619*t363-t596*t11+t623*t129+t553*t88-t626*t317-t607*t55+t629*t59-t631*t19+t633*t173-dqdv(4,7)*t368*t75+t619*t298+t638*t191-t605*t11+t641*t311
    t649 = dqdv(5,7)*t38
    t652 = dqdv(2,7)*dqdv(3,4)
    t653 = dqdv(4,2)*t652
    t655 = dqdv(5,7)*t61
    t660 = dqdv(4,7)*t147
    t664 = dqdv(3,7)*t75
    t668 = dqdv(5,7)*t317
    t670 = -t613*t1+dqdv(5,7)*t345*t1-t638*t293-t629*t52-t649*t345+t555*t122-t653*t314-t655*t311-t542*t7+t575*t78-t585*t2-t660*t314-t623*t300+t631*t17-t664*t289+t615*t324-t567*t293+t668*t16
    t674 = dqdv(5,8)*t36
    t676 = dqdv(5,8)*t38
    t678 = dqdv(3,8)*t327
    t680 = dqdv(4,8)*t40
    t682 = dqdv(5,8)*t78
    t684 = dqdv(4,8)*t191
    t686 = dqdv(4,8)*t147
    t690 = dqdv(2,8)*dqdv(3,4)
    t691 = dqdv(4,2)*t690
    t693 = dqdv(4,8)*t208
    t695 = dqdv(2,8)*dqdv(1,3)
    t696 = dqdv(3,2)*t695
    t699 = dqdv(2,8)*t306
    t703 = dqdv(3,8)*t330
    t705 = dqdv(2,8)*dqdv(3,5)
    t706 = dqdv(5,2)*t705
    t709 = dqdv(5,3)*dqdv(2,8)*dqdv(3,2)
    t711 = dqdv(4,8)*t353
    t713 = -t674*t300-t676*t345-t678*t7+t680*t101-t682*t317-t684*t2-t686*t314-dqdv(4,8)*t368*t75-t691*t314-t693*t306-t696*t19+t684*t324-t699*t1+dqdv(5,8)*t345*t1-t703*t329+t706*t129+t709*t317-t711*t4
    t714 = dqdv(2,8)*t303
    t719 = dqdv(3,8)*t88
    t721 = dqdv(3,8)*t75
    t723 = dqdv(4,8)*t122
    t725 = dqdv(3,8)*t55
    t727 = dqdv(2,8)*t291
    t729 = dqdv(3,8)*t314
    t732 = dqdv(4,3)*dqdv(2,8)*dqdv(3,2)
    t735 = dqdv(2,8)*dqdv(1,4)
    t736 = dqdv(4,2)*t735
    t739 = dqdv(5,8)*t317
    t742 = dqdv(2,8)*dqdv(4,5)
    t743 = dqdv(5,2)*t742
    t746 = dqdv(4,8)*t52
    t748 = -t714*t55+t699*t38-dqdv(3,8)*t348*t277-t719*t19-t721*t289-t723*t11+t725*t345-t727*t52+t729*t61-t732*t319-t706*t300-t736*t11+t723*t288+t739*t16+t696*t17-t743*t293+t714*t35+t746*t345
    t750 = dqdv(5,8)*t8
    t753 = dqdv(4,8)*t319
    t755 = dqdv(5,8)*t311
    t758 = dqdv(3,8)*t129
    t761 = dqdv(5,8)*t53
    t763 = dqdv(5,8)*t22
    t766 = dqdv(5,4)*dqdv(2,8)*dqdv(4,2)
    t768 = dqdv(3,8)*t48
    t771 = dqdv(4,8)*t277
    t775 = dqdv(5,8)*t44
    t777 = dqdv(3,8)*t4
    t779 = dqdv(5,8)*t61
    t781 = t750*t122+t703*t122-t753*t16+t755*t10-t729*t10-t758*t2+t753*t78+t761*t191-t763*t363+t766*t311-t768*t319-t761*t293-t771*t363+t719*t17+t763*t298-t775*t289+t777*t317-t779*t311
    t784 = dqdv(4,8)*t356
    t788 = dqdv(5,8)*t5
    t793 = dqdv(3,8)*t336
    t797 = dqdv(4,8)*t7
    t801 = dqdv(3,8)*t184
    t807 = dqdv(3,8)*t334*t238-t784*t355+t758*t324-t750*t329-t788*t355+t711*t27+t678*t29+t784*t88+t793*t173+t674*t129+t743*t191+t797*t311+t727*t59+t736*t288-t801*t306+dqdv(4,8)*t361*t360+t788*t88+t775*t296
    t814 = dqdv(4,5)*dqdv(5,1)
    t816 = dqdv(2,4)*dqdv(4,1)
    t822 = dqdv(2,3)*dqdv(3,1)
    t824 = dqdv(4,4)*dqdv(2,3)
    t825 = dqdv(3,5)*dqdv(5,1)
    t829 = -t163*t1+t118*t1+t814*dqdv(3,3)*dqdv(2,4)-t816*t288+t150*t360-t68*t59-t225*t35+t121*t238-t822*t17+t825*t824-t814*t208-t150*t44
    t833 = dqdv(3,4)*dqdv(4,1)
    t842 = dqdv(3,4)*dqdv(5,1)
    t844 = t822*t19-t118*t38+t144*t32+t833*t25+t816*t11-t121*t22+t225*t55+t68*t52-t833*t75-t825*t184-t144*t277+t842*t110
    t865 = t17*t170-t19*t170-t11*t97-t14*t220+t5*t211+t8*t94-t319*t261+t317*t251-t8*t161-t314*t87-t317*t121-t311*t150+t311*t189+t319*t144+t11*t74-t5*t132+t327*t100+t14*t140
    t869 = dqdv(1,1)*dqdv(5,5)
    t872 = dqdv(1,1)*dqdv(2,4)
    t874 = dqdv(1,1)*dqdv(2,5)
    t888 = -t75*t97+t25*t97-t184*t869-t41*t94+t330*t872+t55*t874-t178*t17+t132*t238+t140*t32-t129*t163-t330*t71+t178*t19-t140*t277+t105*t125+t41*t161-t132*t22-t25*t74+t129*t118
    t892 = dqdv(3,5)*dqdv(2,4)
    t895 = dqdv(1,1)*dqdv(2,3)
    t906 = dqdv(1,3)*dqdv(5,4)
    t911 = -t32*t170+t277*t170+t892*t211+t208*t869-t356*t895-t52*t874+t71*t288-t161*t360-t91*t25-t71*t11-t892*t132+t191*t163+t91*t75+t32*t115-t906*t158+t356*t178+t161*t44-t191*t118
    t920 = dqdv(3,3)*dqdv(2,1)
    t922 = dqdv(4,4)*dqdv(2,1)
    t930 = dqdv(2,4)*dqdv(3,1)
    t934 = t238*t170-t22*t170-t44*t97+t5*t895+t8*t872-t133*t220-t173*t261+t317*t920+t311*t922-t115*t238-t311*t816+t44*t74+t133*t140-t8*t71+t173*t144+t327*t930-t317*t822-t5*t178
    t936 = dqdv(3,6)*t178
    t940 = dqdv(2,6)*t317
    t942 = dqdv(4,6)*t71
    t948 = dqdv(5,6)*t115
    t952 = dqdv(2,6)*t311
    t956 = dqdv(5,6)*t91
    t958 = dqdv(4,6)*t892
    t962 = -t936*t17-dqdv(2,6)*t17*t170-t940*t251-t942*t288+t470*t920-t442*t178+t442*t895-t436*t163-t948*t238-t470*t822-t409*t170+t952*t150-t529*t816+t404*t895+t956*t110-t958*t211+t520*t874-t488*t874
    t964 = dqdv(3,6)*t132
    t967 = dqdv(3,6)*t41
    t977 = dqdv(3,6)*t25
    t982 = dqdv(2,6)*t11
    t984 = dqdv(3,6)*t140
    t986 = dqdv(2,6)*t319
    t988 = dqdv(4,6)*t32
    t990 = t447*t872-t964*t22-t466*t97+t967*t161-dqdv(2,6)*t14*t140+dqdv(5,6)*t238*t170-t498*t869+t525*t874-t439*t71-t496*t869-t977*t74-t967*t94-t952*t189+t529*t922-t982*t74+t984*t32+t986*t261-t988*t115
    t992 = dqdv(2,6)*t8
    t994 = dqdv(4,6)*t115
    t996 = dqdv(4,6)*t91
    t998 = dqdv(2,6)*t5
    t1011 = dqdv(5,6)*t74
    t1013 = dqdv(4,6)*t161
    t1017 = t992*t161+t994*t277+t996*t25+t998*t132+t964*t238+t940*t121-t996*t75-t992*t94-t483*t170+t988*t170-t404*t178+t942*t11-t445*t163+dqdv(2,6)*t19*t170-t1011*t360+t1013*t360+t439*t872-t984*t277
    t1029 = dqdv(5,6)*t191
    t1031 = dqdv(3,6)*t105
    t1041 = -t1013*t44+t1011*t44-t998*t211+t948*t22+t958*t132+dqdv(2,6)*t314*t87+t436*t118-dqdv(2,6)*t327*t100-t986*t144+t1029*t125+t1031*t125+t445*t118+t982*t97-t481*t97+dqdv(2,6)*t38*t869+t977*t97+t936*t19-t447*t71
    t1045 = dqdv(4,7)*t161
    t1047 = dqdv(5,7)*t74
    t1051 = dqdv(4,7)*t91
    t1053 = dqdv(3,7)*t105
    t1057 = dqdv(2,7)*t319
    t1059 = dqdv(3,7)*t25
    t1061 = dqdv(3,7)*t178
    t1063 = dqdv(4,7)*t32
    t1065 = dqdv(4,7)*t71
    t1070 = dqdv(2,7)*t11
    t1072 = dqdv(2,7)*t317
    t1075 = dqdv(3,7)*t41
    t1077 = -t1045*t44+t1047*t44-t617*t170-t619*t170+t1051*t25+t1053*t125+dqdv(5,7)*t238*t170-t1057*t144+t1059*t97+t1061*t19-t1063*t115+t1065*t11+t553*t895-dqdv(2,7)*t14*t140-t1070*t74+t1072*t121-t544*t71+t1075*t161
    t1078 = dqdv(4,7)*t892
    t1084 = dqdv(3,7)*t132
    t1086 = dqdv(2,7)*t8
    t1088 = dqdv(5,7)*t115
    t1099 = dqdv(2,7)*t311
    t1104 = t1078*t132-dqdv(2,7)*t17*t170-t1075*t94+t559*t895-t1084*t22-t1086*t94+t1088*t22+dqdv(2,7)*t314*t87-t548*t869-t1061*t17+t580*t922-t1088*t238-t571*t97+dqdv(2,7)*t19*t170+t1099*t150-t580*t816+t615*t118+t585*t118
    t1112 = dqdv(2,7)*t5
    t1116 = dqdv(3,7)*t140
    t1119 = dqdv(5,7)*t91
    t1124 = dqdv(4,7)*t115
    t1129 = t555*t872-t603*t869-t1051*t75+dqdv(2,7)*t38*t869-t1065*t288-t1112*t211-t649*t874-t668*t822-t1116*t277-t553*t178+t1119*t110-t1099*t189+t1112*t132+t1116*t32+t1124*t277+t561*t874+t557*t874-t1078*t211
    t1136 = dqdv(5,7)*t191
    t1150 = t1057*t261+t1045*t360-t1072*t251+t1084*t238+t668*t920-t1047*t360+t1136*t125+t1070*t97+t544*t872-t555*t71-dqdv(2,7)*t327*t100-t615*t163-t585*t163+t1086*t161-t664*t97-t1059*t74-t559*t178+t1063*t170
    t1154 = dqdv(2,8)*t311
    t1157 = dqdv(3,8)*t140
    t1162 = dqdv(4,8)*t892
    t1167 = dqdv(4,8)*t32
    t1169 = dqdv(2,8)*t11
    t1173 = dqdv(2,8)*t317
    t1175 = dqdv(4,8)*t91
    t1177 = dqdv(3,8)*t25
    t1181 = -t1154*t189-t721*t97-t1157*t277-t758*t163+t784*t895+t703*t872-t1162*t211+t746*t874-dqdv(5,8)*t133*t220-t1167*t115-t1169*t74-t801*t869+t739*t920-t1173*t251-t1175*t75-t1177*t74-t739*t822-t703*t71
    t1182 = dqdv(3,8)*t178
    t1184 = dqdv(5,8)*t327
    t1189 = dqdv(2,8)*t319
    t1192 = dqdv(3,8)*t132
    t1194 = dqdv(5,8)*t191
    t1198 = dqdv(3,8)*t105
    t1200 = dqdv(4,8)*t161
    t1202 = dqdv(4,8)*t71
    t1204 = dqdv(2,8)*t5
    t1208 = dqdv(2,8)*t8
    t1212 = t1182*t19+t1184*t930+t1162*t132-t684*t163+t1173*t121-t1189*t144+t1154*t150-t1192*t22+t1194*t125+t1157*t32+t758*t118+t1198*t125-t1200*t44+t1202*t11-t1204*t211+dqdv(2,8)*t314*t87-t1208*t94+dqdv(2,8)*t38*t869
    t1219 = dqdv(3,8)*t41
    t1228 = dqdv(5,8)*t91
    t1236 = dqdv(4,8)*t115
    t1238 = t788*t895-t788*t178-t1182*t17+t755*t922-t755*t816+t1219*t161+t1192*t238+t750*t872+t1167*t170+dqdv(5,8)*t238*t170-dqdv(2,8)*t17*t170+t1228*t110-t1202*t288-dqdv(5,8)*t173*t261+t1200*t360+t1177*t97+t1189*t261+t1236*t277
    t1239 = dqdv(5,8)*t115
    t1243 = dqdv(5,8)*t74
    t1262 = -t1239*t238-t693*t869-t1219*t94+t1243*t44+t1175*t25+t684*t118+t725*t874+t1208*t161-t750*t71-dqdv(2,8)*t14*t140-t775*t97+t1169*t97+dqdv(2,8)*t19*t170-t763*t170-t771*t170-t784*t178-dqdv(2,8)*t327*t100+t1204*t132
    t1266 = dqdv(5,5)*dqdv(3,1)
    t1271 = dqdv(2,2)*dqdv(3,5)
    t1274 = dqdv(4,5)*dqdv(3,1)
    t1276 = dqdv(4,4)*dqdv(3,5)
    t1278 = dqdv(3,2)*dqdv(5,1)
    t1282 = dqdv(2,1)*dqdv(3,2)
    t1285 = dqdv(4,2)*dqdv(5,1)
    t1287 = -t1266*t10+t825*t10-t814*dqdv(2,2)*dqdv(3,4)-t150*t1271+t833*t324+t1274*t29-t128*t1276-t1278*t238+dqdv(5,2)*dqdv(3,1)*t238+t1282*t17+t814*t147+t1285*t133
    t1291 = dqdv(3,5)*dqdv(4,2)
    t1292 = dqdv(2,4)*dqdv(5,1)
    t1301 = dqdv(5,5)*dqdv(3,2)
    t1304 = t150*t36-dqdv(4,2)*dqdv(3,1)*t277-t1292*t1291+t816*t301-t104*t14-t833*t2+t128*t5+t1266*t61-t1274*t7-t1282*t19-t816*t1301+t104*t356
    t1309 = dqdv(1,1)*dqdv(3,5)
    t1312 = dqdv(4,5)*dqdv(3,2)
    t1314 = dqdv(1,1)*dqdv(3,4)
    t1317 = dqdv(3,1)*dqdv(1,2)
    t1324 = dqdv(3,5)*dqdv(4,1)
    t1330 = t1301*t97-t301*t97+t46*t1309-t315*t869-t1312*t94+t294*t1314+t183*t1276-t1317*t17-t296*t825+t1317*t19+t1312*t161+t296*t1266-t243*t356+t361*t1324-t183*t5-t294*t91+t243*t14-t1301*t74
    t1340 = dqdv(2,1)*dqdv(1,2)
    t1351 = t17*t65-t19*t65-t43*t874-t32*t207+t53*t94+t22*t230-t319*t181+t317*t235-t1340*t17+t183*t238+t32*t243+t319*t104-t317*t128-t53*t161+t1340*t19-t243*t277-t22*t183+t74*t2
    t1356 = dqdv(5,2)*dqdv(3,4)
    t1365 = dqdv(1,2)*dqdv(3,4)
    t1367 = dqdv(1,2)*dqdv(5,4)
    t1374 = t14*t65-t356*t65+t7*t1309-t1356*t874-t147*t869+t40*t874+t161*t1271-t91*t324-t71*t301-t122*t825+t91*t2-t1365*t118-t1367*t158+t1367*t225+t1365*t163-t161*t36+t71*t1301+t122*t1266
    t1379 = dqdv(1,1)*dqdv(3,2)
    t1383 = dqdv(2,2)*dqdv(3,1)
    t1396 = -t1276*t65+t5*t65+t36*t97-t22*t1379+t892*t207-t53*t1314-t317*t1383+t353*t181-t1317*t238+t1340*t1276-t1340*t5+t243*t133+t317*t1282+t53*t91-t36*t74+t22*t1317-t892*t243-t353*t104
    t1399 = dqdv(2,6)*dqdv(3,1)*dqdv(1,2)
    t1401 = dqdv(5,6)*t71
    t1404 = dqdv(2,6)*dqdv(5,1)*dqdv(1,2)
    t1407 = dqdv(5,6)*dqdv(2,1)*dqdv(1,2)
    t1412 = dqdv(3,6)*t319
    t1414 = dqdv(3,6)*t317
    t1425 = dqdv(3,6)*t22
    t1429 = -t1399*t19-t1401*t1291-t1404*t1276+t1407*t1276+t407*t869-t413*t1266-t500*t1314-t1412*t181+t1414*t235-dqdv(4,6)*t14*t65+t413*t825+t996*t324-t416*t1314-t1013*t1271+t500*t91-dqdv(3,6)*t243*t277-t1425*t183+dqdv(3,6)*t17*t65
    t1439 = dqdv(3,2)*t415
    t1441 = dqdv(3,6)*t1340
    t1446 = dqdv(4,6)*dqdv(2,1)*dqdv(1,2)
    t1450 = dqdv(2,6)*dqdv(5,5)
    t1451 = dqdv(3,2)*t1450
    t1459 = dqdv(3,6)*t183*t238-dqdv(5,6)*t1276*t65+t942*t301-dqdv(3,6)*t61*t869+t404*t65+t1011*t1271+t1439*t94-t1441*t17+dqdv(4,6)*t1367*t158-t1446*t356-dqdv(3,6)*t2*t97+t1451*t74+t449*t1309+t468*t869+t453*t97-t1011*t36-t1439*t161+t470*t1282
    t1461 = dqdv(5,6)*t122
    t1466 = dqdv(2,6)*dqdv(4,1)*dqdv(1,2)
    t1469 = dqdv(5,2)*t455
    t1472 = dqdv(5,6)*t243
    t1480 = dqdv(5,6)*t1317
    t1482 = dqdv(4,6)*t1365
    t1486 = -t1461*t1324-t456*t1266+t1013*t36+t1466*t356-t409*t1379-t1469*t1324+t409*t1317+t1472*t133+t1404*t5-t458*t1309+t1399*t17+t1412*t104-t1407*t5+t1425*t230-t1480*t238-t1482*t163+t1441*t19-t1414*t128
    t1502 = dqdv(3,6)*t53
    t1510 = -t942*t1301+dqdv(3,6)*t32*t243+t416*t91+t1482*t118-t1451*t97+t418*t97-dqdv(3,6)*t19*t65-t517*t1309+t442*t65+t456*t825-t1466*t14-t470*t1383-t429*t874-t1502*t161+t1502*t94+dqdv(4,6)*t1356*t874-t996*t2+dqdv(3,6)*t74*t2
    t1514 = dqdv(3,7)*t317
    t1516 = dqdv(5,2)*t595
    t1529 = dqdv(3,7)*t22
    t1536 = dqdv(2,7)*dqdv(3,1)*dqdv(1,2)
    t1538 = dqdv(3,2)*t566
    t1542 = t1514*t235-t1516*t1324+t596*t825+t638*t91-dqdv(3,7)*t2*t97-dqdv(5,7)*t1276*t65+dqdv(3,7)*t17*t65+t553*t65+t559*t65-t1514*t128-t1529*t183+t605*t825+t1065*t301+dqdv(4,7)*t1367*t158-t1536*t19-t1538*t161+t660*t869-t596*t1266
    t1546 = dqdv(3,7)*t319
    t1551 = dqdv(2,7)*dqdv(4,1)*dqdv(1,2)
    t1557 = dqdv(5,7)*t1317
    t1561 = dqdv(4,7)*t1365
    t1564 = dqdv(2,7)*dqdv(5,1)*dqdv(1,2)
    t1570 = -dqdv(3,7)*t19*t65+t623*t97-t1546*t181-t1051*t2+t1536*t17-t1551*t14+t1045*t36+t619*t1317-dqdv(3,7)*t243*t277-t1557*t238-t619*t1379+t1529*t230+t1561*t118-t1564*t1276-t1561*t163-t1045*t1271+t1051*t324+t1047*t1271
    t1576 = dqdv(5,7)*t71
    t1578 = dqdv(3,7)*t1340
    t1585 = dqdv(5,7)*dqdv(2,1)*dqdv(1,2)
    t1588 = dqdv(3,7)*t53
    t1594 = dqdv(5,7)*t243
    t1596 = dqdv(2,7)*dqdv(5,5)
    t1597 = dqdv(3,2)*t1596
    t1601 = dqdv(3,7)*t32*t243+t593*t97+t1551*t356-t1576*t1291+t1578*t19-t1578*t17+dqdv(3,7)*t183*t238+t668*t1282+t1585*t1276-t564*t1309+t1588*t94+dqdv(4,7)*t1356*t874-t588*t874+t1538*t94+t1594*t133-t1597*t97-t1585*t5+t567*t91
    t1602 = dqdv(5,7)*t122
    t1613 = dqdv(4,7)*dqdv(2,1)*dqdv(1,2)
    t1626 = -t1602*t1324-dqdv(4,7)*t14*t65-dqdv(3,7)*t61*t869+t655*t1309+t653*t869-t567*t1314-t1588*t161-t1613*t356+t1597*t74-t641*t1309-t638*t1314-t668*t1383+dqdv(3,7)*t74*t2+t1564*t5-t1065*t1301-t605*t1266-t1047*t36+t1546*t104
    t1634 = dqdv(2,8)*dqdv(4,1)*dqdv(1,2)
    t1638 = dqdv(3,8)*t53
    t1640 = dqdv(2,8)*dqdv(5,5)
    t1641 = dqdv(3,2)*t1640
    t1647 = dqdv(4,8)*dqdv(2,1)*dqdv(1,2)
    t1658 = t1202*t301-t736*t1266+t1175*t324-t1634*t14-dqdv(3,8)*t61*t869-t1638*t161+t1641*t74-t797*t1309+t788*t65-t723*t1266-t1647*t356+t739*t1282-dqdv(3,8)*t243*t277-t680*t874-t766*t1309+dqdv(4,8)*t1367*t158+t723*t825+t706*t97
    t1664 = dqdv(3,8)*t317
    t1667 = dqdv(2,8)*dqdv(3,1)*dqdv(1,2)
    t1669 = dqdv(3,8)*t319
    t1672 = dqdv(5,8)*t122
    t1680 = dqdv(2,8)*dqdv(5,1)*dqdv(1,2)
    t1688 = dqdv(5,8)*dqdv(2,1)*dqdv(1,2)
    t1690 = dqdv(5,8)*t71
    t1692 = t674*t97-dqdv(3,8)*t19*t65-dqdv(5,8)*t1276*t65+t1664*t235+t1667*t17-t1669*t181-t1200*t1271-t1672*t1324-t739*t1383+dqdv(3,8)*t17*t65-dqdv(4,8)*t14*t65+t1680*t5-dqdv(3,8)*t2*t97-t1641*t97+t1638*t94-t1202*t1301-t1688*t5-t1690*t1291
    t1695 = dqdv(5,2)*t735
    t1698 = dqdv(3,2)*t742
    t1708 = dqdv(3,8)*t22
    t1711 = dqdv(5,8)*t243
    t1713 = dqdv(4,8)*t1365
    t1717 = dqdv(5,8)*t1317
    t1720 = -t1664*t128-t1695*t1324-t1667*t19-t1698*t161-t1175*t2+dqdv(3,8)*t74*t2+dqdv(4,8)*t1356*t874+t1634*t356-t763*t1379+t763*t1317+t1708*t230+t1688*t1276+t1711*t133-t1713*t163+t1713*t118-t1708*t183-t1717*t238-t761*t1314
    t1722 = dqdv(3,8)*t1340
    t1743 = t691*t869-t1722*t17+t736*t825+dqdv(3,8)*t32*t243+t1200*t36-t1680*t1276+t779*t1309+t686*t869-t743*t1314-t1243*t36+t784*t65+dqdv(5,8)*t353*t181+dqdv(3,8)*t183*t238+t761*t91+t1669*t104+t1722*t19+t743*t91+t1698*t94
    t1747 = dqdv(5,5)*dqdv(4,1)
    t1751 = dqdv(2,2)*dqdv(4,5)
    t1762 = t1747*t16-t814*t16+t825*t13+t121*t1751-t144*t324-t1324*t27+t128*t379-dqdv(5,2)*dqdv(4,1)*t360-t104*t288+t1285*t360-t1278*t110-t1747*t78
    t1764 = dqdv(2,3)*dqdv(5,1)
    t1771 = dqdv(5,5)*dqdv(4,3)
    t1776 = dqdv(5,5)*dqdv(4,2)
    t1779 = -t128*t8+t1764*t1312+t144*t2-t1282*t330-t825*t48+t104*t11-t121*t53+t1282*t1771-t822*t294+dqdv(3,2)*dqdv(4,1)*t75+t822*t1776+t1324*t4
    t1786 = dqdv(1,1)*dqdv(4,5)
    t1801 = t294*t170-t1776*t170+t1291*t211+t31*t869-t50*t1786-t301*t220-t183*t379+t243*t288-t298*t1747+t298*t814+t1317*t330+t301*t140+t183*t8-t1317*t1771-t243*t11-t1291*t132-t334*t1274+t1776*t115
    t1808 = dqdv(5,2)*dqdv(4,3)
    t1812 = dqdv(1,2)*dqdv(4,3)
    t1815 = dqdv(1,2)*dqdv(5,3)
    t1824 = t1771*t65-t330*t65+t4*t1786-t48*t869+t336*t874-t1808*t874-t140*t324+t132*t1751+t1812*t163-t1812*t118-t1815*t125+t140*t2-t132*t53-t178*t294-t88*t814+t1815*t68+t178*t1776+t88*t1747
    t1844 = t288*t65-t11*t65-t21*t874-t25*t1379+t44*t230+t36*t211+t311*t235-t314*t1383-t1340*t288+t183*t360-t36*t132+t1340*t11-t311*t128+t314*t1282-t1317*t75+t2*t115+t25*t1317-t44*t183
    t1864 = t379*t65-t8*t65-t53*t170+t36*t220-t41*t1379+t44*t207-t327*t1383+t311*t181-t1340*t379+t243*t360-t1317*t110+t327*t1282+t53*t115+t41*t1317+t1340*t8-t311*t104-t44*t243-t36*t140
    t1872 = dqdv(4,2)*t1450
    t1876 = dqdv(3,6)*t1812
    t1881 = dqdv(4,6)*t311
    t1884 = dqdv(5,2)*t420
    t1890 = t529*t181+t964*t1751-t948*t1751-t421*t814-t418*t140+t1466*t11-t1872*t115-dqdv(3,6)*t1815*t125-t1876*t118-t453*t140-t1404*t8-t1407*t379+t1881*t128+t1446*t288+t1884*t1274+t1872*t170+dqdv(4,6)*t11*t65-t416*t170
    t1912 = dqdv(5,6)*t88
    t1914 = -t500*t170+dqdv(4,6)*t2*t170+t418*t220-t434*t869+t453*t220-t964*t53+t508*t1786+t475*t1786-t529*t104-dqdv(4,6)*t183*t360-t1466*t288+t1404*t379+t1472*t360+dqdv(5,6)*t379*t65-dqdv(4,6)*t288*t65-t439*t65+t402*t1747+t1912*t1274
    t1917 = dqdv(4,2)*t452
    t1919 = dqdv(5,6)*t178
    t1922 = dqdv(4,6)*t44
    t1925 = dqdv(4,6)*t314
    t1928 = dqdv(4,6)*t36
    t1941 = -t1446*t11+t1917*t132+t1919*t1312+t421*t1747+t1922*t183-t984*t324+t1925*t1383-t1881*t235+t1928*t132-t481*t243+dqdv(3,6)*t1771*t65+t948*t53+t1876*t163+t984*t2-t447*t65+t1399*t1771+dqdv(4,6)*t1317*t75-t936*t294
    t1957 = dqdv(4,6)*t25
    t1962 = t1441*t330-t1480*t110-t1925*t1282+t1407*t8-t402*t814+t936*t1776-t994*t2+t461*t874-t505*t869-t1399*t330+t481*t207-t1928*t211-t493*t1786-t1922*t230-t1917*t211+t1957*t1379-dqdv(3,6)*t1808*t874-t1957*t1317
    t1966 = dqdv(4,7)*t44
    t1968 = dqdv(5,7)*t178
    t1973 = dqdv(3,7)*t1812
    t1980 = dqdv(4,7)*t36
    t1989 = dqdv(4,7)*t314
    t1991 = -t1966*t230+t1968*t1312+t538*t1786-t551*t869+t593*t220-t1973*t118+t1088*t53-t1551*t288+t1613*t288+t1564*t379+t583*t1747+t1980*t132-t1613*t11-t623*t140-dqdv(4,7)*t183*t360-t1084*t53+dqdv(4,7)*t2*t170+t1989*t1383
    t2004 = dqdv(4,2)*t622
    t2014 = -t638*t170+t601*t1786-t1980*t211+t571*t207-t1116*t324+dqdv(5,7)*t379*t65-dqdv(4,7)*t288*t65+dqdv(3,7)*t1771*t65+t580*t181+t2004*t132+t1116*t2-t1557*t110-t1124*t2+t1551*t11+t1966*t183-t571*t243+t1061*t1776+t631*t1747
    t2019 = dqdv(4,2)*t1596
    t2037 = -t1564*t8+t1973*t163+t1585*t8-t2019*t115-dqdv(3,7)*t1815*t125-t1989*t1282-t569*t869+t633*t874-t1585*t379-t626*t1786-dqdv(3,7)*t1808*t874+t1084*t1751-t583*t814+t1578*t330-t1061*t294-t1536*t330-t631*t814-t580*t104
    t2038 = dqdv(4,7)*t311
    t2046 = dqdv(5,7)*t88
    t2052 = dqdv(4,7)*t25
    t2054 = dqdv(5,2)*t582
    t2062 = t2038*t128-t2004*t211+t623*t220-t2038*t235-t593*t140-t555*t65+t2019*t170+t2046*t1274-t567*t170+dqdv(4,7)*t1317*t75+t1536*t1771-t2052*t1317+t2054*t1274+t2052*t1379-t544*t65+dqdv(4,7)*t11*t65-t1088*t1751+t1594*t360
    t2073 = dqdv(4,2)*t705
    t2084 = dqdv(3,8)*t1812
    t2086 = dqdv(4,8)*t25
    t2088 = t1192*t1751+t755*t181-t706*t140-t719*t814+t1647*t288-t696*t814+t674*t220+t2073*t132-t768*t869-t1157*t324-dqdv(3,8)*t1808*t874-t1688*t379-t1680*t8-t703*t65+t696*t1747-t743*t170-t2084*t118+t2086*t1379
    t2089 = dqdv(4,2)*t1640
    t2095 = dqdv(4,8)*t44
    t2104 = dqdv(4,8)*t314
    t2110 = dqdv(4,8)*t311
    t2112 = dqdv(5,8)*t88
    t2115 = t2089*t170-t750*t65-t1717*t110+dqdv(3,8)*t1771*t65-t2095*t230-t775*t243-t732*t869+dqdv(5,8)*t379*t65-t1236*t2+t1184*t1282-t761*t170-t2104*t1282-t674*t140+dqdv(4,8)*t1317*t75+t1688*t8-t2110*t235+t2112*t1274+t775*t207
    t2123 = dqdv(5,2)*t695
    t2125 = dqdv(4,8)*t36
    t2139 = t1634*t11-t682*t1786-t2086*t1317-dqdv(4,8)*t288*t65+t2095*t183+t2123*t1274-t2125*t211+t2104*t1383+t777*t1786-t1634*t288-t755*t104+t1667*t1771+t2110*t128+t1182*t1776+t1157*t2+t719*t1747-t2089*t115-dqdv(3,8)*t1815*t125
    t2161 = t706*t220-t1182*t294-t1647*t11+t793*t874+t1239*t53-t1667*t330+t1722*t330+dqdv(4,8)*t11*t65-t1192*t53+t2125*t132-t1184*t1383+t2084*t163+t1680*t379-t2073*t211+t709*t1786-dqdv(4,8)*t183*t360+dqdv(4,8)*t2*t170+t1711*t360
    t2177 = -t189*t16+t150*t16+t121*t10-t144*t29+t842*t13-t833*t27-t104*t59+t1292*t57-t816*t21+t128*t1+t822*t46+t144*t7
    t2191 = t833*t4-t1292*t31-t822*t43-t121*t61+t104*t52-dqdv(2,3)*dqdv(4,1)*t40+t1282*t55-t1282*t35-t842*t48+t816*t50-t128*t38+t1764*t24
    t2206 = t43*t170-t46*t170-t50*t97+t315*t211-t1356*t220+t31*t94+t296*t251-t183*t1+t243*t59-t361*t261-t296*t121-t31*t161
    t2219 = t348*t100+t183*t38+t361*t144+t50*t74-t298*t150-t1317*t55+t1356*t140-t334*t87-t243*t52-t315*t132+t298*t189+t1317*t35
    t2234 = -t35*t65+t55*t65+t244*t207-t184*t230-t48*t94+t4*t97-t140*t29+t132*t10-t129*t235+t105*t181-t178*t43+t178*t46
    t2248 = -t88*t189+t1815*t922+t184*t183+t48*t161-t105*t104-t1812*dqdv(5,4)*dqdv(2,1)+t129*t128+t88*t150-t132*t61-t244*t243+t140*t7-t4*t74
    t2256 = dqdv(2,3)*dqdv(5,4)
    t2264 = t59*t65-t52*t65-t7*t170+t147*t211+t208*t230-t2256*t1379-t161*t16+t191*t235-t906*t1383+t91*t27+t71*t21+t122*t251
    t2278 = -t1367*t920+t1365*dqdv(5,3)*dqdv(2,1)-t208*t183-t91*t4+t161*t78-t191*t128-t147*t132+t906*t1282-t122*t121-t71*t50+t2256*t1317+t7*t115
    t2294 = -t38*t65+t1*t65+t208*t207+t147*t220-t61*t170-t824*t1379-t300*t1383+t91*t13+t191*t181-t74*t16-t208*t243+t1365*dqdv(4,3)*dqdv(2,1)
    t2307 = -t122*t144-t191*t104+t74*t78+t71*t57+t300*t1282+t122*t261+t824*t1317+t61*t115-t91*t48-t71*t31-t147*t140-t289*t920
    t2318 = dqdv(3,6)*t244
    t2323 = t496*t183+t936*t46-t956*t48+t402*t150+t434*t161+t421*t150+t1399*t55+t1919*t24-t2318*t243-t475*t74+t458*t170-t1469*t144
    t2326 = dqdv(5,2)*t406
    t2330 = dqdv(5,6)*t147
    t2340 = -t994*t7-t936*t43-t2326*t140+t413*t121+t1446*t59-t2330*t140-dqdv(4,6)*t88*t100+t1884*t87+t505*t161+dqdv(5,6)*t1*t65-t505*t94+t468*t132
    t2355 = t493*t74+t956*t13+t456*t121-t488*t65-t402*t189+t948*t61-t493*t97-dqdv(4,6)*t178*t40-t1446*t52+t1031*t181-t407*t211+t2326*t220
    t2359 = dqdv(4,6)*t78
    t2369 = t2318*t207+t1441*t35+t508*t97+t2359*t94-t984*t29+t1029*t181+t984*t7+t498*t183-t1466*t59-t942*t21+t1404*t1-t456*t251
    t2387 = -t1011*t16-dqdv(4,4)*dqdv(2,6)*dqdv(5,2)*t170-t964*t61+t436*t128+t994*t29-t2359*t161+t1461*t261+t445*t128-dqdv(4,6)*t59*t65+t525*t65+t1469*t261-t1461*t144
    t2400 = t996*t4-t496*t230-t1399*t35+t1912*t87-t508*t74+t407*t132+t1401*t57-t1407*t1-t413*t251+t520*t65+t1013*t16-t498*t230
    t2411 = dqdv(5,6)*t208
    t2416 = t517*t170-t1441*t55+t2330*t220-t1029*t104-dqdv(3,6)*t35*t65-t421*t189-t996*t27-t436*t235-t2411*t243-t1404*t38+t1407*t38-t1401*t31
    t2430 = -dqdv(4,2)*t420*t100+t1466*t52+t2411*t207-t468*t211-t445*t235+t964*t10+t475*t97-t948*t10-t1031*t104-t449*t170+t942*t50-t434*t94
    t2435 = dqdv(5,7)*t147
    t2439 = dqdv(5,7)*t208
    t2447 = dqdv(4,7)*t78
    t2452 = t2435*t220+t626*t74+t615*t128-t2439*t243-t660*t211-t1116*t29+t585*t128+t1585*t38+t631*t150+t641*t170-t2447*t161-dqdv(4,4)*dqdv(2,7)*dqdv(5,2)*t170
    t2465 = t596*t121+t2439*t207-t601*t74-t603*t230+t1051*t4+t660*t132-t1136*t104-t1053*t104+t1061*t46-t1578*t55-t1119*t48+t1088*t61
    t2483 = -t615*t235+dqdv(5,7)*t1*t65+t583*t150-t1088*t10-t1051*t27-dqdv(4,2)*t582*t100+t653*t132+t1564*t1-t605*t251+t1516*t261-dqdv(3,7)*t35*t65-dqdv(4,7)*t59*t65
    t2484 = dqdv(5,2)*t652
    t2498 = t2484*t220+t1576*t57+t605*t121+t548*t183-dqdv(4,7)*t178*t40+t1065*t50-t585*t235+t1578*t35-t551*t94+t569*t161-t649*t65+t557*t65
    t2507 = dqdv(3,7)*t244
    t2515 = -t631*t189-t1124*t7+t1613*t59-t653*t211-dqdv(4,7)*t88*t100+t2507*t207+t561*t65-t548*t230-t538*t74-t569*t94-t1536*t35-t2435*t140
    t2528 = t1053*t181-t1047*t16+t1124*t29+t1968*t24-t1613*t52+t1136*t181-t2484*t140+t1119*t13+t1602*t261+t1116*t7-t2507*t243+t551*t161
    t2542 = t2046*t87+t1084*t10-t626*t97-t1084*t61-t1061*t43-t1576*t31-t1564*t38+t2447*t94-t1602*t144+t538*t97-t1065*t21+t2054*t87
    t2555 = t603*t183-t1516*t144+t1551*t52-t1585*t1+t1045*t16-t583*t189+t564*t170-t655*t170-t596*t251+t601*t97-t1551*t59+t1536*t55
    t2565 = dqdv(3,8)*t244
    t2576 = -t1672*t144-dqdv(4,2)*t695*t100-dqdv(4,8)*t178*t40-t2565*t243-t1198*t104-t1647*t52+t1688*t38-dqdv(4,8)*t88*t100-t1202*t21+t797*t170-t777*t74+t686*t132
    t2580 = dqdv(5,2)*t690
    t2590 = t1228*t13-t1688*t1+t1175*t4+t2580*t220+t732*t161+t1722*t35+t1680*t1-t686*t211-t691*t211-t709*t74-t1182*t43+t1182*t46
    t2600 = dqdv(5,8)*t147
    t2607 = t2123*t87-t732*t94+t1200*t16-t682*t97+t723*t121+t1239*t61-t1667*t35-t736*t251+t2600*t220+t691*t132+dqdv(5,8)*t178*t24-dqdv(4,8)*t59*t65
    t2620 = t1647*t59+t684*t128-t801*t230+t696*t150+t682*t74-t768*t94+t709*t97+t766*t170-t723*t251+t758*t128-t1695*t144+t1634*t52
    t2623 = dqdv(4,8)*t78
    t2633 = dqdv(5,8)*t208
    t2637 = t2623*t94-t1722*t55+t2565*t207+t719*t150+t768*t161+t1667*t55+t1690*t57+t801*t183-t696*t189+t2633*t207-t779*t170+t1157*t7
    t2651 = -t1680*t38-t1634*t59-t684*t235-t2580*t140+t777*t97-t1175*t27+t1236*t29-t2600*t140+t746*t65+t2112*t87+dqdv(5,8)*t1*t65-t719*t189
    t2668 = t1672*t261-t1157*t29-dqdv(3,8)*t35*t65+t725*t65-dqdv(4,4)*dqdv(2,8)*dqdv(5,2)*t170+t1194*t181+t1198*t181+t736*t121-t758*t235-t693*t230-t1239*t10-t1194*t104
    t2681 = t1192*t10+t1695*t261-t1243*t16-t1192*t61-t1228*t48-t1236*t7-t2623*t161-t676*t65+t1202*t50+t693*t183-t1690*t31-t2633*t243
    !
    dvdq(1,1) = t286*(t34+t63)
    dvdq(1,2) = -t286*t322
    dvdq(1,3) = t286*t350
    dvdq(1,4) = -t286*t377
    dvdq(1,5) = t286*t398
    dvdq(1,6) = -t286*(t438+t472+t503+t533)
    dvdq(1,7) = -t286*(t577+t610+t643+t670)
    dvdq(1,8) = -t286*(t713+t748+t781+t807)
    dvdq(1,9) = 0
    dvdq(2,1) = -t286*(t829+t844)
    dvdq(2,2) = t286*t865
    dvdq(2,3) = -t286*t888
    dvdq(2,4) = t286*t911
    dvdq(2,5) = -t286*t934
    dvdq(2,6) = t286*(t962+t990+t1017+t1041)
    dvdq(2,7) = t286*(t1077+t1104+t1129+t1150)
    dvdq(2,8) = t286*(t1181+t1212+t1238+t1262)
    dvdq(2,9) = 0
    dvdq(3,1) = t286*(t1287+t1304)
    dvdq(3,2) = -t286*t1330
    dvdq(3,3) = t286*t1351
    dvdq(3,4) = -t286*t1374
    dvdq(3,5) = t286*t1396
    dvdq(3,6) = -t286*(t1429+t1459+t1486+t1510)
    dvdq(3,7) = -t286*(t1542+t1570+t1601+t1626)
    dvdq(3,8) = -t286*(t1658+t1692+t1720+t1743)
    dvdq(3,9) = 0
    dvdq(4,1) = -t286*(t1762+t1779)
    dvdq(4,2) = t286*t1801
    dvdq(4,3) = -t286*t1824
    dvdq(4,4) = t286*t1844
    dvdq(4,5) = -t286*t1864
    dvdq(4,6) = t286*(t1890+t1914+t1941+t1962)
    dvdq(4,7) = t286*(t1991+t2014+t2037+t2062)
    dvdq(4,8) = t286*(t2088+t2115+t2139+t2161)
    dvdq(4,9) = 0
    dvdq(5,1) = t286*(t2177+t2191)
    dvdq(5,2) = -t286*(t2206+t2219)
    dvdq(5,3) = t286*(t2234+t2248)
    dvdq(5,4) = -t286*(t2264+t2278)
    dvdq(5,5) = t286*(t2294+t2307)
    dvdq(5,6) = -t286*(t2323+t2340+t2355+t2369+t2387+t2400+t2416+t2430)
    dvdq(5,7) = -t286*(t2452+t2465+t2483+t2498+t2515+t2528+t2542+t2555)
    dvdq(5,8) = -t286*(t2576+t2590+t2607+t2620+t2637+t2651+t2668+t2681)
    dvdq(5,9) = 0
    dvdq(6,1) = 0
    dvdq(6,2) = 0
    dvdq(6,3) = 0
    dvdq(6,4) = 0
    dvdq(6,5) = 0
    dvdq(6,6) = 1
    dvdq(6,7) = 0
    dvdq(6,8) = 0
    dvdq(6,9) = 0
    dvdq(7,1) = 0
    dvdq(7,2) = 0
    dvdq(7,3) = 0
    dvdq(7,4) = 0
    dvdq(7,5) = 0
    dvdq(7,6) = 0
    dvdq(7,7) = 1
    dvdq(7,8) = 0
    dvdq(7,9) = 0
    dvdq(8,1) = 0
    dvdq(8,2) = 0
    dvdq(8,3) = 0
    dvdq(8,4) = 0
    dvdq(8,5) = 0
    dvdq(8,6) = 0
    dvdq(8,7) = 0
    dvdq(8,8) = 1
    dvdq(8,9) = 0
    dvdq(9,1) = 0
    dvdq(9,2) = 0
    dvdq(9,3) = 0
    dvdq(9,4) = 0
    dvdq(9,5) = 0
    dvdq(9,6) = 0
    dvdq(9,7) = 0
    dvdq(9,8) = 0
    dvdq(9,9) = 1
    !
    !test = MATMUL( dVdQ, dQdV ) 
    !DO i = 1, nVar
    !    test(i,i) = test(i,i)-1.0
    !ENDDO
    !error = MAXVAL(ABS(test)) 
    !
    dVdt(1:5) = MATMUL(dVdQ(1:5,:),dQdt) 
    !
    CONTINUE
    !
#endif     
    !
END SUBROUTINE PDEdVdtPrim 