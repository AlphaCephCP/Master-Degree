!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDEFluxGRMHD(F,Q,par)
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE
    INTEGER, PARAMETER :: nVar = 19 
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam)  
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    INTEGER :: iErr 
    INTEGER :: i,j,k,m
    REAL :: V(nVar) 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3)  
    REAL :: gamma1, rho, vx, vy, vz, bx, by, bz, ex, ey, ez, v2, b2, e2, lf, w, ww, uem, wwx, wwy, wwz 
    REAL :: gp, gm, g_cov(3,3), g_contr(3,3), vxB(3), vxB_contr(3), BV(3), BQ(3), vtr(3), vf(3), lapse, shift(3)    
    REAL :: Fij(3,3), vf_cov(3), Qv_contr(3), QB_contr(3), Bv_contr(3), psi 
    !
    F = 0.0 
    !
  CALL PDECons2PrimGRMHD(V,Q,iErr)
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
END SUBROUTINE PDEFluxGRMHD
!
!
! Conservative part of the PDE ( flux tensor F(Q) ) 
!
SUBROUTINE PDESourceGRMHD(S,Q,par,time)
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    INTEGER, PARAMETER :: nVar = 19 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), time   
    REAL, INTENT(OUT) :: S(nVar) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    REAL :: k1, k2, fff, ggg, e, ds, c, xi, sk, sknl, alpha, fa, k0, beta0(3), eta 
    REAL :: g_cov(3,3), det, g_contr(3,3), Christoffel(3,3,3), dgup(3,3,3), faa, b0(3), dChristoffelSrc(3,3,3,3), RiemannSrc(3,3,3,3)
    REAL :: RicciSrc(3,3), gammaup(3), dtmetric(3,3), dtgup(3,3), dtChristoffelSrc(3,3,3), dtGammaUpSrc(3)
    REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25
    !
    S = 0.0     
    !
END SUBROUTINE PDESourceGRMHD
!
! Nonconservative part of the PDE ( B(Q) * gradQ ) 
!    
SUBROUTINE PDENCPGRMHD(BgradQ,Q,gradQ,par)
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE
    INTEGER, PARAMETER :: nVar = 19 
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d), par(nParam)  
    REAL, INTENT(OUT) :: BgradQ(nVar) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    REAL :: k1,k2,fff,ggg,e,c,ds,xi,sk,sknl,g_cov(3,3),g_contr(3,3),Christoffel(3,3,3),dgup(3,3,3)
    REAL :: det, alpha, fa, k0, dk0(3), beta0(3), b0(3), u(3), ialpha 
    REAL :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 
    REAL :: AQx(nVar), BQy(nVar), CQz(nVar) 
    REAL :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL :: v2,vf(3),uem,b2,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim    
    INTEGER :: i,j,k,l,m,iErr, count    
    !
    BgradQ = 0.0 
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
  CALL PDECons2PrimGRMHD(Vc,Q,iErr)
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
END SUBROUTINE PDENCPGRMHD      

SUBROUTINE PDEMatrixBGRMHD(Bn,Q,nv,par)
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE
    INTEGER, PARAMETER :: nVar = 19 
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d)   
    REAL, INTENT(OUT) :: Bn(nVar,nVar) 
    ! Local variables 
    INTEGER :: i 
    REAL    :: p, irho, lam, mu, ialpha  
    REAL    :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)  
    REAL    :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar)  
    REAL    :: k1, k2, fff, ggg, e, ds, cs, xi, sk, sknl, alpha, fa, k0, beta0(3), b0(3)   
    REAL    :: g_cov(3,3), g_contr(3,3), det, Christoffel(3,3,3), dgup(3,3,3), uv(3), bb2   
    REAL    :: lapse, shift(3), gammaij(6), delta(3,3), bv(3), vxb(3), vxb_contr(3), psi, qv_contr(3), qb_contr(3), bv_contr(3) 
    REAL    :: v2,vf(3),uem,e2,gp,gm,vc(nVar),lf,w,ww,gamma1,rho,vf_cov(3), s_contr(3), w_ij, wim    
    INTEGER :: j,k,l,m,iErr, count    
    !
    Bn = 0.0
    !
  !psi = Q(9)
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
  CALL PDECons2PrimGRMHD(Vc,Q,iErr)
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  rho    = Vc(1)
  vf_cov = Vc(2:4)
  p      = Vc(5)
  !
  BV(1:3) = Vc(6:8)
  psi = Vc(9) 
  !
  vf     = MATMUL(g_contr,vf_cov)
  S_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  vxB_contr = MATMUL(g_contr,vxB(1:3))
  BV_contr = MATMUL(g_contr,BV(1:3))
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  bb2    = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(bb2 + e2) 
  !
  !gv_contr = MATMUL(g_contr,gv)
  lf     = 1.0/sqrt(1.0 - v2)
  w      = rho + gamma1*p   ! rho*hentalpy
  ww     = w*lf**2          ! rho*hentalpy*Lorentz^2 
  !
  !DO j=1,3
  A = 0.
  B = 0.
  C = 0.
    !lapse
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=1
    !------ 
    A(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    A(5,10) =  S_contr(j)     !  Q(10) or lapse
    !
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=2
    !------ 
    B(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    B(5,10) =  S_contr(j)     !  Q(10) or lapse
    ! 
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    j=3
    !------
    C(1+j,10) = + (Q(5)+Q(1))   ! Q(10) or lapse
    C(5,10) =  S_contr(j)     !  Q(10) or lapse
    ! 
    count=0
    DO i=1,3
        ! shift
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=1
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        A(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        A(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        !
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=2
        !------ 
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        B(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        B(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        j=3
        !------
        W_ij = ww*vf_cov(i)*vf(j)-vxB(i)*vxB_contr(j)-BV(i)*BV_contr(j)+(p+uem)*delta(i,j)
        !
        C(1+j,10+i) = - Q(1+i)  ! Q(11:13)  shift(i) or shift_contr(i)
        C(5,10+i) = - gp*W_ij    ! Q(11:13)  shift(i) or shift_contr(i)
        ! 
          DO m=1,3
            IF(m.GE.i) THEN  
                !metric
                count=count+1
                !
                Wim = ww*vf(i)*vf(m)-vxB_contr(i)*vxB_contr(m)-BV_contr(i)*BV_contr(m)+(p+uem)*g_contr(i,m)
                Wim = Wim + (1.0 - delta(i,m))*(ww*vf(m)*vf(i)-vxB_contr(m)*vxB_contr(i)-BV_contr(m)*BV_contr(i)+(p+uem)*g_contr(m,i))  ! account also of the remaining symmetric components of gamma for i.NE.m.
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=1
                !------ 
                A(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                A(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                !
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=2
                !------ 
                B(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                B(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                ! 
                !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                j=3
                !------
                C(1+j,13+count) = - 0.5*gp*lapse*Wim  ! Q(14:19) gammaij(count) or  g_cov(i,m)
                C(5,13+count) = - 0.5*gp*Wim*shift(j)   ! Q(14:19) gammaij(count) or  g_cov(i,m)
                ! 
            ENDIF
          ENDDO
      ENDDO
  !ENDDO 
  !
  Bn = A*nv(1) + B*nv(2) + C*nv(3) 
  !  
END SUBROUTINE PDEMatrixBGRMHD  
    

SUBROUTINE PDECons2PrimGRMHD(V,Q,iErr)
    USE typesDef, ONLY : nParam, EQN 
    IMPLICIT NONE
    INTEGER, PARAMETER :: nVar = 19 
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
    !    
    iErr = 0     
    !
  psi = Q(9)
  lapse = Q(10)
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
  CALL MatrixInverse3x3(g_cov,g_contr,gp)
  !
  gp = SQRT(gp)
  !
  gm = 1./gp
  ! 
  !CALL METRIC(x, lapse, gp, gm, shift, g_cov, g_contr)  
  Qloc(1:8)  = gm*Q(1:8)
  !
  !BV(1:3) = Qloc(6:8) ! !wrong: magnetic field is contravariant
  BV_contr(1:3) = Qloc(6:8) ! !wrong: magnetic field is contravariant
  !
  gamma1 = EQN%gamma/(EQN%gamma - 1.0)
  gam    = 1.0/gamma1
  ! Solve for p
  FAILED  = .FALSE.
  dd      = Qloc(1)
  sm_cov  = Qloc(2:4)
  !
  sm   = MATMUL (g_contr, sm_cov)
  !BV_contr   = MATMUL (g_contr, BV)
  BV   = MATMUL (g_cov, BV_contr)
  s2   = sm_cov(1)*sm(1) + sm_cov(2)*sm(2) + sm_cov(3)*sm(3)
  b2   = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  sb   = sm_cov(1)*BV_contr(1) + sm_cov(2)*BV_contr(2) + sm_cov(3)*BV_contr(3)
  sb2  = sb**2
  eps  = 1.e-10 !8
  !!! RHD
  !!e       = Qloc(5) 
  !!x1   = eps      ! min pressure
  !!x2   = 1.0d+5   ! max pressure
  !!p    = RTSAFE_C2P_RHD1(x1,x2,tol,d,e,s2,FAILED)
  !!IF (FAILED) THEN
  !!   p = 1.0e-20
  !!ENDIF
  !!rho  = d / (e + p + d) * SQRT((e + p + d)**2 - s2)
  !!den  = 1.0 / (e + p + d)
  !!!
  !!vf_cov(1:3) = sm_cov(1:3)*den
  !!
  ! First option [Del Zanna et al. (2007) A&A, 473, 11-30 (method 3)]
  e    = Qloc(5) + dd  ! Q(5) = gamma^1/2 ( U - D )
  x1   = 0.      ! 
  x2   = 1.0-eps ! 
  w=0 
  !
  v2   = RTSAFE_C2P_RMHD1(x1,x2,tol,gam,dd,e,s2,b2,sb2,w,FAILED) 
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
     rho  = dd*sqrt(1.-v2)
     vf_cov(1) = (sm_cov(1) + vb*BV(1))*den
     vf_cov(2) = (sm_cov(2) + vb*BV(2))*den
     vf_cov(3) = (sm_cov(3) + vb*BV(3))*den
     p = max(1.e-15, gam*(w*(1.-v2)-rho))
  ENDIF
  !
  V(1:19) = (/ rho, vf_cov(1:3), p, BV_contr(1:3), psi , lapse, shift(1:3), gammaij(1:6)/)
  !  
END SUBROUTINE PDECons2PrimGRMHD  

    
    
    
SUBROUTINE PDEPrim2ConsGRMHD(Q,V)
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE
    INTEGER, PARAMETER    :: nVar = 19 
    ! Argument list 
    REAL, INTENT(IN)      :: V(nVar)     ! primitive variables 
    REAL, INTENT(OUT)     :: Q(nVar)     ! vector of conserved quantities 
    ! Local variables 
    REAL                  :: rho, vx, vy, vz, p, bx, by, bz, ex, ey, ez, v2, b2, e2    
    REAL                  :: lf, gamma1, w, ww, uem 
    REAL                  :: gp, gm, g_cov(3,3), g_contr(3,3), bv(3), vf_cov(3), psi, lapse, shift(3)
    REAL                  :: vf(3), bv_contr(3), qv_contr(3), qb_contr(3), vxb(3), vb_cov, b2_cov, vxb_contr(3) 
    ! 
  rho     = V(1)
  vf_cov  = V(2:4)
  p       = V(5)
  !
  !BV(1:3) = V(6:8)  !wrong
  BV_contr(1:3) = V(6:8)  !wrong
  psi = V(9)
  lapse = V(10)
  shift = V(11:13)          ! NB: we choose V() and Q() being the shift_controvariant!
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
  vf      = MATMUL(g_contr,vf_cov)
  BV = MATMUL(g_cov,BV_contr(1:3))
  Qv_contr = MATMUL(g_contr,Q(2:4))
  QB_contr = MATMUL(g_contr,Q(6:8))
  vxB(1) = vf_cov(2)*BV(3) - vf_cov(3)*BV(2)
  vxB(2) = vf_cov(3)*BV(1) - vf_cov(1)*BV(3)
  vxB(3) = vf_cov(1)*BV(2) - vf_cov(2)*BV(1)
  !vxB_contr = MATMUL(g_contr,vxB(1:3))   ! WRONG!
  vxB_contr(1) = vf(2)*BV_contr(3) - vf(3)*BV_contr(2)
  vxB_contr(2) = vf(3)*BV_contr(1) - vf(1)*BV_contr(3)
  vxB_contr(3) = vf(1)*BV_contr(2) - vf(2)*BV_contr(1)
  !
  !v2     = vx**2 + vy**2 + vz**2
  !b2     = bx**2 + by**2 + bz**2
  !e2     = ex**2 + ey**2 + ez**2 
  vb_cov      = vf_cov(1)*BV_contr(1) + vf_cov(2)*BV_contr(2) + vf_cov(3)*BV_contr(3)    ! VERIFY! 
  v2     = vf(1)*vf_cov(1) + vf(2)*vf_cov(2) + vf(3)*vf_cov(3)
  e2     = vxB_contr(1)*vxB(1) + vxB_contr(2)*vxB(2) + vxB_contr(3)*vxB(3)
  b2     = BV_contr(1)*BV(1) + BV_contr(2)*BV(2) + BV_contr(3)*BV(3)
  !
  uem    = 0.5*(b2 + e2) 
  !
  !
  b2_cov      = SUM(BV(1:3)**2)
  !b2_contr      = SUM(BV_contr(1:3)**2)
  !
  IF (v2.GT.1.0) THEN
     WRITE(*,*)'Superluminal velocity in PDEPrim2Cons!!'
     STOP
  ENDIF
  lf     = 1.0 / sqrt(1.0 - v2)
  gamma1 = EQN%gamma/(EQN%gamma-1.0)
  w      = rho + gamma1*p
  ww     = w*lf**2
  !
  Q(1)    = rho*lf
  Q(2:4)  = ww*vf_cov(1:3) + b2*vf_cov(1:3) - vb_cov*BV(1:3) ! covariant!!!!
  !Q(2:4)  = ww*vf(1:3) + b2_contr*vf(1:3) - vb_contr*BV_contr(1:3) ! contravariant
  Q(5)    = ww - p + uem - Q(1)     !!!!! we subtract PDE(Q(1))!!!!
  Q(6:8) = V(6:8)
  Q(1:8)    = gp*Q(1:8) 
  !
  Q(9:nVar) = V(9:nVar)
    !
END SUBROUTINE PDEPrim2ConsGRMHD       
            
    
FUNCTION RTSAFE_C2P_RMHD1(X1,X2,XACC,gam,d,e,s2,b2,sb2,w,FAILED)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=200
  INTEGER               :: J
  REAL                  :: RTSAFE_C2P_RMHD1
  REAL                  :: X1,X2,XACC,gam,d,e,s2,b2,sb2,w
  REAL                  :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP
  LOGICAL               :: FAILED
  !
  FAILED = .FALSE.
  CALL FUNC_C2P_RMHD1(X1,FL,DF,gam,d,e,s2,b2,sb2,w)
  IF(FL.EQ.0.) THEN
     RTSAFE_C2P_RMHD1=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RMHD1(X2,FH,DF,gam,d,e,s2,b2,sb2,w)
  IF(FH.EQ.0.) THEN
     RTSAFE_C2P_RMHD1=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     RTSAFE_C2P_RMHD1 = 0.
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  RTSAFE_C2P_RMHD1=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1,F,DF,gam,d,e,s2,b2,sb2,w)
  DO 11 J=1,MAXIT
     IF(((RTSAFE_C2P_RMHD1-XH)*DF-F)*((RTSAFE_C2P_RMHD1-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        RTSAFE_C2P_RMHD1=XL+DX
        IF(XL.EQ.RTSAFE_C2P_RMHD1)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=RTSAFE_C2P_RMHD1
        RTSAFE_C2P_RMHD1=RTSAFE_C2P_RMHD1-DX
        IF(TEMP.EQ.RTSAFE_C2P_RMHD1)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNC_C2P_RMHD1(RTSAFE_C2P_RMHD1,F,DF,gam,d,e,s2,b2,sb2,w)
     IF(F.LT.0.) THEN
        XL=RTSAFE_C2P_RMHD1
        FL=F
     ELSE
        XH=RTSAFE_C2P_RMHD1
        FH=F
     ENDIF
11   CONTINUE
     FAILED = .TRUE.
     RTSAFE_C2P_RMHD1 = 0.
     RETURN
END FUNCTION
   !
FUNCTION RTSAFE_C2P_RMHD2(X1,X2,XACC,g1, D, k2, B2, kB, E, H, FAILED)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=200
  INTEGER               :: J
  REAL                  :: RTSAFE_C2P_RMHD2
  REAL                  :: X1,X2,XACC, g1, D, k2, B2, kB, E, H
  REAL                  :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP
  LOGICAL               :: FAILED
  !
  FAILED = .FALSE.
  CALL FUNC_C2P_RMHD2(X1, FL, DF, g1, D, k2, B2, kB, E, H)
  IF(FL.EQ.0.) THEN
     RTSAFE_C2P_RMHD2=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RMHD2(X2, FH, DF, g1, D, k2, B2, kB, E, H)
  IF(FH.EQ.0.) THEN
     RTSAFE_C2P_RMHD2=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  RTSAFE_C2P_RMHD2=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RMHD2(RTSAFE_C2P_RMHD2, F, DF, g1, D, k2, B2, kB, E, H)
  DO 11 J=1,MAXIT
     IF(((RTSAFE_C2P_RMHD2-XH)*DF-F)*((RTSAFE_C2P_RMHD2-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        RTSAFE_C2P_RMHD2=XL+DX
        IF(XL.EQ.RTSAFE_C2P_RMHD2)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=RTSAFE_C2P_RMHD2
        RTSAFE_C2P_RMHD2=RTSAFE_C2P_RMHD2-DX
        IF(TEMP.EQ.RTSAFE_C2P_RMHD2)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNC_C2P_RMHD2(RTSAFE_C2P_RMHD2,F, DF, g1, D, k2, B2, kB, E, H)
     IF(F.LT.0.) THEN
        XL=RTSAFE_C2P_RMHD2
        FL=F
     ELSE
        XH=RTSAFE_C2P_RMHD2
        FH=F
     ENDIF
11   CONTINUE
     FAILED = .TRUE.
     RETURN
END FUNCTION
!    
FUNCTION RTSAFE_C2P_RHD1(X1,X2,XACC,d,e,s2,FAILED)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=200
  INTEGER               :: J
  REAL                  :: RTSAFE_C2P_RHD1
  REAL                  :: X1,X2,XACC,gam,d,e,s2
  REAL                  :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP
  LOGICAL               :: FAILED
  FAILED = .FALSE.
  CALL FUNC_C2P_RHD1(X1,FL,DF,d,e,s2,FAILED)
  IF (FAILED) RETURN
  IF(FL.EQ.0.) THEN
     RTSAFE_C2P_RHD1=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RHD1(X2,FH,DF,d,e,s2,FAILED)
  IF (FAILED) RETURN
  IF(FH.EQ.0.) THEN
     RTSAFE_C2P_RHD1=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  RTSAFE_C2P_RHD1=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RHD1(RTSAFE_C2P_RHD1,F,DF,d,e,s2,FAILED)
  IF (FAILED) RETURN
  DO 11 J=1,MAXIT
     IF(((RTSAFE_C2P_RHD1-XH)*DF-F)*((RTSAFE_C2P_RHD1-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        RTSAFE_C2P_RHD1=XL+DX
        IF(XL.EQ.RTSAFE_C2P_RHD1)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=RTSAFE_C2P_RHD1
        RTSAFE_C2P_RHD1=RTSAFE_C2P_RHD1-DX
        IF(TEMP.EQ.RTSAFE_C2P_RHD1)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNC_C2P_RHD1(RTSAFE_C2P_RHD1,F,DF,d,e,s2,FAILED)
     IF (FAILED) RETURN
     IF(F.LT.0.) THEN
        XL=RTSAFE_C2P_RHD1
        FL=F
     ELSE
        XH=RTSAFE_C2P_RHD1
        FH=F
     ENDIF
11   CONTINUE
     FAILED = .TRUE.
     write(*,*)'RTSAFE exceeding maximum iterations'
     stop
     RETURN
   END FUNCTION
   !
FUNCTION RTSAFE_C2P_RHD2(X1,X2,XACC,qi,ri,kappa,FAILED)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: MAXIT=200
  INTEGER               :: J
  REAL                  :: RTSAFE_C2P_RHD2
  REAL                  :: X1,X2,XACC,qi,ri,kappa
  REAL                  :: FL,FH,DF,XH,XL,SWAP,DXOLD,DX,F,TEMP
  LOGICAL               :: FAILED
  FAILED = .FALSE.
  CALL FUNC_C2P_RHD2(X1,FL,DF,qi,ri,kappa,FAILED)
  IF (FAILED) RETURN
  IF(FL.EQ.0.) THEN
     RTSAFE_C2P_RHD2=X1
     RETURN
  ENDIF
  CALL FUNC_C2P_RHD2(X2,FH,DF,qi,ri,kappa,FAILED)
  IF (FAILED) RETURN
  IF(FH.EQ.0.) THEN
     RTSAFE_C2P_RHD2=X2
     RETURN
  ENDIF
  IF(FL*FH.GT.0.) THEN
     FAILED = .TRUE.
     RETURN
  ENDIF
  IF(FL.LT.0.)THEN
     XL=X1
     XH=X2
  ELSE
     XH=X1
     XL=X2
     SWAP=FL
     FL=FH
     FH=SWAP
  ENDIF
  RTSAFE_C2P_RHD2=.5*(X1+X2)
  DXOLD=ABS(X2-X1)
  DX=DXOLD
  CALL FUNC_C2P_RHD2(RTSAFE_C2P_RHD2,F,DF,qi,ri,kappa,FAILED)
  IF (FAILED) RETURN
  DO 11 J=1,MAXIT
     IF(((RTSAFE_C2P_RHD2-XH)*DF-F)*((RTSAFE_C2P_RHD2-XL)*DF-F).GE.0. &
          .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
        DXOLD=DX
        DX=0.5*(XH-XL)
        RTSAFE_C2P_RHD2=XL+DX
        IF(XL.EQ.RTSAFE_C2P_RHD2)RETURN
     ELSE
        DXOLD=DX
        DX=F/DF
        TEMP=RTSAFE_C2P_RHD2
        RTSAFE_C2P_RHD2=RTSAFE_C2P_RHD2-DX
        IF(TEMP.EQ.RTSAFE_C2P_RHD2)RETURN
     ENDIF
     IF(ABS(DX).LT.XACC) RETURN
     CALL FUNC_C2P_RHD2(RTSAFE_C2P_RHD2,F,DF,qi,ri,kappa,FAILED)
     IF (FAILED) RETURN
     IF(F.LT.0.) THEN
        XL=RTSAFE_C2P_RHD2
        FL=F
     ELSE
        XH=RTSAFE_C2P_RHD2
        FH=F
     ENDIF
11   CONTINUE
     write(*,*)'RTSAFE exceeding maximum iterations'
     stop
     RETURN
   END FUNCTION
   !
   !--------------------------------------------------------
FUNCTION ZBRENT_C2P_RHD2 ( x1, x2, tol, qi, ri, kappa, FAILED)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: x1,x2,tol,qi,ri,kappa
  LOGICAL, INTENT(OUT):: FAILED
  REAL :: ZBRENT_C2P_RHD2, fo
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL    :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,df
  INTERFACE
     PURE SUBROUTINE FUNC_C2P_RHD2(z,f,df,qi,ri,kappa,FAILED)
     USE typesDef, ONLY: EQN 
     IMPLICIT NONE
     INTEGER    :: iter
     REAL       :: z,f,df
     REAL       :: qi,ri,kappa
     LOGICAL    :: FAILED
     REAL       :: gamma_mo, z2, w, epsilon, h, a
     INTENT(IN) :: z,qi,ri,kappa
     INTENT(OUT):: f,df,FAILED
     END SUBROUTINE FUNC_C2P_RHD2
  END INTERFACE
  FAILED = .FALSE.
  a=x1
  b=x2
  CALL FUNC_C2P_RHD2(a,fo,df,qi,ri,kappa,FAILED)
  fa = fo
  CALL FUNC_C2P_RHD2(b,fo,df,qi,ri,kappa,FAILED)
  fb = fo
  IF ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) THEN
    FAILED = .TRUE.
    RETURN
  ENDIF
  c=b
  fc=fb
  do iter=1,ITMAX
     if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        ZBRENT_C2P_RHD2=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     CALL FUNC_C2P_RHD2(b,fo,df,qi,ri,kappa,FAILED)
     fb = fo
  end do
  !call nrerror('zbrent_rout: exceeded maximum iterations')
  ZBRENT_C2P_RHD2=b
END FUNCTION ZBRENT_C2P_RHD2
   !
PURE SUBROUTINE FUNC_C2P_RMHD1(x,f,df,gam,d,e,s2,b2,sb2,w)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2007) A&A, 473, 11-30
  ! and it corresponds to their choice 3 in Section 3.2
  !
  IMPLICIT NONE
  REAL, PARAMETER :: third=1./3.
  INTEGER    :: iter
  REAL       :: x,f,df,v2,rho,c0,c2,c3,dw,dc2,dc3,dlogw,wb,vb2
  REAL       :: gam,d,e,s2,b2,sb2,w
  INTENT(IN) :: x,gam,d,e,s2,b2,sb2
  INTENT(OUT):: f,df,w
  
  v2=x
  rho=d*sqrt(1.-v2)
  
  c3=1.-gam*(1.-v2)
  c2=gam*rho+.5*b2*(1.+v2)-e
  c0=-0.5*sb2
  
  ! For every x=v2, we solve a cubic in W of the form: 
  ! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)
  ! W=y of the paper. If sb=0 ( meaning c0 = 0), 
  ! w = -c2/c3 > 0 and dw = 0 in the do loop below. 
  ! If -c2/c3 < 0 when sb=0, which makes w=0, 
  ! this is a signature that something was wrong before.
  
  if ( abs ( c0) < 1.0d-20) then
     w = -c2 / c3
  else
     w = max ( - c2 / c3, ( -c0 / c3)**third)
     do iter = 1,100
        dw = -((c3*w + c2)*w**2 + c0)/((3*c3*w + 2*c2)*w)
        if (abs(dw/w)<1.e-10) exit
        w = w + dw
     end do
  endif
  
  dc3   = gam
  dc2   = 0.5 * ( b2 - gam * rho / (1.0 - v2))
  dlogw = -( dc3 * w + dc2 ) / ( 3.0 * c3 * w + 2.0 * c2)
  wb    = w + b2
  vb2   = sb2 / w**2
  f     = wb**2 * v2 - ( 2.0 * w + b2) * vb2 - s2
  df    = wb * ( wb + 2.0 * dlogw * ( w * v2 + vb2))
  
END SUBROUTINE FUNC_C2P_RMHD1

PURE SUBROUTINE FUNC_C2P_RMHD2(x, f, df, g1, D, k2, B2, kB, E, H)
  !
  ! This is the CONS2PRIM strategy adopted by Del Zanna et al. (2003) A&A, 400, 397-413
  !                                    and by Dumbser et al. (2008) JCP, 227, 8209-8253
  !
  IMPLICIT NONE
  REAL       :: x, g1, D, k2, B2, kB, E
  REAL       :: f, df
  REAL       :: a1, b1, c1, d1, p1, rho, v2, T2
  REAL       :: dv2, sqrvc, c0, c2, c3, H, dH, drho
  COMPLEX    :: CH
  INTENT(IN) :: x, g1, D, k2, B2, kB, E
  INTENT(OUT):: f, df, H
  REAL, PARAMETER    :: sqr3 = 1.7320508075688772935274463415059
  
  v2    = x
  sqrvc = SQRT((1-v2))
  rho  = D*sqrvc
  T2   = B2*k2-kB*kB
  a1   = (1-(1-v2)/g1)
  b1   = (-E+rho/g1+0.5*B2+2*(1-(1-v2)/g1)*B2)
  c1   = (2*(-E+rho/g1+0.5*B2)*B2+(1-(1-v2)/g1)*B2**2)
  d1   = (-E+rho/g1+0.5*B2)*B2**2+0.5*T2
  !                    
  p1  = c1/a1 - b1*b1/a1/a1/3.
  CH   = 1/a1*(36*c1*b1*a1-108*d1*a1**2-8*b1**3+12*sqr3*csqrt(CMPLX(4*c1**3*a1-c1**2*b1**2-18*c1*b1*a1*d1+27*d1**2*a1**2+4*d1*b1**3))*a1)**(1.D0/3.D0)/6-2.D0/3.D0*(3*c1*a1-b1**2)/a1/(36*c1*b1*a1-108*d1*a1**2-8*b1**3+12*sqr3*csqrt(CMPLX(4*c1**3*a1-c1**2*b1**2-18*c1*b1*a1*d1+27*d1**2*a1**2+4*d1*b1**3))*a1)**(1.D0/3.D0)-b1/a1/3
  H    = REAL(CH)
  !
  f    = H**2*v2 + (2*H+B2)*T2/(H+B2)**2 - k2           

  drho = -0.5*D/SQRT((1-v2))
  dH   = (H**2+H*B2+drho*H+drho*B2)/(-3*H*g1-2*B2*g1-3*H*v2-v2*B2+3*H+B2+2*E*g1-2*rho)
  
  df  = H**2 + 2*H*dH*(v2 - T2/(H+B2)**3 )
  
END SUBROUTINE FUNC_C2P_RMHD2
    
SUBROUTINE FUNC_C2P_RMHD3(x,gam,sb2,b2,s2,d,e,eps,alpha,beta)

  implicit none
  REAL :: x(2),gam,sb2,b2,s2,d,e,eps
  REAL :: alpha(2,2),beta(2)
  REAL :: v2,w,w2,vv,rho,p,ssbb

  v2=x(1)
  w =x(2)

  w2=w**2
  vv=max(eps,1.-v2)
  rho=d*sqrt(vv)
  p=max(eps,gam*(w*vv-rho))

  ssbb=sb2/w2

  beta(1)=-((w+b2)**2*v2-(2.*w+b2)*ssbb-s2)
  beta(2)=-(w-p+.5*(1.+v2)*b2-.5*ssbb-e)

  alpha(1,2)=2.*(w+b2)*v2+2.*(1.+b2/w)*ssbb
  alpha(1,1)=(w+b2)**2
  alpha(2,2)=1.-gam*vv+ssbb/w
  alpha(2,1)=-gam*(-w+.5*rho/vv)+.5*b2

END SUBROUTINE FUNC_C2P_RMHD3


PURE SUBROUTINE FUNC_C2P_RHD1(x,f,df,d,e,s2,FAILED)
  !
  ! This is the CONS2PRIM strategy adopted in Whisky and originally described in Baiotti's PhD thesis. 
  ! See also Appendix D of the book "Relativistic hydrodynamics" by Rezzolla & Zanotti (2013)
  !
  USE typesDef, ONLY: EQN
  IMPLICIT NONE
  INTEGER    :: iter
  REAL       :: x,f,df
  REAL       :: d,e,s2
  LOGICAL    :: FAILED
  REAL       :: gamma_mo, tmp, Z, rho, epsilon, drhodp, depsilondp
  INTENT(IN) :: x,d,e,s2
  INTENT(OUT):: f,df,FAILED
  
  ! x = pressure, the unknown
  FAILED   = .FALSE.
  gamma_mo = EQN%gamma - 1.0
  Z       = e + x + d
  IF ((Z**2 - s2) .GT. 0.0) THEN
     tmp     = SQRT(Z**2 - s2)
  ELSE
     FAILED = .TRUE.
     RETURN
  ENDIF
  rho     = d / Z * tmp
  epsilon = ( tmp - x*Z/tmp - d) / d
  
  drhodp     = d*s2 / ( Z**2 * tmp) 
  depsilondp = x*s2 / ( rho * tmp**2 * Z)
  
  f     = x - rho*epsilon*gamma_mo
  df    = 1.0 - epsilon*gamma_mo*drhodp - rho*gamma_mo*depsilondp 
  
END SUBROUTINE FUNC_C2P_RHD1

!
PURE SUBROUTINE FUNC_C2P_RHD2(z,f,df,qi,ri,kappa,FAILED)
  !
  ! This is the CONS2PRIM strategy adopted by Galeazzi et al. (2013), Phys Rev D 88 064009, Appendix C
  ! See also the book "Relativistic hydrodynamics" by Rezzolla & Zanotti (2013).
  !
  USE typesDef, ONLY: EQN
  IMPLICIT NONE
  INTEGER    :: iter
  REAL       :: z,f,df
  REAL       :: qi,ri,kappa
  LOGICAL    :: FAILED
  REAL       :: z2, lf, epsilon, h
  INTENT(IN) :: z,qi,ri,kappa
  INTENT(OUT):: f,df,FAILED
  
  ! z = (Lorentz factor)*v , the unknown
  
  FAILED   = .FALSE.
  z2      = z**2  
  lf      = SQRT(1.0 + z2)
  epsilon = lf*(qi+1.0) - z*ri - 1.0
  h       = 1.0 + EQN%gamma*epsilon
  
  f  = z - ri/h
  df = 1.0 + ri/h**2 * EQN%gamma * (z/lf*(qi + 1.0) - ri)
  
END SUBROUTINE FUNC_C2P_RHD2
!
!    
FUNCTION ZBRENT_LF_POLY(x1,x2,tol,gamma,kpol,dd,b22,ss,sb)
  ! This is used in the isentropic case
  !
  IMPLICIT NONE
  REAL             :: x1,x2,tol,dd,b22,ss,sb
  REAL             :: zbrent_lf_poly,fo,gamma,kpol
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER  :: EPS=epsilon(x1)
  INTEGER          :: iter
  REAL             :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  REAL             :: gamma_comb,hh,gamma_mo
  
  INTENT(IN) :: x1,x2,tol,gamma,kpol,dd,b22,ss,sb
  
  a=x1
  b=x2
  gamma_mo   = gamma - 1.0
  gamma_comb = kpol*gamma/gamma_mo
  
  !  hh is just to compute fa, the value of the function at the left 
  !  bracketting. Note that hh = Z/D at W=1
  hh = 1.0d0 + gamma_comb*dd**(gamma - 1.0d0)
  
  fa = ss + (sb/dd/hh)**2*(2.0d0*dd*hh + b22)
  
111 fb = (ss + sb*sb/(dd*(1.0d0 + gamma_comb*(dd/b)**(gamma - 1.0d0))*        &
  & b)**2*(2.0d0*dd*(1.0d0 + gamma_comb*(dd/b)**(gamma - 1.0d0))*      &
  & b + b22))*b**2 - ((1.0d0 + gamma_comb*(dd/b)**(gamma - 1.0d0))*    &
  & dd*b + b22)**2*(b**2 - 1.0d0)
  
  if ( fa*fb > 0.0 ) then
     ! update the bracketing interval: increase b !
     b = b*10.0
     !     write(*,*)'from zbrent_lf_poly',b
     goto 111
  endif
  c=b
  fc=fb
  do iter=1,ITMAX
     if ( fb*fc > 0.0 ) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_lf_poly=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm 
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     
     fb = (ss + sb*sb/(dd*(1.0d0 + gamma_comb*(dd/b)**(gamma - 1.0d0))*       &
     & b)**2*(2.0d0*dd*(1.0d0 + gamma_comb*(dd/b)**(gamma - 1.0d0))* &
     & b + b22))*b**2 - ((1.0d0 + gamma_comb*(dd/b)**                &
     & (gamma - 1.0d0))* dd*b + b22)**2*(b**2 - 1.0d0)
     
     
  end do
  zbrent_lf_poly=b
END FUNCTION ZBRENT_LF_POLY

SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      !
      IMPLICIT NONE
      INTEGER         :: NP, N
      INTEGER, PARAMETER :: NMAX=100
      INTEGER            :: INDX(N)
      REAL,    PARAMETER :: TINY=1.0E-20
      REAL               :: A(NP,NP), VV(NMAX)
      ! Local Variables
      INTEGER            :: I, J, K, IMAX
      REAL               :: D, AAMAX, SUM, DUM
      !
      D=1.
      DO I=1,N
        AAMAX=0.
        DO J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        ENDDO
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
      ENDDO  
      DO J=1,N
        IF (J.GT.1) THEN
          DO I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
              ENDDO
              A(I,J)=SUM
            ENDIF
          ENDDO
        ENDIF
        AAMAX=0.
        DO I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        ENDDO
        IF (J.NE.IMAX)THEN
          DO K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          ENDDO
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO I=J+1,N
            A(I,J)=A(I,J)*DUM
          ENDDO
        ENDIF
      ENDDO
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
END SUBROUTINE LUDCMP

! ***************************************************************************

SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      !
      IMPLICIT NONE
      INTEGER   :: NP, N
      INTEGER   :: INDX(N)
      REAL      :: A(NP,NP), B(N)
      ! Local Variables
      INTEGER   :: II, I, J, LL
      REAL      :: SUM
      !
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO J=II,I-1
            SUM=SUM-A(I,J)*B(J)
          ENDDO
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
      ENDDO
      DO I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
          ENDDO
        ENDIF
        B(I)=SUM/A(I,I)
      ENDDO
      RETURN
END SUBROUTINE

SUBROUTINE exact_quartic(solution,coefficients)
  ! ---------------------------------------------------
  IMPLICIT NONE 
  ! ---------------------------------------------------
  REAL    :: coefficients(5)
  COMPLEX :: solution(4) 
  ! ---------------------------------------------------
  ! Variables associated to the solution of the quartic
  ! ---------------------------------------------------
  REAL smallroot
  REAL, PARAMETER  :: smallnum = 1.0e-8
  COMPLEX ec1, ec2
  REAL  temparr_01, temparr_02, temparr_03, temparr_04
  REAL  temparr_05, temparr_06, temparr_07, temparr_08
  REAL  temparr_09, eta
  REAL ::  cubic_a, cubic_b, cubic_c, cubic_d, cubic_p, cubic_q
  REAL :: cubic_disc
  REAL :: quartic_a, quartic_b, quartic_c, quartic_d
  REAL :: quartic_e, quartic_p, quartic_q, quartic_r
  COMPLEX :: cubic_u, cubic_v, cubic_y1, cubic_y2, cubic_y3
  COMPLEX :: ctemparr_01, ctemparr_02, ctemparr_03, ctemparr_04 
  COMPLEX :: quartic_y1, quartic_y2, quartic_y3, quartic_y4
  ! ---------------------------------------------------
  !  Compute the coefficients of the quartic
  ! ---------------------------------------------------
  ! Coefficient of 4-th degree
  quartic_a  = coefficients(1) 
  ! Coefficient of 3-th degree
  quartic_b  = coefficients(2) 
  ! Coefficient of 2-th degree
  quartic_c  = coefficients(3) 
  ! Coefficient of 1-th degree
  quartic_d  = coefficients(4) 
  ! Coefficient of 0-th degree
  quartic_e  = coefficients(5) 
  !     *****************************************************************
  !     * Beginning of solving the quartic.
  !     *****************************************************************
  ec1 = CMPLX ( -0.5, 0.5 * SQRT ( 3.0) )
  ec2 = CONJG ( ec1)
  smallroot = 1.0e-30
  !     * Load the variables and transform the quartic into the corresponding
  !     * reducing cubic.
  quartic_p  = - 3.0 * quartic_b**2 / ( 8.0 * quartic_a**2)    &
  &    + quartic_c  / quartic_a 
  quartic_q  = quartic_b**3 / ( 8.0 * quartic_a**3)            & 
  &     - quartic_b  * quartic_c / ( 2.0 * quartic_a**2 ) &
  &     + quartic_d  / quartic_a 
  quartic_r  = - 3.0 * quartic_b**4                            &
  &     / ( 256.0 * quartic_a**4 )                        &
  &     + quartic_b**2 * quartic_c                        &
  &     / ( 16.0 * quartic_a**3 )                         &
  &     - quartic_b  * quartic_d                          &
  &     / ( 4.0 * quartic_a**2 )                          &
  &     + quartic_e  / quartic_a 
  cubic_a  = 1.0
  cubic_b  = 0.5 * quartic_p 
  cubic_c  = 0.0625 * quartic_p**2 - 0.25 * quartic_r 
  cubic_d  = - 0.015625 * quartic_q**2
  !
  !     *****************************************************************
  !     * Beginning of solving the cubic.
  !     *****************************************************************
  !
  cubic_p  = ( 3.0 * cubic_a  * cubic_c  - cubic_b**2 )     &
  &     / ( 9.0 * cubic_a**2)
  
  cubic_q  = 2.0 * cubic_b**3                               &
  &     / ( 27.0 * cubic_a**3) - cubic_b  * cubic_c    &
  &     / ( 3.0 * cubic_a**2) + cubic_d  / cubic_a 
  !
  cubic_q  = 0.5 * cubic_q 
  !
  cubic_disc  = - cubic_p**3 - cubic_q**2
  !
  IF ( cubic_disc  .GE. 0.0 ) THEN
     !        * Roots are real
     temparr_01  = SQRT ( - cubic_p**3 )
     temparr_02  = - cubic_q / AMAX1 ( smallroot, temparr_01  )
     temparr_02  = AMAX1 ( temparr_02 , -1.0)
     temparr_02  = AMIN1 ( temparr_02 , 1.0)
     temparr_02  = ACOS ( temparr_02  ) / 3.0
     temparr_01  = AMAX1 ( smallroot, temparr_01  )**0.33333333333333
     !
     cubic_u  = CMPLX ( temparr_01 * COS ( temparr_02  ),     &
     &             temparr_01 * SIN ( temparr_02  ) )
     cubic_v  = CONJG ( cubic_u  )
  ELSE
     !        * Roots are complex
     temparr_01  = SQRT ( - cubic_disc )
     temparr_02  = temparr_01 
     temparr_01  = - cubic_q  + temparr_01 
     temparr_02  = - cubic_q  - temparr_02 
     !    
     cubic_u  = SIGN ( 1.0, temparr_01  )
     cubic_u  = cubic_u * AMAX1 ( smallroot,                    &
     &     ABS ( temparr_01  ) )**0.33333333333333
     !
     cubic_v  = SIGN ( 1.0, temparr_02  )
     cubic_v  = cubic_v * AMAX1 ( smallroot,                    &
     &     ABS ( temparr_02  ) )**0.33333333333333
  END IF
  !
  !     * Note that "cubic_y1" will always carry a real root.
  !
  cubic_y1  = cubic_u  + cubic_v 
  cubic_y2  = ec1 * cubic_u  + ec2 * cubic_v 
  cubic_y3  = ec2 * cubic_u  + ec1 * cubic_v 
  temparr_01  = - cubic_b  / ( 3.0 * cubic_a  )
  ctemparr_01  = CMPLX ( temparr_01 , 0.0)
  cubic_y1  = cubic_y1  + ctemparr_01 
  cubic_y2  = cubic_y2  + ctemparr_01 
  cubic_y3  = cubic_y3  + ctemparr_01 
  !     *****************************************************************
  !     * End of solving the cubic.
  !     *****************************************************************
  !     * The (complex) square roots of the roots of the reducing cubic need
  !     * to be chosen so that they satisfy a consistency condition. Here we
  !     * find that specific choice that satisfies the consistency constraint.
  cubic_y1  = CSQRT ( cubic_y1  )
  cubic_y2  = CSQRT ( cubic_y2  )
  cubic_y3  = CSQRT ( cubic_y3  )
  !
  ctemparr_01  = - cubic_y1 
  ctemparr_02  = - cubic_y2 
  ctemparr_03  = - cubic_y3 
  !
  ctemparr_04  = CMPLX ( - 0.125 * quartic_q , 0.0)
  !    
  IF ( CABS ( ctemparr_01  * cubic_y2 * cubic_y3       &
     &        - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y1  = ctemparr_01 
  ELSE IF ( CABS ( cubic_y1  * ctemparr_02 * cubic_y3     &
     &       - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y2  = ctemparr_02 
  ELSE IF ( CABS ( cubic_y1  * cubic_y2 * ctemparr_03     &
     &       - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y3  = ctemparr_03 
  ELSE IF ( CABS ( ctemparr_01  * ctemparr_02 * cubic_y3  &
     &        - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y1  = ctemparr_01 
     cubic_y2  = ctemparr_02 
  ELSE IF ( CABS ( cubic_y1  * ctemparr_02 * ctemparr_03   &
     &       - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y2  = ctemparr_02 
     cubic_y3  = ctemparr_03 
  ELSE IF ( CABS ( ctemparr_01  * cubic_y2 * ctemparr_03   &
     &        - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y1  = ctemparr_01 
     cubic_y3  = ctemparr_03 
  ELSE IF ( CABS ( ctemparr_01  * ctemparr_02 * ctemparr_03 &
     &       - ctemparr_04  ) .LE. smallnum) THEN
     cubic_y1  = ctemparr_01 
     cubic_y2  = ctemparr_02 
     cubic_y3  = ctemparr_03 
  END IF
  !
  temparr_01  = - quartic_b   / ( 4.0 * quartic_a  )
  ctemparr_01  = CMPLX ( temparr_01 , 0.0)
  quartic_y1  = cubic_y1  + cubic_y2 + cubic_y3  + ctemparr_01 
  quartic_y2  = cubic_y1  - cubic_y2 - cubic_y3  + ctemparr_01 
  quartic_y3  = - cubic_y1  + cubic_y2 - cubic_y3  + ctemparr_01 
  quartic_y4  = - cubic_y1  - cubic_y2 + cubic_y3  + ctemparr_01 
  !     * Now sort the four roots in increasing order.
  IF ( REAL ( quartic_y1 ) .GT. REAL ( quartic_y2 ) ) THEN
     ctemparr_01  = quartic_y1 
     quartic_y1  = quartic_y2 
     quartic_y2  = ctemparr_01 
  END IF
  !
  IF ( REAL ( quartic_y2 ) .GT. REAL ( quartic_y3 ) ) THEN
     ctemparr_01  = quartic_y2 
     quartic_y2  = quartic_y3 
     quartic_y3  = ctemparr_01 
  END IF
  !
  IF ( REAL ( quartic_y3 ) .GT. REAL ( quartic_y4 ) ) THEN
     ctemparr_01  = quartic_y3 
     quartic_y3  = quartic_y4 
     quartic_y4  = ctemparr_01 
  END IF
  !
  IF ( REAL ( quartic_y1 ) .GT. REAL ( quartic_y2 ) ) THEN
     ctemparr_01  = quartic_y1 
     quartic_y1  = quartic_y2 
     quartic_y2  = ctemparr_01 
  END IF
  !
  IF ( REAL ( quartic_y2 ) .GT. REAL ( quartic_y3 ) ) THEN
     ctemparr_01  = quartic_y2 
     quartic_y2  = quartic_y3 
     quartic_y3  = ctemparr_01 
  END IF
  !
  IF ( REAL ( quartic_y1 ) .GT. REAL ( quartic_y2 ) ) THEN
     ctemparr_01  = quartic_y1 
     quartic_y1  = quartic_y2 
     quartic_y2  = ctemparr_01 
  END IF
  !
  solution(1) = quartic_y1 
  solution(2) = quartic_y2 
  solution(3) = quartic_y3 
  solution(4) = quartic_y4 
  !
  
END SUBROUTINE exact_quartic

    
    
    



SUBROUTINE Excision_uh 
    USE typesDef
    IMPLICIT NONE
    INTEGER :: iElem, i, j, k 
    REAL    :: x0(d), xGP(d) 
    REAL    :: u0(nVar), par0(nParam)      
    !
    DO iElem = 1, nElem 
        x0 = x(:,tri(1,iElem)) + 0.5*dx   ! get the coordinate of the element barycenter 
        !IF( SQRT(SUM(x0**2))<ExcisionRadius ) THEN
        IF( MAXVAL(ABS(x0)).LT.ExcisionRadius ) THEN
          DO k = 1, nDOF(3) 
           DO j = 1, nDOF(2) 
            DO i = 1, nDOF(1) 
              xGP = x(:,tri(1,iElem)) + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:) 
              !IF( SQRT(SUM(xGP**2))>ExcisionRadius ) CYCLE                     
              CALL InitialField(u0,par0,xGP,0.0) 
              uh(:,i,j,k,iElem) = u0 
              parh(:,i,j,k,iElem) = par0  
            ENDDO
           ENDDO
          ENDDO
        ELSE
            CONTINUE 
        ENDIF 
    ENDDO 
    ! 
END SUBROUTINE Excision_uh
    
    
SUBROUTINE Excision_qh 
    USE typesDef
    IMPLICIT NONE
    REAL    :: x0(d), xGP(d) 
    REAL    :: u0(nVar), par0(nParam)      
    ! Local variables
    INTEGER :: iElem,i,j,k,l,iVar,iDim, iter
    REAL    :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side
    REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side
    REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
    REAL    :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
    REAL    :: lSh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! nonlinear source vector in each space-time DOF
    REAL    :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))                     ! time-averaged nonlinear flux tensor in each space DOF
    REAL    :: lShi(nVar,nDOF(1),nDOF(2),nDOF(3))                       ! time-averaged nonlinear source vector in each space DOF
    REAL    :: lqbnd(nVar,nDOF(2),nDOF(3),6)                            ! time-averaged space-time degrees of freedom
    REAL    :: lFbnd(nVar,nDOF(2),nDOF(3),6)                            ! time-averaged nonlinear flux tensor in each space-time DOF
    !
    DO iElem = 1, nElem 
        x0 = x(:,tri(1,iElem)) + 0.5*dx   ! get the coordinate of the element barycenter 
        !IF( SQRT(SUM(x0**2))<ExcisionRadius ) THEN
        IF( MAXVAL(ABS(x0)).LT.ExcisionRadius ) THEN
         DO l = 1, nDOF(0)  
          DO k = 1, nDOF(3) 
           DO j = 1, nDOF(2) 
            DO i = 1, nDOF(1) 
              xGP = x(:,tri(1,iElem)) + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:) 
              CALL InitialField(u0,par0,xGP,0.0)               
              lqh(:,i,j,k,l) = u0 
              CALL PDEFlux(lFh(:,:,i,j,k,l),u0,parh(:,i,j,k,iElem))
            ENDDO
           ENDDO
          ENDDO
         ENDDO  
        !          
        !
        ! Immediately compute the time-averaged space-time polynomials
        !
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    qhi(:,i,j,k,iElem) = MATMUL( lqh(:,i,j,k,:), wGPN )
                    DO iDim = 1, nDim
                        Fhi(:,iDim,i,j,k,iElem) = MATMUL( lFh(:,iDim,i,j,k,:), wGPN )
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        !
        ! Compute the bounday-extrapolated values for Q and F*n
        !
        ! x-direction: face 1 (left) and face 2 (right)
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                Qbnd(:,j,k,1,iElem) = MATMUL( qhi(:,:,j,k,iElem),   FLCoeff )   ! left
                Qbnd(:,j,k,2,iElem) = MATMUL( qhi(:,:,j,k,iElem),   FRCoeff )   ! right
                Fbnd(:,j,k,1,iElem) = MATMUL( Fhi(:,1,:,j,k,iElem), FLCoeff )   ! left
                Fbnd(:,j,k,2,iElem) = MATMUL( Fhi(:,1,:,j,k,iElem), FRCoeff )   ! right
            ENDDO
        ENDDO
        ! y-direction: face 3 (left) and face 4 (right)
        IF(nDim>=2) THEN
            DO k = 1, nDOF(3)
                DO i = 1, nDOF(1)
                    Qbnd(:,i,k,3,iElem) = MATMUL( qhi(:,i,:,k,iElem),   FLCoeff )   ! left
                    Qbnd(:,i,k,4,iElem) = MATMUL( qhi(:,i,:,k,iElem),   FRCoeff )   ! right
                    Fbnd(:,i,k,3,iElem) = MATMUL( Fhi(:,2,i,:,k,iElem), FLCoeff )   ! left
                    Fbnd(:,i,k,4,iElem) = MATMUL( Fhi(:,2,i,:,k,iElem), FRCoeff )   ! right
                ENDDO
            ENDDO
        ENDIF
        ! z-direction: face 5 (left) and face 6 (right)
        IF(nDim>=3) THEN
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    Qbnd(:,i,j,5,iElem) = MATMUL( qhi(:,i,j,:,iElem),   FLCoeff )   ! left
                    Qbnd(:,i,j,6,iElem) = MATMUL( qhi(:,i,j,:,iElem),   FRCoeff )   ! right
                    Fbnd(:,i,j,5,iElem) = MATMUL( Fhi(:,3,i,j,:,iElem), FLCoeff )   ! left
                    Fbnd(:,i,j,6,iElem) = MATMUL( Fhi(:,3,i,j,:,iElem), FRCoeff )   ! right
                ENDDO
            ENDDO
        ENDIF
        !                  
        ELSE
            CONTINUE 
        ENDIF 
    ENDDO 
    ! 
END SUBROUTINE Excision_qh
    
        
    