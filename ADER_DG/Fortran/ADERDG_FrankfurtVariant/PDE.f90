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
    REAL :: p, irho, lam, mu 
    !
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
#ifdef ELASTICITY
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
    
    
#endif 
    !            
END SUBROUTINE PDEFlux
!
! Nonconservative part of the PDE ( B(Q) * gradQ ) 
!    
SUBROUTINE PDENCP(BgradQ,Q,gradQ,par)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d), par(nParam)  
    REAL, INTENT(OUT) :: BgradQ(nVar,d) 
    ! Local variables 
    REAL :: p, irho, lam, mu 
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    !
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
    lam  = par(1)   
    mu   = par(2) 
    irho = 1./par(3) 
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2)
    Qz = gradQ(:,3)
    !
    BgradQ(1,1) = - (lam+2*mu)*Qx(7) 
    BgradQ(2,1) = - lam*Qx(7) 
    BgradQ(3,1) = - lam*Qx(7) 
    BgradQ(4,1) = - mu *Qx(8) 
    BgradQ(5,1) = 0. 
    BgradQ(6,1) = - mu *Qx(9) 
    BgradQ(7,1) = - irho *Qx(1) 
    BgradQ(8,1) = - irho *Qx(4) 
    BgradQ(9,1) = - irho *Qx(6) 
    
    BgradQ(1,2) = - lam*Qy(8) 
    BgradQ(2,2) = - (lam+2*mu)*Qy(8) 
    BgradQ(3,2) = - lam*Qy(8) 
    BgradQ(4,2) = - mu *Qy(7) 
    BgradQ(5,2) = - mu *Qy(9) 
    BgradQ(6,2) = 0. 
    BgradQ(7,2) = - irho *Qy(4) 
    BgradQ(8,2) = - irho *Qy(2) 
    BgradQ(9,2) = - irho *Qy(5) 
    
    BgradQ(1,3) = - lam*Qz(9) 
    BgradQ(2,3) = - lam*Qz(9) 
    BgradQ(3,3) = - (lam+2*mu)*Qz(9) 
    BgradQ(4,3) = 0. 
    BgradQ(5,3) = - mu *Qz(8) 
    BgradQ(6,3) = - mu *Qz(7) 
    BgradQ(7,3) = - irho *Qz(4) 
    BgradQ(8,3) = - irho *Qz(2) 
    BgradQ(9,3) = - irho *Qz(5) 

!    BgradQ = 0. 
!    BgradQ(1,1) = Qx(1) 
!    BgradQ(1,2) = Qy(1) 
    
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
    REAL :: p, irho, lam, mu 
    REAL :: B1(nVar,nVar), B2(nVar,nVar), B3(nVar,nVar)  
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
    lam  = par(1)   
    mu   = par(2) 
    irho = 1./par(3) 
    !
    B1 = 0. 
    B2 = 0. 
    B3 = 0. 
    !
    B1(1,7) = - (lam+2*mu)
    B1(2,7) = - lam 
    B1(3,7) = - lam
    B1(4,8) = - mu    
    B1(6,9) = - mu  
    B1(7,1) = - irho 
    B1(8,4) = - irho  
    B1(9,6) = - irho  
    
    B2(1,8) = - lam 
    B2(2,8) = - (lam+2*mu)
    B2(3,8) = - lam
    B2(4,7) = - mu 
    B2(5,9) = - mu 
    B2(7,4) = - irho 
    B2(8,2) = - irho 
    B2(9,5) = - irho 
    
    B3(1,9) = - lam
    B3(2,9) = - lam
    B3(3,9) = - (lam+2*mu)
    B3(5,8) = - mu 
    B3(6,7) = - mu 
    B3(7,4) = - irho
    B3(8,2) = - irho
    B3(9,5) = - irho
    
!    B1 = 0. 
!    B2 = 0. 
!    B3 = 0. 
!    B1(1,1) = 1.0 
!    B2(1,1) = 1.0 

#endif 
    !            
    Bn = B1*nv(1) + B2*nv(2) + B3*nv(3) 
    !
END SUBROUTINE PDEMatrixB 
    

SUBROUTINE PDEEigenvalues(Lambda,Q,par,nv)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), par(nParam), nv(d) 
    REAL, INTENT(OUT) :: Lambda(nVar) 
    ! Local variables 
    REAL :: p, u, c, cp, cs, rho0, lam, mu  
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
    ! get the local Lam? constants from the local material parameters 
    lam  = par(1)   
    mu   = par(2) 
    rho0 = par(3) 
    cp = SQRT((lam+2*mu)/rho0) 
    cs = SQRT(mu/rho0) 
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
END SUBROUTINE PDEEigenvalues

SUBROUTINE PDECons2Prim(V,Q,iErr)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)     :: Q(nVar)     ! vector of conserved quantities 
    REAL, INTENT(OUT)    :: V(nVar)     ! primitive variables 
    INTEGER, INTENT(OUT) :: iErr        ! error flag 
    ! Local variables 
    REAL                 :: p 
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
END SUBROUTINE PDECons2Prim 
        
SUBROUTINE PDEPrim2Cons(Q,V)
    USE typesDef, ONLY : nVar, nParam, d, EQN 
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)      :: V(nVar)     ! primitive variables 
    REAL, INTENT(OUT)     :: Q(nVar)     ! vector of conserved quantities 
    ! Local variables 
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
    Name(1) = 'sxx' 
    Name(2) = 'syy' 
    Name(3) = 'szz'
    Name(4) = 'sxy' 
    Name(5) = 'syz' 
    Name(6) = 'sxz' 
    Name(7) = 'u' 
    Name(8) = 'v' 
    Name(9) = 'w' 
#endif 
    !
END SUBROUTINE PDEVarName
    
SUBROUTINE PDEParName(Name) 
    USE typesDef, ONLY : nParam, d, EQN 
    IMPLICIT NONE       
    ! Argument list 
    CHARACTER(LEN=10), INTENT(OUT) :: Name(nParam)
    !
#ifdef EULER 
    !
#endif 
    !
#ifdef ELASTICITY 
    Name(1) = 'lambda' 
    Name(2) = 'mu' 
    Name(3) = 'rho'
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
END SUBROUTINE PDEAssurePositivity    
    
SUBROUTINE PointSource(sigma,t,waveform)
    USE typesDef 
    IMPLICIT NONE
    ! Argument declaration
    INTEGER, INTENT(IN) :: waveform  
    REAL, INTENT(IN)  :: t
    REAL, INTENT(OUT) :: sigma(nVar)  
    ! Local variables 
    REAL :: Pi 
    Pi = ACOS(-1.0) 
    SELECT CASE(waveform)
    CASE(1)
        sigma = 0. 
#ifdef ELASTICITY 
        !sigma(7) = t*t !
        sigma(7) = SIN(2*Pi*t/0.5) 
#endif 
    END SELECT 
END SUBROUTINE     