SUBROUTINE ADERRiemannSolver(lQbndL,lFbndL,lparbndL,lQbndR,lFbndR,lparbndR,nv,bcflag)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)     :: lQbndL(nVar,nDOF(2),nDOF(3))              ! left space-time degrees of freedom 
    REAL, INTENT(IN)     :: lQbndR(nVar,nDOF(2),nDOF(3))              ! right space-time degrees of freedom 
    REAL, INTENT(IN)     :: lparbndL(nParam,nDOF(2),nDOF(3))          ! left material parameters 
    REAL, INTENT(IN)     :: lparbndR(nParam,nDOF(2),nDOF(3))          ! right material parameters 
    REAL, INTENT(IN)     :: nv(d)                                     ! normal vector 
    INTEGER, INTENT(IN)  :: bcflag 
    REAL, INTENT(INOUT)  :: lFbndL(nVar,nDOF(2),nDOF(3))              ! left flux 
    REAL, INTENT(INOUT)  :: lFbndR(nVar,nDOF(2),nDOF(3))              ! right flux 
    ! Local variables 
    INTEGER           :: i,j,k,l,NormalNonZero, ml(1)  
    REAL              :: aux(d), QavL(nVar), QavR(nVar), Id(nVar,nVar), smax, Qav(nVar), parav(nParam)   
    REAL              :: parL(nParam), parR(nParam), xGP, yGP, xdeb, ydeb, sL, sR 
    REAL              :: LL(nVar), LR(nVar), LM(nVar), Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar)
    REAL              :: gradQ(nVar,d), src(nVar), absA(nVar,nVar)        
    !
    ml = MAXLOC(ABS(nv)) 
    NormalNonZero = ml(1) 
    ncp = 0.0 
    !
    ! Compute the average states from the left and the right, which we need to compute the numerical dissipation 
    QavL = 0. 
    QavR = 0. 
    DO k = 1, nDOF(3)
      DO j = 1, nDOF(2)
            aux = (/ 1., wGPN(j), wGPN(k) /) 
            QavL = QavL + PRODUCT(aux(1:nDim))*lQbndL(:,j,k) 
            QavR = QavR + PRODUCT(aux(1:nDim))*lQbndR(:,j,k) 
        ENDDO
    ENDDO
    Qav   = 0.5*(QavL+QavR) 
    !
    ! If necessary, compute the left and right averages of the material parameters 
    ! 
    IF(nParam.GT.0) THEN
        parL = 0. 
        parR = 0. 
        DO k = 1, nDOF(3)
          DO j = 1, nDOF(2)
                aux = (/ 1., wGPN(j), wGPN(k) /) 
                parL = parL + PRODUCT(aux(1:nDim))*lparbndL(:,j,k) 
                parR = parR + PRODUCT(aux(1:nDim))*lparbndR(:,j,k) 
            ENDDO
        ENDDO        
        parav = 0.5*(parL+parR) 
    ENDIF    
    !
    ! We now compute the numerical flux. Note that the scheme is at the moment written in 
    ! CONSERVATION FORM => no fluctuations, but real fluxes. 
    ! Later, this will be converted into the left and right fluctuations. 
    !
#ifdef LINEAR  
    ! 
    CALL PDEMatrixB(Bn,0.5*(QavL+QavR),nv,0.5*(parL+parR))   ! evaluate the system matrix just once in the averaged state 
    CALL DissipationMatrix(DM,QavL,QavR,parL,parR,nv)        ! according to the RIemann solver, compute the dissipation matrix 
    DO k = 1, nDOF(3)
      DO j = 1, nDOF(2)
            !CALL PDEMatrixB(Bn,0.5*(lQbndL(:,j,k)+lQbndR(:,j,k)),nv,0.5*(lparbndL(:,j,k)+lparbndR(:,j,k))) ! evaluate the system matrix in each Gaussian point again (slow) 
            lFbndL(:,j,k) = 0.5*MATMUL( Bn - DM, lQbndR(:,j,k) - lQbndL(:,j,k) ) 
            lFbndR(:,j,k) = 0.5*MATMUL( Bn + DM, lQbndR(:,j,k) - lQbndL(:,j,k) ) 
        ENDDO
    ENDDO
    !
#else
    !
    ! Here, we implement a very simple HLL scheme with scalar dissipation 
    ! We can change this into a more sophisticated Osher or HLLEM Riemann solver whenever needed! 
    !
    CALL PDEEigenvalues(LL,QavL,parL,nv) 
    CALL PDEEigenvalues(LR,QavR,parR,nv) 
    CALL PDEEigenvalues(LM,Qav,parav,nv) 
    !smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
    sL = MIN(0.0, MINVAL(LL), MINVAL(LM) ) 
    sR = MAX(0.0, MAXVAL(LR), MAXVAL(LM) ) 
    ! 
    !sL = -smax
    !sR = +smax 
    !    
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    !
    !Qav = 0.5*(QavL+QavR) 
    !CALL PDEMatrixB(Bn,Qav,nv,0.5*(parL+parR)) 
    !
    gradQ = 0.0  
    !
    IF(bcflag==0) THEN
        DO k = 1, nDOF(3)
          DO j = 1, nDOF(2)
#ifdef NONCONSERVATIVE 
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ,parav)
#endif 
                ! 
#ifdef NOFLUX
                lFbndL(:,j,k) = sL*sR/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )                                                     ! only numerical dissipation, no flux 
#else
                lFbndL(:,j,k) = ( sR*lFbndL(:,j,k) - sL*lFbndR(:,j,k) )/(sR-sL) + sR*sL/(sR-sL)*( lQbndR(:,j,k) - lQbndL(:,j,k) )   ! purely conservative flux 
#endif
                lFbndR(:,j,k) = lFbndL(:,j,k) - sR/(sR-sL)*ncp(:)                                                                   ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) - sL/(sR-sL)*ncp(:)                                                                   ! add the non-conservative product 
            ENDDO
        ENDDO
    ELSE
        IF(bcflag==1) THEN
            Qav = QavL 
        ELSE
            Qav = QavR
        ENDIF        
        ! Qav = 0.5*(QavL+QavR)  
        CALL PDEAbsA(absA,Qav,nv,parav) 
        !absA = 0.0
        !DO i = 1, nVar
        !    absA(i,i) = smax 
        !ENDDO        
        !CALL PDEMatrixB(Bn,Qav,nv,parav)      
        ! 
        DO k = 1, nDOF(3)
          DO j = 1, nDOF(2)
                !IF(bcflag==1) THEN
                !    Qav = lQbndL(:,j,k) 
                !ELSE
                !    Qav = lQbndR(:,j,k)  
                !ENDIF        
                !CALL PDEAbsA(absA,Qav,nv,parav)
                !ncp = MATMUL( Bn, lQbndR(:,j,k) - lQbndL(:,j,k) ) 
                Qav = 0.5*( lQbndR(:,j,k) + lQbndL(:,j,k) ) 
                gradQ(:,NormalNonZero) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ,parav)
                lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*MATMUL( absA, lQbndR(:,j,k) - lQbndL(:,j,k) )       ! purely conservative flux 
                lFbndR(:,j,k) = lFbndL(:,j,k) - 0.5*ncp(:)                                                                      ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) + 0.5*ncp(:)                                                                      ! add the non-conservative product 
            ENDDO
        ENDDO        
    ENDIF    
    !
#endif 
    !
#ifdef CCZ4EINSTEIN
    ! In the CCZ4 system do not add any numerical dissipation to the ODEs of the 4D metric 
    lFbndL(1:6,:,:)   = 0.0 
    lFbndL(17:20,:,:) = 0.0 
    lFbndL(55,:,:)    = 0.0 
    lFbndR(1:6,:,:)   = 0.0 
    lFbndR(17:20,:,:) = 0.0 
    lFbndR(55,:,:)    = 0.0 
#endif 
    !
END SUBROUTINE ADERRiemannSolver 
    
SUBROUTINE DissipationMatrix(DM,QL,QR,parL,parR,nv)
    USE typesDef, ONLY : d, nVar, nParam 
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)     :: QL(nVar), QR(nVar)                        ! left state, right state 
    REAL, INTENT(IN)     :: parL(nParam), parR(nParam)                ! left and right material parameters 
    REAL, INTENT(IN)     :: nv(d)                                     ! normal vector  
    REAL, INTENT(INOUT)  :: DM(nVar, nVar)                            ! Dissipation Matrix >= 0 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d), QavL(nVar), QavR(nVar), Id(nVar,nVar), smax, tv(d), sv(d), parav(nParam)  
    REAL              :: LL(nVar), LR(nVar), iRM(nVar,nVar), RM(nVar,nVar), TM(nVar,nVar), iTM(nVar,nVar), test(nVar,nVar), LM(nVar,nVar), lambda, mu, cp, cs, rho   
    REAL              :: lam, irho, alpha, ialpha, Qrot(nVar), s11, s12, s13, Q(nVar), uv(3), R(nVar,nVar), iR(nVar,nVar)     
    !
#if defined(GODUNOV) || defined(OSHER)      
    !
    IF(ABS(nv(1)-1.0) < 1e-14) THEN
        sv = (/ 0., 1., 0. /)
        tv = (/ 0., 0., 1. /) 
    ENDIF    
    IF(ABS(nv(2)-1.0) < 1e-14) THEN
        sv = (/ 0., 0., 1. /)
        tv = (/ 1., 0., 0.  /) 
    ENDIF    
    IF(ABS(nv(3)-1.0) < 1e-14) THEN
        sv = (/ 1., 0., 0. /)
        tv = (/ 0., 1., 0.  /) 
    ENDIF  
    !
    ! The exact Riemann solver of Godunov 
    !
    TM = 0.0    ! rotation matrix 
    TM(1,1:9) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
    TM(2,1:9) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
    TM(3,1:9) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
    TM(4,1:9) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
    TM(5,1:9) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
    TM(6,1:9) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
    TM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
    TM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
    TM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)  
    TM(10,10) = 1.0 
    TM(11,11) = 1.0 
    TM(12,12) = 1.0 
    TM(13,13) = 1.0 
    TM(14,14) = 1.0 
    !
    iTM = 0.0   ! inverse of the rotation matrix 
    iTM(1,1:9) = (/ nv(1)**2   , nv(2)**2   , nv(3)**2   , 2*nv(1)*nv(2)          , 2*nv(3)*nv(2)          , 2*nv(1)*nv(3)          , 0., 0., 0. /) 
    iTM(2,1:9) = (/ sv(1)**2   , sv(2)**2   , sv(3)**2   , 2*sv(1)*sv(2)          , 2*sv(3)*sv(2)          , 2*sv(1)*sv(3)          , 0., 0., 0. /) 
    iTM(3,1:9) = (/ tv(1)**2   , tv(2)**2   , tv(3)**2   , 2*tv(1)*sv(2)          , 2*tv(3)*tv(2)          , 2*tv(1)*tv(3)          , 0., 0., 0. /) 
    iTM(4,1:9) = (/ nv(1)*sv(1), nv(2)*sv(2), nv(3)*sv(3), nv(1)*sv(2)+sv(1)*nv(2), nv(3)*sv(2)+sv(3)*nv(2), nv(1)*sv(3)+sv(1)*nv(3), 0., 0., 0. /) 
    iTM(5,1:9) = (/ sv(1)*tv(1), sv(2)*tv(2), sv(3)*tv(3), sv(1)*tv(2)+tv(1)*sv(2), sv(3)*tv(2)+tv(3)*sv(2), sv(1)*tv(3)+tv(1)*sv(3), 0., 0., 0. /) 
    iTM(6,1:9) = (/ nv(1)*tv(1), nv(2)*tv(2), nv(3)*tv(3), nv(1)*tv(2)+tv(1)*nv(2), nv(3)*tv(2)+tv(3)*nv(2), nv(1)*tv(3)+tv(1)*nv(3), 0., 0., 0. /) 
    iTM(7,1:9) = (/ 0.,0., 0., 0., 0., 0., nv(1), nv(2), nv(3) /)     
    iTM(8,1:9) = (/ 0.,0., 0., 0., 0., 0., sv(1), sv(2), sv(3) /)     
    iTM(9,1:9) = (/ 0.,0., 0., 0., 0., 0., tv(1), tv(2), tv(3) /)   
    iTM(10,10) = 1.0 
    iTM(11,11) = 1.0 
    iTM(12,12) = 1.0         
    iTM(13,13) = 1.0         
    iTM(14,14) = 1.0     
    !
    !test = MATMUL( TM, TRANSPOSE(TM) )  
    !
    Q = 0.5*(QL+QR) 
    lam  = Q(10) 
    mu   = Q(11) 
    rho  = Q(12) 
    irho = 1./rho 
    Qrot = MATMUL( iTM, Q) 
    alpha = Q(13) 
    IF(Q(13)<=1e-3) THEN
        ialpha = 0.0
        uv = 0.0 
    ELSE
        ialpha = 1./Q(13)
        uv     = Qrot(7:9)*ialpha 
    ENDIF 
    ! 
    s11 = Qrot(1) 
    s12 = Qrot(4) 
    s13 = Qrot(6) 
        
    !
    cp = SQRT((lam+2*mu)/rho) 
    cs = SQRT(mu/rho) 
    !
    RM = 0.0 
    RM(1,1)   = rho*cp**2
    RM(1,10)  = -2*s11
    RM(1,13)  = rho*cp**2
    RM(2,1)   = -rho*(-cp**2+2*cs**2)
    RM(2,4)   = 1
    RM(2,13)  = -rho*(-cp**2+2*cs**2)
    RM(3,1)   = -rho*(-cp**2+2*cs**2)
    RM(3,5)   = 1
    RM(3,13)  = -rho*(-cp**2+2*cs**2)
    RM(4,2)   = rho*cs**2
    RM(4,10)  = -2*s12
    RM(4,12)  = rho*cs**2
    RM(5,6)   = 1
    RM(6,3)   = rho*cs**2
    RM(6,10)  = -2*s13
    RM(6,11)  = rho*cs**2
    RM(7,1)   = cp
    RM(7,10)  = uv(1) 
    RM(7,13)  = -cp
    RM(8,2)   = cs
    RM(8,10)  = uv(2) 
    RM(8,12)  = -cs
    RM(9,3)   = cs
    RM(9,10)  = uv(3) 
    RM(9,11)  = -cs
    RM(10,9)  = 1
    RM(11,8)  = 1
    RM(12,7)  = 1
    RM(13,10) = 1
    RM(14,14) = 1 
    !
    iRM = 0.0 
    iRM(1,1)   = irho/cp**2/2
    iRM(1,7)   = 1/cp/2
    iRM(1,13)  = -1/cp**2/rho*(cp*rho*uv(1)-2*s11)/2
    iRM(2,4)   = irho/cs**2/2
    iRM(2,8)   = 1/cs/2
    iRM(2,13)  = -1/cs**2/rho*(uv(2)*cs*rho-2*s12)/2
    iRM(3,6)   = irho/cs**2/2
    iRM(3,9)   = 1/cs/2
    iRM(3,13)  = -1/cs**2/rho*(uv(3)*cs*rho-2*s13)/2
    iRM(4,1)   = 1/cp**2*(-cp**2+2*cs**2)
    iRM(4,2)   = 1
    iRM(4,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
    iRM(5,1)   = 1/cp**2*(-cp**2+2*cs**2)
    iRM(5,3)   = 1
    iRM(5,13)  = 1/cp**2*2*s11*(-cp**2+2*cs**2) 
    iRM(6,5)   = 1
    iRM(7,12)  = 1
    iRM(8,11)  = 1
    iRM(9,10)  = 1
    iRM(10,13) = 1. 
    iRM(11,6)  = irho/cs**2/2
    iRM(11,9)  = -1/cs/2
    iRM(11,13) = 1/cs**2/rho*(uv(3)*cs*rho+2*s13)/2
    iRM(12,4)  = irho/cs**2/2
    iRM(12,8)  = -1/cs/2
    iRM(12,13) = 1/cs**2/rho*(uv(2)*cs*rho+2*s12)/2
    iRM(13,1)  = irho/cp**2/2
    iRM(13,7)  = -1/cp/2
    iRM(13,13) = 1/cp**2/rho*(cp*rho*uv(1)+2*s11)/2
    iRM(14,14) = 1 
    !
    R = MATMUL(   TM,  RM ) 
    iR = MATMUL( iRM, iTM ) 
    !
    LM = 0. 
    LM(1,1)   = -cp
    LM(2,2)   = -cs
    LM(3,3)   = -cs
    LM(11,11) = +cs 
    LM(12,12) = +cs
    LM(13,13) = +cp
    LM(14,14) = 0.0 

    !
    DM = MATMUL( R, MATMUL( ABS(LM), iR )  )  
    !
    CONTINUE
    !
#else 
    !
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    !
    ! Otherwise, use the simple Rusanov method 
    !
    CALL PDEEigenvalues(LL,QavL,parL,nv) 
    CALL PDEEigenvalues(LR,QavR,parR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
    !
    DM = smax*Id 
    !
#endif         
    !
END SUBROUTINE DissipationMatrix