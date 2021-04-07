SUBROUTINE ADERRiemannSolver(lQbndL,lFbndL,lparbndL,lQbndR,lFbndR,lparbndR,nv)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)     :: lQbndL(nVar,nDOF(2),nDOF(3))              ! left space-time degrees of freedom 
    REAL, INTENT(IN)     :: lQbndR(nVar,nDOF(2),nDOF(3))              ! right space-time degrees of freedom 
    REAL, INTENT(IN)     :: lparbndL(nParam,nDOF(2),nDOF(3))          ! left material parameters 
    REAL, INTENT(IN)     :: lparbndR(nParam,nDOF(2),nDOF(3))          ! right material parameters 
    REAL, INTENT(IN)     :: nv(d)                                     ! normal vector 
    REAL, INTENT(INOUT)  :: lFbndL(nVar,nDOF(2),nDOF(3))              ! left flux 
    REAL, INTENT(INOUT)  :: lFbndR(nVar,nDOF(2),nDOF(3))              ! right flux 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d), QavL(nVar), QavR(nVar), Id(nVar,nVar), smax
    REAL              :: parL(nParam), parR(nParam), Bn(nVar,nVar) 
    REAL              :: LL(nVar), LR(nVar), BMat(nVar,nVar), DM(nVar,nVar)   
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
    ! Here, we implement a very simple Rusanov scheme with scalar dissipation (smax*Id). 
    ! We can change this into a more sophisticated Osher or HLLEM Riemann solver whenever needed! 
    !
    CALL PDEEigenvalues(LL,QavL,parL,nv) 
    CALL PDEEigenvalues(LR,QavR,parR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
    !    
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    !
    CALL PDEMatrixB(Bn,0.5*(QavL+QavR),nv,0.5*(parL+parR)) 
    ! 
    DO k = 1, nDOF(3)
      DO j = 1, nDOF(2)
            lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) )      ! purely conservative flux
            lFbndR(:,j,k) = lFbndL(:,j,k) - 0.5*MATMUL(Bn, lQbndR(:,j,k) - lQbndL(:,j,k) )                          ! subtract the non-conservative product
            lFbndL(:,j,k) = lFbndL(:,j,k) + 0.5*MATMUL(Bn, lQbndR(:,j,k) - lQbndL(:,j,k) )                          ! add the non-conservative product 
        ENDDO
    ENDDO
    !
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
    !
#ifdef GODUNOV     
    !
    parav = 0.5*(parL+parR) 
    lambda = parav(1) 
    mu     = parav(2)
    rho    = parav(3) 
    cp = SQRT((lambda+2*mu)/rho) 
    cs = SQRT(mu/rho) 
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
    TM(1,:) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
    TM(2,:) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
    TM(3,:) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
    TM(4,:) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
    TM(5,:) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
    TM(6,:) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
    TM(7,:) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
    TM(8,:) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
    TM(9,:) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)     
    !
    !test = MATMUL( TM, TRANSPOSE(TM) )  
    !
    RM(1,:) = (/ lambda + 2*mu, 0., 0., 0., 0., 0., 0., 0., lambda + 2*mu /) 
    RM(2,:) = (/ lambda, 0., 0., 0., 1., 0., 0., 0., lambda /) 
    RM(3,:) = (/ lambda, 0., 0., 0., 0., 1., 0., 0., lambda /) 
    RM(4,:) = (/ 0.,     mu, 0., 0., 1., 0., 0., mu, 0.     /) 
    RM(5,:) = (/ 0.,     0., 0., 1., 0., 0., 0., 0., 0.     /)  
    RM(6,:) = (/ 0.,     0., mu, 0., 0., 0., mu, 0., 0.     /)  
    RM(7,:) = (/ cp,     0., 0., 0., 0., 0., 0., 0., -cp    /)  
    RM(8,:) = (/ 0.,     cs, 0., 0., 0., 0., 0.,-cs, 0.     /)  
    RM(9,:) = (/ 0.,     0., cs, 0., 0., 0.,-cs, 0., 0.     /)  
    !
    iRM = 0.0 
    iRM(1,1) = 1.0/(lambda+2*mu)/2.0
    iRM(1,7) = 1/cp/2
    iRM(2,4) = 1/mu/2
    iRM(2,8) = 1/cs/2
    iRM(3,6) = 1/mu/2
    iRM(3,9) = 1/cs/2
    iRM(4,5) = 1.0
    iRM(5,1) = -lambda/(lambda+2*mu)
    iRM(5,2) = 1.0
    iRM(6,1) = -lambda/(lambda+2*mu)
    iRM(6,3) = 1.0
    iRM(7,6) = 1.0/mu/2.0 
    iRM(7,9) = -1.0/cs/2.0
    iRM(8,4) = 1.0/mu/2.0
    iRM(8,8) = -1.0/cs/2.0
    iRM(9,1) = 1.0/(lambda+2*mu)/2.0
    iRM(9,7) = -1.0/cp/2.0
    !
    LM = 0.0
    LM(1,1) = -cp
    LM(2,2) = -cs
    LM(3,3) = -cs
    LM(7,7) = +cs
    LM(8,8) = +cs
    LM(9,9) = +cp 
    !
    !DO i = 1, 9
    !    LM(i,i) = cp 
    !ENDDO    
    !
    DM = MATMUL( TM, MATMUL(RM, MATMUL( ABS(LM), MATMUL( iRM, TRANSPOSE(TM) ) ) )  ) 
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