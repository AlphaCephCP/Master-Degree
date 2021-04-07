#define THIRD_ORDER_INITIAL_GUESS  
    
    
SUBROUTINE ADERSpaceTimePredictorNonlinear(lqhi,lFhi,lShi,lQbnd,lFbnd,luh,lpar)
    USE typesDef
    IMPLICIT NONE
    ! Argument list
    REAL, INTENT(IN)   ::  luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom
    REAL, INTENT(IN)   ::  lpar(nParam,nDOF(1),nDOF(2),nDOF(3))           ! spatial degrees of freedom
    REAL, INTENT(OUT)  ::  lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged space-time degrees of freedom
    REAL, INTENT(OUT)  ::  lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))           ! time-averaged nonlinear flux tensor in each space DOF
    REAL, INTENT(OUT)  ::  lShi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged nonlinear source vector in each space DOF
    REAL, INTENT(OUT)  ::  lqbnd(nVar,nDOF(2),nDOF(3),6)                  ! time-averaged space-time degrees of freedom
    REAL, INTENT(OUT)  ::  lFbnd(nVar,nDOF(2),nDOF(3),6)                  ! time-averaged nonlinear flux tensor in each space-time DOF
    ! Local variables
    INTEGER    :: i,j,k,l,iVar,iDim, iter
    REAL       :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side
    REAL       :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side
    REAL       :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
#ifndef NOFLUX 
    REAL       :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
#endif 
    REAL       :: lSh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! nonlinear source vector in each space-time DOF
    REAL       :: aux(d), w, ltime                                         ! auxiliary variables
    REAL       :: gradQ(nVar,d), BgradQ(nVar), src(nVar)
    REAL       :: Src_BgradQ(nVar), qqt(nVar)                                
    REAL       :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))             ! old space-time degrees of freedom
    REAL       :: lfx(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative fx
    REAL       :: lfy(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative fy
    REAL       :: lfz(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative fz 
    REAL       :: lqx(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qx of q
    REAL       :: lqy(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qy of q
    REAL       :: lqz(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qz of q
    REAL       :: k1RK(nVar,nDOF(1),nDOF(2),nDOF(3))                       ! aux. k1 for RK initial guess 
    REAL       :: k2RK(nVar,nDOF(1),nDOF(2),nDOF(3))                       ! aux. k2 for RK initial guess 
    REAL       :: k3RK(nVar,nDOF(1),nDOF(2),nDOF(3))                       ! aux. k3 for RK initial guess 
    REAL       :: yt(nVar), ytt(nVar), yttt(nVar) 
    !REAL      ::  lqt_debug(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))         ! time derivative qt of q
    REAL       :: res                                                      ! residual
    REAL, PARAMETER    :: tol = 1e-7                                       ! tolerance
#ifdef FIRST_ORDER_INITIAL_GUESS 
    INTEGER, PARAMETER :: nPicard = MAX(1,N)     
#endif 
#ifdef SECOND_ORDER_INITIAL_GUESS  
    INTEGER, PARAMETER :: nPicard = MAX(1,N-1)
#endif     
#ifdef THIRD_ORDER_INITIAL_GUESS 
    INTEGER, PARAMETER :: nPicard = MAX(1,N-2)
#endif     
#ifdef FOURTH_ORDER_INITIAL_GUESS 
    INTEGER, PARAMETER :: nPicard = MAX(1,N-3)  
#endif     
    !
#ifdef PRIMADER
    CALL ADERSpaceTimePredictorPrim(lqhi,lFhi,lShi,lQbnd,lFbnd,luh,lpar) 
    RETURN 
#endif 
    !
    lqx = 0.0
    lqy = 0.0
    lqz = 0.0 
    lFbnd = 0.0 
    ! --------------------------------------------------------------------------------------------------------------------- 
#ifdef SECOND_ORDER_INITIAL_GUESS 
    ltime = time
    lfx = 0.0
    lfy = 0.0
    lfz = 0.0 
    Src_BgradQ = 0.0 
    ! Second order accurate initial guess (MUSCL-Hancock type) 
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),luh(:,i,j,k),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif     
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( luh(:,:,j,k), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( luh(:,i,:,k), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( luh(:,i,j,:), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,luh(:,i,j,k),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                ! use the time derivative to initialize a first order polynomial in time in each spatial DOF 
                DO l = 1, nDOF(0) 
                    lqh(:,i,j,k,l) = luh(:,i,j,k) + dt*xiGPN(l)*qqt 
                ENDDO      
                ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
                ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
                ! inside the Picard loop (DO iter = 1, N+1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                DO iVar = 1, nEvolve 
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
                ENDDO                
            ENDDO
        ENDDO
    ENDDO
#endif 
    ! --------------------------------------------------------------------------------------------------------------------- 
#ifdef THIRD_ORDER_INITIAL_GUESS 
    ltime = time
    lfx = 0.0
    lfy = 0.0
    lfz = 0.0 
    Src_BgradQ = 0.0 
    ! Third order accurate initial guess (Runge-Kutta type) 
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),luh(:,i,j,k),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif             
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( luh(:,:,j,k), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( luh(:,i,:,k), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( luh(:,i,j,:), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,luh(:,i,j,k),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                k1RK(:,i,j,k) = qqt 
                lqh(:,i,j,k,1) = luh(:,i,j,k) + dt*qqt 
            ENDDO
        ENDDO
    ENDDO  
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),lqh(:,i,j,k,1),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif             
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,lqh(:,i,j,k,1),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                ! use the time derivative to initialize a first order polynomial in time in each spatial DOF 
                DO l = 1, nDOF(0) 
                    lqh(:,i,j,k,l) = luh(:,i,j,k) + dt*xiGPN(l)*k1RK(:,i,j,k) + 0.5*dt*xiGPN(l)**2*( qqt - k1RK(:,i,j,k) )  
                ENDDO      
                ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
                ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
                ! inside the Picard loop (DO iter = 1, N+1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                DO iVar = 1, nEvolve 
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
                ENDDO                
            ENDDO
        ENDDO
    ENDDO
#endif     
    ! --------------------------------------------------------------------------------------------------------------------- 
#ifdef FOURTH_ORDER_INITIAL_GUESS 
    ltime = time
    lfx = 0.0
    lfy = 0.0
    lfz = 0.0 
    Src_BgradQ = 0.0 
    ! Fourth order accurate initial guess (CERK Runge-Kutta scheme, B. Owren, M. Zennaro, Derivation of efficient continuous explicit Runge–Kutta methods, SIAM J. Sci. Stat. Comput. 13 (1992) 1488–1501. ) 
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),luh(:,i,j,k),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif     
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( luh(:,:,j,k), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( luh(:,i,:,k), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( luh(:,i,j,:), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,luh(:,i,j,k),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                k1RK(:,i,j,k) = qqt 
                lqh(:,i,j,k,1) = luh(:,i,j,k) + 12./23.*dt*qqt 
            ENDDO
        ENDDO
    ENDDO  
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),lqh(:,i,j,k,1),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif     
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,lqh(:,i,j,k,1),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                k2RK(:,i,j,k) = qqt 
                lqh(:,i,j,k,2) = luh(:,i,j,k) - 68./375.*dt*k1RK(:,i,j,k) + 368./375.*dt*qqt  
            ENDDO
        ENDDO
    ENDDO
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),lqh(:,i,j,k,2),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif     
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,2), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,2), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,2), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,lqh(:,i,j,k,2),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 
                k3RK(:,i,j,k) = qqt 
                lqh(:,i,j,k,3) = luh(:,i,j,k) + 31./144.*dt*k1RK(:,i,j,k) + 529./1152.*dt*k2RK(:,i,j,k) + 125./384.*dt*qqt 
            ENDDO
        ENDDO
    ENDDO    
#ifndef NOFLUX 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                CALL PDEFlux(lFh(:,:,i,j,k,1),lqh(:,i,j,k,3),lpar(:,i,j,k)) 
            ENDDO
        ENDDO
    ENDDO 
#endif     
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
#ifndef NOFLUX 
            lfx(1:nEvolve,:,j,k) = 1.0/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,1), TRANSPOSE(dudx) )         ! compute the gradient of the flux  
#endif             
#ifdef NONCONSERVATIVE 
            lqx(:,:,j,k,1) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,3), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif             
        ENDDO
    ENDDO
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
#ifndef NOFLUX 
                lfy(1:nEvolve,i,:,k) = 1.0/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux 
#endif             
#ifdef NONCONSERVATIVE 
                lqy(:,i,:,k,1) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,3), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                             
            ENDDO
        ENDDO        
    ENDIF
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
#ifndef NOFLUX 
                lfz(1:nEvolve,i,j,:) = 1.0/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,1), TRANSPOSE(dudx) )     ! compute the gradient of the flux                     
#endif             
#ifdef NONCONSERVATIVE 
                lqz(:,i,j,:,1) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,3), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif                                             
            ENDDO
        ENDDO                
    ENDIF    
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                qqt = -lfx(:,i,j,k) - lfy(:,i,j,k) - lfz(:,i,j,k) 
#ifdef NONCONSERVATIVE 
                gradQ(:,1) = lqx(:,i,j,k,1) 
                gradQ(:,2) = lqy(:,i,j,k,1) 
                gradQ(:,3) = lqz(:,i,j,k,1) 
                CALL PDEFusedSrcNCP(Src_BgradQ,lqh(:,i,j,k,3),gradQ,lpar(:,i,j,k),ltime)    
                qqt = qqt + Src_BgradQ 
#endif 

                !yt   = k1RK(:,i,j,k) 
                !ytt  = 2.0*(-65./48.*k1RK(:,i,j,k)+529./384.*k2RK(:,i,j,k)+125./128.*k3RK(:,i,j,k)-qqt)/dt  
                !yttt = 6.0*( 41./72.*k1RK(:,i,j,k)-529./576.*k2RK(:,i,j,k)-125./192.*k3RK(:,i,j,k)+qqt)/dt**2  

                ! use the time derivative to initialize a first order polynomial in time in each spatial DOF 
                DO l = 1, nDOF(0) 
                    lqh(:,i,j,k,l) = luh(:,i,j,k) + dt*xiGPN(l)*k1RK(:,i,j,k) + dt*xiGPN(l)**2*( -65./48.*k1RK(:,i,j,k) + 529./384.*k2RK(:,i,j,k) + 125./128.*k3RK(:,i,j,k) - qqt ) + dt*xiGPN(l)**3*( 41./72.*k1RK(:,i,j,k) - 529./576.*k2RK(:,i,j,k) - 125./192.*k3RK(:,i,j,k) + qqt  ) 
                    !lqh(:,i,j,k,l) = luh(:,i,j,k) + dt*(41./72.*xiGPN(l)**3-65./48.*xiGPN(l)**2+xiGPN(l))*k1RK(:,i,j,k) + dt*(-529./576.*xiGPN(l)**3+529./384.*xiGPN(l)**2)*k2RK(:,i,j,k) + dt*(-125./192.*xiGPN(l)**3+125./128.*xiGPN(l)**2)*k3RK(:,i,j,k) + dt*(xiGPN(l)**3-xiGPN(l)**2)*qqt 
                ENDDO      
                ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
                ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
                ! inside the Picard loop (DO iter = 1, N+1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                DO iVar = 1, nEvolve 
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
                ENDDO                
            ENDDO
        ENDDO
    ENDDO
    ! 
#endif     
    ! --------------------------------------------------------------------------------------------------------------------- 
#ifdef FIRST_ORDER_INITIAL_GUESS 
   ! Simple first order accurate initial guess (constant extrapolation in time)  
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)                
                ! Trivial initial guess (can be significantly improved)
                DO iVar = 1, nVar
                    lqh(iVar,i,j,k,:) = luh(iVar,i,j,k)
                ENDDO
                ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
                ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
                ! inside the Picard loop (DO iter = 1, N+1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                DO iVar = 1, nEvolve 
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
                ENDDO
                !
            ENDDO
        ENDDO
    ENDDO
#endif 
    !
    ! Discrete Picard iterations. This set of nested loops should (theoretically) be a dream for vectorization, since they are rather independent...
    DO iter = 1, nPicard 
        ! save old space-time DOF
        lqhold = lqh
        DO l = 1, nDOF(0) ! loop over DOF in time
            ltime = time + xiGPN(l)*dt   
            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
#ifndef NOFLUX
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        CALL PDEFlux(lFh(:,:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k))
                        !CALL PDESource(lSh(:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k),ltime)
                    ENDDO
                ENDDO
            ENDDO
#endif 
            ! Compute the "derivatives" (contributions of the stiffness matrix)
            ! x direction (independent from the y and z derivatives)
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    aux = (/ wGPN(l), wGPN(j), wGPN(k) /)
#ifndef NOFLUX 
                    rhs(1:nEvolve,:,j,k,l) = rhs0(1:nEvolve,:,j,k,l) - PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(1:nEvolve,1,:,j,k,l), Kxi )
#else
                    rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) 
#endif 
#ifdef NONCONSERVATIVE 
                    lqx(:,:,j,k,l) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,l), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif 
                ENDDO
            ENDDO
            ! y direction (independent from the x and z derivatives) - should not be used for 1D
            IF(nDim>=2) THEN
                DO k = 1, nDOF(3)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(l), wGPN(i), wGPN(k) /)
#ifndef NOFLUX 
                        rhs(1:nEvolve,i,:,k,l) = rhs(1:nEvolve,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(1:nEvolve,2,i,:,k,l), Kxi )
#endif 
#ifdef NONCONSERVATIVE 
                        lqy(:,i,:,k,l) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif 
                    ENDDO
                ENDDO
            ENDIF
            ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
            IF(nDim>=3) THEN
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(l), wGPN(i), wGPN(j) /)
#ifndef NOFLUX 
                        rhs(1:nEvolve,i,j,:,l) = rhs(1:nEvolve,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(1:nEvolve,3,i,j,:,l), Kxi )
#endif 
#ifdef NONCONSERVATIVE 
                        lqz(:,i,j,:,l) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif 
                    ENDDO
                ENDDO
            ENDIF
#if defined(NONCONSERVATIVE) || defined(SOURCE)             
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(i), wGPN(j), wGPN(k) /) 
                        gradQ(:,1) = lqx(:,i,j,k,l) 
                        gradQ(:,2) = lqy(:,i,j,k,l) 
                        gradQ(:,3) = lqz(:,i,j,k,l) 
                        !CALL PDENCP(BgradQ,lqh(:,i,j,k,l),gradQ,lpar(:,i,j,k))    
                        !lSh(:,i,j,k,l) = lSh(:,i,j,k,l) - BgradQ 
                        CALL PDEFusedSrcNCP(Src_BgradQ,lqh(:,i,j,k,l),gradQ,lpar(:,i,j,k),ltime) 
                        lSh(:,i,j,k,l) = Src_BgradQ 
                        rhs(:,i,j,k,l) = rhs(:,i,j,k,l) + PRODUCT(aux(1:nDim))*wGPN(l)*dt*lSh(:,i,j,k,l)  
                    ENDDO
                ENDDO
            ENDDO
#endif             
            !
        ENDDO ! end loop over time DOF
        !
        ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration
        !
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                    lqh(1:nEvolve,i,j,k,:) = 1./(PRODUCT(aux(1:nDim)))*MATMUL( rhs(1:nEvolve,i,j,k,:), TRANSPOSE(iK1) )
                    !lqt_debug(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
                    continue
                ENDDO
            ENDDO
        ENDDO
        !
        ! We can stop the iterations if a certain tolerance has been reached. If you do not like this unpredictable feature (it depends on the solution of the PDE)
        ! simply comment the lines below, so each element will always do the same number of iterations in the predictor step, i.e. the same number of operations
        !
        res = SQRT(SUM((lqh-lqhold)**2))
        IF(res.LT.tol) THEN
            EXIT
        ENDIF
        !
    ENDDO
    !
    ! Immediately compute the time-averaged space-time polynomials
    !
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lqhi(:,i,j,k) = MATMUL( lqh(:,i,j,k,:), wGPN )
#if defined(NONCONSERVATIVE) || defined(SOURCE)             
                lShi(:,i,j,k) = MATMUL( lSh(:,i,j,k,:), wGPN )
#endif                 
#ifndef NOFLUX 
                DO iDim = 1, nDim
                    lFhi(1:nEvolve,iDim,i,j,k) = MATMUL( lFh(1:nEvolve,iDim,i,j,k,:), wGPN )
                ENDDO
#endif 
            ENDDO
        ENDDO
    ENDDO
    !
    ! Compute the bounday-extrapolated values for Q and F*n
    !
    lQbnd = 0.
#ifndef NOFLUX 
    lFbnd = 0.
#endif 
    ! x-direction: face 1 (left) and face 2 (right)
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
            lQbnd(:,j,k,2) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
#ifndef NOFLUX 
            lFbnd(1:nEvolve,j,k,1) = MATMUL( lFhi(1:nEvolve,1,:,j,k), FLCoeff )   ! left
            lFbnd(1:nEvolve,j,k,2) = MATMUL( lFhi(1:nEvolve,1,:,j,k), FRCoeff )   ! right
#endif 
        ENDDO
    ENDDO
    ! y-direction: face 3 (left) and face 4 (right)
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k),   FLCoeff )   ! left
                lQbnd(:,i,k,4) = MATMUL( lqhi(:,i,:,k),   FRCoeff )   ! right
#ifndef NOFLUX 
                lFbnd(1:nEvolve,i,k,3) = MATMUL( lFhi(1:nEvolve,2,i,:,k), FLCoeff )   ! left
                lFbnd(1:nEvolve,i,k,4) = MATMUL( lFhi(1:nEvolve,2,i,:,k), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    ! z-direction: face 5 (left) and face 6 (right)
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:),   FLCoeff )   ! left
                lQbnd(:,i,j,6) = MATMUL( lqhi(:,i,j,:),   FRCoeff )   ! right
#ifndef NOFLUX 
                lFbnd(1:nEvolve,i,j,5) = MATMUL( lFhi(1:nEvolve,3,i,j,:), FLCoeff )   ! left
                lFbnd(1:nEvolve,i,j,6) = MATMUL( lFhi(1:nEvolve,3,i,j,:), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    !
    CONTINUE
    !
    END SUBROUTINE ADERSpaceTimePredictorNonlinear


    SUBROUTINE ADERSpaceTimePredictorLinear(lqhi,lShi,lQbnd,luh,lpar,iElem)
    USE typesDef
    IMPLICIT NONE
    ! Argument list
    INTEGER, INTENT(IN) :: iElem                                        ! element number (needed for the point sources) 
    REAL, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom
    REAL, INTENT(IN)  :: lpar(nParam,nDOF(1),nDOF(2),nDOF(3))           ! spatial degrees of freedom
    REAL, INTENT(OUT) :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged space-time degrees of freedom
    REAL, INTENT(OUT) :: lShi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged nonlinear flux tensor in each space-time DOF
    REAL, INTENT(OUT) :: lqbnd(nVar,nDOF(2),nDOF(3),6)                  ! time-averaged space-time degrees of freedom
    ! Local variables
    INTEGER :: i,j,k,l,ii,jj,kk,iVar,iDim
    REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! temporary work array
    REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0)+1)              ! space-time degrees of freedom
    REAL    :: BgradQ(nVar)                                             ! non-conservative product 
    REAL    :: aux(d), w ,sigma(nVar)                                   ! auxiliary variables
    REAL    :: gradQ(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))            ! spatial gradient of q
    !REAL    :: lqt_debug(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))         ! time derivative qt of q
    REAL    :: dtavFac                                                  ! integral average of the time power terms of the Taylor series
    REAL, PARAMETER :: tol = 1e-6                                       ! tolerance
    !
    ! Init the time derivatives with the point sources (if applicable) 
    ! 
    lqh = 0. 
    gradQ = 0.0 
    !
!    DO i = 1, nPointSource 
!        IF(iElem.EQ.PointSrc(i)%iElem) THEN
!           DO l = 2, nDOF(0)  
!            DO kk = 1, nDOF(3)
!             DO jj = 1, nDOF(2) 
!              DO ii = 1, nDOF(1) 
!               DO iVar = 1, nVar
!                  !
!                  ! Multiply source with (M)^(-1) to get the next higher time derivative
!                  !
!                  CALL PointSource(sigma,time,PointSrc(i)%waveform) 
!                  aux = (/ wGPN(ii), wGPN(jj), wGPN(kk) /)
!                  lqh(iVar,ii,jj,kk,l) = lqh(iVar,ii,jj,kk,l) + PointSrc(i)%phi(ii,jj,kk)/(PRODUCT(aux(1:nDim))*PRODUCT(dx(1:nDim)))*PointSrc(i)%sigma(l-1,iVar)/(dt**(l-2))
!               ENDDO
!              ENDDO
!             ENDDO 
!            ENDDO 
!           ENDDO 
!        ENDIF
!    ENDDO 

    !
    ! The zeroth time derivative (time dof number 1) is the initial condition
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                DO iVar = 1, nVar
                    lqh(iVar,i,j,k,1) = lqh(iVar,i,j,k,1) + luh(iVar,i,j,k)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !
    ! For linear PDE, the fastest space-time predictor is the good old Cauchy-Kovalewski procedure
    !
    DO l = 1, nDOF(0) 
        ! Compute the derivatives in x direction (independent from the y and z derivatives)
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                aux = (/ 1., wGPN(j), wGPN(k) /)
                gradQ(:,1,:,j,k,l) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,l), TRANSPOSE(dudx) )     
            ENDDO
        ENDDO
        ! y direction (independent from the x and z derivatives) - should not be used for 1D
        IF(nDim>=2) THEN
            DO k = 1, nDOF(3)
                DO i = 1, nDOF(1)
                    aux = (/ 1., wGPN(i), wGPN(k) /)
                    gradQ(:,2,i,:,k,l) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,l), TRANSPOSE(dudx) )   
                ENDDO
            ENDDO
        ENDIF
        ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
        IF(nDim>=3) THEN
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    aux = (/ 1., wGPN(i), wGPN(j) /)
                    gradQ(:,3,i,j,:,l) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,l), TRANSPOSE(dudx) )  
                ENDDO
            ENDDO
        ENDIF
        ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    CALL PDENCP(BgradQ,lqh(:,i,j,k,l),gradQ(:,:,i,j,k,l),lpar(:,i,j,k)) 
                    lqh(:,i,j,k,l+1) = lqh(:,i,j,k,l+1) - BgradQ(:) 
                ENDDO
            ENDDO
        ENDDO
        !
        CONTINUE
        !
        !
        ! Multiply with (M)^(-1) to get the next higher time derivative
        !
!        DO k = 1, nDOF(3)
!            DO j = 1, nDOF(2)
!                DO i = 1, nDOF(1)
!                    aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
!                    lqh(:,i,j,k,l+1) = 1./(PRODUCT(aux(1:nDim)))*rhs(:,i,j,k)
!                    lqt_debug(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
!                ENDDO
!            ENDDO
!        ENDDO
        !
    ENDDO
    !
    ! Immediately compute the time-averaged space-time polynomials
    !
    lqhi = lqh(:,:,:,:,1)
    lShi = lqh(:,:,:,:,2)
    dtavFac = 0.5*dt  
    DO l = 2, nDOF(0)
        lqhi(:,:,:,:) = lqhi(:,:,:,:) + dtavFac*lqh(:,:,:,:,l)
        lShi(:,:,:,:) = lShi(:,:,:,:) + dtavFac*lqh(:,:,:,:,l+1)
        dtavFac = dtavFac*dt/REAL(l+1)
    ENDDO
    !
    ! Compute the bounday-extrapolated values for Q and F*n
    !
    lQbnd = 0.
    ! x-direction: face 1 (left) and face 2 (right)
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
            lQbnd(:,j,k,2) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
        ENDDO
    ENDDO
    ! y-direction: face 3 (left) and face 4 (right)
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k),   FLCoeff )   ! left
                lQbnd(:,i,k,4) = MATMUL( lqhi(:,i,:,k),   FRCoeff )   ! right
            ENDDO
        ENDDO
    ENDIF
    ! z-direction: face 5 (left) and face 6 (right)
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:),   FLCoeff )   ! left
                lQbnd(:,i,j,6) = MATMUL( lqhi(:,i,j,:),   FRCoeff )   ! right
            ENDDO
        ENDDO
    ENDIF
    !
    CONTINUE
    !
    END SUBROUTINE ADERSpaceTimePredictorLinear



    SUBROUTINE ADERSpaceTimePredictorPrim(lqhi,lFhi,lShi,lQbnd,lFbnd,luh,lpar)
    USE typesDef
    IMPLICIT NONE
    ! Argument list
    REAL, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom
    REAL, INTENT(IN)  :: lpar(nParam,nDOF(1),nDOF(2),nDOF(3))           ! spatial degrees of freedom
    REAL, INTENT(OUT) :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged space-time degrees of freedom
    REAL, INTENT(OUT) :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))           ! time-averaged nonlinear flux tensor in each space DOF
    REAL, INTENT(OUT) :: lShi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged nonlinear source vector in each space DOF
    REAL, INTENT(OUT) :: lqbnd(nVar,nDOF(2),nDOF(3),6)                  ! time-averaged space-time degrees of freedom
    REAL, INTENT(OUT) :: lFbnd(nVar,nDOF(2),nDOF(3),6)                  ! time-averaged nonlinear flux tensor in each space-time DOF
    ! Local variables
    INTEGER :: i,j,k,l,iVar,iDim, iter, iErr 
    REAL    :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side
    REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side
    REAL    :: lvh(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! primitive variables at time t^n 
    REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
    REAL    :: lph(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
    REAL    :: lphold(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))             ! old space-time degrees of freedom
#ifndef NOFLUX 
    REAL    :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
#endif 
    REAL    :: lSh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! nonlinear source vector in each space-time DOF
    REAL    :: aux(d), w, ltime                                         ! auxiliary variables
    REAL    :: gradQ(nVar,d), BgradQ(nVar), src(nVar)
    REAL    :: Src_BgradQ(nVar), dVdQ(nVar,nVar), dVdt(nVar)                                   
    REAL    :: lfx(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qx of q
    REAL    :: lgy(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qy of q
    REAL    :: lhz(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qz of q
    REAL    :: lqx(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qx of q
    REAL    :: lqy(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qy of q
    REAL    :: lqz(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qz of q
    !REAL    :: lqt_debug(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))         ! time derivative qt of q
    REAL    :: res                                                      ! residual
    REAL, PARAMETER :: tol = 1e-6                                       ! tolerance
    !
    lqx = 0.0
    lqy = 0.0
    lqz = 0.0 
    lfx = 0.0 
    lgy = 0.0
    lhz = 0.0 
    !
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                CALL PDECons2Prim(lvh(:,i,j,k),luh(:,i,j,k),iErr) 
                ! Trivial initial guess (can be significantly improved)
                DO iVar = 1, nVar
                    lph(iVar,i,j,k,:) = lvh(iVar,i,j,k)
                ENDDO
                ! Compute the contribution of the initial condition uh to the time update. I prefer to compute it once
                ! and store it in rhs0, but if you think it is faster, you can also recompute this contribution
                ! inside the Picard loop (DO iter = 1, N+1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                DO iVar = 1, nVar
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*lvh(iVar,i,j,k)
                ENDDO
                !
            ENDDO
        ENDDO
    ENDDO
    !
    ! Discrete Picard iterations. This set of nested loops should (theoretically) be a dream for vectorization, since they are rather independent...
    DO iter = 1, N+1
        ! save old space-time DOF
        lphold = lph
        DO l = 1, nDOF(0) ! loop over DOF in time
            ltime = time + xiGPN(l)*dt   
            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
#ifndef NOFLUX
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        CALL PDEPrim2Cons(lqh(:,i,j,k,l),lph(:,i,j,k,l)) 
                        CALL PDEFluxPrim(lFh(:,:,i,j,k,l),lph(:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k))
                        !CALL PDESource(lSh(:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k),ltime)
                    ENDDO
                ENDDO
            ENDDO
#endif 
            ! Compute the "derivatives" (contributions of the stiffness matrix)
            ! x direction (independent from the y and z derivatives)
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
#ifndef NOFLUX 
                    lfx(:,:,j,k,l) = 1.0/dx(1)*MATMUL( lFh(:,1,:,j,k,l), TRANSPOSE(dudx) ) 
#endif 
#ifdef NONCONSERVATIVE 
                    lqx(:,:,j,k,l) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,l), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
#endif 
                ENDDO
            ENDDO
            ! y direction (independent from the x and z derivatives) - should not be used for 1D
            IF(nDim>=2) THEN
                DO k = 1, nDOF(3)
                    DO i = 1, nDOF(1)
#ifndef NOFLUX 
                        lgy(:,i,:,k,l) = 1.0/dx(2)*MATMUL( lFh(:,2,i,:,k,l), TRANSPOSE(dudx) ) 
#endif 
#ifdef NONCONSERVATIVE 
                        lqy(:,i,:,k,l) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif 
                    ENDDO
                ENDDO
            ENDIF
            ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
            IF(nDim>=3) THEN
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
#ifndef NOFLUX 
                        lhz(:,i,j,:,l) = 1.0/dx(3)*MATMUL( lFh(:,3,i,j,:,l), TRANSPOSE(dudx) ) 
#endif 
#ifdef NONCONSERVATIVE 
                        lqz(:,i,j,:,l) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
#endif 
                    ENDDO
                ENDDO
            ENDIF
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(i), wGPN(j), wGPN(k) /) 
                        gradQ(:,1) = lqx(:,i,j,k,l) 
                        gradQ(:,2) = lqy(:,i,j,k,l) 
                        gradQ(:,3) = lqz(:,i,j,k,l)                         
                        CALL PDEdVdtPrim( dVdt, -lfx(:,i,j,k,l) - lgy(:,i,j,k,l) - lhz(:,i,j,k,l), lph(:,i,j,k,l) ) 
#if defined(NONCONSERVATIVE) || defined(SOURCE)             
                        CALL PDEFusedSrcNCPPrim(Src_BgradQ,lph(:,i,j,k,l),gradQ,lpar(:,i,j,k),ltime) 
                        lSh(:,i,j,k,l) = lSh(:,i,j,k,l) + Src_BgradQ 
#endif             
                        rhs(:,i,j,k,l) = rhs0(:,i,j,k,l) + PRODUCT(aux(1:nDim))*wGPN(l)*dt*( lSh(:,i,j,k,l) + dVdt ) 
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO ! end loop over time DOF
        !
        ! Multiply with (K1)^(-1) to get the discrete time integral of the discrete Picard iteration
        !
        DO k = 1, nDOF(3)
            DO j = 1, nDOF(2)
                DO i = 1, nDOF(1)
                    aux = (/ wGPN(i), wGPN(j), wGPN(k) /)
                    lph(:,i,j,k,:) = 1./(PRODUCT(aux(1:nDim)))*MATMUL( rhs(:,i,j,k,:), TRANSPOSE(iK1) )
                    !lqt_debug(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
                    continue
                ENDDO
            ENDDO
        ENDDO
        !
        ! We can stop the iterations if a certain tolerance has been reached. If you do not like this unpredictable feature (it depends on the solution of the PDE)
        ! simply comment the lines below, so each element will always do the same number of iterations in the predictor step, i.e. the same number of operations
        !
        res = SQRT(SUM((lph-lphold)**2))
        IF(res.LT.tol) THEN
            EXIT
        ENDIF
        !
    ENDDO
    !
    ! Immediately compute the time-averaged space-time polynomials
    !
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lqhi(:,i,j,k) = MATMUL( lqh(:,i,j,k,:), wGPN )
#if defined(NONCONSERVATIVE) || defined(SOURCE)             
                lShi(:,i,j,k) = MATMUL( lSh(:,i,j,k,:), wGPN )
#endif                 
#ifndef NOFLUX 
                DO iDim = 1, nDim
                    lFhi(:,iDim,i,j,k) = MATMUL( lFh(:,iDim,i,j,k,:), wGPN )
                ENDDO
#endif 
            ENDDO
        ENDDO
    ENDDO
    !
    ! Compute the bounday-extrapolated values for Q and F*n
    !
    lQbnd = 0.
#ifndef NOFLUX 
    lFbnd = 0.
#endif 
    ! x-direction: face 1 (left) and face 2 (right)
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
            lQbnd(:,j,k,2) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
#ifndef NOFLUX 
            lFbnd(:,j,k,1) = MATMUL( lFhi(:,1,:,j,k), FLCoeff )   ! left
            lFbnd(:,j,k,2) = MATMUL( lFhi(:,1,:,j,k), FRCoeff )   ! right
#endif 
        ENDDO
    ENDDO
    ! y-direction: face 3 (left) and face 4 (right)
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k),   FLCoeff )   ! left
                lQbnd(:,i,k,4) = MATMUL( lqhi(:,i,:,k),   FRCoeff )   ! right
#ifndef NOFLUX 
                lFbnd(:,i,k,3) = MATMUL( lFhi(:,2,i,:,k), FLCoeff )   ! left
                lFbnd(:,i,k,4) = MATMUL( lFhi(:,2,i,:,k), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    ! z-direction: face 5 (left) and face 6 (right)
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:),   FLCoeff )   ! left
                lQbnd(:,i,j,6) = MATMUL( lqhi(:,i,j,:),   FRCoeff )   ! right
#ifndef NOFLUX 
                lFbnd(:,i,j,5) = MATMUL( lFhi(:,3,i,j,:), FLCoeff )   ! left
                lFbnd(:,i,j,6) = MATMUL( lFhi(:,3,i,j,:), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    !
    CONTINUE
    !
    END SUBROUTINE ADERSpaceTimePredictorPrim 


