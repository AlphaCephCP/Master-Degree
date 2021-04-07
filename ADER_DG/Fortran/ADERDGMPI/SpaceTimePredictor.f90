    SUBROUTINE ADERSpaceTimePredictorNonlinear(lqhi,lFhi,lShi,lQbnd,lFbnd,luh,lpar)
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
    INTEGER :: i,j,k,l,iVar,iDim, iter
    REAL    :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side
    REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side
    REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom
    REAL    :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF
    REAL    :: lSh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! nonlinear source vector in each space-time DOF
    REAL    :: aux(d), w, ltime                                         ! auxiliary variables
    REAL    :: gradQ(nVar,d), BgradQ(nVar), src(nVar)   
    REAL    :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))             ! old space-time degrees of freedom
    REAL    :: lqx(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qx of q
    REAL    :: lqy(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qy of q
    REAL    :: lqz(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qz of q
    REAL    :: lqt(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! time derivative qt of q
    REAL    :: res                                                      ! residual
    REAL, PARAMETER :: tol = 1e-7                                      ! tolerance
    !
    lqx = 0.0
    lqy = 0.0
    lqz = 0.0 
    !
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
                DO iVar = 1, nVar
                    rhs0(iVar,i,j,k,:) = PRODUCT(aux(1:nDim))*F0(:)*luh(iVar,i,j,k)
                ENDDO
                !
            ENDDO
        ENDDO
    ENDDO
    !
    ! Discrete Picard iterations. This set of nested loops should (theoretically) be a dream for vectorization, since they are rather independent...
    DO iter = 1, N+1
        ! save old space-time DOF
        lqhold = lqh
        DO l = 1, nDOF(0) ! loop over DOF in time
            ltime = time + xiGPN(l)*dt   
            ! Compute the fluxes (once these fluxes are available, the subsequent operations are independent from each other)
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        CALL PDEFlux(lFh(:,:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k))
                        CALL PDESource(lSh(:,i,j,k,l),lqh(:,i,j,k,l),lpar(:,i,j,k),ltime)
                    ENDDO
                ENDDO
            ENDDO
            ! Compute the "derivatives" (contributions of the stiffness matrix)
            ! x direction (independent from the y and z derivatives)
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    aux = (/ wGPN(l), wGPN(j), wGPN(k) /)
                    rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) - PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,1,:,j,k,l), Kxi )
                    lqx(:,:,j,k,l) = 1.0/dx(1)*MATMUL( lqh(:,:,j,k,l), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
                ENDDO
            ENDDO
            ! y direction (independent from the x and z derivatives) - should not be used for 1D
            IF(nDim>=2) THEN
                DO k = 1, nDOF(3)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(l), wGPN(i), wGPN(k) /)
                        rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,2,i,:,k,l), Kxi )
                        lqy(:,i,:,k,l) = 1.0/dx(2)*MATMUL( lqh(:,i,:,k,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
                    ENDDO
                ENDDO
            ENDIF
            ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
            IF(nDim>=3) THEN
                DO j = 1, nDOF(2)
                    DO i = 1, nDOF(1)
                        aux = (/ wGPN(l), wGPN(i), wGPN(j) /)
                        rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,3,i,j,:,l), Kxi )
                        lqz(:,i,j,:,l) = 1.0/dx(3)*MATMUL( lqh(:,i,j,:,l), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
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
                        CALL PDENCP(BgradQ,lqh(:,i,j,k,l),gradQ,lpar(:,i,j,k))    
                        lSh(:,i,j,k,l) = lSh(:,i,j,k,l) - BgradQ 
                        rhs(:,i,j,k,l) = rhs(:,i,j,k,l) + PRODUCT(aux(1:nDim))*wGPN(l)*dt*lSh(:,i,j,k,l)  
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
                    lqh(:,i,j,k,:) = 1./(PRODUCT(aux(1:nDim)))*MATMUL( rhs(:,i,j,k,:), TRANSPOSE(iK1) )
                    !lqt(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
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
                lShi(:,i,j,k) = MATMUL( lSh(:,i,j,k,:), wGPN )
                DO iDim = 1, nDim
                    lFhi(:,iDim,i,j,k) = MATMUL( lFh(:,iDim,i,j,k,:), wGPN )
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !
    ! Compute the bounday-extrapolated values for Q and F*n
    !
    lQbnd = 0.
    lFbnd = 0.
    ! x-direction: face 1 (left) and face 2 (right)
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
            lQbnd(:,j,k,2) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
            lFbnd(:,j,k,1) = MATMUL( lFhi(:,1,:,j,k), FLCoeff )   ! left
            lFbnd(:,j,k,2) = MATMUL( lFhi(:,1,:,j,k), FRCoeff )   ! right
        ENDDO
    ENDDO
    ! y-direction: face 3 (left) and face 4 (right)
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k),   FLCoeff )   ! left
                lQbnd(:,i,k,4) = MATMUL( lqhi(:,i,:,k),   FRCoeff )   ! right
                lFbnd(:,i,k,3) = MATMUL( lFhi(:,2,i,:,k), FLCoeff )   ! left
                lFbnd(:,i,k,4) = MATMUL( lFhi(:,2,i,:,k), FRCoeff )   ! right
            ENDDO
        ENDDO
    ENDIF
    ! z-direction: face 5 (left) and face 6 (right)
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:),   FLCoeff )   ! left
                lQbnd(:,i,j,6) = MATMUL( lqhi(:,i,j,:),   FRCoeff )   ! right
                lFbnd(:,i,j,5) = MATMUL( lFhi(:,3,i,j,:), FLCoeff )   ! left
                lFbnd(:,i,j,6) = MATMUL( lFhi(:,3,i,j,:), FRCoeff )   ! right
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
    REAL    :: lqt(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! time derivative qt of q
    REAL    :: dtavFac                                                  ! integral average of the time power terms of the Taylor series
    REAL, PARAMETER :: tol = 1e-7                                       ! tolerance
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
!                    lqt(:,i,j,k,:) = 1.0/dt*MATMUL( lqh(:,i,j,k,:), TRANSPOSE(dudx) )         ! currently used only for debugging purposes, to check if derivatives are correctly computed
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






