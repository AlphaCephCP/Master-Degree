SUBROUTINE RKDG3D 
    USE typesDef
    IMPLICIT NONE
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif    
    ! Local variables
    INTEGER :: i,j,k,iElem,iFace,nRecompute, narg, bcflag  
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
    REAL          :: tCPU3, lx0(d), ldx(d), Int_h(nObs) 
    LOGICAL       :: dmpresult 
    CHARACTER(LEN=100)  :: argstr 
    ! 
    IF(myrank==0) THEN
        PRINT *, ' ************************************************ ' 
        PRINT *, ' **** RUNNING THE ADERDG3D CODE IN RKDG MODE **** ' 
        PRINT *, ' ************************************************ ' 
    ENDIF 
    !
#ifdef PARALLEL 
    tCPU1 = MPI_WTIME() 
#else
    CALL CPU_TIME(tCPU1) 
#endif 
    ! Main loop in time 
    DO timestep = 1, NMAX
        IF(time >= tend) THEN
            EXIT
        ENDIF 
        ! Compute the time step size according to the CFL condition 
        CALL CalcTimeStep 
        ! 
        ! Limiter stuff 
#ifdef NOLIMITER 
        nRecompute   = 0
        recompute(:) = 0 
#endif         
        !
        SELECT CASE(N)
        CASE(1,2)
            ! Kutta's third order scheme 
            CALL Lh(RKDGk1,uh) 
            CALL Lh(RKDGk2,uh+0.5*dt*RKDGk1) 
            CALL Lh(RKDGk3,uh-1.0*dt*RKDGk1+2.0*dt*RKDGk2) 
            uh = uh + dt/6.0*( RKDGk1 + 4.0*RKDGk2 + 1.0*RKDGk3  )             
        CASE(4)
            ! Fifth order Runge-Kutta Fehlberg method 
            CALL Lh( RKDGk1,uh ) 
            CALL Lh( RKDGk2,uh+dt*1./4.*RKDGk1 ) 
            CALL Lh( RKDGk3,uh+dt*( 3./32.*RKDGk1 +	9./32.*RKDGk2 ) )
            CALL Lh( RKDGk4,uh+dt*( 1932./2197.*RKDGk1 - 7200./2197.*RKDGk2 + 7296./2197.*RKDGk3 ) )
            CALL Lh( RKDGk5,uh+dt*( 439./216.*RKDGk1 - 8.*RKDGk2 + 3680./513. *RKDGk3 - 845./4104.*RKDGk4 ) ) 
            CALL Lh( RKDGk6,uh+dt*(   -8./27.*RKDGk1 + 2.*RKDGk2 - 3544./2565.*RKDGk3 + 1859./4104.*RKDGk4 -11./40.*RKDGk5 ) ) 
            uh = uh + dt*( 16./135.*RKDGk1 + 6656./12825.*RKDGk3 + 28561./56430.*RKDGk4 - 9./50.*RKDGk5 + 2./55.*RKDGk6 ) 
        CASE(5)
            ! Sixth order Runge-Kutta scheme of Butcher 
            CALL Lh( RKDGk1,uh ) 
            CALL Lh( RKDGk2,uh+dt*1./3.*RKDGk1 ) 
            CALL Lh( RKDGk3,uh+dt*(  2./3. *RKDGk2 ) ) 
            CALL Lh( RKDGk4,uh+dt*(  1./12.*RKDGk1 + 1./3.*RKDGk2 - 1./12.*RKDGk3 ) ) 
            CALL Lh( RKDGk5,uh+dt*( -1./16.*RKDGk1 + 9./8.*RKDGk2 - 3./16. *RKDGk3 - 3./8.*RKDGk4 ) ) 
            CALL Lh( RKDGk6,uh+dt*(                + 9./8.*RKDGk2 - 3./8.*RKDGk3 - 3./4.*RKDGk4 + 1./2.*RKDGk5 ) ) 
            CALL Lh( RKDGk7,uh+dt*(  9./44.*RKDGk1 - 9./11.*RKDGk2 + 63./44.*RKDGk3 + 18./11.*RKDGk4 - 16./11.*RKDGk6 ) ) 
            uh = uh + dt*( 11./120.*RKDGk1 + 27./40.*RKDGk3 + 27./40.*RKDGk4 - 4./15.*RKDGk5 - 4./15.*RKDGk6 + 11./120.*RKDGk7 ) 
        CASE DEFAULT            
            ! default choice: classic fourth order Runge-Kutta scheme 
            CALL Lh(RKDGk1,uh)
            CALL Lh(RKDGk2,uh+0.5*dt*RKDGk1)
            CALL Lh(RKDGk3,uh+0.5*dt*RKDGk2)
            CALL Lh(RKDGk4,uh+1.0*dt*RKDGk3)
            uh = uh + dt/6.0*( RKDGk1 + 2.0*RKDGk2 + 2.0*RKDGk3 + RKDGk4 )        
        END SELECT  
        !
        IF(MOD(timestep,1)==0) THEN
            IF(myrank==0) THEN
#ifdef NOLIMITER
#ifdef PARALLEL 
                tCPU3 = MPI_WTIME() 
#else
                CALL CPU_TIME(tCPU3) 
#endif 
                WRITE(*,'(a,i5,a,f14.6,a,e13.6,a,e10.3,a,f10.3)') ' n = ', timestep, '   t = ', time, '   WCT [s] = ', tCPU3-tCPU1,   '   TEU [s] = ', (REAL(tCPU3)-REAL(tCPU1))/(nElem*timestep), '  TDU [mus] = ',   (REAL(tCPU3)-REAL(tCPU1))/(nElem*timestep*(N+1)**nDim)*1e6 
#else
                WRITE(*,'(a,i,a,f,a,i,a,f10.3)') ' n = ', timestep, '   t = ', time, '   nRec = ', nRecompute, '   %troub = ', REAL(nRecompute)/REAL(nElem)*100
#endif 
            ENDIF 
        ENDIF  
        time = time + dt 
        ! For the Einstein equations, print the ADM constraints into a file 
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) || defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD)  || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)  
        !CALL RunTimeAnalyseData
#endif 
        !
        !IF(mod(timestep,1)==0) THEN 
        IF(time >= tio) THEN 
            tio = tio + dtio  
            IF(myrank==0) THEN
                PRINT *, ' Time = ', time
            ENDIF 
            CALL WriteData
        ENDIF                     
        !
    ENDDO    
#ifdef PARALLEL 
    tCPU2 = MPI_WTIME() 
#else
    CALL CPU_TIME(tCPU2)
#endif 
    TEU = timestep*nElem 
    IF(myrank==0) THEN
        PRINT *, ' Total CPU time = ', tCPU2-tCPU1 
        PRINT *, ' Time / element update = ', (tCPU2-tCPU1)/TEU 
        PRINT *, ' Time / DOF update = ', (tCPU2-tCPU1)/TEU/PRODUCT(nDOF(1:nDim))  
    ENDIF

    CALL WriteData
    
    CALL AnalyseDG 

#ifdef PARALLEL
  ! Finalize MPI 
  CALL MPI_FINALIZE(mpiErr)        
#endif

    IF(myrank==0) THEN
        PRINT *, ' ----------------------------------------- ' 
        PRINT *, '  Program terminated. Ciao.                ' 
        PRINT *, ' ----------------------------------------- ' 
    ENDIF 
    
END SUBROUTINE RKDG3D 


SUBROUTINE Lh(outduh,auh) 
    USE typesDef
    IMPLICIT NONE 
    INTEGER :: i, iElem, iFace, bcflag  
    REAL    :: auh(nVar,nDOF(1), nDOF(2), nDOF(3),1:nElem) 
    REAL    :: outduh(nVar,nDOF(1), nDOF(2), nDOF(3),1:nElem) 

        ! 
        ! The old, stupid way of doing things 
        ! 
#ifdef STUPID     
        DO iElem  = 1, nElem
           CALL ComputeNodeValues(Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),auh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
        ENDDO 
#else
        DO i = 1, nOuterElements 
            iElem = OuterElements(i) 
            CALL ComputeNodeValues(Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),auh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
        ENDDO  
#endif 
        !
#ifdef PARALLEL 
        CALL MPIExchangeDataPartI 
#endif 
        ! 
#ifndef STUPID
        DO i = 1, nInnerElements 
            iElem = InnerElements(i) 
            CALL ComputeNodeValues(Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),auh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
        ENDDO  
#endif 
        !
        ! Compute the element volume integral 
        DO iElem  = 1, nElem
            CALL ADERVolumeIntegral(duh(:,:,:,:,iElem),auh(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem))  
        ENDDO        
        ! Set the boundary conditions (could depend on space and time)  
        CALL BoundaryConditions
        ! Solve the Riemann problems for the inner faces 
        DO iFace = 1, nInnerFaces 
            bcflag = 0
            CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%paramL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%paramR,Face(iFace)%nv, bcflag )   
            CONTINUE
        ENDDO    
#ifdef PARALLEL 
        CALL MPIExchangeDataPartII  
#endif 
        ! Solve the rest of the Riemann problems for the outer faces 
        DO iFace = nInnerFaces+1, nFace 
            CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%paramL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%paramR,Face(iFace)%nv, 0 )   
        ENDDO    
        ! Compute the surface integrals of the test function multiplied with the numerical flux 
        DO iElem  = 1, nElem
            CALL ADERSurfaceIntegral(duh(:,:,:,:,iElem),FBnd(:,:,:,:,iElem))
        ENDDO       
        ! Do the element update (compute the candidate solution) 
        DO iElem  = 1, nElem
            CALL ElementUpdateRKDG(outduh(:,:,:,:,iElem),duh(:,:,:,:,iElem))
        ENDDO
        ! 
END SUBROUTINE 


SUBROUTINE ComputeNodeValues(lFh,lSh,lQbnd,lFbnd,luh,lpar) 
    USE typesDef 
    IMPLICIT NONE 
    ! Argument list
    REAL, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom
    REAL, INTENT(IN)  :: lpar(nParam,nDOF(1),nDOF(2),nDOF(3))           ! spatial degrees of freedom
    REAL, INTENT(OUT) :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3))            ! nonlinear flux tensor in each space DOF
    REAL, INTENT(OUT) :: lSh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! nonlinear source vector in each space DOF
    REAL, INTENT(OUT) :: lqbnd(nVar,nDOF(2),nDOF(3),6)                  ! space-time degrees of freedom on the boundaries 
    REAL, INTENT(OUT) :: lFbnd(nVar,nDOF(2),nDOF(3),6)                  ! nonlinear flux tensor in each space-time DOF
    ! Local variables
    INTEGER :: i,j,k,l,iVar,iDim, iter
    REAL    :: aux(d), w, ltime                                         ! auxiliary variables
    REAL    :: gradQ(nVar,d), BgradQ(nVar), src(nVar)
    REAL    :: Src_BgradQ(nVar)                                 
    REAL    :: lqx(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative qx of q
    REAL    :: lqy(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative qy of q
    REAL    :: lqz(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! spatial derivative qz of q
    REAL    :: lqt(nVar,nDOF(1),nDOF(2),nDOF(3))                        ! time derivative qt of q
    REAL    :: res                                                      ! residual
    REAL, PARAMETER :: tol = 1e-6                                       ! tolerance
    !
    lqx = 0.0
    lqy = 0.0
    lqz = 0.0 
    !
#ifndef NOFLUX
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                CALL PDEFlux(lFh(:,:,i,j,k),luh(:,i,j,k),lpar(:,i,j,k))
            ENDDO
        ENDDO
    ENDDO
#else
    lFh = 0.0     
#endif 
    ! Compute the "derivatives" (contributions of the stiffness matrix)
    ! x direction (independent from the y and z derivatives)
#ifdef NONCONSERVATIVE 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            lqx(:,:,j,k) = 1.0/dx(1)*MATMUL( luh(:,:,j,k), TRANSPOSE(dudx) )         ! compute the gradient of the solution 
        ENDDO
    ENDDO
    ! y direction (independent from the x and z derivatives) - should not be used for 1D
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lqy(:,i,:,k) = 1.0/dx(2)*MATMUL( luh(:,i,:,k), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
            ENDDO
        ENDDO
    ENDIF
    ! z direction (independent from the x and y derivatives) - should not be used for 1D and 2D
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lqz(:,i,j,:) = 1.0/dx(3)*MATMUL( luh(:,i,j,:), TRANSPOSE(dudx) )     ! compute the gradient of the solution 
            ENDDO
        ENDDO
    ENDIF
#endif 
    ! 
#if defined(NONCONSERVATIVE) || defined(SOURCE)             
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /) 
                gradQ(:,1) = lqx(:,i,j,k) 
                gradQ(:,2) = lqy(:,i,j,k) 
                gradQ(:,3) = lqz(:,i,j,k) 
                !CALL PDENCP(BgradQ,lqh(:,i,j,k,l),gradQ,lpar(:,i,j,k))    
                !lSh(:,i,j,k,l) = lSh(:,i,j,k,l) - BgradQ 
                CALL PDEFusedSrcNCP(Src_BgradQ,luh(:,i,j,k),gradQ,lpar(:,i,j,k),ltime) 
                lSh(:,i,j,k) = Src_BgradQ 
            ENDDO
        ENDDO
    ENDDO
#else
    lSh = 0.0 
#endif             
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
            lQbnd(:,j,k,1) = MATMUL( luh(:,:,j,k),   FLCoeff )    ! left
            lQbnd(:,j,k,2) = MATMUL( luh(:,:,j,k),   FRCoeff )    ! right
#ifndef NOFLUX 
            lFbnd(:,j,k,1) = MATMUL( lFh(:,1,:,j,k), FLCoeff )   ! left
            lFbnd(:,j,k,2) = MATMUL( lFh(:,1,:,j,k), FRCoeff )   ! right
#endif 
        ENDDO
    ENDDO
    ! y-direction: face 3 (left) and face 4 (right)
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1)
                lQbnd(:,i,k,3) = MATMUL( luh(:,i,:,k),   FLCoeff )    ! left
                lQbnd(:,i,k,4) = MATMUL( luh(:,i,:,k),   FRCoeff )    ! right
#ifndef NOFLUX 
                lFbnd(:,i,k,3) = MATMUL( lFh(:,2,i,:,k), FLCoeff )   ! left
                lFbnd(:,i,k,4) = MATMUL( lFh(:,2,i,:,k), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    ! z-direction: face 5 (left) and face 6 (right)
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                lQbnd(:,i,j,5) = MATMUL( luh(:,i,j,:),   FLCoeff )    ! left
                lQbnd(:,i,j,6) = MATMUL( luh(:,i,j,:),   FRCoeff )    ! right
#ifndef NOFLUX 
                lFbnd(:,i,j,5) = MATMUL( lFh(:,3,i,j,:), FLCoeff )   ! left
                lFbnd(:,i,j,6) = MATMUL( lFh(:,3,i,j,:), FRCoeff )   ! right
#endif 
            ENDDO
        ENDDO
    ENDIF
    !
    CONTINUE
    !
END SUBROUTINE      