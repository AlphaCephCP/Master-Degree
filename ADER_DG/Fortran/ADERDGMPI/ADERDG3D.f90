PROGRAM ADERDG3D 
    USE typesDef
    USE F2KCLI
    IMPLICIT NONE
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif    
    ! Local variables
    INTEGER :: i,j,k,iElem,iFace,nRecompute, narg 
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
    LOGICAL       :: dmpresult 
    CHARACTER(LEN=100)  :: argstr 
    !
#ifdef PARALLEL

    narg = COMMAND_ARGUMENT_COUNT()
    IF(narg<nDim) THEN 
        ! using hard-coded CPU distribution        
        nCPUx = (/ 2, 1, 1 /)  
    ELSE
        ! read CPU distribution from the command line 
        DO i = 1, nDim
            CALL GET_COMMAND_ARGUMENT(i,argstr)
            READ(argstr,*) nCPUx(i) 
        ENDDO 
    ENDIF 
    ! All CPU must initialize the MPI environment 
    CALL MPI_INIT(mpiErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpiErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nCPU,   mpiErr)
    ! Automatic detection of single or double precision  
    SELECT CASE(KIND(mpi_real_test))
    CASE(4)
        MPI_AUTO_REAL = MPI_REAL 
        PRINT *, ' Setting MPI_AUTO_REAL to MPI_REAL. ' 
    CASE(8)
        MPI_AUTO_REAL = MPI_DOUBLE_PRECISION   
        PRINT *, ' Setting MPI_AUTO_REAL to MPI_DOUBLE_PRECISION. ' 
    CASE DEFAULT
        PRINT *, ' ERROR: REAL kind incompatible with MPI!. ' 
        PAUSE 
    END SELECT
#else  
    nCPU   = 1
    myrank = 0
#endif
    !
    IF(myrank==0) THEN
        PRINT *, ' --------------------------------------------- ' 
        PRINT *, '    A simple introduction to 3D ADER-DG        ' 
        PRINT *, '    with a posteriori subcell FV limiter       ' 
        PRINT *, '    ExaHyPE pre-prototype version v3.0         ' 
        PRINT *, ' --------------------------------------------- ' 
        PRINT *, '        Written by Michael Dumbser             ' 
        PRINT *, '        University of Trento, Italy            ' 
        PRINT *, ' --------------------------------------------- ' 
        IF(nCPU>=1000) THEN
            WRITE(*,'(a,i6,a)') '   The code is running on nCPU = ', nCPU, '  :-) ' 
        ELSE
            WRITE(*,'(a,i6,a)') '   The code is running on nCPU = ', nCPU, '  :-( ' 
        ENDIF 
        IF(nCPU==PRODUCT(nCPUx(:))) THEN
            WRITE(*,'(a,i4,a,i4,a,i4)') ' The CPU partition is nCPUx = ', nCPUx(1), ' nCPUy = ', nCPUx(2), ' nCPUz = ', nCPUx(3)  
        ELSE
            IF(nCPU>1) THEN
                PRINT *, ' ERROR. PRODUCT(nCPUx(:)) must be nCPU but is ', PRODUCT(nCPUx) 
                STOP
            ENDIF 
        ENDIF
        PRINT *, ' --------------------------------------------- ' 
    ENDIF 
    !
    ! We first need to compute the relevant matrices, set initial
    ! conditions and prepare some necessary stuff...  
    CALL ADERDGInit 
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
        ! Limiter stuff 
        IF(N > 0) THEN
            CALL GetMinMax 
            nRecompute   = 0 
            recompute(:) = 0 
            CALL Saveolduh
        ENDIF     
        ! Process the point sources
        CALL RunPointSources     
        ! ADER predictor step 
        DO iElem  = 1, nElem
#ifdef LINEAR 
          CALL ADERSpaceTimePredictorLinear(qhi(:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem),iElem)  
#else
          CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
#endif             
        ENDDO    
#ifdef PARALLEL 
        CALL MPIExchangeDataPartI 
#endif 
        ! Compute the element volume integral 
        DO iElem  = 1, nElem
            CALL ADERVolumeIntegral(duh(:,:,:,:,iElem),qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem))  
        ENDDO        
        ! Set the boundary conditions (could depend on space and time)  
        CALL BoundaryConditions 
        ! Solve the Riemann problems for the inner faces 
        DO iFace = 1, nInnerFaces 
            CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%paramL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%paramR,Face(iFace)%nv )   
        ENDDO    
#ifdef PARALLEL 
        CALL MPIExchangeDataPartII  
#endif 
        ! Solve the rest of the Riemann problems for the outer faces 
        DO iFace = nInnerFaces+1, nFace 
            CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%paramL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%paramR,Face(iFace)%nv )   
        ENDDO    
        ! Compute the surface integrals of the test function multiplied with the numerical flux 
        DO iElem  = 1, nElem
            CALL ADERSurfaceIntegral(duh(:,:,:,:,iElem),FBnd(:,:,:,:,iElem))
        ENDDO       
        ! Add source terms
        CALL AddPointSources    
        ! Do the element update (compute the candidate solution) 
        DO iElem  = 1, nElem
            CALL ElementUpdate(uh(:,:,:,:,iElem),duh(:,:,:,:,iElem))
            IF(N > 0) THEN
                CALL DMP(dmpresult,uh(:,:,:,:,iElem),Limiter(iElem),0.0) 
                IF(.NOT.dmpresult) THEN
                    recompute(iElem) = 1
                    nRecompute = nRecompute + 1 
                ENDIF 
            ENDIF       
        ENDDO         
        IF( N > 0 ) THEN
            CALL SpreadRecompute 
            CALL AllocateLimiter             
            CALL SubcellRecompute 
            CALL UpdateLimiter 
        ENDIF         
        IF(MOD(timestep,10)==0) THEN
            IF(myrank==0) THEN
                WRITE(*,'(a,i,a,f,a,i,a,f10.3)') ' n = ', timestep, '   t = ', time, '   nRec = ', nRecompute, '   %troub = ', REAL(nRecompute)/REAL(nElem)*100
            ENDIF 
        ENDIF  
        IF(MOD(timestep,500)==0) THEN
            CALL WriteData
        ENDIF                     
        time = time + dt 
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
    
END PROGRAM ADERDG3D
    