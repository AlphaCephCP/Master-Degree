PROGRAM ADERDG3D 
    USE typesDef
    USE F2KCLI
    IMPLICIT NONE
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif    
    ! Local variables
    INTEGER :: i,j,k,iElem,iFace,nRecompute, narg, bcflag  
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
    REAL          :: tCPU3, lx0(d), ldx(d), Int_h(nObs) 
    LOGICAL       :: dmpresult, rkdg = .FALSE.  
    CHARACTER(LEN=100)  :: argstr, rkdgstr  
    !
#ifdef PARALLEL 
    narg = COMMAND_ARGUMENT_COUNT()
    IF(narg<nDim) THEN 
        ! using hard-coded CPU distribution        
        nCPUx = (/ 1, 1, 1 /)  
    ELSE
        ! read CPU distribution from the command line 
        DO i = 1, nDim
            CALL GET_COMMAND_ARGUMENT(i,argstr)
            READ(argstr,*) nCPUx(i) 
        ENDDO 
    ENDIF 
    CALL GET_COMMAND_ARGUMENT(narg,rkdgstr) 
    SELECT CASE(TRIM(rkdgstr))
    CASE('rkdg','RKDG') 
        rkdg = .TRUE. 
    END SELECT        
    ! All CPU must initialize the MPI environment 
    CALL MPI_INIT(mpiErr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpiErr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nCPU,   mpiErr)
    ! Automatic detection of single or double precision  
    SELECT CASE(KIND(mpi_real_test))
    CASE(4)
        MPI_AUTO_REAL = MPI_REAL 
        IF(myrank==0) THEN
            PRINT *, ' Setting MPI_AUTO_REAL to MPI_REAL. ' 
        ENDIF        
    CASE(8)
        MPI_AUTO_REAL = MPI_DOUBLE_PRECISION   
        IF(myrank==0) THEN
            PRINT *, ' Setting MPI_AUTO_REAL to MPI_DOUBLE_PRECISION. ' 
        ENDIF 
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
        PRINT *, ' --------------------------------------------------- ' 
        PRINT *, '       ExaHyPE Fortran prototype version v4.0        ' 
        PRINT *, '         EU H2020 FET-HPC-1 grant No 671698          ' 
        PRINT *, ' --------------------------------------------------- ' 
        PRINT *, '        A simple introduction to 3D ADER-DG          ' 
        PRINT *, '       with a posteriori subcell FV limiter          ' 
        PRINT *, ' --------------------------------------------------- ' 
        PRINT *, '          Written by Michael Dumbser                 ' 
        PRINT *, '          University of Trento, Italy                ' 
        PRINT *, ' --------------------------------------------------- ' 
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
#ifdef RKDG
    ! For comparison purposes, we can also use an RKDG scheme... 
    ! For convenience, we can activate the RKDG mode from the command line. 
    IF(rkdg) THEN
        CALL RKDG3D 
        STOP  
    ENDIF    
#else
    PRINT *, ' In order to run the RKDG scheme, please compile the code with the -DRKDG option! ' 
    STOP 
#endif   
    !
    if(EvalSphereIntegral) then
        CALL ComputeIntegralSphere(Int_h)
    endif
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
#ifdef Z4EINSTEIN         
        ! Dicke Bertha for the Z4 system 
        !CALL DickeBertha(1) 
#endif 
        ! 
        ! Limiter stuff 
#ifdef NOLIMITER 
        nRecompute   = 0
        recompute(:) = 0 
#else        
        IF(N > 0) THEN
            CALL GetMinMax 
            nRecompute   = 0 
            recompute(:) = 0 
            CALL Saveolduh
        ENDIF     
#endif         
        ! Process the point sources
        CALL RunPointSources     
        ! 
#ifdef STUPID         
        ! The old, stupid way of doing things 
        ! 
        DO iElem  = 1, nElem
#ifdef LINEAR 
          CALL ADERSpaceTimePredictorLinear(qhi(:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem),iElem)  
#else
          CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
#endif             
        ENDDO    

#else                 
        ! The way how things should be done properly  
        ! ADER predictor step I (running only over the boundary elements) 
        DO i = 1, nOuterElements 
            iElem = OuterElements(i) 
#ifdef LINEAR 
          CALL ADERSpaceTimePredictorLinear(qhi(:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem),iElem)  
#else
          CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
#endif                         
        ENDDO      
        
#endif  ! ifdef STUPID         
        !
#if defined(Z4GRMHD) || defined(GRMHD)  
        !CALL Excision_qh  
#endif 
        ! 
#ifdef PARALLEL 
        CALL MPIExchangeDataPartI 
#endif 

#ifndef STUPID 
        ! ADER predictor step II: compute the remaining predictor, thus running only over the inner elements 
        DO i = 1, nInnerElements 
            iElem = InnerElements(i) 
#ifdef LINEAR 
          CALL ADERSpaceTimePredictorLinear(qhi(:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem),iElem)  
#else
          CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
#endif                         
        ENDDO      
#endif 
        ! 
        ! Compute the element volume integral 
        DO iElem  = 1, nElem
            CALL ADERVolumeIntegral(duh(:,:,:,:,iElem),qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),Shi(:,:,:,:,iElem))  
        ENDDO        
        ! Set the boundary conditions (could depend on space and time)  
        CALL BoundaryConditions
        ! Solve the Riemann problems for the inner faces 
        DO iFace = 1, nInnerFaces 
#ifdef GODUNOVBC             
            IF(Face(iFace)%Left == 0 ) THEN
                bcflag = 1
            ELSEIF(Face(iFace)%Right == 0 ) THEN
                bcflag = 2 
            ELSE
                bcflag = 0
            ENDIF
#else
            bcflag = 0
#endif        
            ! bcflag = 0 
            CALL ADERRiemannSolver( Face(iFace)%qL,Face(iFace)%FL,Face(iFace)%paramL,Face(iFace)%qR,Face(iFace)%FR,Face(iFace)%paramR,Face(iFace)%nv, bcflag )   
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
#ifndef NOLIMITER        
        IF( N > 0 ) THEN
            CALL SpreadRecompute 
            CALL AllocateLimiter             
            CALL SubcellRecompute 
            CALL UpdateLimiter 
            CONTINUE
        ENDIF    
#endif       
#if defined(Z4EINSTEIN)  || defined(Z4GRMHD) || defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
!        CALL Excision_uh  
!        CALL SpongeCCZ4_uh  
#elif defined(GRMHD)  
        !CALL Excision_uh 
#endif
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
    
END PROGRAM ADERDG3D
    