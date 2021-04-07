PROGRAM ADERDG3D 
    USE typesDef
    IMPLICIT NONE
    ! Local variables
    INTEGER :: i,j,k,iElem,iFace,nRecompute 
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
    LOGICAL       :: dmpresult 
    !
    PRINT *, ' ----------------------------------------- ' 
    PRINT *, '   A simple introduction to 3D ADER-DG 2.0 ' 
    PRINT *, '    with a posteriori subcell FV limiter   ' 
    PRINT *, ' ----------------------------------------- ' 
    PRINT *, '      Written by Michael Dumbser           ' 
    PRINT *, '      University of Trento, Italy          ' 
    PRINT *, ' ----------------------------------------- ' 
    !
    ! We first need to compute the relevant matrices, set initial
    ! conditions and prepare some necessary stuff...  
    CALL ADERDGInit 
    
    CALL CPU_TIME(tCPU1) 
    ! Main loop in time 
    DO timestep = 1, NMAX
        IF(time >= tend) THEN
            EXIT
        ENDIF 
        ! Compute the time step size according to the CFL condition 
        CALL CalcTimeStep 
#ifdef LIMITER
        ! Limiter stuff 
        IF(N > 0) THEN
            CALL GetMinMax 
            nRecompute   = 0 
            recompute(:) = 0 
            CALL Saveolduh
        ENDIF     
#endif
        ! Process the point sources
        CALL RunPointSources     
        ! ADER predictor step 
        DO iElem  = 1, nElem
#ifdef ELASTICITY 
          CALL ADERSpaceTimePredictorLinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem),iElem)  
#else
          CALL ADERSpaceTimePredictorNonlinear(qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem),qBnd(:,:,:,:,iElem),FBnd(:,:,:,:,iElem),uh(:,:,:,:,iElem),parh(:,:,:,:,iElem))  
#endif             
        ENDDO    
        ! Compute the element volume integral 
        DO iElem  = 1, nElem
            CALL ADERVolumeIntegral(duh(:,:,:,:,iElem),qhi(:,:,:,:,iElem),Fhi(:,:,:,:,:,iElem))  
        ENDDO        
        ! Set the boundary conditions (could depend on space and time)  
        CALL BoundaryConditions 
        ! Solve the Riemann problems and compute the surface integrals 
        DO iFace = 1, nFace
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
#ifdef LIMITER
            IF(N > 0) THEN
                CALL DMP(dmpresult,uh(:,:,:,:,iElem),Limiter(iElem),0.0) 
                IF(.NOT.dmpresult) THEN
                    recompute(iElem) = 1
                    nRecompute = nRecompute + 1 
                ENDIF 
            ENDIF   
#endif     
        ENDDO    
#ifdef LIMITER     
        IF( N > 0 ) THEN
            CALL SpreadRecompute 
            CALL AllocateLimiter             
            CALL SubcellRecompute 
            CALL UpdateLimiter 
        ENDIF     
#endif    
        IF(MOD(timestep,10)==0) THEN
            PRINT *, ' n = ', timestep, ' t = ', time, 'nRec = ', nRecompute, ' %troub = ', REAL(nRecompute)/REAL(nElem)*100
            !CALL WriteData
        ENDIF  
        time = time + dt 
    ENDDO    
    CALL CPU_TIME(tCPU2)
    
    TEU = timestep*nElem 
    PRINT *, ' Total CPU time = ', tCPU2-tCPU1 
    PRINT *, ' Time / element update = ', (tCPU2-tCPU1)/TEU 
    PRINT *, ' Time / DOF update = ', (tCPU2-tCPU1)/TEU/PRODUCT(nDOF(1:nDim))  
    
    CALL WriteData
    
    CALL AnalyseDG 

    PRINT *, ' ----------------------------------------- ' 
    PRINT *, '  Program terminated. Ciao.                ' 
    PRINT *, ' ----------------------------------------- ' 
    
    
END PROGRAM ADERDG3D
    
