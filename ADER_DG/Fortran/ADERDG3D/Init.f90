SUBROUTINE ADERDGInit
    USE typesDef
    USE recipies_mod, ONLY : RG 
#ifdef TWOPUNCTURES 
	USE TwoPunctures_mod, ONLY : TwoPunctures_Setup, TwoPunctures_Run, TwoPunctures_AmendParameters
#endif         
#ifdef RNSID
    USE RNSID_mod, ONLY : RNSID_Setup, RNSID_Run, RNSID_AmendParameters 
#endif     
#ifdef TORUS    
    USE Metric_mod
#endif
    IMPLICIT NONE
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif    
    ! Local variables 
    INTEGER          :: i, j, k, ii, jj, kk, l, c, iGP, iElem, iVar, VMAX(d), count, cnt, iErr, iDim, idxs(d), iper(d) 
    INTEGER          :: Left, Right, cInnerFaces, cOuterFaces, itemp(N+1), CPUN(-1:1,-1:1,-1:1)    
    REAL             :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1), phi_xixi(N+1)
    REAL             :: phi_i(N+1), phi_j(N+1), phi_k(N+1) 
    REAL             :: ux(nVar,N+1,N+1,N+1), uy(nVar,N+1,N+1,N+1), uz(nVar,N+1,N+1,N+1) 
    REAL             :: u0(nVar), xGP(d), x0(d), subxi(N+1), aux(d), nv(d), par0(nParam)  
    REAL             :: lparambnd(nParam,6,nDOF(2),nDOF(3))    
    REAL             :: xi, dxi, xi1, xi2, rtemp(N+1), ImLambda(N+1)  
    REAL             :: TestMatrix(N+1,N+1), TestMatrix2(nSubLim,nSubLim), test(nSubLim)
    REAL             :: r1(9) 
    REAL, POINTER    :: LSQM(:,:), iLSQM(:,:), LSQrhs(:,:) 
    LOGICAL          :: dmpresult 
    REAL, PARAMETER  :: Pi = ACOS(-1.0) 
    !
    ! ----------------- Some important preliminary stuff. Do not touch! ------------------- 
    !
    dn(:) = 0 
    DO i = 1, nDim
        dn(i) = 1
    ENDDO 
    ! According to the number of space dimensions, we set the number of degrees of freedom 
    ! i = 0 is the time dimension 
    nDOF(:) = 1 
    DO i = 0, nDim
        nDOF(i) = N+1
    ENDDO 
    ! ------------------------------------------------------------------------------------- 
    !
    maxdt = 1e4     ! if the user wants, he can set an upper limit for the time step 
    !
    ! Some info about the PDE system 
    !
    EQN%Pi    = ACOS(-1.0)                                          ! Pi 
    EQN%gamma = 2.0                                                 ! ratio of specific heats for compressible Euler 
#ifdef SRMHD
    EQN%gamma = 4./3. 
#endif 
#ifdef GRMHD
    !GRMHD-BondiAccretion
    !EQN%gamma = 4./3. 
    !
    !!AlfvenWave SRMHD limit
    EQN%gamma = 5./3. 
#ifdef TORUS
     EQN%gamma = 4./3.
#endif
#endif 
#ifdef Z4GRMHD
    EQN%gamma = 4./3. 
#endif 
#ifdef CCZ4GRMHD
    EQN%gamma = 4./3. 
#endif 
    ! 
    !
    ! Old parameters for the Z4 system 
    !
    !EQN%CCZ4k1  = 0.0
    !EQN%CCZ4k2  = 0.0 
    !EQN%CCZ4k3  = 0.0 
    !EQN%CCZ4eta = 0.0 
    !EQN%CCZ4f   = 0.0 
    !EQN%CCZ4g   = 1.0 
    !EQN%CCZ4xi  = 0.0 
    !EQN%CCZ4e   = 1.0 
    !EQN%CCZ4c   = 0.0 
    !EQN%CCZ4mu  = 0.0 
    !EQN%CCZ4ds  = 1.0 
    !EQN%CCZ4sk  = 0.0
    !EQN%CCZ4bs   = 0.0     
    !EQN%CCZ4LapseType   = 1 
    !EQN%EinsteinAutoAux = 0  
    ! 
    ! Gamma driver for the black hole simulation 
    !EQN%CCZ4k1   = 0.1      ! CCZ4 damping parameter k1 
    !EQN%CCZ4k2   = 0.0      ! CCZ4 damping parameter k2  
    !EQN%CCZ4k3   = 1.0      ! CCZ4 damping parameter k3 
    !EQN%CCZ4eta  = 0.0      ! CCZ4 damping parameter for the PDE of b^i in the gamma driver 
    !EQN%CCZ4itau = 0.0      ! inverse relaxation time for the relaxed approach in order to enforce the unit determinant of the conformal metric and the trace-free matrix A_ij
    !EQN%CCZ4f    = 1.0      ! set f=0.75 (not hyperbolic!!) or f=1.0 for the gamma driver. Typical BSSNOK value: f=0.75. Set f=0 to avoid shift evolution, i.e. to have d/dt beta^i=0.  
    !EQN%CCZ4g    = 0.0      ! not used at the moment: reserved for the generalized harmonic shift   
    !EQN%CCZ4xi   = 1.0      ! set to zero to switch off the gamma driver and to avoid shift evolution (i.e. to have d/dt b^i = 0) 
    !EQN%CCZ4e    = 1.0      ! cleaning speed e>=1 for the Hamiltonian constraint. Typical value for CCZ4 is e=1. However, e>1 gives better constraints and better hyperbolicity!  
    !EQN%CCZ4c    = 1.0      ! set c=0 to remove the algebraic source terms of the type -2*Theta 
    !EQN%CCZ4mu   = 0.2      ! mu=1 adds second order ordering constraints. Has shown to be important for the gamma driver, but too large values can also become a problem... 
    !EQN%CCZ4ds   = 1.0      ! set this value always to ds=1, unless you know what you are doing. It allows to increase the cleaning speed for the momentum constraints, but in CCZ4 this does not seem to do what it should do...  
    !EQN%CCZ4sk   = 1.0      ! setting sk=0 removes the contribution of the shift also in the auxiliary variables. If you want to evolve the shift, set sk=1. 
    !EQN%CCZ4bs   = 1.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
    !EQN%CCZ4LapseType   = 1 ! LapseType = 0 is harmonic lapse, LapseType = 1 is 1+log slicing.   
    !EQN%EinsteinAutoAux = 0  
    ! Parameters for the Gauge wave 
    EQN%CCZ4k1  = 0.0 
    EQN%CCZ4k2  = 0.0 
    EQN%CCZ4k3  = 0.0 
    EQN%CCZ4eta = 0.0 
    EQN%CCZ4f   = 0.0 
    EQN%CCZ4g   = 0.0 
    EQN%CCZ4xi  = 0.0 
    EQN%CCZ4e   = 2.0 
    EQN%CCZ4c   = 0.0 
    EQN%CCZ4mu  = 0.0 
    EQN%CCZ4ds  = 1.0 
    EQN%CCZ4sk  = 0.0
    EQN%CCZ4bs   = 0.0      ! set bs=1 if you want to activate the shift convection for beta, b and B (standard CCZ4 formulation). set it to bs=0 to switch off shift convection for those quantities 
    EQN%CCZ4LapseType   = 0 ! harmonic lapse 
    EQN%EinsteinAutoAux = 0 
    ! 
    dtio = 0.1              
    !
    ! Here, you need to define the computational domain and the number of cells. 
    ! This is typically read from a parameter file 
    !
    !xL = (/ -1.0, -1.0, -1.0 /)                                      ! lower-left corner of the domain 
    !xR = (/ +1.0, +1.0, +1.0 /)                                      ! upper right corner of the domain 
    !xL = 20*(/ -1.0, -1.0, -1.0 /)                                   ! lower-left corner of the domain 
    !xR = 20*(/ +1.0, +1.0, +1.0 /)                                   ! upper right corner of the domain 
    !xL = (/ -5.0, -5.0, -1.0 /)                                      ! lower-left corner of the domain  (for ShuVortex2D) 
    !xR = (/ 15.0, 15.0, +1.0 /)                                      ! upper right corner of the domain (for ShuVortex2D)  
    !xL = (/  0.0,     0.0,    0.0 /)                                 ! lower-left corner of the domain  (for ShuVortex2D) 
    !xR = (/  4000.0,  2000.0, 2000.0 /)                              ! upper right corner of the domain (for ShuVortex2D)  
    !xL = (/  0.0, 0.0,    0.0 /)                                     ! lower-left corner of the domain 
    !xR = (/ +1.0, 1./3., +1./3. /)                                   ! upper right corner of the domain 
    !xL = (/ 1.0, 1.0, 1.0 /) 
    !xR = (/ 2.0, 2.0, 2.0 /) 
    !xL = (/ 0.0, 0.0, 0.0 /)                                         ! shu vortex 2d 
    !xR = (/ 10., 10., 10. /) 
    !xL = (/ 0.0,    0.0,    0.0    /)                                 ! grmhd alfven wave
    !xR = (/ 2.0*Pi, 2.0*Pi, 2.0*Pi /) 
    !xL = (/ -0.5, -0.1, -0.1 /) 
    !xR = (/ +0.5, +0.1, +0.1 /) 
    !xR = (/ 11.0, 11.0, 11.0 /) 
    !xL = 0.5*(/ -1.0, -0.1, -0.1 /) 
    !xR = 0.5*(/  1.0,  0.1,  0.1 /) 
    xL = 0.5*(/ -1.0, -1.0, -1.0 /) 
    xR = 0.5*(/  1.0,  1.0,  1.0 /) 
#ifndef GRMHD    
    ExcisionRadius = -1.0    
#endif    
    !xL = (/ -1.0, -1.0, -0.1 /) 
    !xR = (/  1.0,  1.0,  0.1 /) 
    !xL = (/ -20., -20., -20. /) 
    !xR = (/ +20., +20., +20. /) 
    !xL = (/ 1.0, 0.5, -1.0 /)                                       ! lower-left corner of the domain 
    !xR = (/ 3.0, 2.5, +1.0 /)                                       ! upper right corner of the domain 
    !xL = (/ 1.0, 1.0, -1.0 /)                                       ! lower-left corner of the domain 
    !xR = (/ 1.1, 1.1, +1.0 /)                                       ! upper right corner of the domain 
    !xL = (/ -100., -100., -100. /) 
    !xR = (/ +100., +100., +100. /) 
    !xL = (/ -10.0, -10.0, -10.0 /)                                  ! lower-left corner of the domain 
    !xR = (/ +10.0, +10.0, +10.0 /)                                  ! upper right corner of the domain 
    

    ! 
    IMAX = 5     ! 81                                               ! Number of elements in x,y,z direction 
    JMAX = 5     ! 27
    KMAX = 5
    !
#ifdef TORUS
#ifdef Spherical
    xL = (/ 2.0 , 0.52359877559829887307710723054658, -0.5 /)                                  ! lower-left corner of the domain 
    xR = (/ 18.0 , 2.6179938779914943653855361527329 , 0.5 /)                                  ! upper right corner of the domain 
    IMAX = 20      ! 81                                               ! Number of elements in x,y,z direction 
    JMAX = 20      ! 27
    KMAX = 1  
#else
    xL = (/ -18.0 ,2, -8/)                                  ! lower-left corner of the domain 
    xR = (/ 18.0 , 18 , 8 /)                                  ! upper right corner of the domain 
    IMAX = 20      ! 81                                               ! Number of elements in x,y,z direction 
    JMAX = 10      ! 27
    KMAX = 10  
#endif
#endif
!0.52359877559829887307710723054658                              ! xmin  PI/6
!2.6179938779914943653855361527329                                 ! xmax PI 5/6

    !IMAX = 100     
    !JMAX = 1 
    !KMAX = 1 
    VMAX = (/ IMAX, JMAX, KMAX /)                                    ! Vector of the number of elements in each space dimension 
    dx = (xR-xL)/VMAX                                                ! Mesh spacing 
    NMAX = 10000000                                                  ! Max. number of time steps 
    timestep = 0                                                     ! initial time step number 
    time = 0.                                                        ! initial time 
    tend = 0.2                                                      ! final time 
    !ICType = 'Sod'                                                  ! type of initial condition 'ShuVortex2D'      
    !Basefile = 'Sod'                                                ! Base filename for writing results 
    tend = 1.0                                                      ! final time 
    ICType = 'ShuVortex2D'                                          ! type of initial condition 'ShuVortex2D'      
    Basefile = 'ShuVortex2D_ADERDG'                                 ! Base filename for writing results 
#ifdef RKDG
    Basefile = 'ShuVortex2D_RKDG'                                   ! Base filename for writing results 
#endif 
    !tend = 0.50 
    !ICType = 'EP2D' 
    !Basefile = 'EP2D' 
    !tend = 0.5                                                      ! final time 
    !ICType = 'LOH1'                                                 ! type of initial condition 'ShuVortex2D'      
    !Basefile = 'DebugLOH3D'                                         ! Base filename for writing results 
    !tend = 0.5
    !ICType = 'RPGaussian'
    !Basefile = 'pwave' 
#ifdef Z4EINSTEIN  
    tend = 1.0                                                        ! final time 
    !ICType = 'Z4GaugeWave'                                           ! type of initial condition 'ShuVortex2D'      
    !Basefile = 'Z4GaugeWave'                                         ! Base filename for writing results 
    ICType = 'Z4LinearWave'                                           ! type of initial condition 'ShuVortex2D'      
    Basefile = 'Z4LinearWave'                                         ! Base filename for writing results 
    !ICType = 'Z4Kerr2D'                                              ! type of initial condition 'ShuVortex2D'      
    !Basefile = 'Z4Kerr2D'                                            ! Base filename for writing results 
#endif     
#ifdef SRMHD  
    tend = 2e-2  ! 2.61803398874990	                                  ! final time 
    ICType = 'RMHDAlfvenWave'                                         ! type of initial condition 'ShuVortex2D'      
    Basefile = 'RMHDAlfvenWave'                                       ! Base filename for writing results 
#endif    
#ifdef GRMHD
    !tend = 10.0 
    !ICType = 'GRMHDAccretion'
    !Basefile = 'GRMHDAccretion' 
    !!
    !ICType2 = 'GRMHD-BondiAccretion'  ! special case in Schwarzschild metric
    !
    ! AlfvenWave: 
    !tend = 22.2733119873
    tend = 1.0 
    dtio = 0.1 
    ICType = 'GRMHDAlfven'
    Basefile = 'GRMHDAlfven' 
    !
    !ICType2 = 'GRMHD-Alfven'  ! special case in Schwarzschild metric
#ifdef TORUS
    tend = 10.0 
    ICType = 'GRMHDTorus'
    Basefile = 'GRMHDTorus' 
    !
    ICType2 = 'GRMHDTorus'  ! special case in Schwarzschild metric
    
#endif
#endif  
#ifdef Z4GRMHD
    tend = 1.0  ! e-6 
    ICType = 'Z4GRMHDAccretion'
    Basefile = 'Z4GRMHDAccretion' 
#endif  
#ifdef CCZ4GRMHD
    tend = 1.0  ! e-6 
    !ICType = 'CCZ4GRMHDAccretion'
    !Basefile = 'CCZ4GRMHDAccretion' 
    !ICType = 'CCZ4Minkowski'
    !Basefile = 'CCZ4Minkowski' 
    !ICType = 'CCZ4GasCloud'
    !Basefile = 'CCZ4GasCloud' 
    ICTYPE   = 'CCZ4RNSID' 
    Basefile = 'CCZ4RNSID' 
#endif  
#ifdef CCZ4EINSTEIN 
    ! 
    time = 0.0 
    tend = 0.1  
    ICType   = 'CCZ4GaugeWave'
    Basefile = 'CCZ4GaugeWave'
    !ICType   = 'CCZ4LinearWave'
    !Basefile = 'CCZ4LinearWave-P5-4' 
    !ICType   = 'CCZ4GowdyWave'
    !Basefile = 'CCZ4GowdyWave' 
    !ICType   = 'CCZ4Minkowski'
    !Basefile = 'CCZ4Minkowski' 
    !ICType   = 'CCZ4Puncture'
    !Basefile = 'CCZ4Puncture' 
    !ICType   = 'CCZ4TwoPunctures'
    !Basefile = 'CCZ4TwoPunctures' 
    !ICType   = 'CCZ4GRMHDAccretion'
    !Basefile = 'CCZ4GRMHDAccretion' 
    !ICType    = 'CCZ4Kerr3D'
    !Basefile  = 'CCZ4Kerr3D'
    !ICType    = 'CCZ4Kerr2D'
    !Basefile  = 'CCZ4Kerr2D'
#endif  
#ifdef BSSZ4EINSTEIN 
    tend = 1000.0  ! e-6 
    !ICType   = 'CCZ4GaugeWave'
    !Basefile = 'BSSZ4GaugeWave'  
    ICType = 'CCZ4Minkowski'
    Basefile = 'CCZ4Minkowski' 
#endif  
#ifdef BSSZ4GRMHD 
    tend = 1.0  ! e-6 
    ICType = 'CCZ4GRMHDAccretion'
    Basefile = 'BSSZ4GRMHDAccretion' 
#endif  
#ifdef GPR3D
    ICType = 'ShockNavierStokes'
    Basefile = 'GPRShock'     
    !Basefile = 'GPRShuVortex'     
    EQN%gamma = 1.4
    EQN%rho0  = 1.0 
    EQN%cs    = 5.0 
    EQN%alpha = 5.0 
    EQN%tau1  = 1e-2 
    EQN%tau2  = 0.3111111112e-2  
    EQN%mu    = 1./6.*EQN%rho0*EQN%cs**2*EQN%tau1 
    EQN%kappa = EQN%alpha**2*EQN%tau2 
    EQN%p0    = 0.0 
    EQN%cv    = 1.0 
    tend = 0.1 
#endif 
    ! 
    SELECT CASE(ICType)
    CASE('CCZ4TwoPunctures','Z4TwoPunctures') 
#ifdef TWOPUNCTURES  
        CALL TwoPunctures_Setup()	
	    ! After Setup, you can change the parameters of TwoPunctures, such as Black Hole
	    ! mass, etc. They are stored as C++ class members, so there is no way of changing
	    ! them from Fortran. However, you can do so if you call the following function
	    ! and changes it's implementation
	    !CALL TwoPunctures_AmendParameters()	
	    PRINT *, ' Starting the TwoPunctures algorithm... '
	    CALL TwoPunctures_Run() 
	    PRINT *, ' Done. '
#ifdef PARALLEL
        CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr) 
#endif
#else
        PRINT *, ' Twopunctures not available. Please compile with -DTWOPUNCTURES. ' 
        STOP 
#endif         
    CASE('CCZ4RNSID') 
#ifdef RNSID 
        CALL RNSID_Setup()	 
        CALL RNSID_AmendParameters()  
        PRINT *, ' Starting the RNSID algorithm... '
        CALL RNSID_Run() 
	    PRINT *, ' Done. '        
#ifdef PARALLEL
        CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr) 
#endif
        ! 
#else
        PRINT *, ' RNSID not available. Please compile with -DRNSID. ' 
        STOP 
#endif 
    CASE('GRMHDTorus')
#ifdef TORUS
        CALL init_disk_isentropic_par(METRIC_KSS)
#endif
    ! 
    END SELECT 
    !
    tio  = time + dtio     
    ! 
    !Periodic(:) = (/ .FALSE., .TRUE., .FALSE. /)                          ! no periodic BC 
    Periodic(:) = (/ .TRUE., .FALSE., .TRUE. /)                           ! periodic BC 
#ifdef GRMHD

    Periodic(:) = (/ .TRUE., .TRUE., .TRUE. /)                           ! periodic BC 
#ifdef TORUS
#ifdef Spherical
    Periodic(:) = (/ .FALSE., .FALSE., .TRUE. /) 
#else
    Periodic(:) = (/ .FALSE., .FALSE., .FALSE. /) 
#endif
#endif 
#endif
    AnalyseType = 2                                                  ! comparison with exact solution at the end of the simulation 
    !
    nElem = IMAX*JMAX*KMAX                                           ! Number of elements 
    nNode = (IMAX+1)*(JMAX+1)*(KMAX+1)                               ! number of nodes 
    ALLOCATE( x(d, nNode) )                                          ! Allocate the nodes 
    ALLOCATE( idxn(IMAX+dn(1),JMAX+dn(2),KMAX+dn(3))  )                                              
    ! Define the node coordinates and the node numbers                                 
    count = 0 
    DO k = 1, KMAX+dn(3)  
        DO j = 1, JMAX+dn(2) 
            DO i = 1, IMAX+dn(1) 
                count = count + 1 
                x(:,count) = xL(:) + (/ i-1, j-1, k-1/)*dx(:) 
                idxn(i,j,k) = count 
            ENDDO
        ENDDO
    ENDDO
    ! define the connectivity between the elements and the nodes. You can do this more compactly via a loop, but I prefer the select case, which makes the explicit 
    ! construction of each element clearer 
    ALLOCATE( tri(nVtx, nElem)     ) 
    ALLOCATE( idxe(IMAX,JMAX,KMAX) ) 
    count = 0  
    DO k = 1, KMAX  
        DO j = 1, JMAX
            DO i = 1, IMAX 
                count = count + 1
                idxe(i,j,k) = count 
                SELECT CASE(nDim)
                CASE(1)                    
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k) /) 
                CASE(2)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k) /) 
                CASE(3)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i,j+1,k+1), idxn(i+1,j+1,k+1) /) 
                END SELECT                
            ENDDO
        ENDDO
    ENDDO
    ! compute the Voronoi neighbors of a cell
    ALLOCATE( neighbor(-1:1,-1:1,-1:1,nElem) ) 
    DO i = 1, nElem
        neighbor(:,:,:,i) = i  
    ENDDO    
    DO k = 1, KMAX  
     DO j = 1, JMAX
      DO i = 1, IMAX 
        count = idxe(i,j,k) 
        ! 
        DO ii = -dn(1), dn(1) 
         DO jj = -dn(2), dn(2) 
          DO kk = -dn(3), dn(3)
              IF( (MIN( i+ii, j+jj, k+kk ).GE.1) .AND. (MAX( i+ii-IMAX, j+jj-JMAX, k+kk-KMAX ).LE.0) ) THEN 
                neighbor(ii,jj,kk,count) = idxe(i+ii,j+jj,k+kk) 
              ENDIF
              idxs = (/ i, j, k /) + (/ ii, jj, kk /) 
              DO iDim = 1, d 
                IF(periodic(iDim)) THEN
                  IF(idxs(iDim).LT.1) THEN
                    idxs(iDim) = MIN(VMAX(iDim),MAX(1,VMAX(iDim) + idxs(iDim))) 
                  ENDIF 
                  IF(idxs(iDim).GT.VMAX(iDim)) THEN
                    idxs(iDim) = MIN(VMAX(iDim),MAX(1,idxs(iDim) - VMAX(iDim))) 
                  ENDIF 
                ENDIF 
                IF(ANY(idxs.LT.1) ) THEN                     
                    CYCLE
                ENDIF
                IF(ANY((idxs-VMAX).GT.0)) THEN 
                    CYCLE 
                ENDIF
                IF(idxe(idxs(1),idxs(2),idxs(3)).NE.0) THEN
                    neighbor(ii,jj,kk,count) = idxe(idxs(1),idxs(2),idxs(3))  
                ELSE
                    neighbor(ii,jj,kk,count) = count
                ENDIF
               ENDDO
               !
          ENDDO
         ENDDO
        ENDDO
        !  
      ENDDO
     ENDDO     
    ENDDO
    !
    ! Compute the mapping from face to 3D neighbor index 
    !
    Face2Neigh(:,1) = (/ -1, 0, 0 /)   
    Face2Neigh(:,2) = (/ +1, 0, 0 /)   
    Face2Neigh(:,3) = (/  0,-1, 0 /)   
    Face2Neigh(:,4) = (/  0,+1, 0 /)   
    Face2Neigh(:,5) = (/  0, 0,-1 /)   
    Face2Neigh(:,6) = (/  0, 0,+1 /)   
    !
    ! We now need to define our basis functions. This is done by choosing a set of distinct 1D nodes in the interval [0,1] 
    ! The basis functions will be the Lagrange interpolation polynomials running through these nodes 
    CALL gauleg(0.,1.,xiGPN,wGPN,N+1)
    xin = xiGPN  ! WE make the following choice: the basis functions run through the Gauss-Legendre nodes (=> really orthogonal basis) 
    CALL gauleg(0.,1.,xiGPMm1,wGPMm1,N) 
    !
    CALL InitPointSources 
    !
    ALLOCATE( CPUDistribution(nElem) ) 
    CALL ComputeCPUDistribution(CPUDistribution,VMAX,nCPUx)
    CALL MPIExtractMesh 
    ! 
    DEALLOCATE( idxn )
    !
    ALLOCATE(  uh(nVar, nDOF(1), nDOF(2), nDOF(3),         1:nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE( duh(nVar, nDOF(1), nDOF(2), nDOF(3),         1:nElem) )         ! Allocate the degrees of freedom for the update 
    ALLOCATE( qhi(nVar, nDOF(1), nDOF(2), nDOF(3),         1:nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( Fhi(nVar, d, nDOF(1), nDOF(2), nDOF(3),      1:nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( Shi(nVar, nDOF(1), nDOF(2), nDOF(3),         1:nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( parh(nParam, nDOF(1), nDOF(2), nDOF(3),      1:nElem) )         ! Allocate the degrees of freedom for the material parameters 
    ALLOCATE( qbnd(nVar, nDOF(2), nDOF(3),      6, -nMPIElem:nElem) )         ! Allocate the time-averaged boundary-extrapolated values for Q 
    ALLOCATE( Fbnd(nVar, nDOF(2), nDOF(3),      6, -nMPIElem:nElem) )         ! Allocate the time-averaged boundary-extrapolated values for the normal flux F * n     
    ALLOCATE( parbnd(nParam, nDOF(2), nDOF(3),  6, -nMPIElem:nElem) )         ! Allocate the boundary-extrapolated material parameters 
#ifdef RKDG
    ALLOCATE(  RKDGk1(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE(  RKDGk2(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE(  RKDGk3(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE(  RKDGk4(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
    IF(N>3) THEN
        ALLOCATE(  RKDGk5(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
        ALLOCATE(  RKDGk6(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
        ALLOCATE(  RKDGk7(nVar, nDOF(1), nDOF(2), nDOF(3),     1:nElem) )         ! Allocate the spatial degrees of freedom 
    ENDIF    
#endif 
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
    ALLOCATE(  ADM(4, nDOF(1), nDOF(2), nDOF(3),           1:nElem) )         ! Allocate the ADM degrees of freedom 
#elif defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    ALLOCATE(  ADM(6, nDOF(1), nDOF(2), nDOF(3),           1:nElem) )         ! Allocate the ADM degrees of freedom 
#endif     
    !
    ! define the connectivity between the faces and the elements 
    ! count how many faces we have, taking into account the MPI stuff  
    DO iDim = 1, d
        IF(Periodic(iDim)) THEN
            iper(iDim) = 0
        ELSE
            iper(iDim) = 1 
        ENDIF 
    ENDDO
!    nFace = KMAX*JMAX*(IMAX+iper(1))
!    IF(nDim>=2) THEN 
!        nFace = nFace + KMAX*(JMAX+iper(2))*IMAX
!    ENDIF
!    IF(nDim>=3) THEN 
!        nFace = nFace + (KMAX+iper(3))*JMAX*IMAX 
!    ENDIF
    !
    nInnerElements = 0 
    nOuterElements = 0 
    DO i = 1, nElem        
        DO kk = -1, 1
         DO jj = -1, 1 
          DO ii = -1, 1 
            CPUN(ii,jj,kk) = LocalCPUDistribution( neighbor(ii,jj,kk,i) ) 
          ENDDO
         ENDDO
        ENDDO
        ! 
        IF( MAXVAL(ABS(CPUN-myrank)).GT.0 ) THEN 
            nOuterElements = nOuterElements + 1 
        ELSE
            nInnerElements = nInnerElements + 1 
        ENDIF        
    ENDDO    
    ALLOCATE( OuterElements(nOuterElements) ) 
    ALLOCATE( InnerElements(nInnerElements) ) 
    nInnerElements = 0 
    nOuterElements = 0 
    DO i = 1, nElem        
        DO kk = -1, 1
         DO jj = -1, 1 
          DO ii = -1, 1 
            CPUN(ii,jj,kk) = LocalCPUDistribution( neighbor(ii,jj,kk,i) )  
          ENDDO
         ENDDO
        ENDDO
        IF( MAXVAL(ABS(CPUN-myrank)).GT.0 ) THEN
            nOuterElements = nOuterElements + 1 
            OuterElements(nOuterElements) = i 
        ELSE
            nInnerElements = nInnerElements + 1 
            InnerElements(nInnerElements) = i 
        ENDIF        
    ENDDO 
    ! 
    PRINT *, ' MPI rank ', myrank, ' nInner = ', nInnerElements, ' nOuter = ', nOuterElements  
    !
    nFace = 0 
    nInnerFaces = 0 
    nOuterFaces = 0 
    ! x faces 
    DO k = 1, KMAX
     DO j = 1, JMAX
      DO i = 1, IMAX+iper(1)  
          IF(i.EQ.1) THEN
            IF(Periodic(1)) THEN
                Left  = idxe(IMAX,j,k) 
                Right = idxe(i,j,k)
                IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                    nFace = nFace + 1 
                    IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                        nInnerFaces = nInnerFaces + 1 
                    ELSE
                        nOuterFaces = nOuterFaces + 1
                    ENDIF 
                ENDIF  
            ELSE
                Left  = 0 
                Right = idxe(i,j,k)
                IF( CPUDistribution(Right)==myrank ) THEN
                    nFace = nFace + 1 
                    nInnerFaces = nInnerFaces + 1 
                ENDIF  
            ENDIF 
          ELSEIF(i.EQ.IMAX+1) THEN 
                Left  = idxe(i-1,j,k) 
                Right = 0 
                IF( CPUDistribution(Left)==myrank ) THEN
                    nFace = nFace + 1 
                    nInnerFaces = nInnerFaces + 1 
                ENDIF  
          ELSE              
                Left  = idxe(i-1,j,k) 
                Right = idxe(i,j,k) 
                IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                    nFace = nFace + 1 
                    IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                        nInnerFaces = nInnerFaces + 1 
                    ELSE
                        nOuterFaces = nOuterFaces + 1
                    ENDIF 
                ENDIF  
          ENDIF 
      ENDDO
     ENDDO
    ENDDO 
    ! y faces 
    IF(nDim>=2) THEN
        DO k = 1, KMAX
         DO j = 1, JMAX+iper(2) 
          DO i = 1, IMAX  
              IF(j.EQ.1) THEN
                IF(Periodic(2)) THEN
                    Left  = idxe(i,JMAX,k)
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                            nInnerFaces = nInnerFaces + 1 
                        ELSE
                            nOuterFaces = nOuterFaces + 1
                        ENDIF 
                    ENDIF  
                ELSE
                    Left  = 0
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        nInnerFaces = nInnerFaces + 1 
                    ENDIF  
                ENDIF
              ELSEIF(j.EQ.JMAX+1) THEN 
                    Left  = idxe(i,j-1,k) 
                    Right = 0 
                    IF( CPUDistribution(Left)==myrank ) THEN
                        nFace = nFace + 1
                        nInnerFaces = nInnerFaces + 1  
                    ENDIF  
              ELSE              
                    Left  = idxe(i,j-1,k) 
                    Right = idxe(i,j,k) 
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                            nInnerFaces = nInnerFaces + 1 
                        ELSE
                            nOuterFaces = nOuterFaces + 1
                        ENDIF 
                    ENDIF  
              ENDIF 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    ! z faces 
    IF(nDim>=3) THEN
        DO k = 1, KMAX+iper(3) 
         DO j = 1, JMAX 
          DO i = 1, IMAX  
              IF(k.EQ.1) THEN
                IF(Periodic(3)) THEN
                    Left  = idxe(i,j,KMAX) 
                    Right = idxe(i,j,k) 
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                            nInnerFaces = nInnerFaces + 1 
                        ELSE
                            nOuterFaces = nOuterFaces + 1
                        ENDIF 
                    ENDIF  
                ELSE
                    Left  = 0 
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        nInnerFaces = nInnerFaces + 1 
                    ENDIF  
                ENDIF 
              ELSEIF(k.EQ.KMAX+1) THEN 
                    Left  = idxe(i,j,k-1) 
                    Right = 0 
                    IF( CPUDistribution(Left)==myrank ) THEN
                        nFace = nFace + 1 
                        nInnerFaces = nInnerFaces + 1 
                    ENDIF  
              ELSE              
                    Left  = idxe(i,j,k-1) 
                    Right = idxe(i,j,k) 
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        nFace = nFace + 1 
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN
                            nInnerFaces = nInnerFaces + 1 
                        ELSE
                            nOuterFaces = nOuterFaces + 1
                        ENDIF 
                    ENDIF  
              ENDIF 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    !
    PRINT *, ' myrank = ', myrank, ' nFace = ', nFace 
    !
    ALLOCATE( Face(nFace) ) 
    ! x faces
    cInnerFaces = 0 
    cOuterFaces = nInnerFaces  
    DO k = 1, KMAX
     DO j = 1, JMAX
      DO i = 1, IMAX+iper(1)  
          IF(i.EQ.1) THEN
            IF(Periodic(1)) THEN
                Left  = idxe(IMAX,j,k) 
                Right = idxe(i,j,k) 
                IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                    IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        cOuterFaces = cOuterFaces + 1 
                        count = cOuterFaces 
                    ENDIF 
                ELSE
                    CYCLE  
                ENDIF  
                Face(count)%Left  = Global2LocalElem(idxe(IMAX,j,k)) 
                Face(count)%Right = Global2LocalElem(idxe(i,j,k))  
                Face(count)%qL => qBnd(:,:,:,2,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,:,:,1,Face(count)%Right) 
                Face(count)%FL => FBnd(:,:,:,2,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,:,:,1,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,:,:,2,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,:,:,1,Face(count)%Right) 
            ELSE
                Left  = 0 
                Right = idxe(i,j,k)
                IF( CPUDistribution(Right)==myrank ) THEN
                    cInnerFaces = cInnerFaces + 1 
                    count = cInnerFaces 
                ELSE
                    CYCLE  
                ENDIF  
                Face(count)%Left  = 0 
                Face(count)%Right = Global2LocalElem(idxe(i,j,k)) 
                ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%qR => qBnd(:,:,:,1,Face(count)%Right) 
                Face(count)%FR => FBnd(:,:,:,1,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,:,:,1,Face(count)%Right) 
                Face(count)%paramR => parbnd(:,:,:,1,Face(count)%Right) 
            ENDIF 
          ELSEIF(i.EQ.IMAX+1) THEN 
                Left  = idxe(i-1,j,k) 
                Right = 0 
                IF( CPUDistribution(Left)==myrank ) THEN
                    cInnerFaces = cInnerFaces + 1 
                    count = cInnerFaces 
                ELSE
                    CYCLE  
                ENDIF  
                Face(count)%Left  = Global2LocalElem(idxe(i-1,j,k))  
                Face(count)%Right = 0 
                Face(count)%qL => qBnd(:,:,:,2,Face(count)%Left ) 
                Face(count)%FL => FBnd(:,:,:,2,Face(count)%Left ) 
                ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%paramL => parbnd(:,:,:,2,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,:,:,2,Face(count)%Left ) 
          ELSE              
            Left  = idxe(i-1,j,k) 
            Right = idxe(i,j,k) 
            IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                    cInnerFaces = cInnerFaces + 1 
                    count = cInnerFaces 
                ELSE
                    cOuterFaces = cOuterFaces + 1 
                    count = cOuterFaces 
                ENDIF 
            ELSE
                CYCLE  
            ENDIF  
            Face(count)%Left  = Global2LocalElem(idxe(i-1,j,k)) 
            Face(count)%Right = Global2LocalElem(idxe(i,j,k)) 
            Face(count)%qL => qBnd(:,:,:,2,Face(count)%Left ) 
            Face(count)%qR => qBnd(:,:,:,1,Face(count)%Right) 
            Face(count)%FL => FBnd(:,:,:,2,Face(count)%Left ) 
            Face(count)%FR => FBnd(:,:,:,1,Face(count)%Right) 
            Face(count)%paramL => parbnd(:,:,:,2,Face(count)%Left ) 
            Face(count)%paramR => parbnd(:,:,:,1,Face(count)%Right) 
          ENDIF 
          Face(count)%nv = (/ 1., 0., 0. /)                               ! set face normal vector
          Face(count)%x0 = xL + (/ REAL(i-1), REAL(j-1), REAL(k-1) /)*dx  ! set face lower left corner 
      ENDDO
     ENDDO
    ENDDO 
    ! y faces 
    IF(nDim>=2) THEN
        DO k = 1, KMAX
         DO j = 1, JMAX+iper(2) 
          DO i = 1, IMAX  
              IF(j.EQ.1) THEN
                IF(Periodic(2)) THEN
                    Left  = idxe(i,JMAX,k)
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                            cInnerFaces = cInnerFaces + 1 
                            count = cInnerFaces 
                        ELSE
                            cOuterFaces = cOuterFaces + 1 
                            count = cOuterFaces 
                        ENDIF 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = Global2LocalElem(idxe(i,JMAX,k))
                    Face(count)%Right = Global2LocalElem(idxe(i,j,k))
                    Face(count)%qL => qBnd(:,:,:,4,Face(count)%Left ) 
                    Face(count)%qR => qBnd(:,:,:,3,Face(count)%Right) 
                    Face(count)%FL => FBnd(:,:,:,4,Face(count)%Left ) 
                    Face(count)%FR => FBnd(:,:,:,3,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,:,:,4,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,:,:,3,Face(count)%Right) 
                ELSE
                    Left  = 0
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Right)==myrank ) THEN
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = 0
                    Face(count)%Right = Global2LocalElem(idxe(i,j,k))
                    ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%qR => qBnd(:,:,:,3,Face(count)%Right) 
                    Face(count)%FR => FBnd(:,:,:,3,Face(count)%Right)
                    Face(count)%paramL => parbnd(:,:,:,3,Face(count)%Right) 
                    Face(count)%paramR => parbnd(:,:,:,3,Face(count)%Right) 
                ENDIF
              ELSEIF(j.EQ.JMAX+1) THEN 
                    Left  = idxe(i,j-1,k) 
                    Right = 0 
                    IF( CPUDistribution(Left)==myrank ) THEN
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = Global2LocalElem(idxe(i,j-1,k)) 
                    Face(count)%Right = 0 
                    Face(count)%qL => qBnd(:,:,:,4,Face(count)%Left ) 
                    Face(count)%FL => FBnd(:,:,:,4,Face(count)%Left ) 
                    ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%paramL => parbnd(:,:,:,4,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,:,:,4,Face(count)%Left ) 
              ELSE              
                Left  = idxe(i,j-1,k) 
                Right = idxe(i,j,k) 
                IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                    IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        cOuterFaces = cOuterFaces + 1 
                        count = cOuterFaces 
                    ENDIF 
                ELSE
                    CYCLE  
                ENDIF  
                Face(count)%Left  = Global2LocalElem(idxe(i,j-1,k)) 
                Face(count)%Right = Global2LocalElem(idxe(i,j,k)) 
                Face(count)%qL => qBnd(:,:,:,4,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,:,:,3,Face(count)%Right) 
                Face(count)%FL => FBnd(:,:,:,4,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,:,:,3,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,:,:,4,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,:,:,3,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 1., 0. /)                               ! set face normal vector 
              Face(count)%x0 = xL + (/ REAL(i-1), REAL(j-1), REAL(k-1) /)*dx  ! set face lower left corner 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    ! z faces 
    IF(nDim>=3) THEN
        DO k = 1, KMAX+iper(3) 
         DO j = 1, JMAX 
          DO i = 1, IMAX  
              IF(k.EQ.1) THEN
                IF(Periodic(3)) THEN
                    Left  = idxe(i,j,KMAX) 
                    Right = idxe(i,j,k) 
                    IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                        IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                            cInnerFaces = cInnerFaces + 1 
                            count = cInnerFaces 
                        ELSE
                            cOuterFaces = cOuterFaces + 1 
                            count = cOuterFaces 
                        ENDIF 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = Global2LocalElem(idxe(i,j,KMAX)) 
                    Face(count)%Right = Global2LocalElem(idxe(i,j,k))  
                    Face(count)%qL => qBnd(:,:,:,6,Face(count)%Left ) 
                    Face(count)%qR => qBnd(:,:,:,5,Face(count)%Right) 
                    Face(count)%FL => FBnd(:,:,:,6,Face(count)%Left ) 
                    Face(count)%FR => FBnd(:,:,:,5,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,:,:,6,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,:,:,5,Face(count)%Right) 
                ELSE
                    Left  = 0 
                    Right = idxe(i,j,k)
                    IF( CPUDistribution(Right)==myrank ) THEN
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = 0 
                    Face(count)%Right = Global2LocalElem(idxe(i,j,k)) 
                    ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%qR => qBnd(:,:,:,5,Face(count)%Right) 
                    Face(count)%FR => FBnd(:,:,:,5,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,:,:,5,Face(count)%Right) 
                    Face(count)%paramR => parbnd(:,:,:,5,Face(count)%Right) 
                ENDIF 
              ELSEIF(k.EQ.KMAX+1) THEN 
                    Left  = idxe(i,j,k-1) 
                    Right = 0 
                    IF( CPUDistribution(Left)==myrank ) THEN
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        CYCLE  
                    ENDIF  
                    Face(count)%Left  = Global2LocalElem(idxe(i,j,k-1)) 
                    Face(count)%Right = 0 
                    Face(count)%qL => qBnd(:,:,:,6,Face(count)%Left ) 
                    Face(count)%FL => FBnd(:,:,:,6,Face(count)%Left ) 
                    ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%paramL => parbnd(:,:,:,6,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,:,:,6,Face(count)%Left ) 
              ELSE              
                Left  = idxe(i,j,k-1) 
                Right = idxe(i,j,k) 
                IF( CPUDistribution(Left)==myrank .OR. CPUDistribution(Right)==myrank ) THEN
                    IF( CPUDistribution(Left) == CPUDistribution(Right) ) THEN 
                        cInnerFaces = cInnerFaces + 1 
                        count = cInnerFaces 
                    ELSE
                        cOuterFaces = cOuterFaces + 1 
                        count = cOuterFaces 
                    ENDIF 
                ELSE
                    CYCLE  
                ENDIF  
                Face(count)%Left  = Global2LocalElem(idxe(i,j,k-1))
                Face(count)%Right = Global2LocalElem(idxe(i,j,k))  
                Face(count)%qL => qBnd(:,:,:,6,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,:,:,5,Face(count)%Right) 
                Face(count)%FL => FBnd(:,:,:,6,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,:,:,5,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,:,:,6,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,:,:,5,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 0., 1. /)                               ! set face normal vector 
              Face(count)%x0 = xL + (/ REAL(i-1), REAL(j-1), REAL(k-1) /)*dx  ! set face lower left corner 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    ! 
    IF(EvalSphereIntegral) THEN
        !
      ALLOCATE(SpherePickPoint(nSpherePickPoint))
      CALL DefineSpherePickpoints
      ! Locate Spherepickpoints in the primary mesh 
      !
      DO i = 1, nSpherePickPoint
         ii = MIN( IMAX, MAX( 1, CEILING( (SpherePickPoint(i)%x(1)-xL(1))/dx(1) ) ) ) 
         jj = MIN( JMAX, MAX( 1, CEILING( (SpherePickPoint(i)%x(2)-xL(2))/dx(2) ) ) ) 
         kk = MIN( KMAX, MAX( 1, CEILING( (SpherePickPoint(i)%x(3)-xL(3))/dx(3) ) ) ) 
         SpherePickPoint(i)%iElem0    = idxe(ii,jj,kk) 
         SpherePickPoint(i)%iElem     = SpherePickPoint(i)%iElem0 
         !CALL GetL0Tri(tmptri, SpherePickPoint(i)%iElem0)  
         SpherePickPoint(i)%xi        = ( SpherePickPoint(i)%x(:)-x(:,tri(1,SpherePickPoint(i)%iElem0 )) )/dx 
         SpherePickPoint(i)%picktime  = 0. 
         SpherePickPoint(i)%currlevel = 0   
         !
         CONTINUE
         !
      ENDDO 
      !
      DO i = 1, nSpherePickPoint
        SpherePickPoint(i)%iElem0 = Global2LocalElem( SpherePickPoint(i)%iElem0 ) 
        SpherePickPoint(i)%iElem  = Global2LocalElem( SpherePickPoint(i)%iElem  ) 
        ! If SpherePickpoint is not in the current CPU, then set element numbers to zero 
        IF(SpherePickPoint(i)%iElem0.LE.0) THEN
            SpherePickPoint(i)%iElem0 = 0 
            SpherePickPoint(i)%iElem  = 0 
        ELSE 
            !PRINT *, ' CPU found SpherePickpoint ', myrank, i 
        ENDIF 
      ENDDO 
      !
      ENDIF
    ! 
    !CALL LocateSpherePickPoints  
    !
    DEALLOCATE( idxe ) 
    !
    ! Now, let us compute some of the important matrices in the ADER-DG method 
    ! 
    MM   = 0.   ! Element mass matrix 
    Kxi  = 0.   ! Element stiffness matrix 
    dudx = 0.   ! discrete derivative operator, which projects the derivatives onto the basis  
    DO iGP = 1, N+1 
        CALL BaseFunc1D(phi,phi_xi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! i) Mass-matrix 
                MM(k,l) = MM(k,l) + wGPN(iGP)*phi(k)*phi(l) 
                ! ii) Stiffness-matrix 
                Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  
            ENDDO
        ENDDO        
    ENDDO    
    CALL MatrixInverse(N+1,MM,iMM) 
    dudx = MATMUL( iMM, TRANSPOSE(Kxi) ) 
    !
    Kxixi  = 0.   ! Element stiffness matrix for the fourth order AV operator 
    DO iGP = 1, N+1 
        CALL BaseFunc1Dxx(phi,phi_xi,phi_xixi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! AV stiffness matrix 
                Kxixi(k,l) = Kxixi(k,l) + wGPN(iGP)*phi_xixi(k)*phi_xixi(l)  
            ENDDO
        ENDDO        
    ENDDO   
    CALL RG(N+1,N+1,Kxixi,Lxixi,ImLambda,1,Rxixi,itemp,rtemp,ierr)
    CALL MatrixInverse(N+1,Rxixi,iRxixi)
    CONTINUE
    ! 
!    OPEN(FILE='dudx4.txt', UNIT=333) 
!    do l = 1, N+1
!     do k = 1, N+1
!        WRITE(333,'(a,i2,a,i2,a,e,a)') ' dudx[4][',k-1,'][',l-1,'] = ', dudx(l,k), ';'  
!     enddo 
!    enddo
!    CLOSE(333) 
!    STOP 

    CALL BaseFunc1D(phi0,phi_xi,0.0) ! Compute the basis functions on the left 
    CALL BaseFunc1D(phi1,phi_xi,1.0) ! Compute the basis function on the right 
    ! The flux matrices are all possible combinations of left and right 
    DO k = 1, N+1
        DO l = 1, N+1
            FLm(k,l) = phi0(k)*phi1(l)   ! Left contribution to the left flux matrix    (m = left  of the interface)  
            FLp(k,l) = phi0(k)*phi0(l)   ! Right contribution to the left flux matrix   (p = right of the interface) 
            FRm(k,l) = phi1(k)*phi1(l)   ! Left contribution to the right flux matrix   (m = left  of the interface) 
            FRp(k,l) = phi1(k)*phi0(l)   ! Right contribution to the right flux matrix  (p = right of the interface) 
        ENDDO
    ENDDO        
    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle) 
    F0 = phi0   ! upwinding in time = information comes from smaller times 
    F1 = FRm    ! upwinding in time = information comes from smaller times  
    K1 = F1 - Kxi 
    CALL MatrixInverse(N+1,K1,iK1)   
    FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary 
    FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary 
    !
    ! For point sources, we also need the time mass matrix of the temporal Taylor series
    !
    MT   = 0.   ! Time mass matrix 
    DO iGP = 1, N+1 
        CALL TimeBaseFunc1D(phi,phi_xi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! i) Mass-matrix 
                MT(k,l) = MT(k,l) + wGPN(iGP)*phi(k)*phi(l) 
            ENDDO
        ENDDO        
    ENDDO    
    CALL MatrixInverse(N+1,MT,iMT)     !
    ! For the fine output of each spatial degree of freedom onto a subgrid... 
    DO i = 1, N+1 
       subxi(i) = REAL(i-1)/REAL(N) 
    ENDDO
    cnt = 0 
    DO k = 1, N+1
     DO j = 1, N+1 
      DO i = 1, N+1  
         cnt = cnt + 1 
         CALL BaseFunc1D(phi_i,phi_xi,subxi(i))
         CALL BaseFunc1D(phi_j,phi_xi,subxi(j))
         CALL BaseFunc1D(phi_k,phi_xi,subxi(k))
         count = 0 
         DO kk = 1, nDOF(3) 
          DO jj = 1, nDOF(2)  
           DO ii = 1, nDOF(1) 
             count = count + 1 
             aux = (/ phi_i(ii), phi_j(jj), phi_k(kk) /) 
             SubOutputMatrix(count,cnt) = PRODUCT( aux(1:nDim) )                      
           ENDDO
          ENDDO
         ENDDO
         !
       ENDDO
      ENDDO
    ENDDO
     ! 
     ! subgrid triangulation for subgrid output 
     !
     IF(N==0) THEN
        allsubxi = 0.5
        subtri   = 0 
     ELSE        
         ALLOCATE( idxn(N+1,N+1,N+1) )
         idxn = 0 
         c = 0 
         DO k = 1, N+1
            DO j = 1, N+1 
               DO i = 1, N+1 
                  c = c + 1 
                  allsubxi(:,c) = (/ REAL(i-1)/REAL(N), REAL(j-1)/REAL(N), REAL(k-1)/REAL(N) /)  
                  idxn(i,j,k) = c
               ENDDO
            ENDDO
         ENDDO
         c = 0 
         DO k = 1, N 
            DO j = 1, N
               DO i = 1, N 
                  c = c + 1 
                  subtri(:,c) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i+1,j+1,k), idxn(i,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i+1,j+1,k+1), idxn(i,j+1,k+1) /)         
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE( idxn )     
     ENDIF 
     ! 
     ! subgrid triangulation for subgrid output of the limiter 
     ! 
     ALLOCATE( idxn(nSubLim+1,nSubLim+1,nSubLim+1) )
     idxn = 0 
     c = 0 
     DO k = 1, nSubLim+1
      DO j = 1, nSubLim+1 
        DO i = 1, nSubLim+1  
            c = c + 1 
            idxn(i,j,k) = c
            subxilim(:,c) = (/ (REAL(i)-1.)/REAL(nSubLim), (REAL(j)-1.)/REAL(nSubLim), (REAL(k)-1.)/REAL(nSubLim) /)  
        ENDDO
      ENDDO
     ENDDO
     c = 0 
     DO k = 1, nSubLim 
      DO j = 1, nSubLim
        DO i = 1, nSubLim  
            c = c + 1 
            subtrilim(:,c) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i+1,j+1,k), idxn(i,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i+1,j+1,k+1), idxn(i,j+1,k+1) /)         
        ENDDO
      ENDDO
     ENDDO
     !
    ! 
    ! Now prepare the subcell limiter 
    !
    ALLOCATE( olduh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )       ! Allocate the old spatial degrees of freedom (needed for the limiter) 
    ALLOCATE( Limiter(nElem), recompute(nElem) ) 
    DO i = 1, nElem
        Limiter(i)%oldstatus = 0 
        Limiter(i)%status = 0
        recompute(i) = 0 
    ENDDO 
    nSubLimV = 1 
    DO i = 1, nDim
        nSubLimV(i) = nSubLim 
    ENDDO    
    !
    ALLOCATE( LSQM(N+2,N+2), iLSQM(N+2,N+2), LSQrhs(N+2,nSubLim) )  
    ALLOCATE( uh2lim(nSubLim,N+1), lim2uh(N+1,nSubLim)           )
    ALLOCATE( uh2lob(N+1,N+1)                                    )   
    uh2lim = 0. 
    dxi = 1./REAL(nSubLim)    ! the limiter uses nSubLim subintervals in each cell 
    DO i = 1, nSubLim
     xi1 = REAL(i-1)*dxi     ! left sub-interval boundary 
     xi2 = REAL(i  )*dxi     ! right sub-interval boundary 
     DO ii = 1, N+1  
        xi = xi1 + dxi*xiGPN(ii)     
        CALL BaseFunc1D(phi,phi_xi,xi) 
        uh2lim(i,:) = uh2lim(i,:) + wGPN(ii)*phi(:)
     ENDDO
    ENDDO
    ! To go back from the limited (piecewise constant sub-cell solution) to the uh solution
    ! we use a constrained least-squares reconstruction, which preserves the average exactly 
    LSQM(1:N+1,1:N+1) = 2*MATMUL( TRANSPOSE(uh2lim), uh2lim ) 
    LSQM(N+2,1:N+1)   =  wGPN(:) 
    LSQM(1:N+1,N+2)   = -wGPN(:) 
    LSQM(N+2,N+2) = 0. 
    CALL MatrixInverse(N+2,LSQM,iLSQM)
    LSQrhs(1:N+1,:)   = 2*TRANSPOSE(uh2lim)  
    LSQrhs(N+2,:)     = dxi    
    lim2uh = MATMUL( iLSQM(1:N+1,:), LSQrhs ) 
    ! 
    phi0 = 1. 
    test = MATMUL(uh2lim,phi0) 
    testmatrix = MATMUL(lim2uh,uh2lim) 
    testmatrix2 = MATMUL(uh2lim,lim2uh)     
    !
    ! Compute the Gauss-Lobatto quadrature nodes 
    IF(N>0) THEN
        CALL gaulob(0.,1.,xiLob,wLob,N+1)
        DO ii = 1, N+1
            CALL BaseFunc1D(phi,phi_xi,xiLob(ii))  
            uh2lob(ii,:) = phi(:) 
        ENDDO 
    ELSE
        xiLob = 0.0
        wLob  = 1.0 
        uh2lob = 1.0 
    ENDIF    
    !
    CONTINUE 
    !
    DEALLOCATE( LSQM, iLSQM, LSQrhs ) 
     !
     DO i = 1, nSubLim
         xilimbary(i) = 0.0 + 0.5/REAL(nSubLim) + (i-1)*1.0/REAL(nSubLim)
     ENDDO     
     !
     ! Reference element 
     !
     ReferenceElement(:,1) = (/ 0., 0., 0. /) 
     ReferenceElement(:,2) = (/ 1., 0., 0. /) 
     ReferenceElement(:,3) = (/ 1., 1., 0. /) 
     ReferenceElement(:,4) = (/ 0., 1., 0. /) 
     ReferenceElement(:,5) = (/ 0., 0., 1. /) 
     ReferenceElement(:,6) = (/ 1., 0., 1. /) 
     ReferenceElement(:,7) = (/ 1., 1., 1. /) 
     ReferenceElement(:,8) = (/ 0., 1., 1. /) 
    !
    ! Set the initial condition. Here, we assume a nodal basis. Otherwise, we would have to do formal L2 projection, 
    ! i.e. integration of the initial condition and multiplication with the inverse mass matrix. 
    !
    DO iElem = 1, nElem
        x0 = x(:,tri(1,iElem)) ! get the coordinate of the lower left node 
        DO k = 1, nDOF(3) 
         DO j = 1, nDOF(2)
          DO i = 1, nDOF(1) 
              IF(nDim==2) THEN
                xGP = x0 + (/ xiGPN(i), xiGPN(j), 0.5 /)*dx(:)                               
              ELSE                   
                xGP = x0 + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:)                               
              ENDIF
              ! 
              CALL InitialField(u0,par0,xGP,time) 
              ! 
              !if(myrank==2 .and. iElem==1729) then
              !    print *, 'problematic zone ', xGP, iElem 
              !    print *, ' u0 = ', u0 
              !endif 
              !if(nCPU==1 .and. iElem==7201) then 
              !    continue
              !endif 
              
              uh(:,i,j,k,iElem) = u0 
              parh(:,i,j,k,iElem) = par0  
          ENDDO
         ENDDO
        ENDDO
#ifndef NOLIMITER         
        IF(N>0) THEN
            DO iVar = 1, nVar 
                Limiter(iElem)%lmin(iVar) = MINVAL(uh(iVar,:,:,:,iElem)) 
                Limiter(iElem)%lmax(iVar) = MAXVAL(uh(iVar,:,:,:,iElem)) 
            ENDDO         
            CALL DMP(dmpresult,uh(:,:,:,:,iElem),Limiter(iElem),1e-1)
            IF(.NOT.dmpresult) THEN
                ! If initial condition projection does not satisfy the DMP, then activate the subcell limiter 
                IF(Limiter(iElem)%oldstatus.EQ.0) THEN 
                   ALLOCATE( Limiter(iElem)%Lh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))    ) 
                   ALLOCATE( Limiter(iElem)%NewLh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) ) 
                ENDIF
                Limiter(iElem)%status    = 1 
                Limiter(iElem)%oldstatus = 1 
                DO k = 1, nSubLimV(3) 
                 DO j = 1, nSubLimV(2)
                  DO i = 1, nSubLimV(1) 
                      IF(nDim==2) THEN
                        xGP = x0 + (/ xilimbary(i), xilimbary(j), 0.5 /)*dx(:) 
                      ELSE
                        xGP = x0 + (/ xilimbary(i), xilimbary(j), xilimbary(k) /)*dx(:) 
                      ENDIF 
                      CALL InitialField(u0,par0,xGP,0.0) 
                      Limiter(iElem)%Lh(:,i,j,k) = u0 
                      Limiter(iElem)%NewLh(:,i,j,k) = u0 
                  ENDDO
                 ENDDO
                ENDDO             
            ENDIF              
        ENDIF
#endif         
        !
#ifdef Z4EINSTEIN 
        !
        IF(EQN%EinsteinAutoAux==1) THEN
            ! x derivatives             
            DO k = 1, nDOF(3) 
             DO j = 1, nDOF(2)
                  ux(:,:,j,k) = 1./dx(1)*MATMUL( uh(:,:,j,k,iElem), TRANSPOSE(dudx) ) 
             ENDDO
            ENDDO
            ! y derivatives 
            IF(nDim>=2) THEN 
             DO k = 1, nDOF(3) 
              DO i = 1, nDOF(1) 
                  uy(:,i,:,k) = 1./dx(2)*MATMUL( uh(:,i,:,k,iElem), TRANSPOSE(dudx) ) 
              ENDDO
             ENDDO
            ELSE
                uy = 0.0 
            ENDIF  
            ! z derivatives 
            IF(nDim>=3) THEN 
             DO j = 1, nDOF(2) 
              DO i = 1, nDOF(1) 
                  uz(:,i,j,:) = 1./dx(3)*MATMUL( uh(:,i,j,:,iElem), TRANSPOSE(dudx) ) 
              ENDDO
             ENDDO
            ELSE
                uz = 0.0 
            ENDIF          
            !
            DO k = 1, nDOF(3) 
             DO j = 1, nDOF(2) 
              DO i = 1, nDOF(1) 
                 uh(24,i,j,k,iElem) = ux(17,i,j,k)/uh(17,i,j,k,iElem)   ! A_1   
                 uh(25,i,j,k,iElem) = uy(17,i,j,k)/uh(17,i,j,k,iElem)   ! A_2  
                 uh(26,i,j,k,iElem) = uz(17,i,j,k)/uh(17,i,j,k,iElem)   ! A_3
                 !
                 uh(27,i,j,k,iElem) = ux(18,i,j,k)   ! B_11   
                 uh(28,i,j,k,iElem) = uy(18,i,j,k)   ! B_21
                 uh(29,i,j,k,iElem) = uz(18,i,j,k)   ! B_31
                 uh(30,i,j,k,iElem) = ux(19,i,j,k)   ! B_12
                 uh(31,i,j,k,iElem) = uy(19,i,j,k)   ! B_22
                 uh(32,i,j,k,iElem) = uz(19,i,j,k)   ! B_32
                 uh(33,i,j,k,iElem) = ux(20,i,j,k)   ! B_13
                 uh(34,i,j,k,iElem) = uy(20,i,j,k)   ! B_23
                 uh(35,i,j,k,iElem) = uz(20,i,j,k)   ! B_33
                 !
                 uh(36,i,j,k,iElem) = 0.5*ux(1,i,j,k)    ! D_111   
                 uh(37,i,j,k,iElem) = 0.5*ux(2,i,j,k)    ! D_112   
                 uh(38,i,j,k,iElem) = 0.5*ux(3,i,j,k)    ! D_113   
                 uh(39,i,j,k,iElem) = 0.5*ux(4,i,j,k)    ! D_122   
                 uh(40,i,j,k,iElem) = 0.5*ux(5,i,j,k)    ! D_123   
                 uh(41,i,j,k,iElem) = 0.5*ux(6,i,j,k)    ! D_133   

                 uh(42,i,j,k,iElem) = 0.5*uy(1,i,j,k)    ! D_211   
                 uh(43,i,j,k,iElem) = 0.5*uy(2,i,j,k)    ! D_212   
                 uh(44,i,j,k,iElem) = 0.5*uy(3,i,j,k)    ! D_213   
                 uh(45,i,j,k,iElem) = 0.5*uy(4,i,j,k)    ! D_222   
                 uh(46,i,j,k,iElem) = 0.5*uy(5,i,j,k)    ! D_223   
                 uh(47,i,j,k,iElem) = 0.5*uy(6,i,j,k)    ! D_233   

                 uh(48,i,j,k,iElem) = 0.5*uz(1,i,j,k)    ! D_311   
                 uh(49,i,j,k,iElem) = 0.5*uz(2,i,j,k)    ! D_312   
                 uh(50,i,j,k,iElem) = 0.5*uz(3,i,j,k)    ! D_313   
                 uh(51,i,j,k,iElem) = 0.5*uz(4,i,j,k)    ! D_322   
                 uh(52,i,j,k,iElem) = 0.5*uz(5,i,j,k)    ! D_323   
                 uh(53,i,j,k,iElem) = 0.5*uz(6,i,j,k)    ! D_333   
                 
              ENDDO
             ENDDO
            ENDDO 
            
            !
        ENDIF 
        !
#endif     
        !
    ENDDO    
    !
    !
    !if(myrank==2) then
    !    print *, ' final values for uh in the troubled zone: ' 
    !    print *, uh(1:6,:,:,:,1729) 
    !    stop 
    !endif
    !if(nCPU==1) then
    !    print *, uh(1:6,:,:,:,7201)
    !    stop 
    !endif    
    !
    IF(nParam.GT.0) THEN
        !
        ! Compute the bounday-extrapolated values for the material parameters 
        !
        DO iElem = 1, nElem 
            IF(iElem==0) CYCLE 
            parbnd(:,:,:,:,iElem) = 0. 
            ! x-direction: face 1 (left) and face 2 (right) 
            DO k = 1, nDOF(3) 
             DO j = 1, nDOF(2) 
                parbnd(:,j,k,1,iElem) = MATMUL( parh(:,:,j,k,iElem),   FLCoeff )   ! left 
                parbnd(:,j,k,2,iElem) = MATMUL( parh(:,:,j,k,iElem),   FRCoeff )   ! right 
             ENDDO
            ENDDO 
            ! y-direction: face 3 (left) and face 4 (right) 
            IF(nDim>=2) THEN
                DO k = 1, nDOF(3) 
                 DO i = 1, nDOF(1) 
                    parbnd(:,i,k,3,iElem) = MATMUL( parh(:,i,:,k,iElem),   FLCoeff )   ! left 
                    parbnd(:,i,k,4,iElem) = MATMUL( parh(:,i,:,k,iElem),   FRCoeff )   ! right 
                 ENDDO
                ENDDO 
            ENDIF    
            ! z-direction: face 5 (left) and face 6 (right) 
            IF(nDim>=3) THEN
                DO j = 1, nDOF(2) 
                 DO i = 1, nDOF(1) 
                    parbnd(:,i,j,5,iElem) = MATMUL( parh(:,i,j,:,iElem),   FLCoeff )   ! left 
                    parbnd(:,i,j,6,iElem) = MATMUL( parh(:,i,j,:,iElem),   FRCoeff )   ! right 
                 ENDDO
                ENDDO 
            ENDIF            
            !
        ENDDO
    ENDIF    
    !
    CALL MPIExchangeParameters
    !
    ! Deploy energy in one single element... 
    !i = PointSrc(1)%iElem
    !uh(8,:,:,:,i) = -10.0 
    !
    !CALL DickeBertha(1)   
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) || defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD)  || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)  
    !    CALL RunTimeAnalyseData
#endif 
    !
    CALL WriteData
    !
    CONTINUE 
    !
END SUBROUTINE ADERDGInit
    
    
SUBROUTINE InitialField(u0,par,xGP,tGP) 
#ifdef GRMHD
    USE typesDef, ONLY : nVar, nDim, nParam, d, ICType, EQN, ICType2,ExcisionRadius, dx  , aom, Mbh   
#else
    USE typesDef, ONLY : nVar, nDim, nParam, d, ICType, EQN, ICType2,ExcisionRadius, dx
#endif  
    USE Bessel_mod 
#ifdef TWOPUNCTURES  
	USE TwoPunctures_mod, ONLY : TwoPunctures_Interpolate
#endif 
#ifdef RNSID 
    USE RNSID_mod, ONLY : RNSID_Interpolate
#endif
#ifdef TORUS    
    USE Metric_mod
#endif
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN ) :: xGP(d), tGP        ! spatial position vector and time 
    REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
    REAL, INTENT(OUT) :: par(nParam)        ! material parameter vector 
    ! Local variables 
    INTEGER :: i,j,k,l,nm,iNewton, Maxnewton = 50    
    REAL :: VBase(nVar), ampl(nVar), sigma(d) 
    REAL :: V0(nVar),r,VLL(nVar),VRR(nVar), VZ4(54),V70(70)  
    REAL :: du,dv,drho,dTemp,dp,epsilon,xc(d) 
    REAL :: omega, parL(nParam), parR(nParam) 
    REAL :: r1(9), lambda, lambdat, lambdax, mu, rho, cs, cp, ICA, HH, dxH, Kxx
#ifndef GRMHD
    REAL :: aom,Mbh
#endif
    REAL :: theta, a,  M, zz, r_hz, delta, g_cov(3,3), g_contr(3,3), Kex(3,3), Aex(3,3) 
    REAL :: traceK, alpha, fa, ggg, fff, k0, sig, AA(3), BB(3,3), BMat(3,3), DD(3,3,3), PP(3)  
    REAL :: rho0, p0, eta, b0, tempaa, tempab, tempac, va2, vax, bv(3), BV_contr(3), vv(3), gm, gp, ddet(3) 
    REAL :: beta(3), betadown(3), dbetadown(3,3), nablabeta(3,3), test(3,3), Gtilde(3)   
    REAL :: Christoffel(3,3,3), dxphi, bbb, dxb, dtb, kyy, kzz, Aref, x0, re0, ms, uu0  
    REAL :: vv_cov(3), shift(3), p, lapse, gammaij(6), ng, up, phi,phi2, dphi(3), rc, vc2, tc, pc, rhoc, lf, vtheta, vphi, vx, vy, vz 
    REAL :: urr, tt, c1, c2, df, vr, ut, f, urc, g_tt, betaru, dtt, detg, psi , xloc(d), minR
    REAL :: BJ(0:1),DJ(0:1),BY(0:1),DY(0:1) 
    REAL :: IJ(0:1),DIJ(0:1),IY(0:1),DIY(0:1) 
    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    REAL :: xGP_sph(3),A_coord(3,3), detA,iA_coord(3,3),A_coord_contr(3,3),iA_coord_contr(3,3)
    !
    u0  = 0.0
    par = 0.0 
    !
    ! no local material parameters for Euler equations 
    ! Gaussian perturbation 
    SELECT CASE(TRIM(ICType)) 
    CASE('SmoothWave')
#ifdef EULER
        sigma = (/ 0.05, 0.05, 0.05 /)       ! half-width
        VBase(:) = (/ 1., 0., 0., 0., 1. /)  ! base-state 
        ampl(:)  = 0.                        ! perturbation amplitude vector 
        ampl(5)   = 1e-3                     ! 
        V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM(xGP(1:nDim)**2/sigma(1:nDim)**2) )
#endif         
    CASE('EulerGauss')
#ifdef EULER
        V0 = 0.0  
        ! Gaussian perturbation 
        sigma = (/ 10., 10., 10. /)           ! half-width
        VBase(:) = (/ 1.4, 0., 0., 0., 1e6 /) ! base-state 
        ampl(:)  = 0.                         ! perturbation amplitude vector 
        ampl(5)   = 100.0                     ! 
        xc = (/ 2000.0, 1950.0, 0.0 /)        ! 1950 
        V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM((xGP(1:nDim)-xc(1:nDim))**2/sigma(1:nDim)**2) )  
#endif         
    CASE('Sod')     
#ifdef EULER
        r = SQRT( SUM(xGP(1:nDim)**2) )     
        VLL = (/ 1.0,   0.0, 0.0, 0.0, 1.0 /) 
        VRR = (/ 0.125, 0.0, 0.0, 0.0, 0.1 /) 
        IF(xGP(1)<0) THEN 
            V0 = VLL
        ELSE
            V0 = VRR        
        ENDIF    
#endif         
    CASE('Const')     
#ifdef EULER
        r = SQRT( SUM(xGP(1:nDim)**2) )     
        VLL = (/ 1.0, 0.0, 0.0, 0.0, 1.0 /) 
        VRR = (/ 1.0, 0.0, 0.0, 0.0, 1.0 /) 
        IF(xGP(1)<0) THEN 
            V0 = VLL
        ELSE
            V0 = VRR        
        ENDIF    
#endif         
    CASE('EP2D')     
#ifdef EULER
        r = SQRT( SUM(xGP(1:nDim)**2) )     
        VLL = (/ 1.0,   0.0, 0.0, 0.0, 1.0 /) 
        VRR = (/ 0.125, 0.0, 0.0, 0.0, 0.1 /) 
        IF(r.LT.0.5) THEN
        !IF(xGP(1)<0) THEN 
            V0 = VLL
        ELSE
            V0 = VRR        
        ENDIF    
#endif         
    CASE('ShuVortex2D')
        !
        continue 
        !
#ifdef EULER 
       epsilon = 5.0
       r = SQRT((xGP(1)-tGP-5.)**2+(xGP(2)-tGP-5.)**2)
       du = epsilon/2./EQN%Pi*exp(0.5*(1.-r*r))*(5. - xGP(2) + tGP)
       dv = epsilon/2./EQN%Pi*exp(0.5*(1.-r*r))*(xGP(1)  - 5.- tGP)
       dTemp = -(EQN%gamma-1.)*epsilon**2/8./EQN%gamma/EQN%Pi**2*exp(1.-r*r)
       drho = (1.+dTemp)**(1./(EQN%gamma-1.))-1.
       dp   = (1.+dTemp)**(EQN%gamma/(EQN%gamma-1.))-1.
       !
       V0(1) = 1. + drho
       V0(2) = 1. + du
       V0(3) = 1. + dv
       V0(4) = 0.0
       V0(5) = 1. + dp
#endif        
       !
#ifdef GPR3D  
       epsilon = 5.0
       r = SQRT((xGP(1)-tGP-5.)**2+(xGP(2)-tGP-5.)**2)
       du = epsilon/2./Pi*exp(0.5*(1.-r*r))*(5. - xGP(2) + tGP)
       dv = epsilon/2./Pi*exp(0.5*(1.-r*r))*(xGP(1)  - 5.- tGP)
       dTemp = -(EQN%gamma-1.)*epsilon**2/8./EQN%gamma/Pi**2*exp(1.-r*r)
       drho = (1.+dTemp)**(1./(EQN%gamma-1.))-1.
       dp   = (1.+dTemp)**(EQN%gamma/(EQN%gamma-1.))-1.
       !
       V0 = 0.0 
       ! Density 
       V0(1) = 1. + drho 
       ! Velocity    
       V0(2) = 1.  + du
       V0(3) = 1.  + dv
       V0(4) = 0.
       ! Pressure  
       V0(5) = 1.  + dp
       ! distortion tensor  
       V0(6:14) = V0(1)**(1./3.)*(/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /)   
       ! thermal impulse
       V0(15:17) = 0.0 
       !
       continue
       !
#endif        
       ! 
    CASE('ShockNavierStokes')
       !         
#ifdef GPR3D    
        V0(:) = 0.0   
        ! 
        x0   = 0.0        
        rho0 = 1.0 
        uu0  = 0.0 
        Ms   = 2.0 
        Re0  = rho0 * Ms / EQN%mu         
        p0   = 1./EQN%gamma 
        CALL NSShock(V0(1),vx,V0(5),xGP(1)-Ms*tGP-x0,EQN%gamma,EQN%mu,Re0,Ms,rho0,uu0,p0)
        V0(2:4) = (/ vx, 0.0, 0.0 /) 
        Aref = V0(1)**(1./3.) 
        V0(6:14) = (/ Aref, 0., 0., 0., Aref, 0., 0., 0., Aref /)         
#endif         
       !
    END SELECT 

#ifdef ELASTICITY
    ! Set the local material parameters 
    ! lambda, mu, rho 
    !parR = (/ 4.0, 2.0, 1.0 /) 
    !parL = (/ 2.0, 1.0, 1.0 /) 
    !omega = 0.5*(1+erf(2*(xGP(3)-0.6)))  
    !par = (1.0-omega)*parL + omega*parR  
    
    !parR = (/ 7.509672500e9, 7.50916375e9, 2200. /) 
    !parL = parR 
    !!parL = (/ 3.509672500e9, 3.50916375e9, 2200. /) 
    !omega = 0.5*(1+erf(0.01*(xGP(2)-1500.0)))  
    !par = (1.0-omega)*parL + omega*parR  
    !
    !V0 = 0.0  
    !! Gaussian perturbation 
    !sigma = (/ 100., 100., 100. /)          ! half-width
    !VBase(:) = 0.0                       ! base-state 
    !ampl(:)  = 0.                        ! perturbation amplitude vector 
    !ampl(8)   = 0.0 ! 1e-9                     ! 
    !xc = (/ 2000.0, 1800.0, 0.0 /) 
    !V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM((xGP(1:nDim)-xc(1:nDim))**2/sigma(1:nDim)**2) )   
    
    !V0(1) = 1.0 + 0.1*xGP(1) + 0.2*xGP(2)  

    ! par = (/ 1.0, 1.0, 1.0 /)
 
    SELECT CASE(TRIM(ICType)) 
    CASE('ElasticSineWave')
        ! Sin wave test 
        lambda = 1.0 
        mu     = 1.0 
        rho    = 1.0 
        cp = SQRT((lambda+2*mu)/rho) 
        cs = SQRT(mu/rho) 
        !par = (/ lambda, mu, rho /) 
        ! 
        r1 = (/ 0., 0., 0., mu, 0., 0., 0., -cs, 0.0 /)   ! Eigenvector for the shear wave
        V0(1:9) = SIN(2*ACOS(-1.0)*xGP(1))*r1
        V0(10:14) = (/ lambda, mu, rho, 1.0, 1.0 /) 
    CASE('RPGaussian')
        ! Sin wave test 
        lambda = 2.0 
        mu     = 1.0 
        rho    = 1.0 
        cp = SQRT((lambda+2*mu)/rho) 
        cs = SQRT(mu/rho) 
        !par = (/ lambda, mu, rho /) 
        ! 
        r1 = (/ lambda+2*mu, lambda, lambda, 0., 0., 0., -cp, 0., 0. /)   ! Eigenvector for the shear wave
        V0(1:9) = 0.1*EXP(-0.5*(xGP(1)+0.5)**2/0.05**2)*r1 
        IF(xGP(1)<0.0) THEN
            V0(10:14) = (/ lambda, mu, rho, 1.0, 1.0 /) 
        ELSE
            V0(10:14) = (/ lambda, mu, rho, 0.0, 1.0 /) 
        ENDIF
        ! 
    CASE DEFAULT
        !
        !par = (/ 7.509672500e9, 7.50916375e9, 2200. /) 
        V0  = 0.0 
        V0(10:14) = (/ 7.509672500e9, 7.50916375e9, 2200., 1.0, 1.0 /) 
        ! 
    END SELECT 
     
#endif 

#ifdef ACOUSTIC
    par = 0.0 

    EQN%c0 = 1.0 
    
    V0 = 0.0  
    ! Gaussian perturbation 
    sigma = (/ 0.1, 0.1, 0.1 /)          ! half-width
    VBase(:) = (/ 1., 0., 0., 0. /)     ! base-state 
    ampl(:)  = 0.                        ! perturbation amplitude vector 
    ampl(1)   = 1e-2                     ! 
    xc = (/ 0.0, 0.0, 0.0 /)       ! 1950 
    V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM((xGP(1:nDim)-xc(1:nDim))**2/sigma(1:nDim)**2) )   
    
    !V0(1) = 1.0 + 0.1*xGP(1) + 0.2*xGP(2)  

#endif 
    ! A simple debug check for the computation of derivatives 
    !u0 = 0. 
    !u0(1) = 0.123 !*xGP(1) 
    !
#ifdef Z4EINSTEIN
    SELECT CASE(ICType)
    CASE('Z4GaugeWave')
        V0 = 0.0
        !
        ICA = 0.1
        !  
        ! Gauge wave propagating along x 
        HH     = 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        dxH    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        Kxx    = - Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))/SQRT( 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))  )  ! extrinsic curvature  
        !
        V0(1)  = HH          ! \gamma_xx
        V0(4)  = 1.0         ! \gamma_yy
        V0(6)  = 1.0         ! \gamma_zz
        V0(7)  = Kxx         ! K_xx
        V0(10) = 0.0         ! K_yy
        V0(12) = 0.0         ! K_zz

        V0(17) = SQRT(HH)
        !
        V0(21:22) = xGP(1:2)        
        !
        ! Auxiliary variables
        V0(24) = 1.0/(2.0*HH)*dxH           ! A_x  
        V0(36) = 0.5*dxH                    ! D_xxx
        !
    CASE('Z4LinearWave') 
        ! 
        V0 = 0.0
        !
        ICA = 1e-8 
        !  
        bbb    = ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        dxb    = +2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        dtb    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))
        Kxx    = 0.0 
        Kyy    = -0.5*dtb  ! extrinsic curvature  
        Kzz    = +0.5*dtb  ! extrinsic curvature  
        !
        !
        V0(1)  = 1.0         ! \gamma_xx
        V0(4)  = 1.0 + bbb   ! \gamma_yy
        V0(6)  = 1.0 - bbb   ! \gamma_zz
        V0(7)  = Kxx         ! K_xx
        V0(10) = Kyy         ! K_yy
        V0(12) = Kzz         ! K_zz

        V0(17) = 1.0 
        !
        ! Auxiliary variables
        V0(24) = 0.0                        ! A_x  
        V0(36) = 0.0                        ! D_xxx
        V0(39) = +0.5*dxb                   ! D_xyy 
        V0(41) = -0.5*dxb                   ! D_xzz 
        !
    CASE('Z4Kerr2D','Z4Kerr3D') 
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, ddet )
        test = matmul(g_cov,g_contr) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables  
        betadown = MATMUL(g_cov, beta) 
        ! Derivative of beta down 
        dbetadown = 0.0 
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            dbetadown(k,i) = dbetadown(k,i) + 2.0*DD(k,i,j)*beta(j) + g_cov(i,j)*BB(k,j) 
          ENDDO
         ENDDO
        ENDDO 
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! covariant derivative of the shift beta 
        DO j = 1, 3
         DO i = 1, 3
           nablabeta(i,j) = dbetadown(i,j) 
           DO k = 1, 3
              nablabeta(i,j) = nablabeta(i,j) - Christoffel(i,j,k)*betadown(k) 
           ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        DO j = 1, 3
         DO i = 1, 3
             Kex(i,j) = 1.0/(2*alpha)*( nablabeta(i,j) + nablabeta(j,i) ) 
         ENDDO
        ENDDO
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = SUM( g_contr*Kex ) - 1.0/(alpha*fa)*SUM(beta*AA)   
        !
        ! The metric tensor 
        V0(1) = g_cov(1,1)    ! \gamma_11 
        V0(2) = g_cov(1,2)    ! \gamma_12
        V0(3) = g_cov(1,3)    ! \gamma_13
        V0(4) = g_cov(2,2)    ! \gamma_22
        V0(5) = g_cov(2,3)    ! \gamma_23
        V0(6) = g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature 
        V0(7)  = Kex(1,1) 
        V0(8)  = Kex(1,2) 
        V0(9)  = Kex(1,3) 
        V0(10) = Kex(2,2) 
        V0(11) = Kex(2,3) 
        V0(12) = Kex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Z1 
        V0(14) = 0.0          ! Z2 
        V0(15) = 0.0          ! Z3  
        V0(16) = 0.0          ! Theta 
        ! The lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! for debug only !! 
        !V0(21:23) = xGP(1:3) 
        !        
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        ! 
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(54) = K0 
        !
    CASE('Z4Kerr')
        !
       r      = xGP(1)   ! SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2)
       theta  = xGP(2)   ! 
       a      = 0.9  
       aom    = a 
       Mbh    = 1.0 
       M      = Mbh  
       r_hz   = Mbh + SQRT(Mbh**2 - aom**2) 
       !
       rho   = SQRT(r**2 + a**2*COS(theta)**2) 
       Delta = r**2 - 2*Mbh*r+a**2 
       Sig   = (r**2+a**2)**2 - a**2*Delta*SIN(theta)**2  
       zz    = 2.*M*r/rho**2
       !
       u0(:)  = 0.0
       !      
       g_cov( 1, 1:3) = (/ 1.0+zz,                     0.0,     -a*SIN(theta)**2*(1.0+zz)     /) 
       g_cov( 2, 1:3) = (/ 0.0,                        rho**2,   0.0                          /) 
       g_cov( 3, 1:3) = (/ -a*SIN(theta)**2*(1.0+zz),  0.0,      sig*SIN(theta)**2/rho**2   /) 
       !
       V0(1) = g_cov(1,1)    ! \gamma_11 
       V0(2) = g_cov(1,2)    ! \gamma_12
       V0(3) = g_cov(1,3)    ! \gamma_13
       V0(4) = g_cov(2,2)    ! \gamma_22
       V0(5) = g_cov(2,3)    ! \gamma_23
       V0(6) = g_cov(3,3)    ! \gamma_33 
       !
       ! the extrinsic curvature. general case for M >= 0, a >= 0  
       !
       Kex(1,1) = -2*(r**4+M*r**3-cos(theta)**2*M*a**2*r-a**4*cos(theta)**4)*M*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**6+3*a**2*cos(theta)**2*r**4+2*r**5*M+3*a**4*cos(theta)**4*r**2+4*r**3*cos(theta)**2*M*a**2+a**6*cos(theta)**6+2*cos(theta)**4*a**4*M*r)
       Kex(1,2) = 2*r*M*sin(theta)*cos(theta)*a**2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**4+2*a**2*cos(theta)**2*r**2+2*M*r**3+a**4*cos(theta)**4+2*cos(theta)**2*M*a**2*r)
       Kex(1,3) = sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*a*sin(theta)**2*M*(r**2-a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       Kex(2,1) = 2*r*M*sin(theta)*cos(theta)*a**2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**4+2*a**2*cos(theta)**2*r**2+2*M*r**3+a**4*cos(theta)**4+2*cos(theta)**2*M*a**2*r)
       Kex(2,2) = 2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*M*r**2/(r**2+a**2*cos(theta)**2+2*M*r)
       Kex(2,3) = -2*cos(theta)*r*M*sin(theta)**3*a**3*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**4+2*a**2*cos(theta)**2*r**2+2*M*r**3+a**4*cos(theta)**4+2*cos(theta)**2*M*a**2*r)
       Kex(3,1) = sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*a*sin(theta)**2*M*(r**2-a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       Kex(3,2) = -2*cos(theta)*r*M*sin(theta)**3*a**3*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**4+2*a**2*cos(theta)**2*r**2+2*M*r**3+a**4*cos(theta)**4+2*cos(theta)**2*M*a**2*r)
       Kex(3,3) = 2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*sin(theta)**2*(r*cos(theta)**4*a**4-cos(theta)**4*M*a**4+r**2*cos(theta)**2*M*a**2+2*cos(theta)**2*a**2*r**3+M*a**4*cos(theta)**2+r**5-M*a**2*r**2)*M*r/(r**6+3*a**2*cos(theta)**2*r**4+2*r**5*M+3*a**4*cos(theta)**4*r**2+4*r**3*cos(theta)**2*M*a**2+a**6*cos(theta)**6+2*cos(theta)**4*a**4*M*r)
       !        
       traceK = 2*(r**4+3*M*r**3+2*a**2*cos(theta)**2*r**2+cos(theta)**2*M*a**2*r+a**4*cos(theta)**4)*M*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**6+3*a**4*cos(theta)**4*r**2+a**6*cos(theta)**6+4*M**2*r**4+8*r**3*cos(theta)**2*M*a**2+4*cos(theta)**4*a**4*M*r+4*cos(theta)**2*M**2*a**2*r**2+4*r**5*M+3*a**2*cos(theta)**2*r**4)
       ! 
       V0(7)  = Kex(1,1) 
       V0(8)  = Kex(1,2) 
       V0(9)  = Kex(1,3) 
       V0(10) = Kex(2,2) 
       V0(11) = Kex(2,3) 
       V0(12) = Kex(3,3) 
       !
       ! The cleaning variables Z and Theta 
       !
       V0(13) = 0.0                                      ! Z1 
       V0(14) = 0.0                                      ! Z2 
       V0(15) = 0.0                                      ! Z3  
       V0(16) = 0.0                                      ! Theta 
      
       alpha = 1.0/SQRT(1+zz)                            ! alpha (lapse) 
       V0(17) = alpha 
       
       V0(18) = zz/(1+zz)                                ! beta_1 (shift) 
       V0(19) = 0.0                                      ! beta_2 
       V0(20) = 0.0                                      ! beta_3 
       !
       !! Auxiliary variables
       !! 
       !! vector A_i = \partial_i \alpha / \alpha 
       !!  
       V0(24) = -M*(-r**2+a**2*cos(theta)**2)/(a**4*cos(theta)**4+2*cos(theta)**2*r**2*a**2+2*cos(theta)**2*a**2*M*r+r**4+2*M*r**3)                     ! A1 = alpha_1/alpha  
       V0(25) = -2*a**2*M*r*cos(theta)*sin(theta)/(a**4*cos(theta)**4+2*cos(theta)**2*r**2*a**2+2*cos(theta)**2*a**2*M*r+r**4+2*M*r**3)                 ! A2 = alpha_2/alpha 
       V0(26) = 0.0                                                                                                                                     ! A3 = alpha_3/alpha 
       !
       ! Matrix B_ik = \partial_i \beta_k 
       !
       BMat = 0.0 
       ! special case a=0 
       !BMat(1,1) = -2.0/(r+2)**2 
       ! general case
       BMat(1,1) = -2*M*(r**2-a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+4*M*r**3+a**4*cos(theta)**4+4*cos(theta)**2*M*a**2*r+4*M**2*r**2)
       BMat(2,1) = 4*M*a**2*r*cos(theta)*sin(theta)/(r**4+2*a**2*cos(theta)**2*r**2+4*M*r**3+a**4*cos(theta)**4+4*cos(theta)**2*M*a**2*r+4*M**2*r**2)
       !
       V0(27) = BMat(1,1) 
       V0(28) = BMat(2,1)
       V0(29) = BMat(3,1)
       V0(30) = BMat(1,2)
       V0(31) = BMat(2,2)
       V0(32) = BMat(3,2)
       V0(33) = BMat(1,3)
       V0(34) = BMat(2,3)
       V0(35) = BMat(3,3) 
       !
       ! tensor D_ijk = \partial i \tilde \gamma_jk / 2 
       ! 
       DD(1,1,1) = M*(-r**2+a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(1,1,2) = 0
       DD(1,1,3) = -a*sin(theta)**2*M*(-r**2+a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(1,2,1) = 0
       DD(1,2,2) = r
       DD(1,2,3) = 0
       DD(1,3,1) = -a*sin(theta)**2*M*(-r**2+a**2*cos(theta)**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(1,3,2) = 0
       DD(1,3,3) = sin(theta)**2*(cos(theta)**4*r*a**4-cos(theta)**4*a**4*M+2*r**3*a**2*cos(theta)**2+a**4*M*cos(theta)**2+cos(theta)**2*a**2*M*r**2+r**5-a**2*M*r**2)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(2,1,1) = 2*a**2*M*r*cos(theta)*sin(theta)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(2,1,2) = 0
       DD(2,1,3) = -(a**4*cos(theta)**4+2*a**2*cos(theta)**2*r**2+r**4+2*M*r**3+2*a**2*M*r)*a*cos(theta)*sin(theta)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(2,2,1) = 0
       DD(2,2,2) = -a**2*cos(theta)*sin(theta)
       DD(2,2,3) = 0
       DD(2,3,1) = -(a**4*cos(theta)**4+2*a**2*cos(theta)**2*r**2+r**4+2*M*r**3+2*a**2*M*r)*a*cos(theta)*sin(theta)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(2,3,2) = 0
       DD(2,3,3) = (2*r**4*a**2*cos(theta)**2+2*r**2*a**4*cos(theta)**2+r**6+r**4*a**2+4*a**2*M*r**3+2*a**4*M*r-2*a**4*cos(theta)**4*M*r+cos(theta)**4*a**6+r**2*a**4*cos(theta)**4-4*r**3*a**2*cos(theta)**2*M)*cos(theta)*sin(theta)/(r**4+2*a**2*cos(theta)**2*r**2+a**4*cos(theta)**4)
       DD(3,1,1) = 0
       DD(3,1,2) = 0
       DD(3,1,3) = 0
       DD(3,2,1) = 0
       DD(3,2,2) = 0
       DD(3,2,3) = 0
       DD(3,3,1) = 0
       DD(3,3,2) = 0
       DD(3,3,3) = 0
       
       !    
       V0(36) = DD(1,1,1) 
       V0(37) = DD(1,1,2) 
       V0(38) = DD(1,1,3) 
       V0(39) = DD(1,2,2) 
       V0(40) = DD(1,2,3) 
       V0(41) = DD(1,3,3) 
       V0(42) = DD(2,1,1) 
       V0(43) = DD(2,1,2) 
       V0(44) = DD(2,1,3) 
       V0(45) = DD(2,2,2) 
       V0(46) = DD(2,2,3) 
       V0(47) = DD(2,3,3) 
       V0(48) = DD(3,1,1) 
       V0(49) = DD(3,1,2) 
       V0(50) = DD(3,1,3) 
       V0(51) = DD(3,2,2) 
       V0(52) = DD(3,2,3) 
       V0(53) = DD(3,3,3) 
       !
      SELECT CASE(EQN%CCZ4LapseType) 
      CASE(0)  ! harmonic 
          fa = 1.0 
      CASE DEFAULT  ! 1 + log 
          fa = 2.0/alpha
      END SELECT   
       !
       ggg = EQN%CCZ4g 
       fff = EQN%CCZ4f 
       !
       !K0 = (2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*r**4-r**3*M+6*r**3*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*M+4*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*a**2*cos(theta)**2*r**2+2*r*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*M*a**2*cos(theta)**2+a**2*cos(theta)**2*M*r+2*sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))*a**4*cos(theta)**4)*M/(a**6*cos(theta)**6+r**6+4*cos(theta)**4*a**4*M*r+8*cos(theta)**2*a**2*r**3*M+4*cos(theta)**2*a**2*M**2*r**2+4*r**4*M**2+3*cos(theta)**4*a**4*r**2+3*r**4*a**2*cos(theta)**2+4*r**5*M)
       ! 
       K0 = 2*(r**4*fa-r**3*M+3*r**3*fa*M+2*fa*a**2*cos(theta)**2*r**2+r*fa*a**2*cos(theta)**2*M+a**2*cos(theta)**2*M*r+fa*a**4*cos(theta)**4)*M/sqrt((r**2+a**2*cos(theta)**2+2*M*r)/(r**2+a**2*cos(theta)**2))/(r**2+a**2*cos(theta)**2+2*M*r)/(r**4+2*r**2*a**2*cos(theta)**2+a**4*cos(theta)**4)/fa        
       !         
       V0(54) = K0 
        !
    CASE('Z4TwoPunctures')
        !
#ifdef TWOPUNCTURES  
        CALL TwoPunctures_Interpolate(xGP, V0) 
#else
        PRINT *, ' TwoPunctures not available. Please compile with -DTWOPUNCTURES flag. '
#endif 
        !
    CASE DEFAULT
        V0(:) = 1.0 
        V0(1) = V0(1) + 0.1*SIN(2*Pi*xGP(1)-tGP)     
    END SELECT 
#endif 
    !
#if defined(CCZ4EINSTEIN) || defined(BSSZ4EINSTEIN) 
    SELECT CASE(ICType)
    CASE('CCZ4GaugeWave')
        V0 = 0.0
        !
        ICA = 0.1 
        !ICA = 0.9  
        ! 
        ! Gauge wave propagating along x 
        HH     = 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        dxH    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        dxphi  = - HH**(-7.0/6.0)*dxH/6.0
        phi    = ( 1.0 / HH)**(1.0/6.0)
        Kxx    = - Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))/SQRT( 1.0 - ICA*SIN( 2.0*Pi*( xGP(1) - tGP))  )  ! extrinsic curvature  
        traceK = Kxx/HH
        !
!        HH = 1.0 
!        phi = 1.0 
!        traceK = 0.0
!        Kxx = 0.0 
!        traceK = 0.0 
!        dxH = 0.0
!        dxphi = 0.0 
        !
        V0(:)  = 0.0
        V0(1)  = phi**2*HH                          ! \tilde\gamma_xx
        V0(4)  = phi**2                             ! \tilde\gamma_yy
        V0(6)  = phi**2                             ! \tilde\gamma_zz
        V0(7)  = phi**2*(Kxx - 1.0/3.0*traceK*HH )  ! \tilde A_{xx}
        V0(10) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_{yy}
        V0(12) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_[zz}
        V0(17) = SQRT(HH)
        V0(14) = 2.0/(3.0*HH**(5.0/3.0))*dxH        ! Gtilde 
        !
        ! Auxiliary variables
        V0(24) = 1.0/(2.0*HH)*dxH                ! A_x  
        V0(36) = HH**(-1.0/3.0)*dxH/3.0          ! D_xxx
        V0(39) = phi*dxphi                       ! D_xyy
        V0(41) = phi*dxphi                       ! D_xzz
        ! 
        V0(54) = traceK
        V0(55) = phi
        V0(56) = dxphi/phi                       ! P_x
    CASE('CCZ4LinearWave')        
        !
        V0 = 0.0
        !
        ICA = 1e-8   
        !        
        ! The amplitude ICA is read from the parameter file: linear wave along x direction        
        bbb    = ICA*SIN( 2.0*Pi*( xGP(1) - tGP)) 
        !
        dxb    = +2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP)) 
        dtb    = -2.0*Pi*ICA*COS( 2.0 * Pi*(xGP(1) - tGP))
        ! 
        Kxx    =  0.0                            ! extrinsic curvature  
        Kyy    = -0.5*dtb                        ! extrinsic curvature  
        Kzz    = +0.5*dtb                        ! extrinsic curvature 
        !
        ! Compute the metric and its derivatives 
        alpha       = 1.0 
        beta        = 0.0
        g_cov       = 0.0 
        g_cov(1,1)  = 1.0
        g_cov(2,2)  = 1.0+bbb
        g_cov(3,3)  = 1.0-bbb 
        g_contr     = 0.0
        g_contr(1,1) = 1./g_cov(1,1) 
        g_contr(2,2) = 1./g_cov(2,2) 
        g_contr(3,3) = 1./g_cov(3,3) 
        !
        AA    = 0.0
        BB    = 0.0 
        DD    = 0.0
        DD(1,1,1) =  0.0 
        DD(1,2,2) =  0.5*dxb 
        DD(1,3,3) = -0.5*dxb 
        !
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3) 
        phi = detg**(-1./6.)   
        PP    = 0.0
        PP(1) = 2.0/3.0*ICA**2*cos(2*pi*(xGP(1)-tGP))*pi*sin(2*pi*(xGP(1)-tGP))/(1-ICA**2+ICA**2*cos(2*pi*(xGP(1)-tGP))**2) 
        !        
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        Kex(1,1) = Kxx 
        Kex(2,2) = Kyy 
        Kex(3,3) = Kzz 
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3.*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        !     
        K0 = 0.0 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(59) = K0 
        !
    CASE('CCZ4GowdyWave')        
        !
        V0 = 0.0 
        !
        ICA = 1e-8   
        !        
        ! The amplitude ICA is read from the parameter file: linear wave along x direction        
        CALL JYNA(1,2*Pi*tGP,nm,BJ,DJ,BY,DY) 
        CALL JYNA(1,2*Pi    ,nm,IJ,DIJ,IY,DIY) 
        !  
        bbb    =       BJ(0)*COS(2*Pi*xGP(1)) 
        dtb    =  2*Pi*DJ(0)*COS(2*Pi*xGP(1))  
        dxb    = -2*Pi*BJ(0)*SIN(2*Pi*xGP(1))  
        !
        lambda  = -2*Pi*tGP*BJ(0)*BJ(1)*COS(2*Pi*xGP(1))**2 + 2*Pi**2*tGP**2*(BJ(0)**2+BJ(1)**2) -0.5*( 4*Pi**2*(IJ(0)**2+IJ(1)**2) - 2*Pi*IJ(0)*IJ(1) ) 
        lambdat = tGP*(dtb**2 + dxb**2) 
        lambdax = 2*tGP*dxb*dtb 
        ! 
        Kxx    =  0.25*tGP**(-0.25)*EXP(+0.25*lambda)*(1./tGP - lambdat)          ! extrinsic curvature  
        Kyy    = -0.50*tGP**(+0.25)*EXP(-0.25*lambda)*EXP(+bbb)*(1.0 + tGP*dtb)   ! extrinsic curvature  
        Kzz    = -0.50*tGP**(+0.25)*EXP(-0.25*lambda)*EXP(-bbb)*(1.0 - tGP*dtb)   ! extrinsic curvature 
        !
        ! Compute the metric and its derivatives 
        alpha       = tGP**(-0.25) * EXP(0.25*lambda) 
        beta        = 0.0
        g_cov       = 0.0 
        g_cov(1,1)  = tGP**(-0.5)*EXP(0.5*lambda)         
        g_cov(2,2)  = tGP*EXP( bbb) 
        g_cov(3,3)  = tGP*EXP(-bbb) 
        g_contr     = 0.0
        g_contr(1,1) = 1./g_cov(1,1) 
        g_contr(2,2) = 1./g_cov(2,2) 
        g_contr(3,3) = 1./g_cov(3,3) 
        !
        AA    = 0.0
        AA(1) = 1./4.*lambdax 
        BB    = 0.0 
        DD    = 0.0
        DD(1,1,1) =  0.0 
        DD(1,2,2) =  0.5*tGP*EXP( bbb)*dxb  
        DD(1,3,3) = -0.5*tGP*EXP(-bbb)*dxb   
        !
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3) 
        phi = detg**(-1./6.)   
        PP    = 0.0
        PP(1) = -1./12.*lambdax 
        !        
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        Kex(1,1) = Kxx 
        Kex(2,2) = Kyy 
        Kex(3,3) = Kzz 
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3.*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        !     
        K0 = 0.0 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(59) = K0         
        !
    CASE('CCZ4Puncture')
        V0 = 0.0
        !
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = 0.0 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(59) = K0 
        ! 
    CASE('CCZ4TwoPunctures') 
        !
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature 
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3)         
        !
        CONTINUE 
        !         
    END SELECT
#endif 
    !
#ifdef SRMHD
    SELECT CASE(ICType) 
    CASE('RMHDAlfvenWave') 
       rho0 = 1.
       p0   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + EQN%gamma / ( EQN%gamma - 1.0) * p0 / rho0
       tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2) 
       !
       BV(1) = B0
       BV(2) = eta * B0 * COS( 2*Pi*( xGP(1) - vax*tGP ) )
       BV(3) = eta * B0 * SIN( 2*Pi*( xGP(1) - vax*tGP ) )
       !
       VV(1)   = 0.0
       VV(2:3) = - vax * BV(2:3) / B0
       !
       ! Now convert to conservative variables 
       !
       V0 = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /)        
       !
       CONTINUE
       !
    END SELECT
    
#endif 
    !
#ifdef GRMHD
    SELECT CASE(ICType)
       !
#ifdef TORUS        
    CASE('GRMHDTorus')
       xloc = xGP 
       !
       minR = 0.8*ExcisionRadius
       !
#ifdef Spherical  
       !
       r=xloc(1)
       IF ( r .LT. minR) THEN ! RETURN
           r = minR
           xloc(1) = r
       ENDIF
       !
#else
       !
       r=SQRT(SUM(xloc(1:nDim)**2))
       !
       IF ( r .LT. minR) THEN ! RETURN
           r=minR
           xloc(1) = r*DSIN(THETA)*DCOS(PHI2)
           xloc(2) = r*DSIN(THETA)*DSIN(PHI2)
           xloc(3) = r*DCOS(THETA) 
       ENDIF
       !
#endif 
       !
       CALL METRIC ( xloc , lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(EQN%gamma - 1.0)

       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
       !
       CALL MatrixInverse3x3(g_cov,g_contr,gp)
       !
#ifdef Spherical  
           !
           CALL init_disk_isentropic(xloc, V0(1:8),METRIC_KSS) 
           V0(9) = 0.
           V0(10:19) = (/  lapse, shift(1:3), gammaij(1:6)/)
           ! 
#else  
       ! The Following is for Kerr-Schild Cartesian coordinates 
       r      = SQRT( xloc(1)**2 + xloc(2)**2 + xloc(3)**2) 
       theta  = ACOS( xloc(3)/r)
       phi    = ATAN2( xloc(2), xloc(1))
       ! 
       minR = 0.8*ExcisionRadius
       IF(r.LT.minR) THEN  ! RETURN
           r=ExcisionRadius
           xloc(3) = r*DCOS(theta)
           xloc(2) = r*DSIN(theta)*DSIN(phi)
           xloc(1) = r*DSIN(theta)*DCOS(phi)
       ENDIF 
       !
       xGP_sph = (/r,theta, phi   /) 
       !
       CALL init_disk_isentropic(xGP_sph, V0(1:8),METRIC_KSS)
       !    
        CALL Cart2SphMatrix_cov(A_coord,xGP_sph)
        ! 
        CALL MatrixInverse3x3(A_coord,iA_coord,detA)      ! iA = Sph2CartMatrix_cov    = Cart2SphMatrix_contr
        A_coord_contr(:,:) = TRANSPOSE(iA_coord(:,:))     ! Cart2SphMatrix_contr = TRANSPOSE( Sph2CartMatrix_cov )
        iA_coord_contr(:,:) = TRANSPOSE(A_coord(:,:))    ! Sph2CartMatrix_contr = TRANSPOSE( Cart2SphMatrix_cov )
        !
        VV_cov(1:3) = MATMUL(iA_coord,V0(2:4))  ! Spherical to Cartesian: Velocity is covariant
        BV_contr(1:3) = MATMUL(iA_coord_contr,V0(6:8))  !  Spherical to Cartesian:  Magnetic field is contravariant  
        !
        V0(2:4)=VV_cov(1:3)
        V0(6:8)=BV_contr(1:3)
        !   
        V0(10:19) = (/  lapse, shift(1:3), gammaij(1:6)/)
        !
        IF(MAXVAL(ABS(V0(2:4))).GT.1e-10) THEN
            continue
        ENDIF
        IF(MAXVAL(ABS(V0(10:19))).LT.1e-10) THEN
            continue
        ENDIF
        ! 
#endif        
#endif 
       
    CASE('GRMHDAlfven')
       rho = 1.
       p   = 1.
       eta  = 1. 
       B0   = 1.0  
       !
       hh = 1.0 + EQN%gamma / ( EQN%gamma - 1.0) * p / rho
       tempaa = rho * hh + B0**2 * ( 1.0 + eta**2)
       tempab = 2.0 * eta * B0**2 / tempaa
       tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
       va2 = b0**2 / ( tempaa * tempac)
       vax = sqrt ( va2)       
       !
       ! flat metric: contr and cov are the same in Cartesian
       BV_contr(1) = B0
       BV_contr(2) = eta * B0 * COS( xGP(1) - vax*tGP)
       BV_contr(3) = eta * B0 * SIN( xGP(1) - vax*tGP)
       !
       VV_cov(1)   = 0.0
       VV_cov(2:3) = - vax * BV_contr(2:3) / B0
       !
       ! Now convert to conservative variables 
       !
       lapse = 1.0
       shift(1:3) = 0.
       gammaij = 0.
       gammaij(1) = 1.
       gammaij(4) = 1.
       gammaij(6) = 1. 
       !
       V0 = (/ rho, VV_cov(1:3), p ,BV_contr(1:3), 0.0, lapse, shift(1:3), gammaij(1:6)/) 
       !V0 = (/ 1.0, 0.0, 0.0, 0.0, 1.0 ,0.0, 0.0, 0.0, 0.0, lapse, shift(1:3), gammaij(1:6)/)  
       ! 
    CASE('GRMHDAccretion')       
       !aom = 0.0
       !Mbh = 1.0  
       rc   = 8.0 
       rhoc = 1./16.
       B0 = 4. 
       xloc = xGP
       IF ( aom > 0.0) THEN 
           WRITE(*,*)'Spherical Accretion solution is only for a = 0!'
           STOP
       ENDIF 
       !
       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2) 
       theta  = ACOS( xGP(3)/r)
       phi2    = ATAN2( xGP(2), xGP(1))
       !
       !IF ( r .LT. 0.0001) THEN ! RETURN
       !    lapse = 1.0
       !    shift(1:3) = 0.
       !    !
       !    rho = 1.0
       !    p = 1.0
       !    VV_cov = 0.
       !    BV(1:3) = 0.
       !    !
       !    !shift_contr = MATMUL(g_contr,shift)  !shift is controvariant. See fluxes....
       !    !
       !    gammaij = 0.
       !    gammaij(1) = 1.0
       !    gammaij(4) = 1.0
       !    gammaij(6) = 1.0
       !    !  
       !    V0 = (/ rho, VV_cov(1:3), p ,BV(1:3), 0.0, lapse, shift(1:3), gammaij(1:6)/)
       !    !
       !    CALL PDEPrim2Cons(u0,V0)
       !    RETURN
       !ELSE
       !minR = MAX(MAXVAL(dx),0.5*ExcisionRadius) 
       minR = 0.8*ExcisionRadius
       IF ( r .LT. minR) THEN ! RETURN
           r = minR
           xloc(1) = r*DSIN(THETA)*DCOS(PHI2)
           xloc(2) = r*DSIN(THETA)*DSIN(PHI2)
           xloc(3) = r*DCOS(THETA) 
       ENDIF
       !
       CALL METRIC ( xloc , lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
       ng     = 1.0/(EQN%gamma - 1.0)  
       !
       !
       !
       !phi2    = ACOS( xGP(1) / (r*SIN(theta)))
       zz     = 2.0/r   ! we are computing the solution at theta=pi/2
       betaru = zz/(1.0 + zz)
       g_tt   = zz - 1.0
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
       
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON   
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dtt  = -f/df
          IF (abs(dtt) < 1.e-14) EXIT
          tt = tt + dtt
       ENDDO
       IF(iNewton.gE.MAXNEWTON) THEN
           continue
       ENDIF
       !
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = DSIN(theta)*DCOS(phi2)*vr
       vy  = DSIN(theta)*DSIN(phi2)*vr
       vz  = DCOS(theta)*vr
       !
       VV(1:3) = (/ vx, vy, vz /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       BV(1:3) = 0.
       !
       lapse = lapse
       !
       !shift_contr = MATMUL(g_contr,shift)  !shift is controvariant. See fluxes....
       !
       gammaij(1) = g_cov(1,1)
       gammaij(2) = g_cov(1,2)
       gammaij(3) = g_cov(1,3)
       gammaij(4) = g_cov(2,2)
       gammaij(5) = g_cov(2,3)
       gammaij(6) = g_cov(3,3) 
       !
       SELECT CASE(TRIM(ICType2))
       CASE('GRMHD-BondiAccretion')
           !
#ifdef Spherical
           PRINT *, 'GRMHD-BondiAccretio only for Cartesian coordinates!'
           continue
           !
           STOP
#endif           
           !
           BV_contr(1:3) = 0.
           !
           IF(B0.GT.0) THEN
                BV_contr(3) = 2.2688/Mbh*SQRT(B0)
           ENDIF
           !
           BV_contr(1) = gm*BV_contr(3)*Mbh**2/r**3*xloc(1)
           BV_contr(2) = gm*BV_contr(3)*Mbh**2/r**3*xloc(2)
           BV_contr(3) = gm*BV_contr(3)*Mbh**2/r**3*xloc(3)
           !
           !BV(1:3) = MATMUL(g_cov,BV_contr(1:3))
           ! 
           !
           ! 
       CASE('GRMHDCBS') 
           !
#ifndef Spherical
           PRINT *, 'GRMHDCBS only for spherical coordinates!'
           continue
           !
           STOP
#endif           
           !
           BV_contr(1:3) = 0.
           BV_contr(1) = lapse*B0*DCOS(xGP(2))
           BV_contr(2) = -  lapse*B0*DSIN(xGP(2))/xGP(1) 
           !
           !BV(1:3) = MATMUL(g_cov,BV(1:3))
           ! 
           !
           CASE DEFAULT
           !
           continue
           !
       END SELECT
       !  
       V0 = (/ rho, VV_cov(1:3), p ,BV_contr(1:3), 0.0, lapse, shift(1:3), gammaij(1:6)/) 
       !
       
       IF ( r .LT. 0.5) THEN ! RETURN
           continue
       ENDIF
       CONTINUE
       !
    END SELECT 
#endif
    !
#ifdef Z4GRMHD
    SELECT CASE(ICType) 
    CASE('Z4GRMHDAccretion')
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, ddet )
        test = matmul(g_cov,g_contr) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables  
        betadown = MATMUL(g_cov, beta) 
        ! Derivative of beta down 
        dbetadown = 0.0 
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            dbetadown(k,i) = dbetadown(k,i) + 2.0*DD(k,i,j)*beta(j) + g_cov(i,j)*BB(k,j) 
          ENDDO
         ENDDO
        ENDDO 
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! covariant derivative of the shift beta 
        DO j = 1, 3
         DO i = 1, 3
           nablabeta(i,j) = dbetadown(i,j) 
           DO k = 1, 3
              nablabeta(i,j) = nablabeta(i,j) - Christoffel(i,j,k)*betadown(k) 
           ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        DO j = 1, 3
         DO i = 1, 3
             Kex(i,j) = 1.0/(2*alpha)*( nablabeta(i,j) + nablabeta(j,i) ) 
         ENDDO
        ENDDO
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = SUM( g_contr*Kex ) - 1.0/(alpha*fa)*SUM(beta*AA)   
        !
        ! The metric tensor 
        V0(1) = g_cov(1,1)    ! \gamma_11 
        V0(2) = g_cov(1,2)    ! \gamma_12
        V0(3) = g_cov(1,3)    ! \gamma_13
        V0(4) = g_cov(2,2)    ! \gamma_22
        V0(5) = g_cov(2,3)    ! \gamma_23
        V0(6) = g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature 
        V0(7)  = Kex(1,1) 
        V0(8)  = Kex(1,2) 
        V0(9)  = Kex(1,3) 
        V0(10) = Kex(2,2) 
        V0(11) = Kex(2,3) 
        V0(12) = Kex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Z1 
        V0(14) = 0.0          ! Z2 
        V0(15) = 0.0          ! Z3  
        V0(16) = 0.0          ! Theta 
        ! The lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        ! 
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(54) = K0 
        !
       aom = 0.0
       Mbh = 1.0  
       rc   = 8.0 
       rhoc = 1./16. 

       IF ( aom > 0.0) THEN 
           WRITE(*,*)'Spherical Accretion solution is only for a = 0!'
           STOP
       ENDIF 
       !
       lapse = alpha
       shift = beta 
       ng     = 1.0/(EQN%gamma - 1.0)

       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = MAX(0.1, SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2) ) 
       !
       theta  = ACOS( xGP(3)/r)
       !
       phi    = ATAN2( xGP(2), xGP(1))
       !phi    = ACOS( xGP(1) / (r*SIN(theta)))
       zz     = 2.0/r   ! we are computing the solution at theta=pi/2
       betaru = zz/(1.0 + zz)
       g_tt   = zz - 1.0
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
       
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON   
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dtt  = -f/df
          IF (abs(dtt) < 1.e-10) EXIT
          tt = tt + dtt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = SIN(theta)*COS(phi)*vr
       vy  = SIN(theta)*SIN(phi)*vr
       vz  = COS(theta)*         vr
       !
       VV(1:3) = (/ vx, vy, vz /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       BV(1:3) = 0.
       !
       V0(55:63) = (/ rho, VV_cov(1:3), p, BV(1:3), 0.0 /)
       !
    END SELECT 
#endif 
    !
#if defined(CCZ4GRMHD) || defined(CCZ4EINSTEIN) || defined(BSSZ4GRMHD) || defined(BSSZ4EINSTEIN) 
    SELECT CASE(ICType) 
    CASE('CCZ4Minkowski')
       CALL RANDOM_NUMBER(V0) 
       !V0 = V0*1e-4   
       V0 = 0e-6 
       !V0(1:6) = 0.0 
       ! 
       V0(1)  = V0(1)  + 1.0 
       V0(4)  = V0(4)  + 1.0 
       V0(6)  = V0(6)  + 1.0 
       V0(17) = V0(17) + 1.0 
       V0(55) = V0(55) + 1.0 
       !
       ! V0(13) = 0.1*EXP(-0.5*( xGP(1)**2+xGP(2)**2)/0.1**2 ) 
       !
#ifdef CCZ4GRMHD 
       V0(60:68) = V0(60:68) + (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 /) 
#endif         
       ! 
    END SELECT 
#endif  
    !
#if defined(CCZ4GRMHD) || defined(BSSZ4GRMHD) || defined(CCZ4EINSTEIN)  
    SELECT CASE(ICType) 
    CASE('CCZ4GRMHDAccretion','CCZ4Kerr2D','CCZ4Kerr3D')
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        betadown = MATMUL(g_cov, beta) 
        ! Derivative of beta down 
        dbetadown = 0.0 
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            dbetadown(k,i) = dbetadown(k,i) + 2.0*DD(k,i,j)*beta(j) + g_cov(i,j)*BB(k,j) 
          ENDDO
         ENDDO
        ENDDO 
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! covariant derivative of the shift beta 
        DO j = 1, 3
         DO i = 1, 3
           nablabeta(i,j) = dbetadown(i,j) 
           DO k = 1, 3
              nablabeta(i,j) = nablabeta(i,j) - Christoffel(i,j,k)*betadown(k) 
           ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        DO j = 1, 3
         DO i = 1, 3
             Kex(i,j) = 1.0/(2*alpha)*( nablabeta(i,j) + nablabeta(j,i) ) 
         ENDDO
        ENDDO
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = traceK - 1.0/(alpha*fa)*SUM(beta*AA)   
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 

#if defined(CCZ4GRMHD) || defined(CCZ4EINSTEIN) 
        V0(59) = K0 
#endif        
#ifdef BSSZ4GRMHD         
        V0(62) = K0 
#endif        
        !
#ifdef CCZ4EINSTEIN 
       CALL PDEPrim2Cons(u0,V0)  
       RETURN
#endif 
     
        !
       aom = 0.0
       Mbh = 1.0  
       rc   = 8.0 
       rhoc = 1./16. 

       IF ( aom > 0.0) THEN 
           WRITE(*,*)'Spherical Accretion solution is only for a = 0!'
           STOP
       ENDIF 
       !
       lapse = alpha
       shift = beta 
       ng     = 1.0/(EQN%gamma - 1.0)

       ! The Following is for Kerr-Schild Cartesian coordinates
       r      = MAX(0.1, SQRT( xGP(1)**2 + xGP(2)**2 + xGP(3)**2) ) 
       !
       theta  = ACOS( xGP(3)/r)
       !
       phi    = ATAN2( xGP(2), xGP(1))
       !phi    = ACOS( xGP(1) / (r*SIN(theta)))
       zz     = 2.0/r   ! we are computing the solution at theta=pi/2
       betaru = zz/(1.0 + zz)
       g_tt   = zz - 1.0
       !
       urc = sqrt(1.0 / (2.0*rc))
       vc2 = urc**2 / (1.0 - 3.0*urc**2)
       tc  = ng*vc2 / ((1.0 + ng)*(1.0 - ng*vc2))
       pc  = rhoc*tc
       
       c1 = urc*tc**ng*rc**2
       c2 = (1.0 + ( 1.0 + ng)*tc)**2*(1.0 - 2.0/rc+urc**2)
       !
       tt = tc
       DO iNewton = 1, MAXNEWTON   
          urr = c1 / (r**2*tt**ng)
          f   = (1.0 + (1.0 + ng)*tt)**2*(1.0 - 2.0/r + urr**2) - c2
          df  = 2.0 * (1.0 + ng)*(1.0 + (1.0 + ng)*tt)*(1.0 - 2.0/r + urr**2) - 2.0*ng*urr**2/tt*(1.0 + (1.0 + ng)*tt)**2
          dtt  = -f/df
          IF (abs(dtt) < 1.e-10) EXIT
          tt = tt + dtt
       ENDDO
       ut     = (-zz*urr + sqrt(urr**2 - zz + 1.0))/(zz - 1.0)
       LF     = lapse*ut
       vr     = ( urr / LF + betaru / lapse)
       vtheta = 0.0
       vphi   = 0.0
       !
       vx  = SIN(theta)*COS(phi)*vr
       vy  = SIN(theta)*SIN(phi)*vr
       vz  = COS(theta)*         vr
       !
       VV(1:3) = (/ vx, vy, vz /)
       ! Convert to covariant velocities
       VV_cov = MATMUL(g_cov, VV)
       !
       rho = rhoc*(tt/tc)**ng
       p   = rho*tt
       !
       BV(1:3) = 0.
       !
#ifdef CCZ4GRMHD         
       V0(60:68) = (/ rho, VV_cov(1:3), p, BV(1:3), 0.0 /)
#endif
#ifdef BSSZ4GRMHD         
       V0(63:71) = (/ rho, VV_cov(1:3), p, BV(1:3), 0.0 /)
#endif
       !
    CASE('CCZ4RNSID')
        ! 
        V0 = 0.0
        ! Compute the metric and its derivatives 
        CALL METRIC(  xGP, alpha,  gp, gm, beta, Kex, g_cov, g_contr, phi )
        CALL DMETRIC( xGP, AA, BB, DD, PP )
        test = matmul(g_cov,g_contr)
        detg = g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)
        !phi = detg**(-1./6.) 
        AA = AA/alpha   ! convert pure derivatives to aux. variables 
        DD = 0.5*DD     ! convert pure derivatives to aux. variables 
        PP = PP/phi     ! convert pure derivatives to aux. variables  
        betadown = MATMUL(g_cov, beta) 
        ! Derivative of beta down 
        dbetadown = 0.0 
        DO j = 1, 3
         DO i = 1, 3
          DO k = 1, 3
            dbetadown(k,i) = dbetadown(k,i) + 2.0*DD(k,i,j)*beta(j) + g_cov(i,j)*BB(k,j) 
          ENDDO
         ENDDO
        ENDDO 
        ! Christoffel symbols 
        Christoffel = 0.0
        DO k = 1, 3
         DO j = 1, 3
          DO i = 1, 3
              DO l = 1, 3 
                Christoffel(i,j,k) = Christoffel(i,j,k) + g_contr(k,l)*( DD(i,j,l) + DD(j,i,l) - DD(l,i,j) ) 
              ENDDO 
          ENDDO
         ENDDO
        ENDDO 
        ! covariant derivative of the shift beta 
        DO j = 1, 3
         DO i = 1, 3
           nablabeta(i,j) = dbetadown(i,j) 
           DO k = 1, 3
              nablabeta(i,j) = nablabeta(i,j) - Christoffel(i,j,k)*betadown(k) 
           ENDDO
         ENDDO
        ENDDO 
        ! Extrinsic curvature 
        Kex = 0.0
        DO j = 1, 3
         DO i = 1, 3
             Kex(i,j) = 1.0/(2*alpha)*( nablabeta(i,j) + nablabeta(j,i) ) 
         ENDDO
        ENDDO
        ! Trace of K 
        traceK = SUM( Kex*g_contr ) 
        ! The conformal traceless part of the extrinsic curvature
        Aex = phi**2*( Kex - 1./3*traceK*g_cov ) 
        ! The contracted connection coefficients 
        Gtilde = 0.0
        DO i = 1, 3
         DO j = 1, 3
          DO k = 1, 3
           DO l = 1, 3
               Gtilde(i) = Gtilde(i) + 1./phi**2*( g_contr(i,j)*g_contr(k,l)*( 2*DD(l,j,k) + 2*PP(l)*g_cov(j,k) )  ) 
           ENDDO
          ENDDO 
         ENDDO
        ENDDO 
        SELECT CASE(EQN%CCZ4LapseType) 
        CASE(0)  ! harmonic 
           fa = 1.0 
        CASE DEFAULT  ! 1 + log 
           fa = 2.0/alpha
        END SELECT 
        ! K0 to make the PDE for alpha stationary 
        K0 = traceK - 1.0/(alpha*fa)*SUM(beta*AA)   
        !
        ! The metric tensor 
        V0(1) = phi**2*g_cov(1,1)    ! \gamma_11 
        V0(2) = phi**2*g_cov(1,2)    ! \gamma_12
        V0(3) = phi**2*g_cov(1,3)    ! \gamma_13
        V0(4) = phi**2*g_cov(2,2)    ! \gamma_22
        V0(5) = phi**2*g_cov(2,3)    ! \gamma_23
        V0(6) = phi**2*g_cov(3,3)    ! \gamma_33 
        ! The extrinsic curvature         
        V0(7)  = Aex(1,1)  
        V0(8)  = Aex(1,2) 
        V0(9)  = Aex(1,3) 
        V0(10) = Aex(2,2) 
        V0(11) = Aex(2,3) 
        V0(12) = Aex(3,3) 
        ! The cleaning variables Z and Theta 
        V0(13) = 0.0          ! Theta 
        V0(14) = Gtilde(1)    ! G1 
        V0(15) = Gtilde(2)    ! G2 
        V0(16) = Gtilde(3)    ! G3  
        ! The logarithm of the lapse 
        V0(17) = alpha 
        ! The shift 
        V0(18) = beta(1)  
        V0(19) = beta(2) 
        V0(20) = beta(3) 
        !
        ! Auxiliary variables
        ! 
        !! vector A_i = \partial_i \alpha / \alpha 
        !!  
        V0(24) = AA(1)     ! A1 = alpha_1/alpha  
        V0(25) = AA(2)     ! A2 = alpha_2/alpha    
        V0(26) = AA(3)     ! A3 = alpha_3/alpha   
        !
        ! Matrix B_ik = \partial_i \beta_k 
        !
        V0(27) = BB(1,1) 
        V0(28) = BB(2,1)
        V0(29) = BB(3,1)
        V0(30) = BB(1,2)
        V0(31) = BB(2,2)
        V0(32) = BB(3,2)
        V0(33) = BB(1,3)
        V0(34) = BB(2,3)
        V0(35) = BB(3,3) 
        !
        ! tensor D_ijk = 0.5 \partial i \tilde \gamma_jk   
        !
        DO j = 1, 3
         DO i = 1, 3 
          DO k = 1, 3 
            DD(k,i,j) = phi**2*( DD(k,i,j) + PP(k)*g_cov(i,j) ) 
          ENDDO
         ENDDO
        ENDDO 
        !
        V0(36) = DD(1,1,1) 
        V0(37) = DD(1,1,2) 
        V0(38) = DD(1,1,3) 
        V0(39) = DD(1,2,2) 
        V0(40) = DD(1,2,3) 
        V0(41) = DD(1,3,3) 
        V0(42) = DD(2,1,1) 
        V0(43) = DD(2,1,2) 
        V0(44) = DD(2,1,3) 
        V0(45) = DD(2,2,2) 
        V0(46) = DD(2,2,3) 
        V0(47) = DD(2,3,3) 
        V0(48) = DD(3,1,1) 
        V0(49) = DD(3,1,2) 
        V0(50) = DD(3,1,3) 
        V0(51) = DD(3,2,2) 
        V0(52) = DD(3,2,3) 
        V0(53) = DD(3,3,3) 
        ! trace of K 
        V0(54) = traceK 
        ! logarithm of the conformal factor phi 
        V0(55) = phi 
        ! derivative of phi 
        V0(56) = PP(1) 
        V0(57) = PP(2) 
        V0(58) = PP(3) 
        ! The trace of the extrinsic curvature at the initial time 
        V0(59) = K0 
        !
#ifdef CCZ4EINSTEIN 
       CALL PDEPrim2Cons(u0,V0)  
       RETURN
#endif      
        !
#ifdef RNSID   
        CALL RNSID_Interpolate(xGP, V70) 
        rho         = V70(60) 
        VV_cov(1:3) = V70(61:63) 
        p           = V70(64) 
        BV          = 0.0 
#else
        PRINT *, ' RNSID not available. Please compile with -DRNSID flag. '
        STOP 
#endif      
        
       !
#ifdef CCZ4GRMHD         
       V0(60:68) = (/ rho, VV_cov(1:3), p, BV(1:3), 0.0 /)
#endif
       ! 
    CASE('CCZ4GasCloud')
       !
       V0 = 0.0
       !
       V0(1)  = 1.0 
       V0(4)  = 1.0 
       V0(6)  = 1.0 
       V0(17) = 1.0 
       V0(55) = 1.0 !- 0.8*EXP(-0.5*SUM(xGP(1:nDim)**2)/0.1**2)   
       !
       V0(60:68) = (/ 2e-3 + 1e-2*EXP(-0.5*SUM(xGP(1:nDim)**2)/0.1**2), 0.0, 0.0, 0.0, 2e-3 + 1e-2*EXP(-0.5*SUM(xGP(1:nDim)**2)/0.1**2), 0.0, 0.0, 0.0, 0.0 /) 
       !
       !
       CONTINUE
       !
    END SELECT 
#endif 
    !
    CALL PDEPrim2Cons(u0,V0) 
    !
    CONTINUE
    !
END SUBROUTINE InitialField
    

  SUBROUTINE NSShock(rho,u,p,x,gamma,mu,Re0,M,rho0,u0,p0) 
      IMPLICIT NONE      
      ! Argument list declaration
      REAL :: x, gamma, mu, Re0, M, rho0, u0, p0
      REAL :: rho, u, p 
      ! Local variable declarations
      INTEGER :: i 
      REAL :: lambda,M2,xf 
      REAL :: u1,u2,f1,f2,ui,fi,ubar,pbar,mf,dudx 
      REAL, PARAMETER :: tol = 1e-12
      INTEGER, PARAMETER :: ShockType = 1 
      !
      M2     = M*M 
      lambda = ( 1 + 0.5*(gamma-1)*M2 ) / (0.5*(gamma+1)*M2 ) 
      !
      ! Use the regula falsi in order to find the root of the transcendental velocity equation 
      !
      u1 = lambda
      u2 = 1 
      !
      IF(ShockType.EQ.0) THEN
        xf = x
      ELSE
        xf = -x
      ENDIF      
      !
      CALL NSSFunction(f1,u1,xf,gamma,lambda,M2,Re0) 
      CALL NSSFunction(f2,u2,xf,gamma,lambda,M2,Re0) 
      !
      IF(f1*f2.GT.0.) THEN 
        PRINT *, ' ERROR: Sign does not change in NSSFunction. '
        PRINT *, ' Info: ', u0,u1,x,lambda,M2 
        PRINT *, '       ', f1,f2 
        STOP  
      ENDIF
      !mf = MIN( ABS(f1), ABS(f2) )
      !
      DO i = 1, 500   
          IF(ABS(f1).LT.tol) THEN
              ubar  = u1 
              pbar  = 1.
              EXIT 
          ENDIF
          IF(ABS(f2).LT.tol) THEN
              ubar  = u2 
              pbar  = 1. 
              EXIT
          ENDIF          
          ui = 0.5*(u1+u2) ! (f1*u2 - f2*u1)/(f1-f2)
          CALL NSSFunction(fi,ui,xf,gamma,lambda,M2,Re0) 
          IF(f1*fi.GT.0.) THEN
            u1 = ui
            f1 = fi 
          ELSE
            u2 = ui 
            f2 = fi
          ENDIF          
      ENDDO 
      !
      !IF( ABS(f1)>tol .AND. ABS(f2)>tol ) THEN
      !  PRINT *, ' Newton did not converge. ' 
      !  STOP 
      !ENDIF
      ! 
      ! The function is very stiff. If 100 iterations do not lead to convergence, then
      ! simply take the best values you can get... 
      !
      IF(ABS(f1).LT.ABS(f2)) THEN
          ubar  = u1           
      ELSE
          ubar  = u2 
      ENDIF          
      !
      dudx = 0.5*(gamma+1)*(ubar-1.)*(ubar-lambda)/( 4./3./Re0*gamma*ubar )
      pbar = 1 - ubar +  4./3./Re0*dudx 
      !
      rho   = rho0 / ubar    
      IF(ShockType.EQ.0) THEN
         u     = M*SQRT(gamma*p0/rho0) * ubar 
      ELSE
         u     = M*SQRT(gamma*p0/rho0) * ( 1. - ubar ) 
      ENDIF
      p     = p0 + pbar*M**2*gamma*p0 
      !
  END SUBROUTINE NSShock
  !
  SUBROUTINE NSSFunction(f,u,x,gamma,lambda,M2,Re0)
      IMPLICIT NONE
      ! Argument list declaration
      REAL :: f,u,x,Re0,M2,gamma,lambda  
      ! Local variable declaration
      REAL :: RHS, exponent 
      !
      exponent = MIN( 0.75*Re0*(M2-1)/gamma/M2 * x, 10. ) 
      RHS      = ABS( 0.5*(1-lambda) )**(1-lambda) * EXP( exponent ) 
      f        = ABS(u-1) - ABS(u-lambda)**lambda * RHS 
      !
  END SUBROUTINE      
    
