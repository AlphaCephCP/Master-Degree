MODULE typesDef
  IMPLICIT NONE 
  PUBLIC 
  !
  ! ================================== This part of the typesDef can be modified by the user.  ================================== 
  !
  INTEGER, PARAMETER :: N = 3                               ! Polynomial degree of our approximation in space and time 
  INTEGER, PARAMETER :: nDim = 3                            ! The number of space dimensions that we actually want to simulate 
  REAL, PARAMETER    :: CFL = 0.9                           ! The Courant-Friedrichs-Lewy number < 1 
#ifdef EULER
  INTEGER, PARAMETER :: nVar = 5                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of static parameters for the PDE system 
#endif 
#ifdef ELASTICITY
  INTEGER, PARAMETER :: nVar = 14                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef Z4EINSTEIN
  INTEGER, PARAMETER :: nVar = 54                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef Z4GRMHD
  INTEGER, PARAMETER :: nVar = 63                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef CCZ4EINSTEIN 
  INTEGER, PARAMETER :: nVar = 59                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef CCZ4GRMHD
  INTEGER, PARAMETER :: nVar = 68                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef BSSZ4EINSTEIN 
  INTEGER, PARAMETER :: nVar = 62                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef BSSZ4GRMHD 
  INTEGER, PARAMETER :: nVar = 72                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef ACOUSTIC 
  INTEGER, PARAMETER :: nVar = 4                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef SRMHD 
  INTEGER, PARAMETER :: nVar = 9                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef GRMHD 
  INTEGER, PARAMETER :: nVar = 19                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
#ifdef GPR3D  
  INTEGER, PARAMETER :: nVar = 17                           ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of material parameters for the PDE system 
#endif   
  !
#ifdef GRMHD
  INTEGER, PARAMETER :: nEvolve = 9  
#else   
  INTEGER, PARAMETER :: nEvolve = nVar 
#endif   
  !  
  ! ==================================           Do NOT change the stuff below                 ==================================
  !
  ! The following variables contain important information about the numerical method. Do NOT change.  
  !
  INTEGER, PARAMETER :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !! 
  REAL, PARAMETER    :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes according to von Neumann stability analysis and experience     
  INTEGER            :: nDOF(0:3)                           ! The number of degrees of freedom in space and time 
  REAL               :: xiGPN(N+1), wGPN(N+1)               ! The Gauss-Legendre quadrature nodes and weights 
  REAL               :: xiGPMm1(N), wGPMm1(N)               ! The points for basis of degree N-1 
  REAL               :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes) 
  REAL               :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse 
  REAL               :: Kxi(N+1,N+1), dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator 
  REAL               :: Kxixi(N+1,N+1), Lxixi(N+1)          ! Element stiffness matrix for an internal fourth order artificial dissipation operator and its eigenvalues 
  REAL               :: Rxixi(N+1,N+1), iRxixi(N+1,N+1)     ! Right and left eigenvectors of the fourth order dissipation operator 
  REAL               :: FLm(N+1,N+1), FLp(N+1,N+1)          ! Left flux matrices 
  REAL               :: FRm(N+1,N+1), FRp(N+1,N+1)          ! Right flux matrices 
  REAL               :: FLcoeff(N+1), FRcoeff(N+1)          ! extrapolation coefficients to the left and right boundary 
  REAL               :: F0(N+1), F1(N+1,N+1)                ! Time flux matrices 
  REAL               :: K1(N+1,N+1), iK1(N+1,N+1)           ! F1 - Ktau 
  INTEGER            :: dn(d)                               ! number of direct neighbors in each dimension 
  REAL               :: MT(N+1,N+1), iMT(N+1,N+1)           ! Time mass matrix (for point source terms) 
  INTEGER            :: Face2Neigh(3,6)                     ! Mapping from face to neighbor index 
  ! Periodic boundary conditions 
  LOGICAL            :: Periodic(d)                         ! periodic BC in x, y, z direction 
  ! Stuff related to the problem setup, mesh and output 
  INTEGER            :: IMAX, JMAX, KMAX, NMAX              ! The number of cells in each space dimension & max. number of time steps 
  INTEGER            :: nElem, nFace, nNode                 ! The number of elements, faces and nodes 
  INTEGER            :: timestep                            ! the number of the current time step 
  REAL               :: xL(d), xR(d)                        ! computational domain 
  REAL               :: dx(d), dt, maxdt                    ! The vector of mesh spacings in each dimension and the time step   
  REAL               :: time, tend                          ! current time  and final time 
  REAL               :: tio, dtio                           ! output time and output time interval 
  REAL, POINTER      :: x(:,:)                              ! the node coordinates (nDim, nNode) 
  INTEGER, PARAMETER :: nVtx = 2**nDim, nFac = 2*nDim       ! number of vertices and faces per element 
  INTEGER, POINTER   :: tri(:,:)                            ! connectivity from the element to the nodes 
  INTEGER, POINTER   :: Element2Face(:,:)                   ! connectivity from each element to its faces 
  CHARACTER(LEN=200) :: BaseFile                            ! Basic filename to write the results 
  REAL               :: SubOutputMatrix((N+1)**d,(N+1)**d)  ! Matrix needed for the plotting of the results on a fine subgrid 
  INTEGER            :: subtri(2**d,N**d)                   ! subcell connectivity (for fine output) 
  REAL               :: allsubxi(d,(N+1)**d)                ! subnodes (for fine output) 
  ! Some diagnostics or data analysis                       ! 
  REAL               :: tCPU1, tCPU2                        ! CPU times 
  INTEGER(8)         :: TEU                                 ! total element updates 
  INTEGER            :: AnalyseType   
  CHARACTER(LEN=200) :: ICType, ICType2  
  ! Data needed for the subcell limiter 
  INTEGER, PARAMETER :: nSubLim = 2*N+1                     ! number of subcells 
  INTEGER            :: nSubLimV(d)                         ! vector for number of subcells in each dimension 
  REAL, POINTER      :: uh2lim(:,:), lim2uh(:,:)            ! mapping from DG polynomials to the subcells and back 
  REAL, POINTER      :: uh2lob(:,:)                         ! mapping from the DG polynomials to the Gauss-Lobatto points 
  REAL               :: xiLob(N+1), wLob(N+1)               ! Gauss lobatto points and weights 
  INTEGER, POINTER   :: neighbor(:,:,:,:)                   ! set of Voronoi neighbors of a cell 
  REAL, POINTER      :: olduh(:,:,:,:,:)                    ! for the a posteriori limiter, we need the possibility to go back 
  INTEGER, POINTER   :: recompute(:)                        ! map containing the troubled zones that need to be recomputed 
  INTEGER            :: subtrilim(2**d,(nSubLim)**d)        ! subcell connectivity (for fine output) 
  REAL               :: subxilim(d,(nSubLim+1)**d)          ! subnodes (for fine output) 
  REAL               :: xilimbary(nSubLim)                  ! barycenters of the subcells 
  ! Definition of the unit reference element 
  REAL               :: ReferenceElement(d,2**d) 
  ! Kreiss-Oliger dissipation coefficient 
  REAL               :: epsilon_ko = 1e-2                   ! use some fourth order artificial dissipation if necessary 
  TYPE tFace
    REAL, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector 
    REAL, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector 
    REAL, POINTER    :: paramL(:,:,:), paramR(:,:,:)        ! pointer to left and right boundary-extrapolated material parameters 
    INTEGER          :: Left, Right                         ! pointer to left and right element 
    REAL             :: nv(d), x0(d)                        ! face normal vector and left lower corner coordinate  
  END TYPE      
  TYPE(tFace), POINTER :: Face(:) 
  !
  ! The main variables of the ADER-DG scheme 
  !  
  REAL, POINTER      :: uh(:,:,:,:,:)                       ! the coefficients of the DG polynomial        (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: duh(:,:,:,:,:)                      ! the update coefficients of the DG polynomial (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qhi(:,:,:,:,:)                      ! the time-averaged coefficients of the space-time predictor (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: Fhi(:,:,:,:,:,:)                    ! the time-averaged coefficients of the flux tensor of the space-time predictor (nVar, d, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: Shi(:,:,:,:,:)                      ! the time-averaged coefficients of the source vector of the space-time predictor (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the state vector Q in the element 
  REAL, POINTER      :: Fbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the normal flux F in the element 
  REAL, POINTER      :: parh(:,:,:,:,:)                     ! the coefficients of the DG polynomial for the material parameters (nParam, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: parbnd(:,:,:,:,:)                   ! the boundary-extrapolated material parameters 
  REAL, POINTER      :: ADM(:,:,:,:,:)                      ! ADM constraints in each point 
  REAL, POINTER      :: RKDGk1(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk2(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk3(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk4(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk5(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk6(:,:,:,:,:) 
  REAL, POINTER      :: RKDGk7(:,:,:,:,:) 
  !
  TYPE tLimiter 
      INTEGER        :: status, oldstatus 
      REAL, POINTER  :: Lh(:,:,:,:) 
      REAL, POINTER  :: NewLh(:,:,:,:) 
      REAL           :: lmin(nVar), lmax(nVar)       
  END TYPE tLimiter  
  TYPE(tLimiter), POINTER :: Limiter(:) 
  !
  TYPE tPointSource
      INTEGER        :: iElem, waveform 
      REAL           :: xP(d), xiP(d)  
      REAL, POINTER  :: phi(:,:,:) 
      REAL, POINTER  :: sigma(:,:), SrcInt(:) 
  END TYPE tPointSource
  INTEGER                     :: nPointSource = 0 
  TYPE(tPointSource), POINTER :: PointSrc(:) 
  !
  ! Important info and parameters concerning the governing PDE system 
  !
  TYPE tEquations 
      REAL    :: gamma, Pi, c0   
      REAL    :: CCZ4k1, CCZ4k2, CCZ4k3, CCZ4eta, CCZ4itau, CCZ4f, CCZ4g, CCZ4xi, CCZ4e, CCZ4c, CCZ4mu, CCZ4ds, CCZ4sk, CCZ4bs  
      REAL    :: cs, alpha, cv, rho0, p0, tau1, tau2, mu, kappa  
      INTEGER :: CCZ4LapseType, EinsteinAutoAux = 0   
      REAL    :: DivCleaning_a = 1.0 
  END TYPE tEquations 
  ! 
  TYPE(tEquations)   :: EQN 
  !  
#ifdef GRMHD  
  REAL, PARAMETER    :: ExcisionRadius = 1.0, p_floor = 1.0e-11, rho_floor = 1.0e-10
  REAL, PARAMETER    :: aom = 0.0, Mbh = 1.0
  REAL, PARAMETER    :: P_eps = 1e-4  
  REAL, PARAMETER    :: b0_torus_par = 0.0, al_torus_par = 3.8   ! magnetic field parameter and angular momentum parameter for the torus  
#else
  REAL    :: ExcisionRadius 
#endif
  !
  ! Stuff required for MPI 
  !
  INTEGER               :: myrank, nCPU, mpiErr, MPI_AUTO_REAL
  INTEGER               :: nCPUx(d) = (/ 1, 1, 1 /) 
  INTEGER               :: nMPIElem, nInnerFaces, nOuterFaces, nInnerElements, nOuterElements 
  INTEGER, ALLOCATABLE  :: InnerElements(:), OuterElements(:) 
  INTEGER, ALLOCATABLE  :: CPUDistribution(:), LocalCPUDistribution(:)  
  INTEGER, ALLOCATABLE  :: Local2GlobalElem(:), Global2LocalElem(:)  
  INTEGER, ALLOCATABLE  :: nMPISendElem(:),  MPISendElem(:,:),  nMPIRecvElem(:),  MPIRecvElem(:,:)
  INTEGER, POINTER      :: idxn(:,:,:), idxe(:,:,:)    
  REAL                  :: mpi_real_test
  TYPE tCommList 
      INTEGER          :: nElements    
      INTEGER, POINTER :: Elements(:), FaceIndex(:)  
  END TYPE tCommList 
  TYPE tRealMessage                                                             ! Defines a vector message of type real
      REAL, POINTER      :: Content(:)
  END TYPE tRealMessage
  TYPE tIntegerMessage                                                          ! Defines a vector message of type integer
      INTEGER, POINTER   :: Content(:)
  END TYPE tIntegerMessage
  TYPE(tCommList), POINTER            :: SendElem(:), RecvElem(:)
  INTEGER                             :: nSendCPU, nRecvCPU 
  INTEGER, ALLOCATABLE                :: SendCPU(:), RecvCPU(:), invSendCPU(:), invRecvCPU(:)   
  TYPE(tRealMessage), ALLOCATABLE     :: send_message(:)
  TYPE(tRealMessage), ALLOCATABLE     :: recv_message(:)
  INTEGER, ALLOCATABLE                :: send_request(:)
  INTEGER, ALLOCATABLE                :: recv_request(:)
  INTEGER, ALLOCATABLE                :: send_status_list(:,:)
  INTEGER, ALLOCATABLE                :: recv_status_list(:,:)
  !   
  LOGICAL, PARAMeTER   :: EvalSphereIntegral =.FALSE.
  INTEGER, PARAMETER  :: nObs = 5   ! number of observables integrated along the spherical surface
  REAL, PARAMETER    :: SphereRadius = 2.0
  INTEGER, PARAMETER :: SphPick_nPhi = 100, SphPick_nTheta = 200
  INTEGER, PARAMETER :: nSpherePickPoint = SphPick_nPhi*SphPick_nTheta*(N+1)**2
  !
  TYPE tPickPoint
     REAL                                   :: x(d), picktime                   ! position of pickpoint 
     REAL                                   :: xi(d)                            ! relative position within the element 0 <= xi(i) <= 1
     INTEGER                                :: iElem0, iElem                    ! level0 element and current element containing the pickpoint (changes dynamically) 
     INTEGER                                :: currlevel                        ! current level of the element containing the pickpoint  
     CHARACTER(LEN=200)                     :: FileName                         ! name of the pickpoint file
     REAL                                   :: r(d)                             ! position of pickpoint in polar coordinates
     INTEGER                                :: iTheta,iPhi
  END TYPE tPickPoint
  TYPE(tPickPoint), POINTER                 :: SpherePickPoint(:)                     ! array containing all pickpoint of one CPU (may also change dynamically in the future) 
  !  
END MODULE typesDef 
    
    
    
    