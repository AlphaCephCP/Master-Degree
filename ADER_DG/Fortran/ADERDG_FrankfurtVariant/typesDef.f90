MODULE typesDef
  IMPLICIT NONE 
  PUBLIC 
  !
  ! ================================== This part of the typesDef can be modified by the user.  ================================== 
  !
  INTEGER, PARAMETER :: N = 4                               ! Polynomial degree of our approximation in space and time 
  INTEGER, PARAMETER :: nDim = 2                            ! The number of space dimensions that we actually want to simulate 
  REAL, PARAMETER    :: CFL = 0.9                           ! The Courant-Friedrichs-Lewy number < 1 
#ifdef EULER
  INTEGER, PARAMETER :: nVar = 5                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 0                          ! The number of static parameters for the PDE system 
#endif 
#ifdef ELASTICITY
  INTEGER, PARAMETER :: nVar = 9                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER :: nParam = 3                          ! The number of material parameters for the PDE system 
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
  REAL               :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes) 
  REAL               :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse 
  REAL               :: Kxi(N+1,N+1), dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator 
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
  REAL               :: dx(d), dt                           ! The vector of mesh spacings in each dimension and the time step   
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
  CHARACTER(LEN=200) :: ICType  
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
  ! 
  TYPE tFace
    REAL, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector 
    REAL, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector 
    REAL, POINTER    :: paramL(:,:,:), paramR(:,:,:)        ! pointer to left and right boundary-extrapolated material parameters 
    INTEGER          :: Left, Right                         ! pointer to left and right element 
    REAL             :: nv(d)                               ! face normal vector 
  END TYPE      
  TYPE(tFace), POINTER :: Face(:) 
  !
  ! The main variables of the ADER-DG scheme 
  !  
  REAL, POINTER      :: uh(:,:,:,:,:)                       ! the coefficients of the DG polynomial        (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: duh(:,:,:,:,:)                      ! the update coefficients of the DG polynomial (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qhi(:,:,:,:,:)                      ! the time-averaged coefficients of the space-time predictor (nVar, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: Fhi(:,:,:,:,:,:)                    ! the time-averaged coefficients of the flux tensor of the space-time predictor (nVar, d, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: qbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the state vector Q in the element 
  REAL, POINTER      :: Fbnd(:,:,:,:,:)                     ! the boundary-extrapolated data for the normal flux F in the element 
  REAL, POINTER      :: parh(:,:,:,:,:)                     ! the coefficients of the DG polynomial for the material parameters (nParam, nDOF(1), nDOF(2), nDOF(3), nElem) 
  REAL, POINTER      :: parbnd(:,:,:,:,:)                   ! the boundary-extrapolated material parameters 
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
  INTEGER                     :: nPointSource
  TYPE(tPointSource), POINTER :: PointSrc(:) 
  !
  ! Important info and parameters concerning the governing PDE system 
  !
  TYPE tEquations 
      REAL           :: gamma, Pi  
  END TYPE tEquations 
  
  TYPE(tEquations)   :: EQN 
  
  !
END MODULE typesDef 
    
    
    
    