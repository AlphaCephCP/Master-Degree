MODULE MainVariables 
    !
    IMPLICIT NONE 
    !
    PUBLIC
    !    
    INTEGER, PARAMETER :: d = 3 
    !
#ifdef TOYPDE
    INTEGER, PARAMETER :: EqnType = 1 
    INTEGER, PARAMETER :: nVar = 1 
#endif 
    !
#ifdef SOS
    INTEGER, PARAMETER :: EqnType = 1 
    INTEGER, PARAMETER :: nVar = 2 
#endif 

#ifdef FOS
    INTEGER, PARAMETER :: EqnType = 2 
    INTEGER, PARAMETER :: nVar = 3 
#endif     
    !
#ifdef ACOUSTIC 
    INTEGER, PARAMETER :: EqnType = 2 
    INTEGER, PARAMETER :: nVar = 4   
#endif     
    !
#ifdef CCZ4EINSTEIN 
    INTEGER, PARAMETER :: EqnType = 2 
    INTEGER, PARAMETER :: nVar = 59  
#endif     
    ! 
    INTEGER :: IMAX, JMAX, KMAX, VMAX(d), NMAX, nDim, iter   
    REAL    :: xL(d), xR(d), dx(d), dt, CFL, time, tend, tio, dtio, amax  
    !
    REAL, POINTER :: uh(:,:,:,:), k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:) 
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
    CHARACTER(LEN=200) :: BaseFile 
    !    
END MODULE MainVariables     