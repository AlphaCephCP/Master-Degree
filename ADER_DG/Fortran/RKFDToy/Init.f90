SUBROUTINE Init
    USE MainVariables
    IMPLICIT NONE 
    INTEGER :: i, j, k 
    REAL    :: xGP(d), u0(nVar)  
    !
    xL = (/ -1.0, -1.0, -1.0 /) 
    xR = (/ +1.0, +1.0, +1.0 /)     
    !    
#ifdef CCZ4EINSTEIN 
    xL = (/ -0.5, -0.1, -1.0 /) 
    xR = (/ +0.5, +0.1, +1.0 /)     
#endif     
    !
    tend = 0.5 
    !
    nDim = 1
    !
    IMAX = 400 !  100 !   
    JMAX = 1   !  100 ! 
    KMAX = 1   !    1 !  
    !
    NMAX = 10000000   
    !
    VMAX = (/ IMAX, JMAX, KMAX /) 
    dx = (xR-xL)/MAX(1.0,REAL(VMAX))   
    time = 0.0 
    dtio = 0.1  
    iter = 0 
    !
#ifdef FOS    
    CFL = 0.9    
    amax = 1.0 + SQRT(2.0) 
    BaseFile = 'Test-fos' 
#endif 
#ifdef SOS
    CFL = 0.16    
    BaseFile = 'Test-sos' 
#endif 
#ifdef ACOUSTIC 
    CFL = 0.9  
    amax = 1.0 
    BaseFile = 'AcousticWave' 
    dtio = 0.01 
#endif 
#ifdef CCZ4EINSTEIN 
    CFL = 0.9 ! 0.45        
    amax = 2.0     
    tend = 1000.0 
    BaseFile = 'GaugeWave-A01-FOCCZ4-fvdiss' 
    dtio = 1.0  
#endif 
    !
    tio  = dtio  
    !
    ALLOCATE( uh(nVar,IMAX,JMAX,KMAX)  ) 
    ALLOCATE( k1(nVar,IMAX,JMAX,KMAX), k2(nVar,IMAX,JMAX,KMAX) ) 
    ALLOCATE( k3(nVar,IMAX,JMAX,KMAX), k4(nVar,IMAX,JMAX,KMAX) ) 
    !
    DO k = 1, KMAX
     DO j = 1, JMAX 
      DO i = 1, IMAX
          xGP = xL + 0.5*dx + (/ i-1, j-1, k-1 /)*dx 
          CALL InitialField(u0,xGP,0.0) 
          uh(:,i,j,k) = u0 
      ENDDO
     ENDDO
    ENDDO  
    !
    CALL WriteData
    !
    CONTINUE 
    !
END SUBROUTINE Init 
    

SUBROUTINE InitialField(u0,xGP,tGP) 
    USE MainVariables
    IMPLICIT NONE 
    REAL :: u0(nVar), xGP(d), tGP 
    REAL :: ICA, HH, dxH, dxphi, phi, g_cov(3,3), g_contr(3,3), Kxx, traceK 
    REAL, PARAMETER :: Pi = ACOS(-1.0) 
    !
    u0 = 0.0 
#ifdef TOYPDE     
    IF( xGP(1) <= 0.0 ) THEN 
        u0(1) = 1.0
    ELSE
        u0(1) = 0.1
    ENDIF    
#endif 

#ifdef SOS 
    u0(1) = 1.0 + 0.5*EXP(-0.5*xGP(1)**2/0.1**2)     
    !IF( xGP(1) <= 0.0) THEN
    !    u0(1) = 1.0 
    !    u0(2) = 1.0 
    !ELSE
    !    u0(1) = 0.1 
    !    u0(2) = 0.5 
    !ENDIF    
#endif 

#ifdef FOS
    u0(1) = 1.0 + 0.1*EXP(-0.5*xGP(1)**2/0.1**2)     
    u0(3) = -0.1*xGP(1)/0.1**2*EXP(-0.5*xGP(1)**2/0.1**2) 
#endif 
    !
#ifdef ACOUSTIC 
    u0 = 0.0 
    u0(4) = 0.1*EXP(-0.5*SUM(xGP(1:nDim)**2)/0.1**2) 
    !u0(4) = SIN( 2.0*Pi*( xGP(1) - tGP)) 
#endif 

    !
#ifdef CCZ4EINSTEIN 
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
        u0(:)  = 0.0
        u0(1)  = phi**2*HH                          ! \tilde\gamma_xx
        u0(4)  = phi**2                             ! \tilde\gamma_yy
        u0(6)  = phi**2                             ! \tilde\gamma_zz
        u0(7)  = phi**2*(Kxx - 1.0/3.0*traceK*HH )  ! \tilde A_{xx}
        u0(10) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_{yy}
        u0(12) = phi**2*(0.0 - 1.0/3.0*traceK*1.0)  ! \tilde A_[zz}
        u0(17) = LOG( SQRT(HH) ) 
        u0(14) = 2.0/(3.0*HH**(5.0/3.0))*dxH        ! Gtilde 
        !
        ! Auxiliary variables
        u0(24) = 1.0/(2.0*HH)*dxH                ! A_x  
        u0(36) = HH**(-1.0/3.0)*dxH/3.0          ! D_xxx
        u0(39) = phi*dxphi                       ! D_xyy
        u0(41) = phi*dxphi                       ! D_xzz
        ! 
        u0(54) = traceK
        u0(55) = LOG( phi ) 
        u0(56) = dxphi/phi                       ! P_x
        !
#endif 
    !
END SUBROUTINE InitialField 
    
    