! Set nDim to the desired number of space dimensions in Init.f90
! To solve the FO-CCZ4 system, define the preprocessor -DCCZ4EINSTEIN 
! To solve the acoustic wave equations, define the preprocessor -DACOUSTIC  
! To suppress IO for speed tests, define -DNOIO 
! 
! Have fun :-) 
! 
PROGRAM EinsteinFD
    USE MainVariables 
    IMPLICIT NONE 
    INTEGER :: i,j,k, iVar 
    REAL    :: t1, t2, t3, TEU, L1error(nVar), L2error(nVar), Linferror(nVar), Qe(nVar), xGP(d), locerror(nVar)      
    INTEGER, PARAMETER  :: RK = 4 
    CHARACTER(LEN=10)  :: VarName(nVar)    
    !
    PRINT *, ' --------------------------------------------- ' 
    PRINT *, '      Welcome to a simple FD toy code          ' 
    PRINT *, '         Written by Michael Dumbser            ' 
    PRINT *, '        University of Trento, Italy            ' 
    PRINT *, ' --------------------------------------------- ' 
    !
    CALL Init 
    !
    CALL CPU_TIME(t1) 
    !
    DO iter = 1, NMAX
        !
        dt = CFL/SUM( amax/dx(1:nDim) ) 
        ! 
        IF(time>=tend) THEN
            EXIT 
        ENDIF        
        IF(time+dt>tend) THEN
            dt = tend - time
        ENDIF
        IF(time+dt>tio) THEN 
            dt = tio - time 
        ENDIF 
        ! Explicit Euler scheme
        SELECT CASE(RK) 
        CASE(1) 
            ! Explicit Euler 
            CALL Lh(k1,uh)                    
            uh = uh + dt*k1 
        CASE(2) 
            ! RK2 scheme
            CALL Lh(k1,uh)
            CALL Lh(k2,uh+dt*k1) 
            uh = uh + dt/2.0*(k1+k2)          
        CASE(3) 
            ! Kutta's third order scheme 
            CALL Lh(k1,uh) 
            CALL Lh(k2,uh+0.5*dt*k1) 
            CALL Lh(k3,uh-1.0*dt*k1+2.0*dt*k2) 
            uh = uh + dt/6.0*( k1 + 4.0*k2 + 1.0*k3  ) 
        CASE(4)             
            ! Classical RK4 scheme 
            CALL Lh(k1,uh) 
            CALL Lh(k2,uh+0.5*dt*k1) 
            CALL Lh(k3,uh+0.5*dt*k2) 
            CALL Lh(k4,uh+1.0*dt*k3) 
            uh = uh + dt/6.0*( k1 + 2.0*k2 + 2.0*k3 + k4 ) 
        END SELECT 
        !
#ifdef CCZ4EINSTEIN         
        DO k = 1, KMAX
         DO j = 1, JMAX 
          DO i = 1, IMAX
            CALL EnforceCCZ4Constraints(uh(:,i,j,k)) 
          ENDDO
         ENDDO
        ENDDO 
#endif         
        !
        time = time + dt 
        !
        IF(MOD(iter,100)==0) THEN
            CALL CPU_TIME(t3) 
            !CALL WriteData
            WRITE(*,'(a,i8,a,f14.6,a,e13.6,a,e10.3,a,f10.3)') ' n = ', iter, '   t = ', time, '   WCT [s] = ', t3-t1, '  TDU [mus] = ',  (t3-t1)/(iter*IMAX*JMAX*KMAX)*1e6   
        ENDIF
        IF( time >= tio ) THEN
            CALL WriteData 
            PRINT *, ' Current time = ', time 
            tio = tio + dtio 
        ENDIF 
        !
#ifdef CCZ4EINSTEIN
#ifndef NOIO 
        CALL RunTimeAnalyseData
#endif         
#endif         
        !
    ENDDO    
    !
    CALL CPU_TIME(t2) 
    !
    CALL WriteData 
    !
    PRINT *, ' --------------------------------------------- ' 
    PRINT *, '   Computing error norms... ' 
    L1error     = 0.0
    L2error     = 0.0
    Linferror   = 0.0 
    !
    DO k = 1, KMAX
     DO j = 1, JMAX 
      DO i = 1, IMAX
        xGP = xL + 0.5*dx + (/ i-1, j-1, k-1 /)*dx 
        CALL InitialField(Qe,xGP,time) 
        locerror = ABS( Qe - uh(:,i,j,k) ) 
        L1error = L1error + PRODUCT(dx(1:nDim))*locError  
        L2error = L2error + PRODUCT(dx(1:nDim))*locError**2   
        DO iVar = 1, nVar 
            Linferror(iVar) = MAX( Linferror(iVar), locerror(iVar) ) 
        ENDDO                 
      ENDDO
     ENDDO
    ENDDO 
    L2error = SQRT(L2error) 
    CALL PDEVarName(VarName)     
    DO iVar = 1, nVar
        PRINT *,  VarName(iVar) 
        PRINT *, ' L1error   = ', L1error(iVar) 
        PRINT *, ' L2error   = ', L2error(iVar) 
        PRINT *, ' Linferror = ', Linferror(iVar) 
    ENDDO    
    !
    PRINT *, ' --------------------------------------------- ' 
    PRINT *, '   Program finished. ' 
    PRINT *, ' --------------------------------------------- ' 
    PRINT *, '   Total CPU time [s] = ', t2 - t1 
    TEU = REAL(iter)*REAL(PRODUCT(VMAX(1:nDim)))   
    PRINT *, '   Time / EU [mu s]   = ', (t2 - t1)/TEU 
    PRINT *, ' --------------------------------------------- ' 
    !
END PROGRAM EinsteinFD
    
    