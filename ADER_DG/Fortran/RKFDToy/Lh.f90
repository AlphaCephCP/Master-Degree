SUBROUTINE Lh(dudt,u)
    USE MainVariables
    IMPLICIT NONE 
    REAL    :: dudt(nVar,IMAX,JMAX,KMAX), u(nVar,IMAX,JMAX,KMAX) 
    INTEGER :: i, j, k
    INTEGER :: ip1, im1, im2, ip2, im3, ip3    
    INTEGER :: jp1, jm1, jm2, jp2, jm3, jp3    
    INTEGER :: kp1, km1, km2, kp2, km3, kp3    
    REAL    :: Q(nVar), ux(nVar), uy(nVar), uz(nVar)
    REAL    :: uxx(nVar), uxy(nVar), uxz(nVar) 
    REAL    :: uyx(nVar), uyy(nVar), uyz(nVar) 
    REAL    :: uzx(nVar), uzy(nVar), uzz(nVar), uxxxx(nVar)  
    REAL    :: gradQ(nVar,d), Hessian(nVar,d,d), dQdt(nVar)   
    REAL    :: IR(nVar,nVar), IL(nVar,nVar), smax  
    INTEGER, PARAMETER :: FD = 6 
    !
    gradQ   = 0.0 
    Hessian = 0.0 
    !
    DO k = 1, KMAX
     km1 = k-1
     kp1 = k+1
     km2 = k-2
     kp2 = k+2
     km3 = k-3
     kp3 = k+3
     IF(km1<=0) km1 = KMAX+km1 
     IF(km2<=0) km2 = KMAX+km2 
     IF(km3<=0) km3 = KMAX+km3 
     IF(kp1>=KMAX+1) kp1 = kp1-KMAX 
     IF(kp2>=KMAX+1) kp2 = kp2-KMAX 
     IF(kp3>=KMAX+1) kp3 = kp3-KMAX 
     ! 
     DO j = 1, JMAX 
      jm1 = j-1
      jp1 = j+1
      jm2 = j-2
      jp2 = j+2
      jm3 = j-3
      jp3 = j+3
      IF(jm1<=0) jm1 = JMAX+jm1 
      IF(jm2<=0) jm2 = JMAX+jm2 
      IF(jm3<=0) jm3 = JMAX+jm3 
      IF(jp1>=IMAX+1) jp1 = jp1-JMAX 
      IF(jp2>=IMAX+1) jp2 = jp2-JMAX 
      IF(jp3>=IMAX+1) jp3 = jp3-JMAX 
      ! 
      DO i = 1, IMAX
        im1 = i-1
        ip1 = i+1
        im2 = i-2
        ip2 = i+2
        im3 = i-3
        ip3 = i+3
        IF(im1<=0) im1 = IMAX+im1 
        IF(im2<=0) im2 = IMAX+im2 
        IF(im3<=0) im3 = IMAX+im3 
        IF(ip1>=IMAX+1) ip1 = ip1-IMAX 
        IF(ip2>=IMAX+1) ip2 = ip2-IMAX 
        IF(ip3>=IMAX+1) ip3 = ip3-IMAX 

        ! define the finite differences here 
        Q = u(:,i,j,k)
#ifdef SOS        
        ux  = 1.0/dx(1)*( u(:,i,j,k) - u(:,im1,j,k) ) 
        uxx = 1./dx(1)**2*( u(:,ip1,j,k) - 2*u(:,i,j,k) + u(:,im1,j,k) ) 
        gradQ(:,1) = ux
        Hessian(:,1,1) = uxx 
        ! now call the PDE operator 
        CALL PDEOperator(dQdt,Q,gradQ,Hessian)  
#endif    
#ifdef FOS
        IR = 0.0
        IL = 0.0 
        !
        IR(1,1) = 1      
        IR(1,2) = 0
        IR(1,3) = 0
        IR(2,1) = 0
        IR(2,2) = 1.D0/3.D0
        IR(2,3) = 1.D0/3.D0
        IR(3,1) = 0
        IR(3,2) = 2.D0/3.D0
        IR(3,3) = 2.D0/3.D0
        !
        IL(1,1) = 0      
        IL(1,2) = 0
        IL(1,3) = 0
        IL(2,1) = 0
        IL(2,2) = 2.D0/3.D0
        IL(2,3) = -1.D0/3.D0
        IL(3,1) = 0
        IL(3,2) = -2.D0/3.D0
        IL(3,3) = 1.D0/3.D0        
        !
        ux  = 1.0/dx(1)*( MATMUL( IR, u(:,i,j,k) - u(:,im1,j,k)) + MATMUL( IL, u(:,ip1,j,k) - u(:,i,j,k))  ) 
        !
        uxx = 1./dx(1)**2*( u(:,ip1,j,k) - 2*u(:,i,j,k) + u(:,im1,j,k) ) 
        gradQ(:,1) = ux
        Hessian(:,1,1) = uxx 
        ! now call the PDE operator 
        CALL PDEOperator(dQdt,Q,gradQ,Hessian)  
        !        
#endif
        !
#if defined(CCZ4EINSTEIN) || defined(ACOUSTIC) 
        ! test a simple second order accurate central finite difference 
        !ux = ( u(:,ip1,j,k) - u(:,im1,j,k) ) / (2 * dx(1) ) 
        ! test a fourth order central finite difference 
        !ux = ( 8.0*u(:,ip1,j,k) - 8.0*u(:,im1,j,k) + u(:,im2,j,k) - u(:,ip2,j,k) )/( 12.0*dx(1) ) 
        ! test a sixth order central finite difference 
        ux = ( 9.0*u(:,im2,j,k) - 45.0*u(:,im1,j,k) + 45.0*u(:,ip1,j,k) - 9.0*u(:,ip2,j,k) - u(:,im3,j,k) + u(:,ip3,j,k)  ) / ( 60.0*dx(1) ) 
        IF(nDim>=2) THEN
            uy = ( 9.0*u(:,i,jm2,k) - 45.0*u(:,i,jm1,k) + 45.0*u(:,i,jp1,k) - 9.0*u(:,i,jp2,k) - u(:,i,jm3,k) + u(:,i,jp3,k)  ) / ( 60.0*dx(2) ) 
        ENDIF
        IF(nDim>=3) THEN
            uz = ( 9.0*u(:,i,j,km2) - 45.0*u(:,i,j,km1) + 45.0*u(:,i,j,kp1) - 9.0*u(:,i,j,kp2) - u(:,i,j,km3) + u(:,i,j,kp3)  ) / ( 60.0*dx(3) ) 
        ENDIF 
        !uxx = 1./dx(1)**2*( u(:,ip1,j,k) - 2*u(:,i,j,k) + u(:,im1,j,k) )
        gradQ(:,1) = ux 
        gradQ(:,2) = uy 
        gradQ(:,3) = uz  
        CALL PDEFusedSrcNCP(dQdt,Q,gradQ) 
        ! fourth order derivative (Kreiss-Oliger dissipation) 
        !uxxxx = ( u(:,im2,j,k) - 4.0*u(:,im1,j,k) + 6.0*u(:,i,j,k) - 4.0*u(:,ip1,j,k) +  u(:,ip2,j,k) ) 
        !dQdt = dQdt - dt/dx(1)*uxxxx    ! it does not seem to work without the factor dt 
        ! numerical dissipation that we would have with a fourth order finite volume scheme and Rusanov flux (signal speed smax) 
        smax = EQN%CCZ4e 
        dQdt = dQdt - 1.0/dx(1)* 3.0/256.0*smax* ( -15.0*u(:,ip1,j,k)-u(:,ip3,j,k)-15.0*u(:,im1,j,k)+6.0*u(:,ip2,j,k)+20.0*u(:,i,j,k)+6.0*u(:,im2,j,k)-u(:,im3,j,k) ) 
        IF(nDim>=2) THEN
            dQdt = dQdt - 1.0/dx(2)* 3.0/256.0*smax* ( -15.0*u(:,i,jp1,k)-u(:,i,jp3,k)-15.0*u(:,i,jm1,k)+6.0*u(:,i,jp2,k)+20.0*u(:,i,j,k)+6.0*u(:,i,jm2,k)-u(:,i,jm3,k) ) 
        ENDIF 
        IF(nDim>=3) THEN
            dQdt = dQdt - 1.0/dx(3)* 3.0/256.0*smax* ( -15.0*u(:,i,j,kp1)-u(:,i,j,kp3)-15.0*u(:,i,j,km1)+6.0*u(:,i,j,kp2)+20.0*u(:,i,j,k)+6.0*u(:,i,j,km2)-u(:,i,j,km3) ) 
        ENDIF 
        ! 
        CONTINUE 
#endif 
        ! 
        dudt(:,i,j,k) = dQdt(:) 
      ENDDO 
     ENDDO
    ENDDO    
    !
END SUBROUTINE Lh 
    