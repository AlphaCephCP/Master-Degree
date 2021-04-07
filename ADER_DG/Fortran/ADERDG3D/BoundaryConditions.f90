SUBROUTINE BoundaryConditions 
    USE typesDef
    ! Local variables 
    INTEGER :: iFace 
    INTEGER :: i,j,k
    REAL    :: Qbc(nVar),Fbc(nVar,d),Vbc(nVar),TM(9,9), RM(9,9), sv(3), tv(3)   
    REAL    :: nv(d), a1, a2, a3, QRrot(9), QLrot(9)  
    REAL    :: lambda, mu, rho, xGP(d), QGP(nVar), FGP(nVar,d), parGP(nParam)   
    INTEGER :: ibc(1)  
    !
    CONTINUE
    !
#ifdef LINEARELASTICITY
    !
    Qbc = 0. 
    DO iFace = 1, nFace
        !
        ! Here, we take more care of the boundary conditions... :-)  
        IF(Face(iFace)%Left.EQ.0) THEN
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2)
                    
                    nv = Face(iFace)%nv    ! we need to flip the sign here, since the left element is missing, hence the outward unit normal points to the left ! 
                    IF(ABS(nv(1)-1.0) < 1e-14) THEN
                        nv = (/-1., 0., 0. /) 
                        sv = (/ 0., 0., 1. /)
                        tv = (/ 0., 1., 0. /) 
                    ENDIF    
                    IF(ABS(nv(2)-1.0) < 1e-14) THEN
                        nv = (/ 0.,-1., 0. /) 
                        sv = (/ 1., 0., 0. /)
                        tv = (/ 0., 0., 1.  /) 
                    ENDIF    
                    IF(ABS(nv(3)-1.0) < 1e-14) THEN
                        nv = (/ 0., 0.,-1. /) 
                        sv = (/ 0., 1., 0. /)
                        tv = (/ 1., 0., 0. /) 
                    ENDIF  
                    !
                    ! Rotation matrix to rotate into local coordinate system aliged with normal vector n  
                    !
                    TM(1,:) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
                    TM(2,:) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
                    TM(3,:) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
                    TM(4,:) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
                    TM(5,:) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
                    TM(6,:) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
                    TM(7,:) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
                    TM(8,:) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
                    TM(9,:) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)     
                    
                    ! The material parameters in this point 
                    lambda = Face(iFace)%paramR(1,j,k) 
                    mu     = Face(iFace)%paramR(2,j,k) 
                    rho    = Face(iFace)%paramR(3,j,k)   
                    cp = SQRT((lambda+2*mu)/rho) 
                    cs = SQRT(mu/rho) 
                    ! Right eigenvector matrix (for the local rotated coordinate system) 
                    RM(1,:) = (/ lambda + 2*mu, 0., 0., 0., 0., 0., 0., 0., lambda + 2*mu /) 
                    RM(2,:) = (/ lambda, 0., 0., 0., 1., 0., 0., 0., lambda /) 
                    RM(3,:) = (/ lambda, 0., 0., 0., 0., 1., 0., 0., lambda /) 
                    RM(4,:) = (/ 0.,     mu, 0., 0., 1., 0., 0., mu, 0.     /) 
                    RM(5,:) = (/ 0.,     0., 0., 1., 0., 0., 0., 0., 0.     /)  
                    RM(6,:) = (/ 0.,     0., mu, 0., 0., 0., mu, 0., 0.     /)  
                    RM(7,:) = (/ cp,     0., 0., 0., 0., 0., 0., 0., -cp    /)  
                    RM(8,:) = (/ 0.,     cs, 0., 0., 0., 0., 0.,-cs, 0.     /)  
                    RM(9,:) = (/ 0.,     0., cs, 0., 0., 0.,-cs, 0., 0.     /)                      
                    ! Rotate the inner state in the local coordinate system 
                    QLrot = MATMUL( TRANSPOSE(TM), Face(iFace)%qR(:,j,k))  
                    ! The coefficients below solve the inverse Riemann problem so that the normal stress on the boundary is zero (for the Godunov flux based on the exact Riemann solver) 
                    a1 = 0.5*(-QLrot(1)+rho*cp*QLrot(7))/rho/cp**2 
                    a2 = 0.5*(-QLrot(4)+rho*cs*QLrot(8))/rho/cs**2  
                    a3 = 0.5*(-QLrot(6)+rho*cs*QLrot(9))/rho/cs**2  
                    ! Compute the outer state, and rotate back to the physical x-y-z coordinate system 
                    QRrot = a1*RM(:,1) + a2*RM(:,2) + a3*RM(:,3)                     
                    ! For the Godunov flux, which is based on the EXACT Riemann solver, we can indeed use the simple and sloppy mirror BC, since the Riemann solver will select the correct incoming and outgoing characteristics by itself !!  
                    ! This can be easily checked using the theory of characteristics 
#ifdef GODUNOV
                    QRrot = QLrot 
                    QRrot((/1,4,6/)) = -QRrot((/1,4,6/)) 
#endif 
                    ! Put free surface boundaries everywhere to check our new stuff 
                    IF(ibc(1)==1) then
                        continue
                    endif
                    !IF( ibc(1)==2 ) THEN
                        Face(iFace)%qL(:,j,k) = MATMUL( TM, QRrot ) 
                    !ELSE 
                    !    Face(iFace)%qL(:,j,k) = qbc 
                    !ENDIF 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    !     
                    nv = +Face(iFace)%nv    ! we do not need to flip the sign here, since the right element is missing, hence the outward unit normal points to the right ! 
                    IF(ABS(nv(1)-1.0) < 1e-14) THEN
                        sv = (/ 0., 1., 0. /)
                        tv = (/ 0., 0., 1. /) 
                    ENDIF    
                    IF(ABS(nv(2)-1.0) < 1e-14) THEN
                        sv = (/ 0., 0., 1. /)
                        tv = (/ 1., 0., 0.  /) 
                    ENDIF    
                    IF(ABS(nv(3)-1.0) < 1e-14) THEN
                        sv = (/ 1., 0., 0. /)
                        tv = (/ 0., 1., 0.  /) 
                    ENDIF  
                    !
                    ! Rotation matrix 
                    !
                    TM(1,:) = (/ nv(1)**2   , sv(1)**2   , tv(1)**2   , 2*nv(1)*sv(1)          , 2*sv(1)*tv(1)          , 2*nv(1)*tv(1)          , 0., 0., 0. /) 
                    TM(2,:) = (/ nv(2)**2   , sv(2)**2   , tv(2)**2   , 2*nv(2)*sv(2)          , 2*sv(2)*tv(2)          , 2*nv(2)*tv(2)          , 0., 0., 0. /) 
                    TM(3,:) = (/ nv(3)**2   , sv(3)**2   , tv(3)**2   , 2*nv(3)*sv(3)          , 2*sv(3)*tv(3)          , 2*nv(3)*tv(3)          , 0., 0., 0. /) 
                    TM(4,:) = (/ nv(2)*nv(1), sv(2)*sv(1), tv(2)*tv(1), nv(2)*sv(1)+nv(1)*sv(2), sv(2)*tv(1)+sv(1)*tv(2), nv(2)*tv(1)+nv(1)*tv(2), 0., 0., 0. /) 
                    TM(5,:) = (/ nv(3)*nv(2), sv(3)*sv(2), tv(3)*tv(2), nv(3)*sv(2)+nv(2)*sv(3), sv(3)*tv(2)+sv(2)*tv(3), nv(3)*tv(2)+nv(2)*tv(3), 0., 0., 0. /) 
                    TM(6,:) = (/ nv(3)*nv(1), sv(3)*sv(1), tv(3)*tv(1), nv(3)*sv(1)+nv(1)*sv(3), sv(3)*tv(1)+sv(1)*tv(3), nv(3)*tv(1)+nv(1)*tv(3), 0., 0., 0. /) 
                    TM(7,:) = (/ 0.,0., 0., 0., 0., 0., nv(1), sv(1), tv(1) /)     
                    TM(8,:) = (/ 0.,0., 0., 0., 0., 0., nv(2), sv(2), tv(2) /)     
                    TM(9,:) = (/ 0.,0., 0., 0., 0., 0., nv(3), sv(3), tv(3) /)     
                    ! The material parameters in this point 
                    lambda = Face(iFace)%paramL(1,j,k) 
                    mu     = Face(iFace)%paramL(2,j,k) 
                    rho    = Face(iFace)%paramL(3,j,k)   
                    cp = SQRT((lambda+2*mu)/rho) 
                    cs = SQRT(mu/rho) 
                    ! Right eigenvector matrix 
                    RM(1,:) = (/ lambda + 2*mu, 0., 0., 0., 0., 0., 0., 0., lambda + 2*mu /) 
                    RM(2,:) = (/ lambda, 0., 0., 0., 1., 0., 0., 0., lambda /) 
                    RM(3,:) = (/ lambda, 0., 0., 0., 0., 1., 0., 0., lambda /) 
                    RM(4,:) = (/ 0.,     mu, 0., 0., 1., 0., 0., mu, 0.     /) 
                    RM(5,:) = (/ 0.,     0., 0., 1., 0., 0., 0., 0., 0.     /)  
                    RM(6,:) = (/ 0.,     0., mu, 0., 0., 0., mu, 0., 0.     /)  
                    RM(7,:) = (/ cp,     0., 0., 0., 0., 0., 0., 0., -cp    /)  
                    RM(8,:) = (/ 0.,     cs, 0., 0., 0., 0., 0.,-cs, 0.     /)  
                    RM(9,:) = (/ 0.,     0., cs, 0., 0., 0.,-cs, 0., 0.     /)                      
                    ! Rotate the inner state in the local coordinate system 
                    QLrot = MATMUL( TRANSPOSE(TM), Face(iFace)%qL(:,j,k))  
                    ! If the right neighbor is missing, the incoming characteristics have negative speed (lambda < 0), i.e. we need eigenvectors 1, 2 and 3    
                    ! The coefficients below solve the inverse Riemann problem so that the normal stress on the boundary is zero 
                    a1 = 0.5*(-QLrot(1)+rho*cp*QLrot(7))/rho/cp**2 
                    a2 = 0.5*(-QLrot(4)+rho*cs*QLrot(8))/rho/cs**2  
                    a3 = 0.5*(-QLrot(6)+rho*cs*QLrot(9))/rho/cs**2  
                    ! Compute the outer state, and rotate back to the physical x-y-z coordinate system 
                    QRrot = a1*RM(:,1) + a2*RM(:,2) + a3*RM(:,3)                     
#ifdef GODUNOV
                    ! For the exact Riemann solver, we can still use the old and sloppy mirror boundary condition. We now even have a proof for this on the paper. 
                    QRrot = QLrot 
                    QRrot((/1,4,6/)) = -QRrot((/1,4,6/)) 
#endif 
                    ! Put free surface boundaries everywhere 
                    !IF( ibc(1)==2 ) THEN
                        Face(iFace)%qR(:,j,k) = MATMUL( TM, QRrot ) 
                    !ELSE 
                    !    Face(iFace)%qR(:,j,k) = qbc 
                    !ENDIF 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO 
#elif defined(ACOUSTIC) 
    Qbc = (/ 1., 0., 0., 0. /)  
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use either simple extrapolation (copy from inside the domain) 
        ! or impose a constant value 
        !qbc = (/ 1., 0., 0., 0., 2.5 /) 
        IF(Face(iFace)%Left.EQ.0) THEN
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qL(:,j,k) = Face(iFace)%qR(:,j,k)
                    IF( ibc(1)==2 ) THEN
                        Face(iFace)%qL(3,j,k) = -Face(iFace)%qL(3,j,k)
                        !Face(iFace)%qL(:,j,k) = qbc 
                    ELSE 
                        Face(iFace)%qL(:,j,k) = qbc 
                    ENDIF 
                    !Face(iFace)%FL(:,j,k) = Face(iFace)%FR(:,j,k)
                    !Face(iFace)%qL(:,j,k) = qbc  
                    !CALL PDEFlux(Fbc,Face(iFace)%qL(:,j,k),Face(iFace)%paramL(:,j,k)) 
                    !Face(iFace)%FL(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qR(:,j,k) = Face(iFace)%qL(:,j,k) 
                    IF( ibc(1)==2 ) THEN
                        Face(iFace)%qR(3,j,k) = -Face(iFace)%qR(3,j,k) 
                        !Face(iFace)%qR(:,j,k) = qbc 
                    ELSE 
                        Face(iFace)%qR(:,j,k) = qbc 
                    ENDIF 
                    !Face(iFace)%FR(:,j,k) = Face(iFace)%FL(:,j,k)
                    !Face(iFace)%qR(:,j,k) = qbc  
                    !CALL PDEFlux(Fbc,Face(iFace)%qR(:,j,k),Face(iFace)%paramR(:,j,k)) 
                    !Face(iFace)%FR(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO 
#elif defined(Z4EINSTEIN) || defined(SRMHD) || defined(GRMHD) || defined(Z4GRMHD) || defined(CCZ4GRMHD) || defined(CCZ4EINSTEIN) || defined(GPR3D) || defined(EULER)      
    ! Fix boundary data  
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use either simple extrapolation (copy from inside the domain) 
        ! or impose a constant value 
        nv = Face(iFace)%nv    
        IF(Face(iFace)%Left.EQ.0) THEN
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    IF(ABS(nv(1)-1.0) < 1e-14) THEN 
                        ! x direction 
                        xGP = Face(iFace)%x0 + (/ 0., xiGPN(j), xiGPN(k) /)*dx 
                    ENDIF    
                    IF(ABS(nv(2)-1.0) < 1e-14) THEN
                        ! y direction 
                        xGP = Face(iFace)%x0 + (/ xiGPN(j), 0., xiGPN(k) /)*dx 
                    ENDIF    
                    IF(ABS(nv(3)-1.0) < 1e-14) THEN
                        ! z direction 
                        xGP = Face(iFace)%x0 + (/ xiGPN(j), xiGPN(k), 0. /)*dx 
                    ENDIF 
                    ! Compute the time integral of the boundary condition 
                    Face(iFace)%qL(:,j,k) = 0.0 
                    Face(iFace)%FL(:,j,k) = 0.0 
                    DO i = 1, N+1
                        CALL InitialField(QGP,parGP,xGP,time+xiGPN(i)*dt) 
                        !QGP = 2*QGP - Face(iFace)%qR(:,j,k)    ! set boundary state and flux to two times what we want minus what we have. (penalty on boundary jumps) 
                        Face(iFace)%qL(:,j,k) = Face(iFace)%qL(:,j,k) + wGPN(i)*QGP 
                        CALL PDEFlux(Fbc,QGP,parGP) 
                        Face(iFace)%FL(:,j,k) = Face(iFace)%FL(:,j,k) + wGPN(i)*MATMUL(Fbc(:,:), Face(iFace)%nv) 
                    ENDDO 
                    CONTINUE
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    IF(ABS(nv(1)-1.0) < 1e-14) THEN 
                        ! x direction 
                        xGP = Face(iFace)%x0 + (/ 0., xiGPN(j), xiGPN(k) /)*dx 
                    ENDIF    
                    IF(ABS(nv(2)-1.0) < 1e-14) THEN
                        ! y direction 
                        xGP = Face(iFace)%x0 + (/ xiGPN(j), 0., xiGPN(k) /)*dx 
                    ENDIF    
                    IF(ABS(nv(3)-1.0) < 1e-14) THEN
                        ! z direction 
                        xGP = Face(iFace)%x0 + (/ xiGPN(j), xiGPN(k), 0. /)*dx 
                    ENDIF 
                    ! Compute the time integral of the boundary condition 
                    Face(iFace)%qR(:,j,k) = 0.0 
                    Face(iFace)%FR(:,j,k) = 0.0 
                    DO i = 1, N+1
                        CALL InitialField(QGP,parGP,xGP,time+xiGPN(i)*dt) 
                        !QGP = 2*QGP - Face(iFace)%qL(:,j,k)    ! set boundary state and flux to two times what we want minus what we have. (penalty on boundary jumps) 
                        Face(iFace)%qR(:,j,k) = Face(iFace)%qR(:,j,k) + wGPN(i)*QGP 
                        CALL PDEFlux(Fbc,QGP,parGP) 
                        Face(iFace)%FR(:,j,k) = Face(iFace)%FR(:,j,k) + wGPN(i)*MATMUL(Fbc(:,:), Face(iFace)%nv) 
                    ENDDO 
                    CONTINUE
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO    
#else
    ! Fix boundary data  
    !Vbc = (/ 1., 1., 1., 0., 1. /)    ! primitive variables     
    CALL PDEPrim2Cons(qBC,Vbc)        ! convert into conservative variables    
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use either simple extrapolation (copy from inside the domain) 
        ! or impose a constant value 
        IF(Face(iFace)%Left.EQ.0) THEN
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qL(:,j,k) = Face(iFace)%qR(:,j,k)
                    !Face(iFace)%qL(2:4,j,k) = -Face(iFace)%qL(2:4,j,k)             ! mirror the entire velocity vector (which means cheating for Euler, but is correct for Navier-Stokes) 
                    Face(iFace)%qL(2:4,j,k) = Face(iFace)%qR(2:4,j,k) - 2*DOT_PRODUCT( Face(iFace)%nv, Face(iFace)%qR(2:4,j,k) ) * Face(iFace)%nv 
                    !Face(iFace)%FL(:,j,k) = Face(iFace)%FR(:,j,k)
                    !Face(iFace)%qL(:,j,k) = qbc 
                    CALL PDEFlux(Fbc,Face(iFace)%qL(:,j,k),Face(iFace)%paramL(:,j,k)) 
                    Face(iFace)%FL(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qR(:,j,k)   = Face(iFace)%qL(:,j,k)
                    !Face(iFace)%qR(2:4,j,k) = -Face(iFace)%qR(2:4,j,k)             ! mirror the entire velocity vector (which means  cheating for Euler, but is correct for Navier-Stokes) 
                    Face(iFace)%qR(2:4,j,k) = Face(iFace)%qL(2:4,j,k) - 2*DOT_PRODUCT( Face(iFace)%nv, Face(iFace)%qL(2:4,j,k) ) * Face(iFace)%nv 
                    !Face(iFace)%FR(:,j,k) = Face(iFace)%FL(:,j,k)
                    !Face(iFace)%qR(:,j,k) = qbc 
                    CALL PDEFlux(Fbc,Face(iFace)%qR(:,j,k),Face(iFace)%paramR(:,j,k)) 
                    Face(iFace)%FR(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO    

#endif 
    ! 
END SUBROUTINE BoundaryConditions 
    
    