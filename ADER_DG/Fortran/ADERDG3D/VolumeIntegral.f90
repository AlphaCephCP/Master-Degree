SUBROUTINE ADERVolumeIntegral(lduh,lqhi,lFhi,lShi)
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)  :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))      ! space-time degrees of freedom 
    REAL, INTENT(IN)  :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))    ! nonlinear flux tensor in each space DOF 
    REAL, INTENT(IN)  :: lShi(nVar,nDOF(1),nDOF(2),nDOF(3))      ! nonlinear source vector in each space DOF 
    REAL, INTENT(OUT) :: lduh(nVar,nDOF(1),nDOF(2),nDOF(3))      ! spatial degrees of freedom 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d) 
    !
#ifdef LINEAR  
    ! for linear non-conservative PDE, the volume integral is trivial, since it only involves the element mass matrix, which later will cancel 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            DO i = 1, nDOF(1) 
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /) 
                lduh(:,i,j,k) = lShi(:,i,j,k) * PRODUCT(aux(1:nDim)) 
            ENDDO
        ENDDO
    ENDDO 
    CONTINUE
#else
    ! 
    ! Initialize the update DOF 
    ! 
#if defined(NONCONSERVATIVE) || defined(SOURCE) 
    DO k = 1, nDOF(3) 
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ wGPN(i), wGPN(j), wGPN(k) /)  
                lduh(:,i,j,k) = PRODUCT(aux(1:nDim))*lShi(:,i,j,k)  
            ENDDO
        ENDDO
    ENDDO 
#else
    lduh = 0.0 
#endif 
    ! 
#ifndef NOFLUX 
    ! x - direction 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1., wGPN(j), wGPN(k) /) 
            lduh(:,:,j,k) = lduh(:,:,j,k) + MATMUL( lFhi(:,1,:,j,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(1) 
        ENDDO
    ENDDO
    ! y - direction (not needed in 1D) 
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(k) /) 
                lduh(:,i,:,k) = lduh(:,i,:,k) + MATMUL( lFhi(:,2,i,:,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(2)
            ENDDO
        ENDDO
    ENDIF 
    ! z - direction (not needed in 1D and 2D) 
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ 1., wGPN(i), wGPN(j) /)  
                lduh(:,i,j,:) = lduh(:,i,j,:) + MATMUL( lFhi(:,3,i,j,:), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(3)  
            ENDDO
        ENDDO
    ENDIF
#endif 
    !
#ifdef KREISSOLIGER 
    ! x - direction 
    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1., wGPN(j), wGPN(k) /) 
            lduh(:,:,j,k) = lduh(:,:,j,k) - epsilon_ko*MATMUL( lqhi(:,:,j,k), TRANSPOSE(Kxixi) )*PRODUCT(aux(1:nDim)) 
        ENDDO
    ENDDO
    ! y - direction (not needed in 1D) 
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1., wGPN(i), wGPN(k) /) 
                lduh(:,i,:,k) = lduh(:,i,:,k) - epsilon_ko*MATMUL( lqhi(:,i,:,k), TRANSPOSE(Kxixi) )*PRODUCT(aux(1:nDim)) 
            ENDDO
        ENDDO
    ENDIF 
    ! z - direction (not needed in 1D and 2D) 
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ 1., wGPN(i), wGPN(j) /)  
                lduh(:,i,j,:) = lduh(:,i,j,:) - epsilon_ko*MATMUL( lqhi(:,i,j,:), TRANSPOSE(Kxixi) )*PRODUCT(aux(1:nDim))   
            ENDDO
        ENDDO
    ENDIF
#endif 
    !
#endif
    !
    CONTINUE
    !
END SUBROUTINE ADERVolumeIntegral 
    
    