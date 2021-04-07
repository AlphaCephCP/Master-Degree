SUBROUTINE BaseFunc1D(phi,phi_xi,xi)
   USE typesDef
   IMPLICIT NONE
   ! Argument list 
   REAL, INTENT(IN ) :: xi                              ! coordinate in [0,1] where to evaluate the basis 
   REAL, INTENT(OUT) :: phi(N+1), phi_xi(N+1)           ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER     :: i,j,m
   REAL        :: tmp   
   ! 
   ! Initialize variables 
   phi      = 1. 
   phi_xi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO m = 1, N+1
      DO j = 1, N+1
         IF(j.EQ.m) CYCLE 
         phi(m) = phi(m)*(xi-xin(j))/(xin(m)-xin(j))    
      ENDDO 
      DO i = 1, N+1
         IF(i.EQ.m) CYCLE
         tmp = 1. 
         DO j = 1, N+1
            IF(j.EQ.i) CYCLE 
            IF(j.EQ.m) CYCLE 
            tmp = tmp*(xi-xin(j))/(xin(m)-xin(j))    
         ENDDO 
         phi_xi(m) = phi_xi(m) + tmp/(xin(m)-xin(i)) 
      ENDDO 
   ENDDO 
   !
END SUBROUTINE BaseFunc1D

SUBROUTINE BaseFunc1Dxx(phi,phi_xi,phi_xixi,xi)
   USE typesDef
   IMPLICIT NONE
   ! Argument list 
   REAL, INTENT(IN ) :: xi                              ! coordinate in [0,1] where to evaluate the basis 
   REAL, INTENT(OUT) :: phi(N+1), phi_xi(N+1)           ! the basis and its derivative w.r.t. xi 
   REAL, INTENT(OUT) :: phi_xixi(N+1)                   ! the second derivative of the basis function 
   ! Local variables 
   INTEGER     :: i,j,k,m
   REAL        :: tmp   
   ! 
   ! Initialize variables 
   phi        = 1. 
   phi_xi     = 0. 
   phi_xixi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO i = 1, N+1
      DO j = 1, N+1
         IF(j.EQ.i) CYCLE 
         phi(i) = phi(i)*(xi-xin(j))/(xin(i)-xin(j)) 
      ENDDO 
      DO m = 1, N+1
         IF(i.EQ.m) CYCLE 
         tmp = 1.0 
         DO j = 1, N+1
             IF(j==i) CYCLE
             IF(j==m) CYCLE
             tmp = tmp*(xi-xin(j))/(xin(i)-xin(j)) 
         ENDDO 
         phi_xi(i) = phi_xi(i) + tmp/(xin(i)-xin(m)) 
         DO k = 1, N+1
             IF(k==i) CYCLE
             IF(k==m) CYCLE
             tmp = 1.0 
             DO j = 1, N+1
                 IF(j==i) CYCLE
                 IF(j==m) CYCLE
                 IF(j==k) CYCLE
                 tmp = tmp*(xi-xin(j))/(xin(i)-xin(j))  
             ENDDO          
             phi_xixi(i) = phi_xixi(i) + tmp/(xin(i)-xin(m))/(xin(i)-xin(k)) 
         ENDDO         
      ENDDO 
   ENDDO 
   !
END SUBROUTINE BaseFunc1Dxx    
    

    
SUBROUTINE TimeBaseFunc1D(theta,theta_tau,tau)
   USE typesDef
   IMPLICIT NONE
   ! Argument list 
   REAL, INTENT(IN ) :: tau                              ! coordinate in [0,1] where to evaluate the basis 
   REAL, INTENT(OUT) :: theta(N+1), theta_tau(N+1)           ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER     :: i,j,m
   REAL        :: tmp   
   ! 
   ! Simple Taylor basis 
   theta(1) = 1.
   DO m = 2, N+1
     theta(m) = theta(m-1)/REAL(m-1)*tau 
   ENDDO 
   ! 
   theta_tau(1) = 0. 
   IF(N>0) THEN
    theta_tau(2) = 1.
    DO m = 3, N+1
        theta_tau(m) = theta_tau(m-1)/REAL(m-2)*tau 
    ENDDO  
   ENDIF 
   !
END SUBROUTINE TimeBaseFunc1D
    
SUBROUTINE Mm1BaseFunc1D(phi,phi_xi,xi)
   USE typesDef
   IMPLICIT NONE
   ! Argument list 
   REAL, INTENT(IN ) :: xi                              ! coordinate in [0,1] where to evaluate the basis 
   REAL, INTENT(OUT) :: phi(N), phi_xi(N)               ! the basis and its derivative w.r.t. xi 
   ! Local variables 
   INTEGER     :: i,j,m
   REAL        :: tmp   
   ! 
   ! Initialize variables 
   phi      = 1. 
   phi_xi   = 0. 
   ! Lagrange polynomial and its derivative 
   DO m = 1, N
      DO j = 1, N
         IF(j.EQ.m) CYCLE 
         phi(m) = phi(m)*(xi-xiGPMm1(j))/(xiGPMm1(m)-xiGPMm1(j))    
      ENDDO 
      DO i = 1, N
         IF(i.EQ.m) CYCLE
         tmp = 1. 
         DO j = 1, N
            IF(j.EQ.i) CYCLE 
            IF(j.EQ.m) CYCLE 
            tmp = tmp*(xi-xiGPMm1(j))/(xiGPMm1(m)-xiGPMm1(j))    
         ENDDO 
         phi_xi(m) = phi_xi(m) + tmp/(xiGPMm1(m)-xiGPMm1(i)) 
      ENDDO 
   ENDDO 
   !
END SUBROUTINE Mm1BaseFunc1D
    
    