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
    
        