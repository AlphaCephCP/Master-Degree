!
! Small and simple Gauss package
! written by Michael Dumbser
! 27.05.2011 
!
SUBROUTINE LinSolve(N,A,b,x)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), b(N), x(N)  
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: temp(N+1),piv 
  REAL          :: C(N+1,N)  
  !
  C(1:N,:)     = TRANSPOSE(A)
  C(N+1,:)     = b 
  !    
  ! Forward elimination with column pivoting 
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        DO j = 1, N
           PRINT *, A(j,:) 
        ENDDO
        STOP
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  x = C(N+1,:)  
  !
END SUBROUTINE LinSolve

SUBROUTINE LinSolve2(N,M,A,b,x,iErr)
   IMPLICIT NONE
   INTEGER       :: N,M
   INTEGER       :: iErr 
   REAL          :: A(N,N), b(N,M), x(N,M)
   !
   INTEGER       :: i,j,ml(1)
   REAL          :: temp(N+M),piv
   REAL          :: C(N+M,N)
   !
   iErr = 0 
   ! 
   C(1:N,:)      = TRANSPOSE(A)
   C(N+1:N+M,:)  = TRANSPOSE(b)
   !
   ! Forward elimination with column pivoting
   !
   DO i = 1, N
     ! If pivot element is zero, then swap rows
     ml = MAXLOC(ABS(C(i,i:N)))
     j = i - 1 + ml(1)
     temp   = C(:,j)
     C(:,j) = C(:,i)
     C(:,i) = temp
     IF(C(i,i).EQ.0.) THEN
        !IF(PRESENT(iErr)) THEN
            PRINT *, 'ERROR. Matrix is singular!'
            iErr = 1 
            RETURN 
        !ELSE
        !    PRINT *, 'ERROR. Matrix is singular!'
        !    DO j = 1, N
        !      PRINT *, A(j,:)
        !    ENDDO
        !    STOP
        !ENDIF 
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv
     DO j = i+1, N
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
   ENDDO
   !
   ! Back substitution
   !
   DO i = N,1,-1
    DO j = i-1,1,-1
      C(:,j) = C(:,j) - C(i,j)*C(:,i)
    ENDDO
   ENDDO
   !
   x = TRANSPOSE(C(N+1:N+M,:))
   !
END SUBROUTINE LinSolve2 
    
SUBROUTINE MatrixInverse(N,A,iA)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), iA(N,N)
  !
  INTEGER       :: i,j,flag,ml(1) 
  REAL          :: piv
  REAL          :: temp(2*N)
  REAL          :: C(2*N,N)  
  !
  C(1:N,:)     = TRANSPOSE(A)
  C(N+1:2*N,:) = 0. 
  DO i = 1, N
     C(N+i,i) = 1.
  ENDDO
  !    
  ! Forward elimination and row swapping (if necessary)
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        PRINT *, 'ERROR. Matrix is singular!'
        DO j = 1, N
           PRINT *, A(j,:) 
        ENDDO
        STOP
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  iA = TRANSPOSE( C(N+1:2*N,:) ) 
  !
END SUBROUTINE MatrixInverse


!
! Constrained L2 projection   
! 
! This subroutine solves the quadratic optimization problem associated with the 
! constrained L2 projection. 
!    A x = A b    (1) 
! and a set of linear constraints 
!   C x = R b     (2) 
! then the subroutine CLSQ will solve equation (1) in a least squares sense, respecting the 
! linear constraints (2) exactly. 
! This leads to a simple quadratic optimization problem with constraint, which is solved 
! using a standard Lagrangian multiplier technique. See the above reference for details. 
! The solution of the constrained system is given by x = M b, where the matrix iGH is 
! the output of CLSQ. 
! 
! Input: nEqn   = number of equations ( > nDOF ) 
!        nDOF   = number of unknowns  
!        nConst = number of constraints 
!        A      = the (nDOF,nDOF)   system matrix of (1) 
!        C      = the (nConst,nDOF) matrix of the linear constraints (2) 
!        R      = the (nConst,nEqn) matrix of the linear constraints (2) 
! 
! Output: M    = the (nDOF,nEqn) matrix, which yields the solution x when multiplied with the right hand side vector b 
!
SUBROUTINE CL2(M,nEqn,nDOF,nConst,A,C,R,iErr) 
   !USE operators_mod 
   IMPLICIT NONE 
   ! Argument list 
   INTEGER :: nEqn, nDOF, nConst, iErr  
   REAL    :: M(nDOF,nEqn), A(nDOF,nDOF), C(nConst,nDOF), R(nConst,nEqn) 
   ! Local variables 
   INTEGER :: i 
   REAL    :: H(nDOF+nConst,nEqn), iGH(nDOF+nConst,nEqn), G(nDOF+nConst,nDOF+nConst) ! iG(nDOF+nConst,nDOF+nConst), Id(nDOF+nConst,nDOF+nConst) 
   ! 
   iErr = 0 
   ! Build the right hand side matrix of the associated optimization problem 
   H = 0. 
   H(1:nDOF,1:nDOF)             = 2.0*A              ! right hand side matrix of the associated unconstrained system of normal equations  
   H(nDOF+1:nDOF+nConst,1:nEqn) = R                  ! right hand side matrix of the linear constraints 
   ! 
   ! Build the left hand side matrix of the associated optimization problem 
   G = 0. 
   G(1:nDOF,1:nDOF) = 2.0*A                          ! left hand side matrix of the associated unconstrained system of normal equations   
   G(1:nDOF,nDOF+1:nDOF+nConst) = -TRANSPOSE(C)      ! transpose of the left hand side matrix of the linear constraints 
   G(nDOF+1:nDOF+nConst,1:nDOF) = C                  ! left hand side matrix of the linear constraints  
   !
   CALL LinSolve2(nDOF + nConst,nEqn,G,H,iGH,iErr)        ! Solve the final saddle point problem 
   ! 
   M = iGH(1:nDOF,1:nEqn)                            ! We are only interested in the solution vector x and NOT in the solution for the Lagrange multiplier
   ! 
END SUBROUTINE CL2 

    
    
    