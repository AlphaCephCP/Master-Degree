SUBROUTINE WriteData 
    USE MainVariables
    IMPLICIT NONE
    CHARACTER(LEN=500) :: FileName, VarString
    CHARACTER(LEN=10)  :: VarName(nVar)    
    INTEGER :: i,j,k 
    REAL    :: xGP(d) 
    !
#ifdef NOIO
    PRINT *, ' No output into file! ' 
    RETURN
#endif     
    !    
    WRITE(FileName,'(a,a,i7.7,a)') TRIM(BaseFile),'-',iter,'.dat' 
    PRINT *, ' Writing data to file ', TRIM(FileName) 
    !
    OPEN(UNIT=333,FILE=TRIM(FileName),RECL=1000)
    CALL PDEVarName(VarName) 
    WRITE(VarString,'(a)') ' VARIABLES = "x" "y" "z" '  
    ! 
    DO i = 1, nVar
        WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' "', TRIM(VarName(i)), '" ' 
    ENDDO    
    WRITE(333,*) TRIM(VarString) 
    WRITE(333,*) ' ZONE I = ', IMAX, ' J = ', JMAX, ' K = ', KMAX, ' DATAPACKING=POINT ' 
    
    DO k = 1, KMAX 
     DO j = 1, JMAX 
      DO i = 1, IMAX
         xGP = xL + (/ i-1, j-1, k-1 /)*dx 
         WRITE(333,*) xGP(1), xGP(2), xGP(3), uh(:,i,j,k) 
      ENDDO
     ENDDO 
    ENDDO    
    ! 
    CLOSE(333) 
    !
END SUBROUTINE WriteData
    