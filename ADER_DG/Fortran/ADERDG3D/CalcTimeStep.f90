SUBROUTINE CalcTimeStep
    USE typesDef
    IMPLICIT NONE
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif    
    ! Local variables
    INTEGER :: iElem, iDim, i, j, k  
    REAL    :: lmax(d), Lambda(nVar), dtE, MPIdt  
    REAL    :: nv(d,d), denom, lx0 !(d)
    ! Normal vectors pointing into the three space dimensions 
    nv = 0. 
    DO i = 1, d
        nv(i,i) = 1. 
    ENDDO
    !
    ! This function computes the maximum admissible time step according to the CFL condition for a generic nonlinear PDE 
    !
    dt = 1e20 
    DO iElem = 1, nElem
#ifdef GRMHD
      lx0 = MAXVAL(ABS(x(:,tri(:,iElem))))
      IF( lx0.LT.ExcisionRadius ) THEN
          continue
          CYCLE
      ENDIF 
#endif
        IF(Limiter(iElem)%status==0) THEN
            DO k = 1, nDOF(3) 
             DO j = 1, nDOF(2) 
              DO i = 1, nDOF(1) 
                denom = 0.   
                DO iDim = 1, nDim
                    CALL PDEEigenvalues(Lambda,uh(:,i,j,k,iElem),parh(:,i,j,k,iElem),nv(:,iDim)) 
                    denom = denom + MAXVAL(ABS(Lambda))/dx(iDim) 
                ENDDO        
                dt = MIN(dt, CFL*PNPMTable(N)/denom )    
              ENDDO
             ENDDO
            ENDDO        
        ELSE
            DO k = 1, nSubLimV(3) 
             DO j = 1, nSubLimV(2) 
              DO i = 1, nSubLimV(1) 
                denom = 0.   
                DO iDim = 1, nDim
                    CALL PDEEigenvalues(Lambda,Limiter(iElem)%Lh(:,i,j,k),parh(:,i,j,k,iElem),nv(:,iDim)) 
                    denom = denom + MAXVAL(ABS(Lambda))/dx(iDim) 
                ENDDO        
                dt = MIN(dt, CFL*PNPMTable(N)/denom )    
              ENDDO
             ENDDO
            ENDDO        
        ENDIF        
    ENDDO
    !
#ifdef PARALLEL 
  CALL MPI_ALLREDUCE(dt,MPIdt,1,MPI_AUTO_REAL,MPI_MIN,MPI_COMM_WORLD,mpiErr) 
  dt = MPIdt 
#endif 
    ! 
    dt = MIN(dt, maxdt) 
    IF(time+dt>tend) THEN
        dt = tend - time
    ENDIF 
    IF(time+dt>tio) THEN
        dt = tio - time
    ENDIF 
    !    !
    CONTINUE 
    !
END SUBROUTINE CalcTimeStep 
    