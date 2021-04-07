SUBROUTINE InitPointSources 
    USE typesDef
    IMPLICIT NONE
    ! Argument list 
    ! Local variables 
    INTEGER :: i, ii,jj,kk,idxs(d) 
    REAL    :: aux(d)
    REAL    :: phi1(N+1), phi1_xi(N+1)  
    REAL    :: phi2(N+1), phi2_xi(N+1)  
    REAL    :: phi3(N+1), phi3_xi(N+1)  

#ifdef ELASTICITY 
    nPointSource = 1 
    ALLOCATE( PointSrc(nPointSource) ) 
    IF(nPointSource>0) THEN
        PointSrc(1)%xP = (/ 2000., 1950.0, 1000. /) 
!        PointSrc(1)%xP = (/ 0.0, 0.0, 0.0 /) 
        PointSrc(1)%waveform = 3  
    ENDIF 
    DO i = 1, nPointSource
        ! Locate the sources on the Cartesian :-) grid
        idxs = CEILING( (PointSrc(i)%xP - xL)/dx ) 
        PointSrc(i)%iElem = idxe( idxs(1), idxs(2), idxs(3) ) 
        PointSrc(i)%xiP =  (PointSrc(i)%xP - xL)/dx - REAL( idxs - 1 ) 
        ALLOCATE( PointSrc(i)%phi(nDOF(1),nDOF(2),nDOF(3)) )   
        ALLOCATE( PointSrc(i)%sigma(1:N+1,1:nVar), PointSrc(i)%SrcInt(1:nVar) )    
        CALL BaseFunc1D(phi1,phi1_xi,PointSrc(i)%xiP(1))
        CALL BaseFunc1D(phi2,phi2_xi,PointSrc(i)%xiP(2))
        CALL BaseFunc1D(phi3,phi3_xi,PointSrc(i)%xiP(3))
        DO kk = 1, nDOF(3)
         DO jj = 1, nDOF(2) 
          DO ii = 1, nDOF(1)
            aux = (/ phi1(ii), phi2(jj), phi3(kk) /) 
            PointSrc(i)%phi(ii,jj,kk) = PRODUCT(aux(1:nDim)) 
          ENDDO
         ENDDO
        ENDDO 
        CONTINUE
    ENDDO 
#else
    nPointSource = 0    
#endif 
END SUBROUTINE InitPointSources 

SUBROUTINE RunPointSources 
    USE typesDef
    IMPLICIT NONE
    ! Local variables 
    INTEGER :: i, l, iGP, iVar, ii,jj,kk,idxs(d) 
    REAL    :: aux(d), sigma(nVar), tGP  
    REAL    :: theta(N+1), theta_tau(N+1)  
    REAL    :: thetasigma(N+1,nVar), dtavFac  

    DO i = 1, nPointSource
        IF(PointSrc(i)%iElem<=0) THEN
            CYCLE
        ENDIF 
        thetasigma = 0. 
        PointSrc(i)%SrcInt(:) = 0.0 
        DO iGP = 1, N+1
            tGP = time + dt*xiGPN(iGP) 
            CALL SourceTimeFunction(sigma,tGP,PointSrc(i)%waveform) 
            CALL TimeBaseFunc1D(theta,theta_tau,xiGPN(iGP)) 
            DO iVar = 1, nVar
                thetasigma(:,iVar) = thetasigma(:,iVar) + wGPN(iGP)*theta(:)*sigma(iVar)
            ENDDO  
            PointSrc(i)%SrcInt(:) = PointSrc(i)%SrcInt(:) + wGPN(iGP)*sigma(:)
        ENDDO 
        PointSrc(i)%sigma(:,:) = MATMUL( iMT, thetasigma(:,:) ) 

!        PointSrc(i)%SrcInt(:) = PointSrc(i)%sigma(1,:) 
!        dtavFac = 0.5*1  
!        DO l = 2, N+1 
!            PointSrc(i)%SrcInt(:) = PointSrc(i)%SrcInt(:)  + dtavFac*PointSrc(i)%sigma(l,:)
!            dtavFac = dtavFac*1/REAL(l+1)
!        ENDDO
        CONTINUE
    ENDDO 

END SUBROUTINE RunPointSources 

SUBROUTINE AddPointSources  
    USE typesDef
    IMPLICIT NONE
    ! Local variables 
    INTEGER :: i, iVar, iElem, l, iGP, ii,jj,kk,idxs(d) 
    !
    DO i = 1, nPointSource
        IF(PointSrc(i)%iElem<=0) THEN
            CYCLE
        ENDIF 
        iElem = PointSrc(i)%iElem 
        DO kk = 1, nDOF(3) 
         DO jj = 1, nDOF(2) 
          DO ii = 1, nDOF(1) 
            DO iVar = 1, nVar
                duh(iVar,ii,jj,kk,iElem) = duh(iVar,ii,jj,kk,iElem) + PointSrc(i)%phi(ii,jj,kk)*PointSrc(i)%SrcInt(iVar)/PRODUCT(dx(1:nDim))  
            ENDDO
          ENDDO
         ENDDO  
        ENDDO 
        !
        CONTINUE
        !
    ENDDO 
    !
END SUBROUTINE AddPointSources 