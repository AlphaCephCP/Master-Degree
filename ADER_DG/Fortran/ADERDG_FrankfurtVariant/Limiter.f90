!
! Compute the subcell data lim from a given DG polynomial luh 
!
SUBROUTINE GetSubcellData(lim,luh) 
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   REAL :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)), luh(nVar,nDOF(1),nDOF(2),nDOF(3)) 
   ! Local variables 
   INTEGER :: iii, jjj, kkk
   REAL :: limz(nVar,nDOF(1),nDOF(2),nSubLimV(3)), limy(nVar,nDOF(1),nSubLimV(2),nSubLimV(3))
   !
    IF (nDim == 3) THEN
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                limz(:,iii,jjj,:) = MATMUL( luh(:,iii,jjj,:), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        limz = 0. 
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                limz(:,iii,jjj,1) = luh(:,iii,jjj,1) 
            ENDDO
        ENDDO
    ENDIF
    ! mapping of uh to the sub-cell limiter along y direction
    IF( nDim >= 2 ) THEN
        limy = 0. 
        DO jjj = 1, nSubLimV(3)  ! z 
            DO iii = 1, nDOF(1) ! x
                limy(:,iii,:,jjj) = MATMUL( limz(:,iii,:,jjj), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        limy = 0. 
        DO iii = 1, nDOF(1) ! x
            limy(:,iii,1,1) = luh(:,iii,1,1) 
        ENDDO
    ENDIF    
    ! mapping of uh to the sub-cell limiter along x direction 
    lim = 0. 
    DO jjj = 1, nSubLimV(3)  ! z 
        DO iii = 1, nSubLimV(2) ! y
            lim(:,:,iii,jjj) = MATMUL( limy(:,:,iii,jjj), TRANSPOSE(uh2lim) )  
        ENDDO
    ENDDO   
    !
END SUBROUTINE GetSubcellData  
!
!
! Compute the subcell data of the parameters from a given DG polynomial lparh 
!
SUBROUTINE GetSubcellParam(subpar,lparh) 
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   REAL :: subpar(nParam,nSubLimV(1),nSubLimV(2),nSubLimV(3)), lparh(nParam,nDOF(1),nDOF(2),nDOF(3)) 
   ! Local variables 
   INTEGER :: iii, jjj, kkk
   REAL :: parz(nParam,nDOF(1),nDOF(2),nSubLimV(3)), pary(nParam,nDOF(1),nSubLimV(2),nSubLimV(3))
   !
    IF (nDim == 3) THEN
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                parz(:,iii,jjj,:) = MATMUL( lparh(:,iii,jjj,:), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        parz = 0. 
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                parz(:,iii,jjj,1) = lparh(:,iii,jjj,1) 
            ENDDO
        ENDDO
    ENDIF
    ! mapping of uh to the sub-cell limiter along y direction
    IF( nDim >= 2 ) THEN
        pary = 0. 
        DO jjj = 1, nSubLimV(3)  ! z 
            DO iii = 1, nDOF(1) ! x
                pary(:,iii,:,jjj) = MATMUL( parz(:,iii,:,jjj), TRANSPOSE(uh2lim) ) 
            ENDDO
        ENDDO
    ELSE
        pary = 0. 
        DO iii = 1, nDOF(1) ! x
            pary(:,iii,1,1) = lparh(:,iii,1,1) 
        ENDDO
    ENDIF    
    ! mapping of uh to the sub-cell limiter along x direction 
    subpar = 0. 
    DO jjj = 1, nSubLimV(3)  ! z 
        DO iii = 1, nSubLimV(2) ! y
            subpar(:,:,iii,jjj) = MATMUL( pary(:,:,iii,jjj), TRANSPOSE(uh2lim) )  
        ENDDO
    ENDDO   
    !
END SUBROUTINE GetSubcellParam
!
! Compute the point values lob in the Gauss-Lobatto points from a given DG polynomial luh 
!
SUBROUTINE GetLobattoData(lob,luh) 
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   REAL :: lob(nVar,nDOF(1),nDOF(2),nDOF(3)), luh(nVar,nDOF(1),nDOF(2),nDOF(3)) 
   ! Local variables 
   INTEGER :: iii, jjj, kkk
   REAL :: lobz(nVar,nDOF(1),nDOF(2),nDOF(3)), loby(nVar,nDOF(1),nDOF(2),nDOF(3))
   !
    IF (nDim == 3) THEN
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                lobz(:,iii,jjj,:) = MATMUL( luh(:,iii,jjj,:), TRANSPOSE(uh2lob) ) 
            ENDDO
        ENDDO
    ELSE
        lobz = 0. 
        DO jjj = 1, nDOF(2)  ! y
            DO iii = 1, nDOF(1) ! x
                lobz(:,iii,jjj,1) = luh(:,iii,jjj,1) 
            ENDDO
        ENDDO
    ENDIF
    ! mapping of uh to the Gauss-Lobatto points along y direction
    IF( nDim >= 2 ) THEN
        loby = 0. 
        DO jjj = 1, nDOF(3)  ! z 
            DO iii = 1, nDOF(1) ! x
                loby(:,iii,:,jjj) = MATMUL( lobz(:,iii,:,jjj), TRANSPOSE(uh2lob) ) 
            ENDDO
        ENDDO
    ELSE
        loby = 0. 
        DO iii = 1, nDOF(1) ! x
            loby(:,iii,1,1) = luh(:,iii,1,1) 
        ENDDO
    ENDIF    
    ! mapping of uh to the Gauss-Lobatto points along x direction 
    lob = 0. 
    DO jjj = 1, nDOF(3)  ! z 
        DO iii = 1, nDOF(2) ! y
            lob(:,:,iii,jjj) = MATMUL( loby(:,:,iii,jjj), TRANSPOSE(uh2lob) )  
        ENDDO
    ENDDO   
    !
END SUBROUTINE GetLobattoData  
!
SUBROUTINE PutSubcellData(luh,lim) 
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   REAL :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)), luh(nVar,nDOF(1),nDOF(2),nDOF(3)) 
   ! Local variables 
   INTEGER :: iii, jjj, kkk
   REAL :: luhz(nVar,nSubLimV(1),nSubLimV(2),N+1), luhy(nVar,nSubLimV(1),N+1,N+1)
   !
   IF (nDim >= 3) THEN
        ! mapping of the sub-cell limiter to uh along z direction
        DO jjj = 1, nSubLimV(2)  ! y
            DO iii = 1, nSubLimV(1) ! x
                luhz(:,iii,jjj,:) = MATMUL( lim(:,iii,jjj,:), TRANSPOSE(lim2uh) ) 
            ENDDO
        ENDDO
   ELSE
        DO jjj = 1, nSubLimV(2)  ! y
            DO iii = 1, nSubLimV(1) ! x
                luhz(:,iii,jjj,1) = lim(:,iii,jjj,1) 
            ENDDO
        ENDDO
   ENDIF   
   IF(nDim >= 2) THEN
        ! mapping of the sub-cell limiter to uh along y direction
        luhy = 0. 
        DO jjj = 1, nDOF(3)  ! z 
            DO iii = 1, nSubLimV(1) ! x
                luhy(:,iii,:,jjj) = MATMUL( luhz(:,iii,:,jjj), TRANSPOSE(lim2uh) ) 
            ENDDO
        ENDDO
   ELSE
        luhy = 0. 
        DO iii = 1, nSubLimV(1) ! x
            luhy(:,iii,1,1) = luhz(:,iii,1,1) 
        ENDDO
   ENDIF
   ! mapping of the sub-cell limiter to uh along x direction
   luh = 0. 
   DO jjj = 1, nDOF(3)  ! z 
        DO iii = 1, nDOF(2)  ! y
            luh(:,:,iii,jjj) = MATMUL( luhy(:,:,iii,jjj), TRANSPOSE(lim2uh) )  
        ENDDO
   ENDDO   
   !
END SUBROUTINE PutSubcellData  
!    
SUBROUTINE GetMinMax  
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   INTEGER :: reflev 
   ! Local variables 
   INTEGER :: i, j, ii, jj, kk, iVar 
   REAL :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) 
   REAL :: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
   REAL :: Qmin(nVar), Qmax(nVar) 
   !
   ! Step 1) Loop over all elements and get the min/max for each variable 
   !
   DO i = 1, nElem
    IF(Limiter(i)%status.EQ.0) THEN 
        CALL GetSubcellData( lim, uh(:,:,:,:,i) ) 
        CALL GetLobattoData( lob, uh(:,:,:,:,i) )
        DO iVar = 1, nVar 
            !Limiter(i)%lmin(iVar) = MIN( MINVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MINVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))) ) 
            !Limiter(i)%lmax(iVar) = MAX( MAXVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MAXVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))) ) 
            Limiter(i)%lmin(iVar) = MIN( MINVAL(uh(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3),i)), MINVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MINVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))) ) 
            Limiter(i)%lmax(iVar) = MAX( MAXVAL(uh(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3),i)), MAXVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MAXVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))) ) 
        ENDDO  
    ELSE
        lim = Limiter(i)%Lh 
        DO iVar = 1, nVar 
            Limiter(i)%lmin(iVar) = MINVAL( lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3)) ) 
            Limiter(i)%lmax(iVar) = MAXVAL( lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3)) )  
        ENDDO  
    ENDIF  
   ENDDO
   !
   ! Step 2) Loop over all Voronoi neighbors and get the final min/max for each variable (requires step 1 to be complete) 
   ! 
   DO i = 1, nElem 
    Qmin = Limiter(i)%lmin(:) 
    Qmax = Limiter(i)%lmax(:)
    DO kk = -dn(3), dn(3) 
     DO jj = -dn(2), dn(2) 
      DO ii = -dn(1), dn(1) 
        j = neighbor(ii,jj,kk,i) 
        DO iVar = 1, nVar 
            Qmin(iVar) = MIN( Qmin(iVar), Limiter(j)%lmin(iVar) ) 
            Qmax(iVar) = MAX( Qmax(iVar), Limiter(j)%lmax(iVar) ) 
        ENDDO  
       ENDDO
      ENDDO
     ENDDO
     Limiter(i)%lmin(:) = Qmin 
     Limiter(i)%lmax(:) = Qmax 
    ENDDO
    !
    CONTINUE 
    !
END SUBROUTINE GetMinMax   
!    
SUBROUTINE DMP(dmpresult,arguh,argLimiter,argmax)  
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   LOGICAL        :: dmpresult 
   REAL           :: arguh(nVar,nDOF(1),nDOF(2),nDOF(3)),argmax  
   TYPE(tLimiter) :: argLimiter
   ! Local variables 
   INTEGER :: i, iii, jjj, kkk, iVar, iErr 
   REAL    :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) 
   REAL    :: lob(nVar,nDOF(1),nDOF(2),nDOF(3)) 
   REAL    :: lmin(nVar), lmax(nVar), ldiff, Qout(nVar) 
   !
#ifdef ELASTICITY
   dmpresult = .TRUE. 
   RETURN
#endif    
   !dmpresult = .FALSE. 
   !return 
   dmpresult = .TRUE. 
   !return
   CALL GetSubcellData( lim, arguh )  
   CALL GetLobattoData( lob, arguh )  
   DO iVar = 1, nVar 
        lmin(iVar) = MIN( MINVAL(arguh(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MINVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))), MINVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))) )  
        lmax(iVar) = MAX( MAXVAL(arguh(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))), MAXVAL(lim(iVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3))), MAXVAL(lob(iVar,1:nDOF(1),1:nDOF(2),1:nDOF(3))) )  
        ldiff = MAX( 1e-4, 1e-3*(argLimiter%lmax(iVar) - argLimiter%lmin(iVar)), argmax )   ! *(argLimiter%lmax(iVar) - argLimiter%lmin(iVar))
        IF(lmin(iVar).LT.argLimiter%lmin(iVar)-ldiff) THEN       
            dmpresult = .FALSE. 
            RETURN
        ENDIF 
        IF(lmax(iVar).GT.argLimiter%lmax(iVar)+ldiff) THEN
            dmpresult = .FALSE. 
            RETURN
        ENDIF 
   ENDDO  
   ! Physical detection criteria 
   DO kkk = 1, nSubLimV(3)
    DO jjj = 1, nSubLimV(2)
     DO iii = 1, nSubLimV(1)
        CALL PDEAssurePositivity(iErr, lim(:,iii,jjj,kkk), Qout ) 
        IF(iErr.LT.0) THEN 
            dmpresult = .FALSE. 
            RETURN
        ENDIF         
        DO iVar = 1, nVar 
           IF(ISNAN(lim(iVar,iii,jjj,kkk))) THEN
                dmpresult = .FALSE. 
                RETURN
           ENDIF
        ENDDO  
     ENDDO
    ENDDO
   ENDDO 
   !
   CONTINUE 
   !
END SUBROUTINE DMP   
   
SUBROUTINE Saveolduh 
   USE typesDef     
   IMPLICIT NONE                                                               
   ! Local variables
   INTEGER  :: i 
   !
   DO i = 1, nElem 
      olduh(:,:,:,:,i) = uh(:,:,:,:,i)
   ENDDO 
   !
   CONTINUE
   !
END SUBROUTINE Saveolduh  

SUBROUTINE UpdateLimiter 
   USE typesDef
   IMPLICIT NONE                                                               
   ! Local variables
   INTEGER :: i 
   !
   DO i = 1, nElem
    IF(Limiter(i)%status.EQ.0) THEN
        IF(Limiter(i)%oldstatus.NE.0) THEN 
            DEALLOCATE( Limiter(i)%Lh    ) 
            DEALLOCATE( Limiter(i)%NewLh ) 
        ENDIF        
    ELSE
        Limiter(i)%Lh = Limiter(i)%NewLh 
    ENDIF 
    Limiter(i)%oldstatus = Limiter(i)%status  
    CONTINUE 
   ENDDO 
   !
END SUBROUTINE UpdateLimiter 

SUBROUTINE AllocateLimiter
   USE typesDef     
   IMPLICIT NONE                                                               
   ! Local variables
   INTEGER    :: i  
   !
   DO i = 1, nElem
       IF(recompute(i).EQ.1) THEN
           IF(Limiter(i)%oldstatus.EQ.0) THEN 
               ALLOCATE( Limiter(i)%Lh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))    ) 
               ALLOCATE( Limiter(i)%NewLh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) ) 
           ENDIF
           Limiter(i)%status = 1 
       ELSE
           Limiter(i)%status = 0 
       ENDIF       
   ENDDO 
   !
END SUBROUTINE AllocateLimiter    
    
SUBROUTINE SpreadRecompute 
   USE typesDef     
   IMPLICIT NONE                                                               
   ! Local variables
   INTEGER  :: i, j, ii, jj, kk
   INTEGER  :: iFace
   !
   DO i = 1, nElem 
    IF(recompute(i)==1) THEN
      DO iFace = 1, 2*nDim 
        ii = Face2Neigh(1,iFace) 
        jj = Face2Neigh(2,iFace) 
        kk = Face2Neigh(3,iFace) 
        j = neighbor(ii,jj,kk,i) 
        IF(recompute(j)==0) THEN
            recompute(j) = 2 
        ENDIF 
       ENDDO
     ENDIF 
    ENDDO
   !
   CONTINUE
   !
END SUBROUTINE SpreadRecompute      
    
SUBROUTINE Subcellrecompute    
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   INTEGER :: i 
   ! Local variables 
   INTEGER :: j, ii, jj, kk, iVar, reflev, iDim, iFace   
   INTEGER :: postatus
   REAL    :: ldx(d), xg(d), Wout(nVar), nv(d), nvi(d,d) 
   REAL    :: subpar(nParam,(1-nSubLimV(1)):2*nSubLimV(1),(1-nSubLimV(2)):2*nSubLimV(2),(1-nSubLimV(3)):2*nSubLimV(3)) 
   REAL    :: subuh0(nVar,(1-nSubLimV(1)):2*nSubLimV(1),(1-nSubLimV(2)):2*nSubLimV(2),(1-nSubLimV(3)):2*nSubLimV(3)) 
   REAL    :: subuh1(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) 
   REAL    :: slopex(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: slopey(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: slopez(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wLx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wLy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wLz(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wRx(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wRy(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: wRz(nVar,1-dn(1):nSubLimV(1)+dn(1),1-dn(2):nSubLimV(2)+dn(2),1-dn(3):nSubLimV(3)+dn(3)) 
   REAL    :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)), lpar(nParam,nSubLimV(1),nSubLimV(2),nSubLimV(3)), x00(d) 
   REAL    :: LL(nVar),LR(nVar),Qt(nVar),smax
   REAL    :: FLx(nVar,d), FRx(nVar,d), FLy(nVar,d), FRy(nVar,d), FLz(nVar,d), FRz(nVar,d) 
   REAL    :: Fx(nVar,1:nSubLimV(1)+1,1:nSubLimV(2),1:nSubLimV(3)) 
   REAL    :: Fy(nVar,1:nSubLimV(1),1:nSubLimV(2)+1,1:nSubLimV(3)) 
   REAL    :: Fz(nVar,1:nSubLimV(1),1:nSubLimV(2),1:nSubLimV(3)+1) 
   REAL    :: ARoep(nVar,nVar), ARoem(nVar,nVar), absADT(nVar,nVar), fp(nVar), fm(nVar), Dp(nVar), Dm(nVar), gradQ(nVar,3)  
   REAL    :: AQx(nVar), BQy(nVar), CQz(nVar) 
   REAL    :: Eh(nVar,nVar), Q0(nVar), Q1(nVar), dQ(nVar), Id(nVar,nVar), ff(nVar), src(nVar), Jac(nVar,nVar), Q2(nVar), delta, res, res2, tol = 1e-6    
   INTEGER :: iNewton, inner, MaxNewton = 100  
   LOGICAL :: firstorder, dmpresult 
   !
   Id = 0. 
   DO iVar = 1, nVar
       Id(iVar,iVar) = 1. 
   ENDDO 
   nvi = 0
   DO ii = 1, d 
       nvi(ii,ii) = 1
   ENDDO 
   ldx = dx/REAL(nSubLim)   
   Fx = 0. 
   Fy = 0. 
   Fz = 0.
   slopex = 0. 
   slopey = 0. 
   slopez = 0. 
   wLx = 0. 
   wRx = 0. 
   wLy = 0. 
   wRy = 0. 
   wLz = 0. 
   wRz = 0. 
   !
   DO i = 1, nElem
       IF(recompute(i) > 0) THEN
           ! Get the initial data, either from a previous limiter stage, or from healthy uh initial data 
           ! note that here we only need the values of the face neighbors. For simplicity, we run over all the Voronoi neighbors             
           DO kk = -dn(3), dn(3) 
            DO jj = -dn(2), dn(2) 
             DO ii = -dn(1), dn(1) 
 !          DO iFace = 1, 2*nDim 
!                ii = Face2Neigh(1,iFace) 
!                jj = Face2Neigh(2,iFace) 
!                kk = Face2Neigh(3,iFace) 
                j = neighbor(ii,jj,kk,i)
                IF(Limiter(j)%oldstatus.EQ.0) THEN
                    CALL GetSubcellData(lim,olduh(:,:,:,:,j)) 
                ELSE
                    lim = Limiter(j)%Lh 
                ENDIF                
                subuh0(:,1+ii*nSubLimV(1):(ii+1)*nSublimV(1),1+jj*nSubLimV(2):(jj+1)*nSublimV(2),1+kk*nSubLimV(3):(kk+1)*nSublimV(3)) = lim(:,:,:,:)    
                IF(nParam > 0) THEN
                    CALL GetSubcellParam(lpar,parh(:,:,:,:,j)) 
                    subpar(:,1+ii*nSubLimV(1):(ii+1)*nSublimV(1),1+jj*nSubLimV(2):(jj+1)*nSublimV(2),1+kk*nSubLimV(3):(kk+1)*nSublimV(3)) = lpar(:,:,:,:)    
                ENDIF       
             ENDDO 
            ENDDO         
           ENDDO 
           ! -------------------------------------------------------------- 
           ! Second order TVD MUSCL reconstruction in space and time 
           ! -------------------------------------------------------------- 
           DO kk = 1, nSubLimV(3)
            DO jj = 1, nSubLimV(2)                 
             DO ii = 1, nSubLimV(1) 
                 ! x slopes
                 CALL minmod( slopex(:,ii,jj,kk), subuh0(:,ii+1,jj,kk)-subuh0(:,ii,jj,kk), subuh0(:,ii,jj,kk)-subuh0(:,ii-1,jj,kk), nVar ) 
                 IF(nDim>=2) THEN
                    ! y slopes
                    CALL minmod( slopey(:,ii,jj,kk), subuh0(:,ii,jj+1,kk)-subuh0(:,ii,jj,kk), subuh0(:,ii,jj,kk)-subuh0(:,ii,jj-1,kk), nVar ) 
                 ENDIF
                 IF(nDim>=3) THEN
                    ! z slopes
                    CALL minmod( slopez(:,ii,jj,kk), subuh0(:,ii,jj,kk+1)-subuh0(:,ii,jj,kk), subuh0(:,ii,jj,kk)-subuh0(:,ii,jj,kk-1), nVar ) 
                 ENDIF                 
             ENDDO
            ENDDO
           ENDDO 
           ! Time evolution and boundary extrapolated data 
           DO kk = 1-dn(3), nSubLimV(3)+dn(3)
            DO jj = 1-dn(2), nSubLimV(2)+dn(2)                 
              DO ii = 1-dn(1), nSubLimV(1)+dn(1) 
                  wLx(:,ii,jj,kk) = subuh0(:,ii,jj,kk) - 0.5*slopex(:,ii,jj,kk) 
                  wRx(:,ii,jj,kk) = subuh0(:,ii,jj,kk) + 0.5*slopex(:,ii,jj,kk) 
                  CALL PDEFlux(FLx,wLx(:,ii,jj,kk),subpar(:,ii,jj,kk))
                  CALL PDEFlux(FRx,wRx(:,ii,jj,kk),subpar(:,ii,jj,kk))
                  Qt = - (FRx(:,1)-FLx(:,1))/ldx(1) 
                  IF(nDim>=2) THEN
                      wLy(:,ii,jj,kk) = subuh0(:,ii,jj,kk) - 0.5*slopey(:,ii,jj,kk) 
                      wRy(:,ii,jj,kk) = subuh0(:,ii,jj,kk) + 0.5*slopey(:,ii,jj,kk) 
                      CALL PDEFlux(FLy,wLy(:,ii,jj,kk),subpar(:,ii,jj,kk))
                      CALL PDEFlux(FRy,wRy(:,ii,jj,kk),subpar(:,ii,jj,kk))
                      Qt = Qt - (FRy(:,2)-FLy(:,2))/ldx(2) 
                  ENDIF 
                  IF(nDim>=3) THEN
                      wLz(:,ii,jj,kk) = subuh0(:,ii,jj,kk) - 0.5*slopez(:,ii,jj,kk) 
                      wRz(:,ii,jj,kk) = subuh0(:,ii,jj,kk) + 0.5*slopez(:,ii,jj,kk) 
                      CALL PDEFlux(FLz,wLz(:,ii,jj,kk),subpar(:,ii,jj,kk))
                      CALL PDEFlux(FRz,wRz(:,ii,jj,kk),subpar(:,ii,jj,kk))
                      Qt = Qt - (FRz(:,3)-FLz(:,3))/ldx(3) 
                  ENDIF 
                  wLx(:,ii,jj,kk) = wLx(:,ii,jj,kk) + 0.5*dt*Qt 
                  wRx(:,ii,jj,kk) = wRx(:,ii,jj,kk) + 0.5*dt*Qt 
                  wLy(:,ii,jj,kk) = wLy(:,ii,jj,kk) + 0.5*dt*Qt 
                  wRy(:,ii,jj,kk) = wRy(:,ii,jj,kk) + 0.5*dt*Qt 
                  wLz(:,ii,jj,kk) = wLz(:,ii,jj,kk) + 0.5*dt*Qt 
                  wRz(:,ii,jj,kk) = wRz(:,ii,jj,kk) + 0.5*dt*Qt 
              ENDDO
             ENDDO
           ENDDO           
           
           nv = nvi(:,1) 
           ! Solve the Riemann problems on the x-edges of the subgrid 
           DO kk = 1, nSubLimV(3)
            DO jj = 1, nSubLimV(2)                 
             DO ii = 1, nSubLimV(1)+dn(1)
                 CALL RusanovFluxS(Fx(:,ii,jj,kk),wRx(:,ii-1,jj,kk),wLx(:,ii,jj,kk),subpar(:,ii-1,jj,kk),subpar(:,ii,jj,kk),nv) 
             ENDDO
            ENDDO
           ENDDO    
           IF(nDim>=2) THEN
               nv = nvi(:,2) 
               ! Solve the Riemann problems on the y-edges of the subgrid 
               DO kk = 1, nSubLimV(3)  
                DO jj = 1, nSubLimV(2)+dn(2)                 
                 DO ii = 1, nSubLimV(1)
                     CALL RusanovFluxS(Fy(:,ii,jj,kk),wRy(:,ii,jj-1,kk),wLy(:,ii,jj,kk),subpar(:,ii,jj-1,kk),subpar(:,ii,jj,kk),nv) 
                 ENDDO
                ENDDO
               ENDDO 
           ENDIF           
           IF(nDim>=3) THEN
               nv = nvi(:,3) 
               ! Solve the Riemann problems on the z-edges of the subgrid 
               DO kk = 1, nSubLimV(3)+dn(3) 
                DO jj = 1, nSubLimV(2)
                 DO ii = 1, nSubLimV(1)
                     CALL RusanovFluxS(Fz(:,ii,jj,kk),wRz(:,ii,jj,kk-1),wLz(:,ii,jj,kk),subpar(:,ii,jj,kk-1),subpar(:,ii,jj,kk),nv)
                 ENDDO
                ENDDO
               ENDDO 
           ENDIF           
           
           ! Evolve the subcell data inside cell i using a simple first order finite volume scheme 
           DO kk = 1, nSubLimV(3)
            DO jj = 1, nSubLimV(2) 
             DO ii = 1, nSubLimV(1)
                 subuh1(:,ii,jj,kk) = subuh0(:,ii,jj,kk) - dt/ldx(1)*(Fx(:,ii+1,jj,kk)-Fx(:,ii,jj,kk)) - dt/ldx(2)*(Fy(:,ii,jj+1,kk)-Fy(:,ii,jj,kk)) - dt/ldx(3)*(Fz(:,ii,jj,kk+1)-Fz(:,ii,jj,kk)) 
             ENDDO
            ENDDO
           ENDDO 
           !
           IF(recompute(i).EQ.1) THEN
                ! If a cell is really troubled, then set the limiter status to 1 and set the new subcell data 
                Limiter(i)%status = 1 
                ! Save the evolved data in the limiter
                DO kk = 1, nSubLimV(3) 
                 DO jj = 1, nSubLimV(2) 
                  DO ii = 1, nSubLimV(1)  
                      Limiter(i)%NewLh(:,ii,jj,kk) = subuh1(:,ii,jj,kk)  
                  ENDDO
                 ENDDO
                ENDDO    
                ! Put the evolved cell-averaged data back into the DG polynomial 
                CALL PutSubCellData(uh(:,:,:,:,i),Limiter(i)%NewLh) 
           ELSEIF(recompute(i).EQ.2) THEN
                ! Put the evolved cell-averaged data back into the DG polynomial 
                CALL PutSubCellData(uh(:,:,:,:,i),subuh1) 
                ! If we recompute neighbors of troubled cells, and the resulting reconstructed DG polynomial is troubled, then activate the limiter also in this case 
                CALL DMP(dmpresult,uh(:,:,:,:,i),Limiter(i),0.0)  
                IF(dmpresult) THEN 
                    Limiter(i)%status = 0 
                ELSE
                    Limiter(i)%status = 1  
                    IF(Limiter(i)%oldstatus.EQ.0) THEN 
                       ALLOCATE( Limiter(i)%Lh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))    ) 
                       ALLOCATE( Limiter(i)%NewLh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) ) 
                    ENDIF
                    ! Save the evolved data in the limiter
                    DO kk = 1, nSubLimV(3) 
                     DO jj = 1, nSubLimV(2) 
                      DO ii = 1, nSubLimV(1)  
                          Limiter(i)%NewLh(:,ii,jj,kk) = subuh1(:,ii,jj,kk)  
                      ENDDO
                     ENDDO
                    ENDDO    
                ENDIF            
            ENDIF 
           !
       ENDIF
       !
   ENDDO       
   !
   CONTINUE
   !
END SUBROUTINE Subcellrecompute        
    
SUBROUTINE RusanovFluxS(flux,QL,QR,parL,parR,nv)
    USE typesDef
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: QL(nVar), QR(nVar), parL(nParam), parR(nParam), nv(d) 
    REAL, INTENT(OUT) :: flux(nVar)
    ! Local variables
    REAL :: smax, LL(nVar), LR(nVar) 
    REAL :: FL(nVar,d), FR(nVar,d) 
    ! 
    CALL PDEEigenvalues(LL,QL,parL,nv) 
    CALL PDEEigenvalues(LR,QR,parR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) ) 
    CALL PDEFlux(FL,QL,parL) 
    CALL PDEFlux(FR,QR,parR)     
    flux = 0.5*MATMUL(FL+FR,nv) - 0.5*smax*(QR-QL)   
    ! 
END SUBROUTINE    
    
SUBROUTINE minmod(c,a,b,nnn)
    Implicit none
    ! Argument list 
    INTEGER, INTENT(IN)   :: nnn
    REAL, INTENT(IN)      :: a(nnn), b(nnn) 
    REAL, INTENT(OUT)     :: c(nnn)
    ! Local variables 
    INTEGER :: i 
    !
    DO i = 1, nnn 
        IF(a(i)*b(i).LE.0.) THEN
            c(i) = 0. 
        ELSE
            IF( ABS(a(i)).LT.ABS(b(i)) ) THEN
                c(i) = a(i) 
            ELSE
                c(i) = b(i)
            ENDIF 
        ENDIF 
    ENDDO 
    ! 
END SUBROUTINE minmod 
    
SUBROUTINE GetSubcell_uh(outuh,i) 
   USE typesDef
   IMPLICIT NONE 
   ! Argument list
   REAL    :: outuh(nVar,1:nSubLimV(1)+1,1:nSubLimV(2)+1,1:nSubLimV(3)+1) 
   INTEGER :: i 
   ! Local variables 
   INTEGER :: k, j, kj, ii, jj, kk, iVar, reflev, iDim, nv(d), iii, jjj, kkk, lll, pp, qq, rr, aa, bb, cc, iErr   
   REAL    :: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
   INTEGER :: NodeCounter(1:nSubLimV(1)+1,1:nSubLimV(2)+1,1:nSubLimV(3)+1) 
   REAL    :: subuh(nVar,(1-nSubLimV(1)):2*nSubLimV(1),(1-nSubLimV(2)):2*nSubLimV(2),(1-nSubLimV(3)):2*nSubLimV(3)) 
   !
   ! Get the initial data, either from a previous limiter stage, or from healthy uh initial data 
   DO kk = -dn(3), dn(3) 
    DO jj = -dn(2), dn(2) 
     DO ii = -dn(1), dn(1) 
        j = neighbor(ii,jj,kk,i)
        IF(Limiter(j)%oldstatus.EQ.0) THEN
            CALL GetSubcellData(lim,olduh(:,:,:,:,j)) 
        ELSE
            lim = Limiter(j)%Lh 
        ENDIF
        subuh(:,1+ii*nSubLimV(1):(ii+1)*nSublimV(1),1+jj*nSubLimV(2):(jj+1)*nSublimV(2),1+kk*nSubLimV(3):(kk+1)*nSublimV(3)) = lim(:,:,:,:)    
     ENDDO
    ENDDO
   ENDDO
   !
   NodeCounter = 0 
   outuh = 0. 
   DO kk = 1, nSubLimV(3)
    DO jj = 1, nSubLimV(2)
     DO ii = 1, nSubLimV(1)
        DO aa = 1, nVtx 
          iii = ReferenceElement(1,aa) 
          jjj = ReferenceElement(2,aa) 
          kkk = ReferenceElement(3,aa) 
          outuh(:,ii+iii,jj+jjj,kk+kkk) = outuh(:,ii+iii,jj+jjj,kk+kkk) + subuh(:,ii,jj,kk) 
          NodeCounter(ii+iii,jj+jjj,kk+kkk) = NodeCounter(ii+iii,jj+jjj,kk+kkk) + 1 
        ENDDO 
     ENDDO
    ENDDO
   ENDDO
   !
   DO kk = 1, nSubLimV(3)+dn(3) 
    DO jj = 1, nSubLimV(2)+dn(2)
     DO ii = 1, nSubLimV(1)+dn(1)
         outuh(:,ii,jj,kk) = outuh(:,ii,jj,kk)/NodeCounter(ii,jj,kk) 
     ENDDO
    ENDDO
   ENDDO
   ! 
END SUBROUTINE GetSubcell_uh          