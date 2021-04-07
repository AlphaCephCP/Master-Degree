
    
  SUBROUTINE LocateSpherePickpoints    
   USE typesDef
   IMPLICIT NONE  
   ! Local variables 
   INTEGER  :: i, reflev, iElem, ii, jj, kk, oldlevel  
   REAL     :: lx0(d), ldx(d) 
   LOGICAL  :: found 
   !
   DO i = 1, nSpherePickPoint
     IF(SpherePickPoint(i)%iElem0.EQ.0) THEN
        ! SpherePickpoint is not in the current CPU, so skip it. 
        CYCLE
     ENDIF 
     iElem = SpherePickPoint(i)%iElem 
     iElem = SpherePickPoint(i)%iElem0
        lx0 = x(:,tri(1,iElem))
        SpherePickPoint(i)%iElem     = iElem 
        SpherePickPoint(i)%xi        = ( SpherePickPoint(i)%x-lx0 )/dx  
   ENDDO  
   !
END SUBROUTINE LocateSpherePickpoints 

    
SUBROUTINE DefineSpherePickpoints
    USE typesDef
    IMPLICIT NONE
    !INTEGER, PARAMETER :: nPhi=100,  nTheta=100
    REAL :: theta(N+1,SphPick_nTheta), phi(N+1,SphPick_nPhi) 
    !REAL :: Obs(nDOF,nDOF,nTheta,nPhi)
    REAL :: x1,y1,z1,dtheta,dphi,xiGP(N+1),wGP(N+1),Qloc(nVar)
    INTEGER :: nElemB, count,i,j,ii,jj,iii,jjj,kkk
    REAL :: Int, R
    !nElemB = nPhi*nTheta  ! 2d mesh 
    !nPhi,nTheta,SphereRadius
    R=SphereRadius
    dPhi = 2*EQN%PI/SphPick_nPhi
    dTheta = EQN%PI/SphPick_nTheta
    !  
    CALL gauleg(0.,1.,xiGP,wGP,N+1)
    !
    !phi domain
    DO i=1,SphPick_nPhi
        phi(:,i) = dPhi*real(i-1)+xiGP(:)*dPhi
    ENDDO
    !
    !theta domain
    DO i=1,SphPick_nTheta
        theta(:,i) = dTheta*real(i-1)+xiGP(:)*dTheta
    ENDDO
    !
    !dx(1) = (xR(1)-xL(1))/IMAX
    !dx(2) = (xR(2)-xL(2))/JMAX 
    !dx(3) = (xR(3)-xL(3))/KMAX 
    !  
    count = 0
    DO i=1,SphPick_nTheta
        DO j=1,SphPick_nPhi
            DO ii=1,N+1
                DO jj=1,N+1
                        count = count + 1
                        ! Cartesian coordinates 
                        SpherePickPoint(count)%iTheta = ii
                        SpherePickPoint(count)%iPhi = jj
                        SpherePickPoint(count)%x(1) = R * dsin(theta(ii,i))*dcos(phi(jj,j))
                        SpherePickPoint(count)%x(2) = R * dsin(theta(ii,i))*dsin(phi(jj,j))
                        SpherePickPoint(count)%x(3) = R * dcos(theta(ii,i))
                        ! 
                        SpherePickPoint(count)%r(1) = R
                        SpherePickPoint(count)%r(2) = theta(ii,i)
                        SpherePickPoint(count)%r(3) = phi(jj,j)
                        !
                        !! Locate GL point on the main grid (reflev. zero)
                        !iii = MIN( IMAX, MAX( 1, CEILING( (x1-xL(1))/dx(1) ) ) ) 
                        !jjj = MIN( JMAX, MAX( 1, CEILING( (x1-xL(2))/dx(2) ) ) ) 
                        !kkk = MIN( KMAX, MAX( 1, CEILING( (z1-xL(3))/dx(3) ) ) )  
                        !!
                        !SpherePickPoint(count)%iElem0    = idxe0(iii,jjj,kkk) 
                        !SpherePickPoint(count)%iElem     = idxe0(iii,jjj,kkk)  !tmp
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    !
    !
    CONTINUE
    !
END SUBROUTINE DefineSpherePickpoints
   
    
    
    
SUBROUTINE ComputeIntegralSphere(Int_h) !,Obervable)
    USE typesDef
    IMPLICIT NONE
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif
    ! input variables       
    !INTERFACE
    !  SUBROUTINE Obervable(O,V,x)  
    !     REAL, INTENT(IN) :: V,x
    !     REAL    :: O 
    !  END SUBROUTINE Obervable 
    !END INTERFACE       
    !INTERFACE
    !  SUBROUTINE ObervableVec(O,V,x)  
    !     USE MainVariables, ONLY : d
    !     REAL, INTENT(IN) :: V,x
    !     REAL    :: O(d) 
    !  END SUBROUTINE ObervableVec 
    !END INTERFACE 
    ! local variables
    !INTEGER, PARAMETER :: nPhi=100,  nTheta=100
    REAL :: theta(N+1,SphPick_nTheta), phi(N+1,SphPick_nPhi)
    REAL :: Obs(nObs) !,ObsVec(d)
    REAL :: x1,y1,z1,dtheta,dphi,xiGP(N+1),wGP(N+1),Qloc(nVar),Vloc(nVar) 
    INTEGER :: count,ii,jj,i,iElem
    REAL :: Int_h(nObs),mpiInt(nObs),rloc(d),xloc(d),R,dS,norm(d) !IntVec(d),mpiIntVec(d),
    !  
    !nElemB = nPhi*nTheta  ! 2d mesh 
    !nPhi,nTheta,SphereRadius
    R=SphereRadius
    dPhi = 2.0*EQN%PI/SphPick_nPhi
    dTheta = EQN%PI/SphPick_nTheta
    !
    CALL gauleg(0.,1.,xiGP,wGP,N+1)
    !
    ! evaluate the observable degrees of freedom,
    ! and compute the integral
    !
    IF(R.GT.MINVAL(ABS((/xL, xR/)))) THEN
        PRINT *, 'The selected sphere is outside the chosen spatial domain!'
        WRITE(*,*) 'Radius',R
        WRITE(*,*) xL
        WRITE(*,*) xR
        STOP
    ENDIF
    !
    !
    Int_h = 0.
    !IntVec = 0.
    count = 0
    Obs = 0.
    !
    DO i=1,nSpherePickPoint
        !
        IF(SpherePickPoint(i)%iElem0.EQ.0) THEN
            continue
            CYCLE
        ENDIF 
        !
        ii = SpherePickPoint(i)%iTheta
        jj = SpherePickPoint(i)%iPhi
        iElem = SpherePickPoint(i)%iElem
        ! 
        ! Cartesian coordinates
        xloc = SpherePickPoint(i)%x
        !
        rloc = SpherePickPoint(i)%r
        !      
        ! Evaluate conserved variables at pickpoints:
        CALL DoPick_uh(Qloc,Vloc,iElem,SpherePickPoint(i)%xi)
        ! 
        ! Evaluate primitive variables
        !CALL Cons2Prim(Vloc,Qloc,....)
        !
        ! Evaluate Observable at the Gauss-Legendre q. points of the spherical surface (i.e. at Cartesian pickpoints)
        CALL Observable(Obs,Vloc,xloc,rloc)
        !Obs = 1.0
        !
        ! surface element (may be it is different in curved space)
        dS = R**2*dsin(rloc(2))
        norm = (/ dsin(rloc(2))*dcos(rloc(3)), dsin(rloc(2))*dsin(rloc(3)), dcos(rloc(2)) /)
        !
        Int_h = Int_h + wGP(jj)*wGP(ii)*Obs*dS  !*DOT_PRODUCT(norm,(/1.0, 0. , 0. /) )
        !IntVec = IntVec + wGP(jj)*wGP(ii)*Obs*dS*norm
        ! 
    ENDDO
    !
    ! constant factor for uniform grids
    !
#ifdef PARALLEL
  CALL MPI_ALLREDUCE(Int_h,mpiInt,nObs,MPI_AUTO_REAL,MPI_SUM,MPI_COMM_WORLD,mpiErr) 
  Int_h = mpiInt*dtheta*dphi  
  !CALL MPI_ALLREDUCE(IntVec, mpiIntVec, d,MPI_AUTO_REAL,MPI_SUM,MPI_COMM_WORLD,mpiErr)
  !IntVec  = mpiIntVec*dtheta*dphi 
#else
    Int_h = Int_h*dtheta*dphi 
    !IntVec = IntVec*dtheta*dphi 
#endif 
    
    !
    !
    continue
    !
    END SUBROUTINE ComputeIntegralSphere
    
  
SUBROUTINE Observable(Obs,Vloc,xloc,rloc)
    USE typesDef, ONLY : nVar,d,nObs
    IMPLICIT NONE
    ! input variable
    REAL :: Obs(nObs),Vloc(nVar),xloc(d),rloc(d)
    ! local variables
    REAL :: g_cov(3,3),g_contr(3,3),b2,v2,gp,v_cov(3),v_contr(3),b_cov(3),b_contr(3),norm(d)
    !
    !put here the definition of your observable O = O(V,x)
    !
    ! e.g. Compute the area
    Obs(1) = 1.0
    !
    ! e.g. Compute the flux of a constant vecotr field (should be zero) 
    v_cov(1:3) = (/ 0.1, 0.1, -0.6 /) 
    norm(1:3) = (/ dsin(rloc(2))*dcos(rloc(3)), dsin(rloc(2))*dsin(rloc(3)), dcos(rloc(2)) /) !this is wrong in curved space
    !
    v2 = DOT_PRODUCT(v_cov, norm)  ! becareful here: choose the proper norm vector (contravariant)
    Obs(2) = v2  
    ! 
    ! e.g. Compute the mean vector of the norm... i.e. zero   
    Obs(3:5) = norm
    ! 
    ! e.g. Compute v**2 or b**2
    !g_cov(1,1) = V(14)
    !g_cov(1,2) = V(15)
    !g_cov(1,3) = V(16)
    !g_cov(2,2) = V(17)
    !g_cov(2,3) = V(18)
    !g_cov(3,3) = V(19)
    !g_cov(2,1) = V(15)
    !g_cov(3,1) = V(16)
    !g_cov(3,2) = V(17)
    !!
    !CALL MatrixInverse3x3(g_cov,g_contr,gp)
    !!
    !v_cov(1:3) = V(2:4)
    !b_contr(1:3) = V(6:8)
    !!
    !v_contr = MATMUL(g_contr,v_cov)
    !b_cov = MATMUL(g_cov,b_contr)
    !!
    !v2 = DOT_PRODUCT(v_cov,v_contr)
    !b2 = DOT_PRODUCT(b_cov,b_contr)
    !!
    !Obs(2) = v2
    !Obs(3) = b2
    !
    ! e.g. Compute the mean velocity vector field  
    !v_cov(1:3) = V(2:4)  
    !Obs(1:3) = v_cov(1:3)   
    ! 
    !
    ! e.g. Compute the flux of the velocity 
    !v_cov(1:3) = V(2:4)  
    !norm(1:3) = (/ dsin(rloc(2))*dcos(rloc(3)), dsin(rloc(2))*dsin(rloc(3)), dcos(rloc(2)) /) !this is wrong in curved space
    !
    !v2 = DOT_PRODUCT(v_cov, norm)  ! becareful here: choose the proper norm vector (contravariant)
    !Obs(1) = v2  
    ! 
    continue
    !
END SUBROUTINE Observable
    
SUBROUTINE DoPick_uh(Q,V,iElem,xi) 
    !--------------------------------------------------------------------------
    USE typesDef
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    INTEGER                        :: iElem  
    REAL                           :: Q(nVar), V(nVar)
    REAL                           :: xi(d) 
    ! local Variables
    INTEGER                        :: ii,jj,kk,ll,c,iErr   
    !REAL                           :: lqh(nVar,M+1,M+1,M+1,0:M+1) 
    REAL                           :: theta, aux(d)  
    REAL                           :: psii(N+1), psij(N+1), psik(N+1), psi_xi(N+1) 
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: iElem,xi  
    INTENT(OUT)                    :: Q,V
    !--------------------------------------------------------------------------
    !
    CALL BaseFunc1D(psii,psi_xi,xi(1))  
    CALL BaseFunc1D(psij,psi_xi,xi(2))  
    CALL BaseFunc1D(psik,psi_xi,xi(3))  
    Q = 0.
    c = 0 
     DO kk = 1, N+1
      DO jj = 1, N+1
       DO ii = 1, N+1
        c = c + 1
        aux = (/ psii(ii), psij(jj), psik(kk) /) 
        theta = PRODUCT(aux(1:nDim))
        Q = Q + theta*uh(:,ii,jj,kk,iElem)
       ENDDO
      ENDDO
     ENDDO
    ! 
    !Q = MATMUL( theta, TRANSPOSE(uh(:,:,iElem)) )  
    CALL PDECons2Prim(V,Q,iErr)
    ! 
END SUBROUTINE DoPick_uh

