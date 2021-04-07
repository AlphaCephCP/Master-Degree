SUBROUTINE ADERDGInit
    USE typesDef
    IMPLICIT NONE
    ! Local variables 
    INTEGER          :: i, j, k, ii, jj, kk, l, c, iGP, iElem, iVar, VMAX(d), count, cnt, iDim, idxs(d), iper(d) 
    REAL             :: phi(N+1), phi0(N+1), phi1(N+1), phi_xi(N+1)
    REAL             :: phi_i(N+1), phi_j(N+1), phi_k(N+1) 
    REAL             :: u0(nVar), xGP(d), x0(d), subxi(N+1), aux(d), nv(d), par0(nParam)  
    REAL             :: lparambnd(nParam,6,nDOF(2),nDOF(3))    
    REAL             :: xi, dxi, xi1, xi2 
    REAL             :: TestMatrix(N+1,N+1), TestMatrix2(nSubLim,nSubLim), test(nSubLim)
    REAL, POINTER    :: LSQM(:,:), iLSQM(:,:), LSQrhs(:,:) 
    INTEGER, POINTER :: idxn(:,:,:), idxe(:,:,:)  
    LOGICAL          :: dmpresult 
    ! 
    ! ----------------- Some important preliminary stuff. Do not touch! ------------------- 
    dn(:) = 0 
    DO i = 1, nDim
        dn(i) = 1
    ENDDO 
    ! According to the number of space dimensions, we set the number of degrees of freedom 
    ! i = 0 is the time dimension 
    nDOF(:) = 1 
    DO i = 0, nDim
        nDOF(i) = N+1
    ENDDO 
    ! ------------------------------------------------------------------------------------- 
    !
    ! Some info about the PDE system 
    !
    EQN%Pi    = ACOS(-1.0)                                          ! Pi 
    EQN%gamma = 1.4                                                 ! ratio of specific heats for compressible Euler 
    !
    ! Here, you need to define the computational domain and the number of cells. 
    ! This is typically read from a parameter file 
    !
    xL = (/ -1.0, -1.0, -1.0 /)                                     ! lower-left corner of the domain 
    xR = (/ +1.0, +1.0, +1.0 /)                                     ! upper right corner of the domain 
    !xL = (/  0.0,  0.0, -1.0 /)                                     ! lower-left corner of the domain 
    !xR = (/ 10.0, 10.0, +1.0 /)                                     ! upper right corner of the domain 
    IMAX = 50                                                       ! Number of elements in x,y,z direction 
    JMAX = 50  
    KMAX = 1  
    VMAX = (/ IMAX, JMAX, KMAX /)                                   ! Vector of the number of elements in each space dimension 
    dx = (xR-xL)/VMAX                                               ! Mesh spacing 
    NMAX = 100000                                                   ! Max. number of time steps 
    timestep = 0                                                    ! initial time step number 
    time = 0.                                                       ! initial time 
    tend = 1.0                                                      ! final time 
    !ICType = 'Sod'                                                  ! type of initial condition 'ShuVortex2D'      
    Basefile = 'Test'                                                ! Base filename for writing results 
    !tend = 1.0                                                     ! final time 
    !ICType = 'ShuVortex2D'                                         ! type of initial condition 'ShuVortex2D'      
    !Basefile = 'ShuVortex2D'                                       ! Base filename for writing results 
    !Periodic(:) = (/ .FALSE., .TRUE., .FALSE. /)                     ! periodic BC 
    !AnalyseType = 0                                                 ! comparison with exact solution
    !
    nElem = IMAX*JMAX*KMAX                                          ! Number of elements 
    ALLOCATE(  uh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the spatial degrees of freedom 
    ALLOCATE( duh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the degrees of freedom for the update 
    ALLOCATE( qhi(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )         ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( Fhi(nVar, d, nDOF(1), nDOF(2), nDOF(3), nElem) )      ! Allocate the time-averaged space-time degrees of freedom 
    ALLOCATE( parh(nParam, nDOF(1), nDOF(2), nDOF(3), nElem) )      ! Allocate the degrees of freedom for the material parameters 
    ALLOCATE( qbnd(nVar, 6, nDOF(2), nDOF(3), nElem) )              ! Allocate the time-averaged boundary-extrapolated values for Q 
    ALLOCATE( parbnd(nParam, 6, nDOF(2), nDOF(3), nElem) )          ! Allocate the boundary-extrapolated material parameters 
    ALLOCATE( Fbnd(nVar, 6, nDOF(2), nDOF(3), nElem) )              ! Allocate the time-averaged boundary-extrapolated values for the normal flux F * n 
    nNode = (IMAX+1)*(JMAX+1)*(KMAX+1)                              ! number of nodes 
    ALLOCATE( x(d, nNode) )                                         ! Allocate the nodes 
    ALLOCATE( idxn(IMAX+dn(1),JMAX+dn(2),KMAX+dn(3))  )                                              
    ! Define the node coordinates and the node numbers                                 
    count = 0 
    DO k = 1, KMAX+dn(3)  
        DO j = 1, JMAX+dn(2) 
            DO i = 1, IMAX+dn(1) 
                count = count + 1 
                x(:,count) = xL(:) + (/ i-1, j-1, k-1/)*dx(:) 
                idxn(i,j,k) = count 
            ENDDO
        ENDDO
    ENDDO
    ! define the connectivity between the elements and the nodes. You can do this more compactly via a loop, but I prefer the select case, which makes the explicit 
    ! construction of each element clearer 
    ALLOCATE( tri(nVtx, nElem)     ) 
    ALLOCATE( idxe(IMAX,JMAX,KMAX) ) 
    count = 0  
    DO k = 1, KMAX  
        DO j = 1, JMAX
            DO i = 1, IMAX 
                count = count + 1
                idxe(i,j,k) = count 
                SELECT CASE(nDim)
                CASE(1)                    
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k) /) 
                CASE(2)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k) /) 
                CASE(3)
                    tri(:,count) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i,j+1,k), idxn(i+1,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i,j+1,k+1), idxn(i+1,j+1,k+1) /) 
                END SELECT                
            ENDDO
        ENDDO
    ENDDO
    ! compute the Voronoi neighbors of a cell
    ALLOCATE( neighbor(-1:1,-1:1,-1:1,nElem) ) 
    DO i = 1, nElem
        neighbor(:,:,:,i) = i  
    ENDDO    
    DO k = 1, KMAX  
     DO j = 1, JMAX
      DO i = 1, IMAX 
        count = idxe(i,j,k) 
        ! 
        DO ii = -dn(1), dn(1) 
         DO jj = -dn(2), dn(2) 
          DO kk = -dn(3), dn(3)
              IF( (MIN( i+ii, j+jj, k+kk ).GE.1) .AND. (MAX( i+ii-IMAX, j+jj-JMAX, k+kk-KMAX ).LE.0) ) THEN 
                neighbor(ii,jj,kk,count) = idxe(i+ii,j+jj,k+kk) 
              ENDIF
              idxs = (/ i, j, k /) + (/ ii, jj, kk /) 
              DO iDim = 1, d 
                IF(periodic(iDim)) THEN
                  IF(idxs(iDim).LT.1) THEN
                    idxs(iDim) = MIN(VMAX(iDim),MAX(1,VMAX(iDim) + idxs(iDim))) 
                  ENDIF 
                  IF(idxs(iDim).GT.VMAX(iDim)) THEN
                    idxs(iDim) = MIN(VMAX(iDim),MAX(1,idxs(iDim) - VMAX(iDim))) 
                  ENDIF 
                ENDIF 
                IF(ANY(idxs.LT.1) ) THEN                     
                    CYCLE
                ENDIF
                IF(ANY((idxs-VMAX).GT.0)) THEN 
                    CYCLE 
                ENDIF
                IF(idxe(idxs(1),idxs(2),idxs(3)).NE.0) THEN
                    neighbor(ii,jj,kk,count) = idxe(idxs(1),idxs(2),idxs(3))  
                ELSE
                    neighbor(ii,jj,kk,count) = count
                ENDIF
               ENDDO
               !
          ENDDO
         ENDDO
        ENDDO
        !  
      ENDDO
     ENDDO     
    ENDDO
    !
    ! Compute the mapping from face to 3D neighbor index 
    !
    Face2Neigh(:,1) = (/ -1, 0, 0 /)   
    Face2Neigh(:,2) = (/ +1, 0, 0 /)   
    Face2Neigh(:,3) = (/  0,-1, 0 /)   
    Face2Neigh(:,4) = (/  0,+1, 0 /)   
    Face2Neigh(:,5) = (/  0, 0,-1 /)   
    Face2Neigh(:,6) = (/  0, 0,+1 /)   
    ! 
    DEALLOCATE( idxn )
    ! define the connectivity between the faces and the elements 
    ! count how many faces we have 
    DO iDim = 1, d
        IF(Periodic(iDim)) THEN
            iper(iDim) = 0
        ELSE
            iper(iDim) = 1 
        ENDIF 
    ENDDO 
    nFace = KMAX*JMAX*(IMAX+iper(1))
    IF(nDim>=2) THEN 
        nFace = nFace + KMAX*(JMAX+iper(2))*IMAX
    ENDIF
    IF(nDim>=3) THEN 
        nFace = nFace + (KMAX+iper(3))*JMAX*IMAX 
    ENDIF
    ALLOCATE( Face(nFace) ) 
    ! x faces 
    count = 0 
    DO k = 1, KMAX
     DO j = 1, JMAX
      DO i = 1, IMAX+iper(1)  
          count = count + 1 
          IF(i.EQ.1) THEN
            IF(Periodic(1)) THEN
                Face(count)%Left  = idxe(IMAX,j,k) 
                Face(count)%Right = idxe(i,j,k) 
                Face(count)%qL => qBnd(:,2,:,:,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,1,:,:,Face(count)%Right) 
                Face(count)%FL => FBnd(:,2,:,:,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,1,:,:,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,2,:,:,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,1,:,:,Face(count)%Right) 
            ELSE
                Face(count)%Left  = 0 
                Face(count)%Right = idxe(i,j,k)
                ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%qR => qBnd(:,1,:,:,Face(count)%Right) 
                Face(count)%FR => FBnd(:,1,:,:,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,1,:,:,Face(count)%Right) 
                Face(count)%paramR => parbnd(:,1,:,:,Face(count)%Right) 
            ENDIF 
          ELSEIF(i.EQ.IMAX+1) THEN 
                Face(count)%Left  = idxe(i-1,j,k) 
                Face(count)%Right = 0 
                Face(count)%qL => qBnd(:,2,:,:,Face(count)%Left ) 
                Face(count)%FL => FBnd(:,2,:,:,Face(count)%Left ) 
                ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                Face(count)%paramL => parbnd(:,2,:,:,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,2,:,:,Face(count)%Left ) 
          ELSE              
            Face(count)%Left  = idxe(i-1,j,k) 
            Face(count)%Right = idxe(i,j,k) 
            Face(count)%qL => qBnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%qR => qBnd(:,1,:,:,Face(count)%Right) 
            Face(count)%FL => FBnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%FR => FBnd(:,1,:,:,Face(count)%Right) 
            Face(count)%paramL => parbnd(:,2,:,:,Face(count)%Left ) 
            Face(count)%paramR => parbnd(:,1,:,:,Face(count)%Right) 
          ENDIF 
          Face(count)%nv = (/ 1., 0., 0. /) ! set face normal vector 
      ENDDO
     ENDDO
    ENDDO 
    ! y faces 
    IF(nDim>=2) THEN
        DO k = 1, KMAX
         DO j = 1, JMAX+iper(2) 
          DO i = 1, IMAX  
              count = count + 1 
              IF(j.EQ.1) THEN
                IF(Periodic(2)) THEN
                    Face(count)%Left  = idxe(i,JMAX,k)
                    Face(count)%Right = idxe(i,j,k)
                    Face(count)%qL => qBnd(:,4,:,:,Face(count)%Left ) 
                    Face(count)%qR => qBnd(:,3,:,:,Face(count)%Right) 
                    Face(count)%FL => FBnd(:,4,:,:,Face(count)%Left ) 
                    Face(count)%FR => FBnd(:,3,:,:,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,4,:,:,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,3,:,:,Face(count)%Right) 
                ELSE
                    Face(count)%Left  = 0
                    Face(count)%Right = idxe(i,j,k)
                    ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%qR => qBnd(:,3,:,:,Face(count)%Right) 
                    Face(count)%FR => FBnd(:,3,:,:,Face(count)%Right)
                    Face(count)%paramL => parbnd(:,3,:,:,Face(count)%Right) 
                    Face(count)%paramR => parbnd(:,3,:,:,Face(count)%Right) 
                ENDIF
              ELSEIF(j.EQ.JMAX+1) THEN 
                    Face(count)%Left  = idxe(i,j-1,k) 
                    Face(count)%Right = 0 
                    Face(count)%qL => qBnd(:,4,:,:,Face(count)%Left ) 
                    Face(count)%FL => FBnd(:,4,:,:,Face(count)%Left ) 
                    ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%paramL => parbnd(:,4,:,:,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,4,:,:,Face(count)%Left ) 
              ELSE              
                Face(count)%Left  = idxe(i,j-1,k) 
                Face(count)%Right = idxe(i,j,k) 
                Face(count)%qL => qBnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,3,:,:,Face(count)%Right) 
                Face(count)%FL => FBnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,3,:,:,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,4,:,:,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,3,:,:,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 1., 0. /) ! set face normal vector 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    ! z faces 
    IF(nDim>=3) THEN
        DO k = 1, KMAX+iper(3) 
         DO j = 1, JMAX 
          DO i = 1, IMAX  
              count = count + 1 
              IF(k.EQ.1) THEN
                IF(Periodic(3)) THEN
                    Face(count)%Left  = idxe(i,j,KMAX) 
                    Face(count)%Right = idxe(i,j,k) 
                    Face(count)%qL => qBnd(:,6,:,:,Face(count)%Left ) 
                    Face(count)%qR => qBnd(:,5,:,:,Face(count)%Right) 
                    Face(count)%FL => FBnd(:,6,:,:,Face(count)%Left ) 
                    Face(count)%FR => FBnd(:,5,:,:,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,6,:,:,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,5,:,:,Face(count)%Right) 
                ELSE
                    Face(count)%Left  = 0 
                    Face(count)%Right = idxe(i,j,k)
                    ALLOCATE( Face(count)%qL(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FL(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%qR => qBnd(:,5,:,:,Face(count)%Right) 
                    Face(count)%FR => FBnd(:,5,:,:,Face(count)%Right) 
                    Face(count)%paramL => parbnd(:,5,:,:,Face(count)%Right) 
                    Face(count)%paramR => parbnd(:,5,:,:,Face(count)%Right) 
                ENDIF 
              ELSEIF(k.EQ.KMAX+1) THEN 
                    Face(count)%Left  = idxe(i,j,k-1) 
                    Face(count)%Right = 0 
                    Face(count)%qL => qBnd(:,6,:,:,Face(count)%Left ) 
                    Face(count)%FL => FBnd(:,6,:,:,Face(count)%Left ) 
                    ALLOCATE( Face(count)%qR(nVar,nDOF(2),nDOF(3)) ) 
                    ALLOCATE( Face(count)%FR(nVar,nDOF(2),nDOF(3)) ) 
                    Face(count)%paramL => parbnd(:,6,:,:,Face(count)%Left ) 
                    Face(count)%paramR => parbnd(:,6,:,:,Face(count)%Left ) 
              ELSE              
                Face(count)%Left  = idxe(i,j,k-1) 
                Face(count)%Right = idxe(i,j,k) 
                Face(count)%qL => qBnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%qR => qBnd(:,5,:,:,Face(count)%Right) 
                Face(count)%FL => FBnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%FR => FBnd(:,5,:,:,Face(count)%Right) 
                Face(count)%paramL => parbnd(:,6,:,:,Face(count)%Left ) 
                Face(count)%paramR => parbnd(:,5,:,:,Face(count)%Right) 
              ENDIF 
              Face(count)%nv = (/ 0., 0., 1. /) ! set face normal vector 
          ENDDO
         ENDDO
        ENDDO 
    ENDIF 
    !
    ! We now need to define our basis functions. This is done by choosing a set of distinct 1D nodes in the interval [0,1] 
    ! The basis functions will be the Lagrange interpolation polynomials running through these nodes 
    CALL gauleg(0.,1.,xiGPN,wGPN,N+1)
    xin = xiGPN  ! WE make the following choice: the basis functions run through the Gauss-Legendre nodes (=> really orthogonal basis) 
    !
#ifdef LIMITER
    ! Okay this is not really related to limiter but nevertheless we don't support here (@SK)
    CALL InitPointSources(idxe) 
#endif
    DEALLOCATE( idxe ) 
    !
    ! Now, let us compute some of the important matrices in the ADER-DG method 
    ! 
    MM   = 0.   ! Element mass matrix 
    Kxi  = 0.   ! Element stiffness matrix 
    dudx = 0.   ! discrete derivative operator, which projects the derivatives onto the basis  
    DO iGP = 1, N+1 
        CALL BaseFunc1D(phi,phi_xi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! i) Mass-matrix 
                MM(k,l) = MM(k,l) + wGPN(iGP)*phi(k)*phi(l) 
                ! ii) Stiffness-matrix 
                Kxi(k,l) = Kxi(k,l) + wGPN(iGP)*phi_xi(k)*phi(l)  
            ENDDO
        ENDDO        
    ENDDO    
    CALL MatrixInverse(N+1,MM,iMM) 
    dudx = MATMUL( iMM, TRANSPOSE(Kxi) ) 
    CALL BaseFunc1D(phi0,phi_xi,0.0) ! Compute the basis functions on the left 
    CALL BaseFunc1D(phi1,phi_xi,1.0) ! Compute the basis function on the right 
    ! The flux matrices are all possible combinations of left and right 
    DO k = 1, N+1
        DO l = 1, N+1
            FLm(k,l) = phi0(k)*phi1(l)   ! Left contribution to the left flux matrix    (m = left  of the interface)  
            FLp(k,l) = phi0(k)*phi0(l)   ! Right contribution to the left flux matrix   (p = right of the interface) 
            FRm(k,l) = phi1(k)*phi1(l)   ! Left contribution to the right flux matrix   (m = left  of the interface) 
            FRp(k,l) = phi1(k)*phi0(l)   ! Right contribution to the right flux matrix  (p = right of the interface) 
        ENDDO
    ENDDO        
    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle) 
    F0 = phi0   ! upwinding in time = information comes from smaller times 
    F1 = FRm    ! upwinding in time = information comes from smaller times  
    K1 = F1 - Kxi 
    CALL MatrixInverse(N+1,K1,iK1)   
    FLcoeff = phi0  ! coefficients needed to extrapolate data onto the left  boundary 
    FRcoeff = phi1  ! coefficients needed to extrapolate data onto the right boundary 
    !
    ! For point sources, we also need the time mass matrix of the temporal Taylor series
    !
    MT   = 0.   ! Time mass matrix 
    DO iGP = 1, N+1 
        CALL TimeBaseFunc1D(phi,phi_xi,xiGPN(iGP))
        DO k = 1, N+1
            DO l = 1, N+1
                ! i) Mass-matrix 
                MT(k,l) = MT(k,l) + wGPN(iGP)*phi(k)*phi(l) 
            ENDDO
        ENDDO        
    ENDDO    
    CALL MatrixInverse(N+1,MT,iMT)     !
    ! For the fine output of each spatial degree of freedom onto a subgrid... 
    DO i = 1, N+1 
       subxi(i) = REAL(i-1)/REAL(N) 
    ENDDO
    cnt = 0 
    DO k = 1, N+1
     DO j = 1, N+1 
      DO i = 1, N+1  
         cnt = cnt + 1 
         CALL BaseFunc1D(phi_i,phi_xi,subxi(i))
         CALL BaseFunc1D(phi_j,phi_xi,subxi(j))
         CALL BaseFunc1D(phi_k,phi_xi,subxi(k))
         count = 0 
         DO kk = 1, nDOF(3) 
          DO jj = 1, nDOF(2)  
           DO ii = 1, nDOF(1) 
             count = count + 1 
             aux = (/ phi_i(ii), phi_j(jj), phi_k(kk) /) 
             SubOutputMatrix(count,cnt) = PRODUCT( aux(1:nDim) )                      
           ENDDO
          ENDDO
         ENDDO
         !
       ENDDO
      ENDDO
     ENDDO
     ! 
     ! subgrid triangulation for subgrid output 
     ! 
     ALLOCATE( idxn(N+1,N+1,N+1) )
     idxn = 0 
     c = 0 
     DO k = 1, N+1
        DO j = 1, N+1 
           DO i = 1, N+1 
              c = c + 1 
              allsubxi(:,c) = (/ REAL(i-1)/REAL(N), REAL(j-1)/REAL(N), REAL(k-1)/REAL(N) /)  
              idxn(i,j,k) = c
           ENDDO
        ENDDO
     ENDDO
     c = 0 
     DO k = 1, N 
        DO j = 1, N
           DO i = 1, N 
              c = c + 1 
              subtri(:,c) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i+1,j+1,k), idxn(i,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i+1,j+1,k+1), idxn(i,j+1,k+1) /)         
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE( idxn )     
     ! 
     ! subgrid triangulation for subgrid output of the limiter 
     ! 
     ALLOCATE( idxn(nSubLim+1,nSubLim+1,nSubLim+1) )
     idxn = 0 
     c = 0 
     DO k = 1, nSubLim+1
      DO j = 1, nSubLim+1 
        DO i = 1, nSubLim+1  
            c = c + 1 
            idxn(i,j,k) = c
            subxilim(:,c) = (/ (REAL(i)-1.)/REAL(nSubLim), (REAL(j)-1.)/REAL(nSubLim), (REAL(k)-1.)/REAL(nSubLim) /)  
        ENDDO
      ENDDO
     ENDDO
     c = 0 
     DO k = 1, nSubLim 
      DO j = 1, nSubLim
        DO i = 1, nSubLim  
            c = c + 1 
            subtrilim(:,c) = (/ idxn(i,j,k), idxn(i+1,j,k), idxn(i+1,j+1,k), idxn(i,j+1,k), idxn(i,j,k+1), idxn(i+1,j,k+1), idxn(i+1,j+1,k+1), idxn(i,j+1,k+1) /)         
        ENDDO
      ENDDO
     ENDDO
     !
    ! 
    ! Now prepare the subcell limiter 
    !
    ALLOCATE( olduh(nVar, nDOF(1), nDOF(2), nDOF(3), nElem) )       ! Allocate the old spatial degrees of freedom (needed for the limiter) 
    ALLOCATE( Limiter(nElem), recompute(nElem) ) 
    DO i = 1, nElem
        Limiter(i)%oldstatus = 0 
        Limiter(i)%status = 0
    ENDDO 
    nSubLimV = 1 
    DO i = 1, nDim
        nSubLimV(i) = nSubLim 
    ENDDO    
    !
    ALLOCATE( LSQM(N+2,N+2), iLSQM(N+2,N+2), LSQrhs(N+2,nSubLim) )  
    ALLOCATE( uh2lim(nSubLim,N+1), lim2uh(N+1,nSubLim)           )
    ALLOCATE( uh2lob(N+1,N+1)                                    )   
    uh2lim = 0. 
    dxi = 1./REAL(nSubLim)    ! the limiter uses nSubLim subintervals in each cell 
    DO i = 1, nSubLim
     xi1 = REAL(i-1)*dxi     ! left sub-interval boundary 
     xi2 = REAL(i  )*dxi     ! right sub-interval boundary 
     DO ii = 1, N+1  
        xi = xi1 + dxi*xiGPN(ii)     
        CALL BaseFunc1D(phi,phi_xi,xi) 
        uh2lim(i,:) = uh2lim(i,:) + wGPN(ii)*phi(:)
     ENDDO
    ENDDO
    ! To go back from the limited (piecewise constant sub-cell solution) to the uh solution
    ! we use a constrained least-squares reconstruction, which preserves the average exactly 
    LSQM(1:N+1,1:N+1) = 2*MATMUL( TRANSPOSE(uh2lim), uh2lim ) 
    LSQM(N+2,1:N+1)   =  wGPN(:) 
    LSQM(1:N+1,N+2)   = -wGPN(:) 
    LSQM(N+2,N+2) = 0. 
    CALL MatrixInverse(N+2,LSQM,iLSQM)
    LSQrhs(1:N+1,:)   = 2*TRANSPOSE(uh2lim)  
    LSQrhs(N+2,:)     = dxi    
    lim2uh = MATMUL( iLSQM(1:N+1,:), LSQrhs ) 
    ! 
    phi0 = 1. 
    test = MATMUL(uh2lim,phi0) 
    testmatrix = MATMUL(lim2uh,uh2lim) 
    testmatrix2 = MATMUL(uh2lim,lim2uh)     
    !
    ! Compute the Gauss-Lobatto quadrature nodes 
    CALL gaulob(0.,1.,xiLob,wLob,N+1)
    DO ii = 1, N+1
        CALL BaseFunc1D(phi,phi_xi,xiLob(ii))  
        uh2lob(ii,:) = phi(:) 
    ENDDO 
    !
    CONTINUE 
    !
    DEALLOCATE( LSQM, iLSQM, LSQrhs ) 
     !
     DO i = 1, nSubLim
         xilimbary(i) = 0.0 + 0.5/REAL(nSubLim) + (i-1)*1.0/REAL(nSubLim)
     ENDDO     
     !
     ! Reference element 
     !
     ReferenceElement(:,1) = (/ 0., 0., 0. /) 
     ReferenceElement(:,2) = (/ 1., 0., 0. /) 
     ReferenceElement(:,3) = (/ 1., 1., 0. /) 
     ReferenceElement(:,4) = (/ 0., 1., 0. /) 
     ReferenceElement(:,5) = (/ 0., 0., 1. /) 
     ReferenceElement(:,6) = (/ 1., 0., 1. /) 
     ReferenceElement(:,7) = (/ 1., 1., 1. /) 
     ReferenceElement(:,8) = (/ 0., 1., 1. /) 
    !
    ! Set the initial condition. Here, we assume a nodal basis. Otherwise, we would have to do formal L2 projection, 
    ! i.e. integration of the initial condition and multiplication with the inverse mass matrix. 
    !
    DO iElem = 1, nElem
        x0 = x(:,tri(1,iElem)) ! get the coordinate of the lower left node 
        DO k = 1, nDOF(3) 
         DO j = 1, nDOF(2)
          DO i = 1, nDOF(1) 
              xGP = x0 + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:) 
              CALL InitialField(u0,par0,xGP,0.0) 
              uh(:,i,j,k,iElem) = u0 
              parh(:,i,j,k,iElem) = par0  
          ENDDO
         ENDDO
        ENDDO
#ifdef LIMITER
        IF(N>0) THEN
            DO iVar = 1, nVar 
                Limiter(iElem)%lmin(iVar) = MINVAL(uh(iVar,:,:,:,iElem)) 
                Limiter(iElem)%lmax(iVar) = MAXVAL(uh(iVar,:,:,:,iElem)) 
            ENDDO         
            CALL DMP(dmpresult,uh(:,:,:,:,iElem),Limiter(iElem),1e-1)
            IF(.NOT.dmpresult) THEN
                ! If initial condition projection does not satisfy the DMP, then activate the subcell limiter 
                IF(Limiter(iElem)%oldstatus.EQ.0) THEN 
                   ALLOCATE( Limiter(iElem)%Lh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))    ) 
                   ALLOCATE( Limiter(iElem)%NewLh(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3)) ) 
                ENDIF
                Limiter(iElem)%status    = 1 
                Limiter(iElem)%oldstatus = 1 
                DO k = 1, nSubLimV(3) 
                 DO j = 1, nSubLimV(2)
                  DO i = 1, nSubLimV(1) 
                      xGP = x0 + (/ xilimbary(i), xilimbary(j), xilimbary(k) /)*dx(:) 
                      CALL InitialField(u0,par0,xGP,0.0) 
                      Limiter(iElem)%Lh(:,i,j,k) = u0 
                      Limiter(iElem)%NewLh(:,i,j,k) = u0 
                  ENDDO
                 ENDDO
                ENDDO             
            ENDIF              
        ENDIF
#endif 
    ENDDO    
    !
    IF(nParam.GT.0) THEN
        !
        ! Compute the bounday-extrapolated values for the material parameters 
        !
        DO iElem = 1, nElem 
            !
            parbnd(:,:,:,:,iElem) = 0. 
            ! x-direction: face 1 (left) and face 2 (right) 
            DO k = 1, nDOF(3) 
             DO j = 1, nDOF(2) 
                parbnd(:,1,j,k,iElem) = MATMUL( parh(:,:,j,k,iElem),   FLCoeff )   ! left 
                parbnd(:,2,j,k,iElem) = MATMUL( parh(:,:,j,k,iElem),   FRCoeff )   ! right 
             ENDDO
            ENDDO 
            ! y-direction: face 3 (left) and face 4 (right) 
            IF(nDim>=2) THEN
                DO k = 1, nDOF(3) 
                 DO i = 1, nDOF(1) 
                    parbnd(:,3,i,k,iElem) = MATMUL( parh(:,i,:,k,iElem),   FLCoeff )   ! left 
                    parbnd(:,4,i,k,iElem) = MATMUL( parh(:,i,:,k,iElem),   FRCoeff )   ! right 
                 ENDDO
                ENDDO 
            ENDIF    
            ! z-direction: face 5 (left) and face 6 (right) 
            IF(nDim>=3) THEN
                DO j = 1, nDOF(2) 
                 DO i = 1, nDOF(1) 
                    parbnd(:,5,i,j,iElem) = MATMUL( parh(:,i,j,:,iElem),   FLCoeff )   ! left 
                    parbnd(:,6,i,j,iElem) = MATMUL( parh(:,i,j,:,iElem),   FRCoeff )   ! right 
                 ENDDO
                ENDDO 
            ENDIF            
            !
        ENDDO
    ENDIF    
    !
    CALL WriteData
    !
    CONTINUE 
    !
END SUBROUTINE ADERDGInit
    
    
SUBROUTINE InitialField(u0,par,xGP,tGP) 
    USE typesDef
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN ) :: xGP(d), tGP        ! spatial position vector and time 
    REAL, INTENT(OUT) :: u0(nVar)           ! initial data vector in terms of conserved variables 
    REAL, INTENT(OUT) :: par(nParam)        ! material parameter vector 
    ! Local variables 
    REAL :: VBase(nVar), ampl(nVar), sigma(d) 
    REAL :: V0(nVar),r,VLL(nVar),VRR(nVar)   
    REAL :: du,dv,drho,dTemp,dp,epsilon
    ! 
#ifdef EULER
    ! no local material parameters for Euler equations 
    ! Gaussian perturbation 
    SELECT CASE(TRIM(ICType)) 
    CASE('SmoothWave')
        sigma = (/ 0.05, 0.05, 0.05 /)       ! half-width
        VBase(:) = (/ 1., 0., 0., 0., 1. /)  ! base-state 
        ampl(:)  = 0.                        ! perturbation amplitude vector 
        ampl(5)   = 1e-3                     ! 
        V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM(xGP(1:nDim)**2/sigma(1:nDim)**2) )
    CASE('Sod')     
        r = SQRT( SUM(xGP(1:nDim)**2) )     
        VLL = (/ 1.0,   0.0, 0.0, 0.0, 1.0 /) 
        VRR = (/ 0.125, 0.0, 0.0, 0.0, 0.1 /) 
        IF(xGP(1)<0) THEN 
            V0 = VLL
        ELSE
            V0 = VRR        
        ENDIF    
    CASE('EP2D')     
        r = SQRT( SUM(xGP(1:nDim)**2) )     
        VLL = (/ 1.0,   0.0, 0.0, 0.0, 1.0 /) 
        VRR = (/ 0.125, 0.0, 0.0, 0.0, 0.1 /) 
        IF(r.LT.0.5) THEN
        !IF(xGP(1)<0) THEN 
            V0 = VLL
        ELSE
            V0 = VRR        
        ENDIF    
    CASE('ShuVortex2D')
       epsilon = 5.0
       r = SQRT((xGP(1)-tGP-5.)**2+(xGP(2)-tGP-5.)**2)
       du = epsilon/2./EQN%Pi*exp(0.5*(1.-r*r))*(5. - xGP(2) + tGP)
       dv = epsilon/2./EQN%Pi*exp(0.5*(1.-r*r))*(xGP(1)  - 5.- tGP)
       dTemp = -(EQN%gamma-1.)*epsilon**2/8./EQN%gamma/EQN%Pi**2*exp(1.-r*r)
       drho = (1.+dTemp)**(1./(EQN%gamma-1.))-1.
       dp   = (1.+dTemp)**(EQN%gamma/(EQN%gamma-1.))-1.
       !
       V0(1) = 1. + drho
       V0(2) = 1. + du
       V0(3) = 1. + dv
       V0(4) = 0.0
       V0(5) = 1. + dp
       !
    END SELECT 
#endif
#ifdef ELASTICITY
    ! Set the local material parameters 
    par(1) = 2.0  ! lambda 
    par(2) = 1.0  ! mu 
    par(3) = 1.0  ! rho 
    ! Gaussian perturbation 
    sigma = (/ 0.05, 0.05, 0.05 /)       ! half-width
    VBase(:) = 0.0                       ! base-state 
    ampl(:)  = 0.                        ! perturbation amplitude vector 
    ampl(7)   = 1e-3                     ! 
    V0(:) = VBase(:) + ampl(:)*EXP( -0.5*SUM(xGP(1:nDim)**2/sigma(1:nDim)**2) )   
!    V0 = 0.  
!    V0(1) = 1.0 + exp(-0.5*(xGP(1)**2+xGP(2)**2)/0.1**2)  
#endif 
    ! A simple debug check for the computation of derivatives 
    !u0 = 0. 
    !u0(1) = 0.123 !*xGP(1) 
    CALL PDEPrim2Cons(u0,V0) 
    !
END SUBROUTINE InitialField
    
    
