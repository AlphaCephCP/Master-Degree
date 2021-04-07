SUBROUTINE EnforceCCZ4Constraints(luh)
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) 
    USE typesDef, ONLY : nVar, nDOF 
    IMPLICIT NONE
    ! Argument list
    REAL, INTENT(INOUT)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom 
    ! Local variables
    INTEGER :: i,j,k,l,iVar,iDim, iter
    REAL    :: Q(nVar) 
    REAL    :: g_cov(3,3), det, g_contr(3,3), Aex(3,3), traceA, phi  
    REAL    :: DD(3,3,3), traceDk 
    !
    DO k = 1, nDOF(3)
     DO j = 1, nDOF(2)
       DO i = 1, nDOF(1)
            Q = luh(:,i,j,k) 
            ! 
            g_cov(1,1) = Q(1)
            g_cov(1,2) = Q(2)
            g_cov(1,3) = Q(3)
            g_cov(2,1) = Q(2)
            g_cov(2,2) = Q(4)
            g_cov(2,3) = Q(5)
            g_cov(3,1) = Q(3)
            g_cov(3,2) = Q(5)
            g_cov(3,3) = Q(6)
            ! This determinant should be close to unity, since we use the conformal decomposition 
            det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
            g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
            g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
            g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
            g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
            g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
            g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
            g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
            g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
            g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
            !
            phi = det**(-1./6.) 
            g_cov = phi**2*g_cov
            det = (g_cov(1,1)*g_cov(2,2)*g_cov(3,3)-g_cov(1,1)*g_cov(2,3)*g_cov(3,2)-g_cov(2,1)*g_cov(1,2)*g_cov(3,3)+g_cov(2,1)*g_cov(1,3)*g_cov(3,2)+g_cov(3,1)*g_cov(1,2)*g_cov(2,3)-g_cov(3,1)*g_cov(1,3)*g_cov(2,2)) 
            g_contr(1,1) =  (g_cov(2,2)*g_cov(3,3)-g_cov(2,3)*g_cov(3,2)) / det 
            g_contr(1,2) = -(g_cov(1,2)*g_cov(3,3)-g_cov(1,3)*g_cov(3,2)) / det
            g_contr(1,3) = -(-g_cov(1,2)*g_cov(2,3)+g_cov(1,3)*g_cov(2,2))/ det 
            g_contr(2,1) = -(g_cov(2,1)*g_cov(3,3)-g_cov(2,3)*g_cov(3,1)) / det 
            g_contr(2,2) = (g_cov(1,1)*g_cov(3,3)-g_cov(1,3)*g_cov(3,1))  / det 
            g_contr(2,3) = -(g_cov(1,1)*g_cov(2,3)-g_cov(1,3)*g_cov(2,1)) / det 
            g_contr(3,1) = -(-g_cov(2,1)*g_cov(3,2)+g_cov(2,2)*g_cov(3,1))/ det 
            g_contr(3,2) = -(g_cov(1,1)*g_cov(3,2)-g_cov(1,2)*g_cov(3,1)) / det 
            g_contr(3,3) = (g_cov(1,1)*g_cov(2,2)-g_cov(1,2)*g_cov(2,1))  / det 
            !
            Aex(1,1) = Q(7) 
            Aex(1,2) = Q(8) 
            Aex(1,3) = Q(9) 
            Aex(2,1) = Q(8) 
            Aex(2,2) = Q(10) 
            Aex(2,3) = Q(11) 
            Aex(3,1) = Q(9) 
            Aex(3,2) = Q(11) 
            Aex(3,3) = Q(12)             
            !
            traceA = SUM(g_contr*Aex) 
            !
            Aex = Aex - 1./3.*g_cov*traceA 
            !
            luh( 1,i,j,k) = g_cov(1,1) 
            luh( 2,i,j,k) = g_cov(1,2) 
            luh( 3,i,j,k) = g_cov(1,3) 
            luh( 4,i,j,k) = g_cov(2,2) 
            luh( 5,i,j,k) = g_cov(2,3) 
            luh( 6,i,j,k) = g_cov(3,3) 
            !
            luh( 7,i,j,k) = Aex(1,1) 
            luh( 8,i,j,k) = Aex(1,2) 
            luh( 9,i,j,k) = Aex(1,3) 
            luh(10,i,j,k) = Aex(2,2) 
            luh(11,i,j,k) = Aex(2,3) 
            luh(12,i,j,k) = Aex(3,3)             
            !
            ! As suggested by our PRD referee, we also enforce the algebraic constraint that results from the first spatial derivative of the constraint
            ! det \tilde{\gamma}_ij = 0, which leads to
            !
            ! \tilde{\gamma}^{ij} D_kij = 0
            !
            ! and is thus a condition of trace-freeness on all submatrices D_kij for k=1,2,3.
            !
            DD(1,1,1)=Q(36)
            DD(1,1,2)=Q(37)
            DD(1,1,3)=Q(38)
            DD(1,2,1)=Q(37)
            DD(1,2,2)=Q(39)
            DD(1,2,3)=Q(40)
            DD(1,3,1)=Q(38)
            DD(1,3,2)=Q(40)
            DD(1,3,3)=Q(41)
            !
            DD(2,1,1)=Q(42)
            DD(2,1,2)=Q(43)
            DD(2,1,3)=Q(44)
            DD(2,2,1)=Q(43)
            DD(2,2,2)=Q(45)
            DD(2,2,3)=Q(46)
            DD(2,3,1)=Q(44)
            DD(2,3,2)=Q(46)
            DD(2,3,3)=Q(47)
            !
            DD(3,1,1)=Q(48)
            DD(3,1,2)=Q(49)
            DD(3,1,3)=Q(50)
            DD(3,2,1)=Q(49)
            DD(3,2,2)=Q(51)
            DD(3,2,3)=Q(52)
            DD(3,3,1)=Q(50)
            DD(3,3,2)=Q(52)
            DD(3,3,3)=Q(53)
            !
            DO l = 1, 3
                traceDk = SUM(g_contr*DD(l,:,:))
                DD(l,:,:) = DD(l,:,:) - 1./3.*g_cov*traceDk
            ENDDO
            !
            luh(36,i,j,k) = DD(1,1,1)
            luh(37,i,j,k) = DD(1,1,2)
            luh(38,i,j,k) = DD(1,1,3)
            luh(39,i,j,k) = DD(1,2,2)
            luh(40,i,j,k) = DD(1,2,3)
            luh(41,i,j,k) = DD(1,3,3)
            !
            luh(42,i,j,k) = DD(2,1,1)
            luh(43,i,j,k) = DD(2,1,2)
            luh(44,i,j,k) = DD(2,1,3)
            luh(45,i,j,k) = DD(2,2,2)
            luh(46,i,j,k) = DD(2,2,3)
            luh(47,i,j,k) = DD(2,3,3)
            !
            luh(48,i,j,k) = DD(3,1,1)
            luh(49,i,j,k) = DD(3,1,2)
            luh(50,i,j,k) = DD(3,1,3)
            luh(51,i,j,k) = DD(3,2,2)
            luh(52,i,j,k) = DD(3,2,3)
            luh(53,i,j,k) = DD(3,3,3)            
            !
       ENDDO
     ENDDO
    ENDDO    
#endif 
END SUBROUTINE EnforceCCZ4Constraints 
                
!
! The "Dicke Bertha" was an ultra-large caliber cannon developed by the company Krupp and used by the 
! German imperial army during the first world war. 
! 
! For the Einstein equations of General Relativity, we definitely need new ultra large caliber numerical methods 
! to cope with the complex linear and nonlinear momentum and Hamiltonian constraints, respectively.   
! The hyperbolic Z4 cleaning system is very nice and simple, but in the end, it does not do its job properly for
! general spacetimes!
! 
! Instead, our new nonlinear (!) constrained L2 projection algorithm is supposed to perform better. 
! It should guarantee the satisfaction of the Einstein constraints at least locally inside each cell in a weak form.    
!
SUBROUTINE DickeBertha(MaxNewton) 
    USE typesDef
    IMPLICIT NONE 
    INTEGER :: MaxNewton 
    INTEGER :: i,ii,jj,kk,count 
    REAL    :: lx0(d), ldx(d) 
    REAL    :: Cij(6,(N+1)**nDim), Kij(6,(N+1)**nDim) 
    REAL    :: R((N+1)**nDim), Christoffel(3,3,3,(N+1)**nDim), g_contr(3,3,(N+1)**nDim) 
    REAL    :: Q(nVar), gradQ(nVar,d), lwh(nVar,nDOF(1),nDOF(2),nDOF(3))
    REAL    :: lwx(nVar,nDOF(1),nDOF(2),nDOF(3))    
    REAL    :: lwy(nVar,nDOF(1),nDOF(2),nDOF(3))    
    REAL    :: lwz(nVar,nDOF(1),nDOF(2),nDOF(3))    
    !
    DO i = 1, nElem
      !
      lx0 = xL 
      ldx = dx 
      !
      DO kk = 1, nDOF(3) 
         DO jj = 1, nDOF(2)
            DO ii = 1, nDOF(1) 
               lwh(:,ii,jj,kk) = uh(:,ii,jj,kk,i)
            ENDDO
         ENDDO
      ENDDO   
      ! Compute the state derivatives in each point   
      DO kk = 1, nDOF(3)  
         DO jj = 1, nDOF(2) 
            lwx(:,:,jj,kk) = MATMUL( lwh(:,:,jj,kk), TRANSPOSE(dudx) )/ldx(1)  
         ENDDO
      ENDDO 
      ! Compute the state derivatives in each point   
      IF(nDim>=2) THEN
          DO kk = 1, nDOF(3)  
             DO jj = 1, nDOF(1) 
                lwy(:,jj,:,kk) = MATMUL( lwh(:,jj,:,kk), TRANSPOSE(dudx) )/ldx(2)  
             ENDDO
          ENDDO 
      ELSE
          lwy = 0.0
      ENDIF 
      ! Compute the state derivatives in each point   
      IF(nDim>=3) THEN
          DO kk = 1, nDOF(2)  
             DO jj = 1, nDOF(1) 
                lwz(:,jj,kk,:) = MATMUL( lwh(:,jj,kk,:), TRANSPOSE(dudx) )/ldx(3)  
             ENDDO
          ENDDO       
      ELSE
          lwz = 0.0
     ENDIF 
     !
     count = 0 
     DO kk = 1, nDOF(3)  
      DO jj = 1, nDOF(2)  
       DO ii = 1, nDOF(1) 
         count = count + 1 
         Q(:) = lwh(:,ii,jj,kk)
         gradQ(:,1) = lwx(:,ii,jj,kk)  
         gradQ(:,2) = lwy(:,ii,jj,kk)  
         gradQ(:,3) = lwz(:,ii,jj,kk)  
         Kij(:,count) = Q(7:12) 
         CALL ComputeMetricTerms(g_contr(:,:,count),Christoffel(:,:,:,count),R(count),Q,gradQ)             
         CONTINUE  
       ENDDO
      ENDDO
     ENDDO
     !
     ! We enforce the Einstein constraints by using constrained L2 projection.
     ! This is quasi exact (in all quadrature points) for the nonlinear Hamiltonian constraint.    
     ! It is only exact in a weak integral sense for the linear momentum constraints. 
     !  
     CALL EinsteinConstraints(Cij,Kij,R,g_contr,Christoffel,ldx,MaxNewton) 
     count = 0 
     DO kk = 1, nDOF(3)  
      DO jj = 1, nDOF(2)  
       DO ii = 1, nDOF(1) 
         count = count + 1 
         Q(7:12) = Cij(:,count) 
         uh(7:12,ii,jj,kk,i) = Q(7:12)  
       ENDDO
      ENDDO
     ENDDO 
     !
    ENDDO
    !  
END SUBROUTINE DickeBertha 

SUBROUTINE EinsteinConstraints(Cij,Kij,R,g_contr,Christoffel,ldx,MaxNewtonIn)  
    USE typesDef, ONLY : N, d, nDim, nDOF, wGPN, xiGPN      
    IMPLICIT NONE  
    INTEGER       :: nDegF, nEqn, nConstr, iErr, count, nDegFrRec, c1, c2, ii, jj, kk, ll, mm, pp, qq, rr, iii, jjj, kkk 
    INTEGER       :: idxc(3,(N)**nDim), idxh((N+1)**nDim), idxd(6,(N+1)**nDim), idxk(3,3), idx1, idx2, nGPVNm1(d)   
    INTEGER       :: iNewton, MaxNewtonIn, MaxNewton      
    REAL          :: res, tol = 1e-13  
    REAL          :: ldx(d)   
    REAL, POINTER :: MMat(:,:), AMat(:,:), CMat(:,:), RMat(:,:)
    REAL          :: Cij(6,(N+1)**nDim), Kij(6,(N+1)**nDim), R((N+1)**nDim), Christoffel(3,3,3,(N+1)**nDim), g_contr(3,3,(N+1)**nDim)
    REAL          :: lChristoffel(3,3,3), lg_contr(3,3), lK(3,3), lKup(3,3)    
    REAL          :: MassMatrix((N+1)**nDim,(N+1)**nDim) 
    REAL          :: xi, eta, zeta, wGP, K2((N+1)**nDim), traceK((N+1)**nDim)    
    REAL          :: psi_i(N+1), psi_j(N+1), psi_k(N+1), psi_i_xi(N+1), psi_j_xi(N+1), psi_k_xi(N+1), aux(d)
    REAL          :: phi_i(N), phi_j(N), phi_k(N), phi_xi(N), phi(N**nDim)   
    REAL          :: psi((N+1)**nDim), Dpsi(d,(N+1)**nDim), check(3,(N+1)**nDim), checkHam, test(3*N**nDim+(N+1)**nDim)      
    !
    ! Input: nEqn   = number of equations ( > nDegF ) 
    !        nDegF  = number of unknowns  
    !        nConst = number of constraints 
    !        A      = the (nDegF,nDegF)   system matrix of (1) 
    !        C      = the (nConst,nDegF) matrix of the linear constraints (2) 
    !        R      = the (nConst,nEqn) matrix of the linear constraints (2) 
    ! 
    ! Output: M    = the (nDegF,nEqn) matrix, which yields the solution x when multiplied with the right hand side vector b 
    !
    nGPVNm1 = MAX(1, nDOF(1:d)-1)  
    !
    !CALL CheckHamiltonianConstraint(checkHam,Kij,R,g_contr,ldx)
    !
    !IF(checkHam>1e-7) THEN
    !    MaxNewton = 6 
    !ELSE
    !    MaxNewton = 1 
    !ENDIF 
    !
    MaxNewton = MaxNewtonIn 
    !
    ! Symmetry of the tensor K_ij 
    ! 
    idxk(:,1) = (/ 1, 2, 3 /) 
    idxk(:,2) = (/ 2, 4, 5 /) 
    idxk(:,3) = (/ 3, 5, 6 /) 
    !
    !
    !
    ! Numbering of the momentum constraints (Lagrange multiplier mu_ki) 
    ! 
    c1 = 0
    DO mm = 1, 3    
     c2 = 0 
     DO kk = 1, nGPVNm1(3)   
      DO jj = 1, nGPVNm1(2)   
       DO ii = 1, nGPVNm1(1)   
         c1 = c1 + 1 
         c2 = c2 + 1 
         idxc(mm,c2) = c1 
       ENDDO
      ENDDO
     ENDDO
    ENDDO 
    !
    ! Numbering of the Hamiltonian constraint 
    !  
    c2 = 0 
    DO kk = 1, nDOF(3)   
     DO jj = 1, nDOF(2)   
      DO ii = 1, nDOF(1)   
         c2 = c2 + 1 
         c1 = c1 + 1 
         idxh(c2) = c1 
      ENDDO
     ENDDO
    ENDDO
    !
    nDegF    = 6*(N+1)**nDim 
    IF(MaxNewton>1) THEN 
        nEqn    = 6*(N+1)**nDim + 1*(N+1)**nDim 
        nConstr = 3*N**nDim     + 1  
    ELSE
        nEqn    = 6*(N+1)**nDim   
        nConstr = 3*N**nDim       
    ENDIF 
    ALLOCATE( MMat(nDegF,nEqn), AMat(nDegF,nDegF), CMat(nConstr,nDegF), RMat(nConstr,nEqn) )  
    !
    MassMatrix = 0.0 
    count = 0   
    DO kk = 1, nDOF(3)  
     DO jj = 1, nDOF(2)  
      DO ii = 1, nDOF(1) 
        count = count + 1 
        MassMatrix(count,count) = wGPN(ii)*wGPN(jj)*wGPN(kk) 
      ENDDO
     ENDDO
    ENDDO
    !
    nDegFrRec = (N+1)**nDim 
    AMat = 0.0
    count = 0 
    DO kk = 1, 6 
     DO ii = 1, nDegFrRec 
       count = count + 1 
       idxd(kk,ii) = count 
       AMat(count,count) = MassMatrix(ii,ii)     
     ENDDO
    ENDDO        
    !
    Cij = Kij       ! initial guess 
    !
    DO iNewton = 1, MaxNewton  
    !
    CMat = 0.0 
    RMat = 0.0 
    !
    c1   = 0 
    DPsi = 0.0 
    DO kkk = 1, nDOF(3)   
     DO jjj = 1, nDOF(2)   
      DO iii = 1, nDOF(1)  
        c1 = c1 + 1 
        CALL Mm1BaseFunc1D(phi_i,phi_xi,xiGPN(iii))    
        CALL Mm1BaseFunc1D(phi_j,phi_xi,xiGPN(jjj))    
        CALL Mm1BaseFunc1D(phi_k,phi_xi,xiGPN(kkk))
        count = 0 
        DO kk = 1, nGPVNm1(3)
         DO jj = 1, nGPVNm1(2)
          DO ii = 1, nGPVNm1(1) 
            count = count + 1 
            aux = (/ phi_i(ii), phi_j(jj), phi_k(kk) /) 
            phi(count) = PRODUCT( aux(1:nDim) )
          ENDDO
         ENDDO
        ENDDO 
        ! 
        CALL BaseFunc1D(psi_i,psi_i_xi,xiGPN(iii))    
        CALL BaseFunc1D(psi_j,psi_j_xi,xiGPN(jjj))    
        CALL BaseFunc1D(psi_k,psi_k_xi,xiGPN(kkk)) 
        count = 0 
        DO kk = 1, nDOF(3)
         DO jj = 1, nDOF(2)
          DO ii = 1, nDOF(1) 
            count = count + 1 
            aux = (/ psi_i(ii), psi_j(jj), psi_k(kk) /) 
            psi(count) = PRODUCT( aux(1:nDim) )
            aux = (/ psi_i_xi(ii)/ldx(1), psi_j(jj), psi_k(kk) /)  
            DPsi(1,count) = PRODUCT( aux(1:nDim) ) 
            IF(nDim>=2) THEN
                aux = (/ psi_i(ii), psi_j_xi(jj)/ldx(2), psi_k(kk) /)  
                DPsi(2,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(2,count) = 0.0 
            ENDIF 
            IF(nDim>=3) THEN
                aux = (/ psi_i(ii), psi_j(jj), psi_k_xi(kk)/ldx(3) /)  
                DPsi(3,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(3,count) = 0.0 
            ENDIF 
          ENDDO
         ENDDO
        ENDDO
        !
        ! Compute the contravariant metric tensor and the Christoffel symbol in the current Gaussian quadrature point 
        ! 
        lChristoffel = 0.0
        lg_contr     = 0.0 
        lK(1,1)      = SUM( psi(:)*Cij(1,:) );  lK(1,2)      = SUM( psi(:)*Cij(2,:) );  lK(1,3)      = SUM( psi(:)*Cij(3,:) ) 
        lK(2,1)      = SUM( psi(:)*Cij(2,:) );  lK(2,2)      = SUM( psi(:)*Cij(4,:) );  lK(2,3)      = SUM( psi(:)*Cij(5,:) ) 
        lK(3,1)      = SUM( psi(:)*Cij(3,:) );  lK(3,2)      = SUM( psi(:)*Cij(5,:) );  lK(3,3)      = SUM( psi(:)*Cij(6,:) ) 
        !
        DO ii = 1, 3 
         DO jj = 1, 3 
          lg_contr(ii,jj) = SUM( psi(:)*g_contr(ii,jj,:) )   
          DO kk = 1, 3 
            lChristoffel(ii,jj,kk) = SUM( psi(:)*Christoffel(ii,jj,kk,:) )   
          ENDDO
         ENDDO
        ENDDO 
        !
        lKup = 0.0 
        DO ii = 1, 3 
         DO jj = 1, 3 
          DO pp = 1, 3
           DO qq = 1, 3
              lKup(ii,jj) = lKup(ii,jj) + lg_contr(ii,pp)*lg_contr(jj,qq)*lK(pp,qq) 
           ENDDO
          ENDDO
         ENDDO
        ENDDO 
        traceK(c1) = 0.0
        K2(c1)     = 0.0
        DO ii = 1, 3
         DO jj = 1, 3 
            traceK(c1) = traceK(c1) + lg_contr(ii,jj)*lK(ii,jj) 
            K2(c1)     = K2(c1)     + lKup(ii,jj)*lK(ii,jj) 
         ENDDO
        ENDDO           
        !
        CONTINUE 
        !
        ! Gaussian weight in the respective quadrature point 
        aux = (/ wGPN(iii), wGPN(jjj), wGPN(kkk) /) 
        wGP = PRODUCT( aux(1:nDim) ) 
        !
        ! Momentum constraints 
        !
        DO kk = 1, N**nDim 
         DO rr = 1, (N+1)**nDim
          DO ii = 1, 3 
           DO jj = 1, 3
            DO ll = 1, 3
             idx1 = idxc(ii,kk) ! mu_ki 
             idx2 = idxd(idxk(ii,jj),rr) 
             CMat(idx1,idx2) = CMat(idx1,idx2) + wGP*phi(kk)*lg_contr(jj,ll)*DPsi(ll,rr) ! * C_ij,r    
             idx2 = idxd(idxk(jj,ll),rr) 
             CMat(idx1,idx2) = CMat(idx1,idx2) - wGP*phi(kk)*lg_contr(jj,ll)*DPsi(ii,rr) ! * C_jl,r 
             DO mm = 1, 3
                idx2 = idxd(idxk(mm,ii),rr) 
                CMat(idx1,idx2) = CMat(idx1,idx2) - wGP*phi(kk)*lg_contr(jj,ll)*lChristoffel(jj,ll,mm)*Psi(rr) ! * C_mi,r 
                idx2 = idxd(idxk(mm,ll),rr)    
                CMat(idx1,idx2) = CMat(idx1,idx2) + wGP*phi(kk)*lg_contr(jj,ll)*lChristoffel(jj,ii,mm)*Psi(rr) ! * C_ml,r
             ENDDO  
             ! 
            ENDDO 
           ENDDO 
          ENDDO
         ENDDO 
        ENDDO     
        !
        !
        IF(MaxNewton > 1) THEN
         !
         ! Linearized Hamiltonian constraint (Newton method) 
         !
         ! H = R - K_ij K^ij + K^2 = 0 
         ! HLin = H^n + dH/dK (K-K^n) = 0 
         ! HLin = R - C_ij^n C^ij^n + (C^n)^2 - 2* g^ip g^jq C_pq^n (C_ij - C_ij^n) + 2 g^ij g^pq C_pq^n ( C_ij - C_ij^n ) 
         !
         kk = 1                         ! only enforce the constraint in the integral mean sense 
         !DO kk = 1, (N+1)**nDim 
           idx1 = idxh(kk)                 ! lambda_k  
           !RMat(idx1,nDOF+kk) = -wGP   
           DO rr = 1, (N+1)**nDim
            RMat(idx1,nDegF+rr) = RMat(idx1,nDegF+rr) - wGP*Psi(rr) ! * ( R + K_ij K^ij - K^2 )    
            DO ii = 1, 3 
             DO jj = 1, 3
              DO pp = 1, 3
               DO qq = 1, 3 
               !              
               idx2 = idxd(idxk(ii,jj),rr) 
               ! Newton
               CMat(idx1,idx2) = CMat(idx1,idx2) - wGP*2.0*lg_contr(ii,pp)*lg_contr(jj,qq)*lK(pp,qq)*Psi(rr) ! * C_ij,r     
               CMat(idx1,idx2) = CMat(idx1,idx2) + wGP*2.0*lg_contr(ii,jj)*lg_contr(pp,qq)*lK(pp,qq)*Psi(rr) ! * C_ij,r   
               !
               ENDDO   
              ENDDO 
             ENDDO 
            ENDDO
           ENDDO           
          !ENDDO 
          !
          ENDIF 
        !
      ENDDO
     ENDDO
    ENDDO 
    ! 
    CALL CL2(MMat,nEqn,nDegF,nConstr,AMat,CMat,RMat,iErr) 
    !
    Cij = 0.0 
    c1 = 0 
    DO mm = 1, 6  
     DO ii = 1, nDegFrRec 
      c1 = c1 + 1 
      count = 0 
      DO kk = 1, 6
       DO jj = 1, (N+1)**nDim
        count = count + 1  
        Cij(mm,ii) = Cij(mm,ii) + MMat(c1,count)*Kij(kk,jj) 
       ENDDO 
      ENDDO
      IF(MaxNewton>1) THEN
          DO jj = 1, (N+1)**nDim
            count = count + 1 
            ! Newton 
            Cij(mm,ii) = Cij(mm,ii) + MMat(c1,count)*( R(jj) - traceK(jj)**2 + K2(jj) )    
          ENDDO 
      ENDIF 
      !   
     ENDDO 
    ENDDO  
    !
    test = 0.0
    DO ii = 1, nConstr
     count = 0 
     DO kk = 1, 6
      DO jj = 1, (N+1)**nDim
        count = count + 1  
        test(ii) = test(ii) + CMat(ii,count)*Cij(kk,jj) 
      ENDDO  
     ENDDO
     IF(MaxNewton > 1) THEN
         DO jj = 1, (N+1)**nDim
            count = count + 1 
            ! Newton 
            test(ii) = test(ii) - RMat(ii,count)*( R(jj) - traceK(jj)**2 + K2(jj) )  
         ENDDO       
     ENDIF 
    ENDDO  
    !
    CALL CheckMomentumConstraints(check,Cij,g_contr,Christoffel,ldx)
    CALL CheckHamiltonianConstraint(checkHam,Cij,R,g_contr,ldx)
    !
    res = ABS(checkHam) 
    IF( res < tol ) THEN
        EXIT
    ENDIF 
    !
    CONTINUE
    !
    ENDDO 
    !
    IF(MaxNewton>1) THEN
        IF(res > tol) THEN
            PRINT *, ' Newton did not converge for Hamiltonian constraint ', res, tol
            STOP 
        ENDIF    
    ENDIF 
    !
    DEALLOCATE( MMat, AMat, CMat, RMat )  
    !
    CONTINUE
    ! 
END SUBROUTINE EinsteinConstraints 

SUBROUTINE ComputeMetricTerms(g_contr,Christoffel,R,Q,gradQ)
    USE typesDef, ONLY : nVar, d  
    IMPLICIT NONE
    INTEGER :: ii,jj,kk,ll,mm,nn,qq,rr 
    REAL :: Q(nVar), gradQ(nVar,d) 
    REAL :: R, Christoffel(3,3,3), g_cov(3,3), g_contr(3,3), Riemann(3,3,3,3), Ricci(3,3),det  
    REAL :: dgup(3,3,3), dDD(3,3,3,3), dChristoffel(3,3,3,3)  
    REAL :: K(3,3), DD(3,3,3), dZ(3,3)    
    REAL :: Qx(nVar), Qy(nVar), Qz(nVar) 
    !
#ifdef Z4EINSTEIN     
    !
    Qx = gradQ(:,1) 
    Qy = gradQ(:,2) 
    Qz = gradQ(:,3) 
    ! 
    g_cov( 1, 1:3) = (/ Q(1),  Q(2),  Q(3)   /) 
    g_cov( 2, 1:3) = (/ Q(2),  Q(4),  Q(5)   /) 
    g_cov( 3, 1:3) = (/ Q(3),  Q(5),  Q(6)   /)    
    det = (Q(1)*Q(4)*Q(6)-Q(1)*Q(5)**2-Q(2)**2*Q(6)+2*Q(2)*Q(3)*Q(5)-Q(3)**2*Q(4)) 
    g_contr(1,1) = (Q(4)*Q(6)-Q(5)**2)   / det 
    g_contr(1,2) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(1,3) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(2,1) = -(Q(2)*Q(6)-Q(3)*Q(5))/ det 
    g_contr(2,2) = (Q(1)*Q(6)-Q(3)**2)   / det 
    g_contr(2,3) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,1) = (Q(2)*Q(5)-Q(3)*Q(4)) / det 
    g_contr(3,2) = -(Q(1)*Q(5)-Q(2)*Q(3))/ det 
    g_contr(3,3) = (Q(1)*Q(4)-Q(2)**2)   / det 
    DD(1,:,:) = RESHAPE ( (/Q(36), Q(37), Q(38), Q(37), Q(39), Q(40), Q(38), Q(40), Q(41)  /)  , (/3,3/) )
    DD(2,:,:) = RESHAPE ( (/Q(42), Q(43), Q(44), Q(43), Q(45), Q(46), Q(44), Q(46), Q(47)  /)  , (/3,3/) )
    DD(3,:,:) = RESHAPE ( (/Q(48), Q(49), Q(50), Q(49), Q(51), Q(52), Q(50), Q(52), Q(53)  /)  , (/3,3/) )   
    !
    dDD(1,1,:,:)   = RESHAPE ( (/gradQ(36,1), gradQ(37,1), gradQ(38,1), gradQ(37,1), gradQ(39,1), gradQ(40,1), gradQ(38,1), gradQ(40,1), gradQ(41,1)  /)  , (/3,3/) )
    dDD(2,1,:,:)   = RESHAPE ( (/gradQ(36,2), gradQ(37,2), gradQ(38,2), gradQ(37,2), gradQ(39,2), gradQ(40,2), gradQ(38,2), gradQ(40,2), gradQ(41,2)  /)  , (/3,3/) )
    dDD(3,1,:,:)   = RESHAPE ( (/gradQ(36,3), gradQ(37,3), gradQ(38,3), gradQ(37,3), gradQ(39,3), gradQ(40,3), gradQ(38,3), gradQ(40,3), gradQ(41,3)  /)  , (/3,3/) )   
    ! 
    dDD(1,2,:,:)   = RESHAPE ( (/gradQ(42,1), gradQ(43,1), gradQ(44,1), gradQ(43,1), gradQ(45,1), gradQ(46,1), gradQ(44,1), gradQ(46,1), gradQ(47,1)  /)  , (/3,3/) )
    dDD(2,2,:,:)   = RESHAPE ( (/gradQ(42,2), gradQ(43,2), gradQ(44,2), gradQ(43,2), gradQ(45,2), gradQ(46,2), gradQ(44,2), gradQ(46,2), gradQ(47,2)  /)  , (/3,3/) )
    dDD(3,2,:,:)   = RESHAPE ( (/gradQ(42,3), gradQ(43,3), gradQ(44,3), gradQ(43,3), gradQ(45,3), gradQ(46,3), gradQ(44,3), gradQ(46,3), gradQ(47,3)  /)  , (/3,3/) )
    ! 
    dDD(1,3,:,:)   = RESHAPE ( (/gradQ(48,1), gradQ(49,1), gradQ(50,1), gradQ(49,1), gradQ(51,1), gradQ(52,1), gradQ(50,1), gradQ(52,1), gradQ(53,1)  /)  , (/3,3/) )
    dDD(2,3,:,:)   = RESHAPE ( (/gradQ(48,2), gradQ(49,2), gradQ(50,2), gradQ(49,2), gradQ(51,2), gradQ(52,2), gradQ(50,2), gradQ(52,2), gradQ(53,2)  /)  , (/3,3/) )
    dDD(3,3,:,:)   = RESHAPE ( (/gradQ(48,3), gradQ(49,3), gradQ(50,3), gradQ(49,3), gradQ(51,3), gradQ(52,3), gradQ(50,3), gradQ(52,3), gradQ(53,3)  /)  , (/3,3/) )  
    !
    dZ(1,:) = gradQ(13:15,1) 
    dZ(2,:) = gradQ(13:15,2) 
    dZ(3,:) = gradQ(13:15,3) 
    !
    ! Compute the Christoffel symbol 
    Christoffel   = 0.0 
    DO ii = 1, 3 
     DO jj = 1, 3  
      DO kk = 1, 3 
       DO ll = 1, 3 
            Christoffel(ii,jj,kk)   = Christoffel(ii,jj,kk)   + g_contr(kk,ll)*(DD(ii,jj,ll) + DD(jj,ii,ll) - DD(ll,ii,jj))
       ENDDO
      ENDDO
     ENDDO
    ENDDO   
    !
    dgup   = 0.0 
    DO ii = 1, 3 
     DO jj = 1, 3  
      DO kk = 1, 3 
        DO ll = 1, 3 
          dgup(ii,jj,kk) = dgup(ii,jj,kk) - Christoffel(ii,ll,jj)*g_contr(kk,ll) - Christoffel(ii,ll,kk)*g_contr(jj,ll) 
        ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
   dChristoffel = 0.0 
   DO ii = 1, 3 
    DO qq = 1, 3 
     DO mm = 1, 3 
      DO kk = 1, 3
       DO ll = 1, 3            
         !dChristoffel(kk,ii,qq,mm) = dChristoffel(kk,ii,qq,mm) + dgup(kk,mm,ll)*0.5*( dg_cov(ii,qq,ll) + dg_cov(qq,ii,ll) - dg_cov(ll,ii,qq))                                      & 
         !                                                      + 0.5*g_contr(mm,ll)*( dDD(kk,ii,qq,ll) + dDD(kk,qq,ii,ll) - dDD(kk,ll,ii,qq))                                      & 
         !                                                      + 0.5*g_contr(mm,ll)*( dDD(ii,kk,qq,ll) + dDD(qq,kk,ii,ll) - dDD(ll,kk,ii,qq)) 
         dChristoffel(kk,ii,qq,mm) = dChristoffel(kk,ii,qq,mm) + dgup(kk,mm,ll)*( DD(ii,qq,ll) + DD(qq,ii,ll) - DD(ll,ii,qq))                                                      & 
                                                               + 0.5*g_contr(mm,ll)*( dDD(kk,ii,qq,ll) + dDD(kk,qq,ii,ll) - dDD(kk,ll,ii,qq))                                      & 
                                                               + 0.5*g_contr(mm,ll)*( dDD(ii,kk,qq,ll) + dDD(qq,kk,ii,ll) - dDD(ll,kk,ii,qq)) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO    
   !
   ! Compute the Riemann tensor 
   !
   Riemann = 0.0 
   DO ii = 1, 3 
    DO qq = 1, 3 
     DO mm = 1, 3 
      DO kk = 1, 3         
       Riemann(ii,kk,qq,mm) = dChristoffel(kk,ii,qq,mm) - dChristoffel(qq,ii,kk,mm) 
       DO jj = 1, 3 
          Riemann(ii,kk,qq,mm) = Riemann(ii,kk,qq,mm) + Christoffel(ii,qq,jj)*Christoffel(jj,kk,mm) - Christoffel(ii,kk,jj)*Christoffel(jj,qq,mm) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   !
   ! Compute the 3-Ricci tensor 
   !
   Ricci = 0.0 
   DO mm = 1, 3
    DO nn = 1, 3 
     DO ll = 1, 3 
        Ricci(mm,nn) = Ricci(mm,nn) + Riemann(mm,ll,nn,ll)  
     ENDDO
    ENDDO
   ENDDO   
   ! 
   ! Compute the 3-Ricci scalar
   ! 
   R = 0.0
   DO ii = 1, 3 
       DO jj = 1, 3
           R = R + g_contr(ii,jj)*( Ricci(ii,jj) + 2*dZ(ii,jj) ) 
       ENDDO
   ENDDO   
   !
#endif 
   !
END SUBROUTINE ComputeMetricTerms 


SUBROUTINE CheckMomentumConstraints(check,Cij,g_contr,Christoffel,ldx)  
    USE typesDef, ONLY : N, d, nDim, nDOF, xiGPN, wGPN, xiGPMm1, wGPMm1    
    IMPLICIT NONE  
    INTEGER       :: nDegF, nEqn, nConstr, iErr, count, nDegFrRec, c1, c2, ii, jj, kk, ll, mm, rr, iii, jjj, kkk    
    INTEGER       :: idxc(3,(N)**nDim), idxd(6,(N+1)**nDim), idxk(3,3), idx1, idx2, nGPVNm1(d)   
    REAL          :: check(3,(N+1)**nDim), ldx(d) 
    REAL, POINTER :: MMat(:,:), AMat(:,:), CMat(:,:), RMat(:,:)
    REAL          :: Cij(6,(N+1)**nDim), Kij(6,(N+1)**nDim), R((N+1)**nDim), Christoffel(3,3,3,(N+1)**nDim), g_contr(3,3,(N+1)**nDim)
    REAL          :: lChristoffel(3,3,3), lg_contr(3,3)  
    REAL          :: MassMatrix((N+1)**nDim,(N+1)**nDim) 
    REAL          :: xi, eta, zeta, wGP   
    REAL          :: psi_i(N+1), psi_j(N+1), psi_k(N+1), psi_i_xi(N+1), psi_j_xi(N+1), psi_k_xi(N+1), aux(d)
    REAL          :: psi((N+1)**nDim), Dpsi(d,(N+1)**nDim)    
    !
#ifdef Z4EINSTEIN     
    !
    nGPVNm1 = MAX(1, nDOF(1:d)-1)  
    !
    ! Symmetry of the tensor K_ij 
    ! 
    idxk(:,1) = (/ 1, 2, 3 /) 
    idxk(:,2) = (/ 2, 4, 5 /) 
    idxk(:,3) = (/ 3, 5, 6 /) 
    !
    check = 0.0  
    !
    c1   = 0 
    DPsi = 0.0 
    DO kkk = 1, nDOF(3)   
     DO jjj = 1, nDOF(2)   
      DO iii = 1, nDOF(1)  
        c1 = c1 + 1 
        CALL BaseFunc1D(psi_i,psi_i_xi,xiGPN(iii))    
        CALL BaseFunc1D(psi_j,psi_j_xi,xiGPN(jjj))    
        CALL BaseFunc1D(psi_k,psi_k_xi,xiGPN(kkk)) 
        count = 0 
        DO kk = 1, nDOF(3)
         DO jj = 1, nDOF(2)
          DO ii = 1, nDOF(1) 
            count = count + 1 
            aux = (/ psi_i(ii), psi_j(jj), psi_k(kk) /) 
            psi(count) = PRODUCT( aux(1:nDim) )
            aux = (/ psi_i_xi(ii)/ldx(1), psi_j(jj), psi_k(kk) /)  
            DPsi(1,count) = PRODUCT( aux(1:nDim) ) 
            IF(nDim>=2) THEN
                aux = (/ psi_i(ii), psi_j_xi(jj)/ldx(2), psi_k(kk) /)  
                DPsi(2,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(2,count) = 0.0 
            ENDIF            
            IF(nDim>=3) THEN
                aux = (/ psi_i(ii), psi_j(jj), psi_k_xi(kk)/ldx(3) /)  
                DPsi(3,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(3,count) = 0.0 
            ENDIF
          ENDDO
         ENDDO
        ENDDO
        !
        ! Compute the contravariant metric tensor and the Christoffel symbol in the current Gaussian quadrature point 
        ! 
        lChristoffel = 0.0
        lg_contr     = 0.0 
        DO ii = 1, 3 
         DO jj = 1, 3 
          lg_contr(ii,jj) = SUM( psi(:)*g_contr(ii,jj,:) )   
          DO kk = 1, 3 
            lChristoffel(ii,jj,kk) = SUM( psi(:)*Christoffel(ii,jj,kk,:) )   
          ENDDO
         ENDDO
        ENDDO 
        !
        CONTINUE 
        !
        DO rr = 1, (N+1)**nDim
         DO ii = 1, 3 
          DO jj = 1, 3
           DO ll = 1, 3
             check(ii,c1) = check(ii,c1) + lg_contr(jj,ll)*DPsi(ll,rr)*Cij(idxk(ii,jj),rr)  ! * C_ij,r    
             check(ii,c1) = check(ii,c1) - lg_contr(jj,ll)*DPsi(ii,rr)*Cij(idxk(jj,ll),rr)  ! * C_jl,r 
             DO mm = 1, 3
                check(ii,c1) = check(ii,c1) - lg_contr(jj,ll)*lChristoffel(jj,ll,mm)*Psi(rr)*Cij(idxk(mm,ii),rr) ! * C_mi,r 
                check(ii,c1) = check(ii,c1) + lg_contr(jj,ll)*lChristoffel(jj,ii,mm)*Psi(rr)*Cij(idxk(mm,ll),rr) ! * C_ml,r
             ENDDO  
             ! 
           ENDDO 
          ENDDO 
         ENDDO
        ENDDO     
        !
      ENDDO
     ENDDO
    ENDDO 
    !
#endif     
    !
END SUBROUTINE CheckMomentumConstraints 

SUBROUTINE CheckHamiltonianConstraint(check,Cij,R,g_contr,ldx)  
    USE typesDef, ONLY : N, d, nDim, nDOF, xiGPN, wGPN, xiGPMm1, wGPMm1    
    IMPLICIT NONE  
    INTEGER       :: nDegF, nEqn, nConstr, iErr, count, nDegFrRec, c1, c2, ii, jj, kk, ll, mm, pp, qq, rr, iii, jjj, kkk    
    INTEGER       :: idxc(3,(N)**nDim), idxd(6,(N+1)**nDim), idxk(3,3), idx1, idx2, nGPVNm1(d)   
    REAL          :: check, ldx(d), aux(d) 
    REAL, POINTER :: MMat(:,:), AMat(:,:), CMat(:,:), RMat(:,:)
    REAL          :: Cij(6,(N+1)**nDim), Kij(6,(N+1)**nDim), R((N+1)**nDim), g_contr(3,3,(N+1)**nDim)
    REAL          :: lChristoffel(3,3,3), lg_contr(3,3), K(3,3), Kup(3,3), traceK, K2  
    REAL          :: MassMatrix((N+1)**nDim,(N+1)**nDim) 
    REAL          :: xi, eta, zeta, wGP   
    REAL          :: psi_i(N+1), psi_j(N+1), psi_k(N+1), psi_i_xi(N+1), psi_j_xi(N+1), psi_k_xi(N+1) 
    REAL          :: psi((N+1)**nDim), Dpsi(d,(N+1)**nDim)    
    !
#ifdef Z4EINSTEIN 
    !
    nGPVNm1 = MAX(1, nDOF(1:d)-1)  
    !
    ! Symmetry of the tensor K_ij 
    ! 
    idxk(:,1) = (/ 1, 2, 3 /) 
    idxk(:,2) = (/ 2, 4, 5 /) 
    idxk(:,3) = (/ 3, 5, 6 /) 
    !
    check = 0.0  
    !
    c1   = 0 
    DPsi = 0.0 
    DO kkk = 1, nDOF(3)   
     DO jjj = 1, nDOF(2)   
      DO iii = 1, nDOF(1)  
        c1 = c1 + 1 
        CALL BaseFunc1D(psi_i,psi_i_xi,xiGPN(iii))    
        CALL BaseFunc1D(psi_j,psi_j_xi,xiGPN(jjj))    
        CALL BaseFunc1D(psi_k,psi_k_xi,xiGPN(kkk)) 
        count = 0 
        DO kk = 1, nDOF(3)
         DO jj = 1, nDOF(2)
          DO ii = 1, nDOF(1) 
            count = count + 1 
            aux = (/ psi_i(ii), psi_j(jj), psi_k(kk) /) 
            psi(count) = PRODUCT( aux(1:nDim) )
            aux = (/ psi_i_xi(ii)/ldx(1), psi_j(jj), psi_k(kk) /)  
            DPsi(1,count) = PRODUCT( aux(1:nDim) ) 
            IF(nDim>=2) THEN
                aux = (/ psi_i(ii), psi_j_xi(jj)/ldx(2), psi_k(kk) /)  
                DPsi(2,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(2,count) = 0.0 
            ENDIF            
            IF(nDim>=3) THEN
                aux = (/ psi_i(ii), psi_j(jj), psi_k_xi(kk)/ldx(3) /)  
                DPsi(3,count) = PRODUCT( aux(1:nDim) )  
            ELSE
                DPsi(3,count) = 0.0 
            ENDIF
          ENDDO
         ENDDO
        ENDDO
        !
        CONTINUE 
        !
        Kup = 0.0 
        traceK = 0.0 
        K(1,:) = (/ Cij(1,c1), Cij(2,c1), Cij(3,c1) /) 
        K(2,:) = (/ Cij(2,c1), Cij(4,c1), Cij(5,c1) /) 
        K(3,:) = (/ Cij(3,c1), Cij(5,c1), Cij(6,c1) /) 
        !
        DO ii = 1, 3 
         DO jj = 1, 3
          traceK = traceK + g_contr(ii,jj,c1)*K(ii,jj)  
          DO pp = 1, 3
           DO qq = 1, 3 
              Kup(ii,jj) = Kup(ii,jj) + g_contr(ii,pp,c1)*g_contr(jj,qq,c1)*K(pp,qq)  
           ENDDO 
          ENDDO 
         ENDDO
        ENDDO 
        K2 = 0.0
        DO ii = 1, 3 
         DO jj = 1, 3
            K2 = K2 + Kup(ii,jj)*K(ii,jj) 
         ENDDO
        ENDDO     
        !
        aux = (/ wGPN(iii), wGPN(jjj), wGPN(kkk) /) 
        wGP = PRODUCT( aux(1:nDim) ) 
        ! 
        check = check + wGP * ( R(c1) - K2 + traceK**2 )      
        !
      ENDDO
     ENDDO
    ENDDO 
    !
#endif     
    !
END SUBROUTINE CheckHamiltonianConstraint  

    
    