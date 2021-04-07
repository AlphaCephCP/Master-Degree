SUBROUTINE RunTimeAnalyseData 
    !-------------------------------------------------------------------------!
    USE MainVariables
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER, PARAMETER      :: nConstraints = 6 
    INTEGER                 :: i,j,k,iVar 
    INTEGER                 :: im1,im2,ip1,ip2,im3,ip3 
    REAL                    :: Constraints(nConstraints)
    REAL                    :: L1norm(nConstraints), L2norm(nConstraints), Linfnorm(nConstraints) 
    REAL                    :: Q(nVar), ux(nVar), uy(nVar), uz(nVar), gradQ(nVar,d) 
    CHARACTER(LEN=10)       :: OutName(50) 
    CHARACTER(LEN=200)      :: Filename, format_string
    INTEGER, PARAMETER      :: RunTimeAnalyseType = 1 
    !-------------------------------------------------------------------------!
    !
    L1norm(:)   = 0.
    L2norm(:)   = 0.
    Linfnorm(:) = 0.
    !
#ifdef CCZ4EINSTEIN 
    ! 
    DO k = 1, KMAX
     DO j = 1, JMAX
      DO i = 1, IMAX
        im1 = i-1
        ip1 = i+1
        im2 = i-2
        ip2 = i+2
        im3 = i-3
        ip3 = i+3
        IF(im1<=0) im1 = IMAX+im1 
        IF(im2<=0) im2 = IMAX+im2 
        IF(im3<=0) im3 = IMAX+im3 
        IF(ip1>=IMAX+1) ip1 = ip1-IMAX 
        IF(ip2>=IMAX+1) ip2 = ip2-IMAX 
        IF(ip3>=IMAX+1) ip3 = ip3-IMAX 
        Q = uh(:,i,j,k) 
        !ux = ( 8.0*uh(:,ip1,j,k) - 8.0*uh(:,im1,j,k) + uh(:,im2,j,k) - uh(:,ip2,j,k) )/( 12.0*dx(1) ) 
        ux = ( 9.0*uh(:,im2,j,k) - 45.0*uh(:,im1,j,k) + 45.0*uh(:,ip1,j,k) - 9.0*uh(:,ip2,j,k) - uh(:,im3,j,k) + uh(:,ip3,j,k)  ) / ( 60.0*dx(1) ) 
        !uxx = 1./dx(1)**2*( u(:,ip1,j,k) - 2*u(:,i,j,k) + u(:,im1,j,k) )
        gradQ(:,1) = ux 
        CALL ADMConstraints(Constraints,Q,gradQ)
        ! 
        L1norm(:)   = L1norm(:) + ABS(Constraints(:))    * PRODUCT(dx(1:nDim))
        L2norm(:)   = L2norm(:) + ABS(Constraints(:))**2 * PRODUCT(dx(1:nDim))
        !
        DO iVar = 1, nConstraints 
            Linfnorm(iVar) = MAX( Linfnorm(iVar), ABS(Constraints(iVar)) ) 
        ENDDO 
        !      
      ENDDO
     ENDDO 
    ENDDO
    ! Take square-root to obtain the L2-norm
    L2norm(:) = SQRT(L2norm(:))    
    !
    !
    OPEN(95,FILE='ADM-L1.dat',status='UNKNOWN',position='APPEND')  
    OPEN(96,FILE='ADM-L2.dat',status='UNKNOWN',position='APPEND')
    OPEN(97,FILE='ADM-Linf.dat',status='UNKNOWN',position='APPEND')
    WRITE(95,326) time, L1norm(:) 
    WRITE(96,326) time, L2norm(:)  
    WRITE(97,326) time, Linfnorm(:) 
    CLOSE(95)
    CLOSE(96)
    CLOSE(97)
    !
326 FORMAT(1x,7(E13.6,1x))    
    ! 
#endif    
    !    
END SUBROUTINE RunTimeAnalyseData 

    
SUBROUTINE ADMConstraints( Constraints, Q, gradQ )
   USE MainVariables, ONLY : EQN, nVar, d    
   IMPLICIT NONE
   INTENT(IN)  :: Q, gradQ 
   INTENT(OUT) :: Constraints
   INTEGER, PARAMETER :: nConstraints = 6 
   INTEGER :: i, ip, j, k, l, m, n, iErr, qq, ii, jj, kk, ll, mm, nn  
   REAL :: xGP(d), Constraints(nConstraints), Q(nVar), gradQ(nVar,d), gradQT(d,Nvar)  
   REAL :: traceK, R, phi, KK2
   REAL :: g_contr(3,3), g_cov(3,3), Ricci(3,3)
   REAL :: DD(3,3,3), Atilde(3,3), PP(3), GG(3), dP(3,3)
   REAL :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
   REAL :: dDD(3,3,3,3), Christoffel(3,3,3), ChristoffelNC(3,3,3), Id(3,3), dgup(3,3,3), Riemann(3,3,3,3), dChristoffel(3,3,3,3)
   REAL :: Christoffel_diff(3,3,3,3), Ham, Mom(3), dK(3,3,3), dAtilde(3,3,3), dtraceK(3), Qx(nVar), Qy(nVar), Qz(nVar)
   REAL :: dg_cov(3,3,3), g_covx(3,3), g_covy(3,3), g_covz(3,3), det
   REAL :: alpha, Aex(3,3), Kex(3,3), traceA, k0, dAex(3,3,3), dKex(3,3,3), Amix(3,3), Aup(3,3), Kmix(3,3), Kup(3,3)
   REAL :: ghat(3), theta, dtheta(3), dghat(3,3), AA(3), dAA(3,3), BB(3,3), dBB(3,3,3), dphi(3), dPP(3,3), beta(3) 
   REAL :: Christoffel_tilde(3,3,3), Gtilde(3), Christoffel_kind1(3,3,3), Z(3), Zup(3), Kupdown 

   Constraints(:) = 0.0

   Qx = gradQ(:,1) 
   Qy = gradQ(:,2) 
   Qz = gradQ(:,3) 
   gradQT = TRANSPOSE(gradQ) 

#if defined(CCZ4EINSTEIN) 

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
      
    alpha = EXP(Q(17)) 
    ! 
    K0    = Q(59)
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
    Aex = Aex - 1./3.*g_cov*traceA 
    !
    dAex(:,1,1) = gradQ(7,:) 
    dAex(:,1,2) = gradQ(8,:) 
    dAex(:,1,3) = gradQ(9,:) 
    dAex(:,2,1) = gradQ(8,:) 
    dAex(:,2,2) = gradQ(10,:) 
    dAex(:,2,3) = gradQ(11,:) 
    dAex(:,3,1) = gradQ(9,:) 
    dAex(:,3,2) = gradQ(11,:) 
    dAex(:,3,3) = gradQ(12,:) 
    !
    Amix = MATMUL(g_contr, Aex)
    Aup  = MATMUL(g_contr, TRANSPOSE(Amix)) 
    !
    Theta = Q(13)
    dTheta = gradQ(13,:) 
    ! 
    Ghat = (/ Q(14), Q(15), Q(16) /)
    dGhat(:,1) = gradQ(14,:)
    dGhat(:,2) = gradQ(15,:)
    dGhat(:,3) = gradQ(16,:)
    !
    AA    = (/ Q(24), Q(25), Q(26) /) 
    dAA(:,1) = gradQ(24,:) 
    dAA(:,2) = gradQ(25,:) 
    dAA(:,3) = gradQ(26,:) 
    !
    traceK = Q(54) 
    dtraceK = gradQ(54,:) 
    !
    phi   = EXP(Q(55)) 
    dphi  = gradQ(55,:) 
    PP    = Q(56:58) 
    dPP(:,1) = gradQ(56,:) 
    dPP(:,2) = gradQ(57,:) 
    dPP(:,3) = gradQ(58,:) 
    !
    beta = (/ Q(18), Q(19), Q(20) /) 
    BB(1,1) = Q(27) 
    BB(2,1) = Q(28) 
    BB(3,1) = Q(29) 
    BB(1,2) = Q(30) 
    BB(2,2) = Q(31) 
    BB(3,2) = Q(32) 
    BB(1,3) = Q(33) 
    BB(2,3) = Q(34) 
    BB(3,3) = Q(35) 
    !
    dBB(:,1,1) = gradQ(27,:) 
    dBB(:,2,1) = gradQ(28,:) 
    dBB(:,3,1) = gradQ(29,:) 
    dBB(:,1,2) = gradQ(30,:) 
    dBB(:,2,2) = gradQ(31,:) 
    dBB(:,3,2) = gradQ(32,:) 
    dBB(:,1,3) = gradQ(33,:) 
    dBB(:,2,3) = gradQ(34,:) 
    dBB(:,3,3) = gradQ(35,:) 
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
    dDD(:,1,1,1)=gradQ(36,:) 
    dDD(:,1,1,2)=gradQ(37,:) 
    dDD(:,1,1,3)=gradQ(38,:) 
    dDD(:,1,2,1)=gradQ(37,:) 
    dDD(:,1,2,2)=gradQ(39,:) 
    dDD(:,1,2,3)=gradQ(40,:)
    dDD(:,1,3,1)=gradQ(38,:) 
    dDD(:,1,3,2)=gradQ(40,:) 
    dDD(:,1,3,3)=gradQ(41,:)
    dDD(:,2,1,1)=gradQ(42,:) 
    dDD(:,2,1,2)=gradQ(43,:) 
    dDD(:,2,1,3)=gradQ(44,:) 
    dDD(:,2,2,1)=gradQ(43,:) 
    dDD(:,2,2,2)=gradQ(45,:) 
    dDD(:,2,2,3)=gradQ(46,:)
    dDD(:,2,3,1)=gradQ(44,:) 
    dDD(:,2,3,2)=gradQ(46,:) 
    dDD(:,2,3,3)=gradQ(47,:) 
    dDD(:,3,1,1)=gradQ(48,:) 
    dDD(:,3,1,2)=gradQ(49,:) 
    dDD(:,3,1,3)=gradQ(50,:) 
    dDD(:,3,2,1)=gradQ(49,:) 
    dDD(:,3,2,2)=gradQ(51,:) 
    dDD(:,3,2,3)=gradQ(52,:)
    dDD(:,3,3,1)=gradQ(50,:) 
    dDD(:,3,3,2)=gradQ(52,:) 
    dDD(:,3,3,3)=gradQ(53,:)
    !
    dgup = 0.0 
    DO k = 1, 3 
     DO m = 1, 3 
      DO l = 1, 3 
       DO n = 1, 3
        DO j = 1, 3 
           dgup(k,m,l) = dgup(k,m,l)-g_contr(m,n)*g_contr(j,l)*2*DD(k,n,j) 
        ENDDO
       ENDDO 
      ENDDO 
     ENDDO 
    ENDDO         
    !
    Kex  = Aex/phi**2 + 1./3.*traceK*g_cov/phi**2 
    Kmix = MATMUL( phi**2*g_contr, Kex  ) 
    Kup  = MATMUL( phi**2*g_contr, TRANSPOSE(Kmix) ) 
    !
    Christoffel_tilde = 0.0  
    Christoffel       = 0.0 
    Gtilde = 0.0 
    !
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
       !Christoffel_kind1(i,j,k) =  DD(i,j,k)+DD(j,i,k)-DD(k,i,j)     ! this definition does not work ! 
       Christoffel_kind1(i,j,k) = DD(k,i,j)+DD(j,i,k)-DD(i,j,k)      ! this definition seems to work ! 
       DO l = 1, 3
          Christoffel_tilde(i,j,k) = Christoffel_tilde(i,j,k) + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) 
          Christoffel(i,j,k)       = Christoffel(i,j,k)       + g_contr(k,l)*( DD(i,j,l)+DD(j,i,l)-DD(l,i,j) ) -g_contr(k,l)*( g_cov(j,l)*PP(i)+g_cov(i,l)*PP(j)-g_cov(i,j)*PP(l) ) 
          Gtilde(i)                = Gtilde(i)+2*g_contr(i,j)*g_contr(k,l)*DD(l,j,k) 
        ENDDO 
      ENDDO
     ENDDO
    ENDDO
    Z   = 0.5*MATMUL( g_cov, Ghat - Gtilde ) 
    Zup = MATMUL(phi**2*g_contr, Z) 
    !
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          dChristoffel(k,i,ip,m) = 0 
          DO l = 1, 3 
            dChristoffel(k,i,ip,m) = dChristoffel(k,i,ip,m) + g_contr(m,l)*(0.5*(dDD(k,i,ip,l)+dDD(i,k,ip,l))+0.5*(dDD(k,ip,i,l)+dDD(ip,k,i,l))-0.5*(dDD(k,l,i,ip)+dDD(l,k,i,ip)))         & 
                                                            - g_contr(m,l)*(g_cov(ip,l)*0.5*(dPP(k,i)+dPP(i,k))+g_cov(i,l)*0.5*(dPP(k,ip)+dPP(ip,k))-g_cov(i,ip)*0.5*(dPP(k,l)+dPP(l,k)))  & 
                                                             +dgup(k,m,l)*(DD(i,ip,l)+DD(ip,i,l)-DD(l,i,ip)) - dgup(k,m,l)*(g_cov(ip,l)*PP(i)+g_cov(i,l)*PP(ip)-g_cov(i,ip)*PP(l)) - g_contr(m,l)*( 2*DD(k,ip,l)*PP(i)+2*DD(k,i,l)*PP(ip)-2*DD(k,i,ip)*PP(l) ) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    Riemann = 0.0 
    DO i = 1, 3 
     DO ip = 1, 3 
      DO m = 1, 3 
       DO k = 1, 3
          Riemann(i,k,ip,m)    = dChristoffel(k,i,ip,m)-dChristoffel(ip,i,k,m)
          DO j = 1, 3
           Riemann(i,k,ip,m)    = Riemann(i,k,ip,m)    + Christoffel(i,ip,j)*Christoffel(j,k,m) - Christoffel(i,k,j)*Christoffel(j,ip,m) 
          ENDDO 
       ENDDO
      ENDDO
     ENDDO
    ENDDO    
    ! 
    Ricci = 0.0 
    DO m = 1, 3 
     DO n = 1, 3
      DO l = 1, 3    
         Ricci(m,n) = Ricci(m,n) + Riemann(m,l,n,l)  
      ENDDO
     ENDDO
    ENDDO    
    !
    R = SUM(phi**2*g_contr*Ricci) 
    !
    Kupdown = SUM(Kex*Kup) 
    Ham = R - KupDown + traceK**2 
    !
    dKex = 0.0
    DO j = 1, 3
     DO i = 1, 3
      DO k = 1, 3
          dKex(k,i,j) = 1.0/phi**2*( dAex(k,i,j) - 2.0*Aex(i,j)*PP(k) + 1./3.*dtraceK(k)*g_cov(i,j) + 2./3.*traceK*DD(k,i,j) - 2./3.*traceK*g_cov(i,j)*PP(k) ) 
      ENDDO
     ENDDO
    ENDDO 
    !
    Mom(:) = 0.0
    DO ii = 1, 3
        DO jj = 1, 3
            DO ll = 1, 3
                Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*(dKex(ll,ii,jj) - dKex(ii,jj,ll))
                DO mm = 1, 3
                    Mom(ii) = Mom(ii) + phi**2*g_contr(jj,ll)*( - Christoffel(jj,ll,mm)*Kex(mm,ii) + Christoffel(jj,ii,mm)*Kex(mm,ll)) 
                ENDDO
            ENDDO     
        ENDDO
    ENDDO   
   
    Constraints(1)   = Ham 
    Constraints(2:4) = Mom(1:3)
    Constraints(5)   = det - 1.0 
    Constraints(6)   = traceA 
    !
    CONTINUE
    !
#endif 

END SUBROUTINE ADMConstraints 
    
    
    
    
    