SUBROUTINE AnalyseDG 
    USE typesDef
    IMPLICIT NONE
#ifdef PARALLEL 
    INCLUDE 'mpif.h'  
#endif 
    INTEGER :: iElem, iVar, ii, jj, kk, iErr, iInfElem(nVar)  
    REAL    :: ldx(d), lx0(d), xGP(d), par0(nParam), aux(d) 
    REAL    :: u0(nVar), state(nVar), pstate(nVar), pexact(nVar)   
    REAL    :: L1norm(nVar), L2norm(nVar), Linfnorm(nVar), locError(nVar)
    REAL    :: MPIL1(nVar), MPIL2(nVar), MPILinf(nVar),rExc
    CHARACTER(LEN=10)  :: varName(nVar) 
    !
    IF(AnalyseType==0) THEN
        RETURN
    ENDIF  
    !
    L1norm = 0. 
    L2norm = 0. 
    Linfnorm = 0. 
    ! 
    ldx = dx 
    DO iElem = 1, nElem
      lx0 = x(:,tri(1,iElem))
      !      
#ifndef GRMHD      
      IF( SQRT(SUM((lx0+0.5*dx)**2))<ExcisionRadius ) THEN
          CYCLE
      ENDIF 
#else
      rExc = MAXVAL(ABS(x(:,tri(:,iElem))))
      IF( rExc.LT.ExcisionRadius ) THEN
          CYCLE
      ENDIF 
#endif
      !  
      DO kk = 1, nDOF(3) 
       DO jj = 1, nDOF(2)
        DO ii = 1, nDOF(1) 
            xGP = lx0 + (/ xiGPN(ii), xiGPN(jj), xiGPN(kk) /)*ldx
            SELECT CASE(AnalyseType)
            CASE(1) ! Initialcondition
                CALL InitialField(u0,par0,xGP,0.0)
                CALL PDECons2Prim(pexact,u0,iErr)
            CASE(2) ! Time-shifted Initialcondition
                CALL InitialField(u0,par0,xGP,time)
                CALL PDECons2Prim(pexact,u0,iErr)
            END SELECT 
            state(:) = uh(:,ii,jj,kk,iElem)
            CALL PDECons2Prim(pstate,state,iErr)

            locError(:) = ABS(pexact(:)-pstate(:))
            aux = (/ wGPN(ii), wGPN(jj), wGPN(kk) /)
            ! 
            L1norm(:)   = L1norm(:) + locError(:)    * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim)) 
            L2norm(:)   = L2norm(:) + locError(:)**2 * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim))   
            !
            DO iVar = 1, nVar
                IF(locError(iVar).GT.Linfnorm(iVar)) THEN
                    Linfnorm(iVar) = locError(iVar)
                    iInfElem(iVar) = iElem
                ENDIF
            ENDDO
            ! 
       ENDDO
      ENDDO
     ENDDO
     ! 
    ENDDO 
    !
#ifdef PARALLEL
    DO iVar = 1, nVar
        CALL MPI_ALLREDUCE(L1norm(iVar),  MPIL1(iVar),  1,MPI_AUTO_REAL,MPI_SUM,MPI_COMM_WORLD,mpiErr) 
        CALL MPI_ALLREDUCE(L2norm(iVar),  MPIL2(iVar),  1,MPI_AUTO_REAL,MPI_SUM,MPI_COMM_WORLD,mpiErr) 
        CALL MPI_ALLREDUCE(Linfnorm(iVar),MPILinf(iVar),1,MPI_AUTO_REAL,MPI_MAX,MPI_COMM_WORLD,mpiErr) 
        L1norm(iVar)   = MPIL1(iVar) 
        L2norm(iVar)   = MPIL2(iVar) 
        Linfnorm(iVar) = MPILinf(iVar) 
    ENDDO 
#endif 
    !
    L2norm = SQRT(L2norm) 
    ! 
    CALL PDEVarName(VarName) 
    IF(myrank==0) THEN
        DO iVar = 1, nVar                                                         
            WRITE(*,*) '| Error analysis of variable ', TRIM(VarName(iVar)) 
            WRITE(*,*) '| ======================================'     !
            WRITE(*,*) '|   L1_norm   : ', L1norm(iVar)               !
            WRITE(*,*) '|   L2_norm   : ', L2norm(iVar)               !
            WRITE(*,*) '|   Linf_norm : ', Linfnorm(iVar), iInfElem(iVar)
            WRITE(*,*) '| ======================================'     !
        END DO
    ENDIF 
    !
END SUBROUTINE AnalyseDG 

    
SUBROUTINE RunTimeAnalyseData 
    !-------------------------------------------------------------------------!
    USE typesDef
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif
    !-------------------------------------------------------------------------!
    ! Local variable declaration
#if defined(Z4EINSTEIN) || defined(Z4GRMHD)
    INTEGER, PARAMETER :: nConstraints = 4   
#elif defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)
    INTEGER, PARAMETER :: nConstraints = 6 
#else
    INTEGER, PARAMETER :: nConstraints = 0 
#endif    
    INTEGER                 :: iElem,ii,jj,kk,reflev, iii        ! Loop variables         !
    INTEGER                 :: iVar                  ! Variable to analyse    !
    INTEGER                 :: iIntGP                ! Index of internal GP   !
    INTEGER                 :: iInfElem(nVar)        ! Element of max. error  !
    INTEGER                 :: iErr,allocstat                  ! 
    INTEGER                 :: rse(N+1,N+1,N+1), c, nMonitor
    REAL                    :: xGP(d), LF, Constraints(nConstraints), lwh(nVar,N+1,N+1,N+1) 
    REAL                    :: lwx(nVar,N+1,N+1,N+1), lwy(nVar,N+1,N+1,N+1), lwz(nVar,N+1,N+1,N+1)
    REAL                    :: state(nVar), pstate(nVar), gradstate(nVar,d)    
    REAL                    :: locQ(nConstraints)        ! local error in GP      !
    REAL                    :: locC(nConstraints)  
    REAL                    :: L1norm(nConstraints)      ! L1norm                 !
    REAL                    :: L2norm(nConstraints)      ! L2norm                 !
    REAL                    :: Linfnorm(nConstraints) 
    REAL                    :: MPIL1norm(nConstraints)   ! Total L1norm (MPI)     !
    REAL                    :: MPIL2norm(nConstraints)   ! Total L2norm (MPI)     !
    REAL                    :: MPILinfnorm(nConstraints) ! Total L_inf norm  (MPI)!
    REAL                    :: lx0(d),ldx(d), aux(d), vel(3), B(3)
    CHARACTER(LEN=10)       :: OutName(50) 
    CHARACTER(LEN=200)      :: Filename, format_string
    INTEGER, PARAMETER      :: RunTimeAnalyseType = 1 
    !-------------------------------------------------------------------------!
    !
    L1norm(:)   = 0.
    L2norm(:)   = 0.
    Linfnorm(:) = 0.
    !
    DO iElem = 1, nElem
      ! 
      ldx = dx   
      lx0 = x(:,tri(1,iElem)) + 0.5*ldx   
#ifndef GRMHD
      IF( SQRT(SUM(lx0**2))<ExcisionRadius ) THEN
#else
      IF( MAXVAL(ABS(x(:,tri(:,iElem)))).LT.ExcisionRadius ) THEN
#endif
          ADM(:,:,:,:,iElem) = 0.0 
          CYCLE
      ENDIF 
      !
      DO kk = 1, nDOF(3) 
         DO jj = 1, nDOF(2)
            DO ii = 1, nDOF(1) 
               lwh(:,ii,jj,kk) = uh(:,ii,jj,kk,iElem) 
            ENDDO
         ENDDO
      ENDDO   
      ! Compute the x derivatives in each point   
      DO kk = 1, nDOF(3)  
         DO jj = 1, nDOF(2) 
            lwx(:,:,jj,kk) = MATMUL( lwh(:,:,jj,kk), TRANSPOSE(dudx) )/dx(1)  
         ENDDO
      ENDDO 
      IF(nDim>=2) THEN
          ! Compute the y derivatives in each point   
          DO kk = 1, nDOF(3)  
             DO jj = 1, nDOF(1) 
                lwy(:,jj,:,kk) = MATMUL( lwh(:,jj,:,kk), TRANSPOSE(dudx) )/dx(2)  
             ENDDO
          ENDDO 
      ELSE
          lwy = 0.0
      ENDIF 
      IF(nDim>=3) THEN
          ! Compute the z derivatives in each point 
          DO kk = 1, nDOF(2)  
            DO jj = 1, nDOF(1) 
                lwz(:,jj,kk,:) = MATMUL( lwh(:,jj,kk,:), TRANSPOSE(dudx) )/dx(3)  
            ENDDO
          ENDDO         
      ELSE
          lwz = 0.0 
      ENDIF     
      ! 
      DO kk = 1, nDOF(3) 
         DO jj = 1, nDOF(2)
            DO ii = 1, nDOF(1) 
               !
               xGP = x(:,tri(1,iElem)) + (/ xiGPN(ii), xiGPN(jj), xiGPN(kk) /)*dx
               !IF( SQRT(SUM(xGP**2))<ExcisionRadius ) THEN 
               !     ADM(:,ii,jj,kk,iElem) = 0.0
               !     CYCLE
               !ENDIF 
               !
               state(:)   = uh(:,ii,jj,kk,iElem)
               CALL PDECons2Prim(pstate,state,iErr)
               !
               locQ(:) = 0.0
               !
               SELECT CASE(RunTimeAnalyseType)
               CASE(1) 
                    nMonitor  = nConstraints 
                    gradstate(:,1) = lwx(:,ii,jj,kk) 
                    gradstate(:,2) = lwy(:,ii,jj,kk) 
                    gradstate(:,3) = lwz(:,ii,jj,kk)  
                    CALL ADMConstraints(Constraints,state,gradstate)
                    locQ    = 0.
                    locQ(:) = Constraints(:)
                    ADM(:,ii,jj,kk,iElem) = Constraints(:) 
               CASE DEFAULT
                    IF (myrank == 0) THEN                                                              
                        WRITE(*,*)'ERROR: RunTime Analysis of Type ',RunTimeAnalyseType,' unknown!'
                    ENDIF
                    STOP 
               END SELECT 
               !
               aux = (/ wGPN(ii), wGPN(jj), wGPN(kk) /)
               ! 
               locQ(:)     = ABS(locQ(:))
               L1norm(:)   = L1norm(:) + locQ(:)    * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim)) 
               L2norm(:)   = L2norm(:) + locQ(:)**2 * PRODUCT(aux(1:nDim))*PRODUCT(ldx(1:nDim))   
               !
               DO iVar = 1, nMonitor
                  IF(locQ(iVar).GT.Linfnorm(iVar)) THEN
                     Linfnorm(iVar) = locQ(iVar)
                     iInfElem(iVar) = iElem
                  ENDIF
               ENDDO               
               
            ENDDO
         ENDDO
      ENDDO
      ! 
    ENDDO
    !
    !
#ifdef PARALLEL
    DO iVar = 1, nMonitor                                                         
       CALL MPI_REDUCE(L1norm(iVar),  MPIL1norm(iVar),   1, MPI_AUTO_REAL, MPI_SUM, 0, MPI_COMM_WORLD, iErr) 
       CALL MPI_REDUCE(L2norm(iVar),  MPIL2norm(iVar),   1, MPI_AUTO_REAL, MPI_SUM, 0, MPI_COMM_WORLD, iErr) 
       CALL MPI_REDUCE(Linfnorm(iVar),MPILinfnorm(iVar), 1, MPI_AUTO_REAL, MPI_MAX, 0, MPI_COMM_WORLD, iErr)
    ENDDO 
    L1norm   = MPIL1norm 
    L2norm   = MPIL2norm 
    Linfnorm = MPILinfnorm
#endif 
    ! 
    ! Take square-root to obtain the L2-norm
    L2norm(:) = SQRT(L2norm(:))    
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD) || defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)     
    !
    IF(myrank.EQ.0) THEN 
       OPEN(95,FILE='ADM-L1.dat',status='UNKNOWN',position='APPEND')  
       OPEN(96,FILE='ADM-L2.dat',status='UNKNOWN',position='APPEND')
       OPEN(97,FILE='ADM-Linf.dat',status='UNKNOWN',position='APPEND')
       WRITE(95,326) time, L1norm(:) 
       WRITE(96,326) time, L2norm(:)  
       WRITE(97,326) time, Linfnorm(:) 
       CLOSE(95)
       CLOSE(96)
       CLOSE(97)
    ENDIF          
    !
#endif     
    !
#if defined(Z4EINSTEIN) || defined(Z4GRMHD)
326 FORMAT(1x,5(E13.6,1x))    
#endif    
#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)
326 FORMAT(1x,7(E13.6,1x))    
#endif    
    !    
END SUBROUTINE RunTimeAnalyseData 

SUBROUTINE ADMConstraints( Constraints, Q, gradQ )
   USE typesDef, ONLY : EQN, nVar, d    
   IMPLICIT NONE
   INTENT(IN)  :: Q, gradQ 
   INTENT(OUT) :: Constraints
#if defined(Z4EINSTEIN) || defined(Z4GRMHD)
    INTEGER, PARAMETER :: nConstraints = 4 
#elif defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)
    INTEGER, PARAMETER :: nConstraints = 6 
#else
    INTEGER, PARAMETER :: nConstraints = 0 
#endif    
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

#if defined(Z4EINSTEIN) || defined(Z4GRMHD) 
   
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
   !
   g_covx( 1, 1:3) = (/ Qx(1),  Qx(2),  Qx(3)   /) 
   g_covx( 2, 1:3) = (/ Qx(2),  Qx(4),  Qx(5)   /) 
   g_covx( 3, 1:3) = (/ Qx(3),  Qx(5),  Qx(6)   /) 
   !
   g_covy( 1, 1:3) = (/ Qy(1),  Qy(2),  Qy(3)   /) 
   g_covy( 2, 1:3) = (/ Qy(2),  Qy(4),  Qy(5)   /) 
   g_covy( 3, 1:3) = (/ Qy(3),  Qy(5),  Qy(6)   /) 
   !
   g_covz( 1, 1:3) = (/ Qz(1),  Qz(2),  Qz(3)   /) 
   g_covz( 2, 1:3) = (/ Qz(2),  Qz(4),  Qz(5)   /) 
   g_covz( 3, 1:3) = (/ Qz(3),  Qz(5),  Qz(6)   /) 
   !
   dg_cov(1,:,:) = g_covx 
   dg_cov(2,:,:) = g_covy 
   dg_cov(3,:,:) = g_covz 
   !
   Kex = RESHAPE( (/ Q(7), Q(8), Q(9), Q(8), Q(10), Q(11), Q(9), Q(11), Q(12) /), (/3,3/) )   
   !
   DD(1,:,:) = RESHAPE ( (/Q(36), Q(37), Q(38), Q(37), Q(39), Q(40), Q(38), Q(40), Q(41)  /)  , (/3,3/) )
   DD(2,:,:) = RESHAPE ( (/Q(42), Q(43), Q(44), Q(43), Q(45), Q(46), Q(44), Q(46), Q(47)  /)  , (/3,3/) )
   DD(3,:,:) = RESHAPE ( (/Q(48), Q(49), Q(50), Q(49), Q(51), Q(52), Q(50), Q(52), Q(53)  /)  , (/3,3/) )
   
   ! Initialize the necessary derivatives
   
   dK(1,:,:) = RESHAPE( (/ gradQT(1,7), gradQT(1,8), gradQT(1,9), gradQT(1,8), gradQT(1,10), gradQT(1,11), gradQT(1,9), gradQT(1,11), gradQT(1,12) /), (/3,3/) )   
   dK(2,:,:) = RESHAPE( (/ gradQT(2,7), gradQT(2,8), gradQT(2,9), gradQT(2,8), gradQT(2,10), gradQT(2,11), gradQT(2,9), gradQT(2,11), gradQT(2,12) /), (/3,3/) )   
   dK(3,:,:) = RESHAPE( (/ gradQT(3,7), gradQT(3,8), gradQT(3,9), gradQT(3,8), gradQT(3,10), gradQT(3,11), gradQT(3,9), gradQT(3,11), gradQT(3,12) /), (/3,3/) ) 
   !
   dDD(1,1,:,:)   = RESHAPE ( (/gradQT(1,36), gradQT(1,37), gradQT(1,38), gradQT(1,37), gradQT(1,39), gradQT(1,40), gradQT(1,38), gradQT(1,40), gradQT(1,41)  /)  , (/3,3/) )
   dDD(2,1,:,:)   = RESHAPE ( (/gradQT(2,36), gradQT(2,37), gradQT(2,38), gradQT(2,37), gradQT(2,39), gradQT(2,40), gradQT(2,38), gradQT(2,40), gradQT(2,41)  /)  , (/3,3/) )
   dDD(3,1,:,:)   = RESHAPE ( (/gradQT(3,36), gradQT(3,37), gradQT(3,38), gradQT(3,37), gradQT(3,39), gradQT(3,40), gradQT(3,38), gradQT(3,40), gradQT(3,41)  /)  , (/3,3/) )   
   
   dDD(1,2,:,:)   = RESHAPE ( (/gradQT(1,42), gradQT(1,43), gradQT(1,44), gradQT(1,43), gradQT(1,45), gradQT(1,46), gradQT(1,44), gradQT(1,46), gradQT(1,47)  /)  , (/3,3/) )
   dDD(2,2,:,:)   = RESHAPE ( (/gradQT(2,42), gradQT(2,43), gradQT(2,44), gradQT(2,43), gradQT(2,45), gradQT(2,46), gradQT(2,44), gradQT(2,46), gradQT(2,47)  /)  , (/3,3/) )
   dDD(3,2,:,:)   = RESHAPE ( (/gradQT(3,42), gradQT(3,43), gradQT(3,44), gradQT(3,43), gradQT(3,45), gradQT(3,46), gradQT(3,44), gradQT(3,46), gradQT(3,47)  /)  , (/3,3/) )

   dDD(1,3,:,:)   = RESHAPE ( (/gradQT(1,48), gradQT(1,49), gradQT(1,50), gradQT(1,49), gradQT(1,51), gradQT(1,52), gradQT(1,50), gradQT(1,52), gradQT(1,53)  /)  , (/3,3/) )
   dDD(2,3,:,:)   = RESHAPE ( (/gradQT(2,48), gradQT(2,49), gradQT(2,50), gradQT(2,49), gradQT(2,51), gradQT(2,52), gradQT(2,50), gradQT(2,52), gradQT(2,53)  /)  , (/3,3/) )
   dDD(3,3,:,:)   = RESHAPE ( (/gradQT(3,48), gradQT(3,49), gradQT(3,50), gradQT(3,49), gradQT(3,51), gradQT(3,52), gradQT(3,50), gradQT(3,52), gradQT(3,53)  /)  , (/3,3/) )  
   !
   
   !!! Compute Christoffel in the original metric
   Christoffel   = 0.0 
   ChristoffelNC = 0.0 
   DO ii = 1, 3 
    DO jj = 1, 3  
     DO kk = 1, 3 
       DO ll = 1, 3 
         Christoffel(ii,jj,kk)   = Christoffel(ii,jj,kk)   + g_contr(kk,ll)*(DD(ii,jj,ll) + DD(jj,ii,ll) - DD(ll,ii,jj))
         ChristoffelNC(ii,jj,kk) = ChristoffelNC(ii,jj,kk) + g_contr(kk,ll)/2*(dg_cov(ii,jj,ll) + dg_cov(jj,ii,ll) - dg_cov(ll,ii,jj))
       ENDDO
     ENDDO
    ENDDO
   ENDDO   
   !
   !!! Compute dgup via Christoffel
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
   
   Riemann = 0.0 
   DO ii = 1, 3 
    DO qq = 1, 3 
     DO mm = 1, 3 
      DO kk = 1, 3         
       Riemann(ii,kk,qq,mm) = dChristoffel(kk,ii,qq,mm) - dChristoffel(qq,ii,kk,mm) 
       DO jj = 1, 3 
          !Riemann(ii,kk,qq,mm) = Riemann(ii,kk,qq,mm) + 0.5*( Christoffel(ii,qq,jj)*ChristoffelNC(jj,kk,mm) - Christoffel(ii,kk,jj)*ChristoffelNC(jj,qq,mm) )                     & 
          !                                            + 0.5*( ChristoffelNC(ii,qq,jj)*Christoffel(jj,kk,mm) - ChristoffelNC(ii,kk,jj)*Christoffel(jj,qq,mm) ) 
          Riemann(ii,kk,qq,mm) = Riemann(ii,kk,qq,mm) + Christoffel(ii,qq,jj)*Christoffel(jj,kk,mm) - Christoffel(ii,kk,jj)*Christoffel(jj,qq,mm) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
   ENDDO
   !
   Ricci = 0.0 
   DO mm = 1, 3
    DO nn = 1, 3 
     DO ll = 1, 3 
        Ricci(mm,nn) = Ricci(mm,nn) + Riemann(mm,ll,nn,ll)  
     ENDDO
    ENDDO
   ENDDO   
   
   ! Compute the 3-Ricci scalar
   R = 0.0
   DO ii = 1, 3 
       DO jj = 1, 3
           R = R + g_contr(ii,jj)*Ricci(ii,jj)
       ENDDO
   ENDDO   

   ! Compute the trace of K
   traceK = 0.0
   DO ii = 1, 3 
       DO jj = 1, 3
           traceK = traceK + g_contr(ii,jj)*Kex(ii,jj)
       ENDDO
   ENDDO   
   
   ! Compute K*K
   KK2 = 0.0
   DO ii = 1, 3 
       DO jj = 1, 3
           DO ll = 1, 3
               DO mm = 1, 3
                   KK2 = KK2 + g_contr(ii,ll)*g_contr(jj,mm)*Kex(ii,jj)*Kex(ll,mm)
               ENDDO
           ENDDO     
       ENDDO
   ENDDO
   ! 
   !
   ! Hamiltonian Constraint 
   Ham = R - KK2 + traceK**2
   
   ! Momentum Constrains
   Mom(:) = 0.0
   DO ii = 1, 3
      DO jj = 1, 3
         DO ll = 1, 3
             Mom(ii) = Mom(ii) + g_contr(jj,ll)*(dK(ll,ii,jj) - dK(ii,jj,ll))
             DO mm = 1, 3
                !Mom(ii) = Mom(ii) + g_contr(jj,ll)*( - ChristoffelNC(jj,ll,mm)*K(mm,ii) + ChristoffelNC(jj,ii,mm)*K(mm,ll)) 
                Mom(ii) = Mom(ii) + g_contr(jj,ll)*( - Christoffel(jj,ll,mm)*Kex(mm,ii) + Christoffel(jj,ii,mm)*Kex(mm,ll)) 
             ENDDO
         ENDDO     
      ENDDO
   ENDDO   
   
   Constraints(1)   = Ham 
   Constraints(2:4) = Mom(1:3)
   

   CONTINUE

#endif 

#if defined(CCZ4EINSTEIN) || defined(CCZ4GRMHD) || defined(BSSZ4EINSTEIN) || defined(BSSZ4GRMHD)

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
    
    
    