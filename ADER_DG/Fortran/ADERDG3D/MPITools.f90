SUBROUTINE ComputeCPUDistribution(CPUDistribution,VMAX,nCPUx)  
    USE typesDef, ONLY : d, nElem, idxe, idxn, nCPU    
    IMPLICIT NONE 
    ! Argument list  
    INTEGER :: CPUDistribution(nElem)       ! The CPU distribution 
    INTEGER :: VMAX(d), nCPUx(d)            ! number of cells and CPUs per dimension 
    INTEGER :: myrank                       ! myrank 
    INTENT(IN) :: VMAX, nCPUx 
    INTENT(OUT) :: CPUDistribution 
    ! Local variables 
    INTEGER :: i, j, k, iElem, rank_count 
    REAL    :: DCPU(d) 
    INTEGER, ALLOCATABLE :: CPUTable(:,:,:) 
    INTEGER  :: iCPU(d)  
    !
    IF(nCPU==1) THEN
        CPUDistribution(:) = 0
        RETURN 
    ENDIF
    !
    DO i = 1, d
        DCPU(i) = REAL(VMAX(i))/REAL(nCPUx(i))          ! compute the number of elements per CPU per dimension 
    ENDDO 
    ALLOCATE( CPUTable(nCPUx(1),nCPUx(2),nCPUx(3)) )  
    rank_count = 0 
    DO k = 1, nCPUx(3) 
        DO j = 1, nCPUx(2) 
            DO i = 1, nCPUx(1)
                CPUTable(i,j,k) = rank_count
                rank_count = rank_count + 1 
            ENDDO
        ENDDO
    ENDDO
    DO k = 1, VMAX(3)
        DO j = 1, VMAX(2)
            DO i = 1, VMAX(1)
                iCPU(:) = MIN( CEILING( REAL( (/ i, j, k /) )/REAL(DCPU(:)) ), nCPUx ) 
                iElem = idxe(i,j,k) 
                CPUDistribution(iElem) = CPUTable( iCPU(1), iCPU(2), iCPU(3) )  
            ENDDO
        ENDDO
    ENDDO  
    ! 
    DEALLOCATE(CPUTable) 
    ! 
END SUBROUTINE ComputeCPUDistribution

 SUBROUTINE MPIExtractMesh   
    USE typesDef  
    IMPLICIT NONE 
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif    
    ! Local variables 
    INTEGER                              :: i, j, k, ii, jj, kk, iElem, jElem, iNode, iFace  
    INTEGER                              :: iCPU, jCPU, msglength, index   
    INTEGER                              :: nNewElem, nNewNode 
    INTEGER, ALLOCATABLE                 :: CPUTable(:,:,:) 
    INTEGER, ALLOCATABLE                 :: Global2LocalNode(:) 
    INTEGER, ALLOCATABLE                 :: Local2GlobalNode(:)
    INTEGER, POINTER                     :: trinew(:,:), neighbornew(:,:,:,:)  
    REAL, POINTER                        :: xnew(:,:)  
    INTEGER, ALLOCATABLE                 :: nSendElem(:), nRecvElem(:), cumsum(:)  
    INTEGER                              :: OppositeFace(6) 
    TYPE(tIntegerMessage), ALLOCATABLE   :: send_imessage(:)
    TYPE(tIntegerMessage), ALLOCATABLE   :: recv_imessage(:)
    ! 
    IF(nCPU==1) THEN        
        ALLOCATE( Global2LocalElem(nElem), Local2GlobalElem(nElem) ) 
        ALLOCATE( LocalCPUDistribution(nElem) )         
        DO i = 1, nElem
            Global2LocalElem(i) = i 
            Local2GlobalElem(i) = i 
            LocalCPUDistribution(i) = CPUDistribution(i) 
        ENDDO 
        RETURN 
    ENDIF
    !
    OppositeFace = (/ 2, 1, 4, 3, 6, 5 /) 
    !
#ifdef PARALLEL 
    ! 
    ALLOCATE( Global2LocalElem(nElem), Global2LocalNode(nNode) ) 
    ALLOCATE( Local2GlobalElem(nElem), Local2GlobalNode(nNode) ) 
    Global2LocalElem = 0 
    Global2LocalNode = 0 
    Local2GlobalElem = 0 
    Local2GlobalNode = 0 
    nMPIElem  = 0 
    nNewElem  = 0 
    nNewNode  = 0 
    DO i = 1, nElem
        IF(CPUDistribution(i)==myrank) THEN
            nNewElem = nNewElem + 1 
            Global2LocalElem(i) = nNewElem
            Local2GlobalElem(nNewElem) = i  
            DO k = 1, nVtx
                iNode = tri(k,i) 
                IF(Global2LocalNode(iNode)==0) THEN
                    nNewNode = nNewNode + 1 
                    Global2LocalNode(iNode) = nNewNode
                    Local2GlobalNode(nNewNode) = iNode 
                ENDIF
            ENDDO  
            ! 
            DO  iFace = 1, 2*nDim 
                ii = Face2Neigh(1,iFace) 
                jj = Face2Neigh(2,iFace) 
                kk = Face2Neigh(3,iFace)
                j = neighbor(ii,jj,kk,i) 
                IF(CPUDistribution(j).NE.myrank) THEN
                    IF(Global2LocalElem(j)==0) THEN
                        nMPIElem    = nMPIElem + 1 
                        Global2LocalElem(j) = -nMPIElem    ! use negative element numbers for the stuff we need to receive via MPI 
                    ENDIF  
                ENDIF 
            ENDDO 
        ENDIF 
    ENDDO
    !
    !PRINT *, ' Info: ', myrank, nNewElem, nNewNode, nMPIElem 
    !
    ! Now compute the new neighbor, tri and x, and then destroy the old mesh
    !
    ALLOCATE( trinew(nVtx, nNewElem)                ) 
    ALLOCATE( neighbornew(-1:1,-1:1,-1:1, nNewElem) ) 
    ALLOCATE( xnew(d, nNewNode)                     )
    ALLOCATE( LocalCPUDistribution(-nMPIElem:nNewElem) ) 
    LocalCPUDistribution = myrank     
    DO iElem = 1, nElem
        IF( Global2LocalElem(iElem).NE.0 ) THEN
            LocalCPUDistribution(Global2LocalElem(iElem)) = CPUDistribution(iElem) 
        ENDIF        
    ENDDO    
    ! 
    DO iElem = 1, nNewElem
        i = Local2GlobalElem(iElem) 
        trinew(:,iElem) = Global2LocalNode(tri(:,i)) 
        DO kk = -1, 1
         DO jj = -1, 1
          DO ii = -1, 1 
            neighbornew(ii,jj,kk,iElem) = Global2LocalElem(neighbor(ii,jj,kk,i)) 
          ENDDO
         ENDDO
        ENDDO 
    ENDDO  
    DO iNode = 1, nNewNode
        i = Local2GlobalNode(iNode) 
        xnew(:,iNode) = x(:,i) 
    ENDDO 
    !
    ALLOCATE( nSendElem(0:nCPU-1), nRecvElem(0:nCPU-1) ) 
    nSendElem = 0 
    nRecvElem = 0 
    ! Count how many elements must be received based on face neighbors  
    DO iElem = 1, nElem  
        IF(CPUDistribution(iElem)==myrank) THEN 
            DO iFace = 1, 2*nDim   
                ii = Face2Neigh(1,iFace) 
                jj = Face2Neigh(2,iFace) 
                kk = Face2Neigh(3,iFace)
                jElem = neighbor(ii,jj,kk,iElem) 
                jCPU  = CPUDistribution(jElem)
                ! Neighbor is not myrank and is inside a given CPU number jCPU 
                IF(jCPU.NE.myrank) THEN
                    nRecvElem(jCPU) = nRecvElem(jCPU)+1     ! myrank must communicate with this CPU 
                ENDIF
            ENDDO
        ENDIF 
    ENDDO
    ! Compute a small vector that maps the global CPU number to the local number in the list of CPUs with which we need to communicate 
    nRecvCPU = 0 
    DO iCPU=0,nCPU-1
        IF(nRecvElem(iCPU).GT.0) THEN
            nRecvCPU = nRecvCPU + 1 
        ENDIF
    ENDDO 
    ALLOCATE( RecvCPU(nRecvCPU), invRecvCPU(0:nCPU-1) ) 
    nRecvCPU = 0 
    DO iCPU=0,nCPU-1
        IF(nRecvElem(iCPU).GT.0) THEN
            nRecvCPU = nRecvCPU + 1 
            RecvCPU(nRecvCPU) = iCPU 
            invRecvCPU(iCPU) = nRecvCPU
        ENDIF
    ENDDO 
    !
    ! Send this information to the neighbor CPUs  
    !
    ALLOCATE( send_imessage(0:nCPU-1) )
    ALLOCATE( recv_imessage(0:nCPU-1) )
    ALLOCATE( send_request(0:nCPU-1) )
    ALLOCATE( recv_request(0:nCPU-1) )
    ALLOCATE( send_status_list(MPI_STATUS_SIZE,0:nCPU-1) )
    ALLOCATE( recv_status_list(MPI_STATUS_SIZE,0:nCPU-1) )
    !
    DO iCPU = 0, nCPU-1
        MsgLength = 1
        ALLOCATE( send_imessage(iCPU)%Content(MsgLength) )
        ALLOCATE( recv_imessage(iCPU)%Content(MsgLength) )
        send_imessage(iCPU)%Content(1) = nRecvElem(iCPU) 
        ! Post a send request to neighbor CPU
        CALL MPI_ISEND(send_imessage(iCPU)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI_COMM_WORLD,send_request(iCPU), mpiErr)
        ! Post a receive request from neighbor CPU
        CALL MPI_IRECV(recv_imessage(iCPU)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI_COMM_WORLD, recv_request(iCPU), mpiErr)
    ENDDO
    ! Wait until all communication has finished
    CALL MPI_WAITALL(nCPU,send_request,send_status_list,mpiErr)
    CALL MPI_WAITALL(nCPU,recv_request,recv_status_list,mpiErr)
    !
    ! Decode messages into datastructure and deallocate messages
    !
    nSendCPU = 0 
    DO iCPU = 0, nCPU-1
        nSendElem(iCPU) = recv_imessage(iCPU)%Content(1)
        IF(nSendElem(iCPU).GT.0) THEN
            nSendCPU = nSendCPU+1
        ENDIF 
    ENDDO
    DO iCPU = 0, nCPU-1
        DEALLOCATE( send_imessage(iCPU)%Content )
        DEALLOCATE( recv_imessage(iCPU)%Content )
    ENDDO
    DEALLOCATE( send_imessage,recv_imessage,send_request,recv_request,send_status_list,recv_status_list )  
    !
    ! debug
    !PRINT *, myrank, 'nRecvElem = ', nRecvElem(:)
    !PRINT *, myrank, 'nSendElem = ', nSendElem(:)
    !
    ALLOCATE(SendCPU(nSendCPU),invSendCPU(0:nCPU-1)) 
    ALLOCATE(SendElem(nSendCPU)) 
    ! compute a small vector of CPUs to whom I have to send some data.... 
    nSendCPU = 0   
    DO iCPU = 0, nCPU-1 
        IF(nSendElem(iCPU).GT.0) THEN
            nSendCPU = nSendCPU + 1 
            SendCPU(nSendCPU) = iCPU          ! local list of CPUs with whom we are exchanging. iCPU is the REAL CPU number needed for the communication 
            invSendCPU(iCPU) = nSendCPU 
            SendElem(nSendCPU)%nElements = nSendElem(iCPU)  ! 
            ALLOCATE( SendElem(nSendCPU)%elements(  SendElem(nSendCPU)%nElements ) ) 
            ALLOCATE( SendElem(nSendCPU)%FaceIndex( SendElem(nSendCPU)%nElements ) ) 
        ENDIF
    ENDDO
    !
    ! Allocate the necessary space for the elements that we will need to receive during time stepping 
    ALLOCATE( RecvElem(nRecvCPU) ) 
    DO i = 1, nRecvCPU 
        iCPU = RecvCPU(i) 
        RecvElem(i)%nElements = nRecvElem(iCPU) 
        ALLOCATE( RecvElem(i)%Elements(  nRecvElem(iCPU) ) ) 
        ALLOCATE( RecvElem(i)%FaceIndex( nRecvElem(iCPU) ) ) 
    ENDDO
    ! Now put the data into the lists 
    nRecvElem = 0 
    DO iElem = 1, nElem  
        IF(CPUDistribution(iElem)==myrank) THEN
            DO iFace = 1, 2*nDim   
                ii = Face2Neigh(1,iFace) 
                jj = Face2Neigh(2,iFace) 
                kk = Face2Neigh(3,iFace)
                jElem = neighbor(ii,jj,kk,iElem) 
                jCPU  = CPUDistribution(jElem)
                ! Neighbor is not myrank and is inside a given CPU number jCPU 
                IF(jCPU.NE.myrank) THEN
                    nRecvElem(jCPU) = nRecvElem(jCPU)+1     ! myrank must communicate with this CPU 
                    j = invRecvCPU(jCPU) 
                    RecvElem(j)%Elements(nRecvElem(jCPU)) = jElem 
                    RecvElem(j)%FaceIndex(nRecvElem(jCPU)) = OppositeFace(iFace) 
                ENDIF
            ENDDO
        ENDIF 
    ENDDO
    !
    ALLOCATE( send_imessage(nRecvCPU) )
    ALLOCATE( recv_imessage(nSendCPU) )
    ALLOCATE( send_request(nRecvCPU) )
    ALLOCATE( recv_request(nSendCPU) )
    ALLOCATE( send_status_list(MPI_STATUS_SIZE,nRecvCPU) )
    ALLOCATE( recv_status_list(MPI_STATUS_SIZE,nSendCPU) )
    !
    ! Now exchange element numbers which are to be sent to the neighbors
    !
    DO i = 1, nRecvCPU
        iCPU = RecvCPU(i)    ! i is "local" number in the small communication vector, iCPU is the REAL CPU number 
        MsgLength = RecvElem(i)%nElements*2  
        ALLOCATE( send_imessage(i)%Content(MsgLength) ) 
        send_imessage(i)%Content(1:RecvElem(i)%nElements) = RecvElem(i)%elements(1:RecvElem(i)%nElements)                           ! myrank must SEND to iCPU WHICH elements are needed (i.e. that have to be RECEIVED in further timesteps) 
        send_imessage(i)%Content(RecvElem(i)%nElements+1:2*RecvElem(i)%nElements) = RecvElem(i)%FaceIndex(1:RecvElem(i)%nElements)  ! together with the local face index in the neighbor that is needed  
        ! debug 
        !WRITE(*,*) myrank, iCPU, RecvElem(i)%nElements, 'RE=', RecvElem(i)%elements(:) 
        RecvElem(i)%elements(:) = Global2LocalElem( RecvElem(i)%elements(:) )       ! now convert the global element numbers to CPU-local element numbers 
        ! Post a send request to neighbor CPU
        CALL MPI_ISEND(send_imessage(i)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI_COMM_WORLD, send_request(i), mpiErr)
    ENDDO
    !
    DO i = 1, nSendCPU
        iCPU = SendCPU(i)
        MsgLength = SendElem(i)%nElements*2   
        ALLOCATE( recv_imessage(i)%Content(MsgLength) )
        ! Post a receive request from neighbor CPU
        CALL MPI_IRECV(recv_imessage(i)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI_COMM_WORLD, recv_request(i), mpiErr)
    ENDDO
    ! Wait until all communication has finished
    CALL MPI_WAITALL(nRecvCPU,send_request,send_status_list,mpiErr)
    CALL MPI_WAITALL(nSendCPU,recv_request,recv_status_list,mpiErr)
    !
    ! Decode messages into datastructure and deallocate messages
    !
    DO i = 1, nSendCPU
        iCPU = SendCPU(i)
        SendElem(i)%elements(1:SendElem(i)%nElements)  = recv_imessage(i)%Content(1:SendElem(i)%nElements)
        SendElem(i)%FaceIndex(1:SendElem(i)%nElements) = recv_imessage(i)%Content(SendElem(i)%nElements+1:2*SendElem(i)%nElements)
        ! debug 
        !WRITE(*,*) myrank, iCPU, SendElem(i)%nElements, 'SE=', SendElem(i)%elements(:) 
        DO j = 1, SendElem(i)%nElements
            SendElem(i)%elements(j) = Global2LocalElem( SendElem(i)%elements(j) )  ! also here, we convert the global element numbers to CPU-local element numbers 
        ENDDO
    ENDDO
    !
    DO i = 1, nSendCPU 
        DEALLOCATE( recv_imessage(i)%Content )
    ENDDO
    !
    DEALLOCATE(send_imessage,recv_imessage,send_request,recv_request,send_status_list,recv_status_list)
    !
    ! Update point sources
    !
    DO i = 1, nPointSource
        ! Locate the sources on the Cartesian :-) grid
        PointSrc(i)%iElem = MAX(0, Global2LocalElem( PointSrc(i)%iElem ) ) 
    ENDDO     
    !
    DEALLOCATE( Global2LocalNode ) 
    DEALLOCATE( Local2GlobalNode ) 
    DEALLOCATE( x, tri, neighbor ) 
    !
    ! Now our mesh points only to the part contained in myrank 
    ! 
    x         => xnew 
    tri       => trinew
    neighbor  => neighbornew 
    nElem     = nNewElem
    nNode     = nNewNode  
    !
    PRINT *, ' MPIExtractMesh done for myrank = ', myrank
    CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr) 
    !
#endif 
    !
END SUBROUTINE MPIExtractMesh 

SUBROUTINE MPIExchangeData    
  USE typesDef   
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
#ifdef PARALLEL 
  INCLUDE 'mpif.h'  
  INTEGER                             :: i, j, ii, jj, iFace, iCPU, msglength, iElem, iVar, iDOF, counter 
  ! 
  ALLOCATE( send_message(nSendCPU) )
  ALLOCATE( recv_message(nRecvCPU) )
  ALLOCATE( send_request(nSendCPU) )
  ALLOCATE( recv_request(nRecvCPU) )
  ALLOCATE( send_status_list(MPI_STATUS_SIZE,nSendCPU) )
  ALLOCATE( recv_status_list(MPI_STATUS_SIZE,nRecvCPU) )
  !
  ! Now exchange cell averages for reconstruction 
  !
  DO i = 1, nSendCPU
    iCPU = SendCPU(i)    ! i is "local" number in the small communication vector, iCPU is the REAL CPU number 
    MsgLength = SendElem(i)%nElements*nVar*(2*nDim)*nDOF(2)*nDOF(3)*2   
    !MsgLength = SUM( nMPIFace(SendElem(i)%Elements(:)) ) * nVar*nDOF(2)*nDOF(3)*2  
    ALLOCATE( send_message(i)%Content(MsgLength) )
    counter = 0 
    DO j = 1, SendElem(i)%nElements 
      iElem = SendElem(i)%elements(j) 
      ! debug 
      !print *, ' send info ', myrank, iElem 
      DO iFace = 1, 2*nDim   
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nVar
          counter = counter + 1  
          send_message(i)%Content(counter) = qbnd(iVar, ii, jj, iFace, iElem )  
          counter = counter + 1  
          send_message(i)%Content(counter) = Fbnd(iVar, ii, jj, iFace, iElem )  
         ENDDO
        ENDDO
       ENDDO 
      ENDDO 
    ENDDO 
    ! Post a send request to neighbor CPU
    CALL MPI_ISEND(send_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, send_request(i), mpiErr)
  ENDDO
  !
  DO i = 1, nRecvCPU
     iCPU = RecvCPU(i)
     MsgLength = RecvElem(i)%nElements*nVar*(2*nDim)*nDOF(2)*nDOF(3)*2    
     ALLOCATE( recv_message(i)%Content(MsgLength) )
     ! Post a receive request from neighbor CPU
     CALL MPI_IRECV(recv_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, recv_request(i), mpiErr)
  ENDDO
  ! Wait until all communication has finished
  CALL MPI_WAITALL(nSendCPU,send_request,send_status_list,mpiErr)
  CALL MPI_WAITALL(nRecvCPU,recv_request,recv_status_list,mpiErr)
  !
  ! Decode messages into datastructure and deallocate messages
  !
  DO i = 1, nRecvCPU
     counter = 0 
     DO j = 1, RecvElem(i)%nElements       ! all elements received by this CPU 
      ! 
      iElem = RecvElem(i)%elements(j)  
      ! debug 
      !print *, ' recv info ', myrank, iElem 
      DO iFace = 1, 2*nDim 
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nVar
          counter = counter + 1  
          qbnd(iVar, ii, jj, iFace, iElem ) = recv_message(i)%Content(counter)   
          counter = counter + 1  
          Fbnd(iVar, ii, jj, iFace, iElem ) = recv_message(i)%Content(counter)   
         ENDDO
        ENDDO
       ENDDO 
      ENDDO 
      ! 
     ENDDO 
  ENDDO
  !
  DO i = 1, nRecvCPU 
     DEALLOCATE( recv_message(i)%Content )
  ENDDO
  DO i = 1, nSendCPU
     DEALLOCATE( send_message(i)%Content )
  ENDDO
  !
  DEALLOCATE(send_message,recv_message,send_request,recv_request,send_status_list,recv_status_list)       
  !
#endif 
  !
END SUBROUTINE MPIExchangeData  

SUBROUTINE MPIExchangeParameters     
  USE typesDef   
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
#ifdef PARALLEL 
  INCLUDE 'mpif.h'  
  INTEGER                             :: i, j, ii, jj, iFace, iCPU, msglength, iElem, iVar, iDOF, counter 
  ! 
  ALLOCATE( send_message(nSendCPU) )
  ALLOCATE( recv_message(nRecvCPU) )
  ALLOCATE( send_request(nSendCPU) )
  ALLOCATE( recv_request(nRecvCPU) )
  ALLOCATE( send_status_list(MPI_STATUS_SIZE,nSendCPU) )
  ALLOCATE( recv_status_list(MPI_STATUS_SIZE,nRecvCPU) )
  !
  ! Now exchange cell averages for reconstruction 
  !
  DO i = 1, nSendCPU
    iCPU = SendCPU(i)    ! i is "local" number in the small communication vector, iCPU is the REAL CPU number 
    MsgLength = SendElem(i)%nElements*nParam*(2*nDim)*nDOF(2)*nDOF(3)    
    !MsgLength = SUM( nMPIFace(SendElem(i)%Elements(:)) ) * nVar*nDOF(2)*nDOF(3)*2  
    ALLOCATE( send_message(i)%Content(MsgLength) )
    counter = 0 
    DO j = 1, SendElem(i)%nElements 
      iElem = SendElem(i)%elements(j) 
      ! debug 
      !print *, ' send info ', myrank, iElem 
      DO iFace = 1, 2*nDim   
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nParam
          counter = counter + 1  
          send_message(i)%Content(counter) = parbnd(iVar, ii, jj, iFace, iElem )  
         ENDDO
        ENDDO
       ENDDO 
      ENDDO 
    ENDDO 
    ! Post a send request to neighbor CPU
    CALL MPI_ISEND(send_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, send_request(i), mpiErr)
  ENDDO
  !
  DO i = 1, nRecvCPU
     iCPU = RecvCPU(i)
     MsgLength = RecvElem(i)%nElements*nParam*(2*nDim)*nDOF(2)*nDOF(3)    
     ALLOCATE( recv_message(i)%Content(MsgLength) )
     ! Post a receive request from neighbor CPU
     CALL MPI_IRECV(recv_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, recv_request(i), mpiErr)
  ENDDO
  ! Wait until all communication has finished
  CALL MPI_WAITALL(nSendCPU,send_request,send_status_list,mpiErr)
  CALL MPI_WAITALL(nRecvCPU,recv_request,recv_status_list,mpiErr)
  !
  ! Decode messages into datastructure and deallocate messages
  !
  DO i = 1, nRecvCPU
     counter = 0 
     DO j = 1, RecvElem(i)%nElements       ! all elements received by this CPU 
      ! 
      iElem = RecvElem(i)%elements(j)  
      ! debug 
      !print *, ' recv info ', myrank, iElem 
      DO iFace = 1, 2*nDim 
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nParam
          counter = counter + 1  
          parbnd(iVar, ii, jj, iFace, iElem ) = recv_message(i)%Content(counter)   
         ENDDO
        ENDDO
       ENDDO 
      ENDDO 
      ! 
     ENDDO 
  ENDDO
  !
  DO i = 1, nRecvCPU 
     DEALLOCATE( recv_message(i)%Content )
  ENDDO
  DO i = 1, nSendCPU
     DEALLOCATE( send_message(i)%Content )
  ENDDO
  !
  DEALLOCATE(send_message,recv_message,send_request,recv_request,send_status_list,recv_status_list)       
  !
#endif 
  !
END SUBROUTINE MPIExchangeParameters   


SUBROUTINE MPIExchangeDataPartI     
  USE typesDef   
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
#ifdef PARALLEL 
  INCLUDE 'mpif.h'  
  INTEGER                             :: i, j, ii, jj, iFace, iCPU, msglength, iElem, iVar, iDOF, counter 
  ! 
  ALLOCATE( send_message(nSendCPU) )
  ALLOCATE( recv_message(nRecvCPU) )
  ALLOCATE( send_request(nSendCPU) )
  ALLOCATE( recv_request(nRecvCPU) )
  ALLOCATE( send_status_list(MPI_STATUS_SIZE,nSendCPU) )
  ALLOCATE( recv_status_list(MPI_STATUS_SIZE,nRecvCPU) )
  !
  ! Now exchange cell averages for reconstruction 
  !
  DO i = 1, nSendCPU
    iCPU = SendCPU(i)    ! i is "local" number in the small communication vector, iCPU is the REAL CPU number 
    MsgLength = SendElem(i)%nElements*nVar*nDOF(2)*nDOF(3)*2   
    !MsgLength = SUM( nMPIFace(SendElem(i)%Elements(:)) ) * nVar*nDOF(2)*nDOF(3)*2  
    ALLOCATE( send_message(i)%Content(MsgLength) )
    counter = 0 
    DO j = 1, SendElem(i)%nElements 
      iElem = SendElem(i)%elements(j) 
      ! debug 
      !print *, ' send info ', myrank, iElem 
      iFace = SendElem(i)%FaceIndex(j) 
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nVar
          counter = counter + 1  
          send_message(i)%Content(counter) = qbnd(iVar, ii, jj, iFace, iElem )  
          counter = counter + 1  
          send_message(i)%Content(counter) = Fbnd(iVar, ii, jj, iFace, iElem )  
         ENDDO
        ENDDO
       ENDDO 
    ENDDO 
    ! Post a send request to neighbor CPU
    CALL MPI_ISEND(send_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, send_request(i), mpiErr)
  ENDDO
  !
  DO i = 1, nRecvCPU
     iCPU = RecvCPU(i)
     MsgLength = RecvElem(i)%nElements*nVar*nDOF(2)*nDOF(3)*2    
     ALLOCATE( recv_message(i)%Content(MsgLength) )
     ! Post a receive request from neighbor CPU
     CALL MPI_IRECV(recv_message(i)%Content, MsgLength, MPI_AUTO_REAL, iCPU, 2, MPI_COMM_WORLD, recv_request(i), mpiErr)
  ENDDO
  !
#endif 
  !
END SUBROUTINE MPIExchangeDataPartI   

SUBROUTINE MPIExchangeDataPartII    
  USE typesDef   
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
#ifdef PARALLEL 
  INCLUDE 'mpif.h'  
  INTEGER                             :: i, j, ii, jj, iFace, iCPU, msglength, iElem, iVar, iDOF, counter 
  ! 
  ! Wait until all communication has finished
  CALL MPI_WAITALL(nSendCPU,send_request,send_status_list,mpiErr)
  CALL MPI_WAITALL(nRecvCPU,recv_request,recv_status_list,mpiErr)
  !
  ! Decode messages into datastructure and deallocate messages
  !
  DO i = 1, nRecvCPU
     counter = 0 
     DO j = 1, RecvElem(i)%nElements       ! all elements received by this CPU 
      ! 
      iElem = RecvElem(i)%elements(j)  
      ! debug 
      !print *, ' recv info ', myrank, iElem 
      iFace = RecvElem(i)%FaceIndex(j) 
       DO jj = 1, nDOF(3) 
        DO ii = 1, nDOF(2) 
         DO iVar = 1, nVar
          counter = counter + 1  
          qbnd(iVar, ii, jj, iFace, iElem ) = recv_message(i)%Content(counter)   
          counter = counter + 1  
          Fbnd(iVar, ii, jj, iFace, iElem ) = recv_message(i)%Content(counter)   
         ENDDO
        ENDDO
       ENDDO        
      ! 
     ENDDO 
  ENDDO
  !
  DO i = 1, nRecvCPU 
     DEALLOCATE( recv_message(i)%Content )
  ENDDO
  DO i = 1, nSendCPU
     DEALLOCATE( send_message(i)%Content )
  ENDDO
  !
  DEALLOCATE(send_message,recv_message,send_request,recv_request,send_status_list,recv_status_list)       
  !
#endif 
  !
END SUBROUTINE MPIExchangeDataPartII  