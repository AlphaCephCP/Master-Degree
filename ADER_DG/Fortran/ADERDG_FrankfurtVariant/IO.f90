SUBROUTINE WriteData 
  USE typesDef   
  USE ISO_C_BINDING
  IMPLICIT NONE 
  include 'tecio.f90' 
  CHARACTER(LEN=200) :: Filename,Title,ScratchDir, VarString   
  CHARACTER(LEN=10)  :: VarName(nVar) 
  CHARACTER(LEN=10)  :: ParName(nParam) 
  CHARACTER(LEN=10)  :: cmyrank 
  INTEGER            :: i,j,ii,jj,kk,c,nc,iRet,iDim,iErr  
  REAL               :: QN(nVar),VN(nVar),Vav(nVar),Vmin(nVar),ldx(d),lx0(d),lxb(d),ux(nVar),uy(nVar)
  REAL               :: LocNode(nVar,(N+1)**nDim), GradNode(nVar,(N+1)**nDim,d),  xvec(d) 
  REAL               :: ParNode(nParam,(N+1)**nDim)
  REAL               :: outuh(nVar,1:nSubLimV(1)+1,1:nSubLimV(2)+1,1:nSubLimV(3)+1) 
  INTEGER            :: nSubNodes, nSubPlotElem, nRealNodes, ZoneType, nVertex, nDOFs   
  REAL(8)            :: loctime 
  INTEGER*4, POINTER :: NData(:,:)  
  REAL*4, POINTER    :: DataArray(:,:),TempArray(:) 
  INTEGER*4          :: visdouble
  REAL*4             :: Test 
  POINTER   (NullPtr,Null)
  Integer*4 Null(*)
  !
  visdouble = 0
  nVertex = 2**nDim 
  WRITE(FileName,'(a,a1,i8.8,a)') TRIM(BaseFile),'-', timestep, '.plt'
  PRINT *, ' Writing data to file ', TRIM(FileName) 
  !
  NullPtr = 0 
  nSubPlotElem = 0
  nRealNodes = 0  
  DO i = 1, nElem
     IF(Limiter(i)%status.EQ.0) THEN   
        nSubPlotElem = nSubPlotElem + N**nDim  
        nSubNodes = (N+1)**nDim  
     ELSE
        nSubPlotElem = nSubPlotElem + (nSubLim)**nDim  
        nSubNodes = (nSubLim+1)**nDim  
     ENDIF
     nRealNodes = nRealNodes + nSubNodes 
  ENDDO
  ALLOCATE(NData(nVertex,nSubPlotElem))  
  ALLOCATE(DataArray(nRealNodes,nDim+nVar+nParam+2))  
  c  = 0 
  nc = 0 
  DO i = 1, nElem
    IF(Limiter(i)%status.EQ.0) THEN
        DO j = 1, N**nDim
            c = c + 1  
            NData(:,c) = nc + subtri(1:nVertex,j)
        ENDDO 
        nc = nc + (N+1)**nDim  
    ELSE
        DO j = 1, (nSubLim)**nDim 
            c = c + 1 
            NData(:,c) = nc + subtrilim(1:nVertex,j)
        ENDDO 
        nc = nc + (nSubLim+1)**nDim  
    ENDIF 
    !DO j = 1, N**nDim 
    !    c = c + 1  
    !    NData(:,c) = nc + subtri(1:nVertex,j)
    !ENDDO 
    !nc = nc + (N+1)**nDim   
  ENDDO
  nDOFs = PRODUCT(nDOF(1:nDim)) 
  c = 0 
  DO i = 1, nElem 
    IF(Limiter(i)%status.EQ.0) THEN
        LocNode = MATMUL( RESHAPE( uh(:,:,:,:,i), (/ nVar, nDOFs /) ), SubOutputMatrix(1:nDOFs,1:(N+1)**nDim) ) 
#ifdef ELASTICITY
        ParNode = MATMUL( RESHAPE( parh(:,:,:,:,i), (/ nParam, nDOFs /) ), SubOutputMatrix(1:nDOFs,1:(N+1)**nDim) )  
#endif
        lx0 = x(:,tri(1,i)) 
        DO j = 1, (N+1)**nDim  
            QN(:) = LocNode(:,j) 
            xvec = lx0 + allsubxi(:,j)*dx 
            CALL PDECons2Prim(VN,QN,iErr)
            c = c + 1         
            IF(nParam.EQ.0) THEN
                DataArray(c,:) = (/ xvec(1:nDim), VN, REAL(i), REAL(Limiter(i)%status) /)   
            ELSE
                DataArray(c,:) = (/ xvec(1:nDim), VN, ParNode(:,j), REAL(i), REAL(Limiter(i)%status) /)   
            ENDIF        
        ENDDO
    ELSE
        PRINT *, "LIMITER NOT SUPPORTED IN THIS BUILD"
        EXIT 
#ifdef LIMITER
        lx0 = x(:,tri(1,i)) 
        CALL GetSubcell_uh(outuh,i) 
        nc = 0 
         DO kk = 1, nSubLimV(3)+dn(3) 
           DO jj = 1, nSubLimV(2)+dn(2) 
             DO ii = 1, nSubLimV(1)+dn(1)    
                nc = nc + 1 
                QN(:) = outuh(:,ii,jj,kk) 
                xvec = lx0 + subxilim(:,nc)*dx 
                CALL PDECons2Prim(VN,QN,iErr)
                c = c + 1 
                DataArray(c,:) = (/ xvec(1:nDim), VN, REAL(i), REAL(Limiter(i)%status) /)   
             ENDDO 
           ENDDO
        ENDDO 
#endif
    ENDIF    
  ENDDO

  WRITE(Title,'(a,f9.4,a)') 'Time t = ', time, ''//C_NULL_CHAR  
  WRITE(ScratchDir,'(a)') '.'//C_NULL_CHAR 
  SELECT CASE(nDim)
  CASE(1)
   WRITE(VarString,*) 'x ' 
   ZoneType = 1 ! FEM Line seg   
  CASE(2)
   WRITE(VarString,*) 'x y ' 
   ZoneType = 3 ! FEM Quad  
  CASE(3)
   WRITE(VarString,*) 'x y z ' 
   ZoneType = 5 ! FEM Brick 
  END SELECT 
  CALL PDEVarName(VarName)  
  DO i = 1, nVar        
     WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(VarName(i)) , ' '   
  ENDDO
  CALL PDEParName(ParName)  
  DO i = 1, nParam         
     WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(ParName(i)) , ' '   
  ENDDO
  WRITE(VarString,'(a,a)') TRIM(VarString), ' iE limiter ' 
  
  iret = TecIni112(TRIM(Title)//''//C_NULL_CHAR,TRIM(Varstring)//''//C_NULL_CHAR,TRIM(FileName)//''//C_NULL_CHAR,TRIM(ScratchDir)//''//C_NULL_CHAR,0,0,visdouble) 
  loctime = time 
  iRet = TecZne112('Zone'//C_NULL_CHAR, ZoneType, nRealNodes, nSubPlotElem, 0, 0, 0, 0, loctime, 0, 0, 1, 0, 0, 0, 0, 0, Null, Null, Null, 0) 
  iRet = TecDat112( nRealNodes*(nDim+nVar+nParam+2), DataArray, visdouble )
  iRet = TecNod112(NData)
  iRet = TecEnd112()

  DEALLOCATE(NData,DataArray)  

END SUBROUTINE WriteData 
    

    
