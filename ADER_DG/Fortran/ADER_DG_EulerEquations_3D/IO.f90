SUBROUTINE WriteData 
  USE typesDef   
  USE ISO_C_BINDING
  IMPLICIT NONE 
  include 'tecio.f90' 
  CHARACTER(LEN=200) :: Filename,FilenameCSV,Title,ScratchDir, VarString   
  CHARACTER(LEN=10)  :: Name(nVar) 
  CHARACTER(LEN=10)  :: cmyrank 
  INTEGER            :: v,iElem,k,i,j,ii,jj,c,nc,iRet,iDim,iErr 
  REAL               :: QN(nVar),VN(nVar),Vav(nVar),Vmin(nVar),ldx(d),lx0(d),lxb(d),ux(nVar),uy(nVar)
  REAL               :: LocNode(nVar,(N+1)**nDim), GradNode(nVar,(N+1)**nDim,d),  xvec(d), x0(d), xGP(d)
  INTEGER            :: nSubNodes, nPlotElem, nSubPlotElem, nRealNodes, ZoneType, nVertex, nDOFs   
  REAL(8)            :: loctime 
  INTEGER*4, POINTER :: NData(:,:)  
  REAL*4, POINTER    :: DataArray(:,:),TempArray(:) 
  INTEGER*4          :: visdouble
  REAL*4             :: Test 
  POINTER   (NullPtr,Null)
  Integer*4 Null(*)
  integer, parameter :: out_unit=20
  
  !
  visdouble = 0
  nVertex = 2**nDim 
  WRITE(FileName,'(a,a1,i8.8,a)') TRIM(BaseFile),'-', timestep, '.plt'
  PRINT *, ' Writing data to file ', TRIM(FileName) 
  !
  NullPtr = 0 
  nPlotElem =  0
  nSubPlotElem = 0
  nRealNodes = 0  
  DO i = 1, nElem
     nPlotElem = nPlotElem + 1
     nSubPlotElem = nSubPlotElem + N**nDim  
     nSubNodes = (N+1)**nDim  
     nRealNodes = nRealNodes + nSubNodes 
  ENDDO
  ALLOCATE(NData(nVertex,nSubPlotElem))  
  ALLOCATE(DataArray(nRealNodes,nDim+nVar+1))  
  c  = 0 
  nc = 0 
  DO i = 1, nElem
    DO j = 1, N**nDim 
        c = c + 1  
        NData(:,c) = nc + subtri(1:nVertex,j)
    ENDDO 
    nc = nc + (N+1)**nDim   
  ENDDO
  nDOFs = PRODUCT(nDOF(1:nDim)) 
  c = 0 
  DO i = 1, nElem      
    LocNode = MATMUL( RESHAPE( uh(:,:,:,:,i), (/ nVar, nDOFs /) ), SubOutputMatrix(1:nDOFs,1:(N+1)**nDim) ) 
    lx0 = x(:,tri(1,i)) 
    DO j = 1, (N+1)**nDim  
        QN(:) = LocNode(:,j) 
        xvec = lx0 + allsubxi(:,j)*dx 
        CALL PDECons2Prim(VN,QN,iErr)
        c = c + 1 
        DataArray(c,:) = (/ xvec(1:nDim), VN, REAL(i) /)   
    ENDDO
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
  CALL PDEVarName(Name)  
  DO i = 1, nVar        
     WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(Name(i)) , ' '   
  ENDDO
  WRITE(VarString,'(a,a)') TRIM(VarString), ' iE ' 

  iret = TecIni112(TRIM(Title)//''//C_NULL_CHAR,TRIM(Varstring)//''//C_NULL_CHAR,TRIM(FileName)//''//C_NULL_CHAR,TRIM(ScratchDir)//''//C_NULL_CHAR,0,0,visdouble) 
  loctime = time 
  iRet = TecZne112('Zone'//C_NULL_CHAR, ZoneType, nRealNodes, nSubPlotElem, 0, 0, 0, 0, loctime, 0, 0, 1, 0, 0, 0, 0, 0, Null, Null, Null, 0) 
  iRet = TecDat112( nRealNodes*(nDim+nVar+1), DataArray, visdouble )
  iRet = TecNod112(NData)
  iRet = TecEnd112()

  DEALLOCATE(NData,DataArray)  
  

  if (timestep.ne.0) then
    WRITE(FilenameCSV,'(a,a1,i0.0,a)') TRIM(BaseFile),'-', timestep, '.csv'
  else
    WRITE(FilenameCSV,'(a,a)') TRIM(BaseFile),'-0.csv'
  endif
  
  open (unit=out_unit,file=FilenameCSV)  
  
  ! DOF exporter
  DO iElem = 1, nElem
    x0 = x(:,tri(1,iElem)) ! get the coordinate of the lower left node 
    DO k = 1, nDOF(3)
     DO j = 1, nDOF(2) 
      DO i = 1, nDOF(1) 
      
        xGP = x0 + (/ xiGPN(i), xiGPN(j), xiGPN(k) /)*dx(:) 
        !DO v = 1, nVAR
        !  PRINT *, uh(v,i,j,k, iElem)
        !ENDDO
        !PRINT *, xGP
        !PRINT *, time
        !PRINT *, timestep
        !PRINT *, x0
        !PRINT *, xiGPN(i)
        !PRINT *, xiGPN(j)
        !PRINT *, xiGPN(k)
        !PRINT *, FilenameCSV
        
        
        write (out_unit,'(F,A,F,A,F,A,F,A,F,A,F,A,F,A,F)') xGP(1), ',', xGP(2), ',', time, ',', uh(1,i,j,k, iElem), ',', uh(2,i,j,k, iElem), ',', uh(3,i,j,k, iElem), ',', uh(4,i,j,k, iElem), ',', uh(5,i,j,k, iElem)
        !PRINT *, 
        !PRINT *, '---------------------------------------------'
      ENDDO
     ENDDO 
    ENDDO
  ENDDO
  
  close (out_unit)
END SUBROUTINE WriteData 
    

    