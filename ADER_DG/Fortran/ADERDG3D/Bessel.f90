MODULE Bessel_mod 

   IMPLICIT NONE

   INTERFACE JYNA
     MODULE PROCEDURE JYNA
   END INTERFACE 
   
   PUBLIC :: JYNA 
   
   CONTAINS 
   
    SUBROUTINE JYNA(N,X,NM,BJ,DJ,BY,DY)
!C
!C       ==========================================================
!C       Purpose: Compute Bessel functions Jn(x) & Yn(x) and
!C                their derivatives
!C       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ò 0 )
!C                n --- Order of Jn(x) & Yn(x)
!C       Output:  BJ(n) --- Jn(x)
!C                DJ(n) --- Jn'(x)
!C                BY(n) --- Yn(x)
!C                DY(n) --- Yn'(x)
!C                NM --- Highest order computed
!C       Routines called:
!C            (1) JY01B to calculate J0(x), J1(x), Y0(x) & Y1(x)
!C            (2) MSTA1 and MSTA2 to calculate the starting 
!C                point for backward recurrence
!C       =========================================================
!C
    ! Arguments:
    INTEGER :: N, NM
    REAL :: BJ(0:), BY(0:), DJ(0:), DY(0:)
    REAL :: X
!        DIMENSION BJ(0:N),BY(0:N),DJ(0:N),DY(0:N)
!        DIMENSION BJ(0:*),BY(0:*),DJ(0:*),DY(0:*)
    ! Local Variables:
    REAL :: BJ0, BJ1, BJK
    REAL :: BY0, BY1
    REAL :: DJ0, DJ1
    REAL :: DY0, DY1
    REAL :: CS
    REAL :: F, F0, F1, F2
 
    INTEGER :: K,M

    NM = N

    IF (X.LT.1.0D-100) THEN
       DO K=0,N
          BJ(K)=0.0D0
          DJ(K)=0.0D0
          BY(K)=-1.0+300
          DY(K)=1.0+300
       END DO
       BJ(0)=1.0D0
       DJ(1)=0.5D0
       RETURN
    ENDIF
    CALL JY01B(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    BJ(0)=BJ0
    BJ(1)=BJ1
    BY(0)=BY0
    BY(1)=BY1
    DJ(0)=DJ0
    DJ(1)=DJ1
    DY(0)=DY0
    DY(1)=DY1
    IF (N.LE.1) RETURN
    IF (N.LT.INT(0.9*X)) THEN
       DO K=2,N
          BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
          BJ(K)=BJK
          BJ0=BJ1
          BJ1=BJK
       END DO
    ELSE
       M=MSTA1(X,200)
       IF (M.LT.N) THEN
          NM=M
       ELSE
          M=MSTA2(X,N,15)
       ENDIF
       F2=0.0D0
       F1=1.0D-100
       DO K=M,0,-1
          F=2.0D0*(K+1.0D0)/X*F1-F2
          IF (K.LE.NM) BJ(K)=F
          F2=F1
          F1=F
       END DO
       IF (ABS(BJ0).GT.ABS(BJ1)) THEN
          CS=BJ0/F
       ELSE
          CS=BJ1/F2
       ENDIF
       DO K=0,NM
          BJ(K)=CS*BJ(K)
       END DO
    ENDIF
    DO K=2,NM
       DJ(K)=BJ(K-1)-K/X*BJ(K)
    END DO
    F0=BY(0)
    F1=BY(1)
    DO K=2,NM
       F=2.0D0*(K-1.0D0)/X*F1-F0
       BY(K)=F
       F0=F1
       F1=F
    END DO
    DO K=2,NM
       DY(K)=BY(K-1)-K*BY(K)/X
    END DO
    RETURN
  END SUBROUTINE JYNA


  SUBROUTINE JY01B(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
    !C
    !C       =======================================================
    !C       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
    !C                Y1(x), and their derivatives!
    !C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
    !C       Output:  BJ0 --- J0(x)
    !C                DJ0 --- J0'(x)
    !C                BJ1 --- J1(x)
    !C                DJ1 --- J1'(x)
    !C                BY0 --- Y0(x)
    !C                DY0 --- Y0'(x)
    !C                BY1 --- Y1(x)
    !C                DY1 --- Y1'(x)
    !C       =======================================================
    !C
    ! Arguments:
    REAL :: X
    REAL :: BJ0, BJ1, BY0, BY1
    REAL :: DJ0, DJ1, DY0, DY1
    ! Local Variables:
    REAL    :: A0
    REAL    :: P0, P1
    REAL    :: Q0, Q1
    REAL    :: T, T2, TA0, TA1
    
    ! Constants:
    REAL, PARAMETER :: PI = 3.141592653589793
    
    IF (X.EQ.0.0D0) THEN
       BJ0=1.0D0
       BJ1=0.0D0
       DJ0=0.0D0
       DJ1=0.5D0
       BY0=-1.0+300
       BY1=-1.0+300
       DY0=1.0+300
       DY1=1.0+300
       RETURN
    ELSE IF (X.LE.4.0D0) THEN
       T=X/4.0D0
       T2=T*T
       BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2   &
                -.0709253492D0)*T2+.4443584263D0)*T2     &
               -1.7777560599D0)*T2+3.9999973021D0)      &
             *T2-3.9999998721D0)*T2+1.0D0
       BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)   &
                   *T2-.0236616773D0)*T2+.1777582922D0)*T2  &
                  -.8888839649D0)*T2+2.6666660544D0)*T2    &
                -3.9999999710D0)*T2+1.9999999998D0)
       BY0=(((((((-.567433D-4*T2+.859977D-3)*T2     &
                 -.94855882D-2)*T2+.0772975809D0)*T2      &
                -.4261737419D0)*T2+1.4216421221D0)*T2    &
              -2.3498519931D0)*T2+1.0766115157)*T2     &
            +.3674669052D0
       BY0=2.0D0/PI*DLOG(X/2.0D0)*BJ0+BY0
       BY1=((((((((.6535773D-3*T2-.0108175626D0)*T2 &
                  +.107657606D0)*T2-.7268945577D0)*T2      &
                +3.1261399273D0)*T2-7.3980241381D0)*T2   &
             +6.8529236342D0)*T2+.3932562018D0)*T2    &
           -.6366197726D0)/X
       BY1=2.0D0/PI*DLOG(X/2.0D0)*BJ1+BY1
    ELSE
       T=4.0D0/X
       T2=T*T
       A0=SQRT(2.0D0/(PI*X))
       P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2     &
            +.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
       Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2    &
               -.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
       TA0=X-.25D0*PI
       BJ0=A0*(P0*COS(TA0)-Q0*SIN(TA0))
       BY0=A0*(P0*SIN(TA0)+Q0*COS(TA0))
       P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2     &
            -.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
       Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2   &
               +.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
       TA1=X-.75D0*PI
       BJ1=A0*(P1*COS(TA1)-Q1*SIN(TA1))
       BY1=A0*(P1*SIN(TA1)+Q1*COS(TA1))
    ENDIF
    DJ0=-BJ1
    DJ1=BJ0-BJ1/X
    DY0=-BY1
    DY1=BY0-BY1/X
    RETURN
  END SUBROUTINE JY01B


  INTEGER FUNCTION MSTA1(X,MP)
    !C
    !C       ===================================================
    !C       Purpose: Determine the starting point for backward  
    !C                recurrence such that the magnitude of    
    !C                Jn(x) at that point is about 10^(-MP)
    !C       Input :  x     --- Argument of Jn(x)
    !C                MP    --- Value of magnitude
    !C       Output:  MSTA1 --- Starting point   
    !C       ===================================================
    !C
    ! Arguments:
    REAL    :: X
    INTEGER :: MP
    ! Local Variables:
    REAL    :: A0
    REAL    :: F0, F1, F
    INTEGER :: N0, N1, NN
    INTEGER :: IT

    A0=ABS(X)
    N0=INT(1.1*A0)+1
    F0=ENVJ(N0,A0)-MP
    N1=N0+5
    F1=ENVJ(N1,A0)-MP
    DO IT=1,20             
       NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
       F=ENVJ(NN,A0)-MP
       IF(ABS(NN-N1).LT.1) exit
       N0=N1
       F0=F1
       N1=NN
       F1=F
    END DO
    MSTA1=NN
    RETURN
  END FUNCTION MSTA1


  INTEGER FUNCTION MSTA2(X,N,MP)
    !C
    !C       ===================================================
    !C       Purpose: Determine the starting point for backward
    !C                recurrence such that all Jn(x) has MP
    !C                significant digits
    !C       Input :  x  --- Argument of Jn(x)
    !C                n  --- Order of Jn(x)
    !C                MP --- Significant digit
    !C       Output:  MSTA2 --- Starting point
    !C       ===================================================
    !C
    ! Arguments:
    REAL :: X
    INTEGER :: N, MP
    ! Local Variables:
    REAL :: A0
    REAL :: EJN
    REAL :: F, F0, F1
    REAL :: HMP
    REAL :: OBJ
    INTEGER :: IT
    INTEGER :: N0, N1, NN
    
    A0=ABS(X)
    HMP=0.5D0*MP
    EJN=ENVJ(N,A0)
    IF (EJN.LE.HMP) THEN
       OBJ=MP
       N0=INT(1.1*A0)
    ELSE
       OBJ=HMP+EJN
       N0=N
    ENDIF
    F0=ENVJ(N0,A0)-OBJ
    N1=N0+5
    F1=ENVJ(N1,A0)-OBJ
    DO IT=1,20
       NN=N1-(N1-N0)/(1.0D0-F0/F1)
       F=ENVJ(NN,A0)-OBJ
       IF (ABS(NN-N1).LT.1) EXIT
       N0=N1
       F0=F1
       N1=NN
       F1=F
    END DO
    MSTA2=NN+10
    RETURN
  END FUNCTION MSTA2

  REAL FUNCTION ENVJ(N,X)
    INTEGER ::  N
    REAL :: X

    ENVJ=0.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
    RETURN
  END FUNCTION ENVJ
    
    
END MODULE Bessel_mod 
    
    