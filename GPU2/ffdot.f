      SUBROUTINE FFDOT(I,J,FIRR,FD)
      INCLUDE 'common6.h'
      REAL*8  FIRR(3),FD(3),A(3),DV(3)
*     DO 1 K = 1,3
*         FIRR(K) = 0.0D0
*         FD(K) = 0.0D0
*   1 CONTINUE
      
          DO 5 K = 1,3
              A(K) = X(K,J) - X(K,I)
              DV(K) = XDOT(K,J) - XDOT(K,I)
    5     CONTINUE
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2 + EPS2
          DR3I = BODY(J)/(RIJ2*SQRT(RIJ2))
          DRDV = 3.0*(A(1)*DV(1) + A(2)*DV(2) + A(3)*DV(3))/RIJ2
          DO 10 K = 1,3
              FIRR(K) = FIRR(K) + A(K)*DR3I
              FD(K) = FD(K) + (DV(K) - A(K)*DRDV)*DR3I
   10     CONTINUE
      RETURN
      END
