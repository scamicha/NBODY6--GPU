      SUBROUTINE CMFIRR(I,I1,FIRR,FD)
*
*
*       Irregular force corrections for active KS.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
      REAL*8  FIRR(3),FD(3),DX(3),DV(3),FP(6),FPD(6),FCM(3),FCMD(3)
*
*
*       Set perturber number and KS indices.
      NP = LIST(1,I1)
      I2 = I1 + 1
      J1 = -1
      BODYIN = 1.0/BODY(I)
*
*       Form irregular force components for perturbed KS pair.
      DO 50 LL = 2,NP+1
          J = LIST(LL,I1)
*       See whether perturber has been resolved (at irregular force time).
          IF (TPRED(J).NE.TIME) THEN
              CALL JPRED(J)
          END IF
*       Obtain individual c.m. force with single particle approximation.
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DV(1) = XDOT(1,J) - XDOT(1,I)
          DV(2) = XDOT(2,J) - XDOT(2,I)
          DV(3) = XDOT(3,J) - XDOT(3,I)
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
*       Save c.m. approximation terms for differential correction.
          FCM(1) = A1*DR3I
          FCM(2) = A2*DR3I
          FCM(3) = A3*DR3I
          FCMD(1) = (DV(1) - A1*DRDV)*DR3I
          FCMD(2) = (DV(2) - A2*DRDV)*DR3I
          FCMD(3) = (DV(3) - A3*DRDV)*DR3I
*
*       See whether to resolve binary components.
          IF (J.GT.N) THEN
              J1 = 2*(J - N) - 1
*       Adopt c.m. of unperturbed binary for consistency.
              IF (LIST(1,J1).EQ.0) THEN
                  K = J
              ELSE
                  K = J1
              END IF
          ELSE
              K = J
          END IF
*
*       Evaluate perturbation on first component due to body #K.
   20     dr2 = 0.0
          drdv = 0.0
          DO 25 L = 1,3
              dx(L) = X(L,K) - X(L,I1)
              dv(L) = XDOT(L,K) - XDOT(L,I1)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   25     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 30 L = 1,3
              FP(L) = dx(L)*dr3i
              FPD(L) = (dv(L) - dx(L)*drdv)*dr3i
   30     CONTINUE
*
*       Evaluate perturbation on second component due to body #K.
          dr2 = 0.0
          drdv = 0.0
          DO 32 L = 1,3
              dx(L) = X(L,K) - X(L,I2)
              dv(L) = XDOT(L,K) - XDOT(L,I2)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   32     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(K)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 35 L = 1,3
              FP(L+3) = dx(L)*dr3i
              FPD(L+3) = (dv(L) - dx(L)*drdv)*dr3i
   35     CONTINUE
*
*       Accumulate individual correction terms.
          DO 40 L = 1,3
              FIRR(L) = FIRR(L) + ((BODY(I1)*FP(L) +
     &                              BODY(I2)*FP(L+3))*BODYIN - FCM(L))
              FD(L) = FD(L) + ((BODY(I1)*FPD(L) +
     &                          BODY(I2)*FPD(L+3))*BODYIN - FCMD(L))
   40     CONTINUE
*
*       Include possible second KS component without c.m. contribution.
          IF (K.EQ.J1) THEN
              DO 45 L = 1,3
                  FCM(L) = 0.0
                  FCMD(L) = 0.0
   45         CONTINUE
              K = K + 1
              GO TO 20
          END IF
   50 CONTINUE
*
*       Check force correction due to regularized chain.
      IF (NCH.GT.0) THEN
          CALL KCPERT(I,I1,FIRR,FD)
      END IF 
*
      RETURN
*
      END
