      SUBROUTINE FPCORR2(I,NBLOSS,NBGAIN,JJLIST,XI,XIDOT)
*
*
*       Force polynomial corrections.
*       -----------------------------
*
      INCLUDE 'common6.h'
      INTEGER JJLIST(LMAX)
      REAL*8  XI(3),XIDOT(3),SAVE1(3),SAVE2(3),SAVE3(3),
     &        A(12),F2DOT(3),F3DOT(4)
*
*
*       See whether there has been a change of neighbours.
      NBSUM = NBLOSS + NBGAIN
      IF (NBSUM.EQ.0) GO TO 100
*
*       Initialize the derivative corrections.
      DO 70 K = 1,3
          SAVE1(K) = 0.0
          SAVE2(K) = 0.0
          SAVE3(K) = 0.0
   70 CONTINUE
*
*       Form compact list of NBLOSS & NBGAIN.
      NNB0 = LIST(1,I)
      IF (NBGAIN.GT.0) THEN
          DO 78 L = 1,NBGAIN
              JJLIST(NBLOSS+L) = JJLIST(NNB0+L)
   78     CONTINUE
      END IF
*
*       Accumulate derivative corrections.
      L = 1
   80 J = JJLIST(L)
*
*       Use c.m. values of XDOT, F & FDOT for single KS components.
      IF (J.LT.IFIRST) THEN
          JCM = N + KVEC(J)
*         STEPJ = STEP(JCM)
          S = TIME - T0(JCM)
*
          DO 82 K = 1,3
              A(K) = X(K,J) - XI(K)
              A(K+3) = XDOT(K,J) - XIDOT(K)
              A(K+6) = 2.0*(F(K,JCM) - F(K,I))
              A(K+9) = 6.0*(FDOT(K,JCM) - FDOT(K,I))
   82     CONTINUE
*
          NCMDER = NCMDER + 1
      ELSE
*         STEPJ = STEP(J)
          S = TIME - T0(J)
          S3 = 3.0*S
*
*       Predict F & FDOT of body #J to order FDOT.
          DO 84 K = 1,3
              A(K) = X(K,J) - XI(K)
              A(K+3) = (FDOT(K,J)*S3 + 2.0*F(K,J))*S + X0DOT(K,J) -
     &                                                          XIDOT(K)
              A(K+6) = 2.0*(FDOT(K,J)*S3 + F(K,J) - F(K,I))
              A(K+9) = 6.0*(FDOT(K,J) - FDOT(K,I))
   84     CONTINUE
      END IF
*
      A13 = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
      A14 = BODY(J)*A13*SQRT(A13)
      A15 = (A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A13
      A16 = A15*A15
      A17 = 3.0*A15
      A18 = 6.0*A15
      A19 = 9.0*A15
      A20 = (A(4)*A(4) + A(5)*A(5) + A(6)*A(6) + A(1)*A(7) + A(2)*A(8)
     &                                            + A(3)*A(9))*A13 + A16
      A21 = 9.0*A20
      A20 = 3.0*A20
      A22 = (9.0*(A(4)*A(7) + A(5)*A(8) + A(6)*A(9)) + 3.0*(A(1)*A(10)
     &             + A(2)*A(11) + A(3)*A(12)))*A13 + A17*(A20 - 4.0*A16)
*
      DO 85 K = 1,3
          F1DOTK = A(K+3) - A17*A(K)
          F2DOT(K) = (A(K+6) - A18*F1DOTK - A20*A(K))*A14
          F3DOT(K) = (A(K+9) - A21*F1DOTK - A22*A(K))*A14 - A19*F2DOT(K)
*       Suppress F1DOT terms (already done in REGINT).
*         F1DOT(K) = F1DOTK*A14
   85 CONTINUE
*
*       Check the third derivative in case of small neighbour steps.
*     IF (STEPJ.LT.10.0*SMIN) THEN
*         F3DOT(4) = ABS(F3DOT(1)) + ABS(F3DOT(2)) + ABS(F3DOT(3))
*         A3 = ABS(D3R(1,I)) + ABS(D3R(2,I)) + ABS(D3R(3,I))
*         IF (F3DOT(4).GT.10.0*A3) THEN
*       Ignore large third derivative terms to avoid convergence problems.
*             DO 86 K = 1,3
*                 F3DOT(K) = 0.0
*  86         CONTINUE
*                 NBDER = NBDER + 1
*         END IF
*     END IF
*
*       Change the sign for NBLOSS contributions.
      IF (L.LE.NBLOSS) THEN
          DO 88 K = 1,3
*             F1DOT(K) = -F1DOT(K)
              F2DOT(K) = -F2DOT(K)
              F3DOT(K) = -F3DOT(K)
   88     CONTINUE
      END IF
*
*       Include derivative corrections from losses & gains.
      DO 90 K = 1,3
*         SAVE1(K) = SAVE1(K) + F1DOT(K)
          SAVE2(K) = SAVE2(K) + F2DOT(K)
          SAVE3(K) = SAVE3(K) + F3DOT(K)
   90 CONTINUE
*
      L = L + 1
      IF (L.LE.NBSUM) GO TO 80
*
      NBCORR = NBCORR + 1
*
*       Perform corrections to irregular and regular force derivatives.
      DO 95 K = 1,3
*       Note that corrected value of D1 & D1R already set in routine REGINT.
*         D1(K,I) = D1(K,I) + SAVE1(K)
          D2(K,I) = D2(K,I) + SAVE2(K)
          D3(K,I) = D3(K,I) + SAVE3(K)
*         D1R(K,I) = D1R(K,I) - SAVE1(K)
          D2R(K,I) = D2R(K,I) - SAVE2(K)
          D3R(K,I) = D3R(K,I) - SAVE3(K)
   95 CONTINUE
*
  100 RETURN
*
      END
