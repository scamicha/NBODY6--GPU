      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      PARAMETER (NIMAX=1024,NPMAX=16)
      REAL*8   H2I(NIMAX),XI(3,NIMAX),VI(3,NIMAX),GPUACC(3,NIMAX),
     &         GPUJRK(3,NIMAX),GPUPHI(NIMAX)
      INTEGER  NXTLST(NMAX),IBL(LMAX),NBLIST(NMAX),LISTQ(NMAX),NL(20)
      INTEGER  IRR(NMAX),IREG(NMAX),LISTGP(LMAX,NIMAX)
      LOGICAL LOOP,LSTEPM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE,ISTART
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
      DATA ISAVE,JSAVE,ISTART /0,0,0/
      SAVE TLISTQ,TMIN,NPACT
*     SAVE CPRED
*     DATA CPRED /0.0D0/
*
*
*       Enforce level search on return, except new and terminated KS.
      IF (IPHASE.NE.1.AND.IPHASE.NE.2) LOOP = .TRUE.
*
*       Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = INT((FLOAT(NZERO)/1024.0)**0.333333)
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
      END IF
*
*       Open the GPU library on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = N
          CALL GPUNB_OPEN(NN)
          ISTART = 1
*       Define parameter for predicting active particles in library.
          IF (N.LE.5000) THEN
              NPACT = 20
          ELSE IF (N.LE.50000) THEN
              NPACT = 25
          ELSE IF (N.LE.100000) THEN
              NPACT = 40
          ELSE
              NPACT = 50
          END IF
      END IF
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
      DTM = 1.0
      TPREV = TIME
*       Initialize end-point of integration times and set DTM.
!$omp parallel do private(I)
      DO 1000 I = IFIRST,NTOT
          TNEW(I) = T0(I) + STEP(I)
*         DTM = MIN(DTM,STEP(I))
1000 CONTINUE
!$omp end parallel do
*
*       Take typical small quantized step for lower search index.
      CALL STEPK(DTMIN,DTM)
*       Determine level for the smallest step (ignore extreme values).
      LQS = 20
      DO 1001 L = 6,20
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
1001 CONTINUE
*
*       Specify upper level for optimized membership.
      LQB = LQS - 4
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*       Enforce new block step search initially and on significant change.
      TLISTQ = TIME
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    1 ICALL = ICALL + 1
*       Reset TMIN second & third time after change to catch new chain step.
      IF (TIME.GE.TLISTQ.OR.ICALL.LE.3) THEN
*       Update interval by optimization at major times (sqrt of N-NPAIRS).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,20
                  NL(L) = 0
   10         CONTINUE
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                  DO 12 L = LQB,LQS
                      IF (STEP(I).LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N-N_b).
              NSQ = SQRT(FLOAT(N - NPAIRS))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEW(I).LE.TLISTQ) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEW(I),TMIN)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = LQ - 1
      END IF
*
*       Find all particles in next block (TNEW = TMIN).
      CALL INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = TMIN
      TBLOCK = TIME
*
*     WRITE (6,22)  I, NXTLEN, NSTEPU, NSTEPI, TIME, STEP(I), STEPR(I)
*  22 FORMAT (' INTGRT   I LEN #U #I T S SR  ',2I6,2I11,F9.4,1P,2E10.2)
*     CALL FLUSH(6)
      IF (STEP(I).LT.1.0D-10) THEN
          WRITE (6,24)  I, NAME(I), NXTLEN, NSTEPR, STEP(I), STEPR(I)
   24     FORMAT (' SMALL STEP!!    I NAME LEN #R SI SR ',
     &                              3I6,I11,1P,2E10.2)
      END IF
*
*       Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 1
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
          END IF
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,25)  I, NAME(I), NSTEPI, TIME, STEP(I), TIME/STEP(I)
*  25     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I6,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check for new regularization at end of previous block.
      IF (IKS.GT.0) THEN
          IPHASE = 1
*       Copy the saved KS component indices and time.
          ICOMP = ISAVE
          JCOMP = JSAVE
          TIME = TSAVE
          GO TO 100
      END IF
*
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance ARchain or KS at first new time.
      IF (TIME.GT.TPREV) THEN
          CALL SUBINT(IQ,I10)
          IF (IQ.LT.0) GO TO 999
      END IF
*
*       Define array for regular force condition (STEPR is commensurate).
      NFR = 0
      DO 28 L = 1,NXTLEN
          J = NXTLST(L)
          IF (TNEW(J).GE.T0R(J) + STEPR(J)) THEN
              NFR = NFR + 1
              IREG(NFR) = J
              IRR(L) = J
          ELSE
              IRR(L) = 0
          END IF
   28 CONTINUE
*
*       Decide between merging of neighbour lists or full N prediction.
      IF (NXTLEN.LT.NPACT.AND.NFR.EQ.0) THEN
*
*       Initialize pointers for neighbour lists.
!$omp parallel do private(L)
          DO 30 L = 1,NXTLEN
              IBL(L) = NXTLST(L)
   30     CONTINUE
!$omp end parallel do
*
*       Merge all neighbour lists (with absent members of IBL added).
          CALL NBSORT(NXTLEN,IBL,NNB,NBLIST)
*
*       Predict coordinates & velocities of neighbours and #I to order FDOT.
          NBPRED = NBPRED + NNB
          DO 35 L = 1,NNB
              J = NBLIST(L)
              S = TIME - T0(J)
              S1 = 1.5*S
              S2 = 2.0*S
              X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
              X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
              X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
              XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
              XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
              XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*       Check prediction of perturbed binary.
              IF (J.GT.N) THEN
                  JPAIR = J - N
                  J1 = 2*JPAIR - 1
                  IF (LIST(1,J1).GT.0) THEN
                      ZZ = 1.0
*       Distinguish between low and high-order prediction of U & UDOT.
                      IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
                      CALL KSRES2(JPAIR,J1,J2,ZZ)
                  END IF
              END IF
   35     CONTINUE
      ELSE
          NNPRED = NNPRED + 1
*         TT3 = DBLTIM()
*         DO 40 J = IFIRST,NTOT
*             IF (BODY(J).EQ.0.0D0) GO TO 40
*             S = TIME - T0(J)
*             S1 = 1.5*S
*             S2 = 2.0*S
*             X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
*             X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
*             X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
*             XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
*             XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
*             XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*  40     CONTINUE
*       Perform full-N prediction using Keigo's C++ procedure.
      CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT)
*
*     TT4 = DBLTIM()
*     CPRED = CPRED + (TT4 - TT3)
*     IF (DMOD(TIME,2.0D0).EQ.0.0D0) WRITE (6,540)  CPRED
* 540 FORMAT (' TIMING CXVPRED  ',1P,E10.2)
*       Resolve all perturbed KS pairs at frequent intervals.
          JJ = -1
          DO 45 JPAIR = 1,NPAIRS
              JJ = JJ + 2
              IF (LIST(1,JJ).GT.0) THEN
                  ZZ = 1.0
*       Distinguish between low and high-order prediction of U & UDOT.
                  IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
                  CALL KSRES2(JPAIR,J1,J2,ZZ)
              END IF
   45     CONTINUE
      END IF
*
*       Save new time (output time at TIME > TADJ) and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
      IKS0 = IKS
*
*       Predict chain variables and perturber cordinates at new block-time.
      IF (NCH.GT.0) THEN
          IF (TNEW(ICH).GT.TBLOCK) THEN
              CALL XVPRED(ICH,0)
          END IF
          CALL XCPRED(2)
      END IF
*
*       Choose between standard and parallel irregular integration.
      IF (NXTLEN.LE.NPMAX) THEN
*       Obtain all irregular forces in one loop.
      DO 48 II = 1,NXTLEN
*       Define pointer array and scalar for regular force condition.
          I = NXTLST(II)
*
*       Advance the irregular steps.
          CALL NBINT(I,IKS,IRR(II))
*
*       Save indices and TIME of first KS candidates in the block.
          IF (IKS0.EQ.0.AND.IKS.GT.0) THEN
              ISAVE = ICOMP
              JSAVE = JCOMP
              TSAVE = TIME
              IKS0 = IKS
          END IF
   48 CONTINUE
*
      ELSE
*       Obtain all irregular forces in parallel.
!$omp parallel do private(II, I)
      DO 50 II = 1,NXTLEN
*       Define pointer array and scalar for regular force condition.
          I = NXTLST(II)
*
*       Advance the irregular steps (parallel version).
          CALL NBINTP(I,IRR(II))
   50 CONTINUE
!$omp end parallel do
      END IF
      NSTEPI = NSTEPI + NXTLEN
*
*       Begin regular force loop (general case NFR > 0, counter NBLCKR).
      IF (NFR.GT.0) THEN
      NBLCKR = NBLCKR + 1
*       Include c.m. particles in the neighbour list (new 08/10).
      NN = NTOT - IFIRST + 1
      CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*     CALL GPUNB_REGF(NI,H2I,XI,VI,GPUACC,GPUJRK,LMAX,NNBMAX,LISTGP)
*
*       Perform regular force loop.
      NOFL(1) = 0
      NOFL2 = 0
      JNEXT = 0
      DO 55 II = 1,NFR,NIMAX
  550     NI = MIN(NFR-JNEXT,NIMAX)
*       Copy neighbour radius and state vector for each block.
!$omp parallel do private(LL, I, K)
          DO 52 LL = 1,NI
              I = IREG(JNEXT+LL)
              H2I(LL) = RS(I)**2
              DO 51 K = 1,3
                  XI(K,LL) = X(K,I)
                  VI(K,LL) = XDOT(K,I)
   51         CONTINUE
   52     CONTINUE
!$omp end parallel do
*
          CALL GPUNB_REGF(NI,H2I,XI,VI,GPUACC,GPUJRK,GPUPHI,LMAX,
     &                                               NBMAX,LISTGP)
          DO 54 LL = 1,NI
              I = IREG(JNEXT+LL)
              NNB = LISTGP(1,LL)
*       Repeat last block with reduced RS(I) on NNB < 0 (at end of loop).
              IF (NNB.LT.0) THEN
                  WRITE (41,556)  NSTEPR, NAME(I), LIST(1,I), RS(I)
  556             FORMAT (' OVERFLOW!   #R NAME NB0 RS ',I11,I6,I4,F8.2)
                  CALL FLUSH(41)
                  NOFL(1) = NOFL(1) + 1
                  RS(I) = 0.9*RS(I)
                  H2I(LL) = RS(I)**2
                  NOFL2 = NOFL2 + 1
              END IF
   54     CONTINUE
          IF (NOFL2.GT.0) THEN
              NOFL2 = 0
              GO TO 550
           END IF
!$omp parallel do private(LL, I, ITEMP, NNB, L1, L)
          DO 56 LL = 1, NI
              I = IREG(JNEXT+LL)
              NNB = LISTGP(1,LL)
*       Copy neighbour list but skip self-interaction.
              L1 = 1
              DO 53 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   53         CONTINUE
              LISTGP(1,LL) = L1 - 1
   56     CONTINUE
!$omp end parallel do
*       Obtain irregular force and perform regular force correction.
!$omp parallel do private(LL, I)
          DO 57 LL = 1, NI
              I = IREG(JNEXT+LL)
              CALL GPUCOR(I,XI(1,LL),VI(1,LL),GPUACC(1,LL),
     &                                        GPUJRK(1,LL),LISTGP(1,LL))
*       Note errors in the potential when using predicted coordinates.
*             POT = POT + BODY(I)*GPUPHI(LL)
   57     CONTINUE
!$omp end parallel do
          JNEXT = JNEXT + NI
   55 CONTINUE
*
*       Accumulate the sum of overflows (NOFL(1) holds current number).
      NOFL(2) = NOFL(2) + NOFL(1)
      NSTEPR = NSTEPR + NFR
*
      END IF
*
*       Determine next block time (note STEP may shrink in GPUCOR).
      DO 350 L = 1,NXTLEN
*       Determine next block time (note STEP may shrink in GPUCOR).
          I = NXTLST(L)
          IF (TNEW(I).LT.TMIN) THEN
              TMIN = TNEW(I)
          END IF
  350 CONTINUE
*
*       Copy current coordinates & velocities from corrected values, IQ.NE.0.
      IF (IKS.GT.0.OR.IQ.NE.0.OR.TIME.GE.TADJ.OR.TIME.GE.TNEXT.OR.
     &   (KZ(19).GE.3.AND.DMOD(TIME,STEPX).EQ.0.0D0)) THEN
*       Note: need copy X & XDOT at output (skip in XVPRED); also cf. #31 >0).
!$omp parallel do private(I, L, K)
      DO 60 L = 1,NXTLEN
*       Determine next block time (note STEP may shrink in GPUCOR).
          I = NXTLST(L)
*         TMIN = MIN(TNEW(I),TMIN)
          DO 58 K = 1,3
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   58     CONTINUE
   60 CONTINUE
!$omp end parallel do
      END IF
*
*       Check integration of tidal tail members.
      IF (NTAIL.GT.0) THEN
*       Allow large quantized interval with internal iteration.
          IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
              DO 65 J = ITAIL0,NTTOT
                  IF (TNEW(J).LE.TIME) THEN
                      CALL NTINT(J)
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.NE.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (KZ(37).GT.1.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time at end of block-step.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
              ELSE
                  CALL MLOSS
              END IF
              IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.NMAX) GO TO 1
      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.100*NMAX.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).EQ.3) CALL MYDUMP(1,1)
      END IF
*
*       Check option for general binary search.
      IF (KZ(4).GT.0.AND.TIME - TLASTS.GT.DELTAS) THEN  
          CALL EVOLVE(0,0)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUB.EQ.0)  WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 75 L = 1,NSUB
              STEPS(L) = 0.0D0
   75     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Obtain elapsed wall-clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      SECS = 3600.0*IHOUR + 60.0*IMIN + ISEC
      WTOT = WTOT + SECS - WTOT0
      WTOT0 = SECS
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          WT = WTOT/3600.0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME+TOFF, TCOMP, CPUTOT/60.0, ERRTOT, DETOT, WT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6,'  WTOT =',F7.1)
      END IF
*
*       Close the library and stop.
      CALL GPUNB_CLOSE
      STOP
*
*       Set current global time.
  100 TTOT = TIME + TOFF
*
      RETURN
*
      END
