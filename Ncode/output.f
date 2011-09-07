      SUBROUTINE OUTPUT
*
*
*       Output and data save.
*       ---------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/GALAXY/  GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &                OMEGA,DISK,A,B,V02,RL2
      COMMON/CLOUDS/  XCL(3,MCL),XDOTCL(3,MCL),BODYCL(MCL),RCL2(MCL),
     &                CLM(MCL),CLMDOT(MCL),CLDOT,VCL,SIGMA,RB2,PCL2,
     &                TCL,STEPCL,NCL,NEWCL
      COMMON/ECHAIN/  ECH
      REAL*8  X1(3,4),V1(3,4),UI(4),VI(4),XREL2(3),VREL2(3)
      REAL*4  XS(3,NMAX),VS(3,NMAX),BODYS(NMAX),AS(20)
      REAL*4  XJ(3,6),VJ(3,6),BODYJ(6)
      LOGICAL  FIRST,SECOND,THIRD
      SAVE  FIRST,SECOND,THIRD
      DATA  FIRST,SECOND ,THIRD/.TRUE.,.TRUE.,.TRUE./
*
*
*       Obtain energy error in case routine ADJUST not called recently.
      IF (TIME.GE.TADJ.OR.TIME.LE.0.0D0) GO TO 10
*
*       Predict X & XDOT for all particles (except unperturbed pairs).
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time (resolve all KS pairs).
      CALL ENERGY
*
*       Include KS pairs, triple, quad, mergers, collisions & chain.
      ETOT = ZKIN - POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL + EMDOT
     &                                                         + ECDOT
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      BE(2) = BE(3)
      BE(3) = ETOT
      DE = BE(3) - BE(2)
      DETOT = DETOT + DE
      DE = DE/MAX(ZKIN,ABS(ETOT))
*       Save sum of relative energy error for main output and accumulate DE.
      ERROR = ERROR + DE
      VIR = POT
*
*       Find density centre & core radius (Casertano & Hut, Ap.J. 298, 80).
      IF (N.GE.20.AND.KZ(29).EQ.0) THEN
          CALL CORE
      END IF
*
*       Check optional sorting of Lagrangian radii & half-mass radius.
      IF (KZ(7).GT.0) THEN
          CALL LAGR(RDENS)
      END IF
*
*       Initialize diagnostic variables.
   10 NP = 0
      IUNP = 0
      AMIN = 100.0
      MULT = 0
*
*       Find smallest semi-major axis and count unperturbed KS pairs.
      DO 20 IPAIR = 1,NPAIRS
          NP = NP + LIST(1,2*IPAIR-1)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          IF (SEMI.GT.0.0) AMIN = MIN(AMIN,SEMI)
          IF (LIST(1,2*IPAIR-1).EQ.0) IUNP = IUNP + 1
          IF (NAME(N+IPAIR).LT.-2*NZERO) MULT = MULT + 1
   20 CONTINUE
*
*       Include search of any hierarchical triples.
      DO 25 IM = 1,NMERGE
          ZMB = CM(1,IM) + CM(2,IM)
          SEMI = -0.5*ZMB/HM(IM)
          AMIN = MIN(AMIN,SEMI)
   25 CONTINUE
*
*       Perform time-step & neighbour statistics (NS is # single stars).
      DTI = 0.0
      DTRI = 0.0
      CNNB = 0.0
      CMAX = 0.0
      NNB = 0
      NS = 0
      SUM = 0.0
      DO 30 I = IFIRST,NTOT
          DTI = DTI + 1.0/STEP(I)
          DTRI = DTRI + 1.0/STEPR(I)
          CNNB = CNNB + LIST(1,I)/STEP(I)
          RHO = LIST(1,I)/RS(I)**3
          CMAX = MAX(CMAX,RHO)
          NNB = NNB + LIST(1,I)
          IF (I.LE.N.AND.BODY(I).GT.0.0D0) NS = NS + 1
          SUM = SUM + BODY(I)**2
   30 CONTINUE
      NS = NS - NSUB
*
*       Estimate relative cost & effective neighbour number of AC scheme.
      COST = CNNB/(FLOAT(N - NPAIRS)*DTRI)
      CNNB = CNNB/DTI
*       Scale maximum particle density contrast by the mean value.
      CMAX = 2.0*CMAX*RSCALE**3/FLOAT(N)
*
*       Set average neighbour number & density centre displacement.
      NNB = FLOAT(NNB)/FLOAT(N - NPAIRS)
      RD = SQRT(RDENS(1)**2 + RDENS(2)**2 + RDENS(3)**2)
*
*       Check print frequency indicator & optional model counter.
      NPRINT = NPRINT + 1
      IF (NPRINT.GT.NFIX.OR.TIME.LE.0.0D0) THEN
          NPRINT = 1
          IF (KZ(3).GT.0) MODEL = MODEL + 1
      END IF
*
*       Form binary & merger energy ratios and count escapers (or suppress).
      EB = EBIN/(ZKIN - POT)
      EM = EMERGE/(ZKIN - POT)
      IF (KZ(21).GT.1) THEN
          CALL JACOBI(NESC)
      ELSE
          NESC = 0
      END IF
*
*       Print main output diagnostics.
      I6 = TSCALE*TTOT
*
      WRITE (6,40)  TTOT, N, NNB, NPAIRS, NMERGE, MULT, NS, NSTEPI,
     &              NSTEPB, NSTEPR, NSTEPU, ERROR, BE(3)
   40 FORMAT (//,' T =',F7.0,'  N =',I6,'  <NB> =',I3,'  KS =',I5,
     &           '  NM =',I3,'  MM =',I2,'  NS =',I6,
     &           '  NSTEPS =',I11,I10,2I11,'  DE =',1P,E9.1,
     &           '  E =',0P,F10.6)
*
      IF (KZ(21).GT.0) THEN
          CALL CPUTIM(TCOMP)
          IF (VC.EQ.0.0D0) VC = RSCALE/TCR
          TRC = 1.02*FLOAT(NC)**2*BODYM/(VC**3*LOG(FLOAT(NC)))
          DMIN1 = MIN(DMIN1, DMIN2, DMIN3, DMIN4, DMINC)
          NEFF = ZMASS**2/SUM
*
          WRITE (6,45)  NRUN, MODEL, TCOMP, TRC, DMIN1, DMIN2, DMIN3,
     &                  DMIN4, AMIN, RMAX, RSMIN, NEFF
   45     FORMAT (/,' NRUN =',I3,'  M# =',I3,'  CPU =',F8.1,'  TRC =',
     &                        F5.1, '  DMIN =',1P,4E8.1,'  AMIN =',E8.1,
     &                '  RMAX =',E8.1,'  RSMIN =',0P,F5.2,'  NEFF =',I6)
      END IF
      VRMS = SQRT(0.5*ZMASS/RSCALE)*VSTAR
*
      WRITE (6,50)
   50 FORMAT (/,'    <R>  RTIDE  RDENS   RC    NC   MC   RHOD   RHOM',
     &                 '  CMAX   <Cn>   Ir/R   UN    NP   RCM    VCM',
     &                 '      AZ      EB/E   EM/E   TCR     T6  NESC',
     &                 '  VRMS')
*
      WRITE (6,55)  RSCALE, RTIDE, RD, RC, NC, ZMC, RHOD, RHOM, CMAX,
     &              CNNB, COST, IUNP, NP, CMR(4), CMRDOT(4), AZ, EB, EM,
     &              TCR, I6, NESC, VRMS
   55 FORMAT (' #1',F5.2,F6.1,F6.2,F7.3,I5,F7.3,F6.0,F7.0,F6.0,F6.1,
     &                   F7.3,I5,I6,F8.3,F8.4,F9.4,2F7.3,F6.2,2I6,F6.1)
*
      WRITE (6,60)
   60 FORMAT (/,7X,'NNPRED    NBCORR  NBFULL  NBVOID    NICONV',
     &           '  NBSMIN  NBDIS  NBDIS2  NCMDER  NBDER  NFAST',
     &           '  NBFAST    NBLOCK    NBLCKR     NBPRED     NBFLUX',
     &           '  NIRECT  NRRECT  NURECT  NBRECT')
      WRITE (6,65)  NNPRED, NBCORR, NBFULL, NBVOID, NICONV,
     &              NBSMIN, NBDIS, NBDIS2, NCMDER, NBDER, NFAST,
     &              NBFAST, NBLOCK, NBLCKR, NBPRED, NBFLUX, NIRECT,
     &              NRRECT, NURECT, NBRECT
   65 FORMAT (' #2',I10,I10,I8,I9,I9,I8,I7,2I8,2I7,I8,2I10,2I11,4I8)
*
      WRITE (6,70)
   70 FORMAT (/,5X,'NKSTRY  NKSREG  NKSHYP     NKSPER  NPRECT  NEWKS ',
     &           '  NKSMOD    NTTRY  NTRIP  NQUAD  NCHAIN  NMERG',
     &           '  NEWHI  NSTEPT  NSTEPQ  NSTEPC')
      WRITE (6,75)  NKSTRY, NKSREG,  NKSHYP, NKSPER, NPRECT, NEWKS,
     &              NKSMOD, NTTRY, NTRIP, NQUAD, NCHAIN, NMERG, NEWHI,
     &              NSTEPT, NSTEPQ, NSTEPC
   75 FORMAT (' #3',I9,I7,I8,I11,I8,I7,2I9,2I7,I8,2I7,3I8)
*
*       Check output for mass loss or tidal capture.
      IF (KZ(19).GT.0.OR.KZ(27).GT.0) THEN
          CALL EVENTS
      END IF
*
*       Obtain half-mass radii for two groups (NAME <= NZERO/5 & > NZERO/5).
      IF (KZ(7).GE.2) THEN
          CALL LAGR2(RDENS)
      END IF
*
*       Include diagnostics about cluster orbit in general external field.
      IF (KZ(14).EQ.3) THEN
          GZ = RG(1)*VG(2) - RG(2)*VG(1)
          SX = RBAR/1000.0
          WRITE (6,78)  NTAIL, (RG(K)*SX,K=1,3), (VG(K)*VSTAR,K=1,3),
     &                  GZ, ETIDE
   78     FORMAT (/,5X,'CLUSTER ORBIT    NT RG VG JZ ET ',
     &                                 I5,3F7.2,2X,3F7.1,1P,E16.8,E10.2)
      END IF
      IF (KZ(14).EQ.4) THEN
          WRITE (6,80)  TTOT, N, RSCALE, ZMASS, MP, DETOT
   80     FORMAT (/,5X,'GAS EXPULSION    T N <R> M MP DETOT ',
     &                                   F7.1,I7,3F7.3,1P,E10.2)
      END IF
*
*       Reset minimum encounter distances & maximum apocentre separation.
      DMIN2 = 100.0
      DMIN3 = 100.0
      DMIN4 = 100.0
      DMINC = 100.0
      RSMIN = 100.0
      RMAX = 0.0
*
*       Check integer overflows (2^{32} or 2.1 billion).
      IF (NSTEPI.GT.2000000000.OR.NSTEPI.LT.0) THEN
          NSTEPI = 0
          NIRECT = NIRECT + 1
      END IF
      IF (NSTEPR.GT.2000000000.OR.NSTEPR.LT.0) THEN
          NSTEPR = 0
          NRRECT = NRRECT + 1
      END IF
      IF (NSTEPU.GT.2000000000.OR.NSTEPU.LT.0) THEN
          NSTEPU = 0
          NURECT = NURECT + 1
      END IF
      IF (NBPRED.GT.2000000000.OR.NBPRED.LT.0) THEN
          NBPRED = 0
      END IF
      IF (NBFLUX.GT.2000000000.OR.NBFLUX.LT.0) THEN
          NBFLUX = 0
          NBRECT = NBRECT + 1
      END IF
      IF (NBCORR.GT.2000000000.OR.NBCORR.LT.0) THEN
          NBCORR = 0
      END IF
      IF (NBLOCK.GT.2000000000.OR.NBLOCK.LT.0) THEN
          NBLOCK = 0
      END IF
*
*       Exit if error exceeds restart tolerance (TIME < TADJ means no CHECK).
      IF (ABS(ERROR).GT.5.0*QE.AND.TIME.LT.TADJ) GO TO 100
*
*       Check optional analysis & output of KS binaries.
      IF (KZ(8).GT.0.AND.NPAIRS.GT.0) THEN
          CALL BINOUT
      END IF
*
*       Include optional diagnostics of block-steps.
      IF (KZ(33).GT.0) THEN
          CALL LEVELS
      END IF
*
*       Check option for neighbour list histogram (same as for STEPR).
      IF (KZ(33).GT.1) THEN
          CALL NBHIST
      END IF
*
*       Check optional output of single bodies & binaries.
      IF (KZ(9).GT.0.OR.KZ(6).GT.0) THEN
          CALL BODIES
      END IF
*
*       Include optional diagnostics for interstellar clouds on unit 47.
      IF (KZ(13).GT.0) THEN
          WRITE (47,82)  TSCALE*TTOT, NS, NEWCL, DETOT, E(3),
     &                   RBAR*SQRT(PCL2), RBAR*RSCALE, RBAR*RC
   82     FORMAT (' ',F8.1,I7,I6,2F10.6,F6.2,2F7.3)
          CALL FLUSH(47)
          PCL2 = RB2
      END IF
*
*       See whether to write data bank of binary diagnostics on unit 9.
      IF (KZ(8).GE.2.AND.NPAIRS.GT.0) THEN
          CALL BINDAT
          IF (KZ(8).GT.3) THEN
              CALL HIDAT
          END IF
      END IF
*
*       Check optional diagnostics of evolving stars.
      IF (KZ(12).GT.0.AND.TIME.GE.TPLOT) THEN
          CALL HRPLOT
      END IF
*
*       Check optional writing of data on unit 3 (frequency NFIX). 
      IF (KZ(3).EQ.0.OR.NPRINT.NE.1) GO TO 100
      IF (KZ(3).GT.2.AND.KZ(3).NE.5) GO TO 99
*
      AS(1) = TTOT
      AS(2) = FLOAT(NPAIRS)
      AS(3) = RBAR
      AS(4) = ZMBAR
      AS(5) = RTIDE
      AS(6) = TIDAL(4)
      AS(7) = RDENS(1)
      AS(8) = RDENS(2)
      AS(9) = RDENS(3)
      AS(10) = TTOT/TCR
      AS(11) = TSCALE
      AS(12) = VSTAR
      AS(13) = RC
      AS(14) = NC
      AS(15) = VC
      AS(16) = RHOM
      AS(17) = CMAX
      AS(18) = RSCALE
      AS(19) = RSMIN
      AS(20) = DMIN1
*
*       Convert masses, coordinates & velocities to single precision.
      DO 90 I = 1,NTOT
          BODYS(I) = BODY(I)
          DO 85 K = 1,3
              XS(K,I) = X(K,I)
              VS(K,I) = XDOT(K,I)
   85     CONTINUE
   90 CONTINUE
*
*       Replace any ghosts by actual M, R & V (including 2 binaries).
      DO 95 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
*       Determine merger & ghost index for negative c.m. name.
          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
*       Include possible quartet [[B,S],S] and quintet [[B,S],B] first.
              IF (NAME(ICM).LT.-2*NZERO) THEN
*       Find ghost and merger index at previous hierarchical level.
                  IM = 1
                  DO 92 K = 1,NMERGE
                      IF (NAMEM(K).EQ.NAME(ICM) + 2*NZERO) IM = K
   92             CONTINUE
                  JG = N
                  DO 94 K = IFIRST,NTOT
                      IF (NAMEG(IM).EQ.NAME(K)) JG = K
   94             CONTINUE
*       Determine the current ghost and merger index.
                  CALL FINDJ(J1,JG2,IM2)
*       Distinguish netween quartet and quintet.
                  IF (NAME(J2).LE.NZERO) THEN
                      IF (JG.LE.N) THEN
                          BODYS(JG) = CM(2,IM)
                          BODYS(J1) = CM(1,IM2)
                          BODYS(JG2) = CM(2,IM2)
                      ELSE
                          JP = JG - N
                          BODYS(2*JP-1) = CM(3,IM)
                          BODYS(2*JP) = CM(4,IM)
                          BODYS(JG2) = CM(2,IM2)
                      END IF
                  ELSE
                      IF (JG.LE.N) THEN
                          BODYS(JG) = CM(2,IM)
                          BODYS(JG2) = CM(2,IM2)
                      ELSE
                          JP = JG - N
                          BODYS(2*JP-1) = CM(3,IM2)
                          BODYS(2*JP) = CM(4,IM2)
                      END IF
                  END IF
                  GO TO 95
              END IF
*
              CALL FINDJ(J1,J,IM)
*       Note: J is ghost index and IM is merger index.
              IF (J.LE.0) GO TO 95
              BODYS(J1) = CM(1,IM)
              BODYS(J) = CM(2,IM)
              ZMB = CM(1,IM) + CM(2,IM)
*       Form global coordinates and velocities from c.m. with XREL & VREL.
              DO K = 1,3
                  X1(K,1) = X(K,J1) + CM(2,IM)*XREL(K,IM)/ZMB
                  X1(K,2) = X(K,J1) - CM(1,IM)*XREL(K,IM)/ZMB
                  V1(K,1) = XDOT(K,J1) + CM(2,IM)*VREL(K,IM)/ZMB
                  V1(K,2) = XDOT(K,J1) - CM(1,IM)*VREL(K,IM)/ZMB
*
                  XS(K,J1) = X1(K,1)
                  XS(K,J)  = X1(K,2)
                  VS(K,J1) = V1(K,1)
                  VS(K,J)  = V1(K,2)
              END DO
*       See whether ghost is second merged binary (QUAD) instead of triple.
              IF (J.GT.N) THEN
                  ICM2 = J
*       Save outer binary masses as well as second KS component & ghost c.m.
                  IPAIR = ICM2 - N
                  I1 = 2*IPAIR - 1
                  I2 = I1 + 1
                  BODYS(I1) = CM(3,IM)
                  BODYS(I2) = CM(4,IM)
                  BODYS(J2) = CM(2,IM)
                  BODYS(J) = CM(3,IM) + CM(4,IM)
*       Copy KS variables to local scalars.
                  DO K = 1,4
                      UI(K) = U(K,IPAIR)
                      VI(K) = UDOT(K,IPAIR)
                  END DO
*       Transform to physical variables and multiply by 4 (momentum formula).
                  CALL KSPHYS(UI,VI,XREL2,VREL2)
                  ZM = CM(3,IM) + CM(4,IM)
                  DO K = 1,3
                      VREL2(K) = 4.0*VREL2(K)
                      X1(K,3) = X(K,J2) + CM(4,IM)*XREL2(K)/ZM
                      X1(K,4) = X(K,J2) - CM(3,IM)*XREL2(K)/ZM
                      V1(K,3) = XDOT(K,J2) + CM(4,IM)*VREL2(K)/ZM
                      V1(K,4) = XDOT(K,J2) - CM(3,IM)*VREL2(K)/ZM

                      XS(K,I1) = X1(K,3)
                      XS(K,I2)  = X1(K,4)
                      VS(K,I1) = V1(K,3)
                      VS(K,I2)  = V1(K,4)
                      XS(K,ICM2) = X(K,J2)
                      VS(K,ICM2) = XDOT(K,J2)
                  END DO
              END IF
          END IF
   95 CONTINUE
*
*       Check modification for chain regularization (case NAME(ICM) = 0).
      IF (NCH.GT.0) THEN
          CALL CHDATA(XJ,VJ,BODYJ)
          DO 98 L = 1,NCH
*       Copy global address from common JLIST (set in CHDATA).
              J = JLIST(L)
              BODYS(J) = BODYJ(L)
              DO 97 K = 1,3
                  XS(K,J) = XJ(K,L)
                  VS(K,J) = VJ(K,L)
   97         CONTINUE
   98     CONTINUE
      END IF
*
*       Split into WRITE (3) NTOT & WRITE (3) ..  if disc instead of tape.
      IF (FIRST) THEN
          OPEN (UNIT=3,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT3')
          FIRST = .FALSE.
      END IF
      NK = 20
      WRITE (3)  NTOT, MODEL, NRUN, NK
      WRITE (3)  (AS(K),K=1,NK), (BODYS(J),J=1,NTOT),
     &           ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &           (NAME(J),J=1,NTOT)
*     CLOSE (UNIT=3)
*
*       Produce output file for tidal tail members.
   99 IF (KZ(3).LE.3.AND.NTAIL.GT.0) THEN
          IF (SECOND) THEN
             OPEN (UNIT=33,STATUS='NEW',FORM='UNFORMATTED',FILE='OUT33')
             SECOND = .FALSE.
          END IF
          DO 110 I = ITAIL0,NTTOT
              BODYS(I) = BODY(I)
              DO 105 K = 1,3
                  XS(K,I) = X(K,I) - RG(K)
                  VS(K,I) = XDOT(K,I) - VG(K)
  105         CONTINUE
  110     CONTINUE
*       Include cluster centre just in case.
          DO 115 K = 1,3
              AS(K) = RG(K)
              AS(K+3) = VG(K)
              AS(K+6) = RDENS(K)
  115     CONTINUE
          AS(10) = TTOT
          AS(11) = RBAR
          AS(12) = TSCALE
          AS(13) = VSTAR
          NK = 13
          WRITE (33)  NTAIL, NK
          WRITE (33)  (AS(K),K=1,NK), (BODYS(J),J=ITAIL0,NTTOT),
     &                ((XS(K,J),K=1,3),J=ITAIL0,NTTOT),
     &                ((VS(K,J),K=1,3),J=ITAIL0,NTTOT),
     &                (NAME(J),J=ITAIL0,NTTOT)
      END IF
*
*       Include all stars in same file (KZ(3) > 3; astrophysical units). 
      IF (KZ(3).GT.3.AND.NTAIL.GT.0) THEN
          IF (THIRD) THEN
              OPEN (UNIT=34,STATUS='NEW',FORM='FORMATTED',FILE='OUT34')
              THIRD = .FALSE.
          END IF
          NP = 0
*       Copy cluster members with respect to density centre.
          DO 120 I = IFIRST,NTOT
              IF (BODY(I).LE.0.0D0) GO TO 120
              NP = NP + 1
              DO 118 K = 1,3
                  XS(K,NP) = (X(K,I) - RDENS(K))*RBAR
                  VS(K,NP) = XDOT(K,I)*VSTAR
  118         CONTINUE
              BODYS(NP) = BODY(I)*SMU
  120     CONTINUE
          N1 = NP
*       Add tidal tail in the same frame.
          DO 130 I = ITAIL0,NTTOT
              NP = NP + 1
              DO 125 K = 1,3
                  XS(K,NP) = (X(K,I) - RG(K) - RDENS(K))*RBAR
                  VS(K,NP) = (XDOT(K,I) - VG(K))*VSTAR
  125         CONTINUE
              BODYS(NP) = BODY(I)*SMU
  130     CONTINUE
          WRITE (34,140)  NP, N1, (TIME+TOFF)*TSCALE, RBAR, VSTAR,
     &                    (RDENS(K),K=1,3), (RG(K),K=1,3), (VG(K),K=1,3)
  140     FORMAT (' ',2I6,F8.1,2F6.2,3F7.3,1P,6E10.2)
          DO 150 I = 1,NP
              WRITE (34,145) (XS(K,I),K=1,3), (VS(K,I),K=1,3), BODYS(I)
  145         FORMAT (' ',3F10.3,3F8.1,F7.2)
  150     CONTINUE
      END IF
*
*       Update next output interval and initialize the corresponding error.
  100 TNEXT = TNEXT + DELTAT
      ERROR = 0.0D0
*
      RETURN
*
      END
