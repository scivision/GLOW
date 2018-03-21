C Subroutine GCHEM
C
C This software is part of the GLOW model.  Use is governed by the Open Source
C Academic Research License Agreement contained in the file glowlicense.txt.
C For more information see the file glow.txt.
C
C Stan Solomon, 1988, 1989, 1992, 1999, 2005
C
C Includes quartic solution to electron density equation.
C Electron density must be supplied above 200 km in array ZE, a priori
C values are expected but not strictly required at and below 200 km,
C depending on the value of KCHEM.  An initial guess of the N(2D)
C density in array ZND is also expected but not strictly required.
C Other neutral species (O, N2, O2, NO, N(4S)) must be supplied,
C see subroutine GLOW.
C
C Chemical calculations are controlled by switch KCHEM:
C 0 = no calculations at all are performed.
C 1 = electron density, O+(4S), N+, N2+, O2+, NO+ supplied at all
C     altitudes; O+(2P), O+(2D), excited neutrals, and emission rates
C     are calculated.
C 2 = electron density, O+(4S), O2+, NO+ supplied at all altitudes;
C     O+(2P), O+(2D), N+, N2+, excited neutrals, emissions calculated.
C 3 = electron density supplied at all altitudes; everything else
C     calculated.
C 4 = electron density supplied above 200 km; electron density below
C     200 km is calculated, everything else calculated at all altitudes.
C     An initial guess of the electron density below 200 km should be
C     supplied in array ZE.  Electron density for the next two levels
C     above J200 is log interpolated between E(J200) and E(J200+3).
C
C For definitions of common block /CGLOW/ see subroutine GLOW
C
C Other definitions:
C A        Einstein coefficients; s-1
C B        Branching ratios
C BZ       Altitude-dependent branching ratios
C G        Resonant scattering g-factors at each altitude; s-1
C KZ       Temperature dependent rate coeffs at each altitude; cm3s-1
C OEI      O electron impact ionization rates; cm-3s-1
C O2EI     O2   "       "        "       "      "
C RN2EI    N2   "       "        "       "      "
C O2PI     O2 photoionization rates; cm-3s-1
C RN2PI    N2       "          "      "
C RN2ED    N2 electron impact dissociation rate; cm-3s-1
C SRCED    O2     "      "         "         " (SR continuum); cm-3s-1
C P        Volume production rate for each species, altitude; cm-3s-1
C L        Loss rate for each species, altitude; s-1
C T1       Effective temperature divided by 300 for O+ + N2; K
C T2          "          "        "                 O+ + O2; K
C T3          "          "        "                 N2+ + O; K
C T4          "          "        "                 N2+ + O2; K
C T5          "          "        "                 O+ + NO; K
C QQ, RR, SS, TT, UU, VV, WW, XX:  Combined terms for calculation of
C                              O+(4S) given e
C AA, BB, CC, DD, EE, FF, GG, HH: Combined terms for solution of
C                        electron density equation, Roble & Ridley, 1988
C COEF     Coefficients of quartic electron density equation, see "
C ROOT     Roots of          "       "         "       "
C
C
C References for rate coefficients, transition coefficients, branching
C ratios, and g-factors:
C
C      k1  O+ + N2                 Pavlov, 2014
C      k2  O+ + O2                 Pavlov, 2014
C      k3  N2+ + O -> NO+ + N(2D)  Strickland, 1999
C      k4  N2+ + O2                Sinnhuber et al, 2012
C      k5  N(2D) + O2              T and S, 2014
C      k6  N(2D) + O -> N + O      Pandya and Joshipura, 2014
C      k7  N(2D) + et               Thirupathaiah and Singh, 2014
C      k8  O(1D) + N2              Thirupathaiah and Singh, 2014
C      k9  O(1D) + O2              Thirupathaiah and Singh, 2014
C      k10 O(1D) + e -> O + e-     Bhardwaj and Raghuram, 2012
C      k11 O(1S) + O               Strickland, 1999
C      k12 O+(2D) + N2 -> N2+ + O  Pavlov, 2014
C      k13 O+(2D) + O2             Pavlov, 2014
C      k14 O+(2D) + e              Pavlov, 2014
C      k15 O+(2D) + O              Pavlov, 2014
C      k16 O+(2P) + N2 -> O + N2+  Pavlov, 2014
C      k17 O+(2P) + O2             Pavlov, 2014
C      k18 O+(2P) + e              Pavlov, 2014
C      k19 O+(2P) + O              Pavlov, 2014
C      k20 O2(c) + O               Solheim & Llewellyn, 1979?
C      k21 O2(c) + N2              Solheim & Llewellyn, 1979?
C      k22 NO+ + et                Pavlov, 2014
C      k23 N2+ + e                 Pavlov, 2014 - Sinnhuber,2012
C      k24 O2+ + et                Pavlov, 2014
C      k25 N+ + O2                 Pavlov, 2014
C      k26 N2(A) + O               Zettergren, 2009
C      k27 O(1D) + O               Thirupathaiah and Singh, 2014
C      k28 O + et                  Pavlov, 2014
C      k29 N2(A) + O2              Strickland et al, 1999
C      k30 O2+ + NO                Pavlov, 2014
C      k31 N(2D) + NO              Strickland et al, 1999
C      k32 N+ + O                  Pavlov, 2014
C      k33 N(2P) + O               Pandya and Joshipura, 2014
C      k34 N(2P) + O2              Pandya and Joshipura, 2014
C      k35 N(2P) + NO              Strickland et al, 1999
C      k36 O(1S) + O2              Strickland et al, 1999
C      k37 O2+ + N                 Pavlov, 2014
C      k38 O+ + N(2D)              Strickland et al, 1999
C      k39 N2+ + O -> N2 + O+      Strickland et al, 1999
C      k40 O+ + NO                 Sinnhuber, 2012
C      k41 N2+ + NO                Pavlov, 2014
C      k42 N+ + NO                 Pavlov, 2014
C      k43 O+(2D) + NO             Strickland, 1999
C      k44 N(2P) + N2              Pandya and Joshipura, 2014
C      k45 N(2P) + N               Strickland et al, 1999
C      k46 N(2D) + N2              Strickland et al, 1999
C      k47 N + O2                  Strickland et al, 1999**
C      k48 N + NO                  Strickland et al, 1999**
C      k49 N + O                   Strickland et al, 1999**
C      k50 O(1S) + NO              Strickland et al, 1999
C      k51 O(1D) + NO              Strickland et al, 1999
C      k52 N2(A) + NO              Strickland et al, 1999
C      k53 N2(A) + O               Zettergren, 2009
C      k54 N2(A) + N               Strickland et al, 1999
C      k55 N(2P) + e -> N + e      Strickland et al, 1999
C      k56 N(2P) + e -> N(2D) + e  Strickland et al, 1999
C      k57 O(1S) + e -> O + e      Bhardwaj and Raghuram, 2012
C      k58 O(1S) + e -> O(1D) + e  Bhardwaj and Raghuram, 2012
C      k59 O + e -> O- + hv        Pavlov, 2014
C      k60 O+ + e -> O + hv(1304)  Strickland et al, 1999
C      k61 O+ + e -> O + hv(1356)  Strickland et al, 1999
C      k62 O(1D) + e -> O+(4S)+2e  Bhardwaj and Raghuram, 2012
C      k63 O+(4S) + O- -> O + O*   Strickland et al, 1999
C      k64 O- + O -> O2 + e        Pavlov, 2014
C      k65 O+(2D) + N              Sinnhuber, 2012
C      k66 O+(2P) + N              Sinnhuber, 2012
C      k67 O+(2P) + NO             Strickland, 1999
C      k68 O+(2P) + N2 -> N+ + NO  Sinnhuber, 2012
C      k69 N+ + e -> N + hv        Schunk and Nagy, 2009
C      k70 O+(2D) + N2->O+(4S)+N2  Sinnhuber, 2012
C      k71 O2+ + N2 -> NO+ + NO    Rees, 1989
C      k72 N(2D) + O -> NO+ + e    Strickland et al, 1999
C      k73 O2+ + N(2D) -> N+ + O2  Rees, 1989
C      k74 O+(2P) + N2 -> O+ + N2* Dharwan et al, 2014
C      k75 O+ + e                  Schunk and Nagy, 2009

C      A1  5200       N(4S-2D)     Pandya and Joshipura, 2014
C      A2  6300       O(3P-1D)     Bhardwaj and Raghuram, 2012
C      A3  6364       O(3P-1D)     Bhardwaj and Raghuram, 2012
C      A4  2972       O(3P-1S)     T and S, 2014
C      A5  5577       O(1D-1S)     Pandya and Joshipura, 2014
C      A6  3726       O+(4S-2D)    Strickland et al, 1999
C      A7  2470       O+(4S-2P)    Zettergren, 2009
C      A8  7319-30    O+(2D-2P)    Zettergren, 2009
C      A9  (Hertz II) O2(X-c)      Strickland, 1999
C      A10 (Veg-Kap)  N2(X-A)      Strickland, 1999
C      A11 3466       N(4S-2P)     Pandya and Joshipura, 2014
C      A12 10400      N(2D-2P)     Pandya and Joshipura, 2014
C      B1  O(1S) from O2+ + e      Strickland et al, 1999
C      B2  O(1D) from O2+ + e      Strickland et al, 1999
C      B3  N(2D) from NO+ + e      Strickland et al, 1999
C      B4  N(2D) from N2+ + e      Strickland et al, 1999
C      B5  N(2D) from N2+ + O      Strickland et al, 1999
C      B6  O(1D) from N(2D) + O2   Pavlov, 2014
C      B7  O(1D) from O+(2D) + O   Pavlov, 2014
C      B8  O+(2D) from O+(2P) + e  Pavlov, 2014
C      B9  O+(2P) from O + e*      Strickland et al, 1999
C      B10 O+(2D) from O + e*      Strickland et al, 1999
C      B11 O+ from O + e*      Strickland et al, 1999
C      B12 O+(2P) from O2 + e*     Strickland et al, 1999
C      B13 O+(2D) from O2 + e*     Strickland et al, 1999
C      B14 O+ from O2 + e*     Strickland et al, 1999
C      B15 N+ from N2 + e*         Strickland et al, 1999
C      B16 N(2D) from above        Strickland et al, 1999
C      B17 O(1D) from N+ + O2      Pavlov, 2014
C      B18 O(1S) from N2(A) + O    Pavlov, 2014
C      B19 O(1S) from O2(*) + O    ?
C      B20 O(1D) from N(2D) + O    Pandya and Joshipura, 2014
C      B21 NO+ from N+ + O2        Pavlov, 2014
C      B22 O2+ from N+ + O2        Pavlov, 2014
C      B23 N(2P) from N2+ + e      Pavlov, 2014
C      B24 N2 + e* pri. dis/ion    xsect ratio > 250 eV?
C      B25 N(2D) from N2 + e* dis  Zipf et al, 1980?
C      B26 N(2P) from N2 + e* dis  Zipf et al, 1980?
C      B27 N(2D) from N2 + hv      Richards et al, 1981 (add to B28)?
C      B28 N(2P) from N2 + hv      ? (cf Zipf & McGlaughlin, 1978)?
C      B29 N(2D) from N(2P) + O    Pandya and Joshipura, 2014
C      B30 O+(2P) from O2 + hv     ?
C      B31 O+(2D)  "               ?
C      B32 O+  "               ?
C      B33 O(1S) from O2+ + N      Zettergren, 2009
C      B34 O(1S) from N(2D) + NO   .1 approximation
C      B35 O2 + e* pri SRC/ion     xsect ratio > 250 eV?
C      B36 N2+(B) from N2 + e*     Strickland et al, 1999
C      B37 (0,0) (3914) fr. N2+(B) Shemansky & Broadfoot, 1971?
C      B38 (0,1) (4278) fr. N2+(B) Shemansky & Broadfoot, 1971?
C      B39 (0,0) (3371) fr. N2(C)  Conway, 1983; Benesch et al, 1966?
C      B40 (0,9) (3352) fr. N2(A)  Cartwright, 1978; Shemansky, 1969?
C      B41 O+(2Po) fr. O+(2Pe)     Kirby et al, 1979?
C      B42 O+(2Do) fr. O+(2Pe)     Kirby et al, 1979?
C      B43 N2(A) bound fraction    ?
C      B44 7990 fr. O(3s'3D)       appx. fr. Hecht, p.c.?
C      B45 NO+ fr. N+ + NO         Pavlov, 2014
C      B46 O(1D) fr. O(1S) + O2    Strickland, 1999
C      B47 O(1D) fr. O(1S) + NO    Strickland, 1999
C      B48 O(1S) fr. N2(A) + O     Zettergren, 2009
C      B49 N(2P) fr. N2(A) + N     Strickland, 1999
C      B50 O+ fr. N+ + O2      Pavlov, 2014
C      B51 N(2D) fr. N+ + O2       Pavlov, 2014
C      B52 O(1D) fr. N(2P) + O     Pandya and Joshipura, 2014
C      B53 O(3P) fr. O2+ + e       Pavlov, 2014
C      B54 O+ fr. O+(2D) + N2  Strickland et al, 1999
C      G1  N2+B(0,0) (3914)        Broadfoot, 1967?
C      G2  N2+B(0,1) (4278)        Broadfoot, 1967?
C      **All species involved are major, collision not used
C      ? For references not verified during updates
C
C Array dimensions:
C JMAX    number of altitude levels
C NBINS   number of energetic electron energy bins
C LMAX    number of wavelength intervals for solar flux
C NMAJ    number of major species
C NEX     number of ionized/excited species
C NW      number of airglow emission wavelengths
C NC      number of component production terms for each emission
C NST     number of states produced by photoionization/dissociation
C NEI     number of states produced by electron impact
C NR      number of rate coefficients, branching ratios, A and G factors
C NF      number of available types of auroral fluxes
C
C
      SUBROUTINE GCHEM
C
      INCLUDE 'glow.h'
      PARAMETER (NMAJ=3)
      PARAMETER (NEX=20)
      PARAMETER (NW=20)
      PARAMETER (NC=10)
      PARAMETER (NST=6)
      PARAMETER (NEI=10)
      PARAMETER (NF=4)
      PARAMETER (NR=100)
C
      COMMON /CGLOW/
     >    IDATE, UT, GLAT, GLONG, ISCALE, JLOCAL, KCHEM,
     >    F107, F107A, HLYBR, FEXVIR, HLYA, HEIEW, XUVFAC,
     >    ZZ(JMAX), ZO(JMAX), ZN2(JMAX), ZO2(JMAX), ZNO(JMAX),
     >    ZNS(JMAX), ZND(JMAX), ZRHO(JMAX), ZE(JMAX),
     >    ZTN(JMAX), ZTI(JMAX), ZTE(JMAX),
     >    PHITOP(NBINS), EFLUX(NF), EZERO(NF),
     >    SZA, DIP, EFRAC, IERR,
     >    ZMAJ(NMAJ,JMAX), ZCOL(NMAJ,JMAX),
     >    WAVE1(LMAX), WAVE2(LMAX), SFLUX(LMAX),
     >    ENER(NBINS), DEL(NBINS),
     >    PESPEC(NBINS,JMAX), SESPEC(NBINS,JMAX),
     >    PHOTOI(NST,NMAJ,JMAX), PHOTOD(NST,NMAJ,JMAX), PHONO(NST,JMAX),
     >    QTI(JMAX), AURI(NMAJ,JMAX), PIA(NMAJ,JMAX), SION(NMAJ,JMAX),
     >    UFLX(NBINS,JMAX), DFLX(NBINS,JMAX), AGLW(NEI,NMAJ,JMAX),
     >    EHEAT(JMAX), TEZ(JMAX), E(JMAX),
     >    DEN(NEX,JMAX), ZETA(NW,JMAX), ZCETA(NC,NW,JMAX), VCB(NW)
C
      REAL KZ, L
C
      DIMENSION A(NR), B(NR), BZ(NR,JMAX), G(NR,JMAX), KZ(NR,JMAX),
     >          OEI(JMAX), O2EI(JMAX), RN2EI(JMAX), O2PI(JMAX),
     >          RN2PI(JMAX), RN2ED(JMAX), SRCED(JMAX),
     >          P(NEX,JMAX), L(NEX,JMAX),
     >          T1(JMAX), T2(JMAX), T3(JMAX), T4(JMAX), T5(JMAX),
     >          QQ(JMAX), RR(JMAX), SS(JMAX), TT(JMAX), UU(JMAX),
     >          VV(JMAX), WW(JMAX), XX(JMAX),
     >          AA(JMAX), BB(JMAX), CC(JMAX), DD(JMAX), EE(JMAX),
     >          FF(JMAX), GG(JMAX), HH(JMAX)
C
      DOUBLE PRECISION COEF(JMAX,5), ROOT(JMAX)
C
      DATA RE/6.37E8/
      DATA A / 5.765E-06, 6.478E-03, 2.097E-03, 0.134, 1.260,
     >         8.89E-05, 4.70E-02, 0.1740, .660, 0.3520,
     >         6.50E-03, 3.45E-03, 88*0.0 /
      DATA B/.002, 1.389, .76, 1.12, 1.0, .10, .10, .737, .100, .196,
     >       .283, .033, .065, .094, .10, 0.288, .36, .36, 0.00, .10,
     >       0.45, .50, 0.00, 0.60, 0.54, 0.44, 0.80, 0.20, .10, 0.33,
     >       0.33, .34, .15, 0.10, 0.10, .11, 0.65, 0.2, 0.24, 0.02,
     >       0.18, 0.72, 0.75, 0.10, .89, .31, .64, .36, .9, .05,
     >       .15, .10, .609, .041, 46*0./
C
C
      IF (KCHEM .EQ. 0) RETURN
C
C
C Zero airglow and density arrays:
C
      DO 10 I=1,JMAX
      DO 10 IW=1,NW
      ZETA(IW,I) = 0.0
   10 CONTINUE
C
      DO 20 I=1,JMAX
      DO 20 IW=1,NW
      DO 20 IC=1,NC
      ZCETA(IC,IW,I) = 0.0
   20 CONTINUE
C
      IF (KCHEM .GE. 3) THEN
      DO 30 I=1,JMAX
      DO 30 IX=1,NEX
      DEN(IX,I)=0.
   30 CONTINUE
      ENDIF
C
C
C Assign g-factors at altitudes which are sunlit:
C
      DO 40 I=1,JMAX
      DO 40 N=1,NR
      G(N,I) = 0.0
   40 CONTINUE
C
      DO 50 I=1,JMAX
      GH = (RE+ZZ(I)) * SIN(SZA)
      IF (SZA .LT. 1.6 .OR. GH .GT. RE) THEN
        G(1,I) = 0.041
        G(2,I) = 0.013
      ENDIF
   50 CONTINUE
C
C
C Calculate rate coefficients as a function of altitude:
C
      DO 70 I=1,JMAX
      T1(I) = (16.*ZTN(I)+28.*ZTI(I)) / (16.+28.) / 300.
      T2(I) = (16.*ZTN(I)+32.*ZTI(I)) / (16.+32.) / 300.
      T3(I) = (28.*ZTN(I)+16.*ZTI(I)) / (28.+16.) / 300.
      T4(I) = (28.*ZTN(I)+32.*ZTI(I)) / (28.+32.) / 300.
      T5(I) = (16.*ZTN(I)+30.*ZTI(I)) / (16.+30.) / 300.
      KZ(1,I) = 1.0E-12*(2.05 - .00308*ZTN(I))
      IF (ZTN(I) .GE. 300) KZ(1,I)=1.72E-12-7.2E-13*(ZTN(I)/300.)
     >                    +1.33E-13*(ZTN(I)/300.)**2.
     >                    +9.28E-15*(ZTN(I)/300.)**3.
     >                    +6.40E-16*(ZTN(I)/300.)**4.
      KZ(2,I) = 1.6E-11 * (300./ZTN(I))**.52
     >                    + 5.5E-11 * EXP(-6382./ZTN(I))
      KZ(3,I) = 1.4E-10*(T3(I)**(-.44))*.07*(T3(I)**(.21)) *
     >                  (1 - .07*T3(I)**.21)
      IF(T3(I)*300 .GT. 1500) KZ(3,I) = 5.2E-11*(T3(I)**(.2)) *
     >                  (1 - .07*T3(I)**.21)
      KZ(4,I) = 5.0E-11 * (1./T4(I)) ** 0.8
      KZ(5,I) = 9.7E-12*EXP(-185/ZTN(I))
      KZ(6,I) = 5.3E-13
      KZ(7,I) = 3.8E-12 * (ZTE(I)) ** 0.81
      KZ(8,I) = 1.8E-11 * EXP(107.8/ZTN(I))
      KZ(9,I) = 3.2E-11 * EXP(67.0/ZTN(I))
      KZ(10,I) = 8.1E-10 * (300./ZTE(I)) ** 0.5
      KZ(11,I) = 2.0E-14
      KZ(12,I) = 1.5E-10 * (300./ZTN(I)) ** 0.5
      KZ(13,I) = 1.0E-10 * (300./ZTN(I)) ** 0.5
      KZ(14,I) = 4.0E-08 * (300./ZTE(I)) ** 0.5
      KZ(15,I) = 5.0E-12
      KZ(16,I) = 2.0E-10 * (300./ZTN(I)) ** 0.5
      KZ(17,I) = 3.1E-10 * (300./ZTN(I)) ** 0.5
      KZ(18,I) = 9.50E-08 * (300./ZTE(I)) ** 0.5
      KZ(19,I) = 5.2E-11
      KZ(20,I) = 8.54E-11 * EXP(-3408./ZTN(I))
     >                    + 1.7E-11*(ZTN(I)/300.)**(-.7)
      KZ(21,I) = 5.0E-11*(ZTN(I)/300.)**(-.8)
      KZ(22,I) = 3.5E-07 * (300./ZTE(I)) ** 0.69
      IF (ZTE(I) .GE. 1200.) KZ(22,I) = 3.02E-07*(300./ZTE(I))**0.56
      KZ(23,I) = 2.2E-07 * (300./ZTE(I)) ** 0.39
      IF (ZTE(I) .GE. 1200.) KZ(23,I) = 1.95E-07*(300./ZTE(I))**0.57
      KZ(24,I) = 1.95E-07 * (300./ZTE(I)) ** 0.70
      IF (ZTE(I) .GE. 1200.) KZ(24,I) = 1.93E-07*(300./ZTE(I))**0.61
      KZ(25,I) = 5.5E-10
      KZ(26,I) = 2.8E-11
      KZ(27,I) = 2.5E-11
      KZ(28,I) = 1.384E-15*EXP(-1.759E-04*ZTE(I)+8.56E-08*ZTE(I)**2.
     >                  - 1.427E-11 * ZTE(I)**3.)
      KZ(29,I) = 4.0E-12
      KZ(30,I) = 4.1E-10
      KZ(31,I) = 6.7E-11
      KZ(32,I) = 2.2E-12
      KZ(33,I) = 2.7E-11
      KZ(34,I) = 2.5E-12
      KZ(35,I) = 3.0E-11
      KZ(36,I) = 2.32E-12*EXP(-(6750.-.015*(ZTN(I)**2.))/(8.314*ZTN(I)))
      KZ(37,I) = 1.65E-10
      KZ(38,I) = 1.3E-10
      KZ(39,I) = 1.4E-10*(T3(I)**(-.44))*.07*(T3(I)**(.21)) *
     >                  .07*T3(I)**.21
      IF(T3(I)*300 .GT. 1500) KZ(39,I) = 5.2E-11*(T3(I)**(.2)) *
     >                  .07*T3(I)**.21
      KZ(40,I) = 8.0E-13
      KZ(41,I) = 7.5E-09*((1./ZTN(I))**.52)
      KZ(42,I) = 6.44E-09*((1./ZTN(I))**.44)
      KZ(43,I) = 1.2E-09
      KZ(44,I) = 5.0E-17
      KZ(45,I) = 6.0E-13
      KZ(46,I) = 1.0E-13*EXP(-(510./ZTN(I)))
      KZ(47,I) = 1.5E-11*EXP(-(3573./ZTN(I)))
      KZ(48,I) = 2.2E-11*EXP(-(160./ZTN(I)))
      IF(ZTN(I) .GT. 400) KZ(48,I) = 3.3E-11
      KZ(49,I) = 3.33E-16*((1./ZTN(I))**.5)*(1.0-.567*(1./ZTN(I))**.5)
      KZ(50,I) = 8.0E-11
      KZ(51,I) = 1.5E-10
      KZ(52,I) = 8.9E-11
      KZ(53,I) = 2.8E-11
      KZ(54,I) = 4.0E-11
      KZ(55,I) = 1.6E-12*(ZTE(I)**.85)
      KZ(56,I) = 9.5E-09
      KZ(57,I) = 1.56E-07*(ZTE(I)/300.)**.74
      KZ(58,I) = 8.56E-09
      KZ(59,I) = 1.384E-15*EXP(-1.759E-04*ZTE(I)+8.56E-08*(ZTE(I)**2.)
     >                  -1.427E-11*(ZTE(I)**3.))
      KZ(60,I) = 0.
C     ^Removed for now: 3.4E-13*(1160./ZTE(I))**.5
      KZ(61,I) = 0.
C     ^Removed for now: 6.3E-13*(1160./ZTE(I))**.5
      KZ(62,I) = 0.
C     ^Removed due to photo electron: 1.75E-07
      KZ(63,I) = 1.0E-7
      KZ(64,I) = 2.3E-10
      KZ(65,I) = 7.5E-11
      KZ(66,I) = 1.0E-10
      KZ(67,I) = 2.9E-08
      KZ(68,I) = 1.0E-10
      KZ(69,I) = 3.6E-12*((250./ZTE(I))**.7)
      KZ(70,I) = 8.0E-10
      KZ(71,I) = 5.0E-16
      KZ(72,I) = 2.5E-18*(ZTN(I)**.5)*(2205+ZTN(I))*EXP(-4410./ZTN(I))
      KZ(73,I) = 2.5E-10
      KZ(74,I) = 0.
C     ^Removing due to lack of definition of N2*: 4.8E-10
      KZ(75,I) = 3.7E-12*((250/ZTE(I))**.7)
   70 CONTINUE
C
C
C Calculate Electron impact ionization, photoionization, and electron
C impact dissociation rates at each altitude; put a priori electron
C density in calculated electron density array: put rough estimate of
C O+ and a priori N(2D) in DEN array:

C
      DO 100 I=1,JMAX
      OEI(I)   = SION(1,I)
      O2EI(I)  = SION(2,I)
      RN2EI(I) = SION(3,I)
      O2PI(I)  = PHOTOI(1,2,I) + PHOTOI(2,2,I) + PHOTOI(3,2,I)
      RN2PI(I) = PHOTOI(1,3,I) + PHOTOI(2,3,I) + PHOTOI(3,3,I) +
     >           PHOTOI(4,3,I) + PHOTOI(5,3,I)
      RN2ED(I) = AGLW(5,3,I) + AGLW(6,3,I) + AGLW(7,3,I)
      SRCED(I) = AGLW(4,2,I)
      E(I)     = ZE(I)
      DEN(10,I)= ZND(I)
  100 CONTINUE
C
C
C Find level below which electron density will be calculated:
C
      IF (KCHEM .GE. 4) THEN
        DO 120 I=JMAX,1,-1
        IF (ZZ(I) .GT. 2.0001E7) J200=I-1
  120   CONTINUE
      ELSE
        J200=0
      ENDIF
C
C
C Iterative loop assures that feedback reactions (O+(2P,2D)+e,
C O+(4S)+N(2D), N2++O) are correctly computed:
C
      DO 350 ITER=1,3
C
C
C Calculate atomic ion densities at each altitude:
C
      DO 150 I=1,JMAX
C
C
C O+(2P):
C
      P(1,I)= PHOTOI(3,1,I)
     >      + B(41) * PHOTOI(5,1,I)
     >      + B(30) * PHOTOI(4,2,I)
     >      + B(9)  * OEI(I)
     >      + B(12) * O2EI(I)
      L(1,I)= KZ(16,I) * ZN2(I)
     >      + KZ(17,I) * ZO2(I)
     >      + KZ(19,I) * ZO(I)
     >      + KZ(18,I) * E(I)
     >      + KZ(66,I) * ZNS(I)
     >      + KZ(67,I) * ZNO(I)
     >      + KZ(68,I) * ZN2(I)
     >      + KZ(74,I) * ZN2(I)
     >      + A(8)
     >      + A(7)
      DEN(1,I) = P(1,I) / L(1,I)
C
C
C O+(2D):
C
      P(2,I)= PHOTOI(2,1,I)
     >      + B(42) * PHOTOI(5,1,I)
     >      + B(31) * PHOTOI(4,2,I)
     >      + B(10) * OEI(I)
     >      + B(13) * O2EI(I)
     >      + B(8)  * KZ(18,I) * DEN(1,I) * E(I)
     >      + KZ(19,I) * DEN(1,I) * ZO(I)
     >      + A(8)  * DEN(1,I)
      L(2,I)= KZ(12,I) * ZN2(I) 
     >      + KZ(13,I) * ZO2(I) 
     >      + KZ(15,I) * ZO(I)
     >      + KZ(14,I) * E(I)
     >      + KZ(43,I) * ZNO(I)
     >      + KZ(65,I) * ZNS(I)
     >      + KZ(70,I) * ZN2(I)
     >      + A(6)
      DEN(2,I) = P(2,I) / L(2,I)
C
C
C O+(4S):
C
      IF (KCHEM .GE. 3) THEN
        P(3,I)= PHOTOI(1,1,I) + PHOTOI(4,1,I)
     >        + B(32) * PHOTOI(4,2,I)
     >        + B(11) * OEI(I)
     >        + B(14) * O2EI(I) 
     >        + KZ(14,I) * DEN(2,I) * E(I) 
     >        + KZ(15,I) * DEN(2,I) * ZO(I) 
     >        + A(6) * DEN(2,I)
     >        + (1.-B(8)) * KZ(18,I) * DEN(1,I) * E(I)
     >        + KZ(19,I) * DEN(1,I) * ZO(I) 
     >        + A(7) * DEN(1,I)
     >        + KZ(32,I) * DEN(4,I) * ZO(I)
     >        + KZ(39,I) * DEN(5,I) * ZO(I)
     >        + B(50) * KZ(25,I) * DEN(4,I) * ZO2(I)
     >        + B(54) * KZ(12,I) * DEN(2,I) * ZN2(I)
C    Not used in AURIC:
     >        + KZ(62,I) * DEN(12,I) * E(I)
     >        + KZ(70,I) * DEN(1,I) * ZN2(I)
     >        + KZ(74,I) * DEN(1,I) * ZN2(I)
        L(3,I)= KZ(1,I) * ZN2(I)
     >        + KZ(2,I) * ZO2(I)
     >        + KZ(38,I) * DEN(10,I)
     >        + KZ(40,I) * ZNO(I)
     >        + KZ(75,I) * E(I)
C    Not used in AURIC:
C     >        + 1.98 * KZ(60,I) * E(I)
C     >        + 1.93 * KZ(61,I) * E(I)
        DEN(3,I) = P(3,I) / L(3,I)
      ENDIF
C
C
C N+:
C
      IF (KCHEM .GE. 2) THEN
        P(4,I) = PHOTOI(6,3,I)
     >         + B(15) * RN2EI(I)
     >         + KZ(38,I) * DEN(3,I) * DEN(10,I)
     >         + KZ(65,I) * DEN(2,I) * ZNS(I)
     >         + KZ(66,I) * DEN(1,I) * ZNS(I)
     >         + KZ(68,I) * DEN(1,I) * ZN2(I)
     >         + KZ(73,I) * DEN(10,I) * DEN(6,I)
        L(4,I) = KZ(25,I) * ZO2(I)
     >         + KZ(32,I) * ZO(I)
     >         + KZ(42,I) * ZNO(I)
     >         + KZ(69,I) * E(I)
        DEN(4,I) = P(4,I) / L(4,I)
      ENDIF
C
  150 CONTINUE
C
C
C
C Above 200 km, (or at all altitudes if KCHEM=3) use a priori
C electron density to calculate O+(4S):
C
      IF (KCHEM .GE. 3) THEN
C
      DO 170 I=J200+1,JMAX
      P(5,I)= RN2PI(I)
     >      + (1.-B(15)) * RN2EI(I)
     >      + KZ(12,I) * DEN(2,I) * ZN2(I)
     >      + KZ(16,I) * DEN(1,I) * ZN2(I)
     >      + (1.-B(45))*KZ(42,I) * DEN(4,I) * ZNO(I)
      L(5,I)= KZ(3,I)  * ZO(I)
     >      + KZ(4,I)  * ZO2(I)
     >      + KZ(23,I) * E(I)
     >      + KZ(39,I) * ZO(I)
     >      + KZ(41,I) * ZNO(I)
      DEN(5,I) = P(5,I) / L(5,I)
      QQ(I) = PHONO(1,I)
     >      + KZ(3,I)  * DEN(5,I) * ZO(I)
     >      + B(21) * KZ(25,I) * DEN(4,I) * ZO2(I)
      RR(I) = KZ(30,I) * ZNO(I)
     >      + KZ(37,I) * ZNS(I)
      SS(I) = KZ(1,I) * ZN2(I)
     >      + KZ(40,I) * ZNO(I)
      TT(I) = KZ(22,I) * E(I)
      UU(I) = O2PI(I)
     >      + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >      + KZ(13,I) * DEN(2,I) * ZO2(I)
     >      + KZ(17,I) * DEN(1,I) * ZO2(I)
     >      + KZ(4,I)  * DEN(5,I) * ZO2(I)
     >      + B(22) * KZ(25,I) * DEN(4,I) * ZO2(I)
      VV(I) = KZ(2,I) * ZO2(I)
      WW(I) = KZ(24,I) * E(I)
     >      + KZ(30,I) * ZNO(I)
     >      + KZ(37,I) * ZNS(I)
      XX(I) = DEN(1,I) + DEN(2,I) + DEN(4,I) + DEN(5,I)
      DEN(3,I) = (TT(I)*WW(I)*E(I)-TT(I)*WW(I)*XX(I)-TT(I)*UU(I)
     >            - QQ(I)*WW(I) - RR(I)*UU(I) ) /
     >           (TT(I)*WW(I)+TT(I)*VV(I)+RR(I)*VV(I)+SS(I)*WW(I))
  170 CONTINUE
C
      ENDIF
C
C
C If KCHEM=4, calculate electron density using quartic equation method
C below 200 km:
C
      IF (KCHEM .GE. 4) THEN
C
      DO 200 I=1,J200
      AA(I) =  PHONO(1,I)
     >       + KZ(1,I)  * DEN(3,I) * ZN2(I)
     >       + KZ(40,I) * DEN(3,I) * ZNO(I)
     >       + B(21) * KZ(25,I) * DEN(4,I) * ZO2(I)
      BB(I) =  O2PI(I)
     >       + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >       + KZ(2,I)  * DEN(3,I) * ZO2(I)
     >       + KZ(13,I) * DEN(2,I) * ZO2(I)
     >       + KZ(17,I) * DEN(1,I) * ZO2(I)
     >       + B(22) * KZ(25,I) * DEN(4,I) * ZO2(I)
      CC(I) =  KZ(30,I) * ZNO(I)
     >       + KZ(37,I) * ZNS(I)
      DD(I) =  RN2PI(I)
     >       + (1.-B(15)) * RN2EI(I)
     >       + KZ(12,I) * DEN(2,I) * ZN2(I)
     >       + KZ(16,I) * DEN(1,I) * ZN2(I)
      EE(I) = KZ(3,I) * ZO(I)
      FF(I) = DEN(1,I) + DEN(2,I) + DEN(3,I) + DEN(4,I)
      HH(I) = KZ(4,I)  * ZO2(I)
      GG(I) = KZ(39,I) * ZO(I)
      COEF(I,5) = KZ(22,I) * KZ(23,I) * KZ(24,I)
      COEF(I,4) =(KZ(22,I) * (KZ(24,I)*(EE(I)+HH(I)+GG(I))
     >                        + KZ(23,I)*CC(I))
     >            - KZ(22,I) * KZ(23,I) * KZ(24,I) * FF(I)) / 4.0
      COEF(I,3) =(KZ(22,I) * CC(I) * (EE(I)+HH(I)+GG(I))
     >            - KZ(22,I) * FF(I)
     >                       * (KZ(24,I)*(EE(I)+HH(I)+GG(I))
     >                          + KZ(23,I)*CC(I))
     >            - KZ(24,I) * KZ(23,I) * AA(I)
     >            - KZ(22,I) * KZ(23,I) * BB(I)
     >            - KZ(22,I) * KZ(24,I) * DD(I)) / 6.0
      COEF(I,2) =(- KZ(22,I) * (CC(I)*FF(I)*(EE(I)+HH(I)+GG(I))
     >                         + DD(I)*HH(I) + BB(I)*(EE(I)+HH(I)+GG(I))
     >                          + CC(I)*DD(I))
     >            - KZ(24,I) * (AA(I)*(EE(I)+HH(I)+GG(I)) + DD(I)*EE(I))
     >            - KZ(23,I) * CC(I) * (AA(I)+BB(I))) / 4.0
      COEF(I,1) = - CC(I) * (EE(I)+HH(I)+GG(I)) * (AA(I)+BB(I)+DD(I))
      COEF(I,5) = COEF(I,5) * 1.D10
      COEF(I,4) = COEF(I,4) * 1.D10
      COEF(I,3) = COEF(I,3) * 1.D10
      COEF(I,2) = COEF(I,2) * 1.D10
      COEF(I,1) = COEF(I,1) * 1.D10
  200 CONTINUE
C
      CALL VQUART (COEF, ROOT, J200)
C
      DO 250 I=1,J200
      E(I) = ROOT(I)
  250 CONTINUE
C
      E(J200+1) = E(J200) * ( E(J200+3) / E(J200) )
     >            ** ( (ZZ(J200+1)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
      E(J200+2) = E(J200) * (E(J200+3)/E(J200))
     >            ** ( (ZZ(J200+2)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
C
      ENDIF
C
C
C Calculate molecular ion densities and excited species densites:
C
      DO 300 I=1,JMAX
C
C
C N2+:
C
      IF (KCHEM .GE. 2) THEN
        P(5,I)= RN2PI(I)
     >        + (1.-B(15)) * RN2EI(I)
     >        + KZ(16,I) * DEN(1,I) * ZN2(I)
     >        + (1.-B(45))*KZ(42,I) * DEN(4,I) * ZNO(I)
        L(5,I)= KZ(3,I)  * ZO(I)
     >        + KZ(4,I)  * ZO2(I)
     >        + KZ(23,I) * E(I)
     >        + KZ(39,I) * ZO(I)
     >        + KZ(41,I) * ZNO(I)
          DEN(5,I) = P(5,I) / L(5,I)
      ENDIF
C
C
C O2+:
C
      IF (KCHEM .GE. 3) THEN
        P(6,I)= O2PI(I)
     >        + (1.-B(12)-B(13)-B(14)) * O2EI(I)
     >        + KZ(2,I)  * DEN(3,I) * ZO2(I)
     >        + KZ(13,I) * DEN(2,I) * ZO2(I)
     >        + KZ(17,I) * DEN(1,I) * ZO2(I)
     >        + KZ(4,I)  * DEN(5,I) * ZO2(I)
     >        + B(22) * KZ(25,I) * DEN(4,I) * ZO2(I)
        L(6,I)= KZ(24,I) * E(I)
     >        + KZ(30,I) * ZNO(I)
     >        + KZ(37,I) * ZNS(I)
     >        + KZ(71,I) * ZN2(I)
     >        + KZ(73,I) * DEN(10,I)
        DEN(6,I) = P(6,I)/ L(6,I)
      ENDIF
C
C
C NO+:
C
      IF (KCHEM .GE. 3) THEN
        P(7,I)= PHONO(1,I)
     >        + KZ(1,I)  * DEN(3,I) * ZN2(I)
     >        + KZ(40,I) * DEN(3,I) * ZNO(I)
     >        + KZ(3,I)  * DEN(5,I) * ZO(I)
     >        + KZ(41,I) * DEN(5,I) * ZNO(I)
     >        + B(21) * KZ(25,I) * DEN(4,I) * ZO2(I)
     >        + KZ(30,I) * DEN(6,I) * ZNO(I)
     >        + KZ(37,I) * DEN(6,I) * ZNS(I)
     >        + B(45)*KZ(42,I)*DEN(4,I)*ZNO(I)
     >        + KZ(43,I) * DEN(2,I) * ZNO(I)
     >        + KZ(67,I) * DEN(1,I) * ZNO(I)
     >        + KZ(71,I) * DEN(6,I) * ZN2(I)
     >        + KZ(72,I) * DEN(10,I) * ZO(I)
        L(7,I)= KZ(22,I) * E(I)
        DEN(7,I) = P(7,I) / L(7,I)
      ENDIF
C
C
C N2(A):
C
      P(8,I)= AGLW(1,3,I) + AGLW(2,3,I) + B(43)*AGLW(3,3,I)
     >        + KZ(74,I) * DEN(1,I) * ZN2(I)
      L(8,I)= KZ(26,I) * ZO(I)
     >      + KZ(29,I) * ZO2(I)
     >      + KZ(52,I) * ZNO(I)
     >      + KZ(54,I) * ZNS(I)
     >      + A(10)
      DEN(8,I) = P(8,I) / L(8,I)
C
C
C N(2P):
C
      P(9,I)= B(28) * PHOTOD(1,3,I)
     >      + B(28) * PHOTOI(6,3,I)
     >      + B(26) * RN2ED(I)
C    Not used in AURIC:
     >      + B(23) * KZ(23,I) * DEN(5,I) * E(I)
     >      + B(49) * KZ(54,I) * DEN(8,I) * ZNS(I)
      L(9,I)= KZ(33,I) * ZO(I)
     >      + KZ(34,I) * ZO2(I)
     >      + KZ(35,I) * ZNO(I)
     >      + KZ(44,I) * ZN2(I)
     >      + KZ(45,I) * ZNS(I)
     >      + KZ(55,I) * E(I)
     >      + KZ(56,I) * E(I)
     >      + A(11)
     >      + A(12)
      DEN(9,I) = P(9,I) / L(9,I)
C
C
C N(2D):
C
      P(10,I)= B(27) * PHOTOD(1,3,I)
     >       + B(27) * PHOTOI(6,3,I)
     >       + B(25) * RN2ED(I)
     >       + B(16) * B(15) * RN2EI(I)
     >       + B(3)  * KZ(22,I) * DEN(7,I) * E(I)
     >       + B(4)  * KZ(23,I) * DEN(5,I) * E(I)
     >       + B(5)  * KZ(3,I)  * DEN(5,I) * ZO(I)
     >       + B(29) * KZ(33,I) * DEN(9,I) * ZO(I)
     >       + KZ(45,I) * DEN(9,I) * ZNS(I)
     >       + KZ(56,I) * DEN(9,I) * E(I)
     >       + B(51) * KZ(25,I) * DEN(4,I) *ZO2(I)
     >       + A(12) * DEN(9,I)
C    Not used in AURIC:
     >       + (1-B(49)) * KZ(54,I) * DEN(8,I) * ZNS(I)
      L(10,I)= KZ(5,I)  * ZO2(I)
     >       + KZ(6,I)  * ZO(I)
     >       + KZ(7,I)  * E(I)
     >       + KZ(31,I) * ZNO(I)
     >       + KZ(38,I) * DEN(3,I)
     >       + KZ(46,I) * ZN2(I)
     >       + KZ(73,I) * DEN(6,I)
     >       + A(1)
      DEN(10,I) = P(10,I) / L(10,I)
C
C
C O(1S):
C
      P(11,I)= AGLW(2,1,I)
     >       + PHOTOD(2,2,I)
     >       + B(1) * KZ(24,I) * DEN(6,I)  * E(I)
     >       + B(18) * KZ(26,I) * DEN(8,I) * ZO(I)
C    Not used in AURIC:
     >       + B(33) * KZ(37,I) * DEN(6,I) * ZNS(I)
     >       + B(34) * KZ(31,I) * DEN(10,I) * ZNO(I)
      L(11,I)= KZ(11,I) * ZO(I)
     >       + KZ(36,I) * ZO2(I)
     >       + KZ(50,I) * ZNO(I)
     >       + KZ(57,I) * E(I)
     >       + KZ(58,I) * E(I)
     >       + A(5)
     >       + A(4)
      DEN(11,I) = P(11,I) / L(11,I)
C
C
C O(1D):
C
      P(12,I)= AGLW(1,1,I)
     >       + SRCED(I)
     >       + PHOTOD(1,2,I)
     >       + KZ(28,I) * E(I)  * ZO(I)
     >       + B(2)  * KZ(24,I) * DEN(6,I)  * E(I)
     >       + B(6)  * KZ(5,I)  * DEN(10,I) * ZO2(I)
     >       + B(20) * KZ(6,I)  * DEN(10,I) * ZO(I)
     >       + B(17) * KZ(25,I) * DEN(4,I)  * ZO2(I)
     >       + B(7)  * KZ(15,I) * DEN(2,I)  * ZO(I)
     >       + (1 - B(33)) * KZ(37,I) * DEN(6,I) * ZNS(I)
     >       + B(46) * KZ(36,I) * DEN(11,I) * ZO2(I)
     >       + B(47) * KZ(50,I) * DEN(11,I) * ZNO(I)
     >       + KZ(58,I) * DEN(11,I) * E(I)
     >       + A(5)  * DEN(11,I)
C    Not used in AURIC:
     >       + B(52) * KZ(33,I) * DEN(9,I) * ZO(I)
      L(12,I)= KZ(8,I)  * ZN2(I) 
     >       + KZ(9,I)  * ZO2(I)
     >       + KZ(10,I) * E(I)
     >       + KZ(27,I) * ZO(I)
     >       + KZ(51,I) * ZNO(I)
     >       + A(2)
     >       + A(3)
C    Not used in AURIC:
     >       + KZ(62,I) * E(I)
      DEN(12,I) = P(12,I) / L(12,I)
C
C
C O-:
C
      P(13,I)= KZ(59,I) * E(I) * ZO(I)
      L(13,I)= KZ(63,I) * DEN(3,I)
     >       + KZ(64,I) * ZO(I)
      DEN(13,I) = P(13,I) / L(13,I)
  300 CONTINUE
C
  350 CONTINUE
C
C
C Calculate airglow emission rates; fill ZCETA array with partial rates
C from each source; fill ZETA array with total rate for each emission:
C
      DO 400 I=1,JMAX
C
      ZCETA(1,1,I) = B(39) * AGLW(3,3,I)
      ZCETA(2,1,I) = B(40) * A(10) * P(8,I) / L(8,I)
C
      ZCETA(1,2,I) = B(38) * B(36) * RN2EI(I)
      ZCETA(2,2,I) = B(38) * PHOTOI(3,3,I)
      ZCETA(3,2,I) = G(2,I) * DEN(5,I)
C
      ZCETA(1,3,I) = A(1) * B(27) * PHOTOD(1,3,I) / L(10,I)
      ZCETA(2,3,I) = A(1) * B(27) * PHOTOI(6,3,I) / L(10,I)
      ZCETA(3,3,I) = A(1) * B(25) * RN2ED(I) / L(10,I)
      ZCETA(4,3,I) = A(1) * B(16) * B(15) * RN2EI(I) / L(10,I)
      ZCETA(5,3,I) = A(1) * B(3)  * KZ(22,I) * DEN(7,I) * E(I) /L(10,I)
      ZCETA(6,3,I) = A(1) * B(4)  * KZ(23,I) * DEN(5,I) * E(I) /L(10,I)
      ZCETA(7,3,I) = A(1) * B(5)  * KZ(3,I)  * DEN(5,I) * ZO(I) /L(10,I)
      ZCETA(8,3,I) = A(1) * B(29) * KZ(33,I) * DEN(9,I) * ZO(I) /L(10,I)
      ZCETA(9,3,I) = A(1) * A(12) * DEN(9,I) / L(10,I)
      ZCETA(10,3,I)= A(1) * KZ(45,I) * DEN(9,I) * ZNS(I)/L(10,I)
      ZCETA(11,3,I)= A(1) * (1-B(49))*KZ(54,I)*DEN(8,I) * ZNS(I)/L(10,I)
      ZCETA(12,3,I)= A(1) * KZ(56,I) * DEN(9,I) * E(I) / L(10,I)
C
      ZCETA(1,4,I) = A(5) * AGLW(2,1,I) / L(11,I)
      ZCETA(2,4,I) = A(5) * B(1)*KZ(24,I) * DEN(6,I) * E(I)  /L(11,I)
      ZCETA(3,4,I) = A(5) * B(18) * KZ(26,I) * DEN(8,I) * ZO(I) /L(11,I)
      ZCETA(4,4,I) = A(5) * B(33) * KZ(37,I) * DEN(6,I) * ZNS(I)/L(11,I)
      ZCETA(5,4,I) = A(5) * B(34) * KZ(31,I) * DEN(10,I)* ZNO(I)/L(11,I)
      ZCETA(6,4,I) = PHOTOD(2,2,I) / L(11,I)
      ZCETA(7,4,I) = A(5) * B(48) * KZ(53,I) * DEN(8,I) * ZO(I)/L(11,I)
C
      ZCETA(1,5,I) = A(2) * AGLW(1,1,I) / L(12,I)
      ZCETA(2,5,I) = A(2) * KZ(28,I) * E(I)  * ZO(I) / L(12,I)
      ZCETA(3,5,I) = A(2) * B(2)  * KZ(24,I) * DEN(6,I)  * E(I)/L(12,I)
      ZCETA(4,5,I) = A(2) * B(6)  * KZ(5,I)  * DEN(10,I) *ZO2(I)/L(12,I)
      ZCETA(5,5,I) = A(2) * B(20) * KZ(6,I)  * DEN(10,I) * ZO(I)/L(12,I)
      ZCETA(6,5,I) = A(2) * B(17) * KZ(25,I) * DEN(4,I) * ZO2(I)/L(12,I)
      ZCETA(7,5,I) = A(2) * B(7)  * KZ(15,I) * DEN(2,I)  * ZO(I)/L(12,I)
      ZCETA(8,5,I) = A(2) * SRCED(I) / L(12,I)
      ZCETA(9,5,I) = A(2) * PHOTOD(1,2,I) / L(12,I)
      ZCETA(10,5,I)= A(2) * A(5)  * DEN(11,I) / L(12,I)
      ZCETA(11,5,I)= A(2)*(1-B(33))*KZ(37,I) * DEN(6,I) * ZNS(I)/L(12,I)
      ZCETA(12,5,I)= A(2)* B(46)*KZ(36,I)*DEN(11,I)*ZO2(I)/L(12,I)
      ZCETA(13,5,I)= A(2)* B(47) * KZ(50,I) * DEN(11,I) * ZNO(I)/L(12,I)
      ZCETA(14,5,I)= A(2)* KZ(58,I) * DEN(11,I) * E(I)/L(12,I)
      ZCETA(15,5,I)= A(2)* B(52) * KZ(33,I) * DEN(9,I) * ZO(I)/L(12,I)
C
      ZCETA(1,6,I) = A(8) * (PHOTOI(3,1,I)+B(41)*PHOTOI(5,1,I)) / L(1,I)
      ZCETA(2,6,I) = A(8) * B(30) * PHOTOI(4,2,I) / L(1,I)
      ZCETA(3,6,I) = A(8) * B(9) * OEI(I) / L(1,I)
      ZCETA(4,6,I) = A(8) * B(12) * O2EI(I) / L(1,I)
C
      ZCETA(1,7,I) = A(12) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,7,I) = A(12) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,7,I) = A(12) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,7,I) = A(12) * B(23) * KZ(23,I) * DEN(5,I) * E(I) / L(9,I)
      ZCETA(5,7,I) = A(12) * B(49) * KZ(54,I) * DEN(8,I) * ZNS(I)/L(9,I)
C
      ZCETA(1,8,I) = A(11) * B(28) * PHOTOD(1,3,I) / L(9,I)
      ZCETA(2,8,I) = A(11) * B(28) * PHOTOI(6,3,I) / L(9,I)
      ZCETA(3,8,I) = A(11) * B(26) * RN2ED(I) / L(9,I)
      ZCETA(4,8,I) = A(11) * B(23) * KZ(23,I) * DEN(5,I) * E(I) / L(9,I)
      ZCETA(5,8,I) = A(11) * B(49) * KZ(54,I) * DEN(8,I) * ZNS(I)/L(9,I)
C
      ZCETA(1,9,I) = AGLW(5,1,I)
C
      ZCETA(1,10,I) = AGLW(6,1,I)
      ZCETA(2,10,I) = AGLW(7,1,I)
      ZCETA(3,10,I) = B(44) * AGLW(8,1,I)
      ZCETA(4,10,I) = .98 * KZ(60,I) * DEN(3,I) * E(I)
      ZCETA(5,10,I) = B(53) * KZ(24,I) * E(I)
      ZCETA(6,10,I) = .46 * KZ(63,I) * DEN(13,I) * ZO(I)/L(13,I)
C
      ZCETA(1,11,I) = A(6) * (PHOTOI(2,1,I)+B(42)*PHOTOI(5,1,I))/ L(2,I)
      ZCETA(2,11,I) = A(6) * B(31) * PHOTOI(4,2,I) / L(2,I)
      ZCETA(3,11,I) = A(6) * B(10) * OEI(I) / L(2,I)
      ZCETA(4,11,I) = A(6) * B(13) * O2EI(I) / L(2,I)
      ZCETA(5,11,I) = A(6) * B(8)  * KZ(18,I) * DEN(1,I) * E(I) / L(2,I)
      ZCETA(6,11,I) = A(6) * A(8)  * DEN(1,I) / L(2,I)
C
      ZETA(1,I)  = ZCETA(1,1,I)+ZCETA(2,1,I)
      ZETA(2,I)  = ZCETA(1,2,I)+ZCETA(2,2,I)+ZCETA(3,2,I)
      ZETA(3,I)  = ZCETA(1,3,I)+ZCETA(2,3,I)+ZCETA(3,3,I)
     >            +ZCETA(4,3,I)+ZCETA(5,3,I)+ZCETA(6,3,I)
     >            +ZCETA(7,3,I)+ZCETA(8,3,I)+ZCETA(9,3,I)
     >            +ZCETA(10,3,I)+ZCETA(11,3,I)
      ZETA(4,I)  = ZCETA(1,4,I)+ZCETA(2,4,I)+ZCETA(3,4,I)
     >            +ZCETA(4,4,I)+ZCETA(5,4,I)+ZCETA(6,4,I)
C     >            +ZCETA(7,4,I)
      ZETA(5,I)  = ZCETA(1,5,I)+ZCETA(2,5,I)+ZCETA(3,5,I)
     >            +ZCETA(4,5,I)+ZCETA(5,5,I)+ZCETA(6,5,I)
     >            +ZCETA(7,5,I)+ZCETA(8,5,I)+ZCETA(9,5,I)
     >            +ZCETA(10,5,I)+ZCETA(11,5,I)+ZCETA(12,5,I)
     >            +ZCETA(13,5,I)+ZCETA(14,5,I)+ZCETA(15,5,I)
      ZETA(6,I)  = ZCETA(1,6,I)+ZCETA(2,6,I)+ZCETA(3,6,I)+ZCETA(4,6,I)
      ZETA(7,I)  = ZCETA(1,7,I)+ZCETA(2,7,I)+ZCETA(3,7,I)+ZCETA(4,7,I)
     >            +ZCETA(5,7,I)
      ZETA(8,I)  = ZCETA(1,8,I)+ZCETA(2,8,I)+ZCETA(3,8,I)+ZCETA(4,8,I)
     >            +ZCETA(5,8,I)
      ZETA(9,I)  = ZCETA(1,9,I)
      ZETA(10,I) = ZCETA(1,10,I)+ZCETA(2,10,I)
C     +ZCETA(3,10,I)+ZCETA(4,10,I)+ZCETA(5,10,I)+ZCETA(6,10,I)
      ZETA(11,I) = ZCETA(1,11,I)+ZCETA(2,11,I)+ZCETA(3,11,I)
     >            +ZCETA(4,11,I)+ZCETA(5,11,I)+ZCETA(6,11,I)
C
  400 CONTINUE
C
C
C Calculate vertical column brightnesses:
C
      DO 450 IW=1,NW
      VCB(IW) = 0.
  450 CONTINUE
C
      DO 550 I=1,JMAX
      IF (I .EQ. JMAX) THEN
        DZ = (ZZ(I) - ZZ(I-1))
      ELSE
        IF (I .EQ. 1) THEN
          DZ = (ZZ(I+1) - ZZ(I))
        ELSE
          DZ = (ZZ(I+1) - ZZ(I-1)) / 2.0
        ENDIF
      ENDIF
      DO 500 IW=1,NW
      VCB(IW) = VCB(IW) + ZETA(IW,I) * DZ
  500 CONTINUE
  550 CONTINUE
C
C
C Convert brightnesses to Rayleighs:
C
      DO 600 IW=1,NW
        VCB(IW) = VCB(IW) / 1.E6
  600 CONTINUE
C
C
      RETURN
      END

