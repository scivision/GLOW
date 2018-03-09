! Subroutine GCHEM

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1988, 1989, 1992, 1999, 2005, 2016

! 3/2016: Replaced quartic solution to electron density equation
! with iterative method.  Also added constraint O+ > 0 for KCHEM=3.
! 1/2017: Removed 7774 -> 1356 fudge factor (and reduced 7774 cross section in exsect.f)

! Electron density must be supplied above 200 km in array ZE, a priori
! values are expected but not necessarily required at and below 200 km,
! depending on the value of KCHEM.  An initial guess of the N(2D)
! density in array ZND is also expected but not required.
! Other neutral species (O, N2, O2, NO, N(4S)) must be supplied,
! see subroutine GLOW.

! Chemical calculations are controlled by switch KCHEM:
! 0 = no calculations at all are performed.
! 1 = electron density, O+(4S), N+, N2+, O2+, NO+ supplied at all
!     altitudes; O+(2P), O+(2D), excited neutrals, and emission rates
!     are calculated.
! 2 = electron density, O+(4S), O2+, NO+ supplied at all altitudes;
!     O+(2P), O+(2D), N+, N2+, excited neutrals, emissions calculated.
! 3 = electron density supplied at all altitudes; everything else
!     calculated.  Note that this may violate charge neutrality and/or
!     lead to other unrealistic results below 200 km, if the electron
!     density supplied is significantly different from what the model
!     thinks it should be.  If it is desired to use a specified ionosphere,
!     KCHEM=2 is probably a better option.
! 4 = electron density supplied above 200 km; electron density below
!     200 km is calculated, everything else calculated at all altitudes.
!     Electron density for the next two levels above J200 is log interpolated
!     between E(J200) and E(J200+3).

! For definitions of use-associated variables see subroutine GLOW and module CGLOW.

! Other definitions:
! A        Einstein coefficients; s-1
! B        Branching ratios
! BZ       Altitude-dependent branching ratios
! G        Resonant scattering g-factors at each altitude; s-1
! KZ       Temperature dependent rate coeffs at each altitude; cm3s-1
! OEI      O electron impact ionization rates; cm-3s-1
! O2EI     O2   "       "        "       "      "
! RN2EI    N2   "       "        "       "      "
! TEI      Total "      "        "       "      "
! OPI      O photoionization rate, cm-3s-1
! O2PI     O2       "          "      "
! RN2PI    N2       "          "      "
! TPI      Total    "          "      "
! TIR      Total ionization rate; cm-3 s-1
! RN2ED    N2 electron impact dissociation rate; cm-3s-1
! SRCED    O2     "      "         "         " (SR continuum); cm-3s-1
! P        Volume production rate for each species, altitude; cm-3s-1
! L        Loss rate for each species, altitude; s-1
! OMINUS   O- density for mutual neutralization contribution to O*
! T1       Effective temperature divided by 300 for O+ + N2; K
! T2          "          "        "                 O+ + O2; K
! T3          "          "        "                 N2+ + O; K
! T4          "          "        "                 N2+ + O2; K
! T5          "          "        "                 O+ + NO; K
! QQ, RR, SS, TT, UU, VV, WW, XX:  Combined terms for calculation of
!                              O+(4S) given e

! Array dimensions:
! JMAX    number of altitude levels
! NBINS   number of energetic electron energy bins
! LMAX    number of wavelength intervals for solar flux
! NMAJ    number of major species
! NST     number of states produced by photoionization/dissociation
! NEI     number of states produced by electron impact
! NEX     number of ionized/excited species
! NW      number of airglow emission wavelengths
! NC      number of component production terms for each emission
! NR      number of rate coefficients, branching ratios, A and G factors

! References for rate coefficients, transition coefficients, branching ratios, and g-factors:
!      k1  O+(4S) + N2             St. Maurice & Torr, 1978 (from Albritton et al., 1977)
!      k2  O+(4S) + O2             Combination of (Chen et al, 1978 and St. Maurice & Torr, 1978)
!      k3  N2+ + O -> NO+ + O      New fit to McFarland et al, 1974
!      k4  N2+ + O2                McFarland et al, 1973
!      k5  N(2D) + O2              Lin & Kaufman, 1971, cf. Piper et al, 1987
!      k6  N(2D) + O               Fell et al., 1990
!      k7  N(2D) + e               Frederick & Rusch, 1977; Queffelec et al, 1985
!      k8  O(1D) + N2              Streit et al, 1976
!      k9  O(1D) + O2              Streit et al, 1976
!      k10 O(1D) + e               Link, 1982
!      k11 O(1S) + O               Slanger & Black, 1981
!      k12 O+(2D) + N2             Johnsen & Biondi, 1980
!      k13 O+(2D) + O2             Johnsen & Biondi, 1980
!      k14 O+(2D) + e              Henry et al., 1969
!      k15 O+(2D) + O              Torr & Torr, 1980
!      k16 O+(2P) + N2             Rusch et al, 1977
!      k17 O+(2P) + O2             Link, 1982
!      k18 O+(2P) + e              Henry et al., 1969
!      k19 O+(2P) + O              Rusch et al, 1977
!      k20 O2(c) + O               Solheim & Llewellyn, 1979
!      k21 O2(c) + N2              Solheim & Llewellyn, 1979
!      k22 NO+ + e                 Walls & Dunn, 1974; Torr et al, 1977; Alge et al, 1983;
!                                  Dulaney et al,1987; Davidson & Hobson, 1987
!      k23 N2+ + e                 Mehr & Biondi, 1969
!      k24 O2+ + e                 Mehr & Biondi; Walls & Dunn; Torr et al; Alge et al, 1983
!      k25 N+ + O2                 Langford et al, 1985
!      k26 N2(A) + O               Piper et al, 1981b (av. v=1,2)
!      k27 O(1D) + O               Abreu et al, 1986; Yee, pc, 1991
!      k28 O + et                  Link, 1982
!      k29 N2(A) + O2              Piper et al, 1981a
!      k30 O2+ + NO                Lindeger & Ferguson, 1974; G&R, 1979
!      k31 N(2D) + NO              Black et al, 1969; fr. Roble, 1986
!      k32 N+ + O                  Torr, 1985 (Levine, Photochemistry)
!      k33 N(2P) + O               Zipf et al, 1980; Young & Dunn, 1975
!      k34 N(2P) + O2              Zipf et al, 1980; Rawlins, 1988
!                                  (cf Ianuzzi & Kaufman, 1980, 3.5E-12)
!      k35 N(2P) + NO              Rees & Jones, 1973, Zipf et al, 1980
!      k36 O(1S) + O2              Slanger et al, 1972; fr. Bates, 1978
!      k37 O2+ + N                 Fehsenfeld (1977)
!      k38 O+ + N(2D)              Bates, 1989 (PSS 37, 363)
!      k39 N2+ + O -> N2 + O+      Torr, 1985; Torr et al, 1988; Knutsen et al, 1988
!      k40 O+ + NO                 St. Maurice & Torr, 1978
!      k41 O+ + e -> 7774, 1356    Melendez, 1999; Qin, 2015 (cf. Tinsley 1973; Julienne, 1974)
!      k42 O + e -> O-             Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k43 O- + O+ -> O* + O       Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k44 O- + O+ -> O2 + e       Melendez, 1999; Qin, 2015 (cf. Tinsley 1973)
!      k45 O+ + e -> 8446, 1304    Estimate from Tinsley, 1973; Julienne, 1974)
!      A1  5200       N(4S-2D)     Wiese et al, 1966
!      A2  6300       O(3P-1D)     Baluja and Zeippen, 1988
!      A3  6364       O(3P-1D)     Baluja and Zeippen, 1988
!      A4  2972       O(3P-1S)     Kernahan & Pang, 1975
!      A5  5577       O(1D-1S)     Kernahan & Pang, 1975
!      A6  3726       O+(4S-2D)    Kernahan & Pang, 1975
!      A7  2470       O+(4S-2P)    Weise et al, 1966
!      A8  7319-30    O+(2D-2P)    Weise et al, 1966
!      A9  (Hertz II) O2(X-c)      Solheim & Llewellyn, 1978
!      A10 (Veg-Kap)  N2(X-A)      Shemansky, 1969
!      A11 3466       N(4S-2P)     Chamberlain, 1961
!      A12 10400      N(2D-2P)     Chamberlain, 1961
!      B1  O(1S) from O2+ + e      Yee et al, 1988
!      B2  O(1D) from O2+ + e      Abreu et al, 1986
!      B3  N(2D) from NO+ + e      Kley et al, 1976
!      B4  N(2D) from N2+ + e      Queffelec et al, 1985
!      B5  N(2D) from N2+ + O      Frederick & Rusch, 1977
!      B6  O(1D) from N(2D) + O2   Link, 1983; Langford et al, 1985
!      B7  O(1D) from O+(2D) + O   ?
!      B8  O+(2D) from O+(2P) + e  Link, 1982
!      B9  O+(2P) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B10 O+(2D) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B11 O+(4S) from O + e*      Gerard & Rusch, 1979; Jones, 1975
!      B12 O+(2P) from O2 + e*     Link, 1982; guess of .3/3
!      B13 O+(2D) from O2 + e*     Link, 1982; guess of .3/3
!      B14 O+(4S) from O2 + e*     Link, 1982; guess of .3/3
!      B15 N+ from N2 + e*         Richards & Torr, 1985
!      B16 N(2D) from above        Zipf et al, 1980
!      B17 O(1D) from N+ + O2      Langford et al, 1985
!      B18 O(1S) from N2(A) + O    Sharp & Torr, 1979
!      B19 O(1S) from O2(*) + O    ? (= 0 at present)
!      B20 O(1D) from N(2D) + O    ?
!      B21 NO+ from N+ + O2        Langford et al, 1985
!      B22 O2+ from N+ + O2        Langford et al, 1985
!      B23 N(2P) from N2+ + e      Queffelec et al, 1985
!      B24 N2 + protons -> N + N   ?
!      B25 N(2D) from N2 + e* dis  Zipf et al, 1980
!      B26 N(2P) from N2 + e* dis  Zipf et al, 1980
!      B27 N(2D) from N2 + hv      Richards et al, 1981 (add to B28)
!      B28 N(2P) from N2 + hv      ? (cf Zipf & McGlaughlin, 1978)
!      B29 N(2D) from N(2P) + O    ?
!      B30 O+(2P) from O2 + hv     ?
!      B31 O+(2D)  "               ?
!      B32 O+(4S)  "               ?
!      B33 O(1S) from O2+ + N      Frederick et al, 1976; Kopp ea, 1977
!      B34 O(1S) from N(2D) + NO   Frederick et al, 1976; Kopp ea, 1977
!      B35 O2 + protons -> (O1D)   ?
!      B36 N2+(B) from N2 + e*     Borst & Zipf, 1970; Shemansky & Broadfoot, 1971
!      B37 (0,0) (3914) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B38 (0,1) (4278) fr. N2+(B) Shemansky & Broadfoot, 1971
!      B39 (0,0) (3371) fr. N2(C)  Conway, 1983; Benesch et al, 1966
!      B40 (0,9) (3352) fr. N2(A)  Cartwright, 1978; Shemansky, 1969
!      B41 O+(2Po) fr. O+(2Pe)     Kirby et al, 1979
!      B42 O+(2Do) fr. O+(2Pe)     Kirby et al, 1979
!      B43 N2(C) bound fraction    ?
!      B44 7990 fr. O(3s'3D)       appx. fr. Hecht, p.c.
!      B45                         not currently in use 
!      B46 N 1493 fr. N2+hv DI     guess 
!      B47 N 1493 fr. N2+e* DI     guess, cf. Mumma and Zipf (1973), Meier (1991)
!      B48 N2(a) from (a,a',w)     estimate from comparison with Ajello & Shemansky, GUVI data, etc.
!      B49 7774, 1356 fr. O-+O+    Melendez, 1999; Qin, 2015 (cf. Tinsley 1973; Julienne, 1974)
!      G1  N2+B(0,0) (3914)        Broadfoot, 1967
!      G2  N2+B(0,1) (4278)        Broadfoot, 1967


  SUBROUTINE GCHEM

    use cglow,only: jmax,  nex, nw,  kchem, sza, &
                    zz, zo, zn2, zo2, zno, zns, znd, ze, ztn, zti, zte, &
                    photoi, photod, phono, pia, sion, aglw, &
                    tei, tpi, tir, e=>ecalc, den=>zxden, zeta, zceta, vcb
!
    implicit none
    integer,parameter :: nr=50
    real,parameter :: re=6.37E8
!
    real ::   A(NR), BZ(NR,JMAX), G(NR,JMAX), KZ(NR,JMAX), &
              OEI(JMAX), O2EI(JMAX), RN2EI(JMAX), &
              OPI(JMAX), O2PI(JMAX), RN2PI(JMAX), &
              RN2ED(JMAX), SRCED(JMAX), P(NEX,JMAX), L(NEX,JMAX), OMINUS(JMAX), &
              T1(JMAX), T2(JMAX), T3(JMAX), T4(JMAX), T5(JMAX), &
              QQ(JMAX), RR(JMAX), SS(JMAX), TT(JMAX), UU(JMAX), &
              VV(JMAX), WW(JMAX), XX(JMAX)
    real ::   dz,tatomi(jmax),alphaef(jmax)
    integer :: i,j200,iter,j

    real,parameter :: B(*) = [0.07, 1.20, 0.76, 1.85, 1.00, 0.10, 0.50, 0.81, 0.20, 0.32, &
                              0.48, 0.10, 0.10, 0.10, 0.16, 0.50, 0.30, 0.19, 0.00, 0.10, &
                              0.43, 0.51, 0.10, 0.60, 0.54, 0.44, 0.80, 0.20, 1.00, 0.33, &
                              0.33, 0.34, 0.21, 0.20, 0.10, 0.11, 0.65, 0.20, 0.24, 0.02, &
                              0.18, 0.72, 0.75, 0.10, 0.00, 0.05, 0.02, 0.70, 0.54, 0.00]

    A(:) = 0.
    A(:12) = [1.07E-5, 0.00585, 0.00185, 0.0450, 1.0600, 9.7E-5, 0.0479, 0.1712, 0.0010, 0.7700, &
              0.00540, 0.07900]
              
    IF (KCHEM == 0) RETURN

! Zero airglow and density arrays:
!
    zeta(:,:) = 0.
    zceta(:,:,:) = 0.
    vcb(:) = 0.
    if (kchem >= 3) den(:,:) = 0.
    g(:,:) = 0.

! Assign g-factors at altitudes which are sunlit:

    where (SZA < 1.6 .OR. (RE+ZZ(:)) * SIN(SZA) > RE)
        G(1,:) = 0.041
        G(2,:) = 0.013
    endwhere

! Calculate rate coefficients as a function of altitude:

  T1(:) = (16.*ZTN(:)+28.*ZTI(:)) / (16.+28.) / 300.
  T2(:) = (16.*ZTN(:)+32.*ZTI(:)) / (16.+32.) / 300.
  T3(:) = (28.*ZTN(:)+16.*ZTI(:)) / (28.+16.) / 300.
  T4(:) = (28.*ZTN(:)+32.*ZTI(:)) / (28.+32.) / 300.
  T5(:) = (16.*ZTN(:)+30.*ZTI(:)) / (16.+30.) / 300.
  
  where (T1(:) < 5.6667)
    KZ(1,:) = 1.533E-12 - 5.92E-13*T1(:) + 8.6E-14*T1(:)**2
  ELSEwhere
    KZ(1,:) = 2.73E-12 - 1.155E-12*T1(:) + 1.483E-13*T1(:)**2
  ENDwhere
  
  where (T2(:) < 6.6667) 
    KZ(2,:) = 3.53E-11 - 1.84E-11*T2(:) + 4.62E-12*T2(:)**2 &
                       - 4.95E-13*T2(:)**3 + 2.00E-14*T2(:)**4
  ELSEwhere
    KZ(2,:) = 2.82E-11 - 7.74E-12*T2(:) + 1.073E-12*T2(:)**2 &
                       - 5.17E-14*T2(:)**3 + 9.65E-16*T2(:)**4
  ENDwhere
  
  KZ(3,:) = 1.793e-10 - 6.242e-11*T3(:) + 1.225e-11*T3(:)**2 &
                      - 1.016e-12*T3(:)**3 + 3.087e-14*T3(:)**4
  KZ(4,:) = 5.0E-11 * (1./T4(:)) ** 0.8
  KZ(5,:) = 6.0E-12
  KZ(6,:) = 6.9E-13
  KZ(7,:) = 5.5E-10 * (ZTE(:)/300.) ** 0.5
  KZ(8,:) = 2.0E-11 * EXP(107.8/ZTN(:))
  KZ(9,:) = 2.9E-11 * EXP(67.5 /ZTN(:))
  KZ(10,:) = 8.1E-10 * (ZTE(:)/300.) ** 0.5
  KZ(11,:) = 2.0E-14
  KZ(12,:) = 8.0E-10
  KZ(13,:) = 7.0E-10
  KZ(14,:) = 6.6E-08 * (300./ZTE(:)) ** 0.5
  KZ(15,:) = 1.0E-11
  KZ(16,:) = 4.8E-10
  KZ(17,:) = 4.8E-10
  KZ(18,:) = 1.7E-07 * (300./ZTE(:)) ** 0.5
  KZ(19,:) = 5.2E-11
  KZ(20,:) = 2.1E-11 * EXP(-1136./ZTN(:))
  KZ(21,:) = 1.0E-13
  KZ(22,:) = 4.2E-07 * (300./ZTE(:)) ** 0.85
  KZ(23,:) = 1.8E-07 * (300./ZTE(:)) ** 0.39
  KZ(24,:) = 1.95E-07 * (300./ZTE(:)) ** 0.70
  
  where (ZTE(:) >= 1200.)  KZ(24,:) = 1.6E-07 * (300./ZTE(:)) ** 0.55
  
  KZ(25,:) = 6.0E-10
  KZ(26,:) = 3.1E-11
  KZ(27,:) = 3.0E-12
  
  where (ZTE(:) < 500.) 
    KZ(28,:) = 1.0E-29
  ELSEwhere
    KZ(28,:) = 2.6E-11 * ZTE(:)**0.5 * EXP(-22740./ZTE(:))
  ENDwhere
  
  
  KZ(29,:) = 4.1E-12
  KZ(30,:) = 4.4E-10
  KZ(31,:) = 7.0E-11
  KZ(32,:) = 1.0E-12
  KZ(33,:) = 1.2E-11
  KZ(34,:) = 2.0E-12
  KZ(35,:) = 1.8E-10
  KZ(36,:) = 4.0E-12 * EXP(-865./ZTN(:))
  KZ(37,:) = 1.2E-10
  KZ(38,:) = 1.3E-10
  KZ(39,:) = 2.0E-11
  
  where (T5(:) < 5) 
    KZ(40,:) = 8.36E-13 - 2.02E-13*T5(:) + 6.95E-14*T5(:)**2
  ELSEwhere
    KZ(40,:) = 5.33E-13 - 1.64E-14*T5(:) + 4.72E-14*T5(:)**2 &
                        - 7.05E-16*T5(:)**3
  ENDwhere

  KZ(41,:) = 7.3E-13
  KZ(42,:) = 1.3E-15
  KZ(43,:) = 1.0E-7
  KZ(44,:) = 1.4E-10
  KZ(45,:) = 4.0E-13

! Calculate Electron impact ionization, photoionization, and electron
! impact dissociation rates at each altitude; put a priori electron
! density in calculated electron density array: put a priori N(2D) in DEN array:

  OEI(:)   = SION(1,:)+PIA(1,:)
  O2EI(:)  = SION(2,:)+PIA(2,:)
  RN2EI(:) = SION(3,:)+PIA(3,:)
  TEI(:)   = OEI(:)+O2EI(:)+RN2EI(:)
  OPI(:)   = PHOTOI(1,1,:)+PHOTOI(2,1,:)+PHOTOI(3,1,:)+PHOTOI(4,1,:)+PHOTOI(5,1,:)
  O2PI(:)  = PHOTOI(1,2,:)+PHOTOI(2,2,:)+PHOTOI(3,2,:)
  RN2PI(:) = PHOTOI(1,3,:)+PHOTOI(2,3,:)+PHOTOI(3,3,:)+PHOTOI(4,3,:)+PHOTOI(5,3,:)
  TPI(:)   = OPI(:)+O2PI(:)+RN2PI(:)+PHOTOI(4,2,:)+PHOTOI(6,3,:)+PHONO(1,:)
  TIR(:)   = TEI(:)+TPI(:)
  RN2ED(:) = AGLW(5,3,:)+AGLW(6,3,:)+AGLW(7,3,:)+B(24)*PIA(3,:)
  SRCED(:) = AGLW(4,2,:) + B(35)*PIA(2,:)
  E(:)     = ZE(:)
  DEN(10,:)= ZND(:)

! Find level below which electron density will be calculated:
    J200 = 0
    IF (KCHEM >= 4) THEN
      DO I=JMAX,1,-1
        IF (ZZ(I) > 2.0001E7) then
          J200=I-1   ! don't exit, let it find the last one down.
        endif
      ENDDO
    ENDIF

! Iterative loop assures that electron density and feedback reactions
! (O+(2P,2D)+e, O+(4S)+N(2D), N2++O) are correctly computed:

    itr: DO ITER=1,5

! Calculate atomic ion densities at each altitude:

! O+(2P):
!
      P(1,:)= PHOTOI(3,1,:) &
            + B(41) * PHOTOI(5,1,:) &
            + B(30) * PHOTOI(4,2,:) &
            + B(9)  * OEI(:) &
            + B(12) * O2EI(:)
      L(1,:)= KZ(16,:) * ZN2(:) &
            + KZ(17,:) * ZO2(:) &
            + KZ(19,:) * ZO(:) &
            + KZ(18,:) * E(:) &
            + A(8) &
            + A(7)
      DEN(1,:) = P(1,:) / L(1,:)

! O+(2D):

      P(2,:)= PHOTOI(2,1,:) &
            + B(42) * PHOTOI(5,1,:) &
            + B(31) * PHOTOI(4,2,:) &
            + B(10) * OEI(:) &
            + B(13) * O2EI(:) &
            + B(8)  * KZ(18,:) * DEN(1,:) * E(:) &
            + A(8)  * DEN(1,:)
      L(2,:)= KZ(12,:) * ZN2(:) &
            + KZ(13,:) * ZO2(:) &
            + KZ(15,:) * ZO(:) &
            + KZ(14,:) * E(:) &
            + A(6)
      DEN(2,:) = P(2,:) / L(2,:)
! N+:

      IF (KCHEM >= 2) THEN
        P(4,:) = PHOTOI(6,3,:) &
               + B(15) * RN2EI(:) &
               + KZ(38,:) * DEN(3,:) * DEN(10,:)
        L(4,:) = KZ(25,:) * ZO2(:) &
               + KZ(32,:) * ZO(:)
        DEN(4,:) = P(4,:) / L(4,:)
      ENDIF

! O+(4S):
!
      IF (KCHEM >= 3) THEN
        P(3,:)= PHOTOI(1,1,:) + PHOTOI(4,1,:) &
              + B(32) * PHOTOI(4,2,:) &
              + B(11) * OEI(:) &
              + B(14) * O2EI(:)  &
              + KZ(14,:) * DEN(2,:) * E(:) & 
              + KZ(15,:) * DEN(2,:) * ZO(:) & 
              + A(6) * DEN(2,:) &
              + (1.-B(8)) * KZ(18,:) * DEN(1,:) * E(:) &
              + KZ(19,:) * DEN(1,:) * ZO(:) & 
              + A(7) * DEN(1,:) &
              + KZ(32,:) * DEN(4,:) * ZO(:) &
              + KZ(39,:) * DEN(5,:) * ZO(:)
        L(3,:)= KZ(1,:) * ZN2(:) &
              + KZ(2,:) * ZO2(:) &
              + KZ(38,:) * DEN(10,:)
        DEN(3,:) = P(3,:) / L(3,:)
      ENDIF

! Above 200 km, (or at all altitudes if KCHEM=3) use a priori
! electron density to calculate O+(4S):
!
    IF (KCHEM >= 3) THEN
      j = j200+1
      
      P(5,J:JMAX) = RN2PI(J:JMAX) &
        + (1.-B(15)) * RN2EI(J:JMAX) &
        + KZ(12,J:JMAX) * DEN(2,J:JMAX) * ZN2(J:JMAX) &
        + KZ(16,J:JMAX) * DEN(1,J:JMAX) * ZN2(J:JMAX)
        
      L(5,J:JMAX)= KZ(3,J:JMAX)  * ZO(J:JMAX) &
              + KZ(4,J:JMAX)  * ZO2(J:JMAX) &
              + KZ(23,J:JMAX) * E(J:JMAX) &
              + KZ(39,J:JMAX) * ZO(J:JMAX)
      DEN(5,J:JMAX) = P(5,J:JMAX) / L(5,J:JMAX)
      QQ(J:JMAX) = PHONO(1,J:JMAX) &
              + KZ(3,J:JMAX)  * DEN(5,J:JMAX) * ZO(J:JMAX) &
              + B(21) * KZ(25,J:JMAX) * DEN(4,J:JMAX) * ZO2(J:JMAX)
      RR(J:JMAX) = KZ(30,J:JMAX) * ZNO(J:JMAX) &
              + KZ(37,J:JMAX) * ZNS(J:JMAX)
      SS(J:JMAX) = KZ(1,J:JMAX) * ZN2(J:JMAX) &
              + KZ(40,J:JMAX) * ZNO(J:JMAX)
      TT(J:JMAX) = KZ(22,J:JMAX) * E(J:JMAX)
      UU(J:JMAX) = O2PI(J:JMAX) &
              + (1.-B(12)-B(13)-B(14)) * O2EI(J:JMAX) &
              + KZ(13,J:JMAX) * DEN(2,J:JMAX) * ZO2(J:JMAX) &
              + KZ(17,J:JMAX) * DEN(1,J:JMAX) * ZO2(J:JMAX) &
              + KZ(4,J:JMAX)  * DEN(5,J:JMAX) * ZO2(J:JMAX) &
              + B(22) * KZ(25,J:JMAX) * DEN(4,J:JMAX) * ZO2(J:JMAX)
      VV(J:JMAX) = KZ(2,J:JMAX) * ZO2(J:JMAX)
      WW(J:JMAX) = KZ(24,J:JMAX) * E(J:JMAX) &
              + KZ(30,J:JMAX) * ZNO(J:JMAX) &
              + KZ(37,J:JMAX) * ZNS(J:JMAX)
      XX(J:JMAX) = DEN(1,J:JMAX) + DEN(2,J:JMAX) + DEN(4,J:JMAX) + DEN(5,J:JMAX)
      DEN(3,J:JMAX) = (TT(J:JMAX)*WW(J:JMAX)*E(J:JMAX) - TT(J:JMAX)*WW(J:JMAX)*XX(J:JMAX) - TT(J:JMAX)*UU(J:JMAX) &
                 - QQ(J:JMAX)*WW(J:JMAX) - RR(J:JMAX)*UU(J:JMAX) ) / &
                    (TT(J:JMAX)*WW(J:JMAX) + TT(J:JMAX)*VV(J:JMAX) + RR(J:JMAX)*VV(J:JMAX) + SS(J:JMAX)*WW(J:JMAX))

    
      where (den(3,J:JMAX) < 0.)   den(3,J:JMAX)=0.

    ENDIF

! If KCHEM=4, calculate electron density below 200 km using iterative method:
! First time: approximate electron density using effective recombination rate.
! Subsequent iterations: update electron density from sum of ions.
   
   if (kchem >= 4) then
      if (iter == 1) then
          tatomi(:j200) = sum(den(1:4,:j200),1)
          alphaef(:j200) = (kz(22,:j200) + kz(24,:j200))/2.
          e(:j200)=(tatomi(:j200)+sqrt(tatomi(:j200)**2 + 4*tir(:j200)/alphaef(:j200)))/2.
      else
         e(:j200) = (sum(den(:7,:j200),1) + e(:j200)) / 2.
      endif
!
!
! Smoothly transition to electron density above 200 km:
!
      E(J200+1) = E(J200) * ( E(J200+3) / E(J200) ) &
                  ** ( (ZZ(J200+1)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
      E(J200+2) = E(J200) * (E(J200+3)/E(J200)) &
                  ** ( (ZZ(J200+2)-ZZ(J200)) / (ZZ(J200+3)-ZZ(J200)) )
!
    endif
!
!
! Calculate molecular ion densities and excited species densites:

! N2+:
!
  IF (KCHEM >= 2) THEN
    P(5,:)= RN2PI(:) &
          + (1.-B(15)) * RN2EI(:) &
          + KZ(12,:) * DEN(2,:) * ZN2(:) &
          + KZ(16,:) * DEN(1,:) * ZN2(:)
    L(5,:)= KZ(3,:)  * ZO(:) &
          + KZ(4,:)  * ZO2(:) &
          + KZ(23,:) * E(:) &
          + KZ(39,:) * ZO(:)
      DEN(5,:) = P(5,:) / L(5,:)
  ENDIF
  
! O2+:
!
  IF (KCHEM >= 3) THEN
    P(6,:)= O2PI(:) &
          + (1.-B(12)-B(13)-B(14)) * O2EI(:) &
          + KZ(2,:)  * DEN(3,:) * ZO2(:) &
          + KZ(13,:) * DEN(2,:) * ZO2(:) &
          + KZ(17,:) * DEN(1,:) * ZO2(:) &
          + KZ(4,:)  * DEN(5,:) * ZO2(:) &
          + B(22) * KZ(25,:) * DEN(4,:) * ZO2(:)
    L(6,:)= KZ(24,:) * E(:) &
          + KZ(30,:) * ZNO(:) &
          + KZ(37,:) * ZNS(:)
    DEN(6,:) = P(6,:)/ L(6,:)
  ENDIF
  
  
! NO+:
!
  IF (KCHEM >= 3) THEN
    P(7,:)= PHONO(1,:) &
          + KZ(1,:)  * DEN(3,:) * ZN2(:) &
          + KZ(40,:) * DEN(3,:) * ZNO(:) &
          + KZ(3,:)  * DEN(5,:) * ZO(:) &
          + B(21) * KZ(25,:) * DEN(4,:) * ZO2(:) &
          + KZ(30,:) * DEN(6,:) * ZNO(:) &
          + KZ(37,:) * DEN(6,:) * ZNS(:)
    L(7,:)= KZ(22,:) * E(:)
    DEN(7,:) = P(7,:) / L(7,:)
  ENDIF
  
! N2(A):
!
  P(8,:)= AGLW(1,3,:) + AGLW(2,3,:) + B(43)*AGLW(3,3,:)
  L(8,:)= KZ(26,:) * ZO(:) &
        + KZ(29,:) * ZO2(:) &
        + A(10)
  DEN(8,:) = P(8,:) / L(8,:)
  
! N(2P):
!
  P(9,:)= B(28) * PHOTOD(1,3,:) &
        + B(28) * PHOTOI(6,3,:) &
        + B(26) * RN2ED(:) &
        + B(23) * KZ(23,:) * DEN(5,:) * E(:)
  L(9,:)= KZ(33,:) * ZO(:) &
        + KZ(34,:) * ZO2(:) &
        + KZ(35,:) * ZNO(:) &
        + A(11) &
        + A(12)
  DEN(9,:) = P(9,:) / L(9,:)
  
! N(2D):
!
  P(10,:)= B(27) * PHOTOD(1,3,:) &
         + B(27) * PHOTOI(6,3,:) &
         + B(25) * RN2ED(:) &
         + B(16) * B(15) * RN2EI(:) &
         + B(3)  * KZ(22,:) * DEN(7,:) * E(:) &
         + B(4)  * KZ(23,:) * DEN(5,:) * E(:) &
         + B(5)  * KZ(3,:)  * DEN(5,:) * ZO(:) &
         + B(29) * KZ(33,:) * DEN(9,:) * ZO(:) &
         + A(12) * DEN(9,:)
  L(10,:)= KZ(5,:)  * ZO2(:) &
         + KZ(6,:)  * ZO(:) &
         + KZ(7,:)  * E(:) &
         + KZ(31,:) * ZNO(:) &
         + KZ(38,:) * DEN(3,:) &
         + A(1)
  DEN(10,:) = P(10,:) / L(10,:)

! O(1S):
!
  BZ(1,:) = 0.12 + 0.02 * LOG10 (E(:)/ZO(:)*(300./ZTE(:))**0.7)
  
  where (BZ(1,:) < 0.03)  BZ(1,:)=0.03
  
  P(11,:)= AGLW(2,1,:) &
         + BZ(1,:) * KZ(24,:) * DEN(6,:)  * E(:) &
         + B(18) * KZ(26,:) * DEN(8,:) * ZO(:) &
         + B(33) * KZ(37,:) * DEN(6,:) * ZNS(:) &
         + B(34) * KZ(31,:) * DEN(10,:) * ZNO(:) &
         + PHOTOD(2,2,:)
  L(11,:)= KZ(11,:) * ZO(:) &
         + KZ(36,:) * ZO2(:) &
         + A(5) &
         + A(4)
  DEN(11,:) = P(11,:) / L(11,:)
  

! O(1D):
!
  P(12,:)= AGLW(1,1,:) &
         + KZ(28,:) * E(:)  * ZO(:) &
         + B(2)  * KZ(24,:) * DEN(6,:)  * E(:) &
         + B(6)  * KZ(5,:)  * DEN(10,:) * ZO2(:) &
         + B(20) * KZ(6,:)  * DEN(10,:) * ZO(:) &
         + B(17) * KZ(25,:) * DEN(4,:)  * ZO2(:) &
         + B(7)  * KZ(15,:) * DEN(2,:)  * ZO(:) &
         + SRCED(:) &
         + PHOTOD(1,2,:) &
         + A(5)  * DEN(11,:)
  L(12,:)= KZ(8,:)  * ZN2(:)  &
         + KZ(9,:)  * ZO2(:) &
         + KZ(10,:) * E(:) &
         + KZ(27,:) * ZO(:) &
         + A(2) &
         + A(3)
  DEN(12,:) = P(12,:) / L(12,:)

ENDDO itr   ! bottom of iterative looop


! Impose charge neutrality:

    e(:j200) = sum(den(:7,:j200),1)


! Calculate O- for mutual neutralization source of O*

    ominus(:) = (kz(42,:)*zo(:)*e(:)) / (kz(43,:)*den(3,:)+kz(44,:)*zo(:))

! Calculate airglow emission rates; fill ZCETA array with partial rates
! from each source; fill ZETA array with total rate for each emission:

  ZCETA(1,1,:) = B(39) * AGLW(3,3,:)
  ZCETA(2,1,:) = B(40) * A(10) * P(8,:) / L(8,:)

  ZCETA(1,2,:) = B(38) * B(36) * RN2EI(:)
  ZCETA(2,2,:) = B(38) * PHOTOI(3,3,:)
  ZCETA(3,2,:) = G(2,:) * DEN(5,:)

  ZCETA(1,3,:) = A(1) * B(27) * PHOTOD(1,3,:) / L(10,:)
  ZCETA(2,3,:) = A(1) * B(27) * PHOTOI(6,3,:) / L(10,:)
  ZCETA(3,3,:) = A(1) * B(25) * RN2ED(:) / L(10,:)
  ZCETA(4,3,:) = A(1) * B(16) * B(15) * RN2EI(:) / L(10,:)
  ZCETA(5,3,:) = A(1) * B(3)  * KZ(22,:) * DEN(7,:) * E(:) /L(10,:)
  ZCETA(6,3,:) = A(1) * B(4)  * KZ(23,:) * DEN(5,:) * E(:) /L(10,:)
  ZCETA(7,3,:) = A(1) * B(5)  * KZ(3,:)  * DEN(5,:) * ZO(:) /L(10,:)
  ZCETA(8,3,:) = A(1) * B(29) * KZ(33,:) * DEN(9,:) * ZO(:) /L(10,:)
  ZCETA(9,3,:) = A(1) * A(12) * DEN(9,:) / L(10,:)

  ZCETA(1,4,:) = A(5) * AGLW(2,1,:) / L(11,:)
  ZCETA(2,4,:) = A(5) * BZ(1,:)*KZ(24,:) * DEN(6,:) * E(:)  /L(11,:)
  ZCETA(3,4,:) = A(5) * B(18) * KZ(26,:) * DEN(8,:) * ZO(:) /L(11,:)
  ZCETA(4,4,:) = A(5) * B(33) * KZ(37,:) * DEN(6,:) * ZNS(:)/L(11,:)
  ZCETA(5,4,:) = A(5) * B(34) * KZ(31,:) * DEN(10,:)* ZNO(:)/L(11,:)
  ZCETA(6,4,:) = PHOTOD(2,2,:) / L(11,:)

  ZCETA(1,5,:) = A(2) * AGLW(1,1,:) / L(12,:)
  ZCETA(2,5,:) = A(2) * KZ(28,:) * E(:)  * ZO(:) / L(12,:)
  ZCETA(3,5,:) = A(2) * B(2)  * KZ(24,:) * DEN(6,:)  * E(:)/L(12,:)
  ZCETA(4,5,:) = A(2) * B(6)  * KZ(5,:)  * DEN(10,:) *ZO2(:)/L(12,:)
  ZCETA(5,5,:) = A(2) * B(20) * KZ(6,:)  * DEN(10,:) * ZO(:)/L(12,:)
  ZCETA(6,5,:) = A(2) * B(17) * KZ(25,:) * DEN(4,:) * ZO2(:)/L(12,:)
  ZCETA(7,5,:) = A(2) * B(7)  * KZ(15,:) * DEN(2,:)  * ZO(:)/L(12,:)
  ZCETA(8,5,:) = A(2) * SRCED(:) / L(12,:)
  ZCETA(9,5,:) = A(2) * PHOTOD(1,2,:) / L(12,:)
  ZCETA(10,5,:)= A(2) * A(5)  * DEN(11,:) / L(12,:)

  ZCETA(1,6,:) = A(8) * (PHOTOI(3,1,:)+B(41)*PHOTOI(5,1,:)) / L(1,:)
  ZCETA(2,6,:) = A(8) * B(30) * PHOTOI(4,2,:) / L(1,:)
  ZCETA(3,6,:) = A(8) * B(9) * OEI(:) / L(1,:)
  ZCETA(4,6,:) = A(8) * B(12) * O2EI(:) / L(1,:)

  ZCETA(1,7,:) = A(12) * B(28) * PHOTOD(1,3,:) / L(9,:)
  ZCETA(2,7,:) = A(12) * B(28) * PHOTOI(6,3,:) / L(9,:)
  ZCETA(3,7,:) = A(12) * B(26) * RN2ED(:) / L(9,:)
  ZCETA(4,7,:) = A(12) * B(23) * KZ(23,:) * DEN(5,:) * E(:) / L(9,:)

  ZCETA(1,8,:) = A(11) * B(28) * PHOTOD(1,3,:) / L(9,:)
  ZCETA(2,8,:) = A(11) * B(28) * PHOTOI(6,3,:) / L(9,:)
  ZCETA(3,8,:) = A(11) * B(26) * RN2ED(:) / L(9,:)
  ZCETA(4,8,:) = A(11) * B(23) * KZ(23,:) * DEN(5,:) * E(:) / L(9,:)

  ZCETA(1,9,:) = AGLW(5,1,:)
  ZCETA(2,9,:) = KZ(41,:) * DEN(3,:) * E(:)
  ZCETA(3,9,:) = B(49) * KZ(43,:) * OMINUS(:) * DEN(3,:)

  ZCETA(1,10,:) = AGLW(6,1,:)
  ZCETA(2,10,:) = AGLW(7,1,:)
  ZCETA(3,10,:) = B(44) * AGLW(8,1,:)
   
  ZCETA(1,11,:) = A(6) * (PHOTOI(2,1,:)+B(42)*PHOTOI(5,1,:))/ L(2,:)
  ZCETA(2,11,:) = A(6) * B(31) * PHOTOI(4,2,:) / L(2,:)
  ZCETA(3,11,:) = A(6) * B(10) * OEI(:) / L(2,:)
  ZCETA(4,11,:) = A(6) * B(13) * O2EI(:) / L(2,:)
  ZCETA(5,11,:) = A(6) * B(8)  * KZ(18,:) * DEN(1,:) * E(:) / L(2,:)
  ZCETA(6,11,:) = A(6) * A(8)  * DEN(1,:) / L(2,:)

  ZCETA(1,12,:) = AGLW(4,3,:) * B(48)

  ZCETA(1,13,:) = AGLW(3,1,:)
  ZCETA(2,13,:) = AGLW(5,1,:)
  ZCETA(3,13,:) = KZ(41,:) * DEN(3,:) * E(:)
  ZCETA(4,13,:) = B(49) * KZ(43,:) * OMINUS(:) * DEN(3,:)

  ZCETA(1,14,:) = B(46)*PHOTOI(6,3,:)
  ZCETA(2,14,:) = B(47)*B(15)*RN2EI(:)
   
  ZCETA(1,15,:) = AGLW(4,1,:)
  ZCETA(2,15,:) = AGLW(6,1,:)
  ZCETA(3,15,:) = AGLW(7,1,:)
  ZCETA(4,15,:) = B(44) * AGLW(8,1,:)
  ZCETA(5,15,:) = KZ(45,:) * DEN(3,:) * E(:)
   
   

  ZETA(1,:) = ZCETA(1,1,:)+ZCETA(2,1,:)
  ZETA(2,:) = ZCETA(1,2,:)+ZCETA(2,2,:)+ZCETA(3,2,:)
  ZETA(3,:) = ZCETA(1,3,:)+ZCETA(2,3,:)+ZCETA(3,3,:)+ZCETA(4,3,:)+ZCETA(5,3,:)+ZCETA(6,3,:)+ZCETA(7,3,:)+ZCETA(8,3,:)+ZCETA(9,3,:)
  ZETA(4,:) = ZCETA(1,4,:)+ZCETA(2,4,:)+ZCETA(3,4,:)+ZCETA(4,4,:)+ZCETA(5,4,:)+ZCETA(6,4,:)
  ZETA(5,:)  = ZCETA(1,5,:)+ZCETA(2,5,:)+ZCETA(3,5,:)+ZCETA(4,5,:)+ZCETA(5,5,:)+ZCETA(6,5,:) &
              +ZCETA(7,5,:)+ZCETA(8,5,:)+ZCETA(9,5,:)+ZCETA(10,5,:)
              
  ZETA(6:8,:) = sum(ZCETA(1:4,6:8,:),1)
 
  ZETA(9:10,:) = sum(ZCETA(1:3,9:10,:),1)
  
  ZETA(11,:) = ZCETA(1,11,:)+ZCETA(2,11,:)+ZCETA(3,11,:)+ZCETA(4,11,:)+ZCETA(5,11,:)+ZCETA(6,11,:)
  ZETA(12,:) = ZCETA(1,12,:)
  ZETA(13,:) = ZCETA(1,13,:)+ZCETA(2,13,:)+ZCETA(3,13,:)+ZCETA(4,13,:)
  ZETA(14,:) = ZCETA(1,14,:)+ZCETA(2,14,:)
  ZETA(15,:) = ZCETA(1,15,:)+ZCETA(2,15,:)+ZCETA(3,15,:)+ZCETA(4,15,:)+ZCETA(5,15,:)


! Calculate vertical column brightnesses:
!
    DO I=1,JMAX
      IF (I .EQ. JMAX) THEN
        DZ = (ZZ(I) - ZZ(I-1))
      ELSE
        IF (I .EQ. 1) THEN
          DZ = (ZZ(I+1) - ZZ(I))
        ELSE
          DZ = (ZZ(I+1) - ZZ(I-1)) / 2.0
        ENDIF
      ENDIF
  
    VCB(:) = VCB(:) + ZETA(:,I) * DZ
  
    ENDDO

! Convert brightnesses to Rayleighs:
  VCB(:) = VCB(:) / 1.E6
 

  END SUBROUTINE GCHEM

