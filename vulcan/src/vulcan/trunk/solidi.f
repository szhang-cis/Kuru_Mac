      SUBROUTINE SOLIDI(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .                  DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .                  KUNL4,RETEN,RETEP,TEMPG,CCERO,DESIG,
     .                  CCEROM,CCEROP,BULKM,DISTM,
     .                  DSTRA,STRAN,DTEMG,
     .                  A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                  STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,FPCHL,
     .                  SHAPE,DBASE,SHRIN,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  GPCOD,VANIS,STRS0,TRMTX,EBASE)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES & THE CONSTITUTIVE TENSOR
C
C***********************************************************************
C
C     FOR UNCOUPLED PROBLEMS:
C
C     SET INCREMENT OF TEMPERATURE EQUAL TO ZERO (zerote.f)
C
C
C     FOR COUPLED (THERMOMECHANICAL OR SOLIDIFICATION) PROBLEMS:
C
C     EVALUATES TOTAL TEMPERATURE-DEPENDENT STRESSES & CONSTITUTIVE 
C     TENSOR
C
C     Notes:
C
C     The dimensions of the elastic constitutive tensor DMATX are:
C     -NSTRS*NSTRS in the integration of the constiutive model
C     -NSTRE*NSTRE in the solution of the equilibrium equation
C     The same idea is applied to the elasto-plastic/viscoplastic
C     constitutive tensor DMTEP
C
C
C     Free energy models:
C
C     SMALL DISPLACEMENTS & STRAINS (default model=3)
C
C     IFREN=1  psi= 1/2 rho_o (e-e^p):C(T):(e-e^p)-
C                   -1\rho_o beta(T,T_ref):(e-e^p)(T-T_ref)+
C                   +beta(T_o,T_ref):(e-e^p)(T_o-T_ref)+...     (thesis)
C     IFREN=2  psi= 1/2 rho_o (e-e^p-e^th):C(T):(e-e^p-e^th)+... (Ortiz)
C     IFREN=3  psi= 1/2 \rho_o (e-e^p):C:(e-e^p)-
C                   -1\rho_o C(T):alpha_th^s(T,T_ref):(e-e^p)(T-T_ref)+
C                   +1\rho_o C(T):alpha_th^s(T_o,T_ref):(e-e^p)
C                   (T_o-T_ref)+...                        (imp. thesis)
C
C     LARGE DISPLACEMENTS & STRAINS
C
C     IFREN=1  psi=
C                                     hyperelastic linear material model
C
C     IFREN=2  psi=
C                                     hyperelastic linear spatial model
C                                     tau = C : e
C
C     IFREN=3  psi=
C                                     elastic model
C                                     sigma = C : e
C
C     IFREN=4  psi=
C                       idem IFREN=1 with non-linear thermal deformation
C
C     IFREN=5  psi=
C                       idem IFREN=2 with non-linear thermal deformation
C
C     IFREN=6  psi=
C                       idem IFREN=3 with non-linear thermal deformation
C
C     IFREN=7  psi=
C                       Simo model
C
C     IFREN=8  psi=
C                       identical to IFREN=2 computing in a different
C                       but absolutely equivalent form the stress and
C                       constitutive tensors
C
C     IFREN=9  psi=
C                       hyperelastic linear material-spatial model
C
C     IFREN=10 psi=
C                       idem IFREN=9 with non-linear thermal deformation
C
C     IFREN=...
C
C     IFREN=51          Mooney-Rivlin model
C              psi=C1*(I1_C - 3) + C2*(I2_C - 3) + 
C                  C3*(I1_C - 3)*(I2_C - 3) + Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       Clambda=penalty factor to achieve
C                               incompressibility
C                       f(J)=volumetric contribution
C               Note: the Neo-hookean model is obtained for C1=C2=0
C
C     IFREN=52          Yeoh model
C              psi=C1*(I1_C - 3) + C2*(I1_C - 3)^2 + C3*(I1_C - 3)^3 +
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       Clambda=penalty factor to achieve incompressib.
C               Note: the Neo-hookean model is obtained for C1=C2=0
C
C     IFREN=53          Ogden's model
C              psi=C1/E1*(L1^E1+L2^E1+L3^E1 - 3) +
C                  C2/E2*(L1^E2+L2^E2+L3^E2 - 3) +
C                  C3/E3*(L1^E3+L2^E3+L3^E3 - 3) +
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       E1=first exponent
C                       E2=second exponent
C                       E3=third exponent
C                       L1, L2 & L3=eigenvalues of Right Cauchy Tensor C
C                       Clambda=penalty factor to achieve incompressib.
C               Note: the two-constant Mooney-Rivlin model is obtained
C                     for C1=2*C1_MR, E1=2, C2=-2*C2_MR & E2=-2.
C
C     IFREN=54          Delfino's model
C              psi=C1/C2*(exp(C2/2*(I1_C - 3)) - 1)+
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       Clambda=penalty factor to achieve incompressib.
C
C     IFREN=55          Holzapfel's model
C              psi=C1/2*(I1_C - 3) +
C                  C2/(2*C3)*(exp(C3*(I4_C - 1)*(I4_C - 1)) - 1)+
C                  C2/(2*C3)*(exp(C3*(I6_C - 1)*(I6_C - 1)) - 1)+
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       Clambda=penalty factor to achieve incompressib.
C                       I4_C=C:A1 related to 1st anisotropic direction
C                       I6_C=C:A2 related to 2nd anisotropic direction
C
C     IFREN=56          Myocardium model (by Holzapfel & Ogden, 2009)
C              psi=a/(2b)*exp(b(I1_C - 3)) +
C                  af/(2*bf)*(exp(bf*(I4f_C - 1)*(I4f_C - 1)) - 1)+
C                  as/(2*bs)*(exp(bs*(I4s_C - 1)*(I4s_C - 1)) - 1)+
C                  afs/(2*bfs)*(exp(bfs*I8fs*I8fs) - 1)+
C                  Clambda*f(J)
C                       a,b,af,bf,as,bs,afs,bfs=material parameters (8)
C                       Clambda=penalty factor to achieve incompressib.
C                       I4f_C=f0.C.f0 related to fiber direction
C                       I4s_C=s0.C.s0 related to surface direction
C                       I8fs_C=f0.C.s0 related to both directions
C
C     IFREN=57          Gasser's model
C              psi=C1/2*(I1_C - 3) +
C                  C2/(2*C3)*(exp(C3*(I4_C - 1)*(I4_C - 1)) - 1)+
C                  C2/(2*C3)*(exp(C3*(I6_C - 1)*(I6_C - 1)) - 1)+
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       C4=fourth constant
C                       Clambda=penalty factor to achieve incompressib.
C                       I4_C=C:A1 related to 1st anisotropic direction
C                                 & C4
C                       I6_C=C:A2 related to 2nd anisotropic direction
C                                 & C4
c
C     IFREN=58          KAN model
C              psi=C1*(1/C2*exp(C2*(I1_C - 3)) + 
C                  C3*(I1_C - 2)*(1-ln(I1_C - 2)) - 1/C2 - C3)+
C                  Clambda*f(J)
C                       C1=first constant
C                       C2=second constant
C                       C3=third constant
C                       Clambda=penalty factor to achieve incompressib.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      COMMON/PROPSC/POISC1,POISC2,POISC3,STRRA,STRR4,DAMAP1,DAMAP2
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION DMATX(NSTRS,*), PRESG(*),
     .          PROPS(*),       SGTOT(*),
     .          SIGMA(*),       DMTEP(*),
     .          DESIG(*),       DSTRA(*),
     .          STRAN(*),       STRAP(*)
      DIMENSION SIGM0(6),       DMATE(6,6),
     .          DMATP(6,6),     UNOMA(6)
      DIMENSION DMAT1(6,6),     DMAT2(6,6)
      DIMENSION DMAPL(6),       DEETH(6)
      DIMENSION STRAT(6),       FPCHL(NFPCH,*),
     .          SHAPE(*),       DBASE(*),       VANIS(NANIV,*)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*), GPCOD(*)
      DIMENSION RCGTT(6),       RCGTI(6),       ENEYE(6), ALAMB(3)
      DIMENSION DCOFA(6,6),     UNORA(6,6),     STILD(6)
      DIMENSION STRS0(*),       EBASE(*)
      DIMENSION RCGPI(6),       FINET(6),       DCINV(6,6)
      DIMENSION ROMTX(3,3),     TRMTX(NDIME,*), TRMTT(3,3)
      DIMENSION CPOL9(MCHAI),   CPOL10(MCHAI)
C
C**** INITIALISES SOME PARAMETERS
C
      KUNL4=0
C
      DO 20 ISTR1=1,NSTR1
      STRAT(ISTR1)=0.0D0
   20 UNOMA(ISTR1)=0.0D0
      UNOMA(1)=1.0D0
      UNOMA(2)=1.0D0
      UNOMA(4)=1.0D0
C
      CALL DARECO(DBASE,DAMAG,TAUAM,TAUAP)
      DAMAX=DAMAG
C
C**** MATERIAL PROPERTIES (EXPLICIT TEMPERATURE FUNCTIONS)
C
      CALL YOUNGT(TEMPG,PROPS,
     .            YOUNG,YOUNGF,YOUNGA)
      CALL POISMT(TEMPG,PROPS,
     .            POISM,POISMF,POISMA)
      CALL SHEART(TEMPG,PROPS,
     .            SHEAR,SHEARF,SHEARA)
      CALL CPOLYT(TEMPG,PROPS,
     .            CPOL1,CPOL2,CPOL3,EPOL1,EPOL2,EPOL3,
     .            CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .            CPOL9,CPOL10,
     .            CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .            CEMIN,CEMAX,CEETA)
#ifndef restricted
      IF(ITERME.GE.0) THEN                  ! coupled problems
       CALL ALPHAT(TEINI,PROPS,
     .             ALPHI,ALPHIF,ALPHIA)
       TEMPO=TEMPG-DTEMG
       CALL ALPHAT(TEMPO,PROPS,
     .             ALPHO,ALPHOF,ALPHOA)
       CALL ALPHAT(TEMPG,PROPS,
     .             ALPHA,ALPHAF,ALPHAA)
      ELSE                                  ! uncoupled problems
#endif
       TEREF=0.0D0
       TEINI=0.0D0
       ALPHI=0.0D0
       ALPHO=0.0D0
       ALPHA=0.0D0
#ifndef restricted
      ENDIF
#endif
      CALL CCEROT(TEMPG,PROPS,DTEMG,
     .            CCERO,DCCER,CCEROF,CCEROA)
C
C**** SOLIDIFICATION VARIABLE (EXPLICIT TEMPERATURE FUNCTION)
C
      ROTET=1.0D0
      DEROT=0.0D0
      LIQUI=0
      DELLS=0.0D0
C
      ROTES=1.0D0
      DEROS=0.0D0
      LIQUS=0
      DELSS=0.0D0
C
#ifndef restricted
      IF(ITERME.GE.0)                       ! coupled problems
     . CALL ROTETT(ROTET,DEROT,LIQUI,DELLS,
     .             TEMPG,DTEMG,SHAPE,FPCHL,
     .             ROTES,DEROS,LIQUS,DELSS,
     .             STRRB,STRRC,
     .             YOUNG,CCERO,CCEROM,CCEROP)
C
C**** COMPUTES DEFORMATION GRADIENT & THE RIGHT CAUCHY-GREEN TENSOR
C
      IF(LARGE.NE.0) CALL CALCFJ(XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                           SHAPE,GPCOD,RCGTT,RCGTI,TRACC,
     .                           TRACI,ENEYE,DETJC,SEINC,FACTJ,
     .                           DCOFA,UNORA,DCINV,ROMTX,UNOMA)
C
C**** PHASE-CHANGE DEFORMATION
C
      IF(ITERME.GE.0)                       ! coupled problems
     . CALL CALCPD(UNOMA,STRAT,DEETH,STRRB,STRRC,STRRE,
     .             RCGTT)
C
C**** COMPUTES THE THERMAL DEFORMATION
C
      STRRA=0.0D0
      IF(ITERME.GE.0)                       ! coupled problems
     . CALL CALCTD(PROPS,UNOMA,STRAT,DEETH,STRRA,STRRD,
     .             TEINI,TEMPO,TEMPG,DTEMG,
     .             ALPHI,ALPHO,ALPHA,
     .             STRAN,STRAP,RCGTT,YOUNG,POISM)
C
C**** COMPUTES LIQUID & SHRINKAGE DEFORMATIONS
C
      IF(LIQUI.EQ.1) CALL CALCLD(UNOMA,STRAT,STRAN,STRAP,RCGTT,
     .                           XJACM,XJACI,XJA3M,XJA3I,
     .                           SHRIN)
#endif
C
C**** SECANT CONSTITUTIVE TENSOR TO INTEGRATE THE STRESSES
C
      CALL CALCCT(DMATX,YOUNG,POISM,ROTET,NSTRS,
     .            POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .            BULKM,DISTM,DAMAP1,DAMAP2,PROPS,
     .            VANIS,ROMTX,TRMTX,
     .            YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
#ifndef restricted
      IF(ISOTR.NE.0) THEN
       DO IDIME=1,NDIME     ! transpose of TRMTX
        DO JDIME=1,NDIME
         TRMTT(IDIME,JDIME)=TRMTX(JDIME,IDIME)
        ENDDO
       ENDDO
       CALL FOURTH(DMATX,XJACM,TRMTT,XJA3M,1.0D0,1.0D0,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,NSTRS,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAP,RCGPI,FINET,TRABE,DCINV,EBASE)
      ENDIF
#endif
C
C**** COMPUTES THE PRODUCT (BULK * THERMAL DEFORMATION)
C
      BULKA=BULKM*STRRA
C
C**** TRANSFORMS CONSTITUTIVE TENSOR FOR LARGE STRAIN ANALYSIS
C
#ifndef restricted
      IF(LARGE.NE.0)
     . CALL FOURTH(DMATX,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,NSTRS,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAP,RCGPI,FINET,TRABE,DCINV,EBASE)
#endif
C
C**** STORE SOME CONSTITUTIVE TENSORS
C
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        DMATE(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMATP(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMAT1(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** SECANT CONSTITUTIVE TENSOR TO BE USED IN STIFFNESS MATRIX
C
#ifndef restricted
      IF(LIQUI.EQ.1) THEN
       IF(LARGE.EQ.0) THEN
        DMAT1(3,3)=DMATX(1,1)/ESTAB
        IF(NTYPE.EQ.4) THEN
         DMAT1(5,5)=DMATX(1,1)/ESTAB
         DMAT1(6,6)=DMATX(1,1)/ESTAB
        ENDIF
c       CALL CALCLT(DMATX,DMAT1,ESTAB)     ! alternative option
       ELSE
        ROTET1=1.0D0/ESTAB
        CALL CALCCT(DMAT1,YOUNG,POISM,ROTET1,    6,
     .              POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .              BULKM1,DISTM1,DAMAP11,DAMAP21,PROPS,
     .              VANIS,ROMTX,TRMTX,
     .              YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
c       CALL CALCLT(DMATX,DMAT1,ESTAB)     ! alternative option
        IF(LARGE.NE.0)
     .   CALL FOURTH(DMAT1,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .               POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .               BULKM1,DISTM1,CPOL1,CPOL2,CPOL3,
     .               EPOL1,EPOL2,EPOL3,
     .               CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .               CPOL9,CPOL10,
     .               CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .               CEMIN,CEMAX,CEETA,
     .               ALAMB,ILAEQ,FVOL1,
     .               RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .               BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .               STILD,ATILD,VANIS,
     .               STRAP,RCGPI,FINET,TRABE,DCINV,EBASE)
       ENDIF           ! large.eq.0
      ENDIF            ! liqui.eq.1
#endif
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMTEP(IKONT)=(1.0D0-DAMAX)*DMAT1(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** TOTAL & INCREMENTAL STRESS STATE PREDICTION
C
      CALL CALCST(SIGMA,DMATX,STRAN,STRAP,STRAT,DAMAX,
     .            DETJM,BULKM,DISTM,CPENI,
     .            FVOL1,STRRA,STRR4,
     .            RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,NSTRS,
     .            FACTJ,STILD,ATILD,SGTOT,PRESG,DESIG,DSTRA,STRS0,
     .            RCGPI,FINET)
C
C**** TANGENT CONJUGATE OF THERMAL DILATATION TENSOR (TANGENT BETA)
C
#ifndef restricted
      IF(ITERME.GT.0) THEN                      ! bidirectional problems
       IF(ITERMP.GT.0) THEN                     ! mech. coupling terms
        CALL CALCBT(DMAPL,COUTD,DMATX,UNOMA,BULKM,DAMAX,
     .              STRAN,STRAP,STRAT,RCGTT,RCGTI,STRRD)
       ENDIF                   ! itermp.gt.0
      ENDIF                    ! iterme.gt.0
#endif
C
C**** YIED CRITERIA CHOICE
C
      GO TO (31,32,33,34,35,36,37,38,39,40,41,42) (NCRIT-30)
C
C**** TRESCA
C
   31 CONTINUE
      CALL RUNEND('ERROR: NCRIT=31 NOT IMPLEMENTED')
      GO TO 100
C
C**** VON MISES
C
   32 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COE multiplies the term (3*J2) in yield function
C     A2COE multiplies the plastic hard. funct. in yield function
C     A3COE multiplies the thermal hard. funct. in yield function
C
      IFVER=INT(PROPS(35))
      IF(IFVER.EQ.1) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.2) THEN
       A1COE=1.0D0
       A2COE=1.0D0-ROTET
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.3) THEN
       A1COE=2.0D0/3.0D0
       A2COE=DSQRT(2.0D0/3.0D0)
       A3COE=A2COE
      ENDIF
      IF(IFVER.EQ.4) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=0.0D0
      ENDIF
      GO TO 100
C
C**** MOHR-COULOMB
C
   33 CONTINUE
      CALL RUNEND('ERROR: NCRIT=33 NOT IMPLEMENTED')
      GO TO 100
C
C**** DRUCKER-PRAGER
C
   34 CONTINUE
      CALL RUNEND('ERROR: NCRIT=34 NOT IMPLEMENTED')
      GO TO 100
C 
C**** J. LUBLINER'S THEORY
C
   35 CONTINUE
      CALL RUNEND('ERROR: NCRIT=35 NOT IMPLEMENTED')
      GO TO 100
C
C**** ABOUAF'S MODEL
C
   36 CONTINUE
      A1COE=1.0D0/ROTET
      A2COE=A1COE-1.0D0
      A3COE=0.0D0
      GO TO 100
C
C**** WEBER-BROWN'S MODEL
C
   37 CONTINUE
      A1COE=6.0D0*(1.0D0-(ROTET)**2.0D0)/(2.0D0+(ROTET)**2.0D0)+0.001D0
      A2COE=(3.0D0*(ROTET)**2.0D0)/(2.0D0+(ROTET)**2.0D0)+0.001D0
      A3COE=0.0D0
      GO TO 100
C
C**** SG CAST IRON
C
   38 CONTINUE
      A1COE=1.0D0
      A2COE=1.0D0
      A3COE=1.0D0
      GO TO 100
C
C**** GREEN SAND
C
   39 CONTINUE
      A1COE=1.0D0
      A2COE=1.0D0
      A3COE=1.0D0
      GO TO 100
C
C**** GREEN SAND
C
   40 CONTINUE
      A1COE=1.0D0
      A2COE=1.0D0
      A3COE=1.0D0
      GO TO 100
C
C**** GURSON
C
   41 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COE multiplies the term (3*J2) in yield function
C     A2COE multiplies the plastic hard. funct. in yield function
C     A3COE multiplies the thermal hard. funct. in yield function
C
      IFVER=INT(PROPS(35))
      IF(IFVER.EQ.1) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.2) THEN
       A1COE=1.0D0
       A2COE=1.0D0-ROTET
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.3) THEN
       A1COE=2.0D0/3.0D0
       A2COE=DSQRT(2.0D0/3.0D0)
       A3COE=A2COE
      ENDIF
      IF(IFVER.EQ.4) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=0.0D0
      ENDIF
      GO TO 100
C
C**** HILL 48
C
   42 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COE multiplies the term (3*J2) in yield function
C     A2COE multiplies the plastic hard. funct. in yield function
C     A3COE multiplies the thermal hard. funct. in yield function
C
      IFVER=INT(PROPS(35))
      IF(IFVER.EQ.1) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.2) THEN
       A1COE=1.0D0
       A2COE=1.0D0-ROTET
       A3COE=1.0D0
      ENDIF
      IF(IFVER.EQ.3) THEN
       A1COE=2.0D0/3.0D0
       A2COE=DSQRT(2.0D0/3.0D0)
       A3COE=A2COE
      ENDIF
      IF(IFVER.EQ.4) THEN
       A1COE=1.0D0
       A2COE=1.0D0
       A3COE=0.0D0
      ENDIF
      GO TO 100
C
  100 CONTINUE
C
C**** FLOW POTENTIAL CHOICE
C
      GO TO (131,132,133,134,135,136,137,138,139,140,141,142) (NCRIP-30)
C
C**** TRESCA
C
  131 CONTINUE
      CALL RUNEND('ERROR: NCRIP=31 NOT IMPLEMENTED')
      GO TO 200
C
C**** VON MISES
C
  132 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COP idem A1COE for flow potential
C     A2COP idem A2COE for flow potential
C     A3COP idem A3COE for flow potential
C
      IPVER=INT(PROPS(51))
      IF(IPVER.EQ.1) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.2) THEN
       A1COP=1.0D0
       A2COP=1.0D0-ROTET
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.3) THEN
       A1COP=2.0D0/3.0D0
       A2COP=DSQRT(2.0D0/3.0D0)
       A3COP=A2COP
      ENDIF
      IF(IPVER.EQ.4) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=0.0D0
      ENDIF
      GO TO 200
C
C**** MOHR-COULOMB
C
  133 CONTINUE
      CALL RUNEND('ERROR: NCRIT=33 NOT IMPLEMENTED')
      GO TO 200
C
C**** DRUCKER-PRAGER
C
  134 CONTINUE
      CALL RUNEND('ERROR: NCRIT=34 NOT IMPLEMENTED')
      GO TO 200
C 
C**** J. LUBLINER'S THEORY
C
  135 CONTINUE
      CALL RUNEND('ERROR: NCRIT=35 NOT IMPLEMENTED')
      GO TO 200
C
C**** ABOUAF'S MODEL
C
  136 CONTINUE
      A1COP=1.0D0/ROTET
      A2COP=A1COE-1.0D0
      A3COP=0.0D0
      GO TO 200
C
C**** WEBER-BROWN'S MODEL
C
  137 CONTINUE
      A1COP=6.0D0*(1.0D0-(ROTET)**2.0D0)/(2.0D0+(ROTET)**2.0D0)+0.001D0
      A2COP=(3.0D0*(ROTET)**2.0D0)/(2.0D0+(ROTET)**2.0D0)+0.001D0
      A3COP=0.0D0
      GO TO 200
C
C**** SG CAST IRON
C
  138 CONTINUE
      A1COP=1.0D0
      A2COP=1.0D0
      A3COP=1.0D0
      GO TO 200
C
C**** GREEN SAND
C
  139 CONTINUE
      A1COP=1.0D0
      A2COP=1.0D0
      A3COP=1.0D0
      GO TO 200
C
C**** GREEN SAND
C
  140 CONTINUE
      A1COP=1.0D0
      A2COP=1.0D0
      A3COP=1.0D0
      GO TO 200
C
C**** GURSON
C
  141 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COP idem A1COE for flow potential
C     A2COP idem A2COE for flow potential
C     A3COP idem A3COE for flow potential
C
      IPVER=INT(PROPS(51))
      IF(IPVER.EQ.1) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.2) THEN
       A1COP=1.0D0
       A2COP=1.0D0-ROTET
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.3) THEN
       A1COP=2.0D0/3.0D0
       A2COP=DSQRT(2.0D0/3.0D0)
       A3COP=A2COP
      ENDIF
      IF(IPVER.EQ.4) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=0.0D0
      ENDIF
      GO TO 200
C
C**** HILL 48
C
  142 CONTINUE
C
C**** CALCULATES CONSTANTS
C
C     Notes:
C
C     A1COP idem A1COE for flow potential
C     A2COP idem A2COE for flow potential
C     A3COP idem A3COE for flow potential
C
      IPVER=INT(PROPS(51))
      IF(IPVER.EQ.1) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.2) THEN
       A1COP=1.0D0
       A2COP=1.0D0-ROTET
       A3COP=1.0D0
      ENDIF
      IF(IPVER.EQ.3) THEN
       A1COP=2.0D0/3.0D0
       A2COP=DSQRT(2.0D0/3.0D0)
       A3COP=A2COP
      ENDIF
      IF(IPVER.EQ.4) THEN
       A1COP=1.0D0
       A2COP=1.0D0
       A3COP=0.0D0
      ENDIF
      GO TO 200
C
  200 CONTINUE
C
      RETURN
      END
