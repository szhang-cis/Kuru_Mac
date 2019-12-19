      SUBROUTINE SOLIDI1(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .                   DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .                   KUNL4,RETEN,RETEP,TEMPG,CCEROF,
     .                   CCEROA,DESIG,
     .                   DSTRA,STRAN,DTEMG,
     .                   A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                   STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,
     .                   ROTES,NELASF,NELASA,FPCHL,SHAPE,DBASE,SHRIN,
     .                   XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                   GPCOD)
#ifndef restricted
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
C     INDEX OF VARIABLES:
C
C     ROTET: 1-f^ls_pc (f^ls_pc: liquid-solid phase-change function)
C     ROTES: 1-f^ga_pc (f^ga_pc: solid-solid phase change function)
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
C     LARGE DISPLACEMENTS & STRAINS (not implemented yet!)
C
C     IFREN=1  psi=
C                                     hyperelastic linear material model
C
C     IFREN=2  psi=
C                                     hyperelastic linear spatial model
C
C     IFREN=3  psi=
C                                     elastic model
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
C                       Simo model (not used here)
C
C     IFREN=8  psi=
C                       identical to IFREN=2 computing in a different
C                       but absolutely equivalent form the stress and
C                       constitutive tensors (not used here)
C
C     IFREN=9  psi=
C                       hyperelastic linear material-spatial model
C
C     IFREN=10 psi=
C                       idem IFREN=9 with non-linear thermal deformation
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
      DIMENSION STRATF(6),      STRATA(6)
      DIMENSION DMATF(6,6),     DMATA(6,6)
      DIMENSION FPCHL(NFPCH,*), SHAPE(*),       DBASE(*)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*), GPCOD(*)
      DIMENSION RCGTT(6),       RCGTI(6),       ENEYE(6), ALAMB(3)
      DIMENSION SIGMAF(6),      SIGMAA(6)
      DIMENSION DEETHF(6),      DEETHA(6)
      DIMENSION DCOFA(6,6),     UNORA(6,6),     STILD(6)
      DIMENSION DCINV(6,6)
      DIMENSION ROMTX(3,3)
C
C**** INITIALISES SOME PARAMETERS
C
      KUNL4=0
C
      DO 20 ISTR1=1,NSTR1
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
      IF(ITERME.GE.0) THEN                  ! coupled problems
       CALL ALPHAT(TEINI,PROPS,
     .             ALPHI,ALPHIF,ALPHIA)
       TEMPO=TEMPG-DTEMG
       CALL ALPHAT(TEMPO,PROPS,
     .             ALPHO,ALPHOF,ALPHOA)
       CALL ALPHAT(TEMPG,PROPS,
     .             ALPHA,ALPHAF,ALPHAA)
      ELSE                                  ! uncoupled problems
       TEREF=0.0D0
       TEINI=0.0D0
       ALPHIF=0.0D0
       ALPHIA=0.0D0
       ALPHOF=0.0D0
       ALPHOA=0.0D0
       ALPHAF=0.0D0
       ALPHAA=0.0D0
      ENDIF
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
      NELASF=1
      NELASA=0
C
      IF(ITERME.GE.0)                       ! coupled problems
     . CALL ROTETT(ROTET,DEROT,LIQUI,DELLS,
     .             TEMPG,DTEMG,SHAPE,FPCHL,
     .             ROTES,DEROS,LIQUS,DELSS,
     .             STRRB,STRRC,
     .             YOUNG,CCERO,CCEROM,CCEROP)
C
      IF(LIQUS.EQ.0) THEN
       IF(ROTES.GT.0.0.AND.ROTES.LT.1.0) THEN
        NELASF=1
        NELASA=1
       ENDIF
      ELSE
       NELASF=0
       NELASA=1
      ENDIF
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
      CALL CALCPD(UNOMA,STRATF,DEETHF,STRRBF,STRRCF,STRREF,
     .            RCGTT)
      CALL CALCPD(UNOMA,STRATA,DEETHA,STRRBA,STRRCA,STRREA,
     .            RCGTT)
C
C**** COMPUTES THE THERMAL DEFORMATION
C
      CALL CALCTD(PROPS,UNOMA,STRATF,DEETHF,STRRAF,STRRDF,
     .            TEINI,TEMPO,TEMPG,DTEMG,
     .            ALPHIF,ALPHOF,ALPHAF,
     .            STRAN,STRAP,RCGTT,YOUNGF,POISMF)
      CALL CALCTD(PROPS,UNOMA,STRATA,DEETHA,STRRAA,STRRDA,
     .            TEINI,TEMPO,TEMPG,DTEMG,
     .            ALPHIA,ALPHOA,ALPHAA,
     .            STRAN,STRAP,RCGTT,YOUNGA,POISMA)
C
C**** COMPUTES LIQUID & SHRINKAGE DEFORMATIONS
C
      IF(LIQUI.EQ.1) CALL CALCLD(UNOMA,STRAT,STRAN,STRAP,RCGTT,
     .                           XJACM,XJACI,XJA3M,XJA3I,
     .                           SHRIN)
C
C**** SECANT CONSTITUTIVE TENSOR TO INTEGRATE THE STRESSES
C
      CALL CALCCT(DMATF,YOUNGF,POISMF,ROTET,    6,    ! POISCiF not used
     .            POISC1F,POISC2F,POISC3F,D1PLSF,D2PLSF,D3PLSF,
     .            BULKMF,DISTMF,DAMAP1,DAMAP2,PROPS,
     .            DUMMY,DUMMY,DUMMY,
     .            YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
      CALL CALCCT(DMATA,YOUNGA,POISMA,ROTET,    6,    ! POISCiA not used
     .            POISC1A,POISC2A,POISC3A,D1PLSA,D2PLSA,D3PLSA,
     .            BULKMA,DISTMA,DAMAP1,DAMAP2,PROPS,
     .            DUMMY,DUMMY,DUMMY,
     .            YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
C
C**** STORE MIXED CONSTITUTIVE TENSOR
C
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        DMATX(ISTRS,JSTRS)=ROTES*DMATF(ISTRS,JSTRS)+
     .              (1.0D0-ROTES)*DMATA(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** COMPUTES MIXED BULK AND SHEAR MODULI
C
      BULKM=ROTES*BULKMF+(1.0D0-ROTES)*BULKMA
      DISTM=ROTES*DISTMF+(1.0D0-ROTES)*DISTMA
C
C**** COMPUTES THE MIXED PRODUCT (BULK * THERMAL DEFORMATION)
C
      BULKA=ROTES*BULKMF*STRRAF+(1.0D0-ROTES)*BULKMA*STRRAA
C
C**** DEALS WITH PLANE STRESS
C
      IF(NTYPE.EQ.1) THEN
       D1PLS=ROTES*D1PLSF+(1.0D0-ROTES)*D1PLSA
       D2PLS=ROTES*D2PLSF+(1.0D0-ROTES)*D2PLSA
       POISC1=-D2PLS/D1PLS
       POISC2=POISC1
       POISC3=0.0D0
      ENDIF
C
C**** TRANSFORMS CONSTITUTIVE TENSOR FOR LARGE STRAIN ANALYSIS
C
      IF(LARGE.NE.0)
     . CALL FOURTH(DMATX,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             DUMMY,DUMMY,
     .             DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .             DUMMY,DUMMY,DUMMY,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,NSTRS,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,DUMMY,
     .             DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY)
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
      IF(LIQUI.EQ.1) THEN
       if(large.ne.0) then
        ROTET1=0.05D0
        CALL CALCCT(DMAT1,YOUNG,POISM,ROTET1,    6,
     .              POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .              BULKM1,DISTM1,DAMAP11,DAMAP21,PROPS,
     .              DUMMY,DUMMY,DUMMY,
     .              YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
        IF(LARGE.NE.0)
     .   CALL FOURTH(DMAT1,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .               POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .               BULKM1,DISTM1,CPOL1,CPOL2,CPOL3,
     .               EPOL1,EPOL2,EPOL3,
     .               CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .               DUMMY,DUMMY,
     .               DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     .               DUMMY,DUMMY,DUMMY,
     .               ALAMB,ILAEQ,FVOL1,
     .               RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .               BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .               STILD,ATILD,DUMMY,
     .               DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY)
       else
        DMAT1(3,3)=DMATX(1,1)/ESTAB
        IF(NTYPE.EQ.4) THEN
         DMAT1(5,5)=DMATX(1,1)/ESTAB
         DMAT1(6,6)=DMATX(1,1)/ESTAB
        ENDIF
       endif
      ENDIF            ! liqui.eq.1
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
C**** TOTAL STRESS STATE PREDICTION
C
      CALL CALCST(SIGMAF,DMATF,STRAN,STRAP,STRATF,DAMAX,
     .            DETJM,BULKMF,DISTMF,CPENI,
     .            FVOL1,STRRAF,STRR4F,
     .            RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,    6,
     .            FACTJ,STILD,ATILD,SGTOT,PRESG,DESIG,DSTRA,DUMMY,
     .            DYMMY,DUMMY)
      CALL CALCST(SIGMAA,DMATA,STRAN,STRAP,STRATA,DAMAX,
     .            DETJM,BULKMA,DISTMA,CPENI,
     .            FVOL1,STRRAA,STRR4A,
     .            RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,    6,
     .            FACTJ,STILD,ATILD,SGTOT,PRESG,DESIG,DSTRA,DUMMY,
     .            DYMMY,DUMMY)
C
C**** COMPUTES MIXED STRESS
C
      DO ISTRS=1,NSTRS
       SIGMA(ISTRS)=ROTES*SIGMAF(ISTRS)+(1.0D0-ROTES)*SIGMAA(ISTRS)
      ENDDO
C
C**** DEALS WITH PLANE STRESS
C
      IF(NTYPE.EQ.1) STRR4=ROTES*STRR4F+(1.0D0-ROTES)*STRR4A
C
C**** UPDATE TOTAL & INCREMENTAL STRESS
C
      DO ISTR1=1,NSTR1
       SGTOT(ISTR1)=SIGMA(ISTR1)*(1.0D0-DAMAX)
       DESIG(ISTR1)=SGTOT(ISTR1)-PRESG(ISTR1)   ! presg= ^t strsg
      ENDDO
C
C**** TANGENT CONJUGATE OF THERMAL DILATATION TENSOR (TANGENT BETA)
C
      IF(ITERME.GT.0) THEN                      ! bidirectional problems
       IF(ITERMP.GT.0) THEN                     ! mech. coupling terms
C
        call runend('ERROR: ITERME > 0 & ITERMP > 0 NOT IMPLEMENTED')
C
        CALL CALCBT(DMAPL,COUTD,DMATX,UNOMA,BULKM,DAMAX,
     .              STRAN,STRAP,STRAT,RCGTT,RCGTI,STRRD)
       ENDIF                   ! itermp.gt.0
      ENDIF                    ! iterme.gt.0
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
#endif
      RETURN
      END
