      SUBROUTINE CEPL36(SIGMA,DESIG,STRAP,PRESG,
     .                  EBASE,DMATX,PROPS,SGTOT,
     .                  DVOLU,DMTEP,GPCOD,
     .                  STRAN,DSTRA,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  TEMPG,DTEMG,VANIS,
     .                  PWOEL,SHAPE,DBOLU,THICK,TEINI,
     .                  FPCHL,COUTD,DBASE,SHRIN,FATIG,
     .                  STRS0)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR:
C
C     1) ELASTO-PLASTIC MATERIALS
C     2) ELASTO-VISCOPLASTIC MATERIALS
C
C     Notes:
C
C     In an isothermal problem, the plastic/viscoplastic module do
C     nothing for IITER=0 when only external loads are applied. This is
C     not true for thermomechanical problems due to the variation of the
C     mechanical properties with the temperature.
C     Nevertheless, some thermomechanical problems do not converge when
C     the plastic/viscoplastic module is accessed for IITER=0 (very
C     large irreversible effects are produced). Then, a purely elastic
C     behaviour is assumed for this iteration.
C
C     FUTURE WORK:
C
C     Degradation and densification problems can be considered by means
C     of IPEPE variable (as particular plastic or viscoplastic models
C     or as a new class of problems included in PROPS(4))
C
C***********************************************************************
C
C     Yielding criteria/Flow Potentials:
C
C     NCRIT = 31   Tresca
C     NCRIT = 32   Von-Mises
C     NCRIT = 33   Mohr-Coulomb- (Modif. Oller)
C     NCRIT = 34   Drucker-Prager
C     NCRIT = 35   Lubliner-Oller
C     NCRIT = 36   Abouaf (cap-model)
C     NCRIT = 37   Weber-Brown (cap-model)
C     NCRIT = 38   SG Cast iron
C     NCRIT = 39   Green sand (yield function)
C     NCRIT = 40   Green sand (flow potential)
C     NCRIT = 41   Gurson
C     NCRIT = 42   Hill 48
C
C     Notes:
C           A new NCRIT is necessary if the term
C           \partial F / \partial \sigma is not proportional to the
C           same term of the existing NCRIT (this proportionality is
C           given by the A1COE & A2COE (Yield criteria) and A1COP &
C           A2COP (flow potentials).
C
C           Different versions of a particular NCRIT can exist: they
C           differ in the A1COE & A2COE (Yield criteria) and A1COP &
C           A2COP (flow potentials). Therefore, these different versions
C           give associated models.
C
C     More specifically, the Yielding criteria/Flow Potentials are:
C
C     NCRIT=32   F=A1COE*\sqrt(3*J2)-C^th-A2COE*C^p
C     NCRIP=32   G=A1COP*\sqrt(3*J2)-C^th-A2COP*C^p
C
C***********************************************************************
C
C     Input variables:
C      
C     SIGMA(ISTR1):       Current predictor stress tensor
C     DESIG(ISTRE):       Increment of predictor stress tensor
C     PRESG(ISTR1):       Last converged stress tensor
C     DMATX(NSTRS,NSTRS): Elastic constitutive tensor
C     PROPS(NPROP):       Material properties
C     DVOLU(NGAUS):       Volume of the integration point
C
C     TEINI:              Initial Temperature
C
C     Output variables:
C
C     SGTOT(NSTRS): Current stress tensor
C     DMTEP(NKOST): Elasto-plastic constitutive tensor
C
C     Input/Output variables:
C
C      Plastic model:
C      STRAP(ISTR1):    Internal variable: plastic deformation
C      EBASE(1): EFFPL: Effective plastic deformation
C      EBASE(2): CAPAP: Internal variable: plastic hardening
C      EBASE(3:3+NSTR1-1): BACKA: Associate variable to back stress
C      DBASE(1): DAMAG or TAUAM: Damage or tau-
C      DBASE(2): TAUAP: tau+
C      COUTD: Coupling variable (for ITERME > 0)
C      FATIG: Fatigue variables
C
C      Viscoplastic model:
C      STRAP(ISTR1):    Internal variable: viscoplastic deformation
C      EBASE(1): EFFPL: Effective plastic deformation
C      EBASE(2): CAPAP: Internal variable: plastic hardening
C      EBASE(3:3+NSTR1-1): BACKA: Associate variable to back stress
C      DBASE(1): DAMAG or TAUAM: Damage or tau-
C      DBASE(2): TAUAP: tau+
C      COUTD: Coupling variable (for ITERME > 0)
C      FATIG: Fatigue variables
C
C      Additional variable:
C
C      KUNLD: Load indicator to be considered in stif30.f to update the
C             constitutive tensor
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
      DIMENSION DESIG(*),       DMATX(NSTRS,*), EBASE(*),  DSTRA(*),
     .          PRESG(*),       PROPS(*),       SGTOT(*),  STRAN(*),
     .          SIGMA(*),       STRAP(*),       DMTEP(*),  DVOLU(*),
     .          XJACM(NDIME,*), XJACI(NDIME,*), GPCOD(*)
C
      DIMENSION DMATE(6,6),     DMATP(6,6),     SIGM0(6),  DM(6,6)
      DIMENSION DEEPP(6),       DMAPL(6),       DEETH(6)
      DIMENSION PWOEL(NNODE),   SHAPE(NNODE),   DBASE(*),
     .          FPCHL(NFPCH,*), FATIG(*),       VANIS(NANIV,*)
      DIMENSION STRS0(*),       TRMTX(3,3)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)

C
      TWOPI=6.283185307179586D0
C
C**** INITIALISES LOAD INDICATOR
C
      KUNLD=1
      KUNL1=1
      KUNL2=1
      KUNL3=1
      KUNL4=1
C
C**** LOOK FOR MODEL DEFINITION
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IPEP1=100
      IF(ISOTR.EQ.1) IPEP1=200
      IPEP2=INT(PROPS(2))   ! 10=elastic; 20=plastic; 30=viscoplastic;
C                           ! 40=hyperelastic
      IPEP3=INT(PROPS(3))   ! 1=no temp.-dependent; 2=temp.-depend.
      IPEPE=IPEP1+IPEP2+IPEP3
C
C**** COMPUTES BASIC VARIABLES & INITIALISES PARAMETERS
C
      CALL SOLIDI(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .            DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .            KUNL4,RETEN,RETEP,TEMPG,CCERO,DESIG,
     .            CCEROM,CCEROP,BULKM,DISTM,
     .            DSTRA,STRAN,DTEMG,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,FPCHL,
     .            SHAPE,DBASE,SHRIN,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            GPCOD,VANIS,STRS0,TRMTX,EBASE)
C
C**** COMPUTES PARTICULAR VARIABLES OF VISCOPLASTIC MODEL
C
#ifndef restricted
      IF(IPEP2.EQ.30)
     . CALL SOL132(PROPS,TEMPG,
     .             VISCO,EXPON,
     .             VISCOF,VISCOA,EXPONF,EXPONA,
     .             AZABA,AZABB,AZABC,
     .             AZABAF,AZABAA,AZABBF,AZABBA,AZABCF,AZABCA)
#endif
C
C**** COMPUTES VARIABLES ASSOCIATED WITH THE ISOTROPIC & KINEMATIC
C     HARDENING
C
      CALL SOLISO(PROPS,TEMPG,
     .            CCOEF,CCOEB,CCOEQ,
     .            CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA,
     .            CKOEF,CKOEB,CKOEQ,
     .            CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA)
C
C**** COMPUTES PARTICULAR VARIABLES OF DAMAGE MODEL
C
#ifndef restricted
      IF(IPEP2.EQ.20.OR.IPEP2.EQ.30)
     . CALL SOLDAM(PROPS,TEMPG,
     .             CDAMA,CFRAC)
#endif
C
C**** COMPUTES PARTICULAR VARIABLES OF YIELD FUNCTION, FLOW POTENTIAL &
C     DAMAGE MODEL
C
      CALL SOLYFD(PROPS,TEMPG,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,CCOB7,CCOB8,
     .            CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6)
C
C**** COMPUTES VARIABLES ASSOCIATED WITH THE POROSITY
C
#ifndef restricted
      IF(IPEP2.EQ.20.OR.IPEP2.EQ.30)
     . CALL SOLPOR(PROPS,TEMPG,CPOEF,CPOE1,CPOE2)
#endif
C
C**** ELASTO-PLASTIC ALGORITHM
C
#ifndef restricted
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.20)
     .  CALL PLASTC(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
     .              PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
     .              ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .              DMATP,STRAN,ALFAT,BETAT,BULKM,DISTM,
     .              GAMAT,RETEN,RETEP,CCERO,CCEROM,CCEROP,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              DEEPP,DTEMG,DCCER,
     .              CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .                          DBASE,CDAMA,CFRAC,COUTH,
     .              CPOEF,CPOE1,CPOE2,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,NSTRS)
      ENDIF
C
C**** ELASTO-VISCOPLASTIC ALGORITHM
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.30) THEN
        CALL VISCOC(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
     .              PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
     .              ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .              DMATP,STRAN,ALFAT,BETAT,BULKM,DISTM,
     .              GAMAT,RETEN,RETEP,CCERO,CCEROM,CCEROP,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              DEEPP,DTEMG,DCCER,
     .              CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .              VISCO,EXPON,DBASE,CDAMA,CFRAC,COUTH,
     .              CPOEF,CPOE1,CPOE2,
     .              AZABA,AZABB,AZABC,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,NSTRS,YIELD)
       ENDIF
      ENDIF
C
C**** COMPUTES COUPLING TERM
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(ITERME.GT.0) THEN              ! bidirectional coupled problem
        IF(ITERMP.GT.0)
     .   CALL PLHEAT(PWOEL,SHAPE,DSTRA,DEEPP,DEETH,DMAPL,SGTOT,TEMPG,
     .               DBOLU,COUTD,COUTH)
       ENDIF
      ENDIF
#endif
C
C**** CORRECTION OF THE "ZZ" COMPONENT OF THE DEFORMATION TENSOR FOR
C     THE PLANE STRESS CASE
C
      IF(NTYPE.EQ.1) CALL PLANES(STRAN,STRAP)
C
C**** DEFINES THE LOAD INDICATOR
C
      IF((KUNL1.EQ.0).OR.(KUNL2.EQ.0).OR.(KUNL3.EQ.0)
     .               .OR.(KUNL4.EQ.0)) KUNLD=0
C
      RETURN
      END
