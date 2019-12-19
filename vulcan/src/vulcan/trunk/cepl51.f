      SUBROUTINE CEPL51(SIGMA,DESIG,STRAP,PRESG,
     .                  EBASE,DMATX,PROPS,SGTOT,
     .                  DVOLU,DMTEP,GPCOD,
     .                  STRAN,DSTRA,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  TEMPG,DTEMG,
     .                  PWOEL,SHAPE,DBOLU,THICK,TEINI,
     .                  FPCHL,COUTD,DBASE,SHRIN,FATIG,
     .                  STRS0)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR:
C
C     1) ELASTO-PLASTIC SG CAST IRON
C     2) ELASTO-VISCOPLASTIC SG CAST IRON
C
C     Notes:
C
C     In an isothermal problem, the plastic/viscoplastic module do
C     nothing for IITER=0 when only external loads are applied. This is
C     not true for thermomechanical problems due to the variation of the
C     mechanical properties with the temperature.
C     properties with the temperature.
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
C      EBASE(1): CAPAP: Internal variable: plastic hardening
C      EBASE(2): EFFPL: Effective plastic deformation
C      DBASE(1): --
C      DBASE(2): --
C      COUTD: Coupling variable (for ITERME > 0)
C      FATIG: Fatigue variables
C
C      Viscoplastic model:
C      STRAP(ISTR1):    Internal variable: viscoplastic deformation
C      EBASE(1): CAPAP: Internal variable: plastic hardening
C      EBASE(2): EFFPL: Effective plastic deformation
C      DBASE(1): --
C      DBASE(2): --
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
     .          FPCHL(NFPCH,*), FATIG(*)
      DIMENSION STRS0(*)
C
      TWOPI=6.283185307179586
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
      CALL SOLIDI1(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .             DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .             KUNL4,RETEN,RETEP,TEMPG,CCEROF,
     .             CCEROA,DESIG,
     .             DSTRA,STRAN,DTEMG,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,
     .             ROTES,NELASF,NELASA,FPCHL,SHAPE,DBASE,SHRIN,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             GPCOD)
C
C**** COMPUTES PARTICULAR VARIABLES OF VISCOPLASTIC MODEL
C
      IF(IPEP2.EQ.30)
     . CALL SOL132(PROPS,TEMPG,
     .             VISCO,EXPON,
     .             VISCOF,VISCOA,EXPONF,EXPONA,
     .             AZABA,AZABB,AZABC,
     .             AZABAF,AZABAA,AZABBF,AZABBA,AZABCF,AZABCA)
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
C**** COMPUTES VARIABLES ASSOCIATED WITH THE POROSITY
C
      IF(IPEP2.EQ.20.OR.IPEP2.EQ.30)
     . CALL SOLPOR(PROPS,TEMPG,CPOEF,CPOE1,CPOE2)
C
C**** ELASTO-PLASTIC ALGORITHM
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.20)
     .  call runend('error in cepl51')
c    .  CALL PLASTC1(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
c    .               PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
c    .               ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
c    .               DMATP,STRAN,ALFAT,BETAT,
c    .               GAMAT,RETEN,RETEP,CCEROF,CCEROA,
c    .               A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
c    .               DEEPP,DTEMG,DCCER,
c    .               CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA)
c    .               XJACM,XJACI,XJA3M,XJA3I,DETJM,
c    .               FATIG)
      ENDIF
C
C**** ELASTO-VISCOPLASTIC ALGORITHM
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.30)
     .  CALL VISCOC1(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
     .               PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
     .               ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .               DMATP,STRAN,ALFAT,BETAT,
     .               GAMAT,RETEN,RETEP,CCEROF,CCEROA,
     .               A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .               DEEPP,DTEMG,DCCER,
     .               CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA,
     .               CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA,
     .               VISCOF,VISCOA,EXPONF,EXPONA,
     .               ROTES,NELASF,NELASA,
     .               DBASE,CDAMA,CFRAC,
     .               CPOEF,CPOE1,CPOE2,
     .               AZABAF,AZABBF,AZABCF,AZABAA,AZABBA,AZABCA,
     .               XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .               FATIG)
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
#endif
      END
