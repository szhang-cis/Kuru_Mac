      SUBROUTINE CEPL53(SIGMA,DESIG,STRAP,PRESG,
     .                  EBASE,DMATX,PROPS,SGTOT,
     .                  DVOLU,DMTEP,GPCOD,
     .                  STRAN,DSTRA,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  TEMPG,DTEMG,VANIS,
     .                  PWOEL,SHAPE,DBOLU,THICK,TEINI,
     .                  FPCHL,COUTD,DBASE,SHRIN,FATIG,
     .                  STRS0,VDUAL)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR THE DUAL PHASE STEEL
C     MODEL
C
C***********************************************************************
C
C     To be revised!!!
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
      DIMENSION STRS0(*),       VDUAL(*)
      DIMENSION TRMTX(3,3)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
      DIMENSION DMATXF(6,6),    DMATXA(6,6)
      DIMENSION DMATEF(6,6),    DMATEA(6,6)
      DIMENSION DMATPF(6,6),    DMATPA(6,6)
      DIMENSION STRANF(6),      STRANA(6)
      DIMENSION SIGMAF(6),      SIGMAA(6)
      DIMENSION SGTOTF(6),      SGTOTA(6)  ! iterative history variables
      DIMENSION STRAPF(6),      STRAPA(6)  ! iterative history variables
      DIMENSION EBASEF(2),      EBASEA(2)  ! iterative history variables
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
C**** RECOVERS STRESS TENSOR
C
      DO ISTRS=1,NSTRS
       SGTOT(ISTRS)=SIGMA(ISTRS)
      ENDDO
C
C**** POINTERS FOR THE DUAL PHASE STEEL MODEL VARIABLES
C
      IDP1=1                ! ferrite (f) plastic strain
      IDP2=IDP1+NSTR1       ! martensite (m) plastic strain
      IDP3=IDP2+NSTR1       ! ferrite constitutive tensor
      IDP4=IDP3+NKOST       ! martensite constitutive tensor
      IDP5=IDP4+NKOST       ! ferrite stress
      IDP6=IDP5+NSTR1       ! martensite stress
      IDP7=IDP6+NSTR1       ! ferrite effective plastic strain
      IDP8=IDP7+1           ! ferrite hardening function
C                           ! other ferrite variables if necessary
      IDP9=IDP8+1           ! martensite effective plastic strain
      IDP10=IDP9+1          ! martensite hardening function
C                           ! other martensite variables if necessary
      IDP11=IDP10+1         ! ferrite strain
      IDP12=IDP11+NSTR1     ! martensite strain
      IDP13=IDP12+NSTR1     ! strain partition coefficient
C
C**** RECOVERS STRAINS, STRESSES & INTERNAL VARIABLES (f & m)
C
      DO ISTRS=1,NSTRS
       STRAPF(ISTRS)=VDUAL(IDP1-1+ISTRS)      ! plastic strain
       STRAPA(ISTRS)=VDUAL(IDP2-1+ISTRS)
       SGTOTF(ISTRS)=VDUAL(IDP5-1+ISTRS)      ! stress
       SGTOTA(ISTRS)=VDUAL(IDP6-1+ISTRS)
       STRANF(ISTRS)=VDUAL(IDP11-1+ISTRS)     ! strain
       STRANA(ISTRS)=VDUAL(IDP12-1+ISTRS)
      ENDDO
C
      EBASEF(1)=VDUAL(IDP7)                   ! effective plastic strain
      EBASEA(1)=VDUAL(IDP9)
      EBASEF(2)=VDUAL(IDP8)                   ! hardening stress
      EBASEA(2)=VDUAL(IDP10)
C
      DSISF=VDUAL(IDP13)
C
C**** RECOVERS FERRITE & MARTENSITE FRACTIONS
C
      EFFPL=0.0D0
      CALL QPARTT(PROPS,EFFPL,
     .            QPART,FRACTF,FRACTA)
C
C**** FERRITE & MARTENSITE STRAINS DERIVED FROM MIXING LAW & STRAIN
C     PARTITION COEFFICIENT
C
      DO ISTRS=1,NSTRS
       STRANF(ISTRS)=STRANF(ISTRS)+DSTRA(ISTRS)*(1.0D0+DSISF)
       STRANA(ISTRS)=STRANA(ISTRS)+DSTRA(ISTRS)*
     .                            (1.0D0-FRACTF/FRACTA*DSISF)
      ENDDO
C
C**** COMPUTES BASIC VARIABLES & INITIALISES PARAMETERS (f & m)
C
C     Note: EBASE is not used now in solidi3.f
C
      CALL SOLIDI3(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .             DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .             KUNL4,RETEN,RETEP,TEMPG,CCERO,DESIG,
     .             CCEROM,CCEROP,BULKM,DISTM,
     .             DSTRA,STRAN,DTEMG,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,FPCHL,
     .             SHAPE,DBASE,SHRIN,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             GPCOD,VANIS,STRS0,TRMTX,EBASE,
     .             CCEROF,CCEROA,
     .             DMATXF,DMATXA,DMATEF,DMATEA,DMATPF,DMATPA,
     .             STRANF,STRANA,SIGMAF,SIGMAA,
     .             STRAPF,STRAPA,
     .             VDUAL(IDP3),VDUAL(IDP4),              ! DMTEPF,DMTEPA
     .             SGTOTF,SGTOTA)
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
C**** COMPUTES PARTICULAR VARIABLES OF DAMAGE MODEL
C
      IF(IPEP2.EQ.20.OR.IPEP2.EQ.30)
     . CALL SOLDAM(PROPS,TEMPG,
     .             CDAMA,CFRAC)
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
      IF(IPEP2.EQ.20.OR.IPEP2.EQ.30)
     . CALL SOLPOR(PROPS,TEMPG,CPOEF,CPOE1,CPOE2)
C
C**** ELASTO-PLASTIC ALGORITHM
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.20) call runend('ERROR: plastic model not implem.')
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
     .              FATIG,TRMTX,    6)
      ENDIF
C
C**** ELASTO-VISCOPLASTIC ALGORITHM (f & m)
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       IF(IPEP2.EQ.30) THEN
        CALL VISCOC(SIGMAF,DESIG,STRAPF,PRESG,EBASEF,DMATXF,
     .              PROPS,SGTOTF,VDUAL(IDP3),DMATEF,KUNL2,ALENG,
     .              ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .              DMATPF,STRANF,ALFAT,BETAT,BULKM,DISTM,
     .              GAMAT,RETEN,RETEP,CCEROF,CCEROM,CCEROP,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              DEEPP,DTEMG,DCCER,
     .              CCOEFF,CCOEBF,CCOEQF,CKOEF,CKOEB,CKOEQ,    ! revisar
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .              VISCOF,EXPONF,DBASE,CDAMA,CFRAC,COUTH,
     .              CPOEF,CPOE1,CPOE2,
     .              AZABA,AZABB,AZABC,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,    6,YIELDF)
C
        CALL VISCOC(SIGMAA,DESIG,STRAPA,PRESG,EBASEA,DMATXA,
     .              PROPS,SGTOTA,VDUAL(IDP4),DMATEA,KUNL2,ALENG,
     .              ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .              DMATPA,STRANA,ALFAT,BETAT,BULKM,DISTM,
     .              GAMAT,RETEN,RETEP,CCEROA,CCEROM,CCEROP,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              DEEPP,DTEMG,DCCER,
     .              CCOEFA,CCOEBA,CCOEQA,CKOEF,CKOEB,CKOEQ,    ! revisar
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .              VISCOA,EXPONA,DBASE,CDAMA,CFRAC,COUTH,
     .              CPOEF,CPOE1,CPOE2,
     .              AZABA,AZABB,AZABC,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,    6,YIELDA)
       ENDIF
      ENDIF
C
C**** COMPUTES PARTITION PARAMETER & STRAIN PARTITION COEFFICIENT
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       YIELD=FRACTF*YIELDF+FRACTA*YIELDA ! DP equivalent stress
       EBASE(1)=FRACTF*EBASEF(1)+FRACTA*EBASEA(1)
       EFFPL=EBASE(1)                    ! DP effective plastic strain
       CALL QPARTT(PROPS,EFFPL,
     .             QPART,FRACTF,FRACTA)
       DSISF=(YIELD-YIELDF)/QPART
C
       IF(DSISF.LT.0.0D0) DSISF=0.0D0             ! lower bound
       XX=FRACTF/FRACTA*DSISF
       IF(XX.GT.1.0D0) DSISF=FRACTA/FRACTF        ! upper bound
      ENDIF
C
C**** COMPUTES STRESS TENSOR WITH MIXING LAW (ML)
C
      DO ISTRS=1,NSTRS
       SGTOT(ISTRS)=FRACTF*SGTOTF(ISTRS)+FRACTA*SGTOTA(ISTRS)
      ENDDO
C
C**** COMPUTES VIA ML APPROXIMATE ELASTO-PLASTIC MATERIAL TENSOR
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMTEP(IKONT)=FRACTF*VDUAL(IDP3-1+IKONT)+
     .               FRACTA*VDUAL(IDP4-1+IKONT)
       ENDDO
      ENDDO
C
C**** STORES STRAINS, STRESSES & INTERNAL VARIABLES (f & m)
C
      IF(IITER.GT.ITEPL) THEN            ! see note
       DO ISTRS=1,NSTRS
        VDUAL(IDP1-1+ISTRS)=STRAPF(ISTRS)
        VDUAL(IDP2-1+ISTRS)=STRAPA(ISTRS)
        VDUAL(IDP5-1+ISTRS)=SGTOTF(ISTRS)
        VDUAL(IDP6-1+ISTRS)=SGTOTA(ISTRS)
        VDUAL(IDP11-1+ISTRS)=STRANF(ISTRS)
        VDUAL(IDP12-1+ISTRS)=STRANA(ISTRS)
       ENDDO
C
       VDUAL(IDP7) =EBASEF(1)
       VDUAL(IDP9) =EBASEA(1)
       VDUAL(IDP8) =EBASEF(2)
       VDUAL(IDP10)=EBASEA(2)
C
       VDUAL(IDP13)=DSISF
C
C**** COMPUTES VIA ML APPROXIMATE PLASTIC STRAIN TENSOR, EFFECTIVE
C     PLASTIC STRAIN & HARDENING STRESS
C
       DO ISTRS=1,NSTRS
        STRAP(ISTRS)=FRACTF*STRAPF(ISTRS)+FRACTA*STRAPA(ISTRS)
       ENDDO
C
       EBASE(1)=FRACTF*EBASEF(1)+FRACTA*EBASEA(1)
       EBASE(2)=FRACTF*EBASEF(2)+FRACTA*EBASEA(2)
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
