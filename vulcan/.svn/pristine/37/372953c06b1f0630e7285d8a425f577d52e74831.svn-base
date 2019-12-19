      SUBROUTINE VISCOC1(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
     .                   PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
     .                   ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .                   DMATP,STRAN,ALFAT,BETAT,
     .                   GAMAT,RETEN,RETEP,CCEROF,CCEROA,
     .                   A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                   DEEPP,DTEMG,DCCER,
     .                   CCOEFF,CCOEFA,CCOEBF,CCOEBA,CCOEQF,CCOEQA,
     .                   CKOEFF,CKOEFA,CKOEBF,CKOEBA,CKOEQF,CKOEQA,
     .                   VISCOF,VISCOA,EXPONF,EXPONA,
     .                   ROTES,NELASF,NELASA,
     .                   DBASE,CDAMA,CFRAC,
     .                   CPOEF,CPOE1,CPOE2,
     .                   AZABAF,AZABBF,AZABCF,AZABAA,AZABBA,AZABCA,
     .                   XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                   FATIG)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR THE SG CAST IRON MODEL
C
C***********************************************************************
C
C**** MORE SPECIFICALLY, ROUTINE viscoc1.f DOES:
C
C     1) ENTERS VISCOPLASTIC MODULE vimodl.f
C     2) EVALUATES THE ELASTO-VISCOPLASTIC CONSTITUTIVE TENSOR
C
C***********************************************************************
C
C**** THE VISCOPLASTIC MODULE vimodl.f COMPUTES:
C
C     1) EFFECTIVE STRESS (INVARIANTS)
C     2) VISCOPLASTIC FLUX
C     3) INTERNAL VARIABLES (i.e., INTERNAL VARIABLES WITHOUT
C        VISCOPLASTIC DEFORMATION IN WHAT FOLLOWS)
C     4) TOTAL HARDENING FUNCTION
C     5) DENOMINATOR OF VISCOPLASTIC PARAMETER LAMBDA ("ABETA")
C     6) INCREMENT OF VISCOPLASTIC PARAMETER LAMBDA ("DLAMD")
C
C     FIRST PASS:   ANALYSES THE STATE OF THE POINT BY MEANS OF THE
C                   YIELD FUNCTION (ELASTIC OR VISCOPLASTIC) 
C     SECOND PASS:  COMPUTES FLUX, INTERNAL VARIABLES & LAMBDA MULTIP.
C     THIRD PASS:  -CHECKS CONVERGENCE OF VISCOPLASTIC MODEL
C                  -VERIFICATES THE VISCOPLASTIC DISSIPATION
C     FOURTH PASS:  STORES INTERNAL VARIABLES
C
C***********************************************************************
C
C**** PARAMETERS:
C    
C     KFLUI=0  - DOES NOT COMPUTE FLUX
C     KFLUI=1  - COMPUTES FLUX
C
C     INDEX=1  - COMPUTES EFFECTIVE STRESS (INVARIANTS), INTERNAL 
C                VARIABLES & TOTAL HARDENING FUNCTION
C     INDEX=2  - COMPUTES EFFECTIVE STRESS (INVARIANTS) & TOTAL
C                HARDENING FUNCTION
C     INDEX=3  - STORES INTERNAL VARIABLES
C
C     IVIRI=0    DOES NOT EVALUATE PSEUDO CONSTITUTIVE TENSOR
C     IVIRI=1    EVALUATES PSEUDO CONSTITUTIVE TENSOR
C
C     ICOCO=-1   EVALUATES ELASTIC CONSTITUTIVE TENSOR
C     ICOCO= 0   EVALUATES CONTINUOUS CONSTITUTIVE TENSOR
C     ICOCO= 1   EVALUATES CONSISTENT CONSTITUTIVE TENSOR
C
C***********************************************************************
C
C**** YIELD CRITERION/FLOW POTENTIAL CASES: NCRIT= 31-41
C
C***********************************************************************
C
C     Notes:
C           Due to the fact that the evolution of all the internal
C           variables depends linearly on \lambda, all the internal
C           variables reach their converged values when the condition 
C           \Delta\lambda=0 is verified (considering \lambda=0 at the
C           beginning of the integration of the viscoplastic model).
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION DESIG(*), DMATX(NSTRS,*), EBASE(*),
     .          PRESG(*), PROPS(*),       SGTOT(*), 
     .          SIGMA(*), STRAP(*),       DMTEP(*),
     .          DSTRA(*), STRAN(*),       DBASE(*), FATIG(*)
C
      DIMENSION VECTFF(6),   VECTFA(6),  VECTGF(6),   VECTGA(6),
     .          DEVIAF(6),   DEVIAA(6),
     .          DP(6),       DS(6),      DEEPP(6),
     .          DGVECF(6),   DGVECA(6),  VCONVF(200), VCONVA(200),
     .          DMATE(6,6),
     .          DM(6,6),     DMATP(6,6), HARDEF(6),   HARDEA(6),
     .          BACKA(6),    BACKS(6),   SPTOT(6)
      DIMENSION HTRIXF(6,6), HTRIXA(6,6)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION TRMTX(3,3)          ! not transferred yet!!
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
C**** DECIDE IF SUBINCREMENTATION IS NECESSARY - not implemented yet!!
C


C
C**** DEFINES MIXED VARIABLES
C
      ROTESF=ROTES
      ROTESA=1.0D0-ROTES
C
C**** INITIALIZES VISCOPLASTIC PARAMETERS
C
      TLAMDF=0.0D+00
      TLAMDA=0.0D+00
C
C**** VISCOPLASTIC MODULE: FIRST PASS
C
      KFLUI=0      ! does not compute flux
      INDEX=2      ! computes eff. stress & tot. hardening function
C
      IF(NELASF.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFF,VECTGF,DEVIAF,DEEPP,DGVECF,YIELD,
     .             ALENG,ABETAF,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSF,CCEROF,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFF,CCOEBF,CCOEQF,CKOEFF,CKOEBF,CKOEQF,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOF,EXPONF,HTRIXF,TLAMBF,
     .             TLAMDF,DLAMDF,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEF,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAF,AZABBF,AZABCF,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
      IF(NELASA.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFA,VECTGA,DEVIAA,DEEPP,DGVECA,YIELD,
     .             ALENG,ABETAA,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSA,CCEROA,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFA,CCOEBA,CCOEQA,CKOEFA,CKOEBA,CKOEQA,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOA,EXPONA,HTRIXA,TLAMBA,
     .             TLAMDA,DLAMDA,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEA,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAA,AZABBA,AZABCA,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
C**** ANALYSES THE STATE OF THE POINT (ELASTIC OR VISCOPLASTIC)
C
      IF(NELASF.EQ.0.AND.NELASA.EQ.0)
     . CALL RUNEND('ERROR: NO AUSTENITE & NO FERRITE')
C
      IF(NELASF.EQ.1.AND.NELASA.EQ.0) THEN
       IF(PREYSF.LT.0.0)  ! preys can be zero (i.e., viscoelastic model)
     .  CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LE 0.0')
      ENDIF
      IF(NELASF.EQ.0.AND.NELASA.EQ.1) THEN
       IF(PREYSA.LT.0.0D0) !preys can be zero (i.e., viscoelastic model)
     .  CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LE 0.0')
      ENDIF
      IF(YIELD.LT.0.0D0)
     . CALL RUNEND('ERROR: YIELD FUNCTION LE 0.0       ')
C
      IF(NELASF.EQ.1) ESCURF=YIELD-PREYSF
      IF(NELASA.EQ.1) ESCURA=YIELD-PREYSA
C
      IF(NELASF.EQ.1.AND.NELASA.EQ.0) THEN
       TOLINF=TOPLA*PREYSF
       IF(ESCURF.LT.TOLINF) GO TO 60  ! elastic or unloading case
      ENDIF
      IF(NELASF.EQ.0.AND.NELASA.EQ.1) THEN
       TOLINA=TOPLA*PREYSA
       IF(ESCURA.LT.TOLINA) GO TO 60  ! elastic or unloading case
      ENDIF
      IF(NELASF.EQ.1.AND.NELASA.EQ.1) THEN
       TOLINF=TOPLA*PREYSF
       TOLINA=TOPLA*PREYSA
       IF(ESCURF.LT.TOLINF.AND.ESCURA.LT.TOLINA)
     .  GO TO 60                      ! elastic or unloading case
      ENDIF
C
C**** VISCOPLASTIC CASE               ! loading case
C
      KUNL2=0
C
C**** VISCOPLASTIC MODULE: SECOND PASS
C
      KFLUI=1      ! computes flux
      INDEX=1      ! computes the int. variables & lambda multiplicator
C
      IF(NELASF.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFF,VECTGF,DEVIAF,DEEPP,DGVECF,YIELD,
     .             ALENG,ABETAF,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSF,CCEROF,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFF,CCOEBF,CCOEQF,CKOEFF,CKOEBF,CKOEQF,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOF,EXPONF,HTRIXF,TLAMBF,
     .             TLAMDF,DLAMDF,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEF,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAF,AZABBF,AZABCF,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
      IF(NELASA.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFA,VECTGA,DEVIAA,DEEPP,DGVECA,YIELD,
     .             ALENG,ABETAA,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSA,CCEROA,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFA,CCOEBA,CCOEQA,CKOEFA,CKOEBA,CKOEQA,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOA,EXPONA,HTRIXA,TLAMBA,
     .             TLAMDA,DLAMDA,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEA,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAA,AZABBA,AZABCA,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
      NKONT=0
C
  600 CONTINUE
      NKONT=NKONT+1
C
C**** EVALUATES THE TOTAL VISCOPLASTIC PARAMETERS
C
      IF(NELASF.EQ.1) TLAMDF=TLAMDF+DLAMDF
      IF(NELASA.EQ.1) TLAMDA=TLAMDA+DLAMDA
C
C**** COMPUTATION OF: 1) INCREMENTAL & TOTAL PLASTIC DEFORMATION
C                     2) INCREMENTAL & TOTAL STRESSES
C
C     Note: in the second call to plstre.f, SPTOT is passed as argument
C           instead of STRAP in order to obtain the correct plastic
C           deformation (sum of the perlite & ferrite contributions)
C
      DTIMX=DTIME
      DAMAX=DAMAG
      CALL PLSTRE(VECTGF,STRAP,SPTOT,DEEPP,TLAMDF,DLAMDF,
     .            DMATX,SIGMA,SGTOT,DESIG,PRESG,
     .               DS,   DP,
     .            DTIMX,DAMAX,DAMAY,NAUXI)
c    .            STRAN,
c    .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
c    .            BULKM,DISTM,STRRA,STRR4,
c    .            RCGTT,RCGTI,TRACC,TRACI,ENEYE)
C
      CALL PLSTRE(VECTGA,SPTOT,SPTOT,DEEPP,TLAMDA,DLAMDA,
     .            DMATX,SIGMA,SGTOT,DESIG,PRESG,
     .               DS,   DP,
     .            DTIMX,DAMAX,DAMAY,NAUXI)
c    .            STRAN,
c    .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
c    .            BULKM,DISTM,STRRA,STRR4,
c    .            RCGTT,RCGTI,TRACC,TRACI,ENEYE)
C
C**** VISCOPLASTIC MODULE: THIRD PASS
C
      KFLUJ=INT(PROPS(56))
      KFLUI=KFLUJ  ! computes flux according to kfluj
      INDEX=1      ! computes eff. st., int. var., tot. ha. f. & "abeta"
C
      IF(NELASF.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFF,VECTGF,DEVIAF,DEEPP,DGVECF,YIELD,
     .             ALENG,ABETAF,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSF,CCEROF,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFF,CCOEBF,CCOEQF,CKOEFF,CKOEBF,CKOEQF,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOF,EXPONF,HTRIXF,TLAMBF,
     .             TLAMDF,DLAMDF,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEF,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAF,AZABBF,AZABCF,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
      IF(NELASA.EQ.1)
     . CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTFA,VECTGA,DEVIAA,DEEPP,DGVECA,YIELD,
     .             ALENG,ABETAA,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYSA,CCEROA,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEFA,CCOEBA,CCOEQA,CKOEFA,CKOEBA,CKOEQA,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCOA,EXPONA,HTRIXA,TLAMBA,
     .             TLAMDA,DLAMDA,STRAP,SPTOT,DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEA,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAA,AZABBA,AZABCA,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
C**** VISCOPLASTIC CONVERGENCE VERIFICATION
C 
      ICOPL=0
      IF(NELASF.EQ.1.AND.NELASA.EQ.0) THEN
       ESCURF=DLAMDF/ABETAF            ! analogy with plasticity
C
       VCONVF(NKONT)=ESCURF
       TOLINF=TOPLA*PREYSF
       IF(PREYSF.EQ.0.0D0) TOLINF=TOPLA  ! viscoelastic case
       IF(ESCURF.LT.TOLINF) ICOPL=1    ! the viscopl. alg. has converged
      ENDIF
      IF(NELASF.EQ.0.AND.NELASA.EQ.1) THEN
       ESCURA=DLAMDA/ABETAA              ! analogy with plasticity
C
       VCONVA(NKONT)=ESCURA
       TOLINA=TOPLA*PREYSA
       IF(PREYSA.EQ.0.0D0) TOLINA=TOPLA  ! viscoelastic case
       IF(ESCURA.LT.TOLINA) ICOPL=1    ! the viscopl. alg. has converged
      ENDIF
      IF(NELASF.EQ.1.AND.NELASA.EQ.1) THEN
       ESCURF=DLAMDF/ABETAF              ! analogy with plasticity
       ESCURA=DLAMDA/ABETAA              ! analogy with plasticity
C                                        ! to be revised !!!!!!
       VCONVF(NKONT)=ESCURF
       VCONVA(NKONT)=ESCURA
       TOLINF=TOPLA*PREYSF
       TOLINA=TOPLA*PREYSA
       IF(PREYSF.EQ.0.0D0) TOLINF=TOPLA  ! viscoelastic case
       IF(PREYSA.EQ.0.0D0) TOLINA=TOPLA  ! viscoelastic case
       IF(ESCURF.LT.TOLINF.AND.ESCURA.LT.TOLINA)
     .  ICOPL=1                        ! the viscopl. alg. has converged
      ENDIF
      IF(ICOPL.EQ.1) THEN              ! the viscopl. alg. has converged
C
C**** VISCOPLASTIC MODULE: FOURTH PASS
C
       KFLUI=0      ! does not compute flux
       INDEX=3      ! stores internal variables
C
       IF(NELASF.EQ.1)
     .  CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .              VECTFF,VECTGF,DEVIAF,DEEPP,DGVECF,YIELD,
     .              ALENG,ABETAF,DMATE,COEPD,DM,
     .              DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .              RETEN,RETEP,INDEX,KFLUI,PREYSF,CCEROF,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              PMEAN,CCOEFF,CCOEBF,CCOEQF,CKOEFF,CKOEBF,CKOEQF,
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .              VISCOF,EXPONF,HTRIXF,TLAMBF,
     .              TLAMDF,DLAMDF,STRAP,SPTOT,DP,DAMAG,COUTH,
     .              DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEF,
     .              CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAF,AZABBF,AZABCF,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,DSTRA,NAUXI)
C
       IF(NELASA.EQ.1)
     .  CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .              VECTFA,VECTGA,DEVIAA,DEEPP,DGVECA,YIELD,
     .              ALENG,ABETAA,DMATE,COEPD,DM,
     .              DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .              RETEN,RETEP,INDEX,KFLUI,PREYSA,CCEROA,
     .              A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .              PMEAN,CCOEFA,CCOEBA,CCOEQA,CKOEFA,CKOEBA,CKOEQA,
     .              CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .              CCOB7,CCOB8,CCOB9,CLANK,
     .              CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .              VISCOA,EXPONA,HTRIXA,TLAMBA,
     .              TLAMDA,DLAMDA,STRAP,SPTOT,DP,DAMAG,COUTH,
     .              DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDEA,
     .              CAPAP,EFFPL,BACKA,BACKS,POROS,AZABAA,AZABBA,AZABCA,
     .              XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .              FATIG,TRMTX,DSTRA,NAUXI)
C
C**** MODIFIES ABETA PARAMETER
C
       IF(NELASF.EQ.1) ABETAF=ABETAF*ROTESF
       IF(NELASA.EQ.1) ABETAA=ABETAA*ROTESA
C
C**** EVALUATES THE PSEUDO ELASTO-VISCOPLASTIC CONSTITUTIVE TENSOR
C
       IVIRI=1                     ! better as input
       IF(IVIRI.EQ.1) THEN
        IF(NALGO.GT.1.AND.NALGP.EQ.1) THEN
         IF(NELASF.EQ.1)
     .    CALL PLRIGI(ABETAF,NSTRS,DGVECF,PROPS,DMATX,DMTEP,
     .                EBASE,NCRIT,NCRIP,NTYPE,NSTR1,DMATE,VECTFF,KSYMM,
     .                ICOCO,HARDEF,ALFAP,TLAMD,NAUXI)
C
         IF(NELASA.EQ.1)
     .    CALL PLRIGI(ABETAA,NSTRS,DGVECA,PROPS,DMATX,DMTEP,
     .                EBASE,NCRIT,NCRIP,NTYPE,NSTR1,DMATE,VECTFA,KSYMM,
     .                ICOCO,HARDEA,ALFAP,TLAMD,NAUXI)
        ENDIF         ! nalgo.eq.1
       ENDIF          ! iviri.eq.1
C
       IPLDE=0                   ! better as input
       IF(IPLDE.EQ.1) THEN
        WRITE(LUPRI,901) NKONT,IELEM
  901   FORMAT('NUMBER OF ITERATIONS IN THE PLASTIC ALGORITHM=',I5,2X,
     .  'ELEMENT NUMBER=',I7)
       ENDIF                    ! iplde.eq.1
C
       RETURN
C
      ELSE      ! the viscoplastic algorithm has not converged
C
       IF(NKONT.GE.MKONT) THEN
        IF(NELASF.EQ.1) WRITE(LUPRI,900) (VCONVF(IKONT),IKONT=1,NKONT)
        IF(NELASA.EQ.1) WRITE(LUPRI,900) (VCONVA(IKONT),IKONT=1,NKONT)
        CALL RUNEND('ERROR: NO CONVER. INTEG. CONST. EQ.')
       ENDIF
C
      ENDIF
C
      GO TO 600
C
C**** ELASTIC CASE
C 
   60 KUNL2=1
C 
      RETURN
C
  900 FORMAT(1E14.6)
      END
