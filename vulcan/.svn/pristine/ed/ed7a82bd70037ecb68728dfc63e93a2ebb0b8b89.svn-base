      SUBROUTINE PLASTC(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
     .                  PROPS,SGTOT,DMTEP,DMATE,KUNL2,ALENG,
     .                  ENELI,DELMU,DM,DSTRA,ANGFI,ANGSI,
     .                  DMATP,STRAN,ALFAT,BETAT,BULKM,DISTM,
     .                  GAMAT,RETEN,RETEP,CCERO,CCEROM,CCEROP,
     .                  A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                  DEEPP,DTEMG,DCCER,
     .                  CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .                  CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .                  CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .                  CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .                              DBASE,CDAMA,CFRAC,COUTH,
     .                  CPOEF,CPOE1,CPOE2,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  FATIG,TRMTX,NAUXI)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR THE PLASTIC MODEL
C
C***********************************************************************
C
C**** MORE SPECIFICALLY, ROUTINE plastc.f DOES:
C
C     1) ENTERS PLASTIC MODULE plmodl.f
C     2) EVALUATES THE ELASTO-PLASTIC CONSTITUTIVE TENSOR
C
C***********************************************************************
C
C**** THE PLASTIC MODULE plmodl.f COMPUTES: 
C
C     1) EFFECTIVE STRESS (INVARIANTS)
C     2) PLASTIC FLUX
C     3) INTERNAL VARIABLES (i.e., INTERNAL VARIABLES WITHOUT PLASTIC
C        DEFORMATION IN WHAT FOLLOWS)
C     4) TOTAL HARDENING FUNCTION
C     5) DENOMINATOR OF RATE PLASTIC CONSISTENCY PARAMETER (1/"ABETA")
C
C
C     FIRST PASS:   ANALYSES THE STATE OF THE POINT BY MEANS OF THE
C                   YIELD FUNCTION (ELASTIC OR PLASTIC) 
C     SECOND PASS:  COMPUTES FLUX, INTERNAL VARIABLES AND MULTIPLICATOR
C                   "ABETA" (THE LAST RELATED WITH THE PLASTIC STRAINS)
C     THIRD PASS:  -IMPOSES THE KUHN-TUCKER & CONSISTENCY (PRAGER)
C                   CONDITIONS
C                  -VERIFICATES THE PLASTIC DISSIPATION
C                  -VERIFICATES THE "MAXIMUM PLASTIC DISSIPATION
C                   POSTULATE" (FOR COUPLED & UNCOUPLED PROBLEMS)
C                  -VERIFICATES ADDITIONAL PLASTIC CONDITIONS
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
C                VARIABLES, TOTAL HARDENING FUNCTION & MULTIPLICATOR
C                "ABETA"
C     INDEX=2  - COMPUTES EFFECTIVE STRESS (INVARIANTS) & TOTAL 
C                HARDENING FUNCTION
C     INDEX=3  - STORES INTERNAL VARIABLES
C
C     ICONPL=0 - DOES NOT PERFORM SOME PLASTIC CONTROLS
C     ICONPL=1 - PERFORMS SOME PLASTIC CONTROLS
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
C           variables depends linearly on \dot \lambda, all the internal
C           variables reach their converged values when the condition 
C           \Delta(\dot \lambda)=0 is verified. The condition
C           \Delta(\dot \lambda)=0 is equivalent to the condition F=0
C           (numerically, \Delta(\dot \lambda) is proportional to
C           \DeltaF or directly F due to the fact that F=0 for the last
C           converged step).
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
      DIMENSION DESIG(*), DMATX(NAUXI,*), EBASE(*),
     .          PRESG(*), PROPS(*),       SGTOT(*), 
     .          SIGMA(*), STRAP(*)  ,     DMTEP(*),
     .          DSTRA(*), STRAN(*),       DBASE(*), FATIG(*)
C
      DIMENSION VECTF(6), VECTG(6),       DEVIA(6),
     .          DP(6),    DS(6),          DEEPP(6), DEEPT(6),
     .          DGVEC(6), VCONV(200),     DMATE(6,6),
     .          DM(6,6),  DMATP(6,6),     HARDE(6),
     .          BACKA(6), BACKS(6),       SPTOT(6),       DESGT(6)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION TRMTX(NDIME,*)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
C**** DECIDE IF SUBINCREMENTATION IS NECESSARY
C
      CALL PLSUBI(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,             ! up to here idem to plmodl.f
     .            SIGMA,DESGT,NSUBP,NAUXI)
C
      IF(NSUBP.GT.1) THEN         ! to compute the plastic coupling term
       DO ISTR1=1,NSTR1
        DEEPT(ISTR1)=0.0D0
       ENDDO
      ENDIF
      COUTH=0.0D0                 ! only useful for ITERMP > 0
C
C**** SUBINCREMENTS LOOP
C
      DO 1600 ISUBP=1,NSUBP
C
      IF(NSUBP.GT.1) THEN
       DAMAX=DAMAG
       DO ISTRS=1,NSTRS                        ! redefines SGTOT & SIGMA
        SGTOT(ISTRS)=SGTOT(ISTRS)+DESGT(ISTRS)/NSUBP
        SIGMA(ISTRS)=SGTOT(ISTRS)/(1.0D0-DAMAX)
       ENDDO
      ENDIF                             ! nsubp.gt.1
C
C**** INITIALIZES VARIABLES RELATED TO PLASTIC DEFORMATION
C
      DO ISTR1=1,NSTR1
       DP(ISTR1)=0.0D0
       DEEPP(ISTR1)=0.0D0
      ENDDO
C
C**** INITIALIZES RATE PLASTIC CONSISTENCY PARAMETER
C
      TLAMD=0.0D+00
C
C**** PLASTIC MODULE: FIRST PASS
C
      KFLUI=0      ! does not compute flux
      INDEX=2      ! computes eff. stress & tot. hardening function
C
      CALL PLMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,NAUXI)
C
C**** ANALYSES THE STATE OF THE POINT (ELASTIC OR PLASTIC)
C
      IF(PREYS.LE.0.0D0)
     . CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LE 0.0')
      IF(YIELD.LT.0.0D0)
     . CALL RUNEND('ERROR: YIELD FUNCTION LT 0.0       ')
      ESCUR=YIELD-PREYS
C
      TOLIN=TOPLA*PREYS
      IF(ESCUR.LT.TOLIN) GO TO 60            ! elastic or unloading case
C
C**** PLASTIC CASE                           ! loading case
C
      KUNL2=0
C
C**** PLASTIC MODULE: SECOND PASS
C
      KFLUI=1      ! computes flux
      INDEX=1      ! computes the int. variables & "abeta" multiplicator
C
      CALL PLMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,NAUXI)
C
      NKONT=0
C
  600 CONTINUE
      NKONT=NKONT+1
C
C**** EVALUATES THE INCREMENTAL RATE PLASTIC CONSISTENCY PARAMETER
C
      IF(ABETA.LE.0.0D0)
     . CALL RUNEND('ERROR: ABETA MULTIPLICATOR LE 0    ')
      DLAMD=ESCUR*ABETA
C
C**** EVALUATES THE TOTAL RATE PLASTIC CONSISTENCY PARAMETER
C
      TLAMD=TLAMD+DLAMD
C
C**** COMPUTATION OF: 1) INCREMENTAL & TOTAL PLASTIC DEFORMATION
C                     2) INCREMENTAL & TOTAL STRESSES
C
      DTIMX=1.0D0
      DAMAX=DAMAG
      CALL PLSTRE(VECTG,STRAP,SPTOT,DEEPP,TLAMD,DLAMD,
     .            DMATX,SIGMA,SGTOT,DESIG,PRESG,
     .               DS,   DP,
     .            DTIMX,DAMAX,DAMAY,NAUXI)
c    .            STRAN,
c    .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
c    .            BULKM,DISTM,STRRA,STRR4,
c    .            RCGTT,RCGTI,TRACC,TRACI,ENEYE)
C
C**** PLASTIC MODULE: THIRD PASS
C
      KFLUJ=INT(PROPS(56))
      KFLUI=KFLUJ  ! computes flux according to kfluj
      INDEX=1      ! computes eff. st., int. var., tot. ha. f. & "abeta"
C
      CALL PLMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,NAUXI)
C
C**** ITERATIVE PLASTIC CONTROLS
C
      ICONPL=0     ! better as input !!
      IF(ICONPL.EQ.1)
     . CALL CONPLA(SGTOT,DEVIA,VECTF,DESIG,DMATX,ABETA,DTEMG,
     .             DCCER,NAUXI)
C
C**** PLASTIC CONVERGENCE VERIFICATION
C 
      IF(PREYS.LE.0.0D0)
     . CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LE 0.0')
      IF(YIELD.LT.0.0D0)
     . CALL RUNEND('ERROR: YIELD FUNCTION LT 0.0       ')
      ESCUR=YIELD-PREYS
C
      VCONV(NKONT)=ESCUR
      TOLIN=TOPLA*PREYS
      IF(DABS(ESCUR).LT.TOLIN.AND.
     .   DABS(DLAMD).LT.TOPLA) THEN     ! the plastic alg. has converged
C
C**** PLASTIC MODULE: FOURTH PASS
C
       KFLUI=0      ! does not compute flux
       INDEX=3      ! stores internal variables
C
       CALL PLMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .             ALENG,ABETA,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .             CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .             CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,NAUXI)
C
C**** CONVERGED PLASTIC CONTROLS
C
       ICONPL=0     ! better as input !!
       IF(ICONPL.EQ.1)
     .  CALL CONPLA(SGTOT,DEVIA,VECTF,DESIG,DMATX,ABETA,DTEMG,
     .              DCCER,NAUXI)
C
C**** EVALUATES ELASTO-PLASTIC CONSTITUTIVE TENSOR
C
       IF(NALGO.GT.1.AND.NALGP.EQ.1)
     .  CALL PLRIGI(ABETA,NSTRS,DGVEC,PROPS,DMATX,DMTEP,
     .              EBASE,NCRIT,NCRIP,NTYPE,NSTR1,DMATE,VECTF,KSYMM,
     .              ICOCO,HARDE,ALFAP,TLAMD,NAUXI)
C
       IPLDE=0                    ! better as input !!
       IF(IPLDE.EQ.1) THEN
        WRITE(LUPRI,901) NKONT,IELEM
  901   FORMAT('NUMBER OF ITERATIONS IN THE PLASTIC ALGORITHM=',I5,2X,
     .  'ELEMENT NUMBER=',I7)
       ENDIF                      ! iplde.eq.1
C
       IF(NSUBP.GT.1) THEN        ! to compute the plastic coupling term
        DO ISTR1=1,NSTR1
         DEEPT(ISTR1)=DEEPT(ISTR1)+DEEPP(ISTR1)
         DEEPP(ISTR1)=DEEPT(ISTR1)
        ENDDO
       ENDIF
C
       GO TO 1600
C
      ELSE      ! the plastic algorithm has not converged
C
       IF(NKONT.GE.MKONT) THEN
        WRITE(LUPRI,900) (VCONV(IKONT),IKONT=1,NKONT)
        CALL RUNEND('ERROR: NO CONVER. INTEG. CONST. EQ.')
       ENDIF
       GO TO 600
C
      ENDIF                     ! escur.lt.tolin.and....
C
C**** ELASTIC CASE
C
   60 KUNL2=1
C
 1600 CONTINUE          ! isubp=1,nsubp
C
#endif
      RETURN
C
  900 FORMAT(1E14.6)
      END
