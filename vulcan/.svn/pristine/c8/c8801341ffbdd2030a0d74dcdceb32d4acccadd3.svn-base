      SUBROUTINE VISCOC(SIGMA,DESIG,STRAP,PRESG,EBASE,DMATX,
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
     .                  VISCO,EXPON,DBASE,CDAMA,CFRAC,COUTH,
     .                  CPOEF,CPOE1,CPOE2,
     .                  AZABA,AZABB,AZABC,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  FATIG,TRMTX,NAUXI,YIELD)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES FOR THE VISCOPLASTIC MODEL
C
C***********************************************************************
C
C**** MORE SPECIFICALLY, ROUTINE viscoc.f DOES:
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
C
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
      DIMENSION DESIG(*), DMATX(NAUXI,*), EBASE(*),
     .          PRESG(*), PROPS(*),       SGTOT(*), 
     .          SIGMA(*), STRAP(*),       DMTEP(*),
     .          DSTRA(*), STRAN(*),       DBASE(*), FATIG(*)
C
      DIMENSION VECTF(6), VECTG(6),       DEVIA(6),
     .          DP(6),    DS(6),          DEEPP(6), DEEPT(6),
     .          DGVEC(6), VCONV(200),     DMATE(6,6),
     .          DM(6,6),  DMATP(6,6),     HARDE(6),
     .          BACKA(6), BACKS(6),       SPTOT(6),       DESGT(6)
      DIMENSION HTRIX(6,6)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION TRMTX(NDIME,*)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
C**** DECIDE IF SUBINCREMENTATION IS NECESSARY
C
      CALL VISUBI(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            VISCO,EXPON,HTRIX,TLAMB,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,AZABA,AZABB,AZABC,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,             ! up to here idem to vimodl.f
     .            SIGMA,DESGT,NSUBP,DSTRA,NAUXI)
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
C**** INITIALIZES VISCOPLASTIC PARAMETER
C
      TLAMD=0.0D+00
C
C**** VISCOPLASTIC MODULE: FIRST PASS
C
      KFLUI=0      ! does not compute flux
      INDEX=2      ! computes eff. stress & tot. hardening function
C
      CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            VISCO,EXPON,HTRIX,TLAMB,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,AZABA,AZABB,AZABC,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,DSTRA,NAUXI)
C
C**** ANALYSES THE STATE OF THE POINT (ELASTIC OR VISCOPLASTIC)
C
      IF(PREYS.LT.0.0D0)  ! preys can be zero (i.e., viscoelastic model)
     . CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LT 0.0')
      IF(YIELD.LT.0.0D0)
     . CALL RUNEND('ERROR: YIELD FUNCTION LE 0.0       ')
      ESCUR=YIELD-PREYS
C
      TOLIN=TOPLA*PREYS
      IF(ESCUR.LT.TOLIN) GO TO 60            ! elastic or unloading case
C
C**** VISCOPLASTIC CASE                      ! loading case
C
      KUNL2=0
C
C**** VISCOPLASTIC MODULE: SECOND PASS
C
      KFLUI=1      ! computes flux
      INDEX=1      ! computes the int. variables & lambda multiplicator
C
      CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            VISCO,EXPON,HTRIX,TLAMB,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,AZABA,AZABB,AZABC,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,DSTRA,NAUXI)
C
      NKONT=0
C
  600 CONTINUE
      NKONT=NKONT+1
C
C**** EVALUATES THE TOTAL VISCOPLASTIC PARAMETER
C
      TLAMD=TLAMD+DLAMD
C
C**** COMPUTATION OF: 1) INCREMENTAL & TOTAL PLASTIC DEFORMATION
C                     2) INCREMENTAL & TOTAL STRESSES
C
      DTIMX=DTIME
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
C**** VISCOPLASTIC MODULE: THIRD PASS
C
      KFLUJ=INT(PROPS(56))
      KFLUI=KFLUJ  ! computes flux according to kfluj
      INDEX=1      ! computes eff. st., int. var., tot. ha. f. & "abeta"
C
      CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .            VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .            ALENG,ABETA,DMATE,COEPD,DM,
     .            DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .            RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .            CCEROM,CCEROP,
     .            A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .            PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .            CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .            CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .            CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .            VISCO,EXPON,HTRIX,TLAMB,
     .            TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .            DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .            CAPAP,EFFPL,BACKA,BACKS,POROS,AZABA,AZABB,AZABC,
     .            XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .            FATIG,TRMTX,DSTRA,NAUXI)
C
C**** ITERATIVE VISCOPLASTIC CONTROLS
C
c     ICONPL=0     ! better as input !!
c     IF(ICONPL.EQ.1)
c    . CALL CONPLA(SGTOT,DEVIA,VECTF,DESIG,DMATX,ABETA,DTEMG,
c    .             DCCER,NAUXI)
C
C**** VISCOPLASTIC CONVERGENCE VERIFICATION
C
      IF(PREYS.LT.0.0D0)
     . CALL RUNEND('ERROR: TOT. HARDENING FUNCT. LT 0.0')
      IF(YIELD.LT.0.0D0)
     . CALL RUNEND('ERROR: YIELD FUNCTION LT 0.0       ')
      ESCUR=TLAMB                      ! analogy with plasticity
C
      VCONV(NKONT)=ESCUR
      TOLIN=TOPLA*PREYS
      IF(PREYS.EQ.0.0D0) TOLIN=TOPLA   ! viscoelastic case
      IF(DABS(ESCUR).LT.TOLIN.AND.
     .   DABS(DLAMD).LT.TOPLA) THEN    ! the viscopl. alg. has converged
C
C**** VISCOPLASTIC MODULE: FOURTH PASS
C
       KFLUI=0      ! does not compute flux
       INDEX=3      ! stores internal variables
C
       CALL VIMODL(SGTOT,EBASE,DMATX,PROPS,
     .             VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .             ALENG,ABETA,DMATE,COEPD,DM,
     .             DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .             RETEN,RETEP,INDEX,KFLUI,PREYS,CCERO,
     .             CCEROM,CCEROP,
     .             A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .             PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .             CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .             CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .             CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .             VISCO,EXPON,HTRIX,TLAMB,
     .             TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .             DBASE,CDAMA,CFRAC,CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .             CAPAP,EFFPL,BACKA,BACKS,POROS,AZABA,AZABB,AZABC,
     .             XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             FATIG,TRMTX,DSTRA,NAUXI)
C
C**** CONVERGED PLASTIC CONTROLS
C
c      ICONPL=0     ! better as input !!
c      IF(ICONPL.EQ.1)
c    .  CALL CONPLA(SGTOT,DEVIA,VECTF,DESIG,DMATX,ABETA,DTEMG,
c    .              DCCER,NAUXI)
C
C**** EVALUATES THE PSEUDO ELASTO-VISCOPLASTIC CONSTITUTIVE TENSOR
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
      ELSE      ! the viscoplastic algorithm has not converged
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
