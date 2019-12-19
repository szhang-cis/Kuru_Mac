      SUBROUTINE PLSUBI(SGTOT,EBASE,DMATX,PROPS,
     .                  VECTF,VECTG,DEVIA,DEEPP,DGVEC,YIELD,
     .                  ALENG,ABETA,DMATE,COEPD,DM,
     .                  DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,BULKM,DISTM,
     .                  RETEN,RETEP,INDEX,KFLUJ,PREYS,CCERO,
     .                  CCEROM,CCEROP,DBASE,CDAMA,CFRAC,
     .                  A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                  PMEAN,CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .                  CCOB1,CCOB2,CCOB3,CCOB4,CCOB5,CCOB6,
     .                  CCOB7,CCOB8,CCOB9,CCOB10,CCOB11,CCOB12,CLANK,
     .                  CLAN1,CLAN2,CLAN3,CLAN4,CLAN5,CLAN6,
     .                  TLAMD,DLAMD,STRAP,SPTOT,   DP,DAMAG,COUTH,
     .                  CPOEF,CPOE1,CPOE2,DESIG,HARDE,
     .                  CAPAP,EFFPL,BACKA,BACKS,POROS,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  FATIG,TRMTX,          ! up to here idem plmodl.f
     .                  SIGMA,DESGT,NSUBP,NAUXI)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE NECESSARY NUMBER OF SUBINCREMENTS FOR
C     THE PLASTIC ALGORITHM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION DMATX(NAUXI,*), EBASE(*),
     .          PROPS(*),       SGTOT(*),
     .          DBASE(*),       FATIG(*)
C
      DIMENSION VECTF(6),       VECTG(6),
     .          DEVIA(6),       DEEPP(6),
     .          DGVEC(6)
      DIMENSION DMATE(6,6),     DMATP(6,6),
     .             DM(6,6),     DP(6)
      DIMENSION BACKA(6),       BACKS(6),  STRAP(*), SPTOT(*), DESIG(*),
     .          HARDE(6)
      DIMENSION SIGMA(*),       DESGT(*)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION TRMTX(NDIME,*)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
      NSUBP=1
      IF(EXCTP.EQ.0.0D0) RETURN
C
C**** PLASTIC MODULE: ZERO PASS (required for subincrementation)
C
C     Note: this ZERO PASS is similar to the FIRST PASS of plmodl.f
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
     . CALL RUNEND('ERROR: YIELD FUNCTION LE 0.0       ')
      ESCUR=YIELD-PREYS
C
C**** DETERMINES NUMBER OF SUBINCREMENTS OF THE PLASTIC ALGORITHM
C
      IF(PREYS.GT.0.0D0) THEN
       NSUBP=INT(ESCUR/(EXCTP*PREYS))
      ELSE
       NSUBP=INT(ESCUR/EXCTP)                       ! viscoelastic model
      ENDIF
      IF(NSUBP.LT.1)     NSUBP=1
      IF(NSUBP.GT.MSUBP) NSUBP=MSUBP
C
      IF(NSUBP.GT.1) THEN
       DAMAX=DAMAG
       DO ISTRS=1,NSTRS            ! load DESIG, redefines SGTOT & SIGMA
        DESGT(ISTRS)=DESIG(ISTRS)
        SGTOT(ISTRS)=SGTOT(ISTRS)-DESGT(ISTRS)
        SIGMA(ISTRS)=SGTOT(ISTRS)/(1.0D0-DAMAX)
       ENDDO
      ENDIF                            ! nsubp.gt.1
C
      RETURN
      END
