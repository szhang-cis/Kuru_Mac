      SUBROUTINE PLMODL(SGTOT,EBASE,DMATX,PROPS,
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
     .                  FATIG,TRMTX,NAUXI)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES: The stress invariant
C                             The flow vector
C                             The internal variables
C                             The ABETA coefficient
C 
C
C     Notes:
C
C     h_Cp: plastic hardening coefficient
C     H_Cp: plastic hardening modulus (in thesis, H_Cp=h_Cp \sigma:R)
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
      DIMENSION SGTOTL(6),      DESIGL(6), SPTOTL(6),
     .          DEEPPL(6),      BACKAM(6), BACKSM(6),
     .          VECTGL(6),      VECTFL(6)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION TRMTX(NDIME,*)
      DIMENSION CLAN1(2),       CLAN2(2),       CLAN3(2),
     .          CLAN4(2),       CLAN5(2),       CLAN6(2)
C
C**** RECOVER INTERNAL VARIABLES & THEIR CONJUGATES
C
      IF(INDEX.EQ.1.OR.INDEX.EQ.2) THEN
       CALL PLRECO(EBASE,CAPAP,EFFPL,BACKA,BACKS,POROS,
     .             PREYA,PROPS,CCERO,EFFP2,SRECR,ERECR,YIELO)
       CALL DARECO(DBASE,DAMAG,TAUAM,TAUAP)
       CALL FARECO(FATIG,
     .             FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
C
C**** COMPUTATION OF INVARIANTS & EFFECTIVE STRESS
C
       CALL PLINVA(NCRIT,PROPS,NTYPE,SGTOT,SGTOTL,ANGFI,CAPAP,
     .             PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .             DFTEQ,DFAFI,YIELD,NSTR1,ALFAT,BETAT,
     .             GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,
     .             CCERO,POROS,IPORO,CPOE1,CPOE2,CLANK,
     .             CLAN1(1),CLAN2(1),CLAN3(1),
     .             CLAN4(1),CLAN5(1),CLAN6(1),
     .             LARGE,IFREN,IKINE,NDIME,XJACM,XJA3M,DETJM,
     .             TRMTX)
       CALL DAINVA(NCRIT,PROPS,NTYPE,SGTOTL,ANGFI,CAPAP,
     .             PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .             DFTEQ,DFAFI,YIELD,NSTR1,ALFAT,BETAT,
     .             GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,CCERO,
     .             POROS,IDAMG,DAMAG,TAUAM,TAUAP,
     .             CCEROM,CCEROP,GAMMAM,GAMMAP,LARGE)
       CALL DASTOR(DBASE,DAMAG,TAUAM,TAUAP,INDEX)
      ENDIF         ! index.eq.1.or.index.eq.2
C
C**** COMPUTATION OF PLASTIC FLUX & NORMAL VECTORS
C
C     Note:
C
C     If NCRIT=NCRIP, that means that VECTF is proportional to VECTG,
C     but not necessarily equal => it is necessary to compute both
C
      IF(INDEX.EQ.1.OR.INDEX.EQ.2) THEN
       IF(KFLUJ.EQ.1) THEN
        NC42T=INT(PROPS(45))
        CALL PLYIEL(VECTFL,ANGFI,NCRIT,NC42T,RETEN,ALFAT,
     .              GAMAT,PROPS(48),PMEAN,VARJ2,VARJ3,DEVIA,
     .              STEFF,SMAXS,THETA,PROPS,A1COE,A2COE,A3COE,
     .              CCERO,CAPAP,YIELD,
     .              POROS,CPOE1,CPOE2,CLANK,
     .              CLAN1(1),CLAN2(1),CLAN3(1),
     .              CLAN4(1),CLAN5(1),CLAN6(1))            ! NORMAL TO F
C
        NC42P=INT(PROPS(46))
        CALL PLYIEL(VECTGL,ANGSI,NCRIP,NC42P,RETEP,PROPS(50),
     .              PROPS(51),PROPS(49),PMEAN,VARJ2,VARJ3,DEVIA,
     .              STEFF,SMAXS,THETA,PROPS,A1COP,A2COP,A3COP,
     .              CCERO,CAPAP,YIELD,
     .              POROS,CPOE1,CPOE2,CLANK,
     .              CLAN1(2),CLAN2(2),CLAN3(2),
     .              CLAN4(2),CLAN5(2),CLAN6(2))            ! NORMAL TO G
       ENDIF        ! kfluj.eq.1
      ENDIF         ! index.eq.1.or.index.eq.2
C
C**** COMPUTATION OF INTERNAL VARIABLES & THEIR CONJUGATES
C
C     Note: DMATP,DMATX,DMATE,DP not used
C
      IF(INDEX.EQ.1) THEN
       YIELP=YIELD
       IF(NASNA.EQ.1)                              ! non-associate model
     . CALL PLINVA(NCRIP,PROPS,NTYPE,SGTOT,SGTOTL,ANGFI,CAPAP,
     .             PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .             DFTEQ,DFAFI,YIELP,NSTR1,ALFAT,BETAT,
     .             GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,
     .             CCERO,POROS,IPORO,CPOE1,CPOE2,CLANK,
     .             CLAN1(2),CLAN2(2),CLAN3(2),
     .             CLAN4(2),CLAN5(2),CLAN6(2),
     .             LARGE,IFREN,IKINE,NDIME,XJACM,XJA3M,DETJM,
     .             TRMTX)
       CALL PLVARI(PROPS,ALENG,SGTOTL,DEEPP,DEEPPL,CAPAP,
     .             VECTGL,NSTRS,ANGFI,ANGSI,
     .             DFTEQ,DFAFI,YIELD,CAPAD,NSTR1,YIELP,
     .             COEPD,VADEG,DMATP,COUTH,
     .             DMATX,DMATE,RETEN,DM,CCERO,NCRIP,EFFPL,EFFP1,DP,
     .             CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,DISTM,
     .             BACKA,BACKAM,BACKS,BACKSM,SPTOT,SPTOTL,ISOTT,IKINE,
     .             POROS,CPOEF,CPOE1,CPOE2,IPORO,DESIG,DESIGL,
     .             A2COE,A3COE,PREYA,DTIME,DAMAG,EFFP2,
     .             LARGE,IFREN,NDIME,XJACM,XJA3M,XJACI,XJA3I,IREKI,
     .             IRECR,RECRY,SRECR,ERECR,ASREC,ALREC,TRMTX,TLAMD,
     .             NAUXI)
       CALL DAVARI(PROPS,ALENG,SGTOTL,DEEPPL,CAPAP,
     .             VECTGL,NSTRS,ANGFI,ANGSI,
     .             DFTEQ,DFAFI,YIELD,CAPAD,NSTR1,
     .             COEPD,VADEG,DMATP,COUTH,
     .             DMATX,DMATE,RETEN,DM,CCERO,NCRIP,EFFPL,EFFP1,DP,
     .             CCOEF,CCOEB,CCOEQ,CKOEF,CKOEB,CKOEQ,
     .             BACKA,BACKS,SPTOTL,ISOTT,IKINE,NTYPE,
     .             DTIME,TLAMD,DLAMD,DAMAG,CDAMA,CFRAC,CCOB9,
     .             CCOB10,CCOB11,CCOB12,IDAMG,PMEAN,NAUXI)
      ENDIF         ! index.eq.1
C
C**** COMPUTATION OF INVARIANTS, EFF. STRESS, PL. FLUX & NORMAL VECTORS
C     (NECESSARY WHEN KINEMATIC HARDENING IS USED)
C
      IF(INDEX.EQ.1) THEN
       IF(IREKI.EQ.1) THEN
        CALL PLINVA(NCRIT,PROPS,NTYPE,SGTOT,SGTOTL,ANGFI,CAPAP,
     .              PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .              DFTEQ,DFAFI,YIELD,NSTR1,ALFAT,BETAT,
     .              GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,
     .              CCERO,POROS,IPORO,CPOE1,CPOE2,CLANK,
     .              CLAN1(1),CLAN2(1),CLAN3(1),
     .              CLAN4(1),CLAN5(1),CLAN6(1),
     .              LARGE,IFREN,IKINE,NDIME,XJACM,XJA3M,DETJM,
     .              TRMTX)
        CALL DAINVA(NCRIT,PROPS,NTYPE,SGTOTL,ANGFI,CAPAP,
     .              PMEAN,VARJ2,VARJ3,DEVIA,STEFF,SMAXS,THETA,
     .              DFTEQ,DFAFI,YIELD,NSTR1,ALFAT,BETAT,
     .              GAMAT,RETEN,A1COE,A2COE,A3COE,BACKS,CCERO,
     .              POROS,IDAMG,DAMAG,TAUAM,TAUAP,
     .              CCEROM,CCEROP,GAMMAM,GAMMAP,LARGE)
        CALL DASTOR(DBASE,DAMAG,TAUAM,TAUAP,INDEX)
C
        IF(KFLUJ.EQ.1) THEN
         NC42T=INT(PROPS(45))
         CALL PLYIEL(VECTFL,ANGFI,NCRIT,NC42T,RETEN,ALFAT,
     .               GAMAT,PROPS(48),PMEAN,VARJ2,VARJ3,DEVIA,
     .               STEFF,SMAXS,THETA,PROPS,A1COE,A2COE,A3COE,
     .               CCERO,CAPAP,YIELD,
     .               POROS,CPOE1,CPOE2,CLANK,
     .               CLAN1(1),CLAN2(1),CLAN3(1),
     .               CLAN4(1),CLAN5(1),CLAN6(1))           ! NORMAL TO F
C
         NC42P=INT(PROPS(46))
         CALL PLYIEL(VECTGL,ANGSI,NCRIP,NC42P,RETEP,PROPS(50),
     .               PROPS(51),PROPS(49),PMEAN,VARJ2,VARJ3,DEVIA,
     .               STEFF,SMAXS,THETA,PROPS,A1COP,A2COP,A3COP,
     .               CCERO,CAPAP,YIELD,
     .               POROS,CPOE1,CPOE2,CLANK,
     .               CLAN1(2),CLAN2(2),CLAN3(2),
     .               CLAN4(2),CLAN5(2),CLAN6(2))           ! NORMAL TO G
        ENDIF       ! kfluj.eq.1
       ENDIF        ! ireki.eq.1
      ENDIF         ! index.eq.1
C
C**** COMPUTATION OF TOTAL HARDENING FUNCTION (PREYS)
C
C     Note: SGTOT not used
C
      IF(INDEX.EQ.1.OR.INDEX.EQ.2) THEN
       CALL PLTOHA(NCRIT,PROPS,NTYPE,SGTOT,ANGFI,CAPAP,PREYS,DAMAG,
     .             CCERO,A2COE,A3COE,ISOTT)
       CALL FATOHA(PREYS,YIELD,CCERO,
     .             FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
      ENDIF         ! index.eq.1.or.index.eq.2
C
C**** COMPUTATION OF "HARDS" (TERM OF THE "ABETA" MULTIPLICATOR)
C
C     Note: the case NCRIT=41 for LARGE > 0 has to be revised  !!!
C
      IF(INDEX.EQ.1) THEN
       CALL PLHARD(HARDS,SGTOTL,VECTFL,VECTGL,DMATP,DESIGL,HARDE,BACKS,
     .             NSTRS,NSTR1,NCRIT,
     .             ISOTT,IKINE,IPORO,
     .             CCOEF,A2COE,CCOEB,CCOEQ,CPOEF,CPOE1,CPOE2,
     .             CKOEF,CKOEB,CKOEQ,
     .             EFFPL,EFFP1,EFFP2,YIELD,PREYS,PREYA,PMEAN,POROS,
     .             CCERO,DAMAG)
      ENDIF         ! index.eq.1
C
C**** COMPUTATION OF THE "ABETA" MULTIPLICATOR
C
      IF(INDEX.EQ.1) THEN
       CALL PLFLOW(VECTF,VECTFL,ABETA,DMATP,HARDS,NSTRS,NTYPE,
     .             PROPS,VECTG,VECTGL,DGVEC,NSTR1,
     .             LARGE,IFREN,NDIME,XJACM,XJA3M,TRMTX)
      ENDIF         ! index.eq.1
C
C**** STORE INTERNAL VARIABLES (& THEIR CONJUGATES ONLY FOR POSTPROCESS)
C     AND PLASTIC DEFORMATION
C
      IF(INDEX.EQ.3) THEN
       CALL PLSTOR(EBASE,CAPAP,EFFPL,BACKA,BACKS,POROS,
     .             PREYS,STRAP,SPTOT,EFFP2,SRECR,ERECR,YIELD)
       CALL DASTOR(DBASE,DAMAG,TAUAM,TAUAP,INDEX)
       CALL FASTOR(FATIG,
     .             FREDU,ANUMF,STENI,STEN1,STEN2,STMIN,STMAX,REVEF)
      ENDIF
C
      RETURN
      END
