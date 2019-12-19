      SUBROUTINE SOLIDI3(SIGMA,PRESG,DMATX,PROPS,SGTOT,SIGM0,DMATE,
     .                   DMTEP,DMATP,ANGFI,ANGSI,ALFAT,BETAT,GAMAT,
     .                   KUNL4,RETEN,RETEP,TEMPG,CCERO,DESIG,
     .                   CCEROM,CCEROP,BULKM,DISTM,
     .                   DSTRA,STRAN,DTEMG,
     .                   A1COE,A2COE,A3COE,A1COP,A2COP,A3COP,
     .                   STRAP,DMAPL,DEETH,COUTD,DCCER,TEINI,FPCHL,
     .                   SHAPE,DBASE,SHRIN,
     .                   XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                   GPCOD,VANIS,STRS0,TRMTX,EBASE,
     .                   CCEROF,CCEROA,
     .                   DMATXF,DMATXA,DMATEF,DMATEA,DMATPF,DMATPA,
     .                   STRANF,STRANA,SIGMAF,SIGMAA,
     .                   STRAPF,STRAPA,
     .                   DMTEPF,DMTEPA,
     .                   SGTOTF,SGTOTA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TOTAL STRESSES & THE CONSTITUTIVE TENSOR
C     FOR THE DUAL PHASE STEEL MODEL
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
      DIMENSION DMATXF(6,6),    DMATXA(6,6)
      DIMENSION DMATEF(6,6),    DMATEA(6,6)
      DIMENSION DMATPF(6,6),    DMATPA(6,6)
      DIMENSION DMAT1F(6,6),    DMAT1A(6,6)
      DIMENSION STRANF(6),      STRANA(6)
      DIMENSION SIGMAF(6),      SIGMAA(6)
      DIMENSION STRAPF(6),      STRAPA(6)
      DIMENSION DMTEPF(*),      DMTEPA(*)
      DIMENSION SGTOTF(*),      SGTOTA(*)
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
      IF(ITERME.GE.0) THEN                  ! coupled problems
       call runend('ERROR: coupled problem not implemented')
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
       ALPHI=0.0D0
       ALPHO=0.0D0
       ALPHA=0.0D0
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
C
C**** SECANT FERRITE CONSTITUTIVE TENSOR TO INTEGRATE THE STRESSES
C
      CALL CALCCT(DMATXF,YOUNGF,POISMF,ROTET,    6,
     .            POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .            BULKM,DISTM,DAMAP1,DAMAP2,PROPS,
     .            VANIS,ROMTX,TRMTX,
     .            YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IF(ISOTR.NE.0) THEN
       DO IDIME=1,NDIME     ! transpose of TRMTX
        DO JDIME=1,NDIME
         TRMTT(IDIME,JDIME)=TRMTX(JDIME,IDIME)
        ENDDO
       ENDDO
       CALL FOURTH(DMATXF,XJACM,TRMTT,XJA3M,1.0D0,1.0D0,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAPF,RCGPI,FINET,TRABE,DCINV,EBASE)
      ENDIF
C
C**** COMPUTES THE PRODUCT (BULK * THERMAL DEFORMATION)
C
      BULKA=BULKM*STRRA
C
C**** TRANSFORMS CONSTITUTIVE TENSOR FOR LARGE STRAIN ANALYSIS
C
      IF(LARGE.NE.0)
     . CALL FOURTH(DMATXF,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAPF,RCGPI,FINET,TRABE,DCINV,EBASE)
C
C**** STORE SOME CONSTITUTIVE TENSORS
C
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        DMATEF(ISTRS,JSTRS)=DMATXF(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMATPF(ISTRS,JSTRS)=DMATXF(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMAT1F(ISTRS,JSTRS)=DMATXF(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** SECANT CONSTITUTIVE TENSOR TO BE USED IN STIFFNESS MATRIX
C
      IF(LIQUI.EQ.1) THEN
       IF(LARGE.EQ.0) THEN
        DMAT1F(3,3)=DMATXF(1,1)/ESTAB
        IF(NTYPE.EQ.4) THEN
         DMAT1F(5,5)=DMATXF(1,1)/ESTAB
         DMAT1F(6,6)=DMATXF(1,1)/ESTAB
        ENDIF
c       CALL CALCLT(DMATXF,DMAT1F,ESTAB)     ! alternative option
       ELSE
        ROTET1=1.0D0/ESTAB
        CALL CALCCT(DMAT1F,YOUNGF,POISMF,ROTET1,    6,
     .              POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .              BULKM1,DISTM1,DAMAP11,DAMAP21,PROPS,
     .              VANIS,ROMTX,TRMTX,
     .              YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
c       CALL CALCLT(DMATXF,DMAT1F,ESTAB)     ! alternative option
        IF(LARGE.NE.0)
     .   CALL FOURTH(DMAT1F,XJACM,XJACI,XJA3M,XJA3I,DETJM,
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
     .               STRAPF,RCGPI,FINET,TRABE,DCINV,EBASE)
       ENDIF           ! large.eq.0
      ENDIF            ! liqui.eq.1
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMTEPF(IKONT)=(1.0D0-DAMAX)*DMAT1F(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** SECANT MARTENSITE CONSTITUTIVE TENSOR TO INTEGRATE THE STRESSES
C
      CALL CALCCT(DMATXA,YOUNGA,POISMA,ROTET,    6,
     .            POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .            BULKM,DISTM,DAMAP1,DAMAP2,PROPS,
     .            VANIS,ROMTX,TRMTX,
     .            YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
C
      ISOTR=INT(PROPS(1))   ! 0=isotropic; 1=orthotropic
      IF(ISOTR.NE.0) THEN
       DO IDIME=1,NDIME     ! transpose of TRMTX
        DO JDIME=1,NDIME
         TRMTT(IDIME,JDIME)=TRMTX(JDIME,IDIME)
        ENDDO
       ENDDO
       CALL FOURTH(DMATXA,XJACM,TRMTT,XJA3M,1.0D0,1.0D0,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAPA,RCGPI,FINET,TRABE,DCINV,EBASE)
      ENDIF
C
C**** COMPUTES THE PRODUCT (BULK * THERMAL DEFORMATION)
C
      BULKA=BULKM*STRRA
C
C**** TRANSFORMS CONSTITUTIVE TENSOR FOR LARGE STRAIN ANALYSIS
C
      IF(LARGE.NE.0)
     . CALL FOURTH(DMATXA,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .             POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .             BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .             EPOL1,EPOL2,EPOL3,
     .             CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .             CPOL9,CPOL10,
     .             CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .             CEMIN,CEMAX,CEETA,
     .             ALAMB,ILAEQ,FVOL1,
     .             RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .             BULKA,    6,UNOMA,FACTJ,DCOFA,UNORA,
     .             STILD,ATILD,VANIS,
     .             STRAPA,RCGPI,FINET,TRABE,DCINV,EBASE)
C
C**** STORE SOME CONSTITUTIVE TENSORS
C
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        DMATEA(ISTRS,JSTRS)=DMATXA(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMATPA(ISTRS,JSTRS)=DMATXA(ISTRS,JSTRS)*(1.0D0-DAMAX)
        DMAT1A(ISTRS,JSTRS)=DMATXA(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** SECANT CONSTITUTIVE TENSOR TO BE USED IN STIFFNESS MATRIX
C
      IF(LIQUI.EQ.1) THEN
       IF(LARGE.EQ.0) THEN
        DMAT1A(3,3)=DMATXA(1,1)/ESTAB
        IF(NTYPE.EQ.4) THEN
         DMAT1A(5,5)=DMATXA(1,1)/ESTAB
         DMAT1A(6,6)=DMATXA(1,1)/ESTAB
        ENDIF
c       CALL CALCLT(DMATXA,DMAT1A,ESTAB)     ! alternative option
       ELSE
        ROTET1=1.0D0/ESTAB
        CALL CALCCT(DMAT1A,YOUNGA,POISMA,ROTET1,    6,
     .              POISC11,POISC21,POISC31,D1PLS1,D2PLS1,D3PLS1,
     .              BULKM1,DISTM1,DAMAP11,DAMAP21,PROPS,
     .              VANIS,ROMTX,TRMTX,
     .              YOUNGF,YOUNGA,POISMF,POISMA,SHEAR,SHEARF,SHEARA)
c       CALL CALCLT(DMATXF,DMAT1F,ESTAB)     ! alternative option
        IF(LARGE.NE.0)
     .   CALL FOURTH(DMAT1A,XJACM,XJACI,XJA3M,XJA3I,DETJM,
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
     .               STRAPA,RCGPI,FINET,TRABE,DCINV,EBASE)
       ENDIF           ! large.eq.0
      ENDIF            ! liqui.eq.1
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMTEPA(IKONT)=(1.0D0-DAMAX)*DMAT1A(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** TOTAL & INCREMENTAL FERRITE STRESS STATE PREDICTION
C
C     Notes: DESIG & STRS0 are not used now. If they have to be used,
C            DESIGF & STRS0F must be implemented.
C            No thermal deformation is considered (otherwise, STRATF
C            must be implemented).
C            DSTRA is not used now. If it has to be used, DSTRAF must
C            be implemented).
C
      CALL CALCST(SIGMAF,DMATXF,STRANF,STRAPF,STRAT,DAMAX,
     .            DETJM,BULKM,DISTM,CPENI,
     .            FVOL1,STRRA,STRR4,
     .            RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,    6,
     .            FACTJ,STILD,ATILD,SGTOTF,PRESG,DESIG,DSTRA,STRS0,
     .            RCGPI,FINET)
C
C**** TOTAL & INCREMENTAL MARTENSITE STRESS STATE PREDICTION
C
C     Note: same as ferrite stress
C
      CALL CALCST(SIGMAA,DMATXA,STRANA,STRAPA,STRAT,DAMAX,
     .            DETJM,BULKM,DISTM,CPENI,
     .            FVOL1,STRRA,STRR4,
     .            RCGTT,RCGTI,TRACC,TRACI,ENEYE,DETJC,SEINC,    6,
     .            FACTJ,STILD,ATILD,SGTOTA,PRESG,DESIG,DSTRA,STRS0,
     .            RCGPI,FINET)
C
C**** TANGENT CONJUGATE OF THERMAL DILATATION TENSOR (TANGENT BETA)
C
      IF(ITERME.GT.0) THEN                      ! bidirectional problems
       IF(ITERMP.GT.0) THEN                     ! mech. coupling terms
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
      RETURN
      END
