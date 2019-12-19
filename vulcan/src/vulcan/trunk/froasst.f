      SUBROUTINE FROASST(CSTIF,GLOAD,GSTIF,NDEST,NEVAB,NFRON,NSTIF,
     .                   KSYMM,
     .                   NELEM,NNODE,NPREL,NGRUP,NPROP,NMATS,
     .           NDATA,NPREV,NSTAT,NMATX,NTOTV,NPOIN,NDIME,NTOTVM,NFPCH,
     .                   MATNO,LNODS,PROEL,PROPS,ELDAT,ELPRE,
     .                   ELVAR,ELMAT,WORK1,DISTO,COORD,DISIT,
     .                   ADVEL,TEMPI,PREAS,TGAPS,DISPL,FPCHA,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES THE ELEMENT STIFFNESS MATRIX INTO THE 
C     GLOBAL STIFFNESS ARRAY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), GLOAD(*), GSTIF(*), NDEST(*)
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION DISTO(NTOTV,3),     COORD(NDIME,NPOIN),
     .          DISIT(NTOTV)
      DIMENSION ADVEL(NTOTV*NDIME), TEMPI(NPOIN,2)
      DIMENSION PREAS(NPOIN),       TGAPS(NPOIN)
      DIMENSION DISPL(NTOTVM),      FPCHA(NFPCH,NPOIN),
     .          LACTI(NELEM)
C
      NFUNC(I,J)=(J*J-J)/2+I
      NFUNB(I,J)=(J-1)*NFRON+I
C
C**** READ CSTIF FROM DATA BASE
C 
      IF(NMEMO7.EQ.0) THEN
       CALL DATBAST(CSTIF,    6,    2)
      ELSE
       CALL STIFMST(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .              PROPS,WORK1,DISTO(1,2),DISIT,COORD,
     .              ADVEL,TEMPI(1,1),PREAS,TGAPS,DISTO(1,1),
     .              TEMPI(1,2),DISPL,FPCHA,LACTI,CSTIF)
      ENDIF
C
      DO 10 IEVAB=1,NEVAB
      IDEST=NDEST(IEVAB)
      IF(IDEST.EQ.0) GO TO 20
C
C**** ASSEMBLE THE ELEMENT STIFFNESSES-BUT NOT IN RESOLUTION
C
      IF(KSYMM.EQ.0) THEN        ! unsymmetric case
       ILIMM=NEVAB
      ELSE                       ! symmetric case
       ILIMM=IEVAB
      ENDIF
C
      DO 10 JEVAB=1,ILIMM
      JDEST=NDEST(JEVAB)
      IF(JDEST.EQ.0) GO TO 20
C
      IF(KSYMM.EQ.0) THEN        ! unsymmetric case
       NGASH=NFUNB(JDEST,IDEST)
       GSTIF(NGASH)=GSTIF(NGASH)+CSTIF(IEVAB,JEVAB)
      ELSE                       ! symmetric case
       NGASH=NFUNC(IDEST,JDEST)
       NGISH=NFUNC(JDEST,IDEST)
       IF(JDEST.GE.IDEST) GSTIF(NGASH)=GSTIF(NGASH)+CSTIF(IEVAB,JEVAB)
       IF(JDEST.LT.IDEST) GSTIF(NGISH)=GSTIF(NGISH)+CSTIF(IEVAB,JEVAB)
      ENDIF
C
   20 CONTINUE
   10 CONTINUE
C
      RETURN
      END
