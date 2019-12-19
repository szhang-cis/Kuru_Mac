      SUBROUTINE FROASS(CSTIF,GLOAD,GSTIF,NDEST,NEVAB,NFRON,NSTIF,
     .                  KSYMM,
     .                  NELEM,NNODE,NPREL,NGRUP,NPROP,NMATS,
     .                  NDATA,NPREV,NSTAT,NMATX,NTOTV,NPOIN,NDIME,
     .                  NDISR,NDISO,
     .                  MATNO,LNODS,PROEL,PROPS,ELDAT,ELPRE,
     .                  ELVAR,ELMAT,WORK1,TEMPN,DISPR,DTEMP,
     .                  VNORM,DISTO,COORD,INFRI,COFRI,LACTI)
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
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION CSTIF(NEVAB,NEVAB), GLOAD(*), GSTIF(*), NDEST(*)
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(NPOIN,2),     DISPR(NTOTV,NDISR),
     .          DTEMP(NPOIN),       VNORM(NTOTV)
      DIMENSION DISTO(NTOTV,NDISO), COORD(NDIME,NPOIN)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
      NFUNC(I,J)=(J*J-J)/2+I
      NFUNB(I,J)=(J-1)*NFRON+I
C
C**** READ CSTIF FROM DATA BASE
C 
      IF(NMEMO7M.EQ.0) THEN
       CALL DATBAS(CSTIF,    6,    2)
      ELSE
       CALL STIFMS(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .             PROPS,WORK1,TEMPN(1,1),DISPR(1,1),DTEMP,VNORM,COORD,
     .             DISTO(1,1),INFRI,COFRI,LACTI,CSTIF)
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
