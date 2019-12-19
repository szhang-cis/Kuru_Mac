      SUBROUTINE SOLADDT(IFFIXT,LNODST,LPNTNT,PRESCT,
     .                   LEQNST,LNUEQT,LOCELT,LPONTT,NACVAT,NDESTT,
     .                   NDFROT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE ALLOCATES MEMORY FOR SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE  :: WORK1T(:)

C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'addsolt.inc'
C
      DIMENSION IFFIXT(NTOTVT,*), LPNTNT(*), LNODST(NNODET,*),
     .          PRESCT(NTOTVT,2)
      DIMENSION LEQNST(*),        LNUEQT(*), LOCELT(*), 
     .          LPONTT(*),        NACVAT(*), NDESTT(*),
     .          NDFROT(*)
C
      IMPOST=0
C
      DO 10 ITOTVT=1,NTOTVT
   10 IF(IFFIXT(ITOTVT,1).EQ.1.AND.PRESCT(ITOTVT,1).NE.0) IMPOST=1
C
C**** SKYLINE & DPCG SOLUTION
C
      IF(KSOLVT.EQ.0) THEN
       CALL SKYDEST(LUSOLT,IFFIXT(1,1),LNODST,LPNTNT,LNUEQT,LPONTT,
     .              NDOFCT,
     .              NELEMT,NEQNST,     NLASTT,NNODET,NPOINT,NTOTVT,
     .              NWIDTT,
     .              LEQNST,KRENUT)
       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
       IF(KSMUST.GT.0)
     .  CALL SKYDEST(LUREST,IFFIXT(1,2),LNODST,LPNTNT,LNUEQT,LPONTT,
     .               NDOFCT,
     .               NELEMT,NEQNST,     NLASTT,NNODET,NPOINT,NTOTVT,
     .               NWIDTT,
     .               LEQNST,KRENUT)
      ENDIF
C
C**** FRONTAL SOLUTION
C
      IF(KSOLVT.EQ.1) THEN
       CALL FROAPPT(LNODST,NBUFAT,NDFROT,NDOFCT,NELEMT,NFRONT,NNODET,
     .              NPOINT,NSTIFT,KSYMMT)
       CALL FRODEST(LOCELT,LNODST,NACVAT,NBUFAT,NDESTT,NDOFCT,NELEMT,
     .              NEVACT,NFRONT,NNODET,NSTIFT)
       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
      ENDIF
C
C**** PCG SOLUTION
C
      IF(KSOLVT.EQ.2) THEN
       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
      ENDIF
C
C**** GMRES SOLUTION
C
C     Note: similar to SKYLINE but some variables are not necessary
C           (see addsolt.f & gmresst.f)
C
      IF(KSOLVT.EQ.3) THEN
       CALL SKYDEST(LUSOLT,IFFIXT(1,1),LNODST,LPNTNT,LNUEQT,LPONTT,
     .              NDOFCT,
     .              NELEMT,NEQNST,     NLASTT,NNODET,NPOINT,NTOTVT,
     .              NWIDTT,
     .              LEQNST,KRENUT)
       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
      ENDIF
C
C**** EXPLICIT SOLUTION (NO SOLVER)
C
      IF(KSOLVT.EQ.4) THEN
       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT,WORK1T)
      ENDIF
C
      RETURN
      END
