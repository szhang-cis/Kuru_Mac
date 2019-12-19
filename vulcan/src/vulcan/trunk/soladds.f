      SUBROUTINE SOLADDS(IFFIXS,LNODSS,LPNTNS,PRESCS,
     .                   LEQNSS,LNUEQS,LOCELS,LPONTS,NACVAS,NDESTS,
     .                   NDFROS,WORK1S)
C***********************************************************************
C
C**** THIS ROUTINE ALLOCATES MEMORY FOR SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1S(:)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'addsols.inc'
C
      DIMENSION IFFIXS(NTOTVS,*), LPNTNS(*), LNODSS(NNODES,*),
     .          PRESCS(NTOTVS,2)
      DIMENSION LEQNSS(*),        LNUEQS(*), LOCELS(*), 
     .          LPONTS(*),        NACVAS(*), NDESTS(*),
     .          NDFROS(*)
C
      IMPOSS=0
C
      DO 10 ITOTVS=1,NTOTVS
   10 IF(IFFIXS(ITOTVS,1).EQ.1.AND.PRESCS(ITOTVS,1).NE.0) IMPOSS=1
C
C**** SKYLINE & DPCG SOLUTION
C
      IF(KSOLVS.EQ.0) THEN
       call runends('skydest not implemented')
c      CALL SKYDEST(LUSOLT,IFFIXT(1,1),LNODST,LPNTNT,LNUEQT,LPONTT,
c    .              NDOFCT,
c    .              NELEMT,NEQNST,     NLASTT,NNODET,NPOINT,NTOTVT,
c    .              NWIDTT,
c    .              LEQNST,KRENUT)
c       CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT)
c       IF(KSMUST.GT.0)
c    .   CALL SKYDEST(LUREST,IFFIXT(1,2),LNODST,LPNTNT,LNUEQT,LPONTT,
c    .                NDOFCT,
c    .                NELEMT,NEQNST,     NLASTT,NNODET,NPOINT,NTOTVT,
c    .                NWIDTT,
c    .                LEQNST,KRENUT)
      ENDIF
C
C**** FRONTAL SOLUTION
C
      IF(KSOLVS.EQ.1)THEN
       CALL FROAPPS(LNODSS,NBUFAS,NDFROS,NDOFCS,NELEMS,NFRONS,NNODES,
     .              NPOINS,NSTIFS,KSYMMS)
       CALL FRODESS(LOCELS,LNODSS,NACVAS,NBUFAS,NDESTS,NDOFCS,NELEMS,
     .              NEVACS,NFRONS,NNODES,NSTIFS)
       CALL ADDSOLS(NEQNSS,NFRONS,NLASTS,NSTIFS,WORK1S)
      ENDIF
C
C**** PCG SOLUTION
C
      IF(KSOLVT.EQ.2)THEN
       call runends('ksolvt eq 2 not implemented')
c      CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT)
      ENDIF
C
C**** GMRES SOLUTION (not implemented yet)
C

C
C**** EXPLICIT SOLUTION (NO SOLVER)
C
      IF(KSOLVT.EQ.4)THEN
       call runends('ksolvt eq 4 not implemented')
c      CALL ADDSOLT(NEQNST,NFRONT,NLASTT,NSTIFT)
      ENDIF
C
      RETURN
      END
