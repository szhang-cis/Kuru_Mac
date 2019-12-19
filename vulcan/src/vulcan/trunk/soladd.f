      SUBROUTINE SOLADD(IFFIX,LNODS,LPNTN,PRESC,
     .                  LEQNS,LNUEQ,LOCEL,LPONT,NACVA,NDEST,NDFRO,
     .                  NEQNS,WORK1)
C***********************************************************************
C
C**** THIS ROUTINE ALLOCATES MEMORY FOR SOLVER
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
C
      DIMENSION IFFIX(NTOTV,*), LPNTN(*), LNODS(NNODE,*),
     .          PRESC(NTOTV,2)
      DIMENSION LEQNS(*),       LNUEQ(*),       LOCEL(*), 
     .          LPONT(*),       NACVA(*),       NDEST(*),
     .          NDFRO(*)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
      INCLUDE 'addsol.inc'
C
      IMPOS=0
C
      DO 10 ITOTV=1,NTOTV
   10 IF(IFFIX(ITOTV,1).EQ.1.AND.PRESC(ITOTV,1).NE.0) IMPOS=1
C
C**** SKYLINE & DPCG SOLUTION
C
      IF(KSOLV.EQ.0) THEN
       CALL SKYDES(LUSOL,IFFIX(1,1),LNODS,LPNTN,LNUEQ,LPONT,NDOFC,
     .             NELEM,NEQNS,     NLAST,NNODE,NPOIN,NTOTV,NWIDT,
     .             LEQNS,KRENU)
       CALL ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
       IF(KSMUS.GT.0)
     .  CALL SKYDES(LURES,IFFIX(1,2),LNODS,LPNTN,LNUEQ,LPONT,NDOFC,
     .              NELEM,NEQNS,     NLAST,NNODE,NPOIN,NTOTV,NWIDT,
     .              LEQNS,KRENU)
      ENDIF
C
C**** FRONTAL SOLUTION
C
      IF(KSOLV.EQ.1)THEN
       CALL FROAPP(LNODS,NBUFA,NDFRO,NDOFC,NELEM,NFRON,NNODE,
     .             NPOIN,NSTIF,KSYMM)
       CALL FRODES(LOCEL,LNODS,NACVA,NBUFA,NDEST,NDOFC,NELEM,
     .             NEVAC,NFRON,NNODE,NSTIF)
       CALL ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
      ENDIF
C
C**** PCG SOLUTION
C
      IF(KSOLV.EQ.2)THEN
       CALL ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
      ENDIF
C
C**** GMRES SOLUTION
C
C     Note: similar to SKYLINE but some variables are not necessary
C           (see addsol.f & gmress.f)
C
      IF(KSOLV.EQ.3) THEN
       CALL SKYDES(LUSOL,IFFIX(1,1),LNODS,LPNTN,LNUEQ,LPONT,NDOFC,
     .             NELEM,NEQNS,     NLAST,NNODE,NPOIN,NTOTV,NWIDT,
     .             LEQNS,KRENU)
       CALL ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
      ENDIF
      
      IF(KSOLV.EQ.4) THEN
       CALL PARDES(LUSOL,IFFIX(1,1),LNODS,LPNTN,LNUEQ,LPONT,NDOFC,  ! no cambio aparente
     .             NELEM,NEQNS,     NLAST,NNODE,NPOIN,NTOTV,NWIDT,
     .             LEQNS,KRENU)
       CALL ADDSOL(NEQNS,NFRON,NLAST,NSTIF,WORK1)
      ENDIF
C
      RETURN
      END
