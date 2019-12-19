      SUBROUTINE RESIDUT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,ELVART,
     .                   ELMATT,HEADST,IFFIXT,LNODST,MATNOT,PROELT,
     .                   PROPST,REFORT,TLOADT,DISPLT,PWORKT,PREAST,
     .                   TGAPST,COORDT,TEMPIT,ADVELT,FPCHAT,LACTIT,
     .                   HTLODT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESIDUAL FORCES AND REACTIONS
C     IT ALSO PERFORMS THE LINE SEARCH IF DESIRED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          WORK1T(*)
      DIMENSION DISITT(*),             DISPRT(NTOTVT,3), 
     .          DISTOT(NTOTVT,3),      HEADST(NPOINT,*),
     .          IFFIXT(*),             REFORT(NTOTVT,2),
     .          TLOADT(NTOTVT,2)
      DIMENSION DISPLT(NTOTVM),        PWORKT(NPOINT,3),
     .          PREAST(NPOINT),        TGAPST(NPOINT),
     .          COORDT(NDIMET,NPOINT)
      DIMENSION TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),  LACTIT(NELEMT),
     .          HTLODT(NHLODT,NSUBFT,NFUNCT)
C
      CALL CPUTIMT(TIME1T)
C
C**** COMPUTE DISSIPATED ENERGY & IN CASE CORRECT THE LUMPED EXPWP
C
      CALL ENERGYT(DISITT,REFORT)
C
      ALPHAT=1.0D+00                 ! Initialize LINE SEARCH parameters
      ITRILT=0
      LSRCHT=1
C
C**** LINE SEARCH LOOP
C
      DO WHILE (LSRCHT.EQ.1)
C
       ITRILT=ITRILT+1
C
C**** READ LAST CONVERGED DISPLACEMENTS
C
       CALL DATBAST(DISTOT,   12,   2)
C
C**** UPDATE NODAL VARIABLES (BEFORE RESIDUAL EVALUATION)
C
       CALL CORRECT(ALPHAT,DISITT,DISPRT,DISTOT,DTIMET,KINTET,REFORT,
     .              TALFAT,TBETAT,IITERT,     1)
C
C**** COMPUTE THE INTERNAL RESISTING FORCES
C
       CALL FORCINT(DISTOT(1,1),
     .              ELDATT,ELPRET,ELVART,ELMATT,
     .              REFORT(1,1),
     .              HEADST,LNODST,MATNOT,PROELT,PROPST,WORK1T,
     .              REFORT(1,2),
     .              DISTOT(1,1+KDYNAT),DISPRT(1,1+KDYNAT),
     .              DISPLT,PREAST,TGAPST,COORDT,TEMPIT(1,1),ADVELT,
     .              TEMPIT(1,1+NMEMO10),FPCHAT,LACTIT,HTLODT)
C
C**** EVALUATE THE INERTIAL & DAMPING FORCES
C
       IF(KDYNAT.EQ.1) THEN
        CALL THERDYT(DISTOT(1,2),LNODST,MATNOT,PROELT,PROPST,
     .               REFORT(1,1),ELDATT,ELPRET,ELVART,ELMATT,
     .               WORK1T,DISITT,REFORT(1,2),COORDT,
     .               TEMPIT(1,1),ADVELT,FPCHAT,DISTOT(1,1),DISPLT,
     .               LACTIT)
       ENDIF
C
C**** CALCULATE RESIDUAL HEAT
C
       CALL RESLODT(DISITT,HEADST,IFFIXT,REFORT(1,1),TLOADT,PWORKT)
C
C**** PRINT SOME RELEVANT DATA
C
       IF(IITERT.GE.1) THEN
        IF(ITRILT.EQ.1) THEN
         WRITE(LUREST,900) IITERT,GZEROT
         IF(ABS(AACOET)+ABS(BBCOET).NE.0.0D+00) 
     .    WRITE(LUREST,910) AACOET,BBCOET/AACOET
        ENDIF
        RATIOT=0.0D0
        IF(GZEROT.NE.0.0D0) RATIOT=GCURNT/GZEROT
        WRITE(LUREST,905) ITRILT,ALPHAT,RATIOT
       ENDIF
C
C**** FIND A NEW FALUE FOR ALPHA ; IF REQUIRED
C
       LSRCHT=0                             ! Set mark not to loop again
       IF(LINEST.EQ.1.AND.IITERT.GT.1)
     .  CALL LINSCH(ALPHAT,GCURNT,GZEROT,ITRILT,LSRCHT)
C
      END DO                                       ! WHILE (LSRCHT.EQ.1)
C
C**** UPDATE NODAL VARIABLES (AFTER RESIDUAL EVALUATION)
C
      CALL CORRECT(ALPHAT,DISITT,DISPRT,DISTOT,DTIMET,KINTET,REFORT,
     .             TALFAT,TBETAT,IITERT,     2)
C
C**** ALTERNATIVE PHASE-CHANGE FRONT VELOCITY CONTROL ALGORITHM
C
      CALL PHCHFVT(REFORT,DUMMY,     2)
C
      CALL CPUTIMT(TIME2T)
      CPURET=CPURET+(TIME2T-TIME1T)
C
      RETURN 
  900 FORMAT(132('='),//,
     .          10X,'ITERATION NO. ',I5,28X,' GZERO =',E13.5/)
  905 FORMAT(15X,'SEARCH NO. ',I3,3X,'ALPHA =',F10.5,3X,
     .           'GCURN/GZERO =',E13.5)
  910 FORMAT(15X,'SECANT-NEWTON :      A =',F10.5,10X,' B/A =',E13.5)
      END
