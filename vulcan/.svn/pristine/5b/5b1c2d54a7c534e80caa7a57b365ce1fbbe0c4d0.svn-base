      SUBROUTINE TMSTEPT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,ELVART,
     .                   ELMATT,HEADST,HTLODT,IFFIXT,PRESCT,LNODST,
     .                   MATNOT,PROELT,PROPST,REFORT,RLOADT,RLOAHT,
     .                   FICTOT,TFICTT,TLOADT,DISPLT,PWORKT,PREAST,
     .                   TGAPST,COORDT,TEMPIT,ADVELT,FPCHAT,LACTIT,
     .                   LPNTNT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE ADVANCES IN TIME, UPDATE TIME AND PRESCR. VARIABLES,
C     AND PERFORMS THE PREDICTOR STAGE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT)  :: WORK1T(:)
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

      INCLUDE 'soladdt.inc'

C
      DIMENSION DISITT(NTOTVT),               DISPRT(NTOTVT,3), 
     .          DISTOT(NTOTVT,3),             ELDATT(NDATAT),      
     .          ELPRET(NPREVT),               ELVART(NSTATT),          
     .          ELMATT(NMATXT),               HEADST(NPOINT,4),    
     .          HTLODT(NHLODT,NSUBFT,NFUNCT), IFFIXT(NTOTVT,2),    
     .          LNODST(NNODET,NELEMT),        MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT),        PROPST(NPROPT,NMATST),
     .          REFORT(NTOTVT),               RLOADT(NTOTVT),      
     .          TLOADT(NTOTVT,2),             DISPLT(NTOTVM),
     .          PWORKT(NPOINT,3)
      DIMENSION PREAST(NPREAT,NPOINT),        TGAPST(NPOINT),
     .          COORDT(NDIMET,NPOINT)
      DIMENSION PRESCT(NTOTVT,2),             RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),               TFICTT(NFUNCT)
      DIMENSION TEMPIT(NPOINT,2),             ADVELT(NTOTVT*NDIMET),
     .          FPCHAT(NFPCH,NPOINT),         LACTIT(NELEMT),
     .          LPNTNT(*)
C
      IF(KARCLT.EQ.0.OR.ISTEPT.LE.2) THEN
C
C**** UPDATE DTIME AND TTIME
C
       CALL TMARCHT
C
C**** INCREMENT LOAD FACTORS, APPLIED LOADS AND NON-TENSIONAL STRAINS
C
       CALL INCREMT(ELDATT,ELPRET,ELVART,ELMATT,HEADST,HTLODT,
     .              IFFIXT(1,1),  PRESCT,
     .              LNODST,MATNOT,PROELT,PROPST,RLOADT,RLOAHT,
     .              FICTOT,TFICTT,TLOADT)
C
C**** PREDICTOR PHASE
C
       CALL PREDICT(DISITT,DISPRT,DISTOT)
C
C**** DEALS WITH CONTACT NON-COINCIDENT MESH
C
       IF(ITERME.GT.0) THEN                      ! bidirectional coupled
        IF(ITERMD.GT.0) THEN                     ! deformed shape
         IF(NOCOIT.GT.0) THEN
          CALL SETMTXT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .                 PROPST,COORDT,TEMPIT(1,1),FPCHAT,DISPLT,
     .                 WORK1T)
C
          CALL SOLADDT(IFFIXT,LNODST,LPNTNT,PRESCT,
     .                 WORK1T(IFIXYT( 1)),WORK1T(IFIXYT( 2)),
     .                 WORK1T(IFIXYT( 3)),WORK1T(IFIXYT( 4)),
     .                 WORK1T(IFIXYT( 5)),WORK1T(IFIXYT( 6)),
     .                 WORK1T(IFIXYT( 7)),WORK1T)
         ENDIF
        ENDIF
       ENDIF
C
C**** CALCULATE RESIDUAL FORCES
C
       CALL RESIDUT(DISITT,DISPRT,DISTOT,ELDATT,ELPRET,ELVART,ELMATT,
     .              HEADST,IFFIXT,LNODST,MATNOT,PROELT,PROPST,REFORT,
     .            TLOADT,DISPLT,PWORKT,PREAST(1,1),TGAPST,COORDT,TEMPIT,
     .              ADVELT,FPCHAT,LACTIT,HTLODT,WORK1T)
      ENDIF
C
      RETURN
      END 
