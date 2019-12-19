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
      REAL*8, ALLOCATABLE  :: WORK1T(:)
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
      END SUBROUTINE
