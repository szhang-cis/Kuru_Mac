      SUBROUTINE  INTERVT(DISTOT,ELDATT,ELPRET,ELVART,ELMATT,HTLODT,
     .                    IFFIXT,PRESCT,LNODST,LPNTNT,MATNOT,PROELT,
     .                    PROPST,RLOADT,RLOAHT,FICTOT,TFICTT,ADVELT,
     .                    HEADST,TLOADT,COORDT,TEMPIT,FPCHAT,PREAST,
     .                    DISPLT,LACTIT,WORK1T)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
C     FOR THE FIRST INTERVAL, IT ALSO EVALUATES THE CONSTANT ELEMENT
C     MATRICES OF THE PROBLEM
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT)  :: WORK1T(:)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'


C
      DIMENSION MATNOT(NELEMT),          LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT),   PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),          ELPRET(NPREVT),
     .          ELVART(NSTATT),          ELMATT(NMATXT)
      DIMENSION DISTOT(*),
     .          HTLODT(NHLODT,NSUBFT,*), IFFIXT(NTOTVT,*),
     .          LPNTNT(*),               RLOADT(*)
      DIMENSION PRESCT(NTOTVT,2),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),          TFICTT(NFUNCT)
      DIMENSION ADVELT(NTOTVT*NDIMET)
      DIMENSION HEADST(NPOINT,4),        TLOADT(NTOTVT,2),
     .          COORDT(NDIMET,NPOINT),   TEMPIT(NPOINT,2),
     .          FPCHAT(NFPCH,NPOINT),    PREAST(NPREAT,NPOINT),
     .          DISPLT(NTOTVM),          LACTIT(NELEMT)

      END SUBROUTINE
