      SUBROUTINE INCREMT(ELDATT,ELPRET,ELVART,ELMATT,HEADST,HTLODT,
     .                   IFFIXT,PRESCT,LNODST,MATNOT,PROELT,PROPST,
     .                   RLOADT,RLOAHT,FICTOT,TFICTT,TLOADT)
C***********************************************************************
C
C**** THIS ROUTINE INCREMENTS THE APPLIED LOADING
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION MATNOT(NELEMT),          LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT),   PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),          ELPRET(NPREVT),
     .          ELVART(NSTATT),          ELMATT(NMATXT)
      DIMENSION HEADST(NPOINT,*), 
     .          HTLODT(NHLODT,NSUBFT,*), IFFIXT(*),       
     .          RLOADT(*),               TLOADT(*)
      DIMENSION PRESCT(NTOTVT,2),        RLOAHT(NTOTVT,NFUNCT),
     .          FICTOT(NFUNCT),          TFICTT(NFUNCT)
C
C**** DETERMINE THE LOAD/DISPLACEMENT INCREMENT FACTOR
C
      CALL INCFACT(HTLODT,FICTOT,TFICTT)
C
C**** INCREMENT THE APPLIED LOAD
C
      CALL INCFORT(HEADST,IFFIXT,RLOADT,RLOAHT,FICTOT,TLOADT)
C
      RETURN
      END
