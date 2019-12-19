      SUBROUTINE INMICRO1(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 1
C
C***********************************************************************
C
C     Index of variables
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal-microstructural)
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ALPHAM(NHISTM)
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
      IN=INUPC
c     ALPHAM(IN+1)=1.0D+00   ! only to test femview (not necessary here)
      ALPHAM(IN+1)=0.0D+00   ! only to test femview (not necessary here)
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+1
C
      RETURN
      END
