      SUBROUTINE INMICRO8(ALPHAM,INUPC,IPLAT)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 8
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
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      AA= VPLAT(IPLAT, 6)
      AQ= VPLAT(IPLAT, 7)
      AK= VPLAT(IPLAT, 8)
      AR= VPLAT(IPLAT, 9)
      AP= VPLAT(IPLAT,10)
      AAR=VPLAT(IPLAT,11)
      IAQR=INT(VPLAT(IPLAT,12))
      IF(IAQR.EQ.0) THEN
       AQR=VPLAT(IPLAT,13)
       IAUX=1
      ENDIF
      IF(IAQR.EQ.1) THEN
       AQRN=VPLAT(IPLAT,13)
       AQRA=VPLAT(IPLAT,14)
       AQRB=VPLAT(IPLAT,15)
       AQRC=VPLAT(IPLAT,16)
       IAUX=4
      ENDIF
      AKR=VPLAT(IPLAT,13+IAUX)
C
      IV=13+IAUX                 ! number of VPLAT defined in input data
C
C**** INITIALITES MICROSTRUCTURAL PARAMETERS ("ALPHAM" & OTHERS)
C
      IN=INUPC
      ALPHAM(IN+1)=0.0D0
      ALPHAM(IN+2)=0.0D0
      ALPHAM(IN+3)=0.0D0
      ALPHAM(IN+4)=0.0D0
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4
C
      RETURN
      END
