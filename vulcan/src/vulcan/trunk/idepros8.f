      SUBROUTINE IDEPROS8(PROPST,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 8 (IPCMO=8) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION PROPST(*)
C
      IA2=IA1+2         ! 2=ipcfo,ipcmo
C
C**** AlN PRECIPITATION
C
      AA=PROPST(IA2+1)
      AQ=PROPST(IA2+2)
      AK=PROPST(IA2+3)
      AR=PROPST(IA2+4)
      AP=PROPST(IA2+5)
C
C**** RECRYSTALIZATION
C
      AAR=PROPST(IA2+6)
      IAQR=INT(PROPST(IA2+7))
      IF(IAQR.EQ.0) THEN
       AQR=PROPST(IA2+8)
       IAUX=1
      ENDIF
      IF(IAQR.EQ.1) THEN
       AQRN=PROPST(IA2+ 8)        ! total nitrogen content
       AQRA=PROPST(IA2+ 9)        ! activation parameter a
       AQRB=PROPST(IA2+10)        ! activation parameter b
       AQRC=PROPST(IA2+11)        ! activation parameter c
       IAUX=4
      ENDIF
      AKR=PROPST(IA2+8+IAUX)
C
      VPLAT(IPLAT, 6)=AA
      VPLAT(IPLAT, 7)=AQ
      VPLAT(IPLAT, 8)=AK
      VPLAT(IPLAT, 9)=AR
      VPLAT(IPLAT,10)=AP
      VPLAT(IPLAT,11)=AAR
      VPLAT(IPLAT,12)=FLOAT(IAQR)
      IF(IAQR.EQ.0) THEN
       VPLAT(IPLAT,13)=AQR
      ENDIF
      IF(IAQR.EQ.1) THEN
       VPLAT(IPLAT,13)=AQRN
       VPLAT(IPLAT,14)=AQRA
       VPLAT(IPLAT,15)=AQRB
       VPLAT(IPLAT,16)=AQRC
      ENDIF
      VPLAT(IPLAT,13+IAUX)=AKR
C
      IMODE=8+IAUX              ! imode=total number of prop. of model 8
C
      IA1=IA2+IMODE
C
      RETURN
      END
