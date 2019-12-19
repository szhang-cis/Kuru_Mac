      SUBROUTINE IDEPROS5(PROPST,IPLAT,IA1)
C***********************************************************************
C
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 5 (IPCMO=5) OF RATE PHASE-CHANGE FORMULATIONS
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
      CA =  PROPST(IA2+ 1)
      SIO=  PROPST(IA2+ 2)
      PH =  PROPST(IA2+ 3)
      CU =  PROPST(IA2+ 4)
      QM =  PROPST(IA2+ 5)
      QG =  PROPST(IA2+ 6)
C
      AKSI= PROPST(IA2+ 7)
      IHELAC=DINT(PROPST(IA2+ 8))
C
      INUCMX=INT(PROPST(IA2+ 9))
      IF(INUCMX.EQ.1.OR.                ! Su, Boeri & Lacaze nuc. models
     .   INUCMX.EQ.2.OR.INUCMX.EQ.3) THEN
       ANUCB=PROPST(IA2+10)
       ANUCC=PROPST(IA2+11)
      ENDIF
      INUCAX=INT(PROPST(IA2+12))
      INUCON=INT(PROPST(IA2+13))
      ANUCON=PROPST(IA2+14)
C
      IGROMX=INT(PROPST(IA2+15))
      IF(IGROMX.EQ.1.OR.                ! Su (2) & Lacaze growth models
     .   IGROMX.EQ.2.OR.IGROMX.EQ.3) THEN
       DIFCA=PROPST(IA2+16)
       RSHEL=PROPST(IA2+17)
       RNODO=PROPST(IA2+18)
       AUSGR=PROPST(IA2+19)
      ENDIF
C
      IKMICX=INT(PROPST(IA2+20))
      IKAUX=0
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
       IKAUX=3
       BASKS=PROPST(IA2+21)
       BASKM=PROPST(IA2+22)
       BASKL=PROPST(IA2+23)
      ENDIF
C
      IFPCDT=INT(PROPST(IA2+21+IKAUX))
      IAFLOJ=INT(PROPST(IA2+22+IKAUX))
C
      VPLAT(IPLAT, 6)=CA
      VPLAT(IPLAT, 7)=SIO
      VPLAT(IPLAT, 8)=PH
      VPLAT(IPLAT, 9)=CU
      VPLAT(IPLAT,10)=QM
      VPLAT(IPLAT,11)=QG
C
      VPLAT(IPLAT,12)=AKSI
      VPLAT(IPLAT,13)=FLOAT(IHELAC)
C
      VPLAT(IPLAT,14)=FLOAT(INUCMX)
      IF(INUCMX.EQ.1.OR.                ! Su, Boeri & Lacaze nuc. models
     .   INUCMX.EQ.2.OR.INUCMX.EQ.3) THEN
       VPLAT(IPLAT,15)=ANUCB
       VPLAT(IPLAT,16)=ANUCC
      ENDIF
      VPLAT(IPLAT,17)=FLOAT(INUCAX)
      VPLAT(IPLAT,18)=FLOAT(INUCON)
      VPLAT(IPLAT,19)=ANUCON
C
      VPLAT(IPLAT,20)=FLOAT(IGROMX)
      IF(IGROMX.EQ.1.OR.                ! Su (2) & Lacaze growth models
     .   IGROMX.EQ.2.OR.IGROMX.EQ.3) THEN
       VPLAT(IPLAT,21)=DIFCA
       VPLAT(IPLAT,22)=RSHEL
       VPLAT(IPLAT,23)=RNODO
       VPLAT(IPLAT,24)=AUSGR
      ENDIF
C
      VPLAT(IPLAT,25)=FLOAT(IKMICX)
      IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
       VPLAT(IPLAT,26)=BASKS
       VPLAT(IPLAT,27)=BASKM
       VPLAT(IPLAT,28)=BASKL
      ENDIF
C
      VPLAT(IPLAT,26+IKAUX)=FLOAT(IFPCDT)
      VPLAT(IPLAT,27+IKAUX)=FLOAT(IAFLOJ)
C
      IMODE=22+IKAUX            ! imode=total number of prop. of model 5
C
      IA1=IA2+IMODE
C
      RETURN
      END
