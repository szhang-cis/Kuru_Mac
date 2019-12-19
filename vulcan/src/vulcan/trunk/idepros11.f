      SUBROUTINE IDEPROS11(PROPST,IPLAT,IA1,INUPM)
C***********************************************************************
C     
C**** THIS ROUTINE ORDERS THE MICROSTRUCTURAL PROPERTIES OF MODEL
C     NUMBER 11 (IPCMO=11) OF RATE PHASE-CHANGE FORMULATIONS
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'       ! thermal-microestructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION PROPST(*)
C
      IA2=IA1+2                 ! 2=ipcfo,ipcmo
C
      IX=0
      NCOPC=IPLAC(INUPM,1,1)
      IF(NCOPC.EQ.0) THEN
       CA=PROPST(IA2+1)
       SI=PROPST(IA2+2)
       CU=PROPST(IA2+3)
       QM=PROPST(IA2+4)
       CR=PROPST(IA2+5)
       AN=PROPST(IA2+6)
       IX=6
      ENDIF
C
      TEFERI=PROPST(IA2+IX+1)      ! reference ferrite temperature
      HENERF=PROPST(IA2+IX+2)      ! ferrite latent heat
      IFENU=INT(PROPST(IA2+IX+3))  ! ferrite nucleation model
      IF(IFENU.EQ.1) THEN
       OMEGAF=PROPST(IA2+IX+4)     ! omega_a
      ENDIF
      IFEGR=INT(PROPST(IA2+IX+5))  ! ferrite growth model
      IF(IFEGR.EQ.1) THEN
       RFER0=PROPST(IA2+IX+6)      ! initial ferrite radius
      ENDIF
C
      TEPERI=PROPST(IA2+IX+7)      ! reference pearlite temperature
      HENERP=PROPST(IA2+IX+8)      ! pearlite latent heat
      IPENU=INT(PROPST(IA2+IX+9))  ! pearlite nucleation model
      IF(IPENU.EQ.1) THEN
       OMEGAP=PROPST(IA2+IX+10)     ! omega_p
      ENDIF
      IPEGR=INT(PROPST(IA2+IX+11)) ! pearlite growth model
      IF(IPEGR.EQ.1) THEN
       DELTAP=PROPST(IA2+IX+12)    ! interface thickness
       DIFFBC=PROPST(IA2+IX+13)    ! C diffusion at the pearlite bound.
       SIGMAP=PROPST(IA2+IX+14)    ! ferrite/cementite interfacial ener.
      ENDIF
C
      IFPCDT=INT(PROPST(IA2+IX+15))
      IAFLOJ=INT(PROPST(IA2+IX+16))
C
      IF(NCOPC.EQ.0) THEN
       VPLAT(IPLAT, 6)=CA
       VPLAT(IPLAT, 7)=SI
       VPLAT(IPLAT, 8)=CU
       VPLAT(IPLAT, 9)=QM
       VPLAT(IPLAT,10)=CR
       VPLAT(IPLAT,11)=AN
      ENDIF
C
      VPLAT(IPLAT, 6+IX)=TEFERI
      VPLAT(IPLAT, 7+IX)=HENERF
      VPLAT(IPLAT, 8+IX)=FLOAT(IFENU)
      VPLAT(IPLAT, 9+IX)=OMEGAF
      VPLAT(IPLAT,10+IX)=FLOAT(IFEGR)
      VPLAT(IPLAT,11+IX)=RFER0
C
      VPLAT(IPLAT,12+IX)=TEPERI
      VPLAT(IPLAT,13+IX)=HENERP
      VPLAT(IPLAT,14+IX)=FLOAT(IPENU)
      VPLAT(IPLAT,15+IX)=OMEGAP
      VPLAT(IPLAT,16+IX)=FLOAT(IPEGR)
      VPLAT(IPLAT,17+IX)=DELTAP
      VPLAT(IPLAT,18+IX)=DIFFBC
      VPLAT(IPLAT,19+IX)=SIGMAP
C
      VPLAT(IPLAT,20+IX)=FLOAT(IFPCDT)
      VPLAT(IPLAT,21+IX)=FLOAT(IAFLOJ)
C
      IMODE=IX+16              ! imode=total number of prop. of model 11
C
      IA1=IA2+IMODE
C
      RETURN
      END
