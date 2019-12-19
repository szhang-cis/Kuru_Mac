      SUBROUTINE CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                   DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                   DSOURT,COUTDT,BASMI ,TGAUIT,PSEUDO)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE DENSITY AND:
C
C     WITH ILAHET=-1, CALCULATE THE DENSITY
C     WITH ILAHET= 0, CALCULATE THE CAPACITY COEFFICIENT
C     WITH ILAHET= 1, CALCULATE THE LATENT HEAT RELEASED
C     WITH ILAHET= 2, CALCULATE THE TEMPERATURE DERIVATIVE OF THE
C                     PHASE-CHANGE FUNCTION (REGULARIZED OR LINEARIZED)
C
C     ISINRT=1 USED IN THE CALCULATION OF JACOBIAN MATRIX
C     ISINRT=2 USED IN THE CALCULATION OF "THERMAL FORCES"
C
C***********************************************************************
C
C     Index of variables
C
C     BASMM = Density
C     BASCC = Capacity coefficient
C     SOUR1T= L*Phase-change function at time t
C     SOUR2T= L*Phase-change function at time t+dt
C     COUTDT= Coupling term due to thermal deformation only. 
C             Term computed in plheat.f (transferred to the thermal 
C             problem in trasmeg.f) and used in matrix C_th for 
C             the improved staggered scheme
C     BASMI = Initial density
C
C     TGAUAT= Last converged temperature at Gauss point (at time t)
C     TGAUST= Current temperature at Gauss point (at time t+alpha*dt for
C             ILAHET=-1,0)
C     TGAUSX= Current temperature at Gauss point (at time t+dt for
C             ILAHET=1,2)
C     TGAUIT= Initial temperature at Gauss point
C
C
C     Free energy models:
C
C     IFRENT=1,2,3  input data: tangent capacity coefficient
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*)
      DIMENSION TSOE1(5), TSOE2(5)
      DIMENSION TSOE1FI(5), TSOE2FI(5)
C
C**** DEALS WITH LINEARIZED PHASE-CHANGE
C
      IF(NPLAT.GT.0) THEN
       IF(LINPC.EQ.1) THEN
        ILAHEX=ILAHET
        IF(ILAHEX.EQ.1) ILAHET=2
       ENDIF
      ENDIF
C
C**** COMPUTES DENSITY
C
      IF(TGAUST.LE.VDENS(1,2)) BASMM=VDENS(1,1)
      DO IDENS=2,NDENS
       I1=IDENS-1
       I2=IDENS
       IF(TGAUST.GT.VDENS(I1,2).AND.TGAUST.LE.VDENS(I2,2)) 
     .    BASMM=(VDENS(I2,1)-VDENS(I1,1))/(VDENS(I2,2)-VDENS(I1,2))*
     .          (TGAUST-VDENS(I1,2))+VDENS(I1,1)
      ENDDO
      IF(TGAUST.GT.VDENS(NDENS,2)) BASMM=VDENS(NDENS,1)
C
C**** COMPUTES INITIAL DENSITY
C
      IF(TGAUIT.LE.VDENS(1,2)) BASMI=VDENS(1,1)
      DO IDENS=2,NDENS
       I1=IDENS-1
       I2=IDENS
       IF(TGAUIT.GT.VDENS(I1,2).AND.TGAUIT.LE.VDENS(I2,2)) 
     .    BASMI=(VDENS(I2,1)-VDENS(I1,1))/(VDENS(I2,2)-VDENS(I1,2))*
     .          (TGAUIT-VDENS(I1,2))+VDENS(I1,1)
      ENDDO
      IF(TGAUIT.GT.VDENS(NDENS,2)) BASMI=VDENS(NDENS,1)
C
      IF(IFILL.EQ.1) TGAUXX=TGAUST     ! stores temp. at time t+alpha*dt
C
      IF(ILAHET.EQ.0) THEN
C
C**** ILAHET=0, CALCULATES THE CAPACITY COEFFICIENT AT TIME (t+alpha dt)
C
       IF(TGAUST.LE.VCAPA(1,2)) BASCC=VCAPA(1,1)
       DO ICAPA=2,NCAPA
        I1=ICAPA-1
        I2=ICAPA
        IF(TGAUST.GT.VCAPA(I1,2).AND.TGAUST.LE.VCAPA(I2,2)) 
     .     BASCC=(VCAPA(I2,1)-VCAPA(I1,1))/(VCAPA(I2,2)-VCAPA(I1,2))*
     .           (TGAUST-VCAPA(I1,2))+VCAPA(I1,1)
       ENDDO
       IF(TGAUST.GT.VCAPA(NCAPA,2)) BASCC=VCAPA(NCAPA,1)
C
C**** IMPROVED STAGGERED SCHEME FOR THERMOMECHANICAL COUPLED PROBLEMS
C     (isothermal split)
C
C     Note: the consideration of the initial density in COUTD is
C           necessary for large strains (it is also correct for small
C           strains)
C
       IF(ITERME.GT.0) THEN
        IF(ITERMP.GT.0) THEN
         IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
          ISIMP=1                       ! better as input (see plheat.f)
          IF(ISIMP.EQ.0.OR.ISIMP.EQ.1) THEN   ! simplif. ther. & general
           BASCC=BASCC+COUTDT/BASMI
          ENDIF
          IF(ISIMP.EQ.2) THEN ! more general but inconsistent with coute
           BASCC=BASCC+(TGAUST+CENKEL)*COUTDT/BASMI     ! (see plheat.f)
          ENDIF
         ENDIF
        ENDIF
       ENDIF
C
      ENDIF                                 ! ilahet.eq.0
C
      IF(ILAHET.EQ.1) THEN
C
C**** ILAHET=1, CALCULATES THE LATENT HEAT RELEASED
C
       SOUR1T=0.0D0
       SOUR2T=0.0D0
       DSOURT=0.0D0
       IF(NPLAT.EQ.0) GOTO 1000
C
       TGAUST=TGAUSX                       ! at time t+dt
       TGAUAT=TGAUST-DTEMPT*DTIMET         ! at time t
C
       DO IPLAT=1,NPLAT
        TSOE1(IPLAT)=0.0D0
        TSOE2(IPLAT)=0.0D0
        TEINF=VPLAT(IPLAT,1)
        TESUP=VPLAT(IPLAT,2)
        HENER=VPLAT(IPLAT,3)
        IPCFO=INT(VPLAT(IPLAT,4))
        IPCMO=INT(VPLAT(IPLAT,5))
        DELEE=TESUP-TEINF
C
        IF(IPCMO.EQ.-2) THEN      ! isothermal phase-change for residual
         IF(ISINRT.EQ.2) THEN
          TEINF=(TEINF+TESUP)/2.0D0
          TESUP=TEINF
         ENDIF
        ENDIF
C
        IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1.OR.IPCMO.EQ.-2) THEN  ! f_pc linear
C
         IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .    TSOE1(IPLAT)=1.0D0/DELEE*(TGAUST-TEINF)
         IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
         IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .    TSOE2(IPLAT)=1.0D0/DELEE*(TGAUAT-TEINF)
         IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
C
         IF(ICONVT.EQ.1) THEN
          IF(ICUBIC.EQ.0) THEN    ! fpc cubical (coherent with ILAHET=2)
           IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUST
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE1(IPLAT)=1.0D0-ROTET
           ENDIF
           IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D0
C
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUAT
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE2(IPLAT)=1.0D0-ROTET
           ENDIF
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D0
          ENDIF                             ! icubic.eq.0
         ENDIF                              ! iconvt.eq.1
C
         TSOE1(IPLAT)=TSOE1(IPLAT)*HENER
         TSOE2(IPLAT)=TSOE2(IPLAT)*HENER
C
        ELSE                                ! f_pc general
         CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
     .                IPLAT,HENER,TSOE1,TSOE2,
     .                TEINF,TESUP,DELEE,IPCMO,    0)
        ENDIF                      ! ipcmo.eq.0
       ENDDO                       ! iplat=1,nplat
C
       DO IPLAT=1,NPLAT
        SOUR1T=SOUR1T+TSOE2(IPLAT)
        SOUR2T=SOUR2T+TSOE1(IPLAT)
       ENDDO
C
C**** PHASE-CHANGE FUNCTION RATE (*DTIMET)
C
       DSOURT=SOUR2T-SOUR1T
C
      ENDIF                                 ! ilahet.eq.1
C
      IF(ILAHET.EQ.2) THEN
C
C**** ILAHET=2, CALCULATES THE TEMPERATURE DERIVATIVE OF THE
C               PHASE-CHANGE FUNCTION (REGULARIZED OR LINEARIZED)
C
C     Notes: ILAHET=2 is only needed for ICONVT=1 (regularized) or
C            LINPC=1 (linearized).
C            Computations at time (t+dt) are only needed for reg. pc.
C            Computations at time (t)    are only needed for lin. pc.
C
       SOUR1T=0.0D0
       SOUR2T=0.0D0
       DSOURT=0.0D0
       IF(NPLAT.EQ.0) GOTO 1000
C
       TGAUST=TGAUSX                       ! at time t+dt
       TGAUAT=TGAUST-DTEMPT*DTIMET         ! at time t
C
       DO IPLAT=1,NPLAT
        TSOE1(IPLAT)=0.0D0
        TSOE2(IPLAT)=0.0D0
        TEINF=VPLAT(IPLAT,1)
        TESUP=VPLAT(IPLAT,2)
        HENER=VPLAT(IPLAT,3)
        IPCFO=INT(VPLAT(IPLAT,4))
        IPCMO=INT(VPLAT(IPLAT,5))
        DELEE=TESUP-TEINF
C
        IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN  ! f_pc linear (not good!)
C
         IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .    TSOE1(IPLAT)=1.0D0/DELEE
         IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
         IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .    TSOE2(IPLAT)=1.0D0/DELEE
         IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
C
         IF(ICONVT.EQ.1) THEN
          IF(ICUBIC.EQ.0) THEN               ! fpc cubical
           IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUST
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE1(IPLAT)=-DEROT
           ENDIF
           IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0D0
C
           IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUAT
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE2(IPLAT)=-DEROT
           ENDIF
           IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0D0
          ENDIF                             ! icubic.eq.0
         ENDIF                              ! iconvt.eq.1
C
         TSOE1(IPLAT)=TSOE1(IPLAT)*HENER
         TSOE2(IPLAT)=TSOE2(IPLAT)*HENER                 ! not necessary
C
        ELSE                                ! f_pc general
         CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
     .                IPLAT,HENER,TSOE1,TSOE2,
     .                TEINF,TESUP,DELEE,IPCMO,    0)
        ENDIF                               ! ipcmo.eq.0
       ENDDO                                ! iplat=1,nplat
C
       DO IPLAT=1,NPLAT
        SOUR1T=SOUR1T+TSOE2(IPLAT)
        SOUR2T=SOUR2T+TSOE1(IPLAT)
       ENDDO
C
       IF(LINPC.EQ.1) THEN                  ! deals with linearized pc
        IF(ILAHEX.EQ.1) THEN                ! see fruf05t.f & mauf05t.f
         SOUR2T=SOUR1T
         SOUR1T=0.0D0
         DSOURT=SOUR2T-SOUR1T
        ENDIF
        IF(ILAHEX.EQ.2) THEN                ! see fruc05t.f & stuc05t.f
         SOUR2T=SOUR1T
        ENDIF
        ILAHET=ILAHEX
       ENDIF
C
      ENDIF                                 ! ilahet.eq.2
C
 1000 CONTINUE
C
C**** DEALS WITH FILLING MATERIAL
C
      IF(IFILL.EQ.0) RETURN
C
      TGAUST=TGAUXX                  ! recovers temp. at time t+alpha*dt
C
C**** DEALS WITH LINEARIZED PHASE-CHANGE
C
      IF(NPLATFI.GT.0) THEN
       IF(LINPCFI.EQ.1) THEN
        ILAHEX=ILAHET
        IF(ILAHEX.EQ.1) ILAHET=2
       ENDIF
      ENDIF
C
C**** COMPUTES DENSITY
C
      IF(TGAUST.LE.VDENSFI(1,2)) BASMMFI=VDENSFI(1,1)
      DO IDENS=2,NDENSFI
       I1=IDENS-1
       I2=IDENS
       IF(TGAUST.GT.VDENSFI(I1,2).AND.TGAUST.LE.VDENSFI(I2,2)) 
     .  BASMMFI=(VDENSFI(I2,1)-VDENSFI(I1,1))/
     .          (VDENSFI(I2,2)-VDENSFI(I1,2))*
     .          (TGAUST-VDENSFI(I1,2))+VDENSFI(I1,1)
      ENDDO
      IF(TGAUST.GT.VDENSFI(NDENSFI,2)) BASMMFI=VDENSFI(NDENSFI,1)
C
C**** COMPUTES INITIAL DENSITY
C
      IF(TGAUIT.LE.VDENSFI(1,2)) BASMIFI=VDENSFI(1,1)
      DO IDENS=2,NDENSFI
       I1=IDENS-1
       I2=IDENS
       IF(TGAUIT.GT.VDENSFI(I1,2).AND.TGAUIT.LE.VDENSFI(I2,2)) 
     .  BASMIFI=(VDENSFI(I2,1)-VDENSFI(I1,1))/
     .          (VDENSFI(I2,2)-VDENSFI(I1,2))*
     .          (TGAUIT-VDENSFI(I1,2))+VDENSFI(I1,1)
      ENDDO
      IF(TGAUIT.GT.VDENSFI(NDENSFI,2)) BASMIFI=VDENSFI(NDENSFI,1)
C
      IF(ILAHET.EQ.0) THEN
C
C**** ILAHET=0, CALCULATES THE CAPACITY COEFFICIENT AT TIME (t+alpha dt)
C
       IF(TGAUST.LE.VCAPAFI(1,2)) BASCCFI=VCAPAFI(1,1)
       DO ICAPA=2,NCAPAFI
        I1=ICAPA-1
        I2=ICAPA
        IF(TGAUST.GT.VCAPAFI(I1,2).AND.TGAUST.LE.VCAPAFI(I2,2)) 
     .   BASCCFI=(VCAPAFI(I2,1)-VCAPAFI(I1,1))/
     .           (VCAPAFI(I2,2)-VCAPAFI(I1,2))*
     .           (TGAUST-VCAPAFI(I1,2))+VCAPAFI(I1,1)
       ENDDO
       IF(TGAUST.GT.VCAPAFI(NCAPAFI,2)) BASCCFI=VCAPAFI(NCAPAFI,1)
C
C**** IMPROVED STAGGERED SCHEME FOR THERMOMECHANICAL COUPLED PROBLEMS
C     (isothermal split)
C
       IF(ITERME.GT.0) THEN
        IF(ITERMP.GT.0) THEN
         IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
          ISIMP=1                       ! better as input (see plheat.f)
          IF(ISIMP.EQ.0.OR.ISIMP.EQ.1) THEN   ! simplif. ther. & general
           BASCCFI=BASCCFI+COUTDT/BASMM
          ENDIF
          IF(ISIMP.EQ.2) THEN ! more general but inconsistent with coute
           BASCCFI=BASCCFI+(TGAUST+CENKEL)*COUTDT/BASMM ! (see plheat.f)
          ENDIF
         ENDIF
        ENDIF
       ENDIF
C
      ENDIF                                 ! ilahet.eq.0
C
      IF(ILAHET.EQ.1) THEN
C
C**** ILAHET=1, CALCULATES THE LATENT HEAT RELEASED
C
       SOUR1TFI=0.0D0
       SOUR2TFI=0.0D0
       DSOURTFI=0.0D0
       IF(NPLATFI.EQ.0) GOTO 2000
C
       TGAUST=TGAUSX                       ! at time t+dt
       TGAUAT=TGAUST-DTEMPT*DTIMET         ! at time t
C
       DO IPLAT=1,NPLATFI
        TSOE1FI(IPLAT)=0.0D0
        TSOE2FI(IPLAT)=0.0D0
        TEINF=VPLATFI(IPLAT,1)
        TESUP=VPLATFI(IPLAT,2)
        HENER=VPLATFI(IPLAT,3)
        IPCFO=INT(VPLATFI(IPLAT,4))
        IPCMO=INT(VPLATFI(IPLAT,5))
        DELEE=TESUP-TEINF
C
        IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN  ! f_pc linear
C
         IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .    TSOE1FI(IPLAT)=1.0D0/DELEE*(TGAUST-TEINF)
         IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
         IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .    TSOE2FI(IPLAT)=1.0D0/DELEE*(TGAUAT-TEINF)
         IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D0
C
         IF(ICONVT.EQ.1) THEN
          IF(ICUBIC.EQ.0) THEN    ! fpc cubical (coherent with ILAHET=2)
           IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUST
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE1FI(IPLAT)=1.0D0-ROTET
           ENDIF
           IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D0
C
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUAT
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE2FI(IPLAT)=1.0D0-ROTET
           ENDIF
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D+00
          ENDIF                             ! icubic.eq.0
         ENDIF                              ! iconvt.eq.1
C
         TSOE1FI(IPLAT)=TSOE1FI(IPLAT)*HENER
         TSOE2FI(IPLAT)=TSOE2FI(IPLAT)*HENER
C
        ELSE                                ! f_pc general
         CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
     .                IPLAT,HENER,TSOE1FI,TSOE2FI,
     .                TEINF,TESUP,DELEE,IPCMO,    1)
        ENDIF                               ! ipcmo.eq.0
       ENDDO                                ! iplat=1,nplat
C
       DO IPLAT=1,NPLATFI
        SOUR1TFI=SOUR1TFI+TSOE2FI(IPLAT)
        SOUR2TFI=SOUR2TFI+TSOE1FI(IPLAT)
       ENDDO
C
C**** PHASE-CHANGE FUNCTION RATE (*DTIMET)
C
       DSOURTFI=SOUR2TFI-SOUR1TFI
C
      ENDIF                                 ! ilahet.eq.1
C
      IF(ILAHET.EQ.2) THEN
C
C**** ILAHET=2, CALCULATES THE TEMPERATURE DERIVATIVE OF THE
C               PHASE-CHANGE FUNCTION (REGULARIZED OR LINEARIZED)
C
C     Notes: ILAHET=2 is only needed for ICONVT=1 (regularized) or
C            LINPCFI=1 (linearized).
C            Computations at time (t+dt) are only needed for reg. pc.
C            Computations at time (t)    are only needed for lin. pc.
C
       SOUR1TFI=0.0D0
       SOUR2TFI=0.0D0
       DSOURTFI=0.0D0
       IF(NPLATFI.EQ.0) GOTO 2000
C
       TGAUST=TGAUSX                       ! at time t+dt
       TGAUAT=TGAUST-DTEMPT*DTIMET         ! at time t
C
       DO IPLAT=1,NPLATFI
        TSOE1FI(IPLAT)=0.0D0
        TSOE2FI(IPLAT)=0.0D0
        TEINF=VPLATFI(IPLAT,1)
        TESUP=VPLATFI(IPLAT,2)
        HENER=VPLATFI(IPLAT,3)
        IPCFO=INT(VPLATFI(IPLAT,4))
        IPCMO=INT(VPLATFI(IPLAT,5))
        DELEE=TESUP-TEINF
C
        IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN  ! f_pc linear (not good!)
C
         IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     .    TSOE1FI(IPLAT)=1.0D0/DELEE
         IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=0.0D0
C
         IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     .    TSOE2FI(IPLAT)=1.0D0/DELEE
         IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=0.0D0
C
         IF(ICONVT.EQ.1) THEN
          IF(ICUBIC.EQ.0) THEN               ! fpc cubical
           IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0D0
           IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUST
           AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE1FI(IPLAT)=-DEROT
           ENDIF
           IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=0.0D0
C
           IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0D0
           IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
            TEMRO=TEINF
            TEMPL=TESUP
            TEMPG=TGAUAT
            AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-
     .                                                3.0D0*TEMRO*TEMRO*
     .       (TEMPL-TEMRO)-1.5D0*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
     .                                        2.0D0*TEMRO*(TEMPL-TEMRO))
            ALP11=-1.0D0/AAMAY
            ALP22=-1.5D0*(TEMPL+TEMRO)*ALP11
            ALP33=-3.0D0*ALP11*TEMRO*TEMRO-2.0D0*ALP22*TEMRO
            ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
            ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
     .            ALP44
            DEROT=3.0D0*ALP11*TEMPG*TEMPG+2.0D0*ALP22*TEMPG+ALP33
            TSOE2FI(IPLAT)=-DEROT
           ENDIF
           IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=0.0D0
          ENDIF                             ! icubic.eq.0
         ENDIF                              ! iconvt.eq.1
C
         TSOE1FI(IPLAT)=TSOE1FI(IPLAT)*HENER
         TSOE2FI(IPLAT)=TSOE2FI(IPLAT)*HENER
C
        ELSE                                ! f_pc general
         CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
     .                IPLAT,HENER,TSOE1FI,TSOE2FI,
     .                TEINF,TESUP,DELEE,IPCMO,    1)
        ENDIF                               ! ipcmo.eq.0
       ENDDO                                ! iplat=1,nplat
C
       DO IPLAT=1,NPLATFI
        SOUR1TFI=SOUR1TFI+TSOE2FI(IPLAT)
        SOUR2TFI=SOUR2TFI+TSOE1FI(IPLAT)
       ENDDO
C
       IF(LINPCFI.EQ.1) THEN                ! deals with linearized pc
        IF(ILAHEX.EQ.1) THEN                ! see fruf05t.f & mauf05t.f
         SOUR2TFI=SOUR1TFI
         SOUR1TFI=0.0D0
         DSOURTFI=SOUR2TFI-SOUR1TFI
        ENDIF
        IF(ILAHEX.EQ.2) THEN                ! see fruc05t.f & stuc05t.f
         SOUR2TFI=SOUR1TFI
        ENDIF
        ILAHET=ILAHEX
       ENDIF
C
      ENDIF                                 ! ilahet.eq.2
C
 2000 CONTINUE
C
C**** COMPUTES THE WEIGHTED CONDUCTIVITY COEFFICIENT
C
      BASMM=PSEUDO*BASMM+(1.0D0-PSEUDO)*BASMMFI
      BASMI=PSEUDO*BASMI+(1.0D0-PSEUDO)*BASMIFI
      IF(ILAHET.EQ.0) THEN
       BASCC=PSEUDO*BASCC+(1.0D0-PSEUDO)*BASCCFI
      ENDIF
      IF(ILAHET.EQ.1) THEN
       SOUR1T=PSEUDO*SOUR1T+(1.0D0-PSEUDO)*SOUR1TFI
       SOUR2T=PSEUDO*SOUR2T+(1.0D0-PSEUDO)*SOUR2TFI
       DSOURT=PSEUDO*DSOURT+(1.0D0-PSEUDO)*DSOURTFI
      ENDIF
      IF(ILAHET.EQ.2.OR.ILAHET.EQ.3) THEN
       SOUR1T=PSEUDO*SOUR1T+(1.0D0-PSEUDO)*SOUR1TFI
       SOUR2T=PSEUDO*SOUR2T+(1.0D0-PSEUDO)*SOUR2TFI
      ENDIF
C
      RETURN
      END
