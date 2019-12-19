      SUBROUTINE CAPCOFS(BASMM ,BASCC ,PROPSS,TGAUSS,DTEMPS,SOUR1S,
     .                   SOUR2S,ILAHES,ISINRS,DSOURS,COUTDS,BASMI,
     .                   TGAUIS,PSEUDOS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE DENSITY AND:
C
C     WITH ILAHE=-1 CALCULATE THE DENSITY
C     WITH ILAHE= 0, CALCULATE THE CAPACITY COEFFICIENT
C     WITH ILAHE= 1, CALCULATE THE LATENT HEAT RELEASED
C     WITH ILAHE= 2, CALCULATE THE TEMPERATURE DERIVATIVE OF THE
C                    PHASE-CHANGE FUNCTION
C
C     ISINR=1 USED IN THE CALCULATION OF JACOBIAN MATRIX
C     ISINR=2 USED IN THE CALCULATION OF "THERMAL FORCES"
C     (ISINR ACTUALLY NOT USED!!)
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
C     TGAUAT= Last converged temperature at Gauss point
C     TGAUST= Current temperature at Gauss point
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
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION PROPSS(*)
      DIMENSION TSOE1(5), TSOE2(5)
      DIMENSION TSOE1FI(5), TSOE2FI(5)
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
      TGAUXX=TGAUSS
      IF(KINTES.EQ.1)             ! Euler method
     . TGAUSS=TALFAS*TGAUSS+(1.0-TALFAS)*(TGAUSS-DTEMPS*DTIMES)
C
C**** COMPUTES DENSITY
C
      IF(TGAUSS.LE.VDENSS(1,2)) BASMM=VDENSS(1,1)
      DO IDENS=2,NDENSS
       I1=IDENS-1
       I2=IDENS
       IF(TGAUSS.GT.VDENSS(I1,2).AND.TGAUSS.LE.VDENSS(I2,2)) 
     .    BASMM=(VDENSS(I2,1)-VDENSS(I1,1))/(VDENSS(I2,2)-VDENSS(I1,2))*
     .          (TGAUSS-VDENSS(I1,2))+VDENSS(I1,1)
      ENDDO
      IF(TGAUSS.GT.VDENSS(NDENSS,2)) BASMM=VDENSS(NDENSS,1)
C
C**** COMPUTES INITIAL DENSITY
C
      IF(TGAUIS.LE.VDENSS(1,2)) BASMI=VDENSS(1,1)
      DO IDENS=2,NDENSS
       I1=IDENS-1
       I2=IDENS
       IF(TGAUIS.GT.VDENSS(I1,2).AND.TGAUIS.LE.VDENSS(I2,2)) 
     .    BASMI=(VDENSS(I2,1)-VDENSS(I1,1))/(VDENSS(I2,2)-VDENSS(I1,2))*
     .          (TGAUIS-VDENSS(I1,2))+VDENSS(I1,1)
      ENDDO
      IF(TGAUIS.GT.VDENSS(NDENSS,2)) BASMI=VDENSS(NDENSS,1)
C
      IF(ILAHES.EQ.0) THEN
C
C**** ILAHE=0, CALCULATES THE CAPACITY COEFFICIENT AT TIME (t+alpha dt)
C
       IF(TGAUSS.LE.VCAPAS(1,2)) BASCC=VCAPAS(1,1)
       DO ICAPA=2,NCAPAS
        I1=ICAPA-1
        I2=ICAPA
        IF(TGAUSS.GT.VCAPAS(I1,2).AND.TGAUSS.LE.VCAPAS(I2,2)) 
     .    BASCC=(VCAPAS(I2,1)-VCAPAS(I1,1))/(VCAPAS(I2,2)-VCAPAS(I1,2))*
     .          (TGAUSS-VCAPAS(I1,2))+VCAPAS(I1,1)
       ENDDO
       IF(TGAUSS.GT.VCAPAS(NCAPAS,2)) BASCC=VCAPAS(NCAPAS,1)
C
C**** IMPROVED STAGGERED SCHEME FOR THERMOMECHANICAL COUPLED PROBLEMS
C     (isothermal split)
C
C     Note: the consideration of the initial density in COUTD is
C           necessary for large strains (it is also correct for small
C           strains)
C
       IF(ITERME.GT.0) THEN

        call runends('para en capcofs')

c       IF(ITERMP.GT.0) THEN
c        IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
c         ISIMP=1                       ! better as input (see plheat.f)
c         IF(ISIMP.EQ.0.OR.ISIMP.EQ.1) THEN   ! simplif. ther. & general
c          BASCC=BASCC+COUTDT/BASMI
c         ENDIF
c         IF(ISIMP.EQ.2) THEN ! more general but inconsistent with coute
c          BASCC=BASCC+(TGAUST+CENKEL)*COUTDT/BASMI     ! (see plheat.f)
c         ENDIF
c        ENDIF
c       ENDIF
       ENDIF
C
      ENDIF       ! ilahes.eq.0
C
      IF(ILAHES.EQ.1) THEN
C
C**** ILAHE=1, CALCULATE THE LATENT HEAT RELEASED
C

       call runends('ilahes=1 not implemented')

c      SOUR1T=0.0
c      SOUR2T=0.0
c      DSOURT=0.0
c      IF(NPLAT.EQ.0) GOTO 1000
C
c      TGAUST=TGAUXX
c      TGAUAT=TGAUST-DTEMPT*DTIMET
C
c      DO IPLAT=1,NPLAT
c       TSOE1(IPLAT)=0.0
c       TSOE2(IPLAT)=0.0
c       TEINF=VPLAT(IPLAT,1)
c       TESUP=VPLAT(IPLAT,2)
c       HENER=VPLAT(IPLAT,3)
c       IPCFO=INT(VPLAT(IPLAT,4))
c       IPCMO=INT(VPLAT(IPLAT,5))
c       DELEE=TESUP-TEINF
C
c       IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN           ! f_pc linear
C
c        IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0
c        IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
c    .    TSOE1(IPLAT)=1.0D+00/DELEE*(TGAUST-TEINF)
c        IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D+00
C
c        IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0
c        IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
c    .    TSOE2(IPLAT)=1.0D+00/DELEE*(TGAUAT-TEINF)
c        IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D+00
C
c        IF(ICONVT.EQ.1) THEN
c         IF(ICUBIC.EQ.0) THEN    ! fpc cubical (coherent with ILAHET=2)
c          IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0
c          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
c           TEMRO=TEINF
c           TEMPL=TESUP
c           TEMPG=TGAUST
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c           ALP11=-1.00/AAMAY
c           ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c           ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c           ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c           ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .            ALP44
c           DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c           TSOE1(IPLAT)=1.0-ROTET
c          ENDIF
c          IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0D+00
C
c          IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0
c          IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
c           TEMRO=TEINF
c           TEMPL=TESUP
c           TEMPG=TGAUAT
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c           ALP11=-1.00/AAMAY
c           ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c           ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c           ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c           ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .            ALP44
c           DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c           TSOE2(IPLAT)=1.0-ROTET
c          ENDIF
c          IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0D+00
c         ENDIF                              ! icubic.eq.0
c        ENDIF                               ! iconvt.eq.1
C
c        TSOE1(IPLAT)=TSOE1(IPLAT)*HENER
c        TSOE2(IPLAT)=TSOE2(IPLAT)*HENER
C
c       ELSE                                ! f_pc general
c        CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
c    .                IPLAT,HENER,TSOE1,TSOE2,
c    .                TEINF,TESUP,DELEE,IPCMO,    0)
c       ENDIF                      ! ipcmo.eq.0
c      ENDDO                       ! iplat=1,nplat
C
c      DO IPLAT=1,NPLAT
c       SOUR1T=SOUR1T+TSOE2(IPLAT)
c       SOUR2T=SOUR2T+TSOE1(IPLAT)
c      ENDDO
C
C**** PHASE-CHANGE FUNCTION RATE (*DTIMET)
C
c      DSOURT=SOUR2T-SOUR1T
C
      ENDIF                                         ! ilahes.eq.1
C
      IF(ILAHES.EQ.2.or.ilahes.eq.3) THEN

       call runends('ilahes=2,3 not implemented')

C
C**** ILAHE=2, CALCULATE THE TEMPERATURE DERIVATIVE OF THE
C              PHASE-CHANGE FUNCTION
C
C     Note: ICONVT=1 for ILAHET=2,3
C
c      SOUR1T=0.0                                        ! not necessary
c      SOUR2T=0.0
c      IF(NPLAT.EQ.0) GOTO 1000
C
c      TGAUST=TGAUXX
c      TGAUAT=TGAUST-DTEMPT*DTIMET
C
c      DO IPLAT=1,NPLAT
c       TSOE1(IPLAT)=0.0
c       TSOE2(IPLAT)=0.0                                 ! not necessary
c       TEINF=VPLAT(IPLAT,1)
c       TESUP=VPLAT(IPLAT,2)
c       HENER=VPLAT(IPLAT,3)
c       IPCFO=INT(VPLAT(IPLAT,4))
c       IPCMO=INT(VPLAT(IPLAT,5))
c       DELEE=TESUP-TEINF
C
c       IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN     ! f_pc linear (not good!)
C
c        IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0
c        IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
c    .    TSOE1(IPLAT)=1.0D+00/DELEE
c        IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0
C
c        IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0            ! not necessary
c        IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
c    .    TSOE2(IPLAT)=1.0D+00/DELEE
c        IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0
C
c        icubit=0                ! to control df_pc/dT (better as input)
c        if(ilahet.eq.4) icubit=0           ! could be ilahet.eq.3
c        icubib=icubic
c        if(ilahet.eq.3) icubib=1
C
c        IF(ICUBIB.EQ.1) THEN               ! fpc cubical
c         IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0
c         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
c          TEMRO=TEINF
c          TEMPL=TESUP
c          TEMPG=TGAUST
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c          ALP11=-1.00/AAMAY
c          ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c          ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c          ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c          ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .           ALP44
c          DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c          if(icubit.eq.1) then
c           TSOE1(IPLAT)=-DEROT
c          else
c           tsoea=-derot
c           if(tsoea.lt.tsoe1(iplat)) TSOE1(IPLAT)=-DEROT
c          endif
c         ENDIF
c         IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=0.0
C
c         IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0           ! not necessary
c         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
c          TEMRO=TEINF
c          TEMPL=TESUP
c          TEMPG=TGAUAT
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c          ALP11=-1.00/AAMAY
c          ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c          ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c          ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c          ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .           ALP44
c          DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c          TSOE2(IPLAT)=-DEROT
c         ENDIF
c         IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=0.0
c        ENDIF                              ! icubib.eq.0
C
c        TSOE1(IPLAT)=TSOE1(IPLAT)*HENER
c        TSOE2(IPLAT)=TSOE2(IPLAT)*HENER                 ! not necessary
C
c       ELSE                                ! f_pc general
c        CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
c    .                IPLAT,HENER,TSOE1,TSOE2,
c    .                TEINF,TESUP,DELEE,IPCMO,    0)
c       ENDIF                               ! ipcmo.eq.0
c      ENDDO                                ! iplat=1,nplat
C
c      DO IPLAT=1,NPLAT
c       SOUR1T=SOUR1T+TSOE2(IPLAT)                       ! not necessary
c       SOUR2T=SOUR2T+TSOE1(IPLAT)
c      ENDDO
C
      ENDIF                                ! ilahes.eq.2
C
 1000 CONTINUE
C
C**** DEALS WITH FILLING MATERIAL
C
      IF(IFILLS.EQ.0) RETURN
C
C**** COMPUTES DENSITY
C

      call runends('ifills=1 not implemented')

c     IF(TGAUST.LE.VDENSFI(1,2)) BASMMFI=VDENSFI(1,1)
c     DO IDENS=2,NDENSFI
c      I1=IDENS-1
c      I2=IDENS
c      IF(TGAUST.GT.VDENSFI(I1,2).AND.TGAUST.LE.VDENSFI(I2,2)) 
c    .  BASMMFI=(VDENSFI(I2,1)-VDENSFI(I1,1))/
c    .          (VDENSFI(I2,2)-VDENSFI(I1,2))*
c    .          (TGAUST-VDENSFI(I1,2))+VDENSFI(I1,1)
c     ENDDO
c     IF(TGAUST.GT.VDENSFI(NDENSFI,2)) BASMMFI=VDENSFI(NDENSFI,1)
C
C**** COMPUTES INITIAL DENSITY
C
c     IF(TGAUIT.LE.VDENSFI(1,2)) BASMIFI=VDENSFI(1,1)
c     DO IDENS=2,NDENSFI
c      I1=IDENS-1
c      I2=IDENS
c      IF(TGAUIT.GT.VDENSFI(I1,2).AND.TGAUIT.LE.VDENSFI(I2,2)) 
c    .  BASMIFI=(VDENSFI(I2,1)-VDENSFI(I1,1))/
c    .          (VDENSFI(I2,2)-VDENSFI(I1,2))*
c    .          (TGAUIT-VDENSFI(I1,2))+VDENSFI(I1,1)
c     ENDDO
c     IF(TGAUIT.GT.VDENSFI(NDENSFI,2)) BASMIFI=VDENSFI(NDENSFI,1)
C
c     IF(ILAHET.EQ.0) THEN
C
C**** ILAHE=0, CALCULATES THE CAPACITY COEFFICIENT AT TIME (t+alpha dt)
C
c      IF(TGAUST.LE.VCAPAFI(1,2)) BASCCFI=VCAPAFI(1,1)
c      DO ICAPA=2,NCAPAFI
c       I1=ICAPA-1
c       I2=ICAPA
c       IF(TGAUST.GT.VCAPAFI(I1,2).AND.TGAUST.LE.VCAPAFI(I2,2)) 
c    .   BASCCFI=(VCAPAFI(I2,1)-VCAPAFI(I1,1))/
c    .           (VCAPAFI(I2,2)-VCAPAFI(I1,2))*
c    .           (TGAUST-VCAPAFI(I1,2))+VCAPAFI(I1,1)
c      ENDDO
c      IF(TGAUST.GT.VCAPAFI(NCAPAFI,2)) BASCCFI=VCAPAFI(NCAPAFI,1)
C
C**** IMPROVED STAGGERED SCHEME FOR THERMOMECHANICAL COUPLED PROBLEMS
C     (isothermal split)
C
c      IF(ITERME.GT.0) THEN
c       IF(ITERMP.GT.0) THEN
c        IF(NITERC.EQ.1.OR.NITERC.EQ.2.OR.NITERC.EQ.4) THEN
c         ISIMP=1                       ! better as input (see plheat.f)
c         IF(ISIMP.EQ.0.OR.ISIMP.EQ.1) THEN   ! simplif. ther. & general
c          BASCCFI=BASCCFI+COUTDT/BASMM
c         ENDIF
c         IF(ISIMP.EQ.2) THEN ! more general but inconsistent with coute
c          BASCCFI=BASCCFI+(TGAUST+CENKEL)*COUTDT/BASMM ! (see plheat.f)
c         ENDIF
c        ENDIF
c       ENDIF
c      ENDIF
C
c     ENDIF       ! ilahet.eq.0
C
c     IF(ILAHET.EQ.1) THEN
C
C**** ILAHE=1, CALCULATE THE LATENT HEAT RELEASED
C
c      SOUR1TFI=0.0
c      SOUR2TFI=0.0
c      DSOURTFI=0.0
c      IF(NPLATFI.EQ.0) GOTO 2000
C
c      TGAUST=TGAUXX
c      TGAUAT=TGAUST-DTEMPT*DTIMET
C
c      DO IPLAT=1,NPLATFI
c       TSOE1FI(IPLAT)=0.0
c       TSOE2FI(IPLAT)=0.0
c       TEINF=VPLATFI(IPLAT,1)
c       TESUP=VPLATFI(IPLAT,2)
c       HENER=VPLATFI(IPLAT,3)
c       IPCFO=INT(VPLATFI(IPLAT,4))
c       IPCMO=INT(VPLATFI(IPLAT,5))
c       DELEE=TESUP-TEINF
C
c       IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN           ! f_pc linear
C
c        IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0
c        IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
c    .    TSOE1FI(IPLAT)=1.0D+00/DELEE*(TGAUST-TEINF)
c        IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D+00
C
c        IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0
c        IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
c    .    TSOE2FI(IPLAT)=1.0D+00/DELEE*(TGAUAT-TEINF)
c        IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D+00
C
c        IF(ICONVT.EQ.1) THEN
c         IF(ICUBIC.EQ.0) THEN    ! fpc cubical (coherent with ILAHET=2)
c          IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0
c          IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
c           TEMRO=TEINF
c           TEMPL=TESUP
c           TEMPG=TGAUST
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c           ALP11=-1.00/AAMAY
c           ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c           ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c           ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c           ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .            ALP44
c           DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c           TSOE1FI(IPLAT)=1.0-ROTET
c          ENDIF
c          IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=1.0D+00
C
c          IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0
c          IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
c           TEMRO=TEINF
c           TEMPL=TESUP
c           TEMPG=TGAUAT
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c           ALP11=-1.00/AAMAY
c           ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c           ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c           ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c           ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .            ALP44
c           DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c           TSOE2FI(IPLAT)=1.0-ROTET
c          ENDIF
c          IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=1.0D+00
c         ENDIF                              ! icubic.eq.0
c        ENDIF                               ! iconvt.eq.1
C
c        TSOE1FI(IPLAT)=TSOE1FI(IPLAT)*HENER
c        TSOE2FI(IPLAT)=TSOE2FI(IPLAT)*HENER
C
c       ELSE                                ! f_pc general
c        CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
c    .                IPLAT,HENER,TSOE1FI,TSOE2FI,
c    .                TEINF,TESUP,DELEE,IPCMO,    1)
c       ENDIF                      ! ipcmo.eq.0
c      ENDDO                       ! iplat=1,nplat
C
c      DO IPLAT=1,NPLATFI
c       SOUR1TFI=SOUR1TFI+TSOE2FI(IPLAT)
c       SOUR2TFI=SOUR2TFI+TSOE1FI(IPLAT)
c      ENDDO
C
C**** PHASE-CHANGE FUNCTION RATE (*DTIMET)
C
c      DSOURTFI=SOUR2TFI-SOUR1TFI
C
c     ENDIF                                         ! ilahet.eq.1
C
c     IF(ILAHET.EQ.2.or.ilahet.eq.3) THEN
C
C**** ILAHE=2, CALCULATE THE TEMPERATURE DERIVATIVE OF THE
C              PHASE-CHANGE FUNCTION
C
C     Note: ICONVT=1 for ILAHET=2,3
C
c      SOUR1TFI=0.0                                      ! not necessary
c      SOUR2TFI=0.0
c      IF(NPLATFI.EQ.0) GOTO 2000
C
c      TGAUST=TGAUXX
c      TGAUAT=TGAUST-DTEMPT*DTIMET
C
c      DO IPLAT=1,NPLATFI
c       TSOE1FI(IPLAT)=0.0
c       TSOE2FI(IPLAT)=0.0                               ! not necessary
c       TEINF=VPLATFI(IPLAT,1)
c       TESUP=VPLATFI(IPLAT,2)
c       HENER=VPLATFI(IPLAT,3)
c       IPCFO=INT(VPLATFI(IPLAT,4))
c       IPCMO=INT(VPLATFI(IPLAT,5))
c       DELEE=TESUP-TEINF
C
c       IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN     ! f_pc linear (not good!)
C
c        IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0
c        IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
c    .    TSOE1FI(IPLAT)=1.0D+00/DELEE
c        IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=0.0
C
c        IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0          ! not necessary
c        IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
c    .    TSOE2FI(IPLAT)=1.0D+00/DELEE
c        IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=0.0
C
c        icubit=0                ! to control df_pc/dT (better as input)
c        if(ilahet.eq.4) icubit=1           ! could be ilahet.eq.3
c        icubib=icubic
c        if(ilahet.eq.3) icubib=0
C
c        IF(ICUBIB.EQ.1) THEN               ! fpc cubical
c         IF(TGAUST.LE.TEINF) TSOE1FI(IPLAT)=0.0
c         IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP) THEN
c          TEMRO=TEINF
c          TEMPL=TESUP
c          TEMPG=TGAUST
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c          ALP11=-1.00/AAMAY
c          ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c          ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c          ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c          ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .           ALP44
c          DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c          if(icubit.eq.0) then
c           TSOE1FI(IPLAT)=-DEROT
c          else
c           tsoea=-derot
c           if(tsoea.lt.tsoe1(iplat)) TSOE1FI(IPLAT)=-DEROT
c          endif
c         ENDIF
c         IF(TGAUST.GT.TESUP) TSOE1FI(IPLAT)=0.0
C
c         IF(TGAUAT.LE.TEINF) TSOE2FI(IPLAT)=0.0         ! not necessary
c         IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP) THEN
c          TEMRO=TEINF
c          TEMPL=TESUP
c          TEMPG=TGAUAT
c          AAMAY=(TEMPL*TEMPL*TEMPL-TEMRO*TEMRO*TEMRO)-3.00*TEMRO*TEMRO*
c    .        (TEMPL-TEMRO)-1.50*(TEMPL+TEMRO)*(TEMPL*TEMPL-TEMRO*TEMRO-
c    .        2.00*TEMRO*(TEMPL-TEMRO))
c          ALP11=-1.00/AAMAY
c          ALP22=-1.50*(TEMPL+TEMRO)*ALP11
c          ALP33=-3.00*ALP11*TEMRO*TEMRO-2.00*ALP22*TEMRO
c          ALP44=-ALP11*TEMPL*TEMPL*TEMPL-ALP22*TEMPL*TEMPL-ALP33*TEMPL
c          ROTET=ALP11*TEMPG*TEMPG*TEMPG+ALP22*TEMPG*TEMPG+ALP33*TEMPG+
c    .           ALP44
c          DEROT=3.0*ALP11*TEMPG*TEMPG+2.0*ALP22*TEMPG+ALP33
c          TSOE2FI(IPLAT)=-DEROT
c         ENDIF
c         IF(TGAUAT.GT.TESUP) TSOE2FI(IPLAT)=0.0
c        ENDIF                              ! icubib.eq.1
C
c        TSOE1FI(IPLAT)=TSOE1FI(IPLAT)*HENER
c        TSOE2FI(IPLAT)=TSOE2FI(IPLAT)*HENER             ! not necessary
C
c       ELSE                                ! f_pc general
c        CALL MICPHAS(TGAUST,TGAUAT,ILAHET,
c    .                IPLAT,HENER,TSOE1FI,TSOE2FI,
c    .                TEINF,TESUP,DELEE,IPCMO,    1)
c       ENDIF                               ! ipcmo.eq.0
c      ENDDO                                ! iplat=1,nplat
C
c      DO IPLAT=1,NPLATFI
c       SOUR1TFI=SOUR1TFI+TSOE2FI(IPLAT)                 ! not necessary
c       SOUR2TFI=SOUR2TFI+TSOE1FI(IPLAT)
c      ENDDO
C
c     ENDIF                                ! ilahet.eq.2
C
c2000 CONTINUE
C
C**** COMPUTES THE WEIGHTED CONDUCTIVITY COEFFICIENT
C
c     BASMM=PSEUDO*BASMM+(1.0-PSEUDO)*BASMMFI
c     BASMI=PSEUDO*BASMI+(1.0-PSEUDO)*BASMIFI
c     IF(ILAHET.EQ.0) THEN
c      BASCC=PSEUDO*BASCC+(1.0-PSEUDO)*BASCCFI
c     ENDIF
c     IF(ILAHET.EQ.1) THEN
c      SOUR1T=PSEUDO*SOUR1T+(1.0-PSEUDO)*SOUR1TFI
c      SOUR2T=PSEUDO*SOUR2T+(1.0-PSEUDO)*SOUR2TFI
c      DSOURT=PSEUDO*DSOURT+(1.0-PSEUDO)*DSOURTFI
c     ENDIF
c     IF(ILAHET.EQ.2.OR.ILAHET.EQ.3) THEN
c      SOUR1T=PSEUDO*SOUR1T+(1.0-PSEUDO)*SOUR1TFI
c      SOUR2T=PSEUDO*SOUR2T+(1.0-PSEUDO)*SOUR2TFI
c     ENDIF
C
      RETURN
      END
