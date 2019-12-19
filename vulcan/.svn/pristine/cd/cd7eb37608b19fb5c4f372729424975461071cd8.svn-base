      SUBROUTINE MICROS8(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ALUMINIUM NITRIDE PRECIPITATION
C     WITH THE MICROSTRUCTURAL MODEL NUMBER 8 (IPCMO=8) OF RATE
C     PHASE-CHANGE FORMULATIONS
C
C     ALUMINIUM NITRIDE PRECIPITATION & RECRYSTALIZATION
C
C***********************************************************************
C
C     Index of variables
C
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     DTEMPT= Temperature rate
C
C     TSOE2 = L*Phase-change function at time t (=0 in this case)
C     TSOE1 = L*Phase-change function at time t+dt (=0 in this case)
C
C     ALPHAM= array of microstructural (microscopical) variables
C             First indexes of ALPHAM devoted to f_pc functions
C             (see outsmot.f)
C
C     IN=INUPC (=0 in this case)
C     ALPHAM(IN+1): aluminium nitride (AlN) fraction (at time t+dt)
C     ALPHAM(IN+2): parameter B (at time t+dt)
C     ALPHAM(IN+3): recrystalization fraction (at time t+dt)
C     ALPHAM(IN+4): parameter Br (at time t+dt)
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
      DIMENSION BASKK(NDIMETO)
C
      DIMENSION TSOE1(5), TSOE2(5), TSOC1(5)
C
      HENER=0.0D0                         ! no phase-change in this case
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
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES (LOAD LAST CONVERGED
C     VALUES)
C
      IN=INUPC
      XNO= ALPHAM(IN+1)
      BBO= ALPHAM(IN+2)
      XNRO=ALPHAM(IN+3)
      BBRO=ALPHAM(IN+4)
C
C**** SOLVE MODEL 8 (INTEGRATION OF ALUMINIUM NITRIDE PRECIPITATION &
C     RECRYSTALIZATION EQUATIONS)
C
      XNN= XNO                         ! initial guess
      BBN= BBO
      XNRN=XNRO
      BBRN=BBRO
C
      TEMPE=TGAUST+AP
      DBB=1.0D0/AA*EXP(-(AQ/(AR*TEMPE)))*DTIMET
      BBN=BBO+DBB
      XNN=1.0D0-EXP(-(BBN**AK))                      ! AlN precipitation
C
      IF(IAQR.EQ.1) THEN
       AFNI=AQRN*(1.0D0-XNN)                         ! free nitrogen
       AQR=AQRA*AFNI*AFNI+AQRB*AFNI+AQRC             ! variable act. en.
      ENDIF
      DBBR=1.0D0/AAR*EXP(-(AQR/(AR*TEMPE)))*AKR*TTIMET**(AKR-1.0D0)*
     .                                                            DTIMET
      BBRN=BBRO+DBBR
      XNRN=1.0D0-EXP(-BBRN)                          ! recrystalization
C
      TSOE1(IPLAT)=0.0D0              ! f_pc at time t+dt
      TSOE2(IPLAT)=0.0D0              ! f_pc at time t
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENSITY, CAPACITY,
C     CONDUCTIVITY, PHASE-CHANGE FUNCTION AND LATENT HEAT)
C     => NOT NECESSARY IN THIS CASE
C

C
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER/BASMM   ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER/BASMM   ! L*f_pc at time t+dt
C
C**** OPTIONS IN THE EVALUATION OF df_pc/dT
C     => NOT NECESSARY IN THIS CASE
C

C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=XNN
      ALPHAM(IN+2)=BBN
      ALPHAM(IN+3)=XNRN
      ALPHAM(IN+4)=BBRN
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+4
C
      RETURN
      END
