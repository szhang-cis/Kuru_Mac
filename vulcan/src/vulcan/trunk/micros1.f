      SUBROUTINE MICROS1(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM,INUPC)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH THE MICROSTRUCTURAL MODEL
C     NUMBER 1 (IPCMO=1) OF RATE PHASE-CHANGE FORMULATIONS
C
C***********************************************************************
C
C     Index of variables
C
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     TGAUIT= Initial temperature
C     DTEMPT= Temperature rate
C
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C
C     PCFU1 = Phase-change function at time t
C     PCFU2 = Phase-change function at time t+dt
C     ALATEN= Corresponding specific latent heat associated to
C             the phase-change function
C
C     ALPHAM= array of microstructural (microscopical) variables
C     ALPHAM= array of microstructural (microscopical) variables
C
C     IN=INUPC (INUPC: number of microscopical phase-changes)
C     ALPHAM(IN+1): solid-solid phase-change function (at time t+dt)
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
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      TEINF=VPLAT(IPLAT,1)
      TESUP=VPLAT(IPLAT,2)
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF
C
C**** TRANSFERS "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      IN=INUPC
      FPCSO=ALPHAM(IN+1)
C
C**** LOAD LAST CONVERGED VALUES
C

C
C**** SOLVE MODEL 1
C
      IF(TGAUST.LE.TEINF) TSOE1(IPLAT)=0.0
      IF(TGAUST.GT.TEINF.AND.TGAUST.LE.TESUP)
     . TSOE1(IPLAT)=1.0/DELEE*(TGAUST-TEINF)
      IF(TGAUST.GT.TESUP) TSOE1(IPLAT)=1.0
C

c     if(TSOE1(IPLAT).lt.fpcso) tsoe1(iplat)=fpcso   ! non-reversible

      FPCSS=TSOE1(IPLAT)
C
      IF(TGAUAT.LE.TEINF) TSOE2(IPLAT)=0.0
      IF(TGAUAT.GT.TEINF.AND.TGAUAT.LE.TESUP)
     . TSOE2(IPLAT)=1.0/DELEE*(TGAUAT-TEINF)
      IF(TGAUAT.GT.TESUP) TSOE2(IPLAT)=1.0
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENSITY, CAPACITY,
C     CONDUCTIVITY, PHASE-CHANGE FUNCTION AND LATENT HEAT)
C

C
      TSOE2(IPLAT)=TSOE2(IPLAT)*HENER      ! L*f_pc at time t
      TSOE1(IPLAT)=TSOE1(IPLAT)*HENER      ! L*f_pc at time t+dt
C
      ICACOT=1    ! problem: ICACOT is not defined for each phase-change
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
       IF(IITERT.GT.0)
     .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FPCSS
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+1
C
      RETURN
      END
