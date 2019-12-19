      SUBROUTINE FRI101T(PROPST,ELDIST,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELELTT,VELCMT,BOUCHL,FPCHLT,GPCODT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL FORCES" FOR BOUNDARY
C     ELEMENTS
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
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
      DIMENSION PROPST(*),        ELDIST(NDOFCT,*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELELTT(*),        VELCMT(*)
      DIMENSION BOUCHL(NDOFCT,*), FPCHLT(NFPCH,*),
     .          GPCODT(NDIMET,*)
C
C**** INITIALIZE ELEMENT CONTRIBUTION OF ELELTT (K*T)
C
      DO IEVABT=1,NNODLT
       ELELTT(IEVABT)=0.0D0
      END DO
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUST=1,NGAULT
C
C**** EVALUATE TEMPERATURE
C    
      DTEMPT=0.0D0
      TGAUST=0.0D0
      BGAUST=0.0D0
      PSEUDO=0.0D0
      DO INODLT=1,NNODLT
       IF(KDYNAT.EQ.1)
     . DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
       TGAUST=TGAUST+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
       IF(NMEMO10.EQ.1)
     .  BGAUST=BGAUST+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT)
       IF(IFILL.EQ.1) THEN
        IF(IMICR.EQ.0) THEN
         IPSEU=2*NNUPT+1
        ELSE
         IPSEU=2*NNUPT+NNUPO+1
        ENDIF
        PSEUDO=PSEUDO+SHAPET(INODLT,IGAUST)*FPCHLT(IPSEU,INODLT)
       ENDIF
      END DO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
      IF(KDYNAT.EQ.1) THEN
       IF(KINTET.EQ.1)       ! Euler method
     .  TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
      ENDIF
C
C**** UPDATES THE CONVECTION-RADIATION COEFFICIENT
C
      CALL ALTCONT(TGAUSX,TGAUST,GPCODT(1,IGAUST))
      IF(NMEMO3.EQ.0) THEN
       CALL RADCOFT(COEFHT,PROPST,TGAUSX,TGAUSX,PSEUDO,PSEUDO,    2)
      ELSE
       CALL RADCOFT(EHISTT(1,IGAUST),PROPST,TGAUSX,TGAUSX,PSEUDO,
     .                                                 PSEUDO,    2)
       COEFHT=EHISTT(1,IGAUST)
      ENDIF
C
      FACTRO=1.0D0
      IF(NMEMO10.EQ.1) THEN
       IF(NDIMET.EQ.1) THEN
        FACTRO=BGAUST
       ELSE
        FACTRO=BGAUST**(-(NDIMET-1)/NDIMET)
       ENDIF
      ENDIF
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling
       IF(ITERMD.GT.0) THEN                     ! deformed shape
        IF(LARGET.EQ.1) THEN                    ! TLF
         detjj=1.0D0    ! determinant of F (from mechanical computation)
         FACTRO=DETJJ
         COEFHT=COEFHT
        ENDIF
        IF(LARGET.EQ.2) THEN                    ! ULF
         FACTRO=1.0D0
         COEFHT=COEFHT
        ENDIF
       ENDIF
      ENDIF
C
      COEFHT=COEFHT*FACTRO
C
C**** CALCULATE THE EFFECTIVE HEATS 
C
      TENVIT=PROPST(2)                          ! IFILL=PROPST(1)
      DESIHT=COEFHT*(TGAUST-TENVIT)
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEATS FORCES" DUE TO
C     CONVECTION-RADIATION EFFECTS
C     (AS A EXTERNAL HEAT FLUX >> -, AS A RESIDUAL HEAT FLUX >> +)
C
      DO IEVABT=1,NNODLT
       ELELMT(IEVABT)=ELELMT(IEVABT)+SHAPET(IEVABT,IGAUST)*DESIHT*
     .                DVOLUT(IGAUST)
      END DO 
C
  100 CONTINUE
C
C**** CONTRIBUTION TO ELELTT
C
      DO IEVABT=1,NNODLT
       ELELTT(IEVABT)=ELELMT(IEVABT)
      END DO 
C
      RETURN
      END
