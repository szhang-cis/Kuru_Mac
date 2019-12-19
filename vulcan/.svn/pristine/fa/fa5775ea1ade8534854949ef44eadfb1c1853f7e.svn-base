      SUBROUTINE FRI104T(PROPST,ELDIST,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   DISPTT,ELCODT,PREASL,TGAPSL,ELELTT,VELCMT,
     .                   BOUCHL,VNORLT,FPCHLT,LNODST)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL FORCES" FOR BOUNDARY
C     ELEMENTS
C
C     Notes: IR432 is an index indicating the type of contact element
C            used in the mechanical problem (IR432=1 for ITYPE=4 & 
C            IR432=2 for ITYPE=32) that is useful to evaluate the
C            normal pressure to be used in the defintion of the heat
C            transfer value. 
C
C            SHAPET, ELDIST & VELCMT must be dimensioned
C            with NNODNT due to non-coincident contact mesh (NOCOLT=1)
C
C
C     Index of variables:
C
C     EHISTT(   1) = Heat transfer coefficient
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
     .          SHAPET(NNODNT,*), EHISTT(NHISTT,*)
      DIMENSION DISPTT(NDOFCM,*), ELCODT(NDIMET,NNODLT)
      DIMENSION PREASL(*),        TGAPSL(*),
     .          ELELTT(*),        VELCMT(*)
      DIMENSION BOUCHL(NDOFCT,*), VNORLT(NDIMET,*),
     .          FPCHLT(NFPCH,*),  LNODST(*)
C
C**** PREASL: NORMAL FORCE (NMEMO11=0)
C     TGAPSL: NORMAL GAP   (NMEMO11=0)
C
C     DISPTT: DISPLACEMENTS (USED TO COMPUTE THE NORMAL PRESSURE & GAP
C             WHEN NMEMO11=1)
C
      NNOBOT=NNODLT/2
      IF(NOCOLT.EQ.1) NNOBOT=NNODLT                ! non-coincident mesh
C
C**** INITIALIZE ELEMENT CONTRIBUTION OF ELELTT (K*T)
C
      DO IEVABT=1,NNODLT
       ELELTT(IEVABT)=0.0D00
      END DO
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO 100 IGAUST=1,NGAULT
C
      IGABO2=0                                 ! AIRGAP
      IF(NOCOLT.EQ.1) THEN                     ! checks contact
       IAUXYT=(IGAUST-1)*NNODLT
       IF(LNODST(NNOBOT).EQ.LNODST(NNODLT+NNOBOT+IAUXYT)) THEN
        IF(IGABO3.EQ.0) THEN
         GO TO 100
        ELSE
         IGABO2=1                              ! AIRGAP becomes BOUNDARY
        ENDIF
       ENDIF
      ENDIF
C
C**** EVALUATE TEMPERATURE IN BOTH BODIES
C
      DTEMP1=0.0D0
      DTEMP2=0.0D0
      TGAUS1=0.0D0
      TGAUS2=0.0D0
      BGAUS1=0.0D0
      BGAUS2=0.0D0
      PSEUDO1=0.0D0
      PSEUDO2=0.0D0
      DO INODLT=1,NNOBOT
       IF(KDYNAT.EQ.1) THEN
        DTEMP1=DTEMP1+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
        IF(NOCOLT.EQ.0) THEN
         DTEMP2=DTEMP2+SHAPET(INODLT,IGAUST)*VELCMT(INODLT+NNOBOT)
        ELSE
         IF(IGABO2.EQ.0) THEN
          DTEMP2=DTEMP2+SHAPET(INODLT+NNOBOT,IGAUST)*
     .                                      VELCMT(INODLT+NNOBOT+IAUXYT)
         ELSE
          DTEMP2=0.0D0
         ENDIF
        ENDIF
       ENDIF
       TGAUS1=TGAUS1+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
       IF(NOCOLT.EQ.0) THEN
        TGAUS2=TGAUS2+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT+NNOBOT)
       ELSE
        IF(IGABO2.EQ.0) THEN
         TGAUS2=TGAUS2+SHAPET(INODLT+NNOBOT,IGAUST)*
     .                                    ELDIST(1,INODLT+NNOBOT+IAUXYT)
        ELSE
         TGAUS2=TENVI3
        ENDIF
       ENDIF
       IF(NMEMO10.EQ.1) THEN
        BGAUS1=BGAUS1+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT)
        IF(NOCOLT.EQ.0)
     .   BGAUS2=BGAUS2+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT+NNOBOT)
       ENDIF
       IF(IFILL.EQ.1) THEN
        IF(IMICR.EQ.0) THEN
         IPSEU=2*NNUPT+1
        ELSE
         IPSEU=2*NNUPT+NNUPO+1
        ENDIF
        PSEUDO1=PSEUDO1+SHAPET(INODLT,IGAUST)*FPCHLT(IPSEU,INODLT)
        IF(NOCOLT.EQ.0)
     .   PSEUDO2=PSEUDO2+SHAPET(INODLT,IGAUST)*
     .                                       FPCHLT(IPSEU,INODLT+NNOBOT)
       ENDIF
      END DO
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
      IF(KDYNAT.EQ.1) THEN
       IF(KINTET.EQ.1) THEN           ! Euler method
        TGAUS1=TALFAT*TGAUS1+(1.0D0-TALFAT)*(TGAUS1-DTEMP1*DTIMET)
        TGAUS2=TALFAT*TGAUS2+(1.0D0-TALFAT)*(TGAUS2-DTEMP2*DTIMET)
       ENDIF
      ENDIF
C
      IF(IGABO2.EQ.0) THEN                       ! AIRGAP
       IDECIT=0
       IF(ITERME.GT.0) THEN                      ! bidirectional coupled
        IF(ITERMG.GT.0) IDECIT=1                 ! gap dependency
       ENDIF
       IF(IDECIT.EQ.1) THEN
C
C**** EVALUATE NORMAL GAP & NORMAL PRESSURE
C
        GAPTET=0.0D+00
        PRESUT=0.0D+00
        IF(IR432.EQ.1) GAPTEP=0.0D0
        IF(NMEMO11.EQ.0) THEN
         DO INODLT=1,NNOBOT
          GAPTET=GAPTET+SHAPET(INODLT,IGAUST)*TGAPSL(INODLT)
          IF(IR432.EQ.1) THEN
           PRESUT=PRESUT+SHAPET(INODLT,IGAUST)*PREASL(INODLT)/
     .            DVOLUT(IGAUST)
          ELSE
           PRESUT=PRESUT+SHAPET(INODLT,IGAUST)*PREASL(INODLT)
          ENDIF
         END DO
        ELSE
         DO IDOFN=1,NDOFCM          ! NDOFN=NDOFCM
          DO INODLT=1,NNOBOT
           ELDI1=DISPTT(IDOFN,INODLT)
           IF(NOCOLT.EQ.0) THEN
            ELDI2=DISPTT(IDOFN,INODLT+NNOBOT)
            GAPTET=GAPTET+SHAPET(INODLT,IGAUST)*(ELDI1-ELDI2)*
     .                    VNORLT(IDOFN,IGAUST)
           ELSE
            ELDI2=DISPTT(IDOFN,INODLT+NNOBOT+IAUXYT)
            GAPTET=GAPTET+(SHAPET(INODLT,IGAUST)*ELDI1-
     .                     SHAPET(INODLT+NNOBOT,IGAUST)*ELDI2)*
     .                                              VNORLT(IDOFN,IGAUST)
           ENDIF
           IF(IR432.EQ.1) THEN
            GAPTEP=GAPTEP+SHAPET(INODLT,IGAUST)*(ELDI1-ELDI2)*
     .                    VNORLT(IDOFN,IGAUST)/DVOLUT(IGAUST)
           ENDIF
          END DO
         END DO
         IF(GAPTET.GE.0.0D0) THEN   ! absolute values
          IF(IR432.EQ.1) THEN
           PRESUT=RIGINT*GAPTEP
          ELSE
           PRESUT=RIGINT*GAPTET
          ENDIF
          GAPTET=0.0D0
         ELSE
          PRESUT=0.0D0
          GAPTET=-GAPTET
         ENDIF
        ENDIF                       ! nmemo11.eq.0
C
C**** UP DATE THE CONVECTION-RADIATION COEFFICIENT
C     (TEMPERATURE & GAP DEPENDENT)
C
        IF(NMEMO3.EQ.0) THEN
         CALL RADCOFU(COEFHT,PROPST,TGAUS1,TGAUS2,GAPTET,
     .                PRESUT,PSEUDO1,PSEUDO2,     1)
        ELSE
         CALL RADCOFU(EHISTT(1,IGAUST),PROPST,TGAUS1,TGAUS2,GAPTET,
     .                PRESUT,PSEUDO1,PSEUDO2,     1)
         COEFHT=EHISTT(1,IGAUST)
        ENDIF
C
       ELSE                          ! thermal or unidirectional coupled
C
C**** UP DATE THE CONVECTION-RADIATION COEFFICIENT
C     (ONLY TEMPERATURE DEPENDENT)
C
        IF(NMEMO3.EQ.0) THEN
         CALL RADCOFT(COEFHT,PROPST,TGAUS1,TGAUS2,PSEUDO1,PSEUDO2,    1)
        ELSE
         CALL RADCOFT(EHISTT(1,IGAUST),PROPST,TGAUS1,TGAUS2,PSEUDO1,
     .                                                    PSEUDO2,    1)
         COEFHT=EHISTT(1,IGAUST)
        ENDIF
C
       ENDIF              ! idecit.eq.1
      ENDIF               ! igabo2.eq.0
C
C**** EVALUATES IF AIRGAP BECOMES BOUNDARY DUE TO CRITICAL GAP
C
      IF(NOCOLT.EQ.1) THEN
       IF(IGABO2.EQ.0) THEN
        IDECOT=0
        IF(ITERME.GT.0) THEN                     ! bidirectional coupled
         IF(ITERMD.GT.0) IDECOT=1                ! deformed shape
        ENDIF
        IF(IDECOT.EQ.1) THEN
         IF(IGABO3.EQ.1) THEN
          IF(IGABO4.EQ.1) THEN
           IF(IDECIT.EQ.0) THEN
C
C**** EVALUATE NORMAL GAP
C
            GAPTET=0.0D+00
            IF(NMEMO11.EQ.0) THEN
             DO INODLT=1,NNOBOT
              GAPTET=GAPTET+SHAPET(INODLT,IGAUST)*TGAPSL(INODLT)
             END DO
            ELSE
             DO IDOFN=1,NDOFCM          ! NDOFN=NDOFCM
              DO INODLT=1,NNOBOT
               ELDI1=DISPTT(IDOFN,INODLT)
               IF(NOCOLT.EQ.0) THEN
                ELDI2=DISPTT(IDOFN,INODLT+NNOBOT)
                GAPTET=GAPTET+SHAPET(INODLT,IGAUST)*(ELDI1-ELDI2)*
     .                        VNORLT(IDOFN,IGAUST)
               ELSE
                ELDI2=DISPTT(IDOFN,INODLT+NNOBOT+IAUXYT)
                GAPTET=GAPTET+(SHAPET(INODLT,IGAUST)*ELDI1-
     .                         SHAPET(INODLT+NNOBOT,IGAUST)*ELDI2)*
     .                                              VNORLT(IDOFN,IGAUST)
               ENDIF
              END DO
             END DO
             IF(GAPTET.GE.0.0D0) THEN   ! absolute values
              GAPTET=0.0D0
             ELSE
              GAPTET=-GAPTET
             ENDIF
            ENDIF                       ! nmemo11.eq.0
           ENDIF                        ! idecit.eq.0
           IF(GAPTET.GT.AGABO4) THEN    ! AIRGAP becomes BOUNDARY
            IGABO2=2
            DTEMP2=0.0D0
            TGAUS2=TENVI3
C
C**** UP DATE THE CONVECTION-RADIATION COEFFICIENT (AS BOUNDARY)
C     (ONLY TEMPERATURE DEPENDENT)
C
            IF(NMEMO3.EQ.0) THEN
             CALL RADCOFT(COEFHT,PROPST,TGAUS1,TGAUS2,PSEUDO1,PSEUDO2,
     .                                                                3)
            ELSE
             CALL RADCOFT(EHISTT(1,IGAUST),PROPST,TGAUS1,TGAUS2,PSEUDO1,
     .                                                    PSEUDO2,    3)
             COEFHT=EHISTT(1,IGAUST)
            ENDIF
           ENDIF                       ! gaptet.gt.agabo4
          ENDIF                        ! igabo4.eq.1
         ENDIF                         ! igabo3.eq.1
        ENDIF                          ! idecot.eq.1
       ENDIF                           ! igabo2.eq.0
      ENDIF                            ! nocoit.eq.1
C
C**** UP DATE THE CONVECTION-RADIATION COEFFICIENT (AS BOUNDARY)
C     (ONLY TEMPERATURE DEPENDENT)
C
      IF(NOCOLT.EQ.1) THEN
       IF(IGABO2.EQ.1) THEN
        IF(NMEMO3.EQ.0) THEN
         CALL RADCOFT(COEFHT,PROPST,TGAUS1,TGAUS2,PSEUDO1,PSEUDO2,
     .                                                                3)
        ELSE
         CALL RADCOFT(EHISTT(1,IGAUST),PROPST,TGAUS1,TGAUS2,PSEUDO1,
     .                                                    PSEUDO2,    3)
         COEFHT=EHISTT(1,IGAUST)
        ENDIF
       ENDIF                           ! igabo2.eq.0
      ENDIF                            ! nocoit.eq.1
C
      FACTR1=1.0D0
      FACTR2=1.0D0
      IF(NMEMO10.EQ.1) THEN
       IF(NDIMET.EQ.1) THEN
        FACTR1=BGAUS1
        IF(NOCOLT.EQ.0)
     .   FACTR2=BGAUS2
       ELSE
        FACTR1=BGAUS1**(-(NDIMET-1)/NDIMET)
        IF(NOCOLT.EQ.0)
     .   FACTR2=BGAUS2**(-(NDIMET-1)/NDIMET)
       ENDIF
      ENDIF
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling
       IF(ITERMD.GT.0) THEN                     ! deformed shape
        IF(LARGET.EQ.1) THEN                    ! TLF
         detjj1=1.0D0   ! determinant of F (from mechanical computation)
         FACTR1=DETJJ1
         IF(NOCOLT.EQ.0) THEN
          detjj2=1.0D0  ! determinant of F (from mechanical computation)
          FACTR2=DETJJ2
         ENDIF
         COEFHT=COEFHT
        ENDIF
        IF(LARGET.EQ.2) THEN                    ! ULF
         FACTR1=1.0D0
         IF(NOCOLT.EQ.0)
     .    FACTR2=1.0D0
         COEFHT=COEFHT
        ENDIF
       ENDIF
      ENDIF
C
      COEFH1=COEFHT*FACTR1
      COEFH2=COEFHT*FACTR2
C
C**** CALCULATE THE EFFECTIVE HEATS IN BOTH BODIES
C
      DESI1T=COEFH1*TGAUS1
      DESI2T=COEFH2*TGAUS2
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEATS FORCES" DUE TO 
C     CONVECTION-RADIATION EFFECTS (IN BOTH BODIES FOR NOCOLT=0 OR IN
C     ONLY ONE BODY FOR NOCOLT=1)
C     (AS A EXTERNAL HEAT FLUX >> -, AS A RESIDUAL HEAT FLUX >> +)
C
       DO IEVABT=1,NNOBOT
        ELELMT(IEVABT)=ELELMT(IEVABT)+SHAPET(IEVABT,IGAUST)*
     .                                    (DESI1T-DESI2T)*DVOLUT(IGAUST)
       END DO
C
  100 CONTINUE
C
      IF(NOCOLT.EQ.0) THEN                    ! only for coincident mesh
       DO IEVABT=1,NNOBOT
        ELELMT(IEVABT+NNOBOT)=-ELELMT(IEVABT)
       END DO
      ENDIF
C
C**** CONTRIBUTION TO ELELTT (K*T IN CONVERGENCE CRITERION)
C
      DO IEVABT=1,NNOBOT
       ELELTT(IEVABT)=ELELMT(IEVABT)
       IF(NOCOLT.EQ.0) ELELTT(IEVABT+NNOBOT)=ELELMT(IEVABT+NNOBOT)
      END DO
C
      RETURN
      END
