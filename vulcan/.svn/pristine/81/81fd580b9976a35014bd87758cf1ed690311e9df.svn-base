      SUBROUTINE MAS104T(DVOLUT,PROPST,SHAPET,WSTIFT,EHISTT,ELDIST,
     .                   VELCMT,PREASL,TGAPSL,BOUCHL,VNORLT,DISPTT,
     .                   FPCHLT,LNODST,AUXS1T,AUXS2T)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE MASS MATRIX ( ELEMENT NO. 104 )
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
      DIMENSION DVOLUT(*),        PROPST(*),
     .          SHAPET(NNODNT,*), WSTIFT(*),
     .          EHISTT(NHISTT,*), ELDIST(NDOFCT,*)
      DIMENSION VELCMT(*)
      DIMENSION PREASL(*),        TGAPSL(*)
      DIMENSION BOUCHL(NDOFCT,*), VNORLT(NDIMET,*), DISPTT(NDOFCM,*),
     .          FPCHLT(NFPCH,*)
      DIMENSION LNODST(*)
      DIMENSION AUXS1T(NEVABT,*), AUXS2T(NEVABT,*)
C
      NNOBOT=NNODLT/2
      IF(NOCOLT.EQ.1) NNOBOT=NNODLT                ! non-coincident mesh
      NEVBOT=NNOBOT*NDOFNT
C
C**** LOOP ON INTEGRATION POINTS
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
      IF(NMEMO3.EQ.0) THEN
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
     .    BGAUS2=BGAUS2+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT+NNOBOT)
        ENDIF
        IF(IFILL.EQ.1) THEN
         IF(IMICR.EQ.0) THEN
          IPSEU=2*NNUPT+1
         ELSE
          IPSEU=2*NNUPT+NNUPO+1
         ENDIF
         PSEUDO1=PSEUDO1+SHAPET(INODLT,IGAUST)*FPCHLT(IPSEU,INODLT)
         IF(NOCOLT.EQ.0)
     .    PSEUDO2=PSEUDO2+SHAPET(INODLT,IGAUST)*
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
       IF(IGABO2.EQ.0) THEN                      ! AIRGAP
        IDECIT=0
        IF(ITERME.GT.0) THEN                     ! bidirectional coupled
         IF(ITERMG.GT.0) IDECIT=1                ! gap dependency
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
     .             DVOLUT(IGAUST)
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
     .                     VNORLT(IDOFN,IGAUST)
            ELSE
             ELDI2=DISPTT(IDOFN,INODLT+NNOBOT+IAUXYT)
             GAPTET=GAPTET+(SHAPET(INODLT,IGAUST)*ELDI1-
     .                      SHAPET(INODLT+NNOBOT,IGAUST)*ELDI2)*
     .                                              VNORLT(IDOFN,IGAUST)
            ENDIF
            IF(IR432.EQ.1) THEN
             GAPTEP=GAPTEP+SHAPET(INODLT,IGAUST)*(ELDI1-ELDI2)*
     .                     VNORLT(IDOFN,IGAUST)/DVOLUT(IGAUST)
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
         CALL RADCOFU(COEFHT,PROPST,TGAUS1,TGAUS2,GAPTET,
     .                PRESUT,PSEUDO1,PSEUDO2,     1)
C
        ELSE                         ! thermal or unidirectional coupled
C
C**** UP DATE THE CONVECTION-RADIATION COEFFICIENT
C     (ONLY TEMPERATURE DEPENDENT)
C
         CALL RADCOFT(COEFHT,PROPST,TGAUS1,TGAUS2,PSEUDO1,PSEUDO2,    1)
C
        ENDIF                  ! idecit.eq.1
       ENDIF                   ! igabo2.eq.0
C
C**** EVALUATES IF AIRGAP BECOMES BOUNDARY DUE TO CRITICAL GAP
C
       IF(NOCOLT.EQ.1) THEN
        IF(IGABO2.EQ.0) THEN
         IDECOT=0
         IF(ITERME.GT.0) THEN                    ! bidirectional coupled
          IF(ITERMD.GT.0) IDECOT=1               ! deformed shape
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
     .                         VNORLT(IDOFN,IGAUST)
                ELSE
                 ELDI2=DISPTT(IDOFN,INODLT+NNOBOT+IAUXYT)
                 GAPTET=GAPTET+(SHAPET(INODLT,IGAUST)*ELDI1-
     .                          SHAPET(INODLT+NNOBOT,IGAUST)*ELDI2)*
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
             CALL RADCOFT(EHISTT(1,IGAUST),PROPST,TGAUS1,TGAUS2,
     .                                            PSEUDO1,PSEUDO2,    3)
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
      ELSE
       BGAUS1=0.0D0
       BGAUS2=0.0D0
       IF(NMEMO10.EQ.1) THEN
        DO INODLT=1,NNODLT
         BGAUS1=BGAUS1+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT)
         IF(NOCOLT.EQ.0)
     .    BGAUS2=BGAUS2+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT+NNOBOT)
        END DO
       ENDIF
       COEFHT=EHISTT(1,IGAUST)
      ENDIF                    ! nmemo3.eq.0
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
      COEFHT=COEFHT*(FACTR1+FACTR2)/2.0D0      ! approximation (symm. J)
C
C**** ELEMENTAL MASS MATRIX
C
      IF(NOCOLT.EQ.0) THEN
       CALL WBOUNDT(DVOLUT(IGAUST),NDOFNT,NEVBOT,NNOBOT,PROPST,
     .              SHAPET(1,IGAUST),WSTIFT,COEFHT,
     .              ELDIST,NDOFCT,KSYMMT,AUXS1T,  104)
C
       CALL WINTERT(DVOLUT(IGAUST),NEVABT,NNOBOT,WSTIFT,KSYMMT,
     .              AUXS1T,AUXS2T)
      ELSE
       CALL WBOUNCT(DVOLUT(IGAUST),NDOFNT,NEVBOT,NNOBOT,PROPST,
     .              SHAPET(1,IGAUST),WSTIFT,COEFHT,
     .              ELDIST,NDOFCT,KSYMMT,AUXS1T,NEVABT,AUXS2T,
     .              IAUXYT,LNODST,IGABO2)
      ENDIF
C
  100 CONTINUE
C
      RETURN
      END
