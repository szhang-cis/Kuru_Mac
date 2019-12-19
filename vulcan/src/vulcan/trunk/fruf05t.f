      SUBROUTINE FRUF05T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,ELEL1T,ELEL2T,
     .                   WHAPEL,INPCCT,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE FOR UNIPHASED ELEMENTS
C     ( ELEMENT  NO 5 )
C    
C***********************************************************************
C
C     Index of variables:
C
C     EHISTT(   1) = Density
C     EHISTT(   2) = Specific Heat coefficient
C     EHISTT(   3) = Isotropic Conductivity or Conduct. x (orthot. mat.)
C     EHISTT(   4) = Conductivity y (orthotropic material)
C     EHISTT(   5) = Conductivity z (orthotropic material)
C     EHISTT(3:11) = Conductivity for the fully anisotropc mat. (3D)
C     EHISTT(4+IX) = L*Phase-change function
C     EHISTT(5+IX) = L*Phase-change function rate
C     EHISTT(6+IX) = Initial density
C     EHISTT(7+IX) = Coupling coefficient
C     EHISTT(8+IX) = Temperature derivative of phase-change function
C
C     SHAPE=   SHAPE WHEN INTEGRATING AN UNIPHASE ELEMENT
C              SHAPD WHEN INTEGRATING ONE PHASE OF A MULTIPHASE ELEMENT
C     DVOLU=   DVOLU WHEN INTEGRATING AN UNIPHASE ELEMENT
C              DVOLD WHEN INTEGRATING ONE PHASE OF A MULTIPHASE ELEMENT
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES (thermal-microstructural)
C
      INCLUDE 'nuec_om.f'             ! thermal-mechanical
      INCLUDE 'nued_om.f'             ! thermal-microstructural
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/MULPHASE/MULPHT
C
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), DISIMT(*),
     .          ELELTT(*)
      DIMENSION ELEL1T(*),        ELEL2T(*)
      DIMENSION WHAPEL(NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),  DVOLIT(*)
C
C**** INTEGRATES THE INTERNAL ENERGY IN TIME T
C
      IF(INPCCT.EQ.1) THEN
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO IGAUST=1,NGAULT
C
        DDDTPT=0.0D0
        DTEMPT=0.0D0
        TGAUST=0.0D0
        TGAUIT=0.0D0
        PSEUDO=0.0D0
        DO INODLT=1,NNODLT
         DDDTPT=DDDTPT+SHAPET(INODLT,IGAUST)*DISIMT(INODLT)
         DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
         TGAUST=TGAUST+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
         TGAUIT=TGAUIT+SHAPET(INODLT,IGAUST)*TEINIT(1,INODLT)
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
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
        TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
        IF(KDYNAT.EQ.1) THEN
         IF(KINTET.EQ.1)       ! Euler method
     .    TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
        ENDIF
C
C**** UPDATES THE PHASE-CHANGE FUNCTION
C
        ILAHET=1
        ISINRT=2
C
        IF(NMEMO3.EQ.0) THEN
         CALL CAPCOFT(BASMM, BASCC, PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ELSE
         IX=NDIMETO-1
         IF(IMICR.EQ.1) THEN
          BASMM=EHISTT(1,IGAUST)             ! all computed in micr05t.f
          BASMI=EHISTT(6+IX,IGAUST)
          SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
          SOUR2T=EHISTT(4+IX,IGAUST)
         ELSE
          IF(MULPHT.EQ.1) THEN
           CALL CAPCOFT(BASMM, BASCC, PROPST,TGAUST,TGAUSX,
     .                  DTEMPT,SOUR1T,EHISTT(4+IX,IGAUST),ILAHET,ISINRT,
     .                  EHISTT(5+IX,IGAUST),COUTDT,BASMI,TGAUIT,PSEUDO)
           BASMM=EHISTT(1,IGAUST)                ! computed in frdy05t.f
           BASMI=EHISTT(6+IX,IGAUST)             ! computed in frdy05t.f
           SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
           SOUR2T=EHISTT(4+IX,IGAUST)
          ELSE
           CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                  DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                  DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
          ENDIF
         ENDIF
        ENDIF
C
C**** EVALUATES EFFECTIVE HEATS
C
        CONSTT=1.0D0
ctm     IF(KINTET.EQ.4) CONSTT=TALFAT
C
        IF(LINPC.EQ.1) SOUR1T=SOUR1T*DTEMPT*DTIMET
        DESILT=BASMI*SOUR1T/(DTIMET*CONSTT)
C
        IF(ITERME.GT.0) THEN                    ! bidirectional coupling
         IF(ITERMD.GT.0) THEN                   ! deformed shape
          IF(LARGET.EQ.1) THEN                  ! TLF
           DESILT=BASMI*SOUR1T/(DTIMET*CONSTT)
          ENDIF
          IF(LARGET.EQ.2) THEN                  ! ULF
           detjj=1.0d0  ! determinant of F (from mechanical computation)
           BASMM=BASMI/DETJJ
           DESILT=BASMM*SOUR1T/(DTIMET*CONSTT)
          ENDIF
          IF(LARGET.EQ.3) THEN                  ! Eulerian
c
c computation of dvolit for element subdivision should be implemented !!
c
c          BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c          DESILT=BASMM*SOUR1T/(DTIMET*CONSTT)
          ENDIF
         ENDIF
        ENDIF
C      
C**** INTEGRATES THE HEATS INTO THE INTERNAL "HEAT FORCES" DUE TO
C     TRANSIENT EFFECTS
C
        DO IEVABT=1,NNODLT
         ELEL1T(IEVABT)=ELEL1T(IEVABT)+WHAPEL(IEVABT,IGAUST)*DESILT*
     .                  DVOLUT(IGAUST)
        END DO 
C
       END DO
C
C**** INTEGRATES THE INTERNAL ENERGY IN TIME T+DT
C
      ELSE
C
C**** LOOP OVER INTEGRATION POINTS
C
       DO IGAUST=1,NGAULT
C
        DDDTPT=0.0D0
        DTEMPT=0.0D0
        TGAUST=0.0D0
        TGAUIT=0.0D0
        PSEUDO=0.0D0
        DO INODLT=1,NNODLT
         DDDTPT=DDDTPT+SHAPET(INODLT,IGAUST)*DISIMT(INODLT)
         DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
         TGAUST=TGAUST+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
         TGAUIT=TGAUIT+SHAPET(INODLT,IGAUST)*TEINIT(1,INODLT)
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
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
        TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
        IF(KDYNAT.EQ.1) THEN
         IF(KINTET.EQ.1)       ! Euler method
     .    TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
        ENDIF
C
C**** UPDATE THE PHASE-CHANGE FUNCTION
C
        ILAHET=1
        ISINRT=2
C
        IF(NMEMO3.EQ.0) THEN
         CALL CAPCOFT(BASMM, BASCC, PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ELSE
         IX=NDIMETO-1
         IF(IMICR.EQ.1) THEN
          BASMM=EHISTT(1,IGAUST)             ! all computed in micr05t.f
          BASMI=EHISTT(6+IX,IGAUST)
          SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
          SOUR2T=EHISTT(4+IX,IGAUST)
         ELSE
          IF(MULPHT.EQ.1) THEN
           CALL CAPCOFT(BASMM, BASCC, PROPST,TGAUST,TGAUSX,
     .                  DTEMPT,SOUR1T,EHISTT(4+IX,IGAUST),ILAHET,ISINRT,
     .                  EHISTT(5+IX,IGAUST),COUTDT,BASMI,TGAUIT,PSEUDO)
           BASMM=EHISTT(1,IGAUST)                ! computed in frdy05t.f
           BASMI=EHISTT(6+IX,IGAUST)             ! computed in frdy05t.f
           SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
           SOUR2T=EHISTT(4+IX,IGAUST)
          ELSE
           CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                  DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                  DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
          ENDIF
         ENDIF
        ENDIF
C
C**** EVALUATES EFFECTIVE HEATS
C
        CONSTT=1.0D0
ctm     IF(KINTET.EQ.4) CONSTT=TALFAT
C
        IF(LINPC.EQ.1) SOUR2T=SOUR2T*DTEMPT*DTIMET
        DESILT=BASMI*SOUR2T/(DTIMET*CONSTT)
C
        IF(ITERME.GT.0) THEN                    ! bidirectional coupling
         IF(ITERMD.GT.0) THEN                   ! deformed shape
          IF(LARGET.EQ.1) THEN                  ! TLF
           DESILT=BASMI*SOUR2T/(DTIMET*CONSTT)
          ENDIF
          IF(LARGET.EQ.2) THEN                  ! ULF
           detjj=1.0d0  ! determinant of F (from mechanical computation)
           BASMM=BASMI/DETJJ
           DESILT=BASMM*SOUR2T/(DTIMET*CONSTT)
          ENDIF
          IF(LARGET.EQ.3) THEN                  ! Eulerian
c          BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
c          DESILT=BASMM*SOUR2T/(DTIMET*CONSTT)
          ENDIF
         ENDIF
        ENDIF
C
C**** INTEGRATES THE HEATS INTO THE INTERNAL "HEAT FORCES" DUE TO
C     TRANSIENT EFFECTS
C
        DO IEVABT=1,NNODLT
         ELEL2T(IEVABT)=ELEL2T(IEVABT)+WHAPEL(IEVABT,IGAUST)*DESILT*
     .                  DVOLUT(IGAUST)
        END DO
C
       END DO
C
      ENDIF
C
      RETURN
      END
