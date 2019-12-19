      SUBROUTINE MAUF05T(DVOLUT,PROPST,SHAPET,WSTIFT,EHISTT,VELCMT,
     .                   WSTI1T,WSTI2T,ELDIST,DISIMT,WHAPEL,TEINIT,
     .                   FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE LATENT HEAT CONTRIBUTION 
C     TO THE JACOBIAN MATRIX FOR UNIPHASE ELEMENTS
C     ( ELEMENT NO. 5 )
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
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
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
      DIMENSION DVOLUT(*),        PROPST(*),
     .          SHAPET(NNODLT,*)
      DIMENSION WSTIFT(*),        EHISTT(NHISTT,*)
      DIMENSION VELCMT(*),        WSTI1T(*),
     .          WSTI2T(*),        ELDIST(NDOFCT,*),
     .          DISIMT(*)
      DIMENSION WHAPEL(NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),  DVOLIT(*)
      DIMENSION EMATX(27,27)
C
      IF(IMICR.EQ.1) THEN
       IF(IMICO.EQ.0) RETURN
      ENDIF
C
C**** LOOP ON INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
C**** CALCULATION IN TIME T
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
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** UPDATES THE PHASE-CHANGE FUNCTION
C
       ILAHET=1
       ISINRT=1
C
       IF(NMEMO3.EQ.0) THEN
        CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(MULPHT.EQ.1.OR.IMICR.EQ.1) THEN
         BASMM=EHISTT(1,IGAUST)
	 BASMI=EHISTT(6+IX,IGAUST)
         SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
         SOUR2T=EHISTT(4+IX,IGAUST)
         IF(IMICR.EQ.1) THEN
          SOUR1T=0.0D0
          SOUR2T=EHISTT(8+IX,IGAUST)
          IF(ICACOT.EQ.1) CALL RECACOT(SOUR1T,SOUR2T,DTEMPT)
         ENDIF
        ELSE
         CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ENDIF            ! mulpht.eq.1.or.imicr.eq.1
       ENDIF             ! nmemo3.eq.0
C
       IF(IMICR.EQ.0) THEN
        VELTOT=1.D-10
        FACT1T=0.0D0
        VELA1T=DABS(DTEMPT)
        IF(VELA1T.GT.VELTOT) FACT1T=(SOUR1T/DTIMET)/DTEMPT
        IF(LINPC.EQ.1) FACT1T=SOUR1T
       ELSE
        FACT1T=SOUR1T
       ENDIF
C
C**** DETERMINES DENSITY
C
       BASMM=BASMI
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling
        IF(ITERMD.GT.0) THEN                    ! deformed shape
         IF(LARGET.EQ.1) THEN                   ! TLF
          BASMM=BASMI
         ENDIF
         IF(LARGET.EQ.2) THEN                   ! ULF
          detjj=1.0d0   ! determinant of F (from mechanical computation)
          BASMM=BASMI/DETJJ
         ENDIF
         IF(LARGET.EQ.3) THEN                   ! Eulerian
c
c computation of dvolit for element subdivision should be implemented !!
c
c         BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
         ENDIF
        ENDIF
       ENDIF
C
C**** CALCULATE THE MATRIX
C
       CALL WMHEATT(DVOLUT(IGAUST),NDOFNT,NEVABT,NNODLT,PROPST,
     .              SHAPET(1,IGAUST),WSTI1T,FACT1T,
     .              BASMM,COUTTT,
     .              WHAPEL(1,IGAUST),KSYMMT,EMATX)
C
      END DO
C
C**** LOOP ON INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
C**** CALCULATION IN TIME T+DT
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
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** UPDATES THE PHASE-CHANGE FUNCTION
C
       ILAHET=1
       ISINRT=1
C
       IF(NMEMO3.EQ.0) THEN
        CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(MULPHT.EQ.1.OR.IMICR.EQ.1) THEN
         BASMM=EHISTT(1,IGAUST)
         BASMI=EHISTT(6+IX,IGAUST)
         SOUR1T=EHISTT(4+IX,IGAUST)-EHISTT(5+IX,IGAUST)
         SOUR2T=EHISTT(4+IX,IGAUST)
         IF(IMICR.EQ.1) THEN
          SOUR1T=0.0D0
          SOUR2T=EHISTT(8+IX,IGAUST)
          IF(ICACOT.EQ.1) CALL RECACOT(SOUR1T,SOUR2T,DTEMPT)
         ENDIF
        ELSE
         CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ENDIF            ! mulpht.eq.1.or.imicr.eq.1
       ENDIF             ! nmemo3.eq.0
C
       IF(IMICR.EQ.0) THEN
        VELTOT=1.D-10
        FACT1T=0.0D0
        VELA1T=DABS(DTEMPT)
        IF(VELA1T.GT.VELTOT) FACT1T=(SOUR2T/DTIMET)/DTEMPT
        IF(LINPC.EQ.1) FACT1T=SOUR2T
       ELSE
        FACT1T=SOUR2T
       ENDIF
C
C**** DETERMINES DENSITY
C
       BASMM=BASMI
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling
        IF(ITERMD.GT.0) THEN                    ! deformed shape
         IF(LARGET.EQ.1) THEN                   ! TLF
          BASMM=BASMI
         ENDIF
         IF(LARGET.EQ.2) THEN                   ! ULF
          detjj=1.0d0   ! determinant of F (from mechanical computation)
          BASMM=BASMI/DETJJ
         ENDIF
         IF(LARGET.EQ.3) THEN                   ! Eulerian
c         BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
         ENDIF
        ENDIF
       ENDIF
C
C**** CALCULATE THE MATRIX
C
       CALL WMHEATT(DVOLUT(IGAUST),NDOFNT,NEVABT,NNODLT,PROPST,
     .              SHAPET(1,IGAUST),WSTI2T,FACT1T,
     .              BASMM,COUTTT,
     .              WHAPEL(1,IGAUST),KSYMMT,EMATX)
C
      END DO
C
      RETURN
      END
