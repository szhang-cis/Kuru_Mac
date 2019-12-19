      SUBROUTINE STUC05T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,
     .                   VELCMT,WSTI1T,WSTI2T,WHAPEL,ADVEMT,TEINIT,
     .                   FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE JACOBIAN MATRIX DUE TO ADVECTIVE
C     PHASE-CHANGE EFFECTS
C     ( ELEMENT TYPE NO. 5 )
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
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      COMMON/MULPHASE/MULPHT
C
      COMMON/CONVE3/ICONV2
C
      DIMENSION CARTDT(NDIMET,NNODLT,*), DVOLUT(*),
     .          EHISTT(NHISTT,*),        ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*),
     .          PROPST(*),               SHAPET(NNODLT,*),
     .          STRSGT(NSTR1T,*),        ESTIFT(*),
     .          HSTIFT(NEVABT,NNODET)
      DIMENSION BMATXT(NSTR1T,*),        DMATXT(NSTR1T,*),
     .          SIGMAT(*),               XJACMT(NDIMET,*)
      DIMENSION VELCMT(*)
      DIMENSION WSTI1T(*),               WSTI2T(*)
      DIMENSION WHAPEL(NNODLT,*),        ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*),        FPCHLT(NFPCH,*)
      DIMENSION VELOCT(3)
C
      IF(ITERMEF.GT.0) THEN         ! skip isothermal pc for tf problems
       IF(KVILAT.EQ.102.OR.KVILAT.EQ.104.OR.KVILAT.EQ.107.OR.
     .    KVILAT.EQ.108) THEN
        IF(ISOLIF.EQ.0) RETURN      ! solid phase is not advected
       ENDIF
      ENDIF
C
      IF(IMICR.EQ.1) THEN
       IF(IMICO.EQ.0) RETURN
      ENDIF
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION (TIME t)
C
      DO IGAUST=1,NGAULT
C
C**** CALCULATION OF TEMPERATURE AT GAUSS POINT
C
       DTEMPT=0.0D0
       TGAUST=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        IF(KDYNAT.EQ.1)
     .   DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
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
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+SHAPET(INODLT,IGAUST)*
     .                  ADVEMT(IDIMET,INODLT)
        ENDDO
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
       IF(ICONV2.NE.0) THEN               ! to be revised !!! (Ag/99)
C
C**** UP DATE THE DENSITY & THE TEMPERATURE DERIVATIVE OF f_pc
C
        ILAHET=2
        ISINRT=1
C
        IF(NMEMO3.EQ.0) THEN
         CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,COEFLT,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ELSE
         IX=NDIMETO-1
         IF(MULPHT.EQ.1.OR.IMICR.EQ.1) THEN
          BASMM=EHISTT(1,IGAUST)        ! BASMM=Density
          BASMI=EHISTT(6+IX,IGAUST)     ! BASMI=Density
          COEFLT=EHISTT(8+IX,IGAUST)    ! COEFLT=L*(df_pc/dT)
c         COEFLT=EHISTT(9+IX,IGAUST)    ! COEFLT=L*(df_pc/dT)
         ELSE
          CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                 DTEMPT,SOUR1T,COEFLT,ILAHET,ISINRT,
     .                 DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
         ENDIF           ! mulpht.eq.1.or.imicr.eq.1
        ENDIF            ! nmemo3.eq.0
C
C**** COMPUTES DENSITY
C
        BASMM=BASMI
C
C**** COMPUTE DENSITY*LATENT HEAT*(df_pc/dT)*VELOCITY
C
        DO IDIMET=1,NDIMET            ! nstr1t=ndimet
         DMATXT(IDIMET,1)=BASMM*COEFLT*VELOCT(IDIMET)
        ENDDO
C
       ELSE            ! iconv2=0 (default)
C
C**** UP DATE THE DENSITY & THE TEMPERATURE DERIVATIVE OF f_pc
C
        ILAHET=1       ! approximated df_pc/dT to stabilize the solution
        ISINRT=1
C
	IF(NMEMO3.EQ.0) THEN
         CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,COEFLT,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
         COEFLT=DSOURT                  ! COEFLT=L*(\Delta df_pc)
        ELSE
         IX=NDIMETO-1
         IF(MULPHT.EQ.1.OR.IMICR.EQ.1) THEN
          BASMM=EHISTT(1,IGAUST)        ! BASMM=Density
          BASMI=EHISTT(6+IX,IGAUST)     ! BASMM=Density
          COEFLT=EHISTT(8+IX,IGAUST)    ! COEFLT=L*(df_pc/dT)
c         COEFLT=EHISTT(9+IX,IGAUST)    ! COEFLT=L*(df_pc/dT)
          IF(IMICR.EQ.1) THEN
           SOUR1T=0.0D0
           IF(ICACOT.EQ.1) CALL RECACOT(SOUR1T,COEFLT,DTEMPT)
          ENDIF
         ELSE
          CALL CAPCOFT(BASMM ,BASCC ,PROPST,TGAUST,TGAUSX,
     .                 DTEMPT,SOUR1T,COEFLT,ILAHET,ISINRT,
     .                 DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
          COEFLT=DSOURT                 ! COEFLT=L*(\Delta df_pc)
         ENDIF           ! mulpht.eq.1.or.imicr.eq.1
        ENDIF            ! nmemo3.eq.0
C
        IF(IMICR.EQ.0) THEN
         VELTOT=1.0D-10
         FACT1T=0.0D0
         VELA1T=DABS(DTEMPT)
         IF(VELA1T.GT.VELTOT) FACT1T=COEFLT/(DTEMPT*DTIMET)
         IF(LINPC.EQ.1) FACT1T=COEFLT                       ! L*df_pc/dT
        ELSE
         FACT1T=COEFLT
        ENDIF
C
C**** COMPUTES DENSITY
C
        BASMM=BASMI
C
C**** COMPUTE DENSITY*LATENT HEAT*(df_pc/dT)*VELOCITY
C
        DO IDIMET=1,NDIMET            ! nstr1t=ndimet
         DMATXT(IDIMET,1)=BASMM*FACT1T*VELOCT(IDIMET)
        ENDDO
C
       ENDIF    ! iconv2.ne.0
C
C**** CALCULATES THE CONVECTIVE MATRIX
C
       CALL KMATRIC(BMATXT,        CARTDT(1,1,IGAUST),  DMATXT,
     .              DVOLUT(IGAUST),ESTIFT,    GPCODT(1,IGAUST),
     .              LARGET,KSYMMT, NDIMET,NDOFNT,NEVABT,NKOVAT, NNODLT,
     .              NSTRET,NSTRST, NTYPET,WHAPEL(1,IGAUST))
C
      END DO ! IGAUS=1,NGAUL
C
      RETURN
      END
