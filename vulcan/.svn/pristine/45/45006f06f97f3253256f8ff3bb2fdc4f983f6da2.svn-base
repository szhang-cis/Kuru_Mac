      SUBROUTINE STIF05T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,ELCODT,
     .                   BMATXT,DMATXT,SIGMAT,XJACMT,VELCMT,
     .                   WARTDL,TEINIT,FPCHLT,DVOLIT,AUXS1T)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE CONDUCTIVITY MATRIX
C     ( ELEMENT TYPE NO. 5 )
C
C***********************************************************************
C
C     Index of variables:
C
C     IX=NDIMETO-1
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
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION CARTDT(NDIMET,NNODLT,*), DVOLUT(*),
     .          EHISTT(NHISTT,*),        ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*),
     .          PROPST(*),               SHAPET(NNODLT,*),
     .          STRSGT(NSTR1T,*),        ESTIFT(*),
     .          HSTIFT(NEVABT,NNODET),   ELCODT(NDIMET,*)
      DIMENSION BMATXT(NSTR1T,*),        DMATXT(NSTR1T,*),
     .          SIGMAT(*),               XJACMT(NDIMET,*),
     .          VELCMT(*)
      DIMENSION WARTDL(NDIMET,NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),         DVOLIT(*),
     .          AUXS1T(NEVABT,*)
      DIMENSION COEFKT(9)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO IGAUST=1,NGAULT
C
       IF(NMEMO3.EQ.0) THEN
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
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
        ILAHET=-1
        ISINRT=1
C
        IDECIT=0
        IF(NMEMO10.EQ.1) IDECIT=1                ! density changes
        IF(IDECIT.EQ.1)
     .   CALL CAPCOFT(BASMM,BASCC,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
	CALL CONCOFT(COEFKT,PROPST,TGAUST,PSEUDO)
       ELSE
        IX=NDIMETO-1
        BASMM= EHISTT(1,IGAUST)
        BASMI= EHISTT(6+IX,IGAUST)
        DO I=1,NDIMETO
         COEFKT(I)=EHISTT(3+I-1,IGAUST)
        ENDDO
       ENDIF                          ! nmemo3.eq.0
C
C**** COMPUTES CONDUCTIVITY TENSOR    (only isotropic case available !!)
C
       IF(NMEMO10.EQ.1) THEN
        IF(BASMM.GT.0.0D0) THEN
         DO I=1,NDIMETO
          COEFKT(I)=COEFKT(I)*(BASMI/BASMM)
         ENDDO
        ENDIF
       ENDIF
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling
        IF(ITERMD.GT.0) THEN                    ! deformed shape
         IF(LARGET.EQ.1) THEN                   ! TLF
          DO I=1,NDIMETO
           detjj=1.0d0  ! determinant of F (from mechanical computation)
           BASMM=BASMI/DETJJ
           COEFKT(I)=COEFKT(I)*(BASMI/BASMM)
          ENDDO
         ENDIF
         IF(LARGET.EQ.2) THEN                   ! ULF
          DO I=1,NDIMETO
           COEFKT(I)=COEFKT(I)
          ENDDO
         ENDIF
         IF(LARGET.EQ.3) THEN                   ! Eulerian
          DO I=1,NDIMETO
           COEFKT(I)=COEFKT(I)
           IF(IFREKT.EQ.2)
     .      COEFKT(I)=COEFKT(I)*DVOLIT(IGAUST)/DVOLUT(IGAUST)
          ENDDO
         ENDIF
        ENDIF
       ENDIF
C
       I=1
       DO IDIMET=1,NDIMET                ! nstr1t=ndimet
        DO JDIMET=1,NDIMET
         DMATXT(IDIMET,JDIMET)=0.0D0
         IF(ISOTRT.EQ.0.OR.ISOTRT.EQ.1) THEN
          IF(IDIMET.EQ.JDIMET) THEN
           IF(ISOTRT.EQ.1) I=IDIMET
           DMATXT(IDIMET,JDIMET)=COEFKT(I)
          ENDIF
          IF(ISOTRT.EQ.2) THEN
           CALL RUNENDT('ERROR: FULLY ANIST. MAT. NOT IMPLEMENTED')
          ENDIF
         ENDIF
        ENDDO
       ENDDO
C
C**** CALCULATES THE MATERIAL CONDUCTIVITY MATRIX
C
       NSTRET=NDIMET
       NSTRST=NDIMET
       IF(LARGET.EQ.1.OR.LARGET.EQ.2)
     .  CALL RUNENDT('ERROR IN STIF05T: LARGE=1,2 NOT IMPLEMENTED')
C
       CALL KMATRIT(BMATXT,        CARTDT(1,1,IGAUST),  DMATXT,
     .              DVOLUT(IGAUST),ESTIFT,            GPCODT(1,IGAUST),
     .              LARGET,KSYMMT, NDIMET,NDOFNT,NEVABT,NKOVAT, NNODLT,
     .              NSTRET,NSTRST, NTYPET,SHAPET(1,IGAUST),
     .              WARTDL(1,1,IGAUST),AUXS1T)
C
      END DO         ! igaust=1,ngault
C
      RETURN
      END
