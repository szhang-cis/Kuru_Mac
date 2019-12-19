      SUBROUTINE FRIN05T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   WARTDL,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING HEATS 
C               ( FOR ELEMENT NO. 5 )
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
     .          EHISTT(NHISTT,*),        ELCODT(NDIMET,*),
     .          ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*), 
     .          LNODST(*),               PROPST(*),
     .          RMAT1T(NDIMET,*),        SHAPET(NNODLT,*),
     .          STRANT(NSTR1T,*),        STRSGT(NSTR1T,*),
     .          STRA0T(NSTR1T,*),        STRS0T(NSTR1T,*)
      DIMENSION BMATXT(NSTR1T,*),        BMSIGT(*),
     .          DESIGT(*),               DMATXT(NSTR1T,*),
     .          DSTRAT(*),               PRESGT(*),
     .          SGTOTT(*),               SIGMAT(*),
     .          TSTRAT(*),               XJACMT(NDIMET,*),
     .          ELELTT(*),               VELCMT(*)
      DIMENSION WARTDL(NDIMET,NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),         DVOLIT(*)
      DIMENSION COEFKT(9)
C
C**** INITIALIZE ELEMENT CONTRIBUTION OF ELELTT (K*T)
C
      DO IEVABT=1,NNODLT
       ELELTT(IEVABT)=0.0D0
      END DO
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
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
C**** UPDATES CONDUCTIVITY COEFFICIENT (COEFKT=Conductivitity)
C
       ILAHET=-1
       ISINRT=2
C
       IF(NMEMO3.EQ.0) THEN
        IDECIT=0
        IF(NMEMO10.EQ.1) IDECIT=1                ! density changes
        IF(IDECIT.EQ.1)
     .   CALL CAPCOFT(BASMM,BASCC,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        CALL CONCOFT(COEFKT,PROPST,TGAUST,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(IMICR.EQ.1) THEN
         BASMM= EHISTT(1,IGAUST)                 ! computed in micr05t.f
         BASMI= EHISTT(6+IX,IGAUST)              ! computed in micr05t.f
         DO I=1,NDIMETO
          COEFKT(I)=EHISTT(3+I-1,IGAUST)         ! computed in micr05t.f
         ENDDO
        ELSE
         IDECIT=0
         IF(NMEMO10.EQ.1) IDECIT=1               ! density changes
         IF(IDECIT.EQ.1)
     .    CALL CAPCOFT(EHISTT(1,IGAUST),BASCC,PROPST,TGAUST,TGAUSX,
     .                 DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                 DSOURT,COUTDT,EHISTT(6+IX,IGAUST),TGAUIT,PSEUDO)
         CALL CONCOFT(EHISTT(3,IGAUST),PROPST,TGAUST,PSEUDO)
         IF(IDECIT.EQ.1) THEN
          BASMM= EHISTT(1,IGAUST)
          BASMI= EHISTT(6+IX,IGAUST)
         ENDIF
         DO I=1,NDIMETO
          COEFKT(I)=EHISTT(3+I-1,IGAUST)
         ENDDO
        ENDIF
       ENDIF
C
C**** COMPUTE CONDUCTIVITY TENSOR     (only isotropic case available !!)
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
C**** CALCULATE THE TOTAL HEAT FLUX AT GAUSSIAN POINTS
C
       CALL LIHEATT(CARTDT(1,1,IGAUST),DEPSVT,SGTOTT,
     .              DSTRAT,
     .              ELDIST,GPCODT(1,IGAUST),LARGET,NDIMET,NDOFNT,NNODLT,
     .              NSTR1T,
     .              PROPST,SHAPET(1,IGAUST),STRANT(1,IGAUST),
     .              STRA0T(1,IGAUST),
     .              TSTRAT,XJACMT,NDOFCT,SIGMAT,VELCMT,DMATXT,KDYNAT)
C
C**** STORE TOTAL HEAT FLUX IN THE PROPER ARRAY
C
       IF(NMEMO4.EQ.1) THEN
        DO ISTR1T=1,NSTR1T
         STRSGT(ISTR1T,IGAUST)=-SGTOTT(ISTR1T)
         STRANT(ISTR1T,IGAUST)= TSTRAT(ISTR1T)
        ENDDO
       ENDIF
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEAT FORCES"
C
       CALL EQHEATT(BMATXT,BMSIGT,CARTDT(1,1,IGAUST),DVOLUT(IGAUST),
     .              ELDIST,
     .              GPCODT(1,IGAUST),LARGET,NDIMET,NDOFCT,NDOFNT,
     .              NEVABT,NNODLT,
     .              NSTR1T,NTYPET,SHAPET(1,IGAUST),SGTOTT,XJACMT,
     .              SIGMAT,ELELTT,WARTDL(1,1,IGAUST))
C 
      END DO      ! igaust=1,ngault
C
      RETURN
      END
