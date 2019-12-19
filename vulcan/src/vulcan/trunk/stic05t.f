      SUBROUTINE STIC05T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,WHAPEL,ADVEMT,VELCMT,
     .                   TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE JACOBIAN MATRIX DUE TO CONVECTION
C     EFFECTS
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
     .          HSTIFT(NEVABT,NNODET)
      DIMENSION BMATXT(NSTR1T,*),        DMATXT(NSTR1T,*),
     .          SIGMAT(*),               XJACMT(NDIMET,*)
      DIMENSION WHAPEL(NNODLT,*),        ADVEMT(NDIMET,*)
      DIMENSION VELCMT(*),               TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*)
      DIMENSION VELOCT(3)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO IGAUST=1,NGAULT
C
       DO IDIMET=1,NDIMET
        VELOCT(IDIMET)=0.0D0
        DO INODLT=1,NNODLT
         VELOCT(IDIMET)=VELOCT(IDIMET)+SHAPET(INODLT,IGAUST)*
     .                  ADVEMT(IDIMET,INODLT)
        ENDDO
       END DO
C
C**** COMPUTE DENSITY*CAPACITY*VELOCITY
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
        ILAHET=0
        ISINRT=1
C
	CALL CAPCOFT(COEFMT,COEFCT,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,COEFLT,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        COEFMT=EHISTT(1,IGAUST)          ! COEFMT=Density
        BASMI= EHISTT(6+IX,IGAUST)       ! BASMI=Initial density
        COEFCT=EHISTT(2,IGAUST)          ! COEFCT=Capacity
       ENDIF                     ! nmemo3.eq.0
C
C**** COMPUTES DENSITY
C
       COEFMT=BASMI
C
C**** COMPUTE DENSITY*CAPACITY*VELOCITY
C
       DO IDIMET=1,NDIMET               ! nstr1t=ndimet
        DMATXT(IDIMET,1)=COEFMT*COEFCT*VELOCT(IDIMET)
       ENDDO
C
C**** CALCULATES THE CONVECTIVE MATRIX
C
       IF(LARGET.EQ.1.OR.LARGET.EQ.2)       
     .  CALL RUNENDT('ERROR IN STIF05T: LARGE=1,2 NOT IMPLEMENTED')
C
       CALL KMATRIC(BMATXT,        CARTDT(1,1,IGAUST),  DMATXT,
     .              DVOLUT(IGAUST),ESTIFT,            GPCODT(1,IGAUST),
     .              LARGET,KSYMMT, NDIMET,NDOFNT,NEVABT,NKOVAT, NNODLT,
     .              NSTRET,NSTRST, NTYPET,WHAPEL(1,IGAUST))
C
      END DO ! IGAUS=1,NGAUL
C
      RETURN
      END
