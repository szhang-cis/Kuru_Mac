      SUBROUTINE FRIC05T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   WHAPEL,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE INTERNAL RESISTING HEATS DUE TO
C     CONVECTION EFFECTS
C     ( FOR ELEMENT NO. 5 )
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
      DIMENSION WHAPEL(NNODLT,*),        ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*),        FPCHLT(NFPCH,*)
      DIMENSION VELOCT(3)
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
     .  DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
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
        IF(KINTET.EQ.1)       ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** UPDATE THE DENSITY & CAPACITY COEFFICIENT
C
       ILAHET=0
       ISINRT=2
C
       IF(NMEMO3.EQ.0) THEN
        CALL CAPCOFT(COEFMT,COEFCT,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(IMICR.EQ.1) THEN
         COEFMT=EHISTT(1,IGAUST)  ! COEFMT=Density computed in micr05t.f
         BASMI= EHISTT(6+IX,IGAUST) ! Init. density comp. in micr05t.f
         COEFCT=EHISTT(2,IGAUST)  ! COEFCT=Capacity comp. in micr05t.f
        ELSE
         IDECIT=0
         IF(NMEMO10.EQ.1) IDECIT=1               ! density changes
         IF(IDECIT.EQ.0) THEN
          CALL CAPCOFT(EHISTT(1,IGAUST),EHISTT(2,IGAUST),PROPST,TGAUST,
     .                                                          TGAUSX,
     .                 DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                 DSOURT,EHISTT(7+IX,IGAUST),EHISTT(6+IX,IGAUST),
     .                                                    TGAUIT,PSEUDO)
         ELSE
          CALL CAPCOFT(COEFMT,EHISTT(2,IGAUST),PROPST,TGAUST,TGAUSX,
     .                 DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                 DSOURT,EHISTT(7+IX,IGAUST),BASMI,TGAUIT,PSEUDO)
         ENDIF
         COEFMT=EHISTT(1,IGAUST)     ! computed in frin05t.f if IDECIT=1
         BASMI= EHISTT(6+IX,IGAUST)  ! computed in frin05t.f if IDECIT=1
         COEFCT=EHISTT(2,IGAUST)     ! COEFCT=Capacity
        ENDIF
       ENDIF
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
C**** CALCULATE THE CONVECTIVE TERM
C
       CALL LIHEATC(DSTRAT,TSTRAT,DMATXT,CONVET,NDIMET,NSTR1T,
     .              NNODLT,NDOFCT,CARTDT(1,1,IGAUST),
     .              VELCMT,ELDIST,KDYNAT)
C
C**** CALCULATES THE "CONVECTIVE FORCES"
C
       CALL EQHEATC(BMSIGT,WHAPEL(1,IGAUST),DVOLUT(IGAUST),CONVET,
     .              NNODLT,NEVABT)
C
      END DO ! IGAUST=1,NGAULT
C
      RETURN
      END
