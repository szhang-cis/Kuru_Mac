      SUBROUTINE MICR05T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   VEL1MT,ADVEMT,TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH A MICROSTRUCTURAL MODEL
C     ( FOR ELEMENT NO. 5 )
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
C     EHISTT(NHISTT-NHISTM+1:NHISTT)=array of microstructural variables
C
C-----------------------------------------------------------------------
C
C     VELCMT = current temperature rate
C     VEL1MT = predicted temperature rate (at time t)
C
C-----------------------------------------------------------------------
C
C     Notes:
C
C     In thermal problems, if NRULET=11/13 is given in the input data
C     file, only the transient terms are integrated at nodes.
C     In microstructural problems, if NRULET=11/13 is given in the
C     input data file, f_pc must be computed at nodes and, taking into
C     account that f_pc is a internal variable, the integration points
C     must coincide with the nodes, i.e., all the energy terms are
C     integrated at nodes (not only the transient ones) by means of
C     NEGATT=2 (see elm005t.f) for all these terms
C
C     The use of NRULET=5 (integration points at nodes) in the smoothing
C     operation (see elm005t.f) is not absolutely correct when NRULET=1
C     (Gauss integration points) is given in the input data file
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
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
      DIMENSION VEL1MT(*),               ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*),        FPCHLT(NFPCH,*)
C
      IX=NDIMETO-1
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
C**** COMPUTES TEMP. RATE, INITIAL TEMP. & TEMPERATURE AT GAUSSIAN POINT
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
         IPSEU=2*NNUPT+NNUPO+1
         PSEUDO=PSEUDO+SHAPET(INODLT,IGAUST)*FPCHLT(IPSEU,INODLT)
        ENDIF
       END DO
       DTEMPX=DTEMPT
C
       IF(IMICO.EQ.0.OR.LINPC.EQ.1) THEN
        DTEM1T=0.0D0
        DO INODLT=1,NNODLT
         DTEM1T=DTEM1T+SHAPET(INODLT,IGAUST)*VEL1MT(INODLT)
        END DO
        DTEMPX=DTEM1T
       ENDIF
C
C**** TEMPERATURE AT GAUSS POINT IN TIME t
C
       TGAUAT=TGAUST-DTEMPT*DTIMET
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + dt)
C
       TGAUSX=TGAUST
C
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1) THEN
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** COMPUTES CONVECTIVE TERM ASSOCIATED WITH MICROSTRUCTURAL VARIABLES
C     (approximation at time t; otherwise a field equation should be
C     solved for the microstructural evolution equations)
C
       IF(ICONVT.EQ.1) THEN
        IF(NNUPO.GT.0) THEN
         DO INUPO=1,NNUPO
          APLUOT(INUPO)=0.0D0
          DO IDIMET=1,NDIMET
           APLUOX=0.0D0
           VELOCX=0.0D0
           DO INODLT=1,NNODLT
            APLUOX=APLUOX+
     .         CARTDT(IDIMET,INODLT,IGAUST)*FPCHLT(INUPO+2*NNUPT,INODLT)
            VELOCX=VELOCX+SHAPET(INODLT,IGAUST)*ADVEMT(IDIMET,INODLT)
           ENDDO
           APLUOT(INUPO)=APLUOT(INUPO)+APLUOX*VELOCX
          ENDDO
         ENDDO
         TEMPCO=0.0D0
         DO IDIMET=1,NDIMET
          TEMPGR=0.0D0
          VELOCX=0.0D0
          DO INODLT=1,NNODLT
           TEMPGR=TEMPGR+CARTDT(IDIMET,INODLT,IGAUST)*ELDIST(1,INODLT)
           VELOCX=VELOCX+SHAPET(INODLT,IGAUST)*ADVEMT(IDIMET,INODLT)
          ENDDO
          TEMPCO=TEMPCO+TEMPGR*VELOCX
         ENDDO
c        DTEMPX=DTEMPX+TEMPCO                 ! material derivative of T
        ENDIF            ! nnupo.gt.0
       ENDIF             ! iconvt.eq.1
C
C**** UPDATES DENSITY & CAPACITY
C
       ILAHET=0
       ISINRT=2
       CALL CAPCOFT(EHISTT(1,IGAUST),EHISTT(2,IGAUST),PROPST,TGAUST,
     .                                                       TGAUSX,
     .              DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .              DSOURT,COUTDT, EHISTT(6+IX,IGAUST),TGAUIT,PSEUDO)
C
C**** UPDATES CONDUCTIVITY
C
       CALL CONCOFT(EHISTT(3,IGAUST),PROPST,TGAUST,PSEUDO)
C
C**** UPDATES THE PHASE-CHANGE FUNCTION, DENSITY, CAPACITY & 
C     CONDUCTIVITY ACCORDING WITH A MICROSTRUCTURAL MODEL
C
       CALL MICROST(PROPST,TGAUAT,TGAUSX,TGAUST,TGAUIT,
     .              DTEMPX,DTEMPT,
     .              EHISTT(1,IGAUST),EHISTT(2,IGAUST),
     .              EHISTT(3,IGAUST),EHISTT(4+IX,IGAUST),
     .              EHISTT(5+IX,IGAUST),EHISTT(8+IX,IGAUST),
     .              EHISTT(9+IX,IGAUST),
     .              EHISTT(NHISTT-NHISTM+1,IGAUST))
C
      END DO ! IGAUST=1,NGAULT
C
      RETURN
      END
