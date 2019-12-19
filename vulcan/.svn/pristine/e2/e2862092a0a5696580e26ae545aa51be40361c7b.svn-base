      SUBROUTINE LDGR05T(DVOLUT,PROPST,SHAPET,WORK1T,ELDIST,VELCMT,
     .                   EHISTT,ILDGRT,WHAPEL,TEINIT,FPCHLT,DVOLIT,
     .                   ADVEMT,CARTDT,XJACMT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE EQUIVALENT "NODAL FORCES" DUE TO 
C     THE INTERNAL HEAT ( ELEMENT NO. 5 )
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
      COMMON/SUVOHEAT/SVHEAT
C
      DIMENSION DVOLUT(*),        PROPST(*),
     .          SHAPET(NNODLT,*), WORK1T(*),
     .          ELDIST(NDOFCT,*), VELCMT(*),
     .          EHISTT(NHISTT,*)
      DIMENSION WHAPEL(NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),  DVOLIT(*)
      DIMENSION ADVEMT(NDIMET,*), CARTDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*)
      DIMENSION TSTRAV(6),        TSTREV(6)
      DIMENSION PROACT(6),        VELOCT(3)
C
C**** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
C
      DO IGAUST=1,NGAULT
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
C**** UPDATES THE DENSITY
C
       ILAHET=-1
       ISINRT=2
C
       IF(NMEMO3.EQ.0) THEN
        CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(IMICR.EQ.1) THEN
         BASMM=EHISTT(1,IGAUST)      ! density computed in micr05t
         BASMI=EHISTT(6+IX,IGAUST)   ! density computed in micr05t
        ELSE
         IF(ICONVT.EQ.0) THEN
          IDECIT=0
          IF(NMEMO10.EQ.1) IDECIT=1               ! density changes
          IF(ILDGRT.EQ.0) THEN
           IF(KDYNAT.EQ.0) THEN
            IF(IDECIT.EQ.0) THEN
             CALL CAPCOFT(EHISTT(1,IGAUST), BASCC,PROPST,TGAUST,TGAUSX,
     .                    DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                    DSOURT,COUTDT,EHISTT(6+IX,IGAUST),TGAUIT,
     .                                                           PSEUDO)
            ENDIF
           ENDIF
          ELSE
           IF(KDYNAT.EQ.0) THEN
            IF(IDECIT.EQ.0) THEN
             CALL CAPCOFT(EHISTT(1,IGAUST), BASCC,PROPST,TGAUST,TGAUSX,
     .                    DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                    DSOURT,COUTDT,EHISTT(6+IX,IGAUST),TGAUIT,
     .                                                           PSEUDO)
            ENDIF
           ELSE
            IF(IDECIT.EQ.0) THEN ! EHISTT(1) & EHISTT(6+IX) will be
C                                ! recomputed in frdy05t.f
             CALL CAPCOFT(EHISTT(1,IGAUST), BASCC,PROPST,TGAUST,TGAUSX,
     .                    DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                    DSOURT,COUTDT,EHISTT(6+IX,IGAUST),TGAUIT,
     .                                                           PSEUDO)
            ENDIF
           ENDIF
          ENDIF                      ! ildgrt.eq.0
         ENDIF                       ! iconvt.eq.0
         BASMM=EHISTT(1,IGAUST)      ! computed in frin05t.f if IDECIT=1
C   ! or in fric05t.f if ICONVT=1 or in frdy05t.f if KDYNAT=1 & ILDGRT=0
         BASMI=EHISTT(6+IX,IGAUST)   ! computed in frin05t.f if IDECIT=1
C   ! or in fric05t.f if ICONVT=1 or in frdy05t.f if KDYNAT=1 & ILDGRT=0
        ENDIF
       ENDIF
C
C**** COMPUTES DENSITY
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
          BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
         ENDIF
        ENDIF
       ENDIF
C
       DENSET=BASMM
C
C**** REFERENCE SOURCE
C
       IF(ILDGRT.EQ.0) THEN
        SOURCT=SVHEAT
       ELSE
        SOURCT=0.0D0                 ! this should be improved !
C
        IF(ITERMEF.GT.0) THEN
         IF(ITERMPF.GT.0) THEN
C
C**** CALCULATE THE RATE-OF-DEFORMATION TENSOR
C
          CALL PROMA2(XJACMT,ADVEMT,CARTDT(1,1,IGAUST),NDIMET,NNODLT)
C
C**** CALCULATE VELOCITY AT GAUSS POINT
C
          DO IDIMET=1,NDIMET
           VELOCT(IDIMET)=0.0D0
           DO INODLT=1,NNODLT
            VELOCT(IDIMET)=VELOCT(IDIMET)+SHAPET(INODLT,IGAUST)*
     .                     ADVEMT(IDIMET,INODLT)
           ENDDO
          END DO
C
C**** CALCULATE VISCOSITY ACCORDING TO KVILAT
C
C     Note:
C
C     FPCHLT(1) should be the liquid-solid phase-change (see trasterf.f)
C
          VISCOT=VISC1F
          IF(IFILL.EQ.1) VISCOT=PSEUDO*VISC1F+(1.0D0-PSEUDO)*VISC2F
          IF(KCONST.GT.0) THEN
           IF(KVILAT.GT.0) THEN
            CRIPST=0.5D0                   ! to be improved
            PSCONT=1.0D0
            PROACT(2)=VISC1F
C
            FSLPG1=0.0D0
            DO INODLT=1,NNODLT
             FSLPG1=FSLPG1+SHAPET(INODLT,IGAUST)*FPCHLT(1,INODLT)
            END DO
C
            TGAUAT=TGAUST-DTEMPT*DTIMET
            IITER7=10000                 ! better as input
            IF(IITERT.GT.IITER7)
     .       TGAUST=TGAUAT               ! tuning to improve convergence
C
            CALL VISLAW(VISCOT,FPARATT,KVILAT,CRIPST,
     .                  VELOCT,XJACMT,TGAUST,PROACT,PSCONT,NDIMET,
     .                  FSLPG1,IERROR,DUMMY,DUMMY)
            IF(IERROR.EQ.1)
     .       CALL RUNENDT('ERROR: WRONG FREE PARAMETERS FOR THIS MODEL')
           ENDIF
          ENDIF
C
          DENSET=1.0D0              ! density not used in coupling term
C
          GO TO (1,2,3), NDIMET     ! NTYPET should be used
C
    1     NSTR1V=1
          TSTRAV(1)=XJACMT(1,1)
          TSTREV(1)=2.0D0*VISCOT*TSTRAV(1)
          GOTO 6
    2     NSTR1V=3
          TSTRAV(1)=XJACMT(1,1)
          TSTRAV(2)=XJACMT(2,2)
          TSTRAV(3)=XJACMT(1,2)+XJACMT(2,1)
          TSTREV(1)=2.0D0*VISCOT*TSTRAV(1)
          TSTREV(2)=2.0D0*VISCOT*TSTRAV(2)
          TSTREV(3)=      VISCOT*TSTRAV(3)
          GOTO 6
    3     NSTR1V=6
          TSTRAV(1)=XJACMT(1,1)
          TSTRAV(2)=XJACMT(2,2)
          TSTRAV(3)=XJACMT(1,2)+XJACMT(2,1)
          TSTRAV(4)=XJACMT(3,3)
          TSTRAV(5)=XJACMT(1,3)+XJACMT(3,1)
          TSTRAV(6)=XJACMT(2,3)+XJACMT(3,2)
          TSTREV(1)=2.0D0*VISCOT*TSTRAV(1)
          TSTREV(2)=2.0D0*VISCOT*TSTRAV(2)
          TSTREV(3)=      VISCOT*TSTRAV(3)
          TSTREV(4)=2.0D0*VISCOT*TSTRAV(4)
          TSTREV(5)=      VISCOT*TSTRAV(5)
          TSTREV(6)=      VISCOT*TSTRAV(6)
          GOTO 6
    6     CONTINUE
C
C**** COMPUTES THE VISCOUS ENERGY TERM
C
          VISCET=0.0D0
          DO ISTR1V=1,NSTR1V
           VISCET=VISCET+TSTREV(ISTR1V)*TSTRAV(ISTR1V)
          ENDDO
          SOURCT=VISCET
         ENDIF                      ! itermpf.gt.0
        ENDIF                       ! itermef.gt.0
        SOURCT=-SOURCT              ! treated as a residual heat
       ENDIF                        ! ildgrt.eq.1
C
C**** EVALUATES THE SOURCE TERM
C
       CALL SOURCTT(DENSET,SOURCT,QINTET)
C
C**** CALCULATE LOADS AND ASSOCIATE WITH ELEMENT NODAL POINTS
C
       IEVABT=0
       DO INODET=1,NNODLT
        SHAVOT=WHAPEL(INODET,IGAUST)*DVOLUT(IGAUST)
C
        IEVABT=IEVABT+1
        WORK1T(IEVABT)=WORK1T(IEVABT)+QINTET*SHAVOT
       ENDDO
C
      END DO      ! igaust=1,ngault
C
      RETURN
      END
