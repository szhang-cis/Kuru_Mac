      SUBROUTINE FRRA05T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,WHAPEL,TEINIT,FPCHLT,DVOLIT,ADVEMT,
     .                   CARTDT,XJACMT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE HEAT TERM FOR THE IMPROVED
C     THERMO-FLUID ALGORITHM
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
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*)
      DIMENSION WHAPEL(NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),  DVOLIT(*)
C
      DIMENSION CARTDT(NDIMET,NNODLT,*), ADVEMT(NDIMET,*)
C
      DIMENSION XJACMT(NDIMET,*)
      DIMENSION TSTRAV(6),        TSTREV(6)
C
      IF(IGALFA.EQ.0) RETURN
C
      IF(IFILL.EQ.1) THEN
       IF(IMICR.EQ.0) THEN
        IPSEU=2*NNUPT+1
       ELSE
        IPSEU=2*NNUPT+NNUPO+1
       ENDIF
       ICOUP=IPSEU+1
      ELSE
       ICOUP=2*NNUPT+1
      ENDIF
C
C**** LOOP OVER INTEGRATION POINTS
C
      DO IGAUST=1,NGAULT
C
       DO KDIMET=1,NDIMET                           ! momentum equations
C
C**** PRESSURE TERM
C
        DESIHP=0.0D0
        DO INODLT=1,NNODLT
         DESIHP=DESIHP+CARTDT(KDIMET,INODLT,IGAUST)*FPCHLT(ICOUP,INODLT)
        ENDDO
C
C**** BOUSSINESQ TERM
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
     .    TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
        ENDIF
C
C**** UPDATES THE DENSITY
C
        ILAHET=-1
        ISINRT=2
C
        IF(NMEMO3.EQ.0) THEN
         CALL CAPCOFT(BASMM,BASCC,PROPST,TGAUST,TGAUSX,
     .                DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
        ELSE
         IX=NDIMETO-1
         IF(IMICR.EQ.1) THEN
          BASMM=EHISTT(1,IGAUST)                 ! computed in micr05t.f
          BASMI=EHISTT(6+IX,IGAUST)              ! computed in micr05t.f
         ELSE
          IF(ICONVT.EQ.0) THEN
           IDECIT=0
           IF(NMEMO10.EQ.1) IDECIT=1             ! density changes
           IF(IDECIT.EQ.0) THEN
            CALL CAPCOFT(EHISTT(1,IGAUST),EHISTT(2,IGAUST),PROPST,
     .                                                    TGAUST,TGAUSX,
     .                   DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                   DSOURT,EHISTT(7+IX,IGAUST),EHISTT(6+IX,IGAUST),
     .                                                    TGAUIT,PSEUDO)
           ELSE
            CALL CAPCOFT(BASMM,EHISTT(2,IGAUST),PROPST,TGAUST,TGAUSX,
     .                   DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                   DSOURT,EHISTT(7+IX,IGAUST),BASMI,TGAUIT,PSEUDO)
           ENDIF
          ENDIF
          BASMM=EHISTT(1,IGAUST)     ! computed in frin05t.f if IDECIT=1
C                                    ! or in fric05t.f if ICONVT=1
          BASMI=EHISTT(6+IX,IGAUST)  ! computed in frin05t.f if IDECIT=1
C                                    ! or in fric05t.f if ICONVT=1
         ENDIF
        ENDIF
C
C**** EVALUATE EFFECTIVE HEATS
C
        DILAF=DILA1F
        IF(IFILL.EQ.1) DILAF=PSEUDO*DILA1F+(1.0D0-PSEUDO)*DILA2F
C
        DESIHT=BASMI*GRAVYFF-BASMI*GRATEFF*DILAF*(TGAUST-TREFEFF)
        DESIHT=DESIHT*GVECTFF(KDIMET)
C
C**** CONVECTIVE TERM
C
        DESIHC=0.0D0
        DO IDIMET=1,NDIMET
         DESIH1=0.0D0
         DO INODLT=1,NNODLT
          DESIH1=DESIH1+
     .                CARTDT(IDIMET,INODLT,IGAUST)*ADVEMT(KDIMET,INODLT)
         ENDDO
         DESIH2=0.0D0
         DO INODLT=1,NNODLT
          DESIH2=DESIH2+SHAPET(INODLT,IGAUST)*ADVEMT(IDIMET,INODLT)
         ENDDO
         DESIHC=DESIHC+DESIH1*DESIH2
        ENDDO
        DESIHC=BASMI*DESIHC
C
C**** INERTIAL TERM
C
        DESIHI=0.0D0
        IF(KDYNAT.EQ.1) THEN
         DO INODLT=1,NNODLT
          DESIHI=DESIHI+SHAPET(INODLT,IGAUST)*ADVEMT(KDIMET,INODLT)
         ENDDO
         DESIHI=BASMI*DESIHI
        ENDIF
C
C**** DIFFUSIVE TERM
C
        VISCOT=VISC1F
        IF(IFILL.EQ.1) VISCOT=PSEUDO*VISC1F+(1.0-PSEUDO)*VISC2F
C
C**** CALCULATE THE RATE-OF-DEFORMATION TENSOR
C
        CALL PROMA2(XJACMT,ADVEMT,CARTDT(1,1,IGAUST),NDIMET,NNODLT)
C
        GO TO (1,2,3), NDIMET     ! NTYPET should be used
C
    1   NSTR1V=1
        TSTRAV(1)=XJACMT(1,1)
        TSTREV(1)=2.0D0*VISCOT*TSTRAV(1)
        GOTO 6
    2   NSTR1V=3
        TSTRAV(1)=XJACMT(1,1)
        TSTRAV(2)=XJACMT(2,2)
        TSTRAV(3)=XJACMT(1,2)+XJACMT(2,1)
        TSTREV(1)=2.0D0*VISCOT*TSTRAV(1)
        TSTREV(2)=2.0D0*VISCOT*TSTRAV(2)
        TSTREV(3)=      VISCOT*TSTRAV(3)
        GOTO 6
    3   NSTR1V=6
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
    6   CONTINUE
C
C**** EVALUATES RESIDUAL
C
        DESIHF=DESIHI+DESIHP+DESIHC-DESIHT
C
        DO IEVABT=1,NNODLT
         TESTF=SHAPET(IEVABT,IGAUST)*GALFA*GVECTFF(KDIMET)
         ELELMT(IEVABT)=ELELMT(IEVABT)+TESTF*DESIHF*DVOLUT(IGAUST)
C
         DO IDIMET=1,NDIMET
          TESTFK=CARTDT(IDIMET,IEVABT,IGAUST)*GALFA*GVECTFF(KDIMET)
          ELELMT(IEVABT)=ELELMT(IEVABT)+
     .                              TESTFK*TSTREV(IDIMET)*DVOLUT(IGAUST)
         ENDDO
        END DO
C
       ENDDO      ! kdimet=1,ndimet
C
      END DO      ! igaust=1,ngault
C
      RETURN
      END
