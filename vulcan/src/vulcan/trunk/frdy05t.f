      SUBROUTINE FRDY05T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,WHAPEL,TEINIT,FPCHLT,
     .                   DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES"
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
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), DISIMT(*),
     .          ELELTT(*)
      DIMENSION WHAPEL(NNODLT,*), TEINIT(NDOFCT,*),
     .          FPCHLT(NFPCH,*),  DVOLIT(*)
C
C**** INITIALIZE ELEMENT CONTRIBUTION
C
C     Note: this is only necessary if the transient terms are included
C           in the convergence criterion
C
      LHINC=0                          ! better as input (see therdyt.f)
      IF(LHINC.EQ.1) THEN
       DO IEVABT=1,NEVABT
        ELELTT(IEVABT)=0.0D0
       END DO
      ENDIF
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
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
C**** UPDATES THE HEAT CAPACITY COEFFICIENT
C
       ILAHET=0
       ISINRT=2
C
       IF(NMEMO3.EQ.0) THEN
        CALL CAPCOFT(BASMM,BASCC,PROPST,TGAUST,TGAUSX,
     .               DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .               DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
       ELSE
        IX=NDIMETO-1
        IF(IMICR.EQ.1) THEN
         BASMM=EHISTT(1,IGAUST)                  ! computed in micr05t.f
         BASMI=EHISTT(6+IX,IGAUST)               ! computed in micr05t.f
         BASCC=EHISTT(2,IGAUST)                  ! computed in micr05t.f
        ELSE
         IF(ICONVT.EQ.0) THEN
          IDECIT=0
          IF(NMEMO10.EQ.1) IDECIT=1              ! density changes
          IF(IDECIT.EQ.0) THEN
           CALL CAPCOFT(EHISTT(1,IGAUST),EHISTT(2,IGAUST),PROPST,TGAUST,
     .                                                           TGAUSX,
     .                  DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                  DSOURT,EHISTT(7+IX,IGAUST),EHISTT(6+IX,IGAUST),
     .                                                    TGAUIT,PSEUDO)
          ELSE
           CALL CAPCOFT(BASMM,EHISTT(2,IGAUST),PROPST,TGAUST,TGAUSX,
     .                  DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .                  DSOURT,EHISTT(7+IX,IGAUST),BASMI,TGAUIT,PSEUDO)
          ENDIF
         ENDIF
         BASMM=EHISTT(1,IGAUST)      ! computed in frin05t.f if IDECIT=1
C                                    ! or in fric05t.f if ICONVT=1
         BASMI=EHISTT(6+IX,IGAUST)   ! computed in frin05t.f if IDECIT=1
C                                    ! or in fric05t.f if ICONVT=1
         BASCC=EHISTT(2,IGAUST)      ! computed in fric05t.f if ICONVT=1
        ENDIF
       ENDIF
C
C**** EVALUATE EFFECTIVE HEATS
C
       DESIHT=BASMI*BASCC*DTEMPT
C
       IF(ITERME.GT.0) THEN                     ! bidirectional coupling
        IF(ITERMD.GT.0) THEN                    ! deformed shape
         IF(LARGET.EQ.1) THEN                   ! TLF
          DESIHT=BASMI*BASCC*DTEMPT
         ENDIF
         IF(LARGET.EQ.2) THEN                   ! ULF
          detjj=1.0d0   ! determinant of F (from mechanical computation)
          BASMM=BASMI/DETJJ
          DESIHT=BASMM*BASCC*DTEMPT
         ENDIF
         IF(LARGET.EQ.3) THEN                   ! Eulerian
          BASMM=BASMI*DVOLIT(IGAUST)/DVOLUT(IGAUST)
          DESIHT=BASMM*BASCC*DTEMPT
         ENDIF
        ENDIF
       ENDIF
C
C**** INTEGRATE THE HEATS INTO THE INTERNAL "HEATS FORCES" DUE TO
C     TRANSIENT EFFECTS
C
       DO IEVABT=1,NNODLT
        ELELMT(IEVABT)=ELELMT(IEVABT)+WHAPEL(IEVABT,IGAUST)*DESIHT*
     .                 DVOLUT(IGAUST)
        IF(LHINC.EQ.1) THEN
         ELELTT(IEVABT)=ELELTT(IEVABT)+SHAPET(IEVABT,IGAUST)*DESIHT*
     .                  DVOLUT(IGAUST) 
        ENDIF
       END DO 
C
      END DO ! IGAUST=1,NGAULT
C
      RETURN
      END
