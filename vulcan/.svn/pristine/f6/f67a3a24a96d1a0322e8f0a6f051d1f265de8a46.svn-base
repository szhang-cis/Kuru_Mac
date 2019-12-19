      SUBROUTINE MAS101T(DVOLUT,PROPST,SHAPET,WSTIFT,EHISTT,ELDIST,
     .                   VELCMT,BOUCHL,FPCHLT,EMATXT,GPCODT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE MASS MATRIX ( ELEMENT NO. 101 )
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
      DIMENSION DVOLUT(*), PROPST(*),        SHAPET(NNODLT,*)
      DIMENSION WSTIFT(*), EHISTT(NHISTT,*), ELDIST(NDOFCT,*)
      DIMENSION VELCMT(*)
      DIMENSION BOUCHL(NDOFCT,*), FPCHLT(NFPCH,*)
      DIMENSION EMATXT(NEVABT,*), GPCODT(NDIMET,*)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUST=1,NGAULT
C
      IF(NMEMO3.EQ.0) THEN
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
C
       DTEMPT=0.0D0
       TGAUST=0.0D0
       BGAUST=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        IF(KDYNAT.EQ.1)
     .  DTEMPT=DTEMPT+SHAPET(INODLT,IGAUST)*VELCMT(INODLT)
        TGAUST=TGAUST+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
        IF(NMEMO10.EQ.1)
     .   BGAUST=BGAUST+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT)
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
C**** TEMPERATURE AT GAUSS POINT IN TIME (t + alpha dt)
C
       IF(KDYNAT.EQ.1) THEN
        IF(KINTET.EQ.1)              ! Euler method
     .   TGAUST=TALFAT*TGAUST+(1.0D0-TALFAT)*(TGAUST-DTEMPT*DTIMET)
       ENDIF
C
       CALL ALTCONT(TGAUSX,TGAUST,GPCODT(1,IGAUST))
       CALL RADCOFT(COEFHT,PROPST,TGAUSX,TGAUSX,PSEUDO,PSEUDO,    2)
      ELSE
       BGAUST=0.0D0
       IF(NMEMO10.EQ.1) THEN
        DO INODLT=1,NNODLT
         BGAUST=BGAUST+SHAPET(INODLT,IGAUST)*BOUCHL(1,INODLT)
        END DO
       ENDIF
       COEFHT=EHISTT(1,IGAUST)
      ENDIF                          ! nmemo3.eq.0
C
      FACTRO=1.0D0
      IF(NMEMO10.EQ.1) THEN
       IF(NDIMET.EQ.1) THEN
        FACTRO=BGAUST
       ELSE
        FACTRO=BGAUST**(-(NDIMET-1)/NDIMET)
       ENDIF
      ENDIF
C
      IF(ITERME.GT.0) THEN                      ! bidirectional coupling
       IF(ITERMD.GT.0) THEN                     ! deformed shape
        IF(LARGET.EQ.1) THEN                    ! TLF
         detjj=1.0D0    ! determinant of F (from mechanical computation)
         FACTRO=DETJJ
         COEFHT=COEFHT
        ENDIF
        IF(LARGET.EQ.2) THEN                    ! ULF
         FACTRO=1.0D0
         COEFHT=COEFHT
        ENDIF
       ENDIF
      ENDIF
C
      COEFHT=COEFHT*FACTRO
C
C**** ELEMENTAL MASS MATRIX
C
      CALL WBOUNDT(DVOLUT(IGAUST),NDOFNT,NEVABT,NNODLT,PROPST,
     .             SHAPET(1,IGAUST),WSTIFT,COEFHT,ELDIST,
     .             NODFCT,KSYMMT,EMATXT,  101)
C
  100 CONTINUE
C
      RETURN
      END
