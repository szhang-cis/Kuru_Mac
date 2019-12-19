      SUBROUTINE MARX05T(DVOLUT,PROPST,SHAPET,ESTIFT,EHISTT,WHAPEL,
     .                   ELDIST,VELCMT,TEINIT,FPCHLT,DVOLIT,EMATXT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE HEAT CAPACITY MATRIX ( ELEMENT NO. 5 )
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
      DIMENSION DVOLUT(*),        PROPST(*),
     .          SHAPET(NNODLT,*)
      DIMENSION ESTIFT(*),        EHISTT(NHISTT,*)
      DIMENSION WHAPEL(NNODLT,*)
      DIMENSION ELDIST(NDOFCT,*), VELCMT(*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
      DIMENSION EMATXT(NEVABT,*)
C
      IF(ITERMEF.LE.0) RETURN
      IF(KTEM1F.NE.3) RETURN

      isanto=0
      if(isanto.eq.0) then

C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUST=1,NGAULT
C
C**** COMPUTE TEMPERATURE AT GAUSSIAN POINT
C
      IF(NMEMO3.EQ.0) THEN
C
       DTEMPT=0.0D0
       TGAUST=0.0D0
       TGAUIT=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
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
       ILAHET=0
       ISINRT=1
C
       CALL CAPCOFT( BASMM, BASCC,PROPST,TGAUST,TGAUSX,
     .              DTEMPT,SOUR1T,SOUR2T,ILAHET,ISINRT,
     .              DSOURT,COUTDT, BASMI,TGAUIT,PSEUDO)
      ELSE
       IX=NDIMETO-1
       BASMM=EHISTT(1,IGAUST)
       BASMI=EHISTT(6+IX,IGAUST)
       BASCC=EHISTT(2,IGAUST)
      ENDIF                          ! nmemo3.eq.0
C
C**** DETERMINES CAPACITY
C
      BASMM=BASMI*TFACTF
C
C**** ELEMENTAL HEAT CAPACITY MATRIX
C
      CALL WMHEATT(DVOLUT(IGAUST),NDOFNT,NEVABT,NNODLT,PROPST,
     .             SHAPET(1,IGAUST),ESTIFT,BASCC,
     .             BASMM,COUTTT,
     .             WHAPEL(1,IGAUST),KSYMMT,EMATXT)
C
  100 CONTINUE
C
      endif






C
      RETURN
      END
