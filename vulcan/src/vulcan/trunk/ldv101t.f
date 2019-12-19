      SUBROUTINE LDV101T(ELDIST,DVOLUT,PROPST,SHAPET,XJACMT,FORCET,
     .                   EHISTT,HTLODT,BOUCHL,FPCHLT,ELCO1T,GPCO1T,
     .                   ELCOIT,GPCODT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE EQUIVALENT NODAL FORCES DUE TO 
C     THE TRACTION FORCES ACTING ON THE SURFACE 2D/3D AS VOLUME
C     FORCES ON BOUNDARY ELEMENTS 1D/2D ( ELEMENT NO. 101 )
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
      COMMON/SUVOHEAT/SVHEAT
C
      DIMENSION DVOLUT(*),             ELDIST(NDOFCT,*),
     .          PROPST(*),             SHAPET(NNODLT,*),
     .          XJACMT(NDIMET,NDIMLT), FORCET(NEVABT),
     .          EHISTT(NHISTT,*),      HTLODT(NHLODT,NSUBFT,NFUNCT)
      DIMENSION BOUCHL(NDOFCT,*),      FPCHLT(NFPCH,*),
     .          ELCO1T(NDIMET,*),      GPCO1T(NDIMET,*),
     .          ELCOIT(NDIMET,*),      GPCODT(NDIMET,*)
C
C**** ENTER LOOPS FOR LENGHT/AREA NUMERICAL INTEGRATION
C
      DO IGAUST=1,NGAULT
C
C**** TEMPERATURE AT GAUSS POINT IN TIME t+dt
C
       TGAUST=0.0D0
       PSEUDO=0.0D0
       DO INODLT=1,NNODLT
        TGAUST=TGAUST+SHAPET(INODLT,IGAUST)*ELDIST(1,INODLT)
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
C**** EVALUATE BOUNDARY HEAT FLUX
C
C     Note:
C     AS AN EXTERNAL HEAT >>> SIGN +, CONSTANT BUT TIME DEPENDENT LOAD
C     AS A RESIDUAL HEAT  >>> SIGN -, NOT CONSTANT
C     (see elm101t.f)
C
       IF(ILDV1T.EQ.1) FACTAT=SVHEAT              ! as an external heat
       IF(ILDV2T.EQ.1) THEN                       ! as a residual heat
        CALL SOURCB(QHEATT,HTLODT,ELCO1T,ELCOIT,  ! computes heat
     .              GPCO1T(1,IGAUST),GPCODT(1,IGAUST),SHAPET(1,IGAUST))
        FACTAT=-QHEATT                            ! SVHEAT should be 1.0
       ENDIF
C
C**** CALCULATE LOADS AND ASSOCIATE WITH ELEMENT NODAL POINTS
C
       IEVABT=0
       DO INODLT=1,NNODLT
         IEVABT=IEVABT+1
         FORCET(IEVABT)=FORCET(IEVABT)+FACTAT*SHAPET(INODLT,IGAUST)*
     .                  DVOLUT(IGAUST)
       END DO
C
      END DO                         ! igaus=1,ngaul
C
      RETURN
      END
