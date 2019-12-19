      SUBROUTINE ELM004(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,
     .                  EMASS,                                    !ELDAT
     .                  STRA0,STRS0,TEMPC,                        !ELPRE
     .                  ELDIS,EHIST,STRAN,STRSG,                  !ELVAR
     .                  CSTIF,ESTIF,WSTIF,PSTIF,QSTIF,HSTIF,      !ELMAT
     .                                                      ITASK)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 4 :
C
C     2D NODAL CONTACT ELEMENT
C
C     FUNCTIONS:
C
C       ITASK =  3   Evaluate Stiffness Matrix
C       ITASK =  4   Evaluate Internal Resisting Forces
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING4/RITEN,RITEF,TRUPL,TRUPM,TOLGA,TOLGAM
C
      DIMENSION LNODS(*), PROEL(*), PROPS(*), INFRI(*), COFRI(*),
     .          NOPRF(*), PRESF(*)
      DIMENSION ELCOD(*), CARTD(*), DVOLU(*), GPCOD(*), SHAPE(*), 
     .          EPMTX(*), RMAT1(*), EMASS(*)
      DIMENSION STRA0(*), STRS0(*), TEMPC(*)
      DIMENSION ELDIS(*), EHIST(*), STRAN(*), STRSG(*)
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*), PSTIF(*), QSTIF(*),
     .          HSTIF(*)
      DIMENSION WORK1(*)
C
C**** SELECT ELEMENT PROPERTIES
C
      NDIML=1
      NNODL=INT(PROEL( 2))
      NRULE=INT(PROEL( 3))
      NTYPE=INT(PROEL( 6))
C
      RITEN=    PROEL( 8)
      RITEF=    PROEL( 9)
      IF(ICONC.EQ.0) THEN
       TRUPL=   PROEL(10)
       TRUPM=   PROEL(11)
      ENDIF
      TOLGA=    PROEL(12)
      TOLGAM=PROEL(13)
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10,  1, 30, 40,  1,  1,  1,  1,  1,  1,
     .         1,  1,  1,  1,  1,  1,  1,  1,  1), ITASK
    1 RETURN ! Nothing
C
C**** SET UP A CONSTANT MATRIX
C
   10 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      CALL SETI04(WORK1(ISETM( 9)),                        ! VNORL
     .            ELCOD,ELDIS,
     .            WORK1(ISETM(10)),                        ! ELCOD
     .            WORK1(ISETM(17)),                        ! ELDI1
     .            ITASK)
      RETURN
C
C**** EVALUATE STIFFNESS
C
   30 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      CALL SETI04(WORK1(ISTIF( 8)),                        ! VNORL
     .            ELCOD,ELDIS,
     .            WORK1(ISTIF( 9)),                        ! ELCOD
     .            WORK1(ISTIF(27)),                        ! ELDI1
     .            ITASK)
C
      CALL ASIL00(WORK1(ISTIF(25)),WORK1(ISTIF(26)),       ! ESTIF,WSTIF
     .            NKOVA,KDYNA,
     .            WORK1(ISTIF(34)),                        ! EMATX
     .            NEVAB,    1)
C
      IF(NPOIC.EQ.0) THEN
       IF(NRULE.EQ.10) THEN
        CALL STIF04A(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ELSE
        CALL STIF04(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ENDIF
      ELSE          ! npoic > 0
       IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
        CALL STIF04B(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ENDIF
       IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
        CALL STIF04C(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ENDIF
       IF(IAUGM.EQ.4) THEN
        CALL STIF04Y(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ENDIF
       IF(IAUGM.EQ.5) THEN
        CALL STIF04Z(PROPS,
     .              WORK1(ISTIF(25)),                      ! ESTIF
     .              WORK1(ISTIF( 5)),                      ! TENOD
     .              WORK1(ISTIF(27)),                      ! ELDIS
     .              WORK1(ISTIF( 9)),                      ! ELCOD
     .              WORK1(ISTIF( 7)),WORK1(ISTIF( 8)))     ! DISPL,VNORL
       ENDIF
      ENDIF         ! npoic.eq.0
C
C**** ASSEMBLY PROCESS (NMEMO6M=1)
C
      CALL ASAL00(ESTIF,WSTIF,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            WORK1(ISTIF(26)),                        ! WSTIF
     .            NKOVA,KDYNA,NMEMO6M)
C
      CALL ASEL04(CSTIF,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            WORK1(ISTIF(26)),                        ! WSTIF
     .            PROPS,
     .            WORK1(ISTIF(28)),                        ! CSTI1
     .            LNODS,
     .            WORK1(ISTIF(32)),WORK1(ISTIF(33)),       ! AUXMA,TRAMA
     .            INFRI,COFRI)
C
      RETURN
C
C**** EVALUATE INTERNAL RESISTING FORCES
C
   40 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      CALL SETI04(WORK1(IFORC(18)),                        ! VNORL
     .            ELCOD,ELDIS,
     .            WORK1(IFORC(19)),                        ! ELCOD
     .            WORK1(IFORC(37)),                        ! ELDI1
     .            ITASK)
C
      IF(NPOIC.EQ.0) THEN
       IF(NRULE.EQ.10) THEN
        CALL FRIN04A(WORK1(IFORC(37)),                     ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)))     ! TGAPL,VNORL
       ELSE
        CALL FRIN04(WORK1(IFORC(37)),                      ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)))     ! TGAPL,VNORL
       ENDIF
      ELSE          ! npoic > 0
       IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
        CALL FRIN04B(WORK1(IFORC(37)),                     ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)))     ! TGAPL,VNORL
       ENDIF
       IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
        CALL FRIN04C(WORK1(IFORC(37)),                     ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)),     ! TGAPL,VNORL
     .              WORK1(IFORC(38)))                      ! DISIL
       ENDIF
       IF(IAUGM.EQ.4) THEN
        CALL FRIN04Y(WORK1(IFORC(37)),                     ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)))     ! TGAPL,VNORL
       ENDIF
       IF(IAUGM.EQ.5) THEN
        CALL FRIN04Z(WORK1(IFORC(37)),                     ! ELDIS
     .              PROPS,
     .              WORK1(IFORC( 1)),                      ! BMSIG
     .              WORK1(IFORC(11)),                      ! TENOD
     .              WORK1(IFORC(19)),                      ! ELCOD
     .              WORK1(IFORC(14)),WORK1(IFORC(16)),     ! PWOEL,PREAL
     .              WORK1(IFORC(17)),WORK1(IFORC(18)))     ! TGAPL,VNORL
       ENDIF
      ENDIF         ! npoic.eq.0
      RETURN
C
      END
