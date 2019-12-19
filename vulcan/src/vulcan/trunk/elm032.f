      SUBROUTINE ELM032(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,
     .                  EMASS,VNORI,                              !ELDAT
     .                  STRA0,STRS0,TEMPC,                        !ELPRE
     .                  ELDIS,EHIST,STRAN,STRSG,                  !ELVAR
     .                  CSTIF,ESTIF,WSTIF,PSTIF,QSTIF,HSTIF,      !ELMAT
     .                                                      ITASK)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 32 :
C
C     2D/3D ELEMENT CONTACT/LINK ELEMENT
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  3   Evaluate Contact Stiffness Matrix
C       ITASK =  4   Evaluate Internal Contact Resisting Forces
C       ITASK =  9   Check correctness of integration rule
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
     .          EPMTX(*), RMAT1(*), EMASS(*), VNORI(*)
      DIMENSION STRA0(*), STRS0(*), TEMPC(*)
      DIMENSION ELDIS(*), EHIST(*), STRAN(*), STRSG(*)
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*), PSTIF(*), QSTIF(*),
     .          HSTIF(*)
      DIMENSION WORK1(*)
C
C**** SELECT ELEMENT PROPERTIES
C
      NNODL=INT(PROEL( 2))
      NRULE=INT(PROEL( 3))
      NGAUL=INT(PROEL( 4))
      NTYPE=INT(PROEL( 6))
      THICK=    PROEL( 7)
C
      RITEN=PROEL( 8)
      RITEF=PROEL( 9)
      IF(ICONC.EQ.0) THEN
       TRUPL=PROEL(10)
       TRUPM=PROEL(11)
      ENDIF
      TOLGA=PROEL(12)
      TOLGAM=PROEL(13)
C
      NDIML=2
      IF(NTYPE.EQ.5.OR.NTYPE.EQ.3) NDIML=1
C
      NOCOL=INT(PROEL(22))
C
      NNODN=NNODL
      IF(NOCOL.EQ.1) NNODN=INT(PROEL(21))
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10,  1, 30, 40,  1,  1,  1,  1, 90,  1,
     .         1,  1,  1,  1,  1,  1,  1,  1,  1), ITASK
    1 RETURN ! Nothing
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
   10 CONTINUE
      CALL SETI32(DVOLU,ELCOD,GPCOD,LNODS,PROPS,SHAPE,THICK,ELDIS,
     .            VNORI,
     .                             WORK1(ISETM(12)),       !       dvolu
     .            WORK1(ISETM(13)),WORK1(ISETM(14)),       ! gpcod,shape
     .            WORK1(ISETM( 1)),WORK1(ISETM( 2)),       ! deriv,posgp
     .            WORK1(ISETM( 3)),WORK1(ISETM( 4)),       ! weigp,xjacm
     .            WORK1(ISETM( 9)),WORK1(ISETM(10)),       ! VNORL,ELCO1
     .            WORK1(ISETM(17)),                        ! ELDI1
     .            STRSG,
     .            WORK1(ISETM(24)),WORK1(ISETM(25)),       ! DVOL2,ELCO2
     .            WORK1(ISETM(26)),WORK1(ISETM(27)),       ! GPCO2,SHAP2
     .            WORK1(ISETM(28)),WORK1(ISETM(29)),       ! DISI2,DISP2
     .            WORK1(ISETM(30)),                        ! VTANL
     .            ITASK)
      RETURN
C
C**** EVALUATE CONTACT STIFFNESS MATRIX
C
   30 CONTINUE
C
      CALL SETI32(DVOLU,ELCOD,GPCOD,LNODS,PROPS,SHAPE,THICK,ELDIS,
     .            VNORI,
     .                             WORK1(ISTIF(11)),       !       dvolu
     .            WORK1(ISTIF(12)),WORK1(ISTIF(13)),       ! gpcod,shape
     .            WORK1(ISTIF(15)),WORK1(ISTIF(16)),       ! deriv,posgp
     .            WORK1(ISTIF(17)),WORK1(ISTIF(18)),       ! weigp,xjacm
     .            WORK1(ISTIF( 8)),WORK1(ISTIF( 9)),       ! VNORL,ELCO1
     .            WORK1(ISTIF(27)),                        ! ELDI1
     .            STRSG,
     .            WORK1(ISTIF(36)),WORK1(ISTIF(37)),       ! DVOL2,ELCO2
     .            WORK1(ISTIF(38)),WORK1(ISTIF(39)),       ! GPCO2,SHAP2
     .            WORK1(ISTIF(40)),WORK1(ISTIF( 7)),       ! DISI2,DISPL
     .            WORK1(ISTIF(41)),                        ! VTANL
     .            ITASK)
C
      CALL ASIL00(WORK1(ISTIF(25)),WORK1(ISTIF(26)),       ! ESTIF,WSTIF
     .            NKOVA,KDYNA,
     .            WORK1(ISTIF(35)),                        ! FMATX
     .            NEVAB,    1)
C
      IF(NPOIC.EQ.0) THEN
       CALL MASS32(WORK1(ISTIF(11)),                       ! DVOLU
     .             PROPS,LNODS,
     .             WORK1(ISTIF(13)),WORK1(ISTIF(25)),      ! SHAPE,ESTIF
     .             EHIST,
     .             WORK1(ISTIF(27)),                       ! ELDIS
     .             WORK1(ISTIF( 5)),WORK1(ISTIF( 7)),      ! TENOD,DISPL
     .             WORK1(ISTIF( 9)),                       ! ELCOD
     .             WORK1(ISTIF( 8)),WORK1(ISTIF( 2)),      ! VNORL,DMATX
     .             WORK1(ISTIF(34)),WORK1(ISTIF(35)),      ! EMATX,FMATX
     .             WORK1(ISTIF( 1)),WORK1(ISTIF(41)))      ! BMATX,VTANL
      ELSE
       IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
c       CALL MASS32B(WORK1(ISTIF(11)),                     ! DVOLU
        CALL MASS32Z(WORK1(ISTIF(11)),                     ! DVOLU
     .             PROPS,
     .             WORK1(ISTIF(13)),WORK1(ISTIF(25)),      ! SHAPE,ESTIF
     .             EHIST,
     .             WORK1(ISTIF(27)),                       ! ELDIS
     .             WORK1(ISTIF( 5)),WORK1(ISTIF( 7)),      ! TENOD,DISPL
     .             WORK1(ISTIF( 9)),                       ! ELCOD
     .             WORK1(ISTIF( 8)),WORK1(ISTIF( 2)))      ! VNORL,DMATX
       ENDIF
       IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
        CALL RUNEND('ERROR: IAUGM=2,3 WITH ITYPE=32')
       ENDIF
      ENDIF
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
      CALL SETI32(DVOLU,ELCOD,GPCOD,LNODS,PROPS,SHAPE,THICK,ELDIS,
     .            VNORI,
     .                             WORK1(IFORC(22)),       !       dvolu
     .            WORK1(IFORC(23)),WORK1(IFORC(24)),       ! gpcod,shape
     .            WORK1(IFORC(26)),WORK1(IFORC(27)),       ! deriv,posgp
     .            WORK1(IFORC(28)),WORK1(IFORC(29)),       ! weigp,xjacm
     .            WORK1(IFORC(18)),WORK1(IFORC(19)),       ! VNORL,ELCO1
     .            WORK1(IFORC(37)),                        ! ELDI1
     .            STRSG,
     .            WORK1(IFORC(42)),WORK1(IFORC(43)),       ! DVOL2,ELCO2
     .            WORK1(IFORC(44)),WORK1(IFORC(45)),       ! GPCO2,SHAP2
     .            WORK1(IFORC(38)),WORK1(IFORC(15)),       ! DISI2,DISPL
     .            WORK1(IFORC(47)),                        ! VTANL
     .            ITASK)
C
      IF(NPOIC.EQ.0) THEN
       CALL FRIN32(PROPS,LNODS,
     .             WORK1(IFORC(37)),                       ! ELDIS
     .             WORK1(IFORC( 1)),                       ! BMSIG
     .             WORK1(IFORC(22)),WORK1(IFORC(24)),      ! DVOLU,SHAPE
     .             EHIST,
     .             WORK1(IFORC(11)),                       ! TENOD
     .             WORK1(IFORC(19)),                       ! ELCOD
     .             WORK1(IFORC(14)),WORK1(IFORC(16)),      ! PWOEL,PREAL
     .             WORK1(IFORC(17)),WORK1(IFORC(18)),      ! TGAPL,VNORL
     .             WORK1(IFORC(47)),                       ! VTANL
     .             WORK1(IFORC(38)),WORK1(IFORC(15)))      ! DISI2,DISPL
      ELSE
       IF(IAUGM.EQ.0.OR.IAUGM.EQ.1) THEN
c       CALL FRIN32B(PROPS,
        CALL FRIN32Z(PROPS,
     .             WORK1(IFORC(37)),                       ! ELDIS
     .             WORK1(IFORC( 1)),                       ! BMSIG
     .             WORK1(IFORC(22)),WORK1(IFORC(24)),      ! DVOLU,SHAPE
     .             EHIST,
     .             WORK1(IFORC(11)),                       ! TENOD
     .             WORK1(IFORC(19)),                       ! ELCOD
     .             WORK1(IFORC(14)),WORK1(IFORC(16)),      ! PWOEL,PREAL
     .             WORK1(IFORC(17)),WORK1(IFORC(18)))      ! TGAPL,VNORL
       ENDIF
       IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
        CALL RUNEND('ERROR: IAUGM=2,3 WITH ITYPE=32')
       ENDIF
      ENDIF
      RETURN
C
C**** CHECK CORRECTNESS OF INTEGRATION RULE
C
   90 CONTINUE
      NNOBO=NNODL/2
      IF(NPOIC.GT.0) THEN
       IF(NNODC.GT.0) THEN
        NNOBO=NNODL/3
       ENDIF
      ENDIF
      IF(NOCOL.EQ.1) NNOBO=NNODL                   ! non-coincident-mesh
      CALL CHEK01(NDIML,IELEM,NNOBO,NRULE,NGAUL,    1)
      RETURN
C
      END
