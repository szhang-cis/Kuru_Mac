      SUBROUTINE ELM003(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,
     .                  EMASS,                                   !ELDAT
     .                  STRA0,STRS0,TEMPC,                       !ELPRE
     .                  ELDIS,EHIST,STRAN,STRSG,                 !ELVAR
     .                  CSTIF,ESTIF,WSTIF,PSTIF,QSTIF,HSTIF,     !ELMAT
     .                                                      ITASK)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 3 :
C
C     1D/2D/3D LINK ELEMENT
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
      NTYPE=INT(PROEL( 6))
      IF(NNODL.NE.2) CALL RUNEND('ELM003: NNODL MUST BE 2 FOR LINKS  ')
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO (  1,  1, 30, 40,  1,  1,  1,  1,  1,  1,
     .         1,  1,  1,  1,  1,  1,  1,  1,  1), ITASK
    1 RETURN ! Nothing
C
C**** EVALUATE STIFFNESS
C
   30 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      CALL SETI03(ELCOD,ELDIS,
     .            WORK1(ISTIF( 9)),                        ! ELCOD
     .            WORK1(ISTIF(27)))                        ! ELDI1
C
      CALL ASIL00(WORK1(ISTIF(25)),WORK1(ISTIF(26)),       ! ESTIF,WSTIF
     .            NKOVA,KDYNA,
     .            WORK1(ISTIF(34)),                        ! EMATX
     .            NEVAB,    1)
C
      CALL STIF03(PROPS,
     .            WORK1(ISTIF(25)),                        ! ESTIF
     .            WORK1(ISTIF( 5)))                        ! TENOD
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
      RETURN
C
C**** EVALUATE INTERNAL RESISTING FORCES
C
   40 CONTINUE
C
C**** SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
      CALL SETI03(ELCOD,ELDIS,
     .            WORK1(IFORC(19)),                        ! ELCOD
     .            WORK1(IFORC(37)))                        ! ELDI1
C
      CALL FRIN03(WORK1(IFORC(37)),                        ! ELDIS
     .            PROPS,
     .            WORK1(IFORC(1)),WORK1(IFORC(11)),        ! BMSIG,TENOD
     .            WORK1(IFORC(14)))                        ! PWOEL
      RETURN
C
      END
