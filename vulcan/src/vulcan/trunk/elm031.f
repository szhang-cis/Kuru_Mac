      SUBROUTINE ELM031(LNODS,PROEL,PROPS,WORK1,
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
C     2D (only) CONTACT-FRICTION ELEMENT
C
C     FUNCTIONS:
C
C       ITASK =  3   Evaluate Stiffness Matrix
C       ITASK =  4   Evaluate Internal Resisting Forces
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION LNODS(*), PROEL(*), PROPS(*)
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
   30 CALL STIF31(PROPS,ESTIF,WORK1(ISTIF(5)),ELDIS,EHIST)
      RETURN
C
C**** EVALUATE INTERNAL RESISTING FORCES
C
   40 CALL FRIN31(ELDIS,PROPS,WORK1(IFORC(1)), WORK1(IFORC(11)),
     .            EHIST,ELCOD,WORK1(IFORC(14)),WORK1(IFORC(15)),
     .                        WORK1(IFORC(16)),WORK1(IFORC(17)))
      RETURN
C
      END
