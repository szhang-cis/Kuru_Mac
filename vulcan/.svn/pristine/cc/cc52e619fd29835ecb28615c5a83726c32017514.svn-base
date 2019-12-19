      SUBROUTINE ELM002(LNODS,PROEL,PROPS,WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,EMASS,
     .                  STRA0,STRS0,TEMPC,
     .                  ELDIS,EHIST,STRAN,STRSG,
     .                  CSTIF,ESTIF,WSTIF,HSTIF,PSTIF,QSTIF,
     .                                                            ITASK)
C************************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 2 :
C
C    2D/3D LINEAR ELASTIC ELEMENT
C
C    FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   Evaluate Mass Matrix
C       ITASK =  3   Evaluate Stiffness ( & Coupling ) Matrices
C       ITASK =  4   Evaluate Internal Resisting Forces
C
C    FOR COUPLED ELEMENTS:
C
C       ITASK = 21   Evaluate Permeability, Compressib. & Coupling Mx.
C       ITASK = 22   Smoothing                           ( # )
C
C********************************************************************
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
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*), HSTIF(*), PSTIF(*),
     .          QSTIF(*)
      DIMENSION WORK1(*)
C
C***SELECT ELEMENT PROPERTIES
C
      NNODL=INT(PROEL( 2))
      NRULE=INT(PROEL( 3))
      NGAUL=INT(PROEL( 4))
      NTYPE=INT(PROEL( 6))
      NCRIT=0
C
      NSTRE=1
      NSTRS=1
      NKOST=1
C
C***DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO (  1,  1, 30, 40), ITASK
    1 RETURN ! Nothing
C
C***EVALUATE STIFFNESS ( & COUPLING MATRICES )
C
   30 CONTINUE
      CALL STIF02(EHIST,ELCOD,PROPS,ESTIF,HSTIF)
      RETURN
C
C***EVALUATE INTERNAL RESISTING FORCES
C
   40 CONTINUE
      CALL FRIN02(EHIST,ELCOD,ELDIS,PROPS,
     .            STRAN,STRSG,STRA0,STRS0,WORK1)
      RETURN
C
      END
