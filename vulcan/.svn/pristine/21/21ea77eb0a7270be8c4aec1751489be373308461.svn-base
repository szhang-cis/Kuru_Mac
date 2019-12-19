      SUBROUTINE ELM001(LNODS,PROEL,PROPS,WORK1,
     .                  ELCOD,CARTD,DVOLU,GPCOD,SHAPE,EPMTX,RMAT1,EMASS,
     .                  STRA0,STRS0,TEMPC,
     .                  ELDIS,EHIST,STRAN,STRSG,
     .                  CSTIF,ESTIF,WSTIF,HSTIF,PSTIF,QSTIF,
     .                                                            ITASK)
C************************************************************************
C
C**** THIS ROUTINE PERFORMS COMPUTATIONS AT ELEMENT LEVEL FOR
C     ELEMENT NO. 1 :
C
C     1D/2D/3D ELASTO/PLASTIC ELEMENT WITH/WITHOUT COUPLED
C     PORE FLUID FIELD
C
C    FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   Evaluate Mass Matrix
C       ITASK =  3   Evaluate Stiffness ( & Coupling ) Matrices
C       ITASK =  4   Evaluate Internal Resisting Forces
C       ITASK =  6   Evaluate equivalent Gravity Forces
C       ITASK =  7   Evaluate equivalent Surface Forces
C       ITASK =  8   Check correctness of integration rule
C       ITASK =  9   Increment Non-Tensional Strains 
C       ITASK = 10   Read Initial State Variables at Gauss Points
C       ITASK = 11   Output Gaussian Variables
C       ITASK = 12   Evaluate contribution for Nodal stresses
C       ITASK = 13   Evaluate contribution for Nodal strains in cracks
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
      DIMENSION NDATO(2,5)
      DATA NDATO/3,3,   ! 2D PLANE STRESS
     .           3,4,   ! 2D PLANE STRAIN
     .           4,4,   ! 2D AXISYMMETRIC
     .           6,6,   ! 3D
     .           1,1/   ! 1D
C
C***SELECT ELEMENT PROPERTIES
C
      NNODL=INT(PROEL( 2))
      NRULE=INT(PROEL( 3))
      NGAUL=INT(PROEL( 4))
      NTYPE=INT(PROEL( 6))
      THICK=    PROEL( 7)
C
      NCRIT=INT(PROPS(36))
C
      NSTRE=NDATO(1,NTYPE)
      NSTRS=NDATO(2,NTYPE)
      NKOST=(NSTRS-1)*NSTRS/2+NSTRS
C
C***DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO ( 10, 20, 30, 40,  1, 60, 70, 80, 90,100,
     .       110,120,130,  1,  1,  1,  1,  1,  1,  1,
     .         1                                     ), ITASK
    1 RETURN ! Nothing
C
C***SET UP SOME CONSTANT MATRICES AND INTEGRATION VALUES
C
   10 CONTINUE
      CALL SETM01(CARTD,DVOLU,ELCOD,EMASS,EPMTX,GPCOD,LNODS,PROPS,
     .            RMAT1,SHAPE,THICK,
     .            WORK1(ISETM(1)),WORK1(ISETM(2)),WORK1(ISETM(3)),
     .            WORK1(ISETM(4)),WORK1(ISETM(5)))
      RETURN
C
C***EVALUATE MASS MATRIX
C
   20 CONTINUE
      CALL MASS01(DVOLU,PROPS,SHAPE,WSTIF)
      RETURN
C
C***EVALUATE STIFFNESS ( & COUPLING MATRICES )
C
   30 CONTINUE
      CALL STIF01(CARTD,DVOLU,EHIST,ELCOD,ELDIS,EPMTX,GPCOD,
     .            PROPS,SHAPE,STRSG,ESTIF,HSTIF,EMASS,
     .            WORK1(ISTIF(1)),WORK1(ISTIF(2)),WORK1(ISTIF(3)),
     .            WORK1(ISTIF(4)),WORK1(ISTIF(6)),STRAN)
      RETURN
C
C***EVALUATE INTERNAL RESISTING FORCES
C
   40 CONTINUE
      CALL FRIN01(CARTD,DVOLU,EHIST,ELCOD,ELDIS,EMASS,EPMTX,GPCOD,
     .            LNODS,PROPS,RMAT1,SHAPE,STRAN,STRSG,STRA0,STRS0,
     .            THICK,
     .            WORK1(IFORC( 1)),WORK1(IFORC( 2)),WORK1(IFORC( 3)),
     .            WORK1(IFORC( 4)),WORK1(IFORC( 5)),WORK1(IFORC( 6)),
     .            WORK1(IFORC( 7)),WORK1(IFORC( 8)),WORK1(IFORC( 9)),
     .            WORK1(IFORC(10)),WORK1(IFORC(11)),WORK1(IFORC(12)),
     .            WORK1(IFORC(14)),WORK1(IFORC(16)),WORK1(IFORC(17)))
      RETURN
C
C***EVALUATE EQUIVALENT GRAVITY LOADS
C
   60 CONTINUE
      CALL LDGR01(DVOLU,PROPS,SHAPE,
     .            WORK1)
      RETURN
C
C***EVALUATE EQUIVALENT SURFACE LOADS
C
   70 CONTINUE
      CALL LDSF01(ELCOD,LNODS,PROPS,THICK,
     .            WORK1(ILOSU( 1)),WORK1(ILOSU( 2)),WORK1(ILOSU( 3)),
     .                             WORK1(ILOSU( 5)),WORK1(ILOSU( 6)),
     .            WORK1(ILOSU( 7)),WORK1(ILOSU( 8)),WORK1(ILOSU( 9)),
     .            WORK1(ILOSU(10)),                 WORK1(ILOSU(12)),
     .            WORK1(ILOSU(13)))
      RETURN
C
C***CHECK CORRECTNESS OF INTEGRATION RULES
C
   80 CONTINUE
      CALL CHEK01(NDIME,IELEM,NNODL,NRULE,NGAUL)
      RETURN
C
C***INCREMENT NON-TENSIONAL STRAINS
C
   90 CONTINUE
      CALL INCS01(EHIST,PROPS,STRA0,TEMPC)
      RETURN
C
C***READ INITIAL STATE VARIABLES AT GAUSS POINTS
C
  100 CONTINUE
      CALL PRGA01(PROPS,EHIST,STRA0,STRS0,STRSG)
      RETURN
C
C***OUTPUT GAUSSIAN STATE VARIABLES
C
  110 CONTINUE
      IF(LARGE.EQ.1) CALL OUTCAU(CARTD,ELDIS,STRSG,WORK1)
      CALL OUTG01(EHIST,STRAN,STRSG)
      RETURN
C
C***COMPUTE ELEMENTAL CONTRIBUTION FOR NODAL STRESSES
C
  120 CONTINUE
      IPUN1=1
      IPUN2=IPUN1+NSTR1*NNODE
      CALL SMOG01(DVOLU,EMASS,SHAPE,STRSG,
     .            WORK1(IPUN1),WORK1(IPUN2))
      RETURN
C
C***COMPUTE ELEMENTAL CONTRIBUTION FOR NODAL STRAINS IN CRACKS
C
  130 CONTINUE
      IPUN1=1
      IPUN2=IPUN1+NDIME*NNODE
      CALL SMOF01(DVOLU,EMASS,SHAPE,EHIST,
     .            WORK1(IPUN1),WORK1(IPUN2))
      RETURN
C
      END
