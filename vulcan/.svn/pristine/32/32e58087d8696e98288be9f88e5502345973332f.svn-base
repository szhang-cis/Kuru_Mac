      SUBROUTINE ELMLIB(LNODS,PROEL,PROPS,
     .                  INFRI,COFRI,NOPRF,PRESF,VANIS,
     .                  ELDAT,ELPRE,ELVAR,ELMAT,
     .                  WORK1,                     ITASK)
C***********************************************************************
C
C**** THIS ROUTINE CALLS THE APPROPRIATE ELEMENT PROCESSOR
C     TO PERFORM THE OPERATION REQUIRED BY THE VALUE GIVEN TO ITASK
C
C     ELEMENT TYPES:
C
C       LTYPE =   1     Continuum elements ( 1D/2D/3D )  (not used now!)
C       LTYPE =   2     Truss elements ( 1D/2D/3D )      (not used now!)
C       LTYPE =   3     Nodal linking or free surface elements ( 2D/3D )
C       LTYPE =   4     Nodal contact elements ( 2D/3D )
C       LTYPE =  30     B-bar for Continuum elements ( 1D/2D/3D )
C       LTYPE =  31     Contact elements ( 2D only )     (not used now!)
C       LTYPE =  32     Elemental contact/linking elements ( 2D/3D )
C       LTYPE =  33     Discontinuous Galerkin (DG) elements ( 2D/3D )
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   Evaluate Mass Matrix
C       ITASK =  3   Evaluate Stiffness ( & Coupling ) Matrices
C       ITASK =  4   Evaluate Internal Resisting Forces
C       ITASK =  5   ------
C       ITASK =  6   Evaluate Internal Resisting Dynamic Forces
C       ITASK =  7   Evaluate equivalent Volume Forces
C       ITASK =  8   Evaluate equivalent Surface Forces
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   ------ (Increment Non-Tensional Strains) 
C       ITASK = 11   ------ (Read Initial State Variables at Gauss P.)
C       ITASK = 12   Output Gaussian Variables
C       ITASK = 13   Evaluate contribution for Nodal stresses
C       ITASK = 14   Evaluate contribution for Nodal strains
C       ITASK = 15   Evaluate contribution for Nodal internal variables
C       ITASK = 16   Evaluate variables for outpos.f
C       ITASK = 17   ------
C       ITASK = 18   Nothing (see outsmo.f)
C       ITASK = 19   ------
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
C
      DIMENSION LNODS(*), PROEL(*), PROPS(*)
      DIMENSION INFRI(*), COFRI(*), NOPRF(*), PRESF(*), VANIS(*)
      DIMENSION ELDAT(*), ELPRE(*), ELVAR(*), ELMAT(*)
      DIMENSION WORK1(*)
C
C**** ZERO RELEVANT ARRAYS (Send WORK1 as real and integer array)
C
      CALL ELM000(WORK1,WORK1,
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                       ITASK)
C
C**** CALL APPROPRIATE ELEMENT PROCESSOR
C
      LTYPE=INT(PROEL(5))
C
      GOTO ( 1, 2, 3, 4,99,99,99,99,99,99,
     .      99,99,99,99,99,99,99,99,99,99,
     .      99,99,99,99,99,99,99,99,99,30,
     .      31,32,33                      ) LTYPE
C
   99 RETURN ! Not implemented yet
C
    1 CALL ELM001(LNODS,PROEL,PROPS,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
    2 CALL ELM002(LNODS,PROEL,PROPS,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
    3 CALL ELM003(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
    4 CALL ELM004(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
   30 CALL ELM030(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,VANIS,
     .            WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),ELDAT(IDATA(9)),
     .            ELDAT(IDATA(10)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
   31 CALL ELM031(LNODS,PROEL,PROPS,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
   32 CALL ELM032(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),ELDAT(IDATA(11)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
   33 CALL ELM033(LNODS,PROEL,PROPS,INFRI,COFRI,NOPRF,PRESF,WORK1,
     .            ELDAT(IDATA(1)),ELDAT(IDATA(2)),ELDAT(IDATA(3)),
     .            ELDAT(IDATA(4)),ELDAT(IDATA(5)),ELDAT(IDATA(6)),
     .            ELDAT(IDATA(7)),ELDAT(IDATA(8)),ELDAT(IDATA(11)),
     .            ELPRE(IPREV(1)),ELPRE(IPREV(2)),ELPRE(IPREV(3)),
     .            ELVAR(ISTAT(1)),ELVAR(ISTAT(2)),ELVAR(ISTAT(3)),
     .            ELVAR(ISTAT(4)),
     .            ELMAT(IMATX(1)),ELMAT(IMATX(2)),ELMAT(IMATX(3)),
     .            ELMAT(IMATX(4)),ELMAT(IMATX(5)),ELMAT(IMATX(6)),
     .                                                      ITASK)
      RETURN
C
      END
