      SUBROUTINE ELMLIBS(LNODSS,PROELS,PROPSS,
     .                   ELDATS,ELPRES,ELVARS,ELMATS,
     .                   WORK1S,                        ITASKS)
C***********************************************************************
C
C**** THIS ROUTINE CALLS THE APPROPRIATE ELEMENT PROCESSOR
C     TO PERFORM THE OPERATION REQUIRED BY THE VALUE GIVEN TO ITASK
C
C     ELEMENT TYPES:
C
C       LTYPE =   5     Continuum elements ( 1D/2D/3D )
C       LTYPE = 101     Boundary element ( not used yet )
C       LTYPE = 104     Interaction (gap) element ( not used yet )
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   -------
C       ITASK =  3   Evaluate : Conductivity Matrix
C                               Heat Capacity Matrix (only for
C                               Transient Problems)
C       ITASK =  4   Evaluate Internal Resisting Heat
C       ITASK =  5   -------
C       ITASK =  6   Evaluate : Internal Dynamic Resisting Heat
C                               (only for Transient Problems)
C                               Internal Dynamic Phase Change Heat
C       ITASK =  7   Evaluate equivalent Volume Heat (Internal Heat)
C       ITASK =  8   ------- (Evaluate equivalent Surface Forces)
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   ------- (Increment Non-Tensional Strains)
C       ITASK = 11   ------- (Read Initial State Variables at Gauss P.)
C       ITASK = 12   Output Gaussian variables
C       ITASK = 13   Evaluate contribution for nodal heat fluxes
C       ITASK = 14   ------- (Evaluate contribution for nodal strains)
C       ITASK = 15   Evaluate contribution for nodal internal variables
C       ITASK = 16   Evaluate variables for outpost.f
C       ITASK = 17   Evaluate contribution for nodal porosity criteria
C       ITASK = 18   Nothing (see outsmot.f)
C       ITASK = 19   -------
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
C
      DIMENSION LNODSS(*), PROELS(*), PROPSS(*)
      DIMENSION ELDATS(*), ELPRES(*), ELVARS(*), ELMATS(*)
      DIMENSION WORK1S(*)
C
C**** ZERO RELEVANT ARRAYS (Send WORK1 as real and integer array)
C
      CALL ELM000S(WORK1S,WORK1S,
     .             ELMATS(IMATXS(1)),ELMATS(IMATXS(2)),
     .             ELMATS(IMATXS(3)),ELMATS(IMATXS(4)),
     .             ELMATS(IMATXS(5)),ELMATS(IMATXS(6)),
     .                                                 ITASKS)
C
C**** CALL APPROPRIATE ELEMENT PROCESSOR
C
      LTYPES=INT(PROELS(5))
C
      IF(LTYPES.LE.100) THEN 
       GOTO (99,99,99,99, 5) LTYPES
      ELSE
       GOTO (101, 99, 99,104) LTYPES-100
      END IF
   99 RETURN ! Not implemented yet
C
    5 CALL ELM005S(LNODSS,PROELS,PROPSS,WORK1S,
     .             ELDATS(IDATAS( 1)),ELDATS(IDATAS( 2)),
     .             ELDATS(IDATAS( 3)),ELDATS(IDATAS( 4)),
     .             ELDATS(IDATAS( 5)),ELDATS(IDATAS( 6)),
     .             ELDATS(IDATAS( 7)),ELDATS(IDATAS( 8)),
     .             ELPRES(IPREVS( 1)),ELPRES(IPREVS( 2)),
     .             ELPRES(IPREVS( 3)),
     .             ELVARS(ISTATS( 1)),ELVARS(ISTATS( 2)),
     .             ELVARS(ISTATS( 3)),
     .             ELVARS(ISTATS( 4)),                   
     .             ELMATS(IMATXS( 1)),ELMATS(IMATXS( 2)),
     .             ELMATS(IMATXS( 3)),
     .             ELMATS(IMATXS( 4)),ELMATS(IMATXS( 5)),
     .             ELMATS(IMATXS( 6)),
     .                                                        ITASKS)
      RETURN
C
  101 call runends('elm101s not implemented')
c 101 CALL ELM101T(LNODST,PROELT,PROPST,WORK1T,
c    .             ELDATT(IDATAT(1)),ELDATT(IDATAT(2)),
c    .             ELDATT(IDATAT(3)),
c    .             ELDATT(IDATAT(4)),ELDATT(IDATAT(5)),
c    .             ELDATT(IDATAT(6)),
c    .             ELDATT(IDATAT(7)),ELDATT(IDATAT(8)),
c    .             ELPRET(IPREVT(1)),ELPRET(IPREVT(2)),
c    .             ELPRET(IPREVT(3)),
c    .             ELVART(ISTATT(1)),ELVART(ISTATT(2)),
c    .             ELVART(ISTATT(3)),
c    .             ELVART(ISTATT(4)),                   
c    .             ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
c    .             ELMATT(IMATXT(3)),
c    .             ELMATT(IMATXT(4)),ELMATT(IMATXT(5)),
c    .             ELMATT(IMATXT(6)),
c    .                                                      ITASKT)
      RETURN
C
  104 call runends('elm101s not implemented')
c 104 CALL ELM104T(LNODST,PROELT,PROPST,WORK1T,
c    .             ELDATT(IDATAT(1)),ELDATT(IDATAT(2)),
c    .             ELDATT(IDATAT(3)),
c    .             ELDATT(IDATAT(4)),ELDATT(IDATAT(5)),
c    .             ELDATT(IDATAT(6)),
c    .             ELDATT(IDATAT(7)),ELDATT(IDATAT(8)),
c    .             ELDATT(IDATAT(10)),
c    .             ELPRET(IPREVT(1)),ELPRET(IPREVT(2)),
c    .             ELPRET(IPREVT(3)),
c    .             ELVART(ISTATT(1)),ELVART(ISTATT(2)),
c    .             ELVART(ISTATT(3)),
c    .             ELVART(ISTATT(4)),                   
c    .             ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
c    .             ELMATT(IMATXT(3)),
c    .             ELMATT(IMATXT(4)),ELMATT(IMATXT(5)),
c    .             ELMATT(IMATXT(6)),
c    .                                                      ITASKT)
      RETURN
C
      END
