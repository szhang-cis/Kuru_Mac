      SUBROUTINE ELMLIBT(LNODST,PROELT,PROPST,
     .                   ELDATT,ELPRET,ELVART,ELMATT,HTLODT,
     .                   WORK1T,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE CALLS THE APPROPRIATE ELEMENT PROCESSOR
C     TO PERFORM THE OPERATION REQUIRED BY THE VALUE GIVEN TO ITASK
C
C     ELEMENT TYPES:
C
C       LTYPE =   5  Thermal elements (1D/2D/3D) with phase-change
C       LTYPE = 101  Thermal Boundary element
C       LTYPE = 104  Thermal Interaction (gap) element
C
C     FUNCTIONS:
C
C       ITASK =  1   Perform some needed initial computations
C       ITASK =  2   -------
C       ITASK =  3   Evaluate : Conductivity Matrix
C                               Heat Capacity Matrix (only for
C                               Transient Problems)
C                               Phase Change Matrix
C       ITASK =  4   Evaluate Internal Resisting Heat
C       ITASK =  5   Evaluate Internal Resisting Heat
C                    (alternative option for micro advective effects)
C       ITASK =  6   Evaluate : Internal Dynamic Resisting Heat
C                               (only for Transient Problems)
C                               Internal Dynamic Phase Change Heat
C       ITASK =  7   Evaluate equivalent Volume Heat (Internal Heat)
C       ITASK =  8   -------
C       ITASK =  9   Check correctness of integration rule
C       ITASK = 10   -------
C       ITASK = 11   -------
C       ITASK = 12   Output Gaussian variables
C       ITASK = 13   Evaluate contribution for nodal heat fluxes
C       ITASK = 14   -------
C       ITASK = 15   Evaluate contribution for nodal internal variables
C       ITASK = 16   Evaluate variables for outpost.f
C       ITASK = 17   Evaluate contribution for nodal porosity criteria
C       ITASK = 18   Nothing (see outsmot.f)
C       ITASK = 19   Evaluate contribution for nodal L*f_pc
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
C
      DIMENSION LNODST(*), PROELT(*), PROPST(*)
      DIMENSION ELDATT(*), ELPRET(*), ELVART(*), ELMATT(*), HTLODT(*)
      DIMENSION WORK1T(*)
C
C**** ZERO RELEVANT ARRAYS (Send WORK1 as real and integer array)
C
      CALL ELM000T(WORK1T,WORK1T,
     .             ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
     .             ELMATT(IMATXT(3)),ELMATT(IMATXT(4)),
     .             ELMATT(IMATXT(5)),ELMATT(IMATXT(6)),
     .                                                        ITASKT)
C
C**** CALL APPROPRIATE ELEMENT PROCESSOR
C
      LTYPET=INT(PROELT(5))
C
      IF(LTYPET.LE.100) THEN 
       GOTO (99,99,99,99, 5) LTYPET
      ELSE
       GOTO (101, 99, 99,104) LTYPET-100
      END IF
   99 RETURN ! Not implemented yet
C
    5 CALL ELM005T(LNODST,PROELT,PROPST,WORK1T,
     .             ELDATT(IDATAT( 1)),ELDATT(IDATAT( 2)),
     .             ELDATT(IDATAT( 3)),ELDATT(IDATAT( 4)),
     .             ELDATT(IDATAT( 5)),ELDATT(IDATAT( 6)),
     .             ELDATT(IDATAT( 7)),ELDATT(IDATAT( 8)),
     .             ELPRET(IPREVT( 1)),ELPRET(IPREVT( 2)),
     .             ELPRET(IPREVT( 3)),
     .             ELVART(ISTATT( 1)),ELVART(ISTATT( 2)),
     .             ELVART(ISTATT( 3)),
     .             ELVART(ISTATT( 4)),                   
     .             ELMATT(IMATXT( 1)),ELMATT(IMATXT( 2)),
     .             ELMATT(IMATXT( 3)),
     .             ELMATT(IMATXT( 4)),ELMATT(IMATXT( 5)),
     .             ELMATT(IMATXT( 6)),
     .                                                        ITASKT)
      RETURN
C
  101 CALL ELM101T(LNODST,PROELT,PROPST,WORK1T,
     .             ELDATT(IDATAT(1)),ELDATT(IDATAT(2)),
     .             ELDATT(IDATAT(3)),
     .             ELDATT(IDATAT(4)),ELDATT(IDATAT(5)),
     .             ELDATT(IDATAT(6)),
     .             ELDATT(IDATAT(7)),ELDATT(IDATAT(8)),
     .             ELPRET(IPREVT(1)),ELPRET(IPREVT(2)),
     .             ELPRET(IPREVT(3)),
     .             ELVART(ISTATT(1)),ELVART(ISTATT(2)),
     .             ELVART(ISTATT(3)),
     .             ELVART(ISTATT(4)),                   
     .             ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
     .             ELMATT(IMATXT(3)),
     .             ELMATT(IMATXT(4)),ELMATT(IMATXT(5)),
     .             ELMATT(IMATXT(6)),HTLODT,
     .                                                        ITASKT)
      RETURN
C  
  104 CALL ELM104T(LNODST,PROELT,PROPST,WORK1T,
     .             ELDATT(IDATAT(1)),ELDATT(IDATAT(2)),
     .             ELDATT(IDATAT(3)),
     .             ELDATT(IDATAT(4)),ELDATT(IDATAT(5)),
     .             ELDATT(IDATAT(6)),
     .             ELDATT(IDATAT(7)),ELDATT(IDATAT(8)),
     .             ELDATT(IDATAT(10)),
     .             ELPRET(IPREVT(1)),ELPRET(IPREVT(2)),
     .             ELPRET(IPREVT(3)),
     .             ELVART(ISTATT(1)),ELVART(ISTATT(2)),
     .             ELVART(ISTATT(3)),
     .             ELVART(ISTATT(4)),                   
     .             ELMATT(IMATXT(1)),ELMATT(IMATXT(2)),
     .             ELMATT(IMATXT(3)),
     .             ELMATT(IMATXT(4)),ELMATT(IMATXT(5)),
     .             ELMATT(IMATXT(6)),
     .                                                        ITASKT)
      RETURN
C  
      END
