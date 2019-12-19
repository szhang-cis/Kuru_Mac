      SUBROUTINE ELM000S(WORK1S,NWORIS,
     .                   CSTIFS,ESTIFS,WSTIFS,HSTIFS,PSTIFS,
     .                   QSTIFS,                                  !ELMAT
     .                                                   ITASKS)
C***********************************************************************
C
C**** THIS ROUTINE ZEROES RELEVANT ARRAYS FOR DIFFERENT ITASK VALUES 
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
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION CSTIFS(*), ESTIFS(*), WSTIFS(*), HSTIFS(*), PSTIFS(*),
     .          QSTIFS(*)
      DIMENSION WORK1S(*)
      DIMENSION NWORIS(*)
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO (  1, 20, 30, 40,  1, 40, 40, 40,  1,  1,
     .         1,  1,130,140,150,  1,170,  1,  1), ITASKS
    1 RETURN ! Nothing
C
C**** ZERO MASS MATRIX
C
   20 CONTINUE       ! set to zero in elm...
c     DO IKOVAS=1,NKOVAS
c      WSTIFS(IKOVAS)=0.0
c     END DO
      RETURN
C
C**** ZERO STIFFNESS ( & COUPLING MATRICES )
C
   30 CONTINUE       ! set to zero in elm...
c     DO IKOVAS=1,NKOVAS
c      ESTIFS(IKOVAS)=0.0
c      IF(KPROBS.EQ.4) WSTIFS(IKOVAS)=0.0 !Only useful for thermal pr.
c     END DO
c     IF(KFLAGS.EQ.1)THEN
c      DO INDEXS=1,NNODES*NEVABS
c       HSTIFS(INDEXS)=0.0
c      END DO
c     ENDIF
      RETURN
C
C**** ZERO FORCES
C
   40 DO IEVABS=1,NEVABS
       WORK1S(IEVABS)=0.0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL HEAT FLUXES
C
  130 DO INDEXS=IGSMOS(4),IGSMOS(6)-1
       WORK1S(INDEXS)=0.0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL TEMP. GRADIENTS
C
  140 DO INDEXS=IGSMOS(6),IGSMOS(8)-1
       WORK1S(INDEXS)=0.0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL INTERNAL VARIABLES
C
  150 DO INDEXS=IGSMOS(9),IGSMOS(11)-1
       WORK1S(INDEXS)=0.0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL POROSITY CRITERIA
C
  170 DO INDEXS=IGSMOS(23),IGSMOS(25)-1
       WORK1S(INDEXS)=0.0
      END DO
      RETURN
C
      END
