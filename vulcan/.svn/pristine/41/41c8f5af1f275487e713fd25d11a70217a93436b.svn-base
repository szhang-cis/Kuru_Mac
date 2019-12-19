      SUBROUTINE ELM000T(WORK1T,NWORIT,
     .                   CSTIFT,ESTIFT,WSTIFT,HSTIFT,PSTIFT,
     .                   QSTIFT,                                  !ELMAT
     .                                                   ITASKT)
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
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION CSTIFT(*), ESTIFT(*), WSTIFT(*), HSTIFT(*), PSTIFT(*),
     .          QSTIFT(*)
      DIMENSION WORK1T(*)
      DIMENSION NWORIT(*)
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO (  1, 20, 30, 40, 40, 40, 40, 40,  1,  1,
     .         1,  1,130,140,150,  1,170,  1,190), ITASKT
    1 RETURN ! Nothing
C
C**** ZERO MASS MATRIX
C
   20 CONTINUE       ! set to zero in elm...
c     DO IKOVAT=1,NKOVAT
c      WSTIFT(IKOVAT)=0.0D0
c     END DO
      RETURN
C
C**** ZERO STIFFNESS ( & COUPLING MATRICES )
C
   30 CONTINUE       ! set to zero in elm...
c     DO IKOVAT=1,NKOVAT
c      ESTIFT(IKOVAT)=0.0D0
c      IF(KPROBT.EQ.4) WSTIFT(IKOVAT)=0.0 !Only useful for thermal pr.
c     END DO
c     IF(KFLAGT.EQ.1)THEN
c      DO INDEXT=1,NNODET*NEVABT
c       HSTIFT(INDEXT)=0.0D0
c      END DO
c     ENDIF
      RETURN
C
C**** ZERO FORCES
C
   40 DO IEVABT=1,NEVABT
       WORK1T(IEVABT)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL HEAT FLUXES
C
  130 DO INDEXT=IGSMOT(4),IGSMOT(6)-1
       WORK1T(INDEXT)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL TEMP. GRADIENTS
C
  140 DO INDEXT=IGSMOT(6),IGSMOT(8)-1
       WORK1T(INDEXT)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL INTERNAL VARIABLES
C
  150 DO INDEXT=IGSMOT(9),IGSMOT(11)-1
       WORK1T(INDEXT)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL POROSITY CRITERIA
C
  170 DO INDEXT=IGSMOT(23),IGSMOT(25)-1
       WORK1T(INDEXT)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL L*f_pc AND ITS RATE
C
  190 DO INDEXT=IFORCT(53),IFORCT(54)+2*NNODET-1
       WORK1T(INDEXT)=0.0D0
      END DO
      RETURN
C
C**** ZERO PERMEABILITY, COMPRESSIBILITY AND COUPLING MATRICES
C
c 230 DO IKONDT=1,NKONDT
c      PSTIFT(IKONDT)=0.0D0
c      QSTIFT(IKONDT)=0.0D0
c     END DO
c     DO INDEXT=1,NNODET*NEVABT
c      HSTIFT(INDEXT)=0.0D0
c     END DO
c     RETURN
C
      END
