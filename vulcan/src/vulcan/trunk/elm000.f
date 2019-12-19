      SUBROUTINE ELM000(WORK1,NWORI,
     .                  CSTIF,ESTIF,WSTIF,HSTIF,PSTIF,QSTIF,    !ELMAT
     .                                                   ITASK)
C***********************************************************************
C
C**** THIS ROUTINE ZEROES RELEVANT ARRAYS FOR DIFFERENT ITASK VALUES 
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
      INCLUDE 'auxl_om.f'
C
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*), HSTIF(*), PSTIF(*),
     .          QSTIF(*)
      DIMENSION WORK1(*)
      DIMENSION NWORI(*)
C
C**** DIRECT TO APPROPRIATE PROCESSOR
C
      GO TO (  1, 20, 30, 40,  1, 40, 40, 40,  1,  1,
     .         1,  1,130,140,150,  1,  1,  1,  1), ITASK
    1 RETURN ! Nothing
C
C**** ZERO MASS MATRIX
C
   20 CONTINUE        ! set to zero in elm030
c     DO IKOVA=1,NKOVA
c      WSTIF(IKOVA)=0.0
c     END DO
      RETURN
C
C**** ZERO STIFFNESS ( & COUPLING MATRICES )
C
   30 CONTINUE        ! set to zero in elm030
c     DO IKOVA=1,NKOVA
c      ESTIF(IKOVA)=0.0
c     END DO
c     IF(KFLAG.EQ.1)THEN
c      DO INDEX=1,NNODE*NEVAB
c       HSTIF(INDEX)=0.0
c      END DO
c     ENDIF
      RETURN
C
C**** ZERO FORCES
C
   40 DO IEVAB=1,NEVAB
       WORK1(IEVAB)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL STRESSES
C
  130 DO INDEX=IGSMO(4),IGSMO(6)-1
       WORK1(INDEX)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL STRAINS
C
  140 DO INDEX=IGSMO(6),IGSMO(8)-1
       WORK1(INDEX)=0.0D0
      END DO
      RETURN
C
C**** ZERO ELEMENTAL CONTRIBUTION FOR NODAL INTERNAL VARIABLES
C
  150 DO INDEX=IGSMO(9),IGSMO(11)-1
       WORK1(INDEX)=0.0D0
      END DO
      RETURN
C
      END
