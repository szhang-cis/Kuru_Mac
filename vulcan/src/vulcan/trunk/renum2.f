      SUBROUTINE RENUM2(NADJ,NSTART,LEV,NAUX,NODES,NS)
C***********************************************************************
C
C**** COMPUTES A SET OF POSIBLES STARTINGS NODES
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BETTER
      DIMENSION NAUX(*),NADJ(*),LEV(*),NSTART(*)
C
C.... BEGIN ITERATION
C.... SELECT INTIAL ROOT NODE ARBITRARILY AND GENERATE ITS LEVEL
C.... STRUCTURE
C
      IROOT = 1
      BETTER = .TRUE.
      DO WHILE (BETTER)

        CALL RENUM4(NSTART,LEV,IDEPTH,NADJ,IWIDTH,NODES,IROOT,LHW)
C
C.... CREATE A LIST OF NODES WHICH ARE AT MAXIMUM DISTANCE FROM ROOT
C.... NODE AND STORE THE ROOT 
C
        NS = IROOT
C
C.... LOOP OVER NODES AT MAXIMUM DISTANCE FROM ROOT NODE
C.... GENERATE LEVEL STRUCTURE FOR EACH NODE
C.... SET SWITCH IF A LEVEL STRUCTURE OF GREATER DEPTH OCCURS
C
        BETTER = .FALSE.

        DO I = 1,LHW
          NNODED = NSTART(I)
          CALL RENUM4(NAUX,LEV,NDEPTH,NADJ,NWIDTH,NODES,NNODED,LR)
          IF(NDEPTH.GE.IDEPTH) THEN
            IF((NDEPTH.NE.IDEPTH).OR.(NWIDTH.LT.IWIDTH)) THEN
              IROOT  = NNODED
              IDEPTH = NDEPTH
              IWIDTH = NWIDTH
              BETTER = .TRUE.
            END IF
          END IF
        END DO
      END DO
      RETURN
      END
