      SUBROUTINE RENUM4(NDEG,LEV,LSD,NADJ,MLW,NODES,NROOT,NLOC)
C***********************************************************************
C
C**** COMPUTE LEVEL STRUCTURE ROOTED AT NROOT
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL*1 BACK
      DIMENSION LEV(*),NDEG(*),NADJ(*)
C
C.... INITIALIZATION
C
      DO I = 1,NODES
        LEV(I) = 0
      END DO
      LEV(NROOT) = 1    
      NDEG(1) = NROOT
      BACK = .FALSE.
      LSD = 1
      KOUNT = 1
      NLOC = 1
      NLOCN = 0
      MLW = 1
C
C.... ASSIGN LEVELS
C
      DO WHILE(KOUNT.LT.NODES)
        DO IL = 1,NLOC
          IF(.NOT.BACK) IP = NDEG(IL)
          IF(BACK) IP = NDEG(NODES+1-IL)
          NAD = NADJ(IP)
          DO WHILE (NAD.GT.0)
            JP = NADJ(NAD)
            IF(LEV(JP).EQ.0) THEN
              NLOCN = NLOCN + 1
              LEV(JP) = LSD + 1
              KOUNT = KOUNT + 1
              IF(BACK) NDEG(NLOCN) = JP
              IF(.NOT.BACK) NDEG(NODES+1-NLOCN) = JP
            END IF
            NAD = NADJ(NAD+1)
          END DO
        END DO
        NLOC = NLOCN
        NLOCN = 0
        IF(NLOC.GT.MLW) MLW = NLOC
        LSD = LSD + 1
        BACK = .NOT.BACK
      END DO
C
      IF(BACK) THEN
        NHALF = NODES/2
        DO I = 1,NHALF
          NN = NDEG(I)
          NDEG(I) = NDEG(NODES+1-I)
          NDEG(NODES+1-I) = NN
        END DO
      END IF
      RETURN
      END
