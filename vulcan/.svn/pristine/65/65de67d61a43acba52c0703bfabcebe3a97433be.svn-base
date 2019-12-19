      SUBROUTINE RENUM3(NADJ,NACT,NODA,NEWNN,NODES,I)
C***********************************************************************
C
C**** RESEQUENCE NODES FOR MINIMUM PROFILE
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NADJ(*),NEWNN(*),NACT(*),NODA(*)
 
      LARGE = 5**5
C
C.... KING'S SCHEME
C
      DO J = 1,NODES
        NEWNN(J) = 0
        NODA(J) = 0
      END DO
      NEWNN(I) = 1
      NAC = 0
C
C.... NEGATE ALL NDEG ENTRIES FOR NODES WHICH ARE
C.... ADJACENT TO STARTING NODE I
C
      MAXFRT = 0
      NAD = NADJ(I)
      DO WHILE (NAD.GT.0)
        MAXFRT = MAXFRT + 1
        NPJ = NADJ(NAD)
        IF(NODA(NPJ).EQ.0) THEN
          NAC = NAC + 1
          NODA(NPJ) = NAC
          NACT(NAC) = NPJ
        END IF
        NAD = NADJ(NAD+1)
      END DO
      NODA(I) = LARGE
C
C.... LOOP OVER NODES TO BE RENUMBERED
C
      DO K = 2,NODES
        MINNEW = LARGE
        LMIN = LARGE
C
C.... LOOP OVER ACTIVE NODES
C.... SKIP TO NEXT NODE IF OLD NODE IS ALREADY RENUMBERED
C
        DO IACT = 1,NAC
          J = NACT(IACT)
          IF(NEWNN(J).LE.0) THEN
            NEW = -1
            MIN = LARGE
C
C.... COMPUTE THE INCREMENT IN ACTIVE NODES FOR EACH NODE J
C.... COMPUTE WHEN THIS NODE WAS FIRST ACTIVE BY CHECKING FOR RENUMBERED
C.... NEIGHBOURS WITH LOWEST NUMBERS
C
            NAD = NADJ(J)
            DO WHILE (NAD.GT.0)
              N = NADJ(NAD)
              IF(NODA(N).EQ.0) NEW = NEW + 1
              IF(NEWNN(N).NE.0) THEN
                IF(NEWNN(N).LT.MIN)MIN = NEWNN(N)
              END IF
              NAD = NADJ(NAD+1)
            END DO
C
C.... SELECT NODE WITH SMALLEST INCREMENT IN ACTIVE NODES
C.... IN THE CASE OF A TIE, SELECT NODE WHICH HAS BEEN LONGEST ACTIVE
C
            IF(NEW.LE.MINNEW) THEN
              IF((NEW.NE.MINNEW).OR.(MIN.LT.LMIN)) THEN
                MINNEW = NEW
                LMIN = MIN
                NEXT = J
              END IF
            END IF
          END IF
        END DO
C
C.... RENUMBER NODE AND COMPUTE NUMBER OF ACTIVE NODES
C
        NEWNN(NEXT) = K
        NIF = NIF+MINNEW
        IF(NIF.GT.MAXFRT) MAXFRT = NIF
C
C.... SET NODES WHICH ARE ADJACENT TO THE NODE JUST RENUMBERED
C.... AS ACTIVES NODES, DEACTIVATE NEXT
C
        NPOS = NODA(NEXT)
        ILAST = NACT(NAC)
        NODA(ILAST) = NPOS
        NACT(NPOS) = ILAST
        NAC = NAC - 1

        IF(MINNEW .NE. -1) THEN
          NAD = IABS(NADJ(NEXT))
          DO WHILE (NAD.GT.0)
            N = NADJ(NAD)
            IF(NODA(N).EQ.0) THEN
              NAC = NAC + 1
              NODA(N) = NAC
              NACT(NAC) = N
            END IF
            NAD = NADJ(NAD+1)
          END DO
        END IF
      END DO
      RETURN
      END
