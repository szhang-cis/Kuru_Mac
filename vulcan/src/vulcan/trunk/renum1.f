      SUBROUTINE RENUM1(LNODS,NODAD,NNODE,NPOIN,NELEM,NPOSI)
C***********************************************************************
C
C**** COMPUTES THE CONNECTIONS BETWEEN DEGREES OF FREEDOM (STORED AS 
C     A LINKED LIST)
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION LNODS(nnode,*),NODAD(*)
      DO ILIB = 1,NPOIN
         NODAD(ILIB) = 0
      END DO
      NFREE = NPOIN+1

      DO IELEM = 1,NELEM
        DO INODE = 1,NNODE
          ILIB = LNODS(INODE,IELEM)
          IF(ILIB.NE.0) THEN
            DO JNODE = 1,NNODE
              IF(JNODE.NE.INODE) THEN
                JLIB = LNODS(JNODE,IELEM)
                IF(JLIB.NE.0) THEN
                   NADR = NODAD(ILIB)
                   NADRO = ILIB
                   DO WHILE (NADR.GT.0)
                     KLIB = NODAD(NADR)
                     IF(KLIB.EQ.JLIB) GO TO 10
                     NADRO = NADR + 1
                     NADR = NODAD(NADRO)
                   END DO
                   NODAD(NADRO) = NFREE
                   NODAD(NFREE) = JLIB
                   NODAD(NFREE+1) = 0
                   NFREE = NFREE + 2
10                 CONTINUE
                 END IF
              END IF
            END DO
          END IF
        END DO
      END DO
      NPOSI = NFREE
      RETURN
      END
