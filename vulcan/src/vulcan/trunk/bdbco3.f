      SUBROUTINE BDBCO3( BMATX,DMATX,ESTIF,CONST,NVAR,NSTRE,NDIME )
C***********************************************************************
C
C**** THIS ROUTINE COMPUTE AND ADD UPPER TRIANGLE OF THE PRODUCT OF
C     MATRIX
C
C     ESTIF(P) + CONST*BMATX(K,I)*DMATX(K,L)*BMATX(L,J) -> ESTIF(P)
C
C         I = 1..NVAR   J = I..NVAR   UPPER TRIANGLE P(I,J)
C         K,L = 1..NSTRE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ESTIF(*), BMATX(NSTRE,NVAR), DMATX(NSTRE,NSTRE)
C
      NQUE=NDIME-NVAR
      DO ISTRE=1,NSTRE
         DO JSTRE=1,NSTRE
            K=0
            TEMP1=DMATX(ISTRE,JSTRE)*CONST
            DO IVAR=1,NVAR
              TEMP2=BMATX(ISTRE,IVAR)*TEMP1
              DO JVAR=IVAR,NVAR
                 K=K+1
                 ESTIF(K) = ESTIF(K) + TEMP2*BMATX(JSTRE,JVAR)
              END DO
              K=K+NQUE
           END DO
        END DO
      END DO
      RETURN
      END
