
      SUBROUTINE BDBCO2(BMATX,DMATX,EMATX,NCOLU,NROWS)
C**********************************************************************
C                             T 
C****THIS ROUTINE PERFORMS E=B D B OF TWO GENERAL MATRICES
C    WHEN DMATX IS DEFINED AS AN ARRAY AND THE RESULTING 
C    EMATX IS GIVEN AS AN ARRAY TOO.
C
C                 BMATX = NROWS * NCOLU
C                 DMATX = NROWS * NROWS
C                 EMATX = NCOLU * NCOLU
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BMATX(NROWS,*), DMATX(*), EMATX(*)
C
      NKOVA=NCOLU*(NCOLU+1)/2
      DO IKOVA=1,NKOVA
        EMATX(IKOVA)=0.0
      ENDDO
C
      KOUN1=0
      DO ICOLU=1,NCOLU
        DO JCOLU=ICOLU,NCOLU
          KOUN1=KOUN1+1
          DO KROWS=1,NROWS
C
            KOUNT=KROWS
            DO LROW1=1,KROWS
              EMATX(KOUN1)=EMATX(KOUN1)+
     .            BMATX(KROWS,ICOLU)*DMATX(KOUNT)*BMATX(LROW1,JCOLU)
              KOUNT=KOUNT+(NROWS-LROW1)
            ENDDO
C  
            KOUNT=KOUNT-(NROWS-(LROW1-1))
            DO LROW1=LROW1,NROWS
              KOUNT=KOUNT+1
              EMATX(KOUN1)=EMATX(KOUN1)+
     .             BMATX(KROWS,ICOLU)*DMATX(KOUNT)*BMATX(LROW1,JCOLU)
            ENDDO
C
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END
