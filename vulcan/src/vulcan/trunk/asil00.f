      SUBROUTINE ASIL00(ESTII,WSTII,NKOVA,KDYNA,
     .                  EMATX,NEVAB,INDES)
C***********************************************************************
C
C**** THIS ROUTINE MANIPULATES CONTRIBUTIONS OF JACOBIAN MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ESTII(*), WSTII(*), EMATX(NEVAB,*)
C
      DO IKOVA=1,NKOVA
       ESTII(IKOVA)=0.0D+00
       IF(KDYNA.EQ.1) WSTII(IKOVA)=0.0D+00
      ENDDO
C
      IF(INDES.EQ.1) THEN
       DO IEVAB=1,NEVAB
        DO JEVAB=1,NEVAB
         EMATX(IEVAB,JEVAB)=0.0D0
        ENDDO
       ENDDO
      ENDIF
C
      RETURN
      END
