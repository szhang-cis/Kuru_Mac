      SUBROUTINE ASAL00(ESTIF,WSTIF,ESTII,WSTII,NKOVA,KDYNA,NMEMO6)
C***********************************************************************
C
C**** THIS ROUTINE MANIPULATES CONTRIBUTIONS OF JACOBIAN MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ESTIF(NKOVA), WSTIF(NKOVA),
     .          ESTII(NKOVA), WSTII(NKOVA)
C
      IF(NMEMO6.EQ.0) THEN
       DO IKOVA=1,NKOVA
        ESTIF(IKOVA)=ESTII(IKOVA)
        IF(KDYNA.EQ.1) WSTIF(IKOVA)=WSTII(IKOVA)
       ENDDO
      ENDIF
C
      RETURN
      END
