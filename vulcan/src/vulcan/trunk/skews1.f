      SUBROUTINE SKEWS1(JSKE1,ORIGI,ZLOCA,XLOCA,COFRI)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ROTATION MATRIX USING OPTION 1
C     
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION COFRI(NSKEW,NDIME,*)
      DIMENSION ORIGI(3), ZLOCA(3), XLOCA(3)
      DIMENSION YLOCA(3)
C
      XNORM=0.0D+00
      ZNORM=0.0D+00
C
C**** 1-D
C
      IF(NDIME.EQ.1) THEN
       XLOCA(1)=XLOCA(1)-ORIGI(1)
       XNORM=XNORM+XLOCA(1)*XLOCA(1)
       XNORM=DSQRT(XNORM)
       XLOCA(1)=XLOCA(1)/XNORM
C
       COFRI(JSKE1,1,1)=XLOCA(1)
      ENDIF
C
C**** 2-D
C
      IF(NDIME.EQ.2) THEN
       DO IDIME=1,NDIME
        XLOCA(IDIME)=XLOCA(IDIME)-ORIGI(IDIME)
        XNORM=XNORM+XLOCA(IDIME)*XLOCA(IDIME)
       ENDDO
       XNORM=DSQRT(XNORM)
       DO IDIME=1,NDIME
        XLOCA(IDIME)=XLOCA(IDIME)/XNORM
       ENDDO
       YLOCA(1)=-XLOCA(2)
       YLOCA(2)= XLOCA(1)
C
       DO IDIME=1,NDIME
        COFRI(JSKE1,1,IDIME)=XLOCA(IDIME)
        COFRI(JSKE1,2,IDIME)=YLOCA(IDIME)
       ENDDO
      ENDIF
C
C**** 3-D
C
      IF(NDIME.EQ.3) THEN
       DO IDIME=1,NDIME
        XLOCA(IDIME)=XLOCA(IDIME)-ORIGI(IDIME)
        ZLOCA(IDIME)=ZLOCA(IDIME)-ORIGI(IDIME)
        XNORM=XNORM+XLOCA(IDIME)*XLOCA(IDIME)
        ZNORM=ZNORM+ZLOCA(IDIME)*ZLOCA(IDIME)
       ENDDO
       XNORM=DSQRT(XNORM)
       ZNORM=DSQRT(ZNORM)
       DO IDIME=1,NDIME
        XLOCA(IDIME)=XLOCA(IDIME)/XNORM
        ZLOCA(IDIME)=ZLOCA(IDIME)/ZNORM
       ENDDO
       YLOCA(1)=ZLOCA(2)*XLOCA(3)-ZLOCA(3)*XLOCA(2)
       YLOCA(2)=ZLOCA(3)*XLOCA(1)-ZLOCA(1)*XLOCA(3)
       YLOCA(3)=ZLOCA(1)*XLOCA(2)-ZLOCA(2)*XLOCA(1)
C
       DO IDIME=1,NDIME
        COFRI(JSKE1,1,IDIME)=XLOCA(IDIME)
        COFRI(JSKE1,2,IDIME)=YLOCA(IDIME)
        COFRI(JSKE1,3,IDIME)=ZLOCA(IDIME)
       ENDDO
      ENDIF
C
      RETURN
      END