      SUBROUTINE STIF03(PROPS,ESTIF,TENOD)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS MATRIX FOR THE LINK 
C     (FREE SURFACE ELEMENT NO. 3)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION PROPS(NPROP), ESTIF(NKOVA)
      DIMENSION TENOD(NNODL)
C
C**** PROPS(1+IDOFN)=PENALTIES PARAMETERS IN X-Y-Z DIRECTIONS
C     PROPS(5)=LIQUIDUS TEMPERATURE
C
      TEMPL=PROPS(5)
C
      IF(TENOD(1).GT.TEMPL.AND.TENOD(2).GT.TEMPL) THEN
C
C**** EVALUATE ELEMENT CONTRIBUTION
C
       DO IDOFN=1,NDOFN
        I=IDOFN
        J=IDOFN
        I1=(2*NEVAB-I)*(I-1)/2+J
        IF(KSYMM.EQ.0) I1=(I-1)*NEVAB+J  ! unsymmetric
        I=IDOFN
        J=IDOFN+NDOFN
        I2=(2*NEVAB-I)*(I-1)/2+J
        IF(KSYMM.EQ.0) I2=(I-1)*NEVAB+J  ! unsymmetric
        I=IDOFN+NDOFN
        J=IDOFN+NDOFN
        I3=(2*NEVAB-I)*(I-1)/2+J
        IF(KSYMM.EQ.0) I3=(I-1)*NEVAB+J  ! unsymmetric
C
        ESTIF(I1)= PROPS(1+IDOFN)
        ESTIF(I2)=-PROPS(1+IDOFN)
        ESTIF(I3)= PROPS(1+IDOFN)
       END DO ! IDOFN=1,NDOFN
C
C**** LOAD ESTIF IN A SQUARE FORM FOR UNSYMMETRIC SOLVER
C
       IF(KSYMM.EQ.0) THEN
        DO IEVAB=1,NEVAB
         DO JEVAB=IEVAB,NEVAB
          KLOCS=(IEVAB-1)*NEVAB+JEVAB
          KLOCI=(JEVAB-1)*NEVAB+IEVAB
          ESTIF(KLOCI)=ESTIF(KLOCS)
         END DO
        END DO
       END IF
C
      ENDIF
C
      RETURN
      END
