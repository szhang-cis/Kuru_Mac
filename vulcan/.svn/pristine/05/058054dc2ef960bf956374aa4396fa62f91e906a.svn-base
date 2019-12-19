      SUBROUTINE IDEPRHT(PROPST,INDEP)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THERMAL BOUNDARY ELEMENTS
C
C***********************************************************************
C
C     NRADH1
C     VRADH1(1)=Radiation-Convection Coefficient values at
C     VRADH1(2)=temperature of body 1
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION  PROPST(*)
C
C**** TEMPERATURE OF BODY 1-DEPENDENCY
C
      I1=INDEP+1                                       ! IFILL=PROPST(1)
      NRADH1=INT(PROPST(I1))
C
      DO IRADH1=1,NRADH1
       IA=I1-1+2*IRADH1
       IB=IA+1
       VRADH1(IRADH1,1)=PROPST(IA)
       VRADH1(IRADH1,2)=PROPST(IB)
      ENDDO
C
      I2=I1+2*NRADH1+1
      NALTF=INT(PROPST(I2))
      VRADH1(NRADH1+1,1)=PROPST(I2)
      IF(NALTF.EQ.1) THEN
       VRADH1(NRADH1+2,1)=PROPST(I2+1)
       VRADH1(NRADH1+2,2)=PROPST(I2+2)
      ENDIF
C
      RETURN
      END
