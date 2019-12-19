      SUBROUTINE IDEPRFT(PROPST,INDEP)
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
C     NRADH2
C     VRADH2(1)=Radiation-Convection Coefficient values at
C     VRADH2(2)=temperature of body 2
C
C
C     NRADH1FI
C     VRADH1FI(1)=Radiation-Convection Coefficient values at
C     VRADH1FI(2)=temperature of body 1
C
C     NRADH2FI
C     VRADH2FI(1)=Radiation-Convection Coefficient values at
C     VRADH2FI(2)=temperature of body 2
C
C
C     INDEP=1    >>    elm104
C     INDEP=2    >>    elm101 (TENVI)
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
      I1F=I1+2*NRADH1+1
      NRADH1FI=INT(PROPST(I1F))
C
      DO IRADH1=1,NRADH1FI
       IA=I1F-1+2*IRADH1
       IB=IA+1
       VRADH1FI(IRADH1,1)=PROPST(IA)
       VRADH1FI(IRADH1,2)=PROPST(IB)
      ENDDO
C
      IF(INDEP.EQ.2) RETURN
C
C**** TEMPERATURE OF BODY 2-DEPENDENCY (ONLY FOR GAP ELEMENT)
C
      I2=I1F+2*NRADH1FI+1
      NRADH2=INT(PROPST(I2))
C
      DO IRADH2=1,NRADH2
       IA=I2-1+2*IRADH2
       IB=IA+1
       VRADH2(IRADH2,1)=PROPST(IA)
       VRADH2(IRADH2,2)=PROPST(IB)
      ENDDO
C
      I2F=I2+2*NRADH2+1
      NRADH2FI=INT(PROPST(I2F))
C
      DO IRADH2=1,NRADH2FI
       IA=I2F-1+2*IRADH2
       IB=IA+1
       VRADH2FI(IRADH2,1)=PROPST(IA)
       VRADH2FI(IRADH2,2)=PROPST(IB)
      ENDDO
C
      RETURN
      END
