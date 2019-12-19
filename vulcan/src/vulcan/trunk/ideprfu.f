      SUBROUTINE IDEPRFU(PROPST,INDEP)
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
C     NRADU
C     VRADU(1)=Radiation-Convection Coefficient values at
C     VRADU(2)=normal gap
C
C     NRADP
C     VRADP(1)=Radiation-Convection Coefficient values at
C     VRADP(2)=normal pressure
C
C     RIGINT=normal stiffness
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
C     NRADUFI
C     VRADUFI(1)=Radiation-Convection Coefficient values at
C     VRADUFI(2)=normal gap
C
C     NRADPFI
C     VRADPFI(1)=Radiation-Convection Coefficient values at
C     VRADPFI(2)=normal pressure
C
C
C     INDEP=1    >>    elm104t
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*)
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
C**** TEMPERATURE OF BODY 2-DEPENDENCY
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
C**** GAP-DEPENDENCY
C
      I3=I2F+2*NRADH2FI+1
      NRADU=INT(PROPST(I3))
C
      DO IRADU=1,NRADU
       IA=I3-1+2*IRADU
       IB=IA+1
       VRADU(IRADU,1)=PROPST(IA)
       VRADU(IRADU,2)=PROPST(IB)
      ENDDO
C
      I3F=I3+2*NRADU+1
      NRADUFI=INT(PROPST(I3F))
C
      DO IRADU=1,NRADUFI
       IA=I3F-1+2*IRADU
       IB=IA+1
       VRADUFI(IRADU,1)=PROPST(IA)
       VRADUFI(IRADU,2)=PROPST(IB)
      ENDDO
C
C**** PRESSURE-DEPENDENCY
C
      I4=I3F+2*NRADUFI+1
      NRADP=INT(PROPST(I4))
C
      DO IRADP=1,NRADP
       IA=I4-1+2*IRADP
       IB=IA+1
       VRADP(IRADP,1)=PROPST(IA)
       VRADP(IRADP,2)=PROPST(IB)
      ENDDO
C
      I4F=I4+2*NRADP+1
      NRADPFI=INT(PROPST(I4F))
C
      DO IRADP=1,NRADPFI
       IA=I4F-1+2*IRADP
       IB=IA+1
       VRADPFI(IRADP,1)=PROPST(IA)
       VRADPFI(IRADP,2)=PROPST(IB)
      ENDDO
C
C**** NORMAL STIFFNESS COMING FROM THE MECHANICAL PROBLEM (see checkp.f)
C
      I5=I4F+2*NRADPFI+1
      RIGINT=PROPST(I5)
      I6=I5+1
      IR432=INT(PROPST(I6))
C
      RETURN
      END
