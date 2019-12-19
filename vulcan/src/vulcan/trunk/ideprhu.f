      SUBROUTINE IDEPRHU(PROPST,INDEP)
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
C**** TEMPERATURE OF BODY 2-DEPENDENCY
C
      I2=I1+2*NRADH1+1
      NRADH2=INT(PROPST(I2))
C
      DO IRADH2=1,NRADH2
       IA=I2-1+2*IRADH2
       IB=IA+1
       VRADH2(IRADH2,1)=PROPST(IA)
       VRADH2(IRADH2,2)=PROPST(IB)
      ENDDO
C
C**** GAP-DEPENDENCY
C
      I3=I2+2*NRADH2+1
      NRADU=INT(PROPST(I3))
C
      DO IRADU=1,NRADU
       IA=I3-1+2*IRADU
       IB=IA+1
       VRADU(IRADU,1)=PROPST(IA)
       VRADU(IRADU,2)=PROPST(IB)
      ENDDO
C
C**** PRESSURE-DEPENDENCY
C
      I4=I3+2*NRADU+1
      NRADP=INT(PROPST(I4))
C
      DO IRADP=1,NRADP
       IA=I4-1+2*IRADP
       IB=IA+1
       VRADP(IRADP,1)=PROPST(IA)
       VRADP(IRADP,2)=PROPST(IB)
      ENDDO
C
C**** NORMAL STIFFNESS COMING FROM THE MECHANICAL PROBLEM (see checkp.f)
C
      I5=I4+2*NRADP+1
      RIGINT=PROPST(I5)
      I6=I5+1
      IR432=INT(PROPST(I6))
C
C**** FRICTIONAL HEATING
C
      I10=I6+1
      IGAFRT=INT(PROPST(I10))
C
      IF(IGAFRT.EQ.1) THEN
       I10=I10+1
       AGAFPT=PROPST(I10)
      ENDIF
C
C**** BOUNDARY HEATING
C
      I10=I10+1
      IGABO3=INT(PROPST(I10))
C
      IF(IGABO3.EQ.1) THEN
       I10=I10+1
       TENVI3=PROPST(I10)
C
       I10=I10+1
       NRADH3=INT(PROPST(I10))
C
       DO IRADH3=1,NRADH3
        IA=I10-1+2*IRADH3
        IB=IA+1
        VRADH3(IRADH3,1)=PROPST(IA)
        VRADH3(IRADH3,2)=PROPST(IB)
       ENDDO
C
       I10=I10+2*NRADH3+1
       IGABO4=INT(PROPST(I10))
C
       IF(IGABO4.EQ.1) THEN
        I10=I10+1
        AGABO4=PROPST(I10)
       ENDIF
C
      ENDIF
C
      RETURN
      END
