      SUBROUTINE IDENPHT(PROPST,INDEP)
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
C     INDEP=1    >>    elm104
C     INDEP=2    >>    elm101 (TENVI)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION  PROPST(*)
C
C**** LOOKS FOR MODEL DEFINITION
C
      IFILL=INT(PROPST(1))   ! 0=standard; 1=filling material
C
      IF(IFILL.EQ.0) THEN
       IF(INDEP.EQ.1) THEN                       ! gap element
        CALL IDEPRHU(PROPST,    1)
       ELSE                                      ! boundary element
        CALL IDEPRHT(PROPST,    2)
       ENDIF
      ENDIF
      IF(IFILL.EQ.1) THEN
       IF(INDEP.EQ.1) THEN                       ! gap element
        CALL IDEPRFU(PROPST,    1)
       ELSE                                      ! boundary element
        CALL IDEPRFT(PROPST,    2)
       ENDIF
      ENDIF
C
      RETURN
      END
