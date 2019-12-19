      SUBROUTINE IDENPRS(PROPSS)
C***********************************************************************
C
C**** THIS ROUTINE CHOOSE THE PROPER SUBROUTINE TO ORDER THE PROPERTIES
C     FOR THE DIFFERENT MODELS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION PROPSS(*)
C
C**** LOOKS FOR MODEL DEFINITION
C
      ISOTRS=INT(PROPSS(NPROPM+1))   ! 0=isotropic; 1=orthotropic
      IFILLS=INT(PROPSS(NPROPM+2))   ! 0=standard; 1=filling material
C
      IF(ISOTRS.EQ.0) THEN
       NDIMETOS=1
       IF(IFILLS.EQ.0) THEN
        CALL IDEPROSS(PROPSS)
       ENDIF
       IF(IFILLS.EQ.1) THEN
        call runends('ideprofs not implemented - idenprs')
c       CALL IDEPROFS(PROPSS)
       ENDIF
      ENDIF
C
      IF(ISOTRS.EQ.1) THEN
       NDIMETOS=NDIMES
       IF(IFILLS.EQ.0) THEN
        call runends('ideproso not implemented - idenprs')
c       CALL IDEPROSO(PROPSS)
       ENDIF
       IF(IFILLS.EQ.1) THEN
        call runends('ideprfos not implemented - idenprs')
c       CALL IDEPRFOS(PROPSS)
       ENDIF
      ENDIF
C
      IF(ISOTRS.EQ.2) THEN
       NDIMETOS=NDIMES*NDIMES
       CALL RUNENDS('ERROR: FULLY ANISOTROPIC MATERIAL NOT IMPLEM.')
      ENDIF
C
      RETURN
      END
