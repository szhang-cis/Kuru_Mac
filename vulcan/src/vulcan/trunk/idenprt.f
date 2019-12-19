      SUBROUTINE IDENPRT(PROPST)
C***********************************************************************
C
C**** THIS ROUTINE CHOOSE THE PROPER SUBROUTINE TO ORDER THE PROPERTIES
C     FOR THE DIFFERENT MODELS
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
C**** LOOKS FOR MODEL DEFINITION
C
      ISOTRT=INT(PROPST(1))  ! 0=isotropic; 1=orthotropic
      IFILL=INT(PROPST(2))   ! 0=standard; 1=filling material
C
      IF(ISOTRT.EQ.0) THEN
       NDIMETO=1
       IF(IFILL.EQ.0) THEN
        CALL IDEPROT(PROPST)
       ENDIF
       IF(IFILL.EQ.1) THEN
        CALL IDEPROF(PROPST)
       ENDIF
      ENDIF
C
      IF(ISOTRT.EQ.1) THEN
       NDIMETO=NDIMET
       IF(IFILL.EQ.0) THEN
        CALL IDEPROTO(PROPST)
       ENDIF
       IF(IFILL.EQ.1) THEN
        CALL IDEPROFO(PROPST)
       ENDIF
      ENDIF
C
      IF(ISOTRT.EQ.2) THEN
       NDIMETO=NDIMET*NDIMET
       CALL RUNENDT('ERROR: FULLY ANISOTROPIC MATERIAL NOT IMPLEM.')
      ENDIF
C
      RETURN
      END
