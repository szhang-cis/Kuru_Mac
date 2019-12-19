      SUBROUTINE DASTOR(DBASE,DAMAG,TAUAM,TAUAP,INDEX)
C***********************************************************************
C
C**** THIS SUBROUTINE STORES INTERNAL VARIABLES & THEIR CONJUGATES
C
C     Note: for IDAMG=21, tau- (TAUAM) & tau+ (TAUAP) must be stored
C           even in absence of damage effects
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'prob_om.f'
C
      DIMENSION DBASE(*)
C
      IF(IDAMG.EQ.0) RETURN
C
      IF(INDEX.EQ.1.OR.INDEX.EQ.2) THEN
       IF(LARGE.NE.0) RETURN      ! not implemented yet for large strains
       IF(IDAMG.NE.21) RETURN
      ENDIF
C
      IF(IDAMG.LE.20) THEN   ! damage models governed by pl. or viscopl.
c      IF(IDAMG.EQ.1) DBASE(1)=DAMAG
c      IF(IDAMG.EQ.4) DBASE(1)=DAMAG
c      IF(IDAMG.EQ.5) DBASE(1)=DAMAG
c      IF(IDAMG.EQ.6) DBASE(1)=DAMAG
       DBASE(1)=DAMAG
      ELSE                   ! damage models governed by a damage crit.
       IF(IDAMG.EQ.21) THEN  ! concrete damage model
        DBASE(1)=TAUAM
        DBASE(2)=TAUAP
       ENDIF
      ENDIF
C
      RETURN
      END
