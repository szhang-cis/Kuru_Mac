      SUBROUTINE MASS01(DVOLU,PROPS,SHAPE,WSTIF)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE MASS MATRIX ( ELEMENT NO. 1 )
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DVOLU(*),  PROPS(*),  SHAPE(NNODL,*)
      DIMENSION WSTIF(*)
C
C**** LOOP ON INTEGRATION POINTS
C
      DO 100 IGAUS=1,NGAUL
C
C**** ELEMENTAL MASS MATRIX
C
      CALL WMATRI(DVOLU(IGAUS),   NDOFN,NEVAB,NNODL, PROPS,
     .            SHAPE(1,IGAUS), WSTIF,KSYMM)
C
  100 CONTINUE
C
C**** COMPLETE MATRIX FOR RECTANGULAR CASE
C
      IF(KSYMM.EQ.0) THEN
       DO IEVAB=1,NEVAB
        DO JEVAB=IEVAB,NEVAB
         KLOCS=(IEVAB-1)*NEVAB+JEVAB
         KLOCI=(JEVAB-1)*NEVAB+IEVAB
         WSTIF(KLOCI)=WSTIF(KLOCS)
        END DO
       END DO
      END IF
C
      RETURN
      END
