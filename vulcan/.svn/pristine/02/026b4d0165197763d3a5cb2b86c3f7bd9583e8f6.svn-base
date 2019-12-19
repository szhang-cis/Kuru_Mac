      SUBROUTINE ASELMT(MATNO,PROEL,PROPS,LNODS,INFRI,COFRI,
     .                  CSTIF,ESTIF,WSTIF,PSTIF,QSTIF,HSTIF,
     .                  WORK1)
C***********************************************************************
C
C**** THIS ROUTINE ASSEMBLES ELEMENTAL "EFFECTIVE" MATRIX   
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION MATNO(*), PROEL(*), PROPS(*)
      DIMENSION CSTIF(*), ESTIF(*), WSTIF(*),
     .          PSTIF(*), QSTIF(*), HSTIF(*)
C
      DIMENSION LNODS(NNODE,NELEM),   INFRI(NPOIN),
     .          COFRI(NSKEW,NDIME,*), WORK1(*)
C
      CALL CPUTIM(TIME1)
C
      IF(KPROB.EQ.4) THEN
       CALL ASELM3(MATNO,PROEL,PROPS,CSTIF,ESTIF,WSTIF)
      ELSE IF(KPORE.NE.2) THEN
       CALL ASELM1(MATNO,PROEL,PROPS,CSTIF,ESTIF,WSTIF,LNODS,
     .             WORK1(ILOSY(1)),WORK1(ILOSY(2)),INFRI,COFRI)
      ELSE
       CALL ASELM2(MATNO,PROEL,PROPS,CSTIF,ESTIF,WSTIF,PSTIF,QSTIF,
     .             HSTIF)
      ENDIF
C
      CALL CPUTIM(TIME2)
      CPUAS=CPUAS+(TIME2-TIME1)
C
      RETURN
      END
