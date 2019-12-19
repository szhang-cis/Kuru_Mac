      SUBROUTINE RENUMNT(KRENUT,LNODST,LPNTNT,NELEMT,NNODET,NPOINT,
     .                   WORK1T)
C***********************************************************************
C                 
C**** THIS ROUTINE SETS UP ARRAY LPNTN IF NODE RENUMBERING IS DESIRED
C                 
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_omt.f'
C
      DIMENSION LNODST(NNODET,*), LPNTNT(*), WORK1T(*)
C
      CALL CPUTIMT(TIME1T)
C
C**** NODE RENUMBERING FOR PROFILE MINIMIZATION
C
      IF(KRENUT.EQ.1) THEN
       CALL RENUM0T(LNODST,LPNTNT,NELEMT,NNODET,NPOINT,WORK1T)
      ENDIF
C
      CALL CPUTIMT(TIME2T)
      CPURNT=CPURNT+(TIME2T-TIME1T)
C
      RETURN
      END
