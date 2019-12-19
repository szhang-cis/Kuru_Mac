      SUBROUTINE RENUMN(KRENU,LNODS,LPNTN,NELEM,NNODE,NPOIN,WORK1)
C***********************************************************************
C                 
C**** THIS ROUTINE SETS UP ARRAY LPNTN IF NODE RENUMBERING IS DESIRED
C                 
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
C
      DIMENSION LNODS(NNODE,*), LPNTN(*), WORK1(*)
C
      CALL CPUTIM(TIME1)
C
C**** NODE RENUMBERING FOR PROFILE MINIMIZATION
C
      IF(KRENU.EQ.1) THEN
       CALL RENUM0(LNODS,LPNTN,NELEM,NNODE,NPOIN,WORK1)
      ENDIF
C
      CALL CPUTIM(TIME2)
      CPURN=CPURN+(TIME2-TIME1)
C
      RETURN
      END
