      SUBROUTINE RENUMNS(KRENUS,LNODSS,LPNTNS,NELEMS,NNODES,NPOINS,
     .                   WORK1S)
C***********************************************************************
C                 
C**** THIS ROUTINE SETS UP ARRAY LPNTN IF NODE RENUMBERING IS DESIRED
C                 
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_oms.f'
C
      DIMENSION LNODSS(NNODES,*), LPNTNS(*), WORK1S(*)
C
      CALL CPUTIMS(TIME1S)
C
C**** NODE RENUMBERING FOR PROFILE MINIMIZATION
C
      IF(KRENUS.EQ.1) THEN
       call runends('ERROR: renumbering not implemented')
c      CALL RENUM0S(LNODSS,LPNTNS,NELEMS,NNODES,NPOINS,WORK1S)
      ENDIF
C
      CALL CPUTIMS(TIME2S)
      CPURNS=CPURNS+(TIME2S-TIME1S)
C
      RETURN
      END
