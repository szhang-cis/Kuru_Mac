      SUBROUTINE INIMDT(DISTO)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES DISPLACEMENT VECTOR AT FICTITIOUS TIMES
C
C     Only needed for: CENTRAL DIFFERENCE METHOD (KINTE=3)
C                      HOUBOLT METHOD            (KINTE=4)
C
C     Note: the acceleration at time 0 should be computed (is zero
C           only if U_0 is zero and no body forces are present)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION DISTO(NTOTV,*)
C
      IF(KDYNA.EQ.0) RETURN
      IF(KINTE.EQ.1.OR.KINTE.EQ.2.OR.KINTE.EQ.5.OR.KINTE.EQ.6) RETURN
C
      IF(KINTE.EQ.3) THEN                        ! computes (-delta t)^U
       DO ITOTV=1,NTOTV
        DISTO(ITOTV,4)=DISTO(ITOTV,1)-DTIME*DISTO(ITOTV,2)+
     .                                  DTIME*DTIME/2.0D0*DISTO(ITOTV,3)
       ENDDO
      ENDIF
C
      RETURN
      END
