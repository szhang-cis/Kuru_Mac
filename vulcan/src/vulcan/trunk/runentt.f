      SUBROUTINE RUNENTT(MESSAGET)
C***********************************************************************
C
C**** THIS ROUTINE MAKES THE SUMMARY OF THE STEPS
C     PERFOMED UP TO NOW (DOES NOT STOP THE THERMAL PROBLEM)
C
C***********************************************************************
      INCLUDE 'prob_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      CHARACTER*35 MESSAGET
C
      WRITE(LUPRIT,900) MESSAGET
      WRITE(LUREST,900) MESSAGET
C
C**** WRITE CPU TIME TABLE
C
      CALL SUMMART
C
C**** CLOSE WORKING FILES
C
      CLOSE(LUDTST)                       ! Input file
      CLOSE(LUPRIT)                       ! Print-out file
      CLOSE(LUREST)                       ! Output file
C
      IF(KPOSTT.EQ.1) CLOSE(LUPOST)       ! Post-process file
      DO ICURVT=1,NCOLDT                  ! Plotting files
        NFILET=LUCU1T-1+ICURVT
        CLOSE(NFILET)
      ENDDO
C
  900 FORMAT(1H1,///,5X,A35)
      END