      SUBROUTINE RUNENDS(MESSAGES)
C***********************************************************************
C
C**** THIS ROUTINE STOPS THE RUN AND MAKES THE SUMMARY OF THE STEPS
C     PERFOMED UP TO NOW
C
C***********************************************************************
      INCLUDE 'prob_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      CHARACTER*35 MESSAGES
C
      IF(IEVFI.EQ.1) WRITE(LUPRIS,900) MESSAGES
      WRITE(LURESS,900) MESSAGES
C
C**** WRITE CPU TIME TABLE
C
      CALL SUMMARS
C
C**** CLOSE WORKING FILES
C
c     CLOSE(LUDTSS)                       ! Input file
      IF(IEVFI.EQ.1) CLOSE(LUPRIS)        ! Print-out file
      CLOSE(LURESS)                       ! Output file
      IF(KPOSTS.EQ.1) CLOSE(LUPOSS)       ! Post-process file
      DO ICURVS=1,NCOLDS                  ! Plotting files
        NFILES=LUCU1S-1+ICURVS
        CLOSE(NFILES)
      ENDDO
C
      STOP 'VULCAN: END OF RUN'
  900 FORMAT(1H1,///,5X,A35)
      END
