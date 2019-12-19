      SUBROUTINE RUNEND(MESSAGE)
C***********************************************************************
C
C**** THIS ROUTINE STOPS THE RUN AND MAKES THE SUMMARY OF THE STEPS
C     PERFOMED UP TO NOW
C
C***********************************************************************
      INCLUDE 'prob_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      CHARACTER(LEN=*) ::  MESSAGE
C
      WRITE(LUPRI,900) MESSAGE
      WRITE(LURES,900) MESSAGE
C
C**** WRITE CPU TIME TABLE
C
      CALL SUMMAR
C
C**** CLOSE WORKING FILES
C
      CLOSE(LUDTS)                      ! Input file
      CLOSE(LUPRI)                      ! Print-out file
      CLOSE(LURES)                      ! Output file
      IF(KPOST.EQ.1) CLOSE(LUPOS)       ! Post-process file
      DO ICURV=1,NCOLD                  ! Plotting files
        NFILE=LUCU1-1+ICURV
        CLOSE(NFILE)
      ENDDO
C
      STOP 'VULCAN: END OF RUN'
  900 FORMAT(1H1,///,5X,A)
      END
