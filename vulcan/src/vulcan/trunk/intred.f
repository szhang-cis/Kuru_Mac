      SUBROUTINE INTRED(TFICT)
C*********************************************************************
C
C**** THIS ROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
C
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION TFICT(NFUNC)
C
      CHARACTER XCOMM*9,YCOMM*123
C
C**** NEW RUN: LOOK FOR 'INTE'RVAL CARD OR FOR END OF ANALYSIS ('STOP')
C
      IF(IREST.EQ.0) THEN
       CALL INTRD0(DTIME,KSTEP,NSTEP,NSSTEP)
       GO TO 1000
      ENDIF
C
C**** RESTART RUN:
C
      IF(IREST.EQ.1) THEN
C
C**** FIRST OF ALL READ ALL DATA UP TO THE END OF THE GENERAL DATA
C     AND COUNT THE NUMBER OF "INTERVAL DATA BLOCKS" UNTIL
C     THE ONE TO RESTART IS FOUND
C
       KOUNT=0
  200  CONTINUE
       READ(LUDAT,900) XCOMM,YCOMM
       IF(XCOMM.NE.'INTERVAL_') GOTO 200
       KOUNT=KOUNT+1
       IF(KOUNT.LT.ITIME) GOTO 200
C
C**** NOW LOOK FOR THE LAST "END_INTERVAL" CARD
C
  210  CONTINUE
       READ(LUDAT,900) XCOMM,YCOMM
       IF(XCOMM.NE.'END_INTER') GOTO 210
C
C**** IF THE LAST INTERVAL HAS BEEN FINISHED ( NSTEP < KSTEP )
C     IN A "RESTART CONTINUE" THE NEXT CARDS ARE SUPPOSED TO BE 
C     THE DATA FOR A NEW INTERVAL
C
       IF(ISKIP.EQ.0) THEN
        KSTEP=ISTEP+1
        IF(KSTEP.LE.NSTEP) GO TO 1000
C
        IREST=2
C
        CALL INTRD0(DTIME,KSTEP,NSTEP,NSSTEP)
        GO TO 1000
C
       ENDIF
C
C**** IN A "RESTART SKIP" THE NEXT CARDS ARE SUPPOSED TO BE 
C     THE NEW STRATEGY FOR THE NEXT INTERVAL 
C
       IF(ISKIP.EQ.1)THEN
C
        CALL INTRD0(DTIME,KSTEP,NSTEP,NSSTEP)
        GO TO 1000
C
       ENDIF
C
       GOTO 1000
C
      ENDIF
C
C**** END OF SEARCH
C
 1000 CONTINUE
C
C**** IF A NEW INTERVAL IS STARTING UPDATE SOME VALUES
C
      IF(KSTEP.EQ.1) THEN
       ITIME=ITIME+1
       DO IFUNC=1,NFUNC
        TFICT(IFUNC)=0.0
       ENDDO 
       KUNLD=0
      ENDIF
C
C**** SOME CONTROL
C
      IF(ITERME.GE.0) THEN       ! uni or bidirectional coupled problems
       IF(NSSTEP.GT.1.AND.NFURESM.GT.0)
     .  CALL RUNMEN('WARNING: INCREMENTATION AND FUTURE ANALYSIS COULD
     .  NOT BE COMPATIBLE')
      ENDIF
C
      RETURN
  900 FORMAT(A9,A123)
      END
