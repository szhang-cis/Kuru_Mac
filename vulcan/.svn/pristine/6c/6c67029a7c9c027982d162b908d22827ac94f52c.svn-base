      PROGRAM VULCAN
C***********************************************************************
C*                                                                     *
C*    MASTER OF VULCAN-T                                               *
C*                                                                     *
C*                       THERMAL ANALYSIS ONLY                         *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE  :: WORKST(:),WORK1T(:)
      INTEGER                    :: ERROR
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microstructural
      INCLUDE 'nuee_om.f'       ! mechanical-microstructural
      INCLUDE 'nuef_om.f'       ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'

C
C**** Se decalara la interface para pasar variables dinamicas a subrutinas
C      
      INCLUDE 'modult.inc'

C
C**** DATA-BASE VECTOR
C
      DIMENSION IPRINT(50)
C
C
C**** NO MECHANICAL COUPLING
C
      ITERME=-1
C
C**** NO MICROSTRUCTURAL COUPLING
C
      IMICR=0
      IMICRM=0
C
C**** NO FLOW COUPLING
C
      ITERMEF=0
C
C**** SETS DATA
C
      CALL SETDATT
C
C**** THERMAL FILES ASSIGNMENT
C
      CALL ASSIFIT
C
      CALL FINDUST
C
C**** CLOSE UNIT=6
C
C     Note: this CLOSE is only necessary for CONVEX machines; why?
C
      if(nmachi.eq.1) then
       CLOSE(6)
      endif
C
C**** CALL SYSTEM ROUTINES TO INITIALIZE ERRORS AND CPUTIME
C
      CALL CPUTIMT(CPUINT)
C
C**** CHECK IF THIS IS A FRESH OR A RESTART RUN
C
      CALL CHECK0T
C
C**** WRITE DIMENSIONS (in MB) OF PERMANENT & TEMPORARY ARRAYS
C
      WRITE(LUREST,3900) 
     . MPRINT*8/1000000,MBUFFERT*8/1000000,MWORK1T*8/1000000
 3900 FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .             5X,'==============================  ',/,
     . 15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB')
C 
C**** WRITE MAIN DIMENSION PARAMETERS
C
      IF(MWRITE.EQ.1) STOP
      IF(MWRITE.EQ.2) THEN
       WRITE(LUREST,900)
  900  FORMAT(//'* * * MAXIMUM DIMENSIONS * * *',/)
       WRITE(LUREST,901) MDATAT
  901  FORMAT(' MDATAT    =',I10)
       WRITE(LUREST,902) MPREVT
  902  FORMAT(' MPREVT    =',I10)
       WRITE(LUREST,903) MSTATT
  903  FORMAT(' MSTATT    =',I10)
       WRITE(LUREST,904) MMATXT
  904  FORMAT(' MMATXT    =',I10)
       WRITE(LUREST,905) MPRINT
  905  FORMAT(' MPRINT    =',I10)
       WRITE(LUREST,906) MBUFFERT
  906  FORMAT(' MBUFFERT  =',I10)
       WRITE(LUREST,908) MWORK1T
  908  FORMAT(' MWORK1T   =',I10,//)
      ENDIF
C
C**** OPEN RESTART FILE
C
      CALL RSOPENT
C
C**** READ CONTROL DATA FOR THE THERMAL PROBLEM
C
      IF(IRESTT.EQ.0) THEN
C
C**** NEW RUN:
C                                              - READ FROM INPUT FILE
       CALL CONTROT
C                                              - WRITE TO DATA BASE
       CALL RSTAR1T(1)
C
      ELSE
C
C**** RESTART RUN:
C                                              - READ FROM DATA BASE
       CALL RSTAR1T(2)
C                                              - CHECK RESTART FILE
       CALL RSCHEKT
C
      ENDIF
C
C**** SET UP POINTERS FOR DATA BASE FILE
C
      CALL RSSET0T
C
C**** OPEN POSTPROCESS FILE IF NECESSARY
C
      CALL PSOPENT
C
C**** CALCULATE MEMORY POINTER FOR PERMANENT ARRAYS
C
      CALL ADDPRIT(IPRINT)
C
C**** CALCULATE MEMORY POINTER FOR TEMPORARY ARRAYS
C
      CALL ADDWORT
C
C**** RESERVE MEMORY FOR WORKSS AND WORK1T
C
      ERROR=0
      ALLOCATE(WORKST(LPRINT),STAT=ERROR)
      WRITE(LUREST,*) ' WORKST(LPRINT) =WORKST(',LPRINT,')'
      WRITE(LUREST,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRINT/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKST IN vul-t")

      ERROR=0
      ALLOCATE(WORK1T(LWOR1T),STAT=ERROR)
      WRITE(LUREST,*) ' WORK1T(LWOR1T) =WORKST(',LWOR1T,')'
      WRITE(LUREST,*) 'REQUIRED SOLVER MEMORY   (TEMPORARY) ='
     .                 ,LWOR1T/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0)
     .  CALL RUNEND("I CAN'T ALLOCATE WORK1T IN vul-t.f")



C
C**** STARTS PROGRAM COMPUTATION
C
      CALL MODULT(WORKST,WORK1T,IPRINT,MPRINT)
C
      END
