      PROGRAM VULCAN
C***********************************************************************
C*                                                                     *
C*    MASTER OF VULCAN-S                                               *
C*                                                                     *
C*                       SPECIES CONSERVARTION ANALYSIS ONLY           *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE            :: WORKSS(:),WORK1S(:)
      INTEGER                        :: ERROR
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_oms.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microstructural
      INCLUDE 'nuee_om.f'       ! mechanical-microstructural
      INCLUDE 'nuef_om.f'       ! thermal-flow
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'prob_oms.f'

C
C**** Se decalara la interface para pasar variables dinamicas a subrutinas
C      
      INCLUDE 'moduss.inc'


C
C**** DATA-BASE VECTOR
C
C
      DIMENSION IPRINS(50)
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
      CALL SETDATS
C
C**** THERMAL FILES ASSIGNMENT
C
      CALL ASSIFIS
C
      CALL FINDUSS
C
C**** CLOSE UNIT=6
C
C     Note: this CLOSE is only necessary for CONVEX machines; why?
C
      if(nmachis.eq.1) then
       CLOSE(6)
      endif
C
C**** CALL SYSTEM ROUTINES TO INITIALIZE ERRORS AND CPUTIME
C
      CALL CPUTIMS(CPUINS)
C
C**** CHECK IF THIS IS A FRESH OR A RESTART RUN
C
      CALL CHECK0S
C
C**** WRITE DIMENSIONS (in MB) OF PERMANENT & TEMPORARY ARRAYS
C
      WRITE(LURESS,4900) 
     . MPRINS*8/1000000,MBUFFERS*8/1000000,MWORK1S*8/1000000
 4900 FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .             5X,'==============================  ',/,
     . 15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB')
C 
C**** WRITE MAIN DIMENSION PARAMETERS
C
      IF(MWRITES.EQ.1) STOP
      IF(MWRITES.EQ.2) THEN
       WRITE(LURESS,5900)
 5900  FORMAT(//'* * * MAXIMUM DIMENSIONS * * *',/)
       WRITE(LURESS,5901) MDATAS
 5901  FORMAT(' MDATAS    =',I10)
       WRITE(LURESS,5902) MPREVS
 5902  FORMAT(' MPREVS    =',I10)
       WRITE(LURESS,5903) MSTATS
 5903  FORMAT(' MSTATS    =',I10)
       WRITE(LURESS,5904) MMATXS
 5904  FORMAT(' MMATXS    =',I10)
       WRITE(LURESS,5905) MPRINS
 5905  FORMAT(' MPRINS    =',I10)
       WRITE(LURESS,5906) MBUFFERS
 5906  FORMAT(' MBUFFERS  =',I10)
       WRITE(LURESS,5908) MWORK1S
 5908  FORMAT(' MWORK1S   =',I10,//)
      ENDIF
C
C**** OPEN RESTART FILE
C
      CALL RSOPENS
C
C**** READ CONTROL DATA FOR THE THERMAL PROBLEM
C
      IF(IRESTS.EQ.0) THEN
C
C**** NEW RUN:
C                                              - READ FROM INPUT FILE
       CALL CONTROS
C                                              - WRITE TO DATA BASE
       CALL RSTAR1S(1)
C
      ELSE
C
C**** RESTART RUN:
C                                              - READ FROM DATA BASE
       CALL RSTAR1S(2)
C                                              - CHECK RESTART FILE
       CALL RSCHEKS
C
      ENDIF
C
C**** SET UP POINTERS FOR DATA BASE FILE
C
      CALL RSSET0S
C
C**** OPEN POSTPROCESS FILE IF NECESSARY
C
      CALL PSOPENS
C
C**** CALCULATE MEMORY POINTER FOR PERMANENT ARRAYS
C
      CALL ADDPRIS(IPRINS)
C
C**** CALCULATE MEMORY POINTER FOR TEMPORARY ARRAYS
C
      CALL ADDWORS
C
C**** RESERVE MEMORY FOR WORKSS
C
      ERROR=0
      ALLOCATE(WORKSS(LPRINS),STAT=ERROR)
      WRITE(LURESS,*) ' WORKSS(LPRINS) =WORKSS(',LPRINS,')'
      WRITE(LURESS,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRINS/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKSS IN vul-s")

      ERROR=0
      ALLOCATE(WORK1S(LWOR1S),STAT=ERROR)
      WRITE(LURESS,*) ' WORK1S(LWOR1S) =WORK1S(',LWOR1S,')'
      WRITE(LURESS,*) 'REQUIRED SOLVER MEMORY   (TEMPORARY) ='
     .                 ,LWOR1S/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0)
     .  CALL RUNEND("I CAN'T ALLOCATE WORK1S IN vul-s.f")



C
C**** STARTS PROGRAM COMPUTATION
C
      CALL MODUSS(WORKSS,WORK1S,IPRINS,MPRINS)
C
      END
