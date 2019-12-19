      PROGRAM VULCAN
C***********************************************************************
C*                                                                     *
C*    MASTER OF VULCAN-TM                                              *
C*                                                                     *
C*                       THERMOMECHANICAL ANALYSIS                     *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE        :: WORKS(:),WORKST(:)
      REAL*8, ALLOCATABLE        :: WORK1(:)
      INTEGER                    :: ERROR
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
      INCLUDE 'para_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'       ! thermal-mechanical
      INCLUDE 'nued_om.f'       ! thermal-microstructural
      INCLUDE 'nuee_om.f'       ! mechanical-microstructural
      INCLUDE 'nuef_om.f'       ! thermal-flow
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
     
C
C**** Se decalara la interface para pasar variables dinamicas a subrutinas
C      
      INTERFACE
        INCLUDE 'moduls.inc'
        INCLUDE 'moduns.inc'
      END INTERFACE
C
      PARAMETER(
c    . MWORK1TM=max(MWORK1T,MWORK1) )                 ! sg
c    . MWORK1TM=(MWORK1T+MWORK1) )                    ! linux (1)
     . MWORK1TM=MWORK1 )                              ! linux (2)
C
C**** DATA-BASE VECTOR
C
C     Note: to chnage the dimension of IPRIN, change it in vul-*m*.f,
C           the "mod"*.f routines called in them and in addpri.f.
C
C
      DIMENSION IPRIN(50)
      DIMENSION IPRINT(50)
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
      CALL SETDAT       ! mechanical problem
      CALL SETDATT      ! thermal problem
      CALL SETDATC      ! coupling problem
C
C**** COUPLING, MECHANICAL & THERMAL FILES ASSIGNMENT
C
      CALL ASSIFI
      CALL ASSIFIT
      CALL ASSIFIC
C
      CALL FINDUS
      CALL FINDUST
      CALL FINDUSC
C
C**** CALL SYSTEM ROUTINES TO INITIALIZE ERRORS AND CPUTIME
C
      CALL CPUTIM(CPUIN)
      CALL CPUTIMT(CPUINT)
C
C**** CHECK IF THIS IS A FRESH OR A RESTART RUN
C
      CALL CHECK0
      CALL CHECK0T
C
C**** WRITE DIMENSIONS (in MB) OF PERMANENT & TEMPORARY ARRAYS
C
      WRITE(LURES,2900) 
     . MPRIN*8/1000000,MBUFFER*8/1000000,MWORK1*8/1000000
 2900 FORMAT(//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .             5X,'==============================  ',/,
     . 15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB',//,
     . 15X,'TEMPORARY MEMORY OF MECHANICAL PROBLEM SHARED WITH',/,
     . 15X,'TEMPORARY MEMORY OF THERMAL PROBLEM')
C
      WRITE(LUREST,3900) 
     . MPRINT*8/1000000,MBUFFERT*8/1000000,MWORK1T*8/1000000
 3900 FORMAT(//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .             5X,'==============================  ',/,
     . 15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB',//,
     . 15X,'TEMPORARY MEMORY OF THERMAL PROBLEM SHARED WITH',/,
     . 15X,'TEMPORARY MEMORY OF MECHANICAL PROBLEM')
C
C**** WRITE MAIN DIMENSION PARAMETERS (mechanical & thermal problems)
C
      IF(MWRITEM.EQ.1) STOP
      IF(MWRITEM.EQ.2) THEN
      WRITE(LURES,1900)
 1900 FORMAT(//'* * * MAXIMUM DIMENSIONS * * *',/)
      WRITE(LURES,1901) MDATA
 1901 FORMAT(' MDATA    =',I10)
      WRITE(LURES,1902) MPREV
 1902 FORMAT(' MPREV    =',I10)
      WRITE(LURES,1903) MSTAT
 1903 FORMAT(' MSTAT    =',I10)
      WRITE(LURES,1904) MMATX
 1904 FORMAT(' MMATX    =',I10)
      WRITE(LURES,1905) MPRIN
 1905 FORMAT(' MPRIN    =',I10)
      WRITE(LURES,1906) MBUFFER
 1906 FORMAT(' MBUFFER  =',I10)
      WRITE(LURES,1908) MWORK1
 1908 FORMAT(' MWORK1   =',I10,//)
      ENDIF
C
      IF(MWRITE.EQ.1) STOP
      IF(MWRITE.EQ.2) THEN
      WRITE(LUREST,900)
  900 FORMAT(//'* * * MAXIMUM DIMENSIONS * * *',/)
      WRITE(LUREST,901) MDATAT
  901 FORMAT(' MDATAT    =',I10)
      WRITE(LUREST,902) MPREVT
  902 FORMAT(' MPREVT    =',I10)
      WRITE(LUREST,903) MSTATT
  903 FORMAT(' MSTATT    =',I10)
      WRITE(LUREST,904) MMATXT
  904 FORMAT(' MMATXT    =',I10)
      WRITE(LUREST,905) MPRINT
  905 FORMAT(' MPRINT    =',I10)
      WRITE(LUREST,906) MBUFFERT
  906 FORMAT(' MBUFFERT  =',I10)
      WRITE(LUREST,908) MWORK1T
  908 FORMAT(' MWORK1T   =',I10,//)
      ENDIF
C
C**** OPEN RESTART FILE
C
      CALL RSOPEN
      CALL RSOPENT
C
C**** READ CONTROL DATA FOR THE MECHANICAL PROBLEM
C
      IF(IREST.EQ.0) THEN
C
C**** NEW RUN:
C                                              - READ FROM INPUT FILE
       CALL CONTRO
C                                              - WRITE TO DATA BASE
       CALL RSTAR1(1)
C
      ELSE
C
C**** RESTART RUN: (not implemented yet !)
C                                              - READ FROM DATA BASE
       CALL RSTAR1(2)
C                                              - CHECK RESTART FILE
       CALL RSCHEK
C
      ENDIF
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
C**** RESTART RUN: (not implemented yet !!)
C                                              - READ FROM DATA BASE
       CALL RSTAR1T(2)
C                                              - CHECK RESTART FILE
       CALL RSCHEKT
C
      ENDIF
C
C**** ESTABLISH THE COUPLING CONTROL PARAMETERS
C
      IF(IREST.NE.IRESTT)
     . CALL RUNEND('ERROR: IREST NE IRESTT  ')
C
      IF(IREST.EQ.0) THEN
C
C**** NEW RUN:
C                                              - READ FROM INPUT FILE
       CALL COUINC
C
      ELSE
C
C**** RESTART RUN: (not implemented yet !!)
C                                              - READ FROM DATA BASE
       CALL COUINC
C
      ENDIF
C
C**** SET UP POINTERS FOR MECHANICAL DATA BASE FILE
C
      CALL RSSET0
C
C**** SET UP POINTERS FOR THERMAL DATA BASE FILE
C
      CALL RSSET0T
C
C**** OPEN POSTPROCESS FILE IF NECESSARY
C
      CALL PSOPEN
      CALL PSOPENT
C
C**** CALCULATE MEMORY POINTER FOR PERMANENT ARRAYS
C
      CALL ADDPRI(IPRIN)
      CALL ADDPRIT(IPRINT)
C
C**** CALCULATE MEMORY POINTER FOR TEMPORARY ARRAYS
C
      CALL ADDWOR
      CALL ADDWORT



C**** RESERVE MEMORY FOR WORKS,WORKT,WORK1
C
      ERROR=0
      ALLOCATE(WORKST(LPRINT),STAT=ERROR)
      WRITE(LUREST,*) ' WORKST(LPRINT) =WORKST(',LPRINT,')'
      WRITE(LUREST,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRINT/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKST IN vul-tm")

      ALLOCATE(WORKS(LPRIN),STAT=ERROR)
      WRITE(LURES,*) ' WORKS(LPRIN) =WORKS(',LPRIN,')'
      WRITE(LURES,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRIN/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKS IN vul-tm")

      ERROR=0
      LWOR1=MAX(LWOR1,LWOR1T)
      ALLOCATE(WORK1(LWOR1),STAT=ERROR)
      WRITE(LURES,*) ' WORK1(LWOR1) =WORK1(',LWOR1,')'
      WRITE(LURES,*) 'REQUIRED SOLVER MEMORY   (TEMPORARY) ='
     .                 ,LWOR1/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0)
     .  CALL RUNEND("I CAN'T ALLOCATE WORK1 IN vul-t.f")


C
C**** STARTS PROGRAM COMPUTATION
C
C***********************************************************************
C
C     ICOSTA DEFINES AN ASPECT OF THE STAGGERED STRATEGY USED
C
C           =0 :  CONVERGED STAGGERED STRATEGY
C           =1 :  NOT CONVERGED STAGGERED STRATEGY
C
C***********************************************************************
C
      IF(ICOSTA.EQ.0) THEN
       CALL MODULS(WORKS, IPRIN, MPRIN,
     .             WORKST,IPRINT,MPRINT,
     .             WORK1)
      ELSE
       call runend('icosta=1 not implemented yet')
       CALL MODUNS(WORKS, IPRIN, MPRIN,
     .             WORKST,IPRINT,MPRINT,
     .             WORK1)
      ENDIF
C
      END
