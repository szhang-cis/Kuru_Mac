      PROGRAM VULCAN
C***********************************************************************
C*                                                                     *
C*    MASTER OF VULCAN-TS                                              *
C*                                                                     *
C*                       THERMOMICROSTRUCTURAL ANALYSIS                *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE         :: WORKST(:),WORKSS(:)
      REAL*8, ALLOCATABLE         :: WORK1T(:)
      INTEGER                     :: ERROR
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
      INCLUDE 'para_oms.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
      INCLUDE 'addi_oms.f'
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
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'prob_oms.f'
C
C**** Se decalara la interface para pasar variables dinamicas a subrutinas
C      
      INCLUDE 'modults.inc'

C
      PARAMETER(
c    . MWORK1TS=max(MWORK1T,MWORK1S) )        ! sg
     . MWORK1TS=MWORK1T )                     ! PC & linux
C
C**** DATA-BASE VECTOR
C
C
      DIMENSION IPRINT(50)
      DIMENSION IPRINS(50)
C
C**** NO MECHANICAL-MICROSTRUCTURAL COUPLING
C
      IMICRM=0
C
C**** NO MECHANICAL COUPLING
C
      ITERME=-1
C
C**** NO FLOW COUPLING
C
      ITERMEF=0
C
C**** SETS DATA
C
      CALL SETDATD    ! coupling problem
      CALL SETDATT    ! thermal problem
      CALL SETDATS    ! microstructural problem
C
C**** COUPLING (not now), MICROSTRUCTURAL & THERMAL FILES ASSIGNMENT
C
      CALL ASSIFIT
      CALL ASSIFIS
c     CALL ASSIFID
C
      CALL FINDUST
      CALL FINDUSS
c     CALL FINDUSD
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
      CALL CPUTIMS(CPUINS)
C
C**** CHECK IF THIS IS A FRESH OR A RESTART RUN
C
      CALL CHECK0T
      CALL CHECK0S
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
      IF(IEVFI.EQ.1) THEN
       WRITE(LURESS,4900)
     .  MPRINS*8/1000000,MBUFFERS*8/1000000,MWORK1S*8/1000000
 4900  FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .              5X,'==============================  ',/,
     .  15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     .  15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     .  15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB')
C
       IF(MWRITES.EQ.1) STOP
       IF(MWRITES.EQ.2) THEN
        WRITE(LURESS,5900)
 5900   FORMAT(//'* * * MAXIMUM DIMENSIONS * * *',/)
        WRITE(LURESS,5901) MDATAS
 5901   FORMAT(' MDATAS    =',I10)
        WRITE(LURESS,5902) MPREVS
 5902   FORMAT(' MPREVS    =',I10)
        WRITE(LURESS,5903) MSTATS
 5903   FORMAT(' MSTATS    =',I10)
        WRITE(LURESS,5904) MMATXS
 5904   FORMAT(' MMATXS    =',I10)
        WRITE(LURESS,5905) MPRINS
 5905   FORMAT(' MPRINS    =',I10)
        WRITE(LURESS,5906) MBUFFERS
 5906   FORMAT(' MBUFFERS  =',I10)
        WRITE(LURESS,5908) MWORK1S
 5908   FORMAT(' MWORK1S   =',I10,//)
       ENDIF
      ENDIF           ! ievfi.eq.1
C
C**** OPEN RESTART FILE
C
      CALL RSOPENT
      CALL RSOPENS
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
C**** READ CONTROL DATA FOR THE MICROSTRUCTURAL PROBLEM
C
      IF(IEVFI.EQ.1) THEN
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
      ENDIF             ! ievfi.eq.1
C
C**** ESTABLISH THE COUPLING CONTROL PARAMETERS
C
      IF(IRESTT.NE.IRESTS)
     . CALL RUNENDT('ERROR: IRESTT NE IRESTS ')
C
      IF(IRESTT.EQ.0) THEN
C
C**** NEW RUN:
C                                              - READ FROM INPUT FILE
       CALL COUIND
C
      ELSE
C
C**** RESTART RUN: (not implemented yet !!)
C                                              - READ FROM DATA BASE
       CALL COUIND
C
      ENDIF
C
C**** SET UP POINTERS FOR THERMAL DATA BASE FILE
C
      CALL RSSET0T
C
C**** SET UP POINTERS FOR MICROSTRUCTURAL DATA BASE FILE
C
      IF(IEVFI.EQ.1) CALL RSSET0S
C
C**** OPEN POSTPROCESS FILE IF NECESSARY
C
      CALL PSOPENT
      IF(IEVFI.EQ.1) CALL PSOPENS
C
C**** CALCULATE MEMORY POINTER FOR PERMANENT ARRAYS
C
      CALL ADDPRIT(IPRINT)
      IF(IEVFI.EQ.0) NMATSS=NMATST
      CALL ADDPRIS(IPRINS)
C
C**** CALCULATE MEMORY POINTER FOR TEMPORARY ARRAYS
C
      CALL ADDWORT
      IF(IEVFI.EQ.1) CALL ADDWORS
C
C**** RESERVE MEMORY FOR WORKSS
C
      ERROR=0
      ALLOCATE(WORKST(LPRINT),STAT=ERROR)
      WRITE(LUREST,*) ' WORKST(LPRINT) =WORKST(',LPRINT,')'
      WRITE(LUREST,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRINT/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKST IN vul-ts")

      ALLOCATE(WORKSS(LPRINS),STAT=ERROR)
      WRITE(LURESS,*) ' WORKSS(LPRINS) =WORKSS(',LPRINS,')'
      WRITE(LURESS,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRINS/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0) CALL RUNEND("CAN'T ALLOCATE WORKSS IN vul-ts")

      ERROR=0
      LWOR1T=MAX(LWOR1T,LWOR1S)
      ALLOCATE(WORK1T(LWOR1T),STAT=ERROR)
      WRITE(LURESS,*) ' WORK1T(LWOR1T) =WORK1T(',LWOR1T,')'
      WRITE(LURESS,*) 'REQUIRED SOLVER MEMORY   (TEMPORARY) ='
     .                 ,LWOR1T/(1024.0*1024.0)*8.0,'MB'
      IF (ERROR.NE.0)
     .  CALL RUNEND("I CAN'T ALLOCATE WORK1T IN vul-ts.f")

C
C**** STARTS PROGRAM COMPUTATION
C
      CALL MODULTS(WORKST,IPRINT,MPRINT,
     .             WORKSS,IPRINS,MPRINS,
     .             WORK1T)
C
      END
