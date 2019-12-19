C> ESTE2 es un comentario de prueba probosa
C> HOLA!
      PROGRAM VULCAN
C> ESTE es un comentario de prueba probosa
C***********************************************************************
C*                                                                     *
C*    MASTER OF VULCAN-M                                               *
C*                                                                     *
C*                       MECHANICAL ANALYSIS ONLY                      *
C*                                                                     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'      ! thermal-mechanical
      INCLUDE 'nued_om.f'      ! thermal-microstructural
      INCLUDE 'nuee_om.f'      ! mechanical-microstructural
      INCLUDE 'nuef_om.f'      ! thermal-flow
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
C**** DATA-BASE VECTOR
C
C     Note: to chnage the dimension of IPRIN, change it in vul-*m*.f,
C           the "mod"*.f routines called in them and in addpri.f.
C
      REAL*8, ALLOCATABLE :: WORKS(:),WORK1(:)
C
      DIMENSION IPRIN(50)      
C
C**** Se decalara la interface para pasar variables dinamicas a subrutinas
C      
      INCLUDE 'modulm.inc'
C
C
C**** NO THERMAL COUPLING
C
      ITERME=-1
C
C**** NO THERMAL-MICROSTRUCTURAL COUPLING
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
      CALL SETDAT
C
C**** MECHANICAL FILES ASSIGNMENT
C
      CALL ASSIFI
C
      CALL FINDUS
C
C**** CALL SYSTEM ROUTINES TO INITIALIZE ERRORS AND CPUTIME
C
      CALL CPUTIM(CPUIN)
C
C**** CHECK IF THIS IS A FRESH OR A RESTART RUN
C
      CALL CHECK0
C
C**** WRITE DIMENSIONS (in MB) OF PERMANENT & TEMPORARY ARRAYS
C
      WRITE(LURES,3900) 
     . (MPRIN/1000000)*8,(MBUFFER/1000000)*8,(MWORK1/1000000)*8
 3900 FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS (MAXIMUM):',/,
     .             5X,'==============================  ',/,
     . 15X,'REQUIRED ARRAYS MEMORY   (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED DATABASE MEMORY (PERMANENT) =',I8,' MB',/,
     . 15X,'REQUIRED SOLVER MEMORY   (TEMPORARY) =',I8,' MB')
C
C**** WRITE MAIN DIMENSION PARAMETERS
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
C**** OPEN RESTART FILE
C
      CALL RSOPEN
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
C**** RESTART RUN:
C                                              - READ FROM DATA BASE
        CALL RSTAR1(2)
C                                              - CHECK RESTART FILE
        CALL RSCHEK
C
      ENDIF
C
C**** SET UP POINTERS FOR DATA BASE FILE
C
      CALL RSSET0
C
C**** OPEN POSTPROCESS FILE IF NECESSARY
C
      CALL PSOPEN
C
C**** CALCULATE MEMORY POINTER FOR PERMANENT ARRAYS
C
      CALL ADDPRI(IPRIN)
C
C**** CALCULATE MEMORY POINTER FOR TEMPORARY ARRAYS
C
      CALL ADDWOR
C
C**** ALLOCATE MEMORY TO WORKS AND WORK1
C
      IF (.NOT. ALLOCATED(WORKS)) THEN
        IERROR=0
        ALLOCATE(WORKS(LPRIN),STAT=IERROR)
        WRITE(LURES,*) ' WORKS(LPRIN) =WORKS(',LPRIN,')'
       WRITE(LURES,*) 'REQUIRED ARRAYS MEMORY   (PERMANENT) ='
     .                 ,LPRIN/(1024.0*1024.0)*8.0,'MB'
      ENDIF
       IF (IERROR.NE.0) 
     .  CALL RUNEND("I CAN'T ALLOCATE WORKS IN vul-m.f") 
     
      IF (.NOT. ALLOCATED(WORK1)) THEN
        IERROR=0
        ALLOCATE(WORK1(LWOR1),STAT=IERROR)
        WRITE(LURES,*) ' WORK1(LWOR1) =WORKS(',LWOR1,')'
        WRITE(LURES,*) 'REQUIRED SOLVER MEMORY   (TEMPORARY) ='
     .                 ,LWOR1/(1024.0*1024.0)*8.0,'MB'
      ENDIF
       IF (IERROR.NE.0) 
     .  CALL RUNEND("I CAN'T ALLOCATE WORK1 IN vul-m.f")
               
C
C**** STARTS PROGRAM COMPUTATION
C
      CALL MODULM(WORKS,WORK1,IPRIN,MPRIN)
C
      END
