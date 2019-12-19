      SUBROUTINE DATBAST(ARRAYT,IFLAGT,ITASKT)
C***********************************************************************
C
C**** THIS ROUTINE INTERFACES VULCAN-DATA WITH DATABASE-FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE       :: BUFFERT(:), TEMPT(:)
      INTEGER, SAVE             :: INIT=0,NBUFFERT,NBUFFERTTMP
      INTEGER                   :: ERROR

      SAVE BUFFERT
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DIMENSION ARRAYT(*)
C
C**** BEGIN
C
      NUNITT=LUDTST                                ! .dts
C
      IF(IDATPT(IFLAGT,5).LT.0) THEN               ! DATA IS OUT-CORE
C
       NRECDT=(IELEMT-1)*NLENPT+IDATPT(IFLAGT,4)+1
       CALL RESTDKT(NUNITT,ITASKT,ARRAYT,IDATPT(IFLAGT,1)*2,LENRCT,
     .              NRECDT)
C
      ELSE                                         ! DATA IS IN-CORE
C
       IF(IFLAGT.EQ.12) THEN ! ACCESS TO DISTO
        NRECDT=1
       ELSE 
        NRECDT=NTOTVT+(IELEMT-1)*NWORPT+IDATPT(IFLAGT,3)+1
       ENDIF

C
C**** INICIALIZA LA VARIABLE BUFFERS
C

       IF(INIT.EQ.0) THEN
         NBUFFERT = (NRECDT+IDATPT(IFLAGT,1)-1)*64
         ERROR = 0
         ALLOCATE(BUFFERT(NBUFFERT),STAT=ERROR)
         IF (ERROR.NE.0) CALL RUNEND("CANT ALLOCATE BUFFERT IN datbast")
         INIT = 1
       ENDIF

C
C**** SE COMPRUEBA Y SE REDIMENSIONA SI ES NECESARIO LA VARIABLE BUFFER
C
      NBUFFERTTMP=(NRECDT+IDATPT(IFLAGT,1))
      IF(NBUFFERTTMP.GT.NBUFFERT) THEN
       ERROR=0
       IF (1.5*NBUFFERT.GT.NBUFFERTTMP) NBUFFERTTMP=1.5*NBUFFERT
       WRITE(LUREST,*) "REALLOCATE BUFFERT (FROM, TO, NEED)",
     .            NBUFFERT,NBUFFERTTMP,NRECDT+IDATPT(IFLAGT,1)
       WRITE(LUREST,*) 'REQUIRED DATABASE MEMORY (PERMANENT) ='
     .                 ,NBUFFERTTMP/(1024.0*1024.0)*8.0,'MB'
       ALLOCATE(TEMPT(NBUFFERTTMP),STAT=ERROR)
       IF(ERROR.NE.0)CALL RUNEND("CANT REALLOCATE BUFFERT IN datbast.f")
       TEMPT(1:NBUFFERT)=BUFFERT(:)
       CALL MOVE_ALLOC(FROM=TEMPT, TO=BUFFERT)
       NBUFFERT=NBUFFERTTMP
      END IF


       CALL RESTVM(BUFFERT,ITASKT,ARRAYT,IDATPT(IFLAGT,1),NRECDT)
C
      ENDIF
C
      RETURN
      END
