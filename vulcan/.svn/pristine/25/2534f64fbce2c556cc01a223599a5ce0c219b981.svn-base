      SUBROUTINE RSTAR4(DISPR,DISTO,ELPRE,ELVAR,TLOAD,LACTI)
C***********************************************************************
C
C**** THIS ROUTINE REFORMS THE FOLLOWING TASKS:
C
C     NEW RUN OR RESTART: WRITES THE RESULTS OF THE CURRENT TIME STEP
C                         INTO *.fan file (future analysis file)
C
C     OPTIONS IN THE INPUT DATA FILE (see check0.f):
C     1) START,NOINITIAL,FUTURE_ANALYSIS=FORM1                 (new run)
C     2) START,INITIAL,FUTURE_ANALYSIS=FORM1                   (new run)
C     3) START,PREVIOUS_ANALYSIS,time of the previous analysis, \
C        FUTURE_ANALYSIS=FORM1                                 (restart)
C
C     Note: DISPR is not necessary to be written/read
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION DISPR(NTOTV,NDISR), DISTO(NTOTV,NDISO)
      DIMENSION ELPRE(NPREV),       ELVAR(NSTAT)
      DIMENSION TLOAD(NTOTV,2)
      DIMENSION LACTI(NELEM)
C
      CALL CPUTIM(TIME1)
C
C**** OPEN FUTURE ANALYSIS FILE
C
      OPEN(UNIT=LUFAN,FILE=CO,STATUS='OLD',ACCESS='DIRECT',
     .     FORM='UNFORMATTED',RECL=LENRC)
C
C**** READ FROM PROCESS AREA AND WRITE TO CONVERGED AREA ELPRE & ELVAR
C
      DO IELEM=1,NELEM
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
       IF(NMEMOM.EQ.1)
     .  CALL DATBAS(ELPRE,    2,    2)
       CALL DATBAS(ELVAR,    3,    2)
C
C**** WRITE ELPRE & ELVAR
C
       IF(NMEMOM.EQ.1)
     .  CALL DATRST(ELPRE,    2,    1)
       CALL DATRST(ELVAR,    3,    1)
      ENDDO
C
C**** READ OR WRITE ALL MATRICES
C
      CALL DATRST(DISTO,   4,    1)
C
      CALL DATRST(TLOAD,   8,    1)
C
      IF(NACTI.EQ.1) THEN
       CALL DATRST(LACTI,  12,    1)
      ENDIF                      ! nacti.eq.1
C
C**** CLOSE FUTURE ANALYSIS FILE
C
      CLOSE(LUFAN)
C
      CALL CPUTIM(TIME2)
      CPURS=CPURS+(TIME2-TIME1)
      RETURN
C
      END
