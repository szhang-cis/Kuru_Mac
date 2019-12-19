      SUBROUTINE RSTAR4S(TEMPIS,DISTOS,ELPRES,ELVARS)
C***********************************************************************
C
C**** THIS ROUTINE REFORMS THE FOLLOWING TASKS:
C
C     NEW RUN OR RESTART: WRITES THE RESULTS OF THE CURRENT TIME STEP
C                         INTO *.fan file (future analysis file)
C
C     OPTIONS IN THE INPUT DATA FILE (see check0t.f):
C     1) START,NOINITIAL,FUTURE_ANALYSIS=FORM1                 (new run)
C     2) START,INITIAL,FUTURE_ANALYSIS=FORM1                   (new run)
C     3) START,PREVIOUS_ANALYSIS,time of the previous analysis, \
C        FUTURE_ANALYSIS=FORM1                                 (restart)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION TEMPIS(NTOTVS,2),             DISTOS(NTOTVS,3)
      DIMENSION ELPRES(NPREVS),               ELVARS(NSTATS)


      call runends('ERROR: in rstar4s')   ! to be improved-see rstar4t.f

C
c     CALL CPUTIMS(TIME1S)
C
C**** OPEN PREVIOUS ANALYSIS FILE
C
c     OPEN(UNIT=LUPANT,FILE=CNT,STATUS='OLD',ACCESS='DIRECT',
c    .     FORM='UNFORMATTED',RECL=LENRCT)
C
C**** READ FROM PROCESS AREA AND WRITE TO CONVERGED AREA ELPRE & ELVAR
C
c     IF(IFLAGT.EQ.1) THEN
c      DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
c       IF(NMEMO.EQ.1)
c    .   CALL DATBAST(ELPRET,    2,    2)
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .   CALL DATBAST(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR
C
c       IF(NMEMO.EQ.1)
c    .   CALL DATRSTT(ELPRET,    2,    1)
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .   CALL DATRSTT(ELVART,    3,    1)
c      ENDDO
c     ENDIF
C
C**** READ FROM CONVERGED AREA (PREVIOUS ANALYSIS) ELPRE & ELVAR
C
c     IF(IFLAGT.EQ.2) THEN
c      DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
c       IF(NMEMO.EQ.1)
c    .   CALL DATRSTT(ELPRET,    2,    2)
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .   CALL DATRSTT(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR ( CURRENT & LAST CONVERGED )
C
c       IF(NMEMO.EQ.1)
c    .   CALL DATBAST(ELPRET,    2,    1)
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .   CALL DATBAST(ELVART,    3,    1)         ! current
c       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .   CALL DATBAST(ELVART,    5,    1)         ! last converged
c      ENDDO
c     ENDIF
C
C**** READ OR WRITE ALL MATRICES
C
c     CALL DATRSTT(DISTOT,   4,IFLAGT)
C
c     CALL DATRSTT(TEMPIT,   5,IFLAGT)               ! instead of disprt
C
C**** CLOSE PREVIOUS ANALYSIS FILE
C
c     CLOSE(LUPANT)
C
C**** PRINTS PREVIOUS RESULTS (only temperatures)
C
c     IF(INITIT.EQ.2) THEN
c      WRITE(LUREST,921)
c      DO IPOINT=1,NPOINT
c       ITOTVA=(IPOINT-1)*NDOFCT
c       WRITE(LUREST,930) IPOINT,(DISTOT(ITOTVA+IDOFNT,1),
c    .                    IDOFNT=1,NDOFNT)
c      ENDDO
C
c      WRITE(LUREST,940)
c     ENDIF
C
c     CALL CPUTIMS(TIME2S)
c     CPURSS=CPURSS+(TIME2S-TIME1S)
      RETURN
C
  921 FORMAT(//,1X,'NODE',4X,'PREVIOUS TEMPERATURE',/)
  930 FORMAT(I5,10X,3(E15.6,5X))
  940 FORMAT(/)
      END
