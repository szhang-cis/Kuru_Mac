      SUBROUTINE RSTAR5S(TEMPIS,DISTOS,ELPRES,ELVARS)
C***********************************************************************
C
C**** THIS ROUTINE REFORMS THE FOLLOWING TASKS:
C
C     RESTART RUN: READS THE RESULTS FROM THE PREVIOUS ANALYSIS AS
C                  INITIAL CONDITIONS FOR THE RESTART ANALYSIS
C
C     OPTIONS IN THE INPUT DATA FILE:
C     1) START,PREVIOUS_ANALYSIS,time of the previous analysis, \
C        FUTURE_ANALYSIS=FORM1                             (restart run)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION TEMPIS(NTOTVS,2),             DISTOS(NTOTVS,3)
      DIMENSION ELPRES(NPREVS),               ELVARS(NSTATS)
C
c     CALL CPUTIMS(TIME1S)

      call runends('ERROR: in rstar5s')   ! to be improved-see rstar5t.f


C
C**** OPEN PREVIOUS ANALYSIS FILE
C
c     OPEN(UNIT=LUPANT,FILE=CNT,STATUS='OLD',ACCESS='DIRECT',
c    .     FORM='UNFORMATTED',RECL=LENRCT)
C
c     DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
c      IF(NMEMO.EQ.1)
c    .  CALL DATPANT(ELPRET,    2,    2)
c      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .  CALL DATPANT(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR ( CURRENT & LAST CONVERGED )
C
c      IF(NMEMO.EQ.1)
c    .  CALL DATBAST(ELPRET,    2,    1)
c      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .  CALL DATBAST(ELVART,    3,    1)         ! current
c      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .  CALL DATBAST(ELVART,    5,    1)         ! last converged
c     ENDDO
C
C**** READ ALL MATRICES
C
c     CALL DATPANT(DISTOT,   4,     2)
C
c     CALL DATPANT(TEMPIT,   5,     2)               ! instead of disprt
C
C**** CLOSE PREVIOUS ANALYSIS FILE
C
c     CLOSE(LUPANT)
C
C**** PRINTS PREVIOUS RESULTS (temperatures)
C
c     WRITE(LUREST,921)
c     DO IPOINT=1,NPOINT
c      ITOTVA=(IPOINT-1)*NDOFCT
c      WRITE(LUREST,930) IPOINT,(DISTOT(ITOTVA+IDOFNT,1),
c    .                   IDOFNT=1,NDOFNT)
c     ENDDO
C
c     WRITE(LUREST,940)
C
C**** PRINTS PREVIOUS RESULTS (internal variables)  to be improved !!!
C
c     DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
c      IF(NMEMO.EQ.1)
c    .  CALL DATBAST(ELPRET,    2,    1)
c      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .  CALL DATBAST(ELVART,    3,    2)         ! current
C
c      if(ielemt.eq.1) then
c       write(lurest,*) 'elvar  current - ielemt=',ielemt
c       do i=istatt(2),istatt(3)-1
c        write(lurest,*) elvart(i)
c       enddo
c      endif
c     ENDDO
C
C**** PRINTS PREVIOUS RESULTS (internal variables)
C
c     DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( LAST CONVERGED )        to be improved !!!
C
c      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
c    .  CALL DATBAST(ELVART,    3,    2)         ! current
C
c      if(ielemt.eq.1) then
c       write(lurest,*) 'elvar  last converged - ielemt=',ielemt
c       do i=istatt(2),istatt(3)-1
c        write(lurest,*) elvart(i)
c       enddo
c      endif
c     ENDDO
C
c     CALL CPUTIMT(TIME2T)
c     CPURST=CPURST+(TIME2T-TIME1T)
      RETURN
C
  921 FORMAT(//,1X,'NODE',4X,'PREVIOUS TEMPERATURE',/)
  930 FORMAT(I5,10X,3(E15.6,5X))
  940 FORMAT(/)
      END
