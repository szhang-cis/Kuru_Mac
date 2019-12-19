      SUBROUTINE RSTAR5T(TEMPIT,DISTOT,ELPRET,ELVART,TLOADT)
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
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION TEMPIT(NTOTVT,2),             DISTOT(NTOTVT,3)
      DIMENSION ELPRET(NPREVT),               ELVART(NSTATT)
      DIMENSION TLOADT(NTOTVT,2)
C
      CALL CPUTIMT(TIME1T)
C
C**** OPEN PREVIOUS ANALYSIS FILE
C
      OPEN(UNIT=LUPANT,FILE=CNT,STATUS='OLD',ACCESS='DIRECT',
     .     FORM='UNFORMATTED',RECL=LENRCT)
C
C**** READ FROM CONVERGED AREA (PREVIOUS ANALYSIS) ELPRE & ELVAR
C
      DO IELEMT=1,NELEMT
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
       IF(NMEMO.EQ.1)
     .  CALL DATPANT(ELPRET,    2,    2)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATPANT(ELVART,    3,    2)
C
C**** WRITE ELPRE & ELVAR ( CURRENT & LAST CONVERGED )
C
       IF(NMEMO.EQ.1)
     .  CALL DATBAST(ELPRET,    2,    1)
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    3,    1)         ! current
       IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     .  CALL DATBAST(ELVART,    5,    1)         ! last converged
      ENDDO
C
C**** READ ALL MATRICES
C
      CALL DATPANT(DISTOT,   4,     2)
C
      CALL DATPANT(TEMPIT,   5,     2)               ! instead of disprt
C
      CALL DATPANT(TLOADT,   8,     2)

c                      to be improved !!!!!!!!!!!!!!!!
      DO INDEXT=1,NTOTVT
       DO J=1,2
        TLOADT(INDEXT,J)=0.0
       ENDDO
      ENDDO


C
C**** CLOSE PREVIOUS ANALYSIS FILE
C
      CLOSE(LUPANT)
C
C**** PRINTS PREVIOUS RESULTS (temperatures)
C
      WRITE(LUREST,921)
      DO IPOINT=1,NPOINT
       ITOTVA=(IPOINT-1)*NDOFCT
       WRITE(LUREST,930) IPOINT,(DISTOT(ITOTVA+IDOFNT,1),
     .                   IDOFNT=1,NDOFNT)
      ENDDO
C
      WRITE(LUREST,940)
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
      CALL CPUTIMT(TIME2T)
      CPURST=CPURST+(TIME2T-TIME1T)
      RETURN
C
  921 FORMAT(//,1X,'NODE',4X,'PREVIOUS TEMPERATURE',/)
  930 FORMAT(I5,10X,3(E15.6,5X))
  940 FORMAT(/)
      END
