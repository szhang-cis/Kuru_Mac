      SUBROUTINE RSTAR5(DISPR,DISTO,ELPRE,ELVAR,LNODS,TLOAD,LACTI)
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
      DIMENSION LNODS(NNODE,NELEM), TLOAD(NTOTV,2)
      DIMENSION LACTI(NELEM)
C
      CALL CPUTIM(TIME1)
C
C**** OPEN PREVIOUS ANALYSIS FILE
C
      OPEN(UNIT=LUPAN,FILE=CN,STATUS='OLD',ACCESS='DIRECT',
     .     FORM='UNFORMATTED',RECL=LENRC)
C
C**** READ FROM CONVERGED AREA (PREVIOUS ANALYSIS) ELPRE & ELVAR
C
      DO IELEM=1,NELEM
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
       IF(NMEMOM.EQ.1)
     .  CALL DATPAN(ELPRE,    2,    2)
       CALL DATPAN(ELVAR,    3,    2)
C
C**** WRITE ELPRE & ELVAR ( CURRENT & LAST CONVERGED )
C
       IF(NMEMOM.EQ.1)
     .  CALL DATBAS(ELPRE,    2,    1)
       CALL DATBAS(ELVAR,    3,    1)         ! current
       CALL DATBAS(ELVAR,    5,    1)         ! last converged
      ENDDO
C
C**** READ ALL MATRICES
C
      CALL DATPAN(DISTO,   4,    2)
C
      CALL DATPAN(TLOAD,   8,    2)
C
      IF(NACTI.EQ.1) THEN
       CALL DATPAN(LACTI,  12,    2)
      ENDIF                      ! nacti.eq.1
C
C**** CLOSE PREVIOUS ANALYSIS FILE
C
      CLOSE(LUPAN)
C
C**** PRINTS PREVIOUS RESULTS (displacements)
C
      WRITE(LURES,921)
      DO IPOIN=1,NPOIN
       ITOTVA=(IPOIN-1)*NDOFC
       WRITE(LURES,930) IPOIN,(DISTO(ITOTVA+IDOFN,1),IDOFN=1,NDOFN)
      ENDDO
C
      WRITE(LURES,940)
C
C**** PRINTS PREVIOUS RESULTS (internal variables)  to be improved !!!
C
c     DO IELEM=1,NELEM
C
C**** READ ELPRE & ELVAR ( CURRENT )
C
c      IF(NMEMOM.EQ.1)
c    .  CALL DATBAST(ELPRET,    2,    1)
c      CALL DATBAST(ELVART,    3,    2)         ! current
C
c      if(ielem.eq.1) then
c       write(lures,*) 'elvar  current - ielem=',ielem
c       do i=istat(2),istat(3)-1
c        write(lures,*) elvar(i)
c       enddo
c      endif
c     ENDDO
C
C**** PRINTS PREVIOUS RESULTS (internal variables)
C
c     DO IELEM=1,NELEM
C
C**** READ ELPRE & ELVAR ( LAST CONVERGED )        to be improved !!!
C
c      CALL DATBAST(ELVART,    3,    2)         ! current
C
c      if(ielem.eq.1) then
c       write(lures,*) 'elvar  last converged - ielem=',ielem
c       do i=istat(2),istat(3)-1
c        write(lures,*) elvar(i)
c       enddo
c      endif
c     ENDDO
C
      CALL CPUTIM(TIME2)
      CPURS=CPURS+(TIME2-TIME1)
      RETURN
C
  921 FORMAT(//,1X,'NODE',4X,'PREVIOUS DISPLACEMENTS',/)
  930 FORMAT(I5,10X,3(E15.6,5X))
  940 FORMAT(/)
      END
