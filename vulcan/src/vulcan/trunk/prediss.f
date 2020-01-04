      SUBROUTINE PREDISS(DISTOS,NCOMPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE INITIAL VALUES FOR DISTO
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'inpo_oms.f'
C
      COMMON/PRFILES/JTAPES,IFLAGS
C
      DIMENSION DISTOS(NDOFCS,*)
      DIMENSION SUBTAS(12)
C
      IOLDE=1                       ! new reading
C
c     IERR1T=1
c     IF(JTAPET.NE.LUDATT.AND.JTAPET.NE.LUINIT) GO TO 1000
c     IF(JTAPET.EQ.LUINIT) THEN       ! external .ini file
c      OPEN(UNIT=LUINIT,FILE=CD1T,STATUS='OLD',ERR=1000)
c      IERR1T=0
C
c1000  IF(IERR1T.NE.0) THEN
c       IF(IERR1T.EQ.1) WRITE(LUREST,901)
c       CALL RUNENDT('PREDIST: ERROR IN OPENING FILES')
c      ENDIF
c     ENDIF
C
C**** READ FOR NUMBER OF NODES WITH INITIAL CONDITIONS
C
      IF(IOLDE.EQ.0) THEN
C
C**** OLD READING (IOLDE=0)
C
c      READ(JTAPET,900) NPOI1T
c      WRITE(LUREST,905) NPOI1T
c      IF(NPOINT.LT.NPOI1T) THEN
c       WRITE(LUREST,910)
c       CALL RUNENDT('PREDIST:ERROR WHEN READING         ')
c      ENDIF
C
      ELSE
C
C**** NEW READING (IOLDE=1)
C
       NPRINS=1
       CALL LISTENS('PREDISS',NPRINS,JTAPES)
       IF(PARAMS(2).NE.0.0)
     .  CALL RUNENDS('PREDISS: ERROR WHEN READING THE NUMBER OF POINTS
     .                         WITH INITIAL TEMPERATURE')
       NPOI1S=DINT(PARAMS(1))
       WRITE(LURESS,905) NPOI1S
       IF(NPOINS.LT.NPOI1S) THEN
        WRITE(LURESS,910)
        CALL RUNENDS('PREDISS: ERROR WHEN READING THE NUMBER OF POINTS
     .                         WITH INITIAL TEMPERATURE')
       ENDIF
C
      ENDIF                   ! iolde.eq.0
C
C**** READ TITLE (ONLY WITH IOLDE=0)
C
      IF(IOLDE.EQ.0) THEN
       READ(JTAPES,915) SUBTAS
       WRITE(LURESS,916) SUBTAS
      ENDIF
C
C**** READ NODAL VARIABLES
C
      DO IPOINS=1,NPOI1S
       IF(IOLDE.EQ.0) THEN
        READ(JTAPES,920) JPOINS, (DISTOS(IDOFNS,JPOINS),IDOFNS=1,NDOFNS)
       ELSE
        CALL LISTENS('PREDISS',NPRINS,JTAPES)
        JPOINS=INT(PARAMS(1))
        DO IDOFNS=1,NDOFNS
         DISTOS(IDOFNS,JPOINS)=PARAMS(1+IDOFNS)
        ENDDO
       ENDIF
      ENDDO
C
C**** ECHO NODAL VALUES FOR THE WHOLE MESH
C
      IF(NCOMPS.EQ.1) WRITE(LURESS,921)
      IF(NCOMPs.EQ.2) WRITE(LURESS,922)
      IF(NCOMPS.EQ.3) WRITE(LURESS,923)
      DO IPOINS=1,NPOINS
        WRITE(LURESS,930)  IPOINS,(DISTOS(IDOFNS,IPOINS),
     .                             IDOFNS=1,NDOFNS)
      ENDDO
C
      WRITE(LURESS,940)
C
      RETURN
C
  900 FORMAT(I5)
  901 FORMAT(' ERROR IN OPENING INITIAL TEMPERATURE INPUT FILE 243',
     . '(53 LINUX)')
  905 FORMAT(/,5X,'NUMBER OF NODES WITH INITIAL DATA',I5,/)
  910 FORMAT(//,
     .    'ERROR: NUMBER OF INPUT POINTS GT TOTAL NO. POINTS')
  915 FORMAT(12A6)
  916 FORMAT(5X,12A6)
  920 FORMAT(I5,3E15.6)
  921 FORMAT(1X,'NODE',4X,'INITIAL TEMPERATURE',/)
  922 FORMAT(1X,'NODE',8X,'INITIAL VELOCITY',/)
  923 FORMAT(1X,'NODE',4X,'INITIAL ACCELERATION',/)
  930 FORMAT(I5,10X,3(E15.6,5X))
  940 FORMAT(/)
      END