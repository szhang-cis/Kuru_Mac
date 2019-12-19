      SUBROUTINE COUINC
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT CONTROL COUPLING PARAMETERS
C     FOR THIS RUN
C
C     FOR CHANGING THE DEFAULT PARAMETERS
C
C     PARAMETERS         CONTROLLING WORDS
C
C     ICOSTA             'C_COS'
C     ITERME             'C_TYP'
C     ISTAGG             'STAGG'
C     INICOU             'NO_IN'
C     NPOINC             'C_POI'
C     NELEMC             'C_ELE'
C     COUFAC             'C_FAC'
C     CENKEL             'C_FAC'
C     NITERC             'C_FAC'
C     NPLCOU             'C_PLA'
C     FPLCOU             'C_PLA'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'inpo_om.f'
      INCLUDE 'prob_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
C
      PARAMETER (MCOMM=9)
      CHARACTER*5 COMMD(MCOMM)
      CHARACTER*6 NOTESC
C
      CHARACTER*8     TITLEC
      COMMON/TITLESCA/TITLEC(8)
C
      DATA COMMD/'C_COS','C_TYP','STAGG','NO_IN','C_POI','C_ELE',
     .           'C_FAC','C_PLA','END_C'/
C
      IBAND1=0
      IBAND2=0
C
      IERR1C=1
      OPEN(UNIT=LUDATC,FILE=CCOA,STATUS='OLD',ERR=1001)
      IERR1C=2
      OPEN(UNIT=LURESC,FILE=CCOB,STATUS='UNKNOWN',ERR=1001)
      IERR1C=0
C
 1001 IF(IERR1C.NE.0) THEN
       IF(IERR1C.EQ.1) WRITE(LURESC,901)
       IF(IERR1C.EQ.2) WRITE(LURESC,902)
      ENDIF
C
C**** READ HEADING CARD
C
   10 READ(LUDATC,900) NOTESC,(TITLEC(I),I=1,8)
      IF(NOTESC.EQ.'VULCAN') GO TO 99
      GO TO 10
C
   99 WRITE(LURESC,905) NOTESC,(TITLEC(I),I=1,8)
C
C**** READ CONTROL PARAMETERS
C
      IPRINC=0
      ITAPEC=LUDATC
      CALL LISTENC('COUINC',IPRINC,ITAPEC)
      IF(WORDS(1).NE.'COUPL') GO TO 2000
 1000 CALL LISTENC('COUINC',IPRINC,ITAPEC)
C
C**** IDENTIFY COMMAND
C
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GOTO 100
      ENDDO
      GO TO 2000
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GOTO (1,2,3,4,5,6,7,8,9), ICOMM
C
    1 CONTINUE
ctm   ICOSTA=INT(PARAM(1))
      icosta=0
      GO TO 1000
C
    2 CONTINUE
      ITERME=INT(PARAM(1))
C
      ITERMG=0
      ITERMP=0
      ITERMD=0
      ITERMF=0
      IF(ITERME.EQ.0) THEN                     ! unidirectional coupling
       IF(LARGET.NE.0) THEN
        LARGET=0
        CALL RUNMENT('WARNING: LARGET NE 0 IS NOT CONSIDERED')
       ENDIF
      ENDIF
      IF(ITERME.GT.0) THEN                     ! bidirectional coupling
       ITERMG=1                                ! defaults
       ITERMP=0              ! should be 1
       ITERMD=0              ! should be 1
       ITERMF=0              ! should be 1
       DO I=2,5
        IF(WORDS(I).EQ.'GAP_D') THEN         ! gap dependency
         ITERMG=INT(PARAM(I))
        ENDIF
        IF(WORDS(I).EQ.'MECHA') THEN         ! mechanical coupling terms
         ITERMP=INT(PARAM(I))
        ENDIF
        IF(WORDS(I).EQ.'DEFOR') THEN         ! deformed shape
         ITERMD=INT(PARAM(I))
         IF(ITERMD.EQ.0) THEN
          IF(LARGE.NE.0)
     .     CALL RUNMENT
     .     ('WARNING: LARGE STRAIN ANALYSIS WITH UNDEFORMED SHAPE')
          IF(LARGET.NE.0) THEN
           LARGET=0
           CALL RUNMENT('WARNING: LARGET NE 0 IS NOT CONSIDERED')
          ENDIF
         ELSE
          IF(LARGE.EQ.0) THEN
c          ITERMD=0                  ! option 1
           IF(LARGET.NE.0) THEN
            LARGET=0
            CALL RUNMENT('WARNING: LARGET NE 0 IS NOT CONSIDERED')
           ENDIF
           LARGET=3                  ! option 2
           CALL RUNMENT('WARNING: A SPATIAL CONFIGURATION IS USED')
          ELSE
           IF(LARGET.EQ.0) THEN
            LARGET=3                 ! Spatial configuration
            CALL RUNMENT('WARNING: A SPATIAL CONFIGURATION IS USED')
           ENDIF
           IF(LARGET.EQ.1) THEN
            CALL RUNENDT('ERROR: LARGET=1 NOT IMPLEMENTED')
           ENDIF
           IF(LARGET.EQ.2) THEN
            CALL RUNENDT('ERROR: LARGET=2 NOT IMPLEMENTED')
           ENDIF
           IF(LARGET.EQ.3) THEN
            CALL RUNMENT('WARNING: ONLY DENSITY AT THE MATERIAL CONFIGUR
     .ATION IS NEEDED TO BE INPUT')
           ENDIF
          ENDIF                      ! large.eq.0
         ENDIF                       ! itermd.eq.0
        ENDIF                        ! words(i).eq.'defor'
        IF(WORDS(I).EQ.'FRICT') THEN         ! frictional heating
         ITERMF=INT(PARAM(I))
         IF(ITERMF.EQ.1.AND.NOCOIT.NE.1) ! to be improved; see pro104t.f
     .    CALL RUNENDT('ERROR: FRIC. HEATING ONLY IMPLEM. WITH NC MESH')
        ENDIF
       ENDDO                         ! i=2,5
      ENDIF                          ! iterme.eq.2
C
      IF(ITERMD.EQ.1.AND.NMEMO10.EQ.1) THEN
       CALL RUNMENT
     . ('WARNING: ITERMD=1 & NMEMO10=1 INCOMPATIBLE; NMEMO10 SET TO 0')
       NMEMO10=0
      ENDIF
      IF(NMEMO10.EQ.1) THEN
       CALL RUNMEN('WARNING: DENSITY VARIATION IS NOT RIGOROUSLY 
     .  CONSISTENT WITH THE MASS CONSERVATION')
       CALL RUNMENT('WARNING: DENSITY VARIATION IS NOT RIGOROUSLY
     .  CONSISTENT WITH THE MASS CONSERVATION')
      ENDIF
      GO TO 1000
C
    3 CONTINUE
      ISTAGG=INT(PARAM(1))
      GO TO 1000
C
    4 CONTINUE
      INICOU=1
      GO TO 1000
C
    5 CONTINUE
      NPOINC=INT(PARAM(1))
      IF(NPOINC.EQ.0) CALL RUNEND('COUINC: NPOINC=0')
      IBAND1=1
      GO TO 1000
C
    6 CONTINUE
      NELEMC=INT(PARAM(1))
      IF(NELEMC.EQ.0) CALL RUNEND('COUINC: NELEMC=0')
      IBAND2=1
      GO TO 1000
C
    7 CONTINUE
      COUFAC=PARAM(1)
      CENKEL=PARAM(2)
      NITERC=INT(PARAM(3))
C
      IF(NITERC.EQ.1) THEN
       IF(ITERME.GT.0) THEN
        IF(ITERMP.GT.0) THEN
         IF(NMEMO3.EQ.0) THEN
          CALL RUNMENT('WARNING: NMEMO3 IS SET TO 1 BECAUSE NITERC=1')
          NMEMO3=1
          NHISTT=NHISTT+NHIST1                           ! see consett.f
         ENDIF
        ENDIF
       ENDIF
      ENDIF
C
      IF(NITERC.EQ.2.OR.NITERC.EQ.3)
     . CALL RUNENDT('ERROR: NITERC.EQ.2.OR.NITERC.EQ.3 not implemented')
C
      IF(ITERME.LE.0) THEN
       IF(NITERC.NE.0)
     .  CALL RUNMENT('WARNING: NITERC IS SET TO 0 WHEN ITERME LE 0')
        NITERC=0
      ELSE
       IF(ITERMP.EQ.0) THEN
        IF(NITERC.NE.0)
     .   CALL RUNMENT('WARNING: NITERC IS SET TO 0 WHEN ITERMP LE 0')
        NITERC=0
       ENDIF
      ENDIF
C
      IF (NITERC.EQ.4) THEN
       IF(LARGE.EQ.0) THEN
        CALL RUNMENT('WARNING: NITERC IS SET TO 1 WHEN LARGE = 0')
        NITERC=1
       ENDIF
C
       IF(ITERME.GT.0) THEN
        IF(ITERMP.GT.0) THEN
         IF(NMEMO3.EQ.0) THEN
          CALL RUNMENT('WARNING: NMEMO3 IS SET TO 1 BECAUSE NITERC=4')
          NMEMO3=1
         ENDIF
        ENDIF
       ENDIF
      ENDIF
C
ctm   NFCOUC=INT(PARAM(4))
      GO TO 1000
C
    8 CONTINUE                              ! only useful for ITERMP > 0
      NPLCOU=INT(PARAM(1))
      IF(NPLCOU.EQ.1) FPLCOU=PARAM(2)
      IF(ITERMP.EQ.0)
     . CALL RUNMENT('WARNING: NO PLASTIC COUPLING TERM IS CONSIDERED')
      GO TO 1000
C
    9 CONTINUE
C
      IF(IBAND1.NE.1)
     . CALL RUNEND('NPOINC IN COUPLED FILE HAS NOT BEEN FOUND')
      IF(IBAND2.NE.1)
     . CALL RUNEND('NELEMC IN COUPLED FILE HAS NOT BEEN FOUND')
C
      IF(INICOU.EQ.0.AND.KDYNA.EQ.1)
     . CALL RUNMEN('WARNING: INICOU=0 WITH KDYNA=1')
      IF(INICOU.EQ.1.AND.KDYNAT.EQ.1)
     . CALL RUNMEN('WARNING: INICOU=1 WITH KDYNAT=1')
C
      IF(NOCOI.NE.NOCOIT)
     . CALL RUNEND('ERROR: NON-COINCIDENT MESH MUST BE IN M & T PROB.')
C
C**** LOOK FOR 'STOP' CARD
C
      CALL LISTENC('COUINC',IPRINC,ITAPEC)
      IF(WORDS(1).NE.'STOP') GO TO 2001
C
C**** ESTABLISH SOME DIMENSIONS
C
      NTOTVM=NTOTV
      NDOFCM=NDOFC
C
      RETURN 
 2000 CALL RUNEND('COUINC: ERROR IN COUPLED DATA BLOCK')
 2001 CALL RUNEND('COUINC: NO STOP CARD IN COUPLED DATA BLOCK')
C
  900 FORMAT(A6,1X,8A8)
  901 FORMAT(' ERROR IN OPENING OUTPUT  FILE  36 (35 LINUX)')
  902 FORMAT(' ERROR IN OPENING PROCESS FILE 236 (36 LINUX)')
  905 FORMAT(///5X,A6,1X,8A8//)
      END
