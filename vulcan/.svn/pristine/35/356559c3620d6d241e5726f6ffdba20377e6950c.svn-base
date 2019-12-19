      SUBROUTINE LOADPST(ELDATT,ELPRET,ELVART,ELMATT,
     .                   LNODST,MATNOT,PROELT,PROPST,RLOADT,RLOAHT,
     .                   TEMPIT,
     .                   WORK1T,
     .                   ADVELT,COORDT,FPCHAT,DISPLT,LACTIT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE CONSISTENT NODAL FORCES FOR EACH
C     ELEMENT
C
C
C     Notes:
C
C     DINT is used to avoid problems in challenge machines
C     
C     IPRINAUX=variable helping in the printing of the reference heat
C              array
C
C     IENDOFT=variable useful when using .loa external file. Warning: in
C     this case, no blank lines are allowed at the end of the .loa file
C     (routines listent*.f do not work!).
C
C     The 'PRINT' option can be used more than once in the input data.
C     Nevertheless, its use is recommended only at the end of the
C     HEAT_DATA in order to print all the heat contributions to the
C     reference array.
C     Actually is not used (March/99).
C
C     When ACTIVE_ELEMENTS are used, SURFACE_HEAT or VOLUME_HEAT cards
C     must be defined after the activation definition in the input data.
C     In this case, only the heats corresponding to the active elements
C     will be considered in the analysis. Therefore, if new elements are
C     activated (in a new interval) with new heats associated to them, 
C     new SURFACE_HEAT or VOLUME_HEAT must be defined. If new elements
C     are activated but no new heats are associated to them, the use of
C     the option OLD_HEAT is correct.
C
C     When the external .loa file is used, input the RETURN card in this
C     file at the end of each interval load data.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      COMMON/LDFILET/ITAPET
      COMMON/ENDOFFILET/IENDOFT
C
      DIMENSION MATNOT(NELEMT),        LNODST(NNODET,NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT),
     .          WORK1T(*)
      DIMENSION RLOADT(*),             RLOAHT(NTOTVT,NFUNCT),
     .          TEMPIT(NPOINT,2),      ADVELT(NTOTVT*NDIMET),
     .          COORDT(NDIMET,NPOINT), FPCHAT(NFPCH,NPOINT),
     .          DISPLT(NTOTVM),        LACTIT(NELEMT)
C
      PARAMETER (MCOMMT=5)
C
      CHARACTER*5 COMMDT(MCOMMT)
C
      DATA COMMDT/'TITLE','POINT','VOLUM','SURFA','PRINT'/
C
      WRITE(LUREST,900)
C
      IERR1T=1
      IF(ITAPET.NE.LUDATT.AND.ITAPET.NE.LULOAT) GOTO 1000
      IF(ITAPET.EQ.LULOAT) THEN       ! external .loa file
       IF(IOLOAT.EQ.0) OPEN(UNIT=LULOAT,FILE=CE1T,STATUS='OLD',ERR=1000)
       IOLOAT=1
       IERR1T=0
C
 1000  IF(IERR1T.NE.0) THEN
        IF(IERR1T.EQ.1) WRITE(LUREST,901)
        CALL RUNENDT('ERROR IN OPENING FILE ')
       ENDIF
      ENDIF      ! itapet.eq.luloat
C
C**** INITIALISE REFERENCE LOAD ARRAY
C
      DO ITOTVT=1,NTOTVT
       RLOADT(ITOTVT)=0.0
       DO IFUNCT=1,NFUNCT
        RLOAHT(ITOTVT,IFUNCT)=0.0
       ENDDO
      END DO
C
      IPRINAUX=0
C
      IPOINAUX=0
      IFACEAUX=0
      IGRAVAUX=0
C
C**** READ NEW COMMAND
C
  200 NPRINT=0
      CALL LISTENT('LOADPST',NPRINT,ITAPET)
C
      IF(IENDOFT.EQ.1) GO TO 2000                 ! end of file detected
C
C**** IDENTIFY COMMAND
C
      DO ICOMMT=1,MCOMMT
       IF(WORDST(1).EQ.COMMDT(ICOMMT)) GO TO 300
      END DO
      GO TO 2000
C
C**** EXECUTE APROPRIATE COMMAND
C
  300 CONTINUE
      GO TO (1,2,3,4,5), ICOMMT
C
C**** 'TITLE' READ LOAD TITLE
C
    1 READ(ITAPET,902) SUBTIT
      WRITE(LUREST,905) SUBTIT
      GO TO 200
C
C**** 'POINT' NODAL POINT SECTION
C
    2 CONTINUE
      IPOINAUX=IPOINAUX+1
      IF(IPOINAUX.GT.1)
     . CALL RUNENDT('ERROR: POINT_HEAT DEFINED TWICE')
      IF(IPRINAUX.EQ.1)
     . CALL RUNMENT('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(PARAMT(1).EQ.0.0)
     . CALL RUNENDT('ERROR: WRONG NUMBER OF POINT HEATS')
      NLOPOT=DINT(PARAMT(1))
      CALL LOAPOIT(NDIMET,NDOFCT,NDOFNT,NPOINT,KPROBT,NTOTVT,NFUNCT,
     .             RLOADT,RLOAHT,NLOPOT,LUREST)
      GO TO 200
C
C**** 'VOLUM' SOURCE HEAT SECTION
C
    3 CONTINUE
      IGRAVAUX=IGRAVAUX+1
      IF(IGRAVAUX.GT.1)
     . CALL RUNENDT('ERROR: VOLUME_LOAD DEFINED TWICE')
      IF(IPRINAUX.EQ.1)
     . CALL RUNMENT('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(NACTIT.EQ.1) THEN
       IF(NEWACT.EQ.-1)
     .  CALL RUNENDT('ERROR: VOLUME CARD MUST BE DEFINED AFTER ACTIV.')
       ISUVOACT=1
      ENDIF
      IF(PARAMT(1).EQ.0.0)
     . CALL RUNENDT('ERROR: WRONG NUMBER OF VOLUME HEATS')
      NLOEL1=DINT(PARAMT(1))
      CALL LOAGRAT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .             PROPST,RLOADT,RLOAHT,TEMPIT(1,1),3,WORK1T,NLOEL1,
     .             ADVELT,COORDT,TEMPIT(1,2),FPCHAT,DISPLT,LACTIT)
      GO TO 200
C
C**** 'SURFA' DISTRIBUTED EDGE HEAT SECTION
C
    4 CONTINUE
      IFACEAUX=IFACEAUX+1
      IF(IFACEAUX.GT.1)
     . CALL RUNENDT('ERROR: SURFACE_LOAD DEFINED TWICE')
      IF(IPRINAUX.EQ.1)
     . CALL RUNMENT('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(NACTIT.EQ.1) THEN
       IF(NEWACT.EQ.-1)
     .  CALL RUNENDT('ERROR: SURFACE CARD MUST BE DEFINED AFTER ACTIV.')
       ISUVOACT=1
      ENDIF
      IF(PARAMT(1).EQ.0.0)
     . CALL RUNENDT('ERROR: WRONG NUMBER OF SURFACE HEATS')
      NLOEL2=DINT(PARAMT(1))
      CALL LOAGRAT(ELDATT,ELPRET,ELVART,ELMATT,LNODST,MATNOT,PROELT,
     .             PROPST,RLOADT,RLOAHT,TEMPIT(1,1),4,WORK1T,NLOEL2,
     .             ADVELT,COORDT,TEMPIT(1,2),FPCHAT,DISPLT,LACTIT)
      GO TO 200
C
C**** 'PRINT' PRINTOUT TOTAL NODAL FORCES
C
    5 WRITE(LUREST,920)
C
      IPRINAUX=1
C
      WRITE(LUREST,925)
      DO IPOINT=1,NPOINT
       ITOT0T=(IPOINT-1)*NDOFCT+1
       ITOTFT=IPOINT*NDOFCT
       WRITE(LUREST,940) IPOINT,(RLOADT(ITOTVT),ITOTVT=ITOT0T,ITOTFT)
      END DO
C
      GO TO 200
C
 2000 CONTINUE
      RETURN
C
  900 FORMAT(1H1,///,10X,'ECHO OF HEAT & FIXITY DATA FOLLOWS :',/
     .               10X,'=================================='/)
  901 FORMAT('ERROR IN OPENING HEAT INPUT FILE')
  902 FORMAT(8A8)  
  905 FORMAT(//5X,8A8)
c 906 FORMAT(//,4X,'HEAT INCREMENTATION DATA',
c    .          '  ( NO. OF FUNCTIONS :',I3,')'//,
c    .          4X,'NUMBER',6X,'TYPE',35X,'PARAMETERS')
c 907 FORMAT(I5,8E15.6)
c 908 FORMAT(5X,I5,5X,I5,12(5X,E15.6))
c 920 FORMAT(5X,' REFERENCE ARRAY ( HEATS & PRESCRIBED TEMP. )',/)
  920 FORMAT(5X,' REFERENCE ARRAY ( HEAT FLUXES )',/)
  925 FORMAT(5H NODE,6X,5H HEAT)
  940 FORMAT(1X,I5,5X,10E12.4)
      END
