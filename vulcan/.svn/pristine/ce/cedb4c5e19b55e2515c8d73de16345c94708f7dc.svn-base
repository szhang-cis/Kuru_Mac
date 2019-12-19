      SUBROUTINE LOADPS(ELDAT,ELPRE,ELVAR,ELMAT,
     .                  LNODS,MATNO,PROEL,PROPS,RLOAD,RLOAH,COORD,
     .                  LACTI,NOPRF,PREHF,WORK1)
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
C     IENDOF=variable useful when using .loa external file. Warning: in
C     this case, no blank lines are allowed at the end of the .loa file
C     (routines listen*.f do not work!).
C     Actually is not used (March/99).
C
C     The 'PRINT' option can be used more than once in the input data.
C     Nevertheless, its use is recommended only at the end of the
C     LOAD_DATA in order to print all the heat contributions to the
C     reference array.
C
C     When ACTIVE_ELEMENTS are used, SURFACE_LOAD or GRAVITY_LOAD cards
C     must be defined after the activation definition in the input data.
C     In this case, only the loads corresponding to the active elements
C     will be considered in the analysis. Therefore, if new elements are
C     activated (in a new interval) with new loads associated to them,
C     new SURFACE_LOAD or GRAVITY_LOAD must be defined. If new elements
C     are activated but no new loads are associated to them, the use of
C     the option OLD_LOAD is correct.
C     Note that GRAVITY_LOAD must be defined for each element with,
C     in general, different associated functions (as VOLUME_HEAT in the
C     thermal problem) in order to consider it properly when active
C     elements are used (see loagra.f). In this case, the function
C     input in the card following the GRAVITY_LOAD is not used. 
C     This procedure is different from the standard one (no element
C     activation) in which the GRAVITY_LOAD is defined for all the
C     elements with the same associated function.
C
C     When the external .loa file is used, input the RETURN card in this
C     file at the end of each interval load data.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/LDFILE/ITAPE
      COMMON/LDFILE1/IAXES
      COMMON/ENDOFFILE/IENDOF
C
      DIMENSION MATNO(NELEM),         LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP),   PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),         ELPRE(NPREV),
     .          ELVAR(NSTAT),         ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION RLOAD(NTOTV),         RLOAH(NTOTV,NFUNC)
      DIMENSION COORD(NDIME,NPOIN),   LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM),   PREHF(NNODE,NDIME,NELEM,NFUNC)
C
      PARAMETER (MCOMM=5)
C
      CHARACTER*5 COMMD(MCOMM)
C
      DATA COMMD/'TITLE','POINT','GRAVI','FACE_','PRINT'/
C
      WRITE(LURES,900)
C
      IERR1=1
      IF(ITAPE.NE.LUDAT) THEN
       ITAPE=LULOA     ! external .loa file
       IF(IOLOA.EQ.0) OPEN(UNIT=LULOA,FILE=CE1,STATUS='OLD',ERR=1000)
       IOLOA=1
       IERR1=0
C
 1000  IF(IERR1.NE.0) THEN
        IF(IERR1.EQ.1) WRITE(LURES,901)
        CALL RUNEND('ERROR IN OPENING FILE ')
       ENDIF
      ENDIF      ! itapet.eq.44
C
C**** INITIALISE REFERENCE LOAD ARRAY
C
      DO ITOTV=1,NTOTV
       RLOAD(ITOTV)=0.0D0
       DO IFUNC=1,NFUNC
        RLOAH(ITOTV,IFUNC)=0.0D0
       END DO
      END DO
C
C**** INITIALISE DEFORMATION-DEPENDENT LOAD
C
      IF(NLDSF.EQ.1) THEN
       DO IELEM=1,NELEM
        DO INODE=1,NNODE
         NOPRF(INODE,IELEM)=0
         DO IDIME=1,NDIME
          DO IFUNC=1,NFUNC
           PREHF(INODE,IDIME,IELEM,IFUNC)=0.0D0
          END DO
         END DO
        END DO
       END DO
      END IF
C
      IPRINAUX=0
C
      IPOINAUX=0
      IFACEAUX=0
      IGRAVAUX=0
C
      ILDSF=0        ! global index for deformation-dependent face loads
C
C**** COMPUTES NUMBER OF POINTS WITHOUT CONTACT (AL METHOD)
C
      NPOI1=NPOIN-NPOIC
C
C**** READ NEW COMMAND
C
  200 NPRIN=0
      CALL LISTEN('LOADPS',NPRIN,ITAPE)
C
      IF(IENDOF.EQ.1) GO TO 2000                  ! end of file detected
C
C**** IDENTIFY COMMAND
C
      DO ICOMM=1,MCOMM
       IF(WORDS(1).EQ.COMMD(ICOMM)) GO TO 300
      END DO
      GO TO 2000
C
C**** EXECUTE APROPRIATE COMMAND
C
  300 CONTINUE
      GO TO (1,2,3,4,5), ICOMM
C
C**** 'TITLE' READ LOAD TITLE
C
    1 READ(ITAPE,902) SUBTI
      WRITE(LURES,905) SUBTI
      GO TO 200
C
C**** 'POINT' NODAL POINT SECTION
C
    2 CONTINUE
      IPOINAUX=IPOINAUX+1
      IF(IPOINAUX.GT.1)
     . CALL RUNEND('ERROR: POINT_LOAD DEFINED TWICE')
      IF(IPRINAUX.EQ.1)
     . CALL RUNMEN('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(PARAM(1).EQ.0)
     . CALL RUNEND('ERROR: WRONG NUMBER OF POINT LOADS')
      NLOPO=DINT(PARAM(1))
      CALL LOAPOI(NDIME,NDOFC,NDOFN,NPOIN,KPROB,NTOTV,NFUNC,
     .            RLOAD,RLOAH,NLOPO,LURES)
      GO TO 200
C
C**** 'GRAVI' GRAVITY LOAD SECTION
C
    3 CONTINUE
      IGRAVAUX=IGRAVAUX+1
      IF(IGRAVAUX.GT.1)
     . CALL RUNEND('ERROR: GRAVITY_LOAD DEFINED TWICE')
      IF(IPRINAUX.EQ.1)
     . CALL RUNMEN('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(NACTI.EQ.1) THEN
       IF(NEWAC.EQ.-1)
     .  CALL RUNEND('ERROR: VOLUME CARD MUST BE DEFINED AFTER ACTIV.')
       ISUVOAC=1
      ENDIF
      CALL LOAGRA(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .            PROPS,RLOAD,RLOAH,COORD,LACTI,WORK1)
      GO TO 200
C
C**** 'FACE_' DISTRIBUTED EDGE LOADS SECTION
C
    4 CONTINUE
      IFACEAUX=IFACEAUX+1      ! face load can be defined more than once
      IF(IPRINAUX.EQ.1)
     . CALL RUNMEN('WARNING: PRINT OF REF.ARRAY MAY BE NOT COMPLETE')
      IF(NACTI.EQ.1) THEN
       IF(NEWAC.EQ.-1)
     .  CALL RUNEND('ERROR: SURFACE CARD MUST BE DEFINED AFTER ACTIV.')
       ISUVOAC=1
      ENDIF
C
      IAXES=-1
      IF(WORDS(2).EQ.'LOCAL') IAXES=0
      IF(WORDS(2).EQ.'GLOBA') IAXES=1
      IF(IAXES.LT.0)
     . CALL RUNEND('ERROR: LOCAL OR GLOBAL LABEL IS LACKING FOR LOADS')
      IF(PARAM(1).EQ.0.0)
     . CALL RUNEND('ERROR: WRONG NUMBER OF SURFACE LOADS')
      NEDGE=DINT(PARAM(1))
C
      KLDSF=0                              ! load set index
      IF(WORDS(3).EQ.'DEFOR') KLDSF=1      ! deformation-dependent loads
      IF(KLDSF.EQ.1) ILDSF=1
      IF(ILDSF.EQ.1.AND.NLDSF.EQ.0)
     . CALL RUNEND('ERROR: NLDSF CANNOT BE 0 FOR DEF-DEP. FACE LOAD')
C
      CALL LOASUR(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .            PROPS,RLOAD,RLOAH,COORD,LACTI,NOPRF,PREHF,
     .            NEDGE,WORK1)
      GO TO 200
C
C**** 'PRINT' PRINTOUT TOTAL NODAL FORCES
C
    5 WRITE(LURES,920)
C
      IPRINAUX=1
C
      IF(NDIME.EQ.2.AND.NTYPE.NE.3) WRITE(LURES,925)
      IF(NDIME.EQ.2.AND.NTYPE.EQ.3) WRITE(LURES,930)
      IF(NDIME.EQ.3)                WRITE(LURES,935)
C
      DO IPOIN=1,NPOI1
       ITOT0=(IPOIN-1)*NDOFC+1
       ITOTF=IPOIN*NDOFC
       WRITE(LURES,940) IPOIN,(RLOAD(ITOTV),ITOTV=ITOT0,ITOTF)
      END DO
C
      GO TO 200
C
 2000 CONTINUE
C
C**** CONTROL
C
      IF(NLDSF.EQ.1) THEN
       IF(IFACEAUX.EQ.0) THEN
        CALL RUNMEN('WARNING: NO DEF-DEP. FACE LOAD FOR THIS INTERVAL')
       ELSE
        IF(ILDSF.EQ.0)
     .   CALL RUNMEN('WARNING: NO DEF-DEP. FACE LOAD FOR THIS INTERVAL')
       ENDIF
      ENDIF
C
      RETURN
C
  900 FORMAT(1H1,///,10X,'ECHO OF LOAD & FIXITY DATA FOLLOWS :',/
     .               10X,'=================================='/)
  901 FORMAT('ERROR IN OPENING LOAD INPUT FILE')
  902 FORMAT(8A8)  
  905 FORMAT(/,5X,8A8,/)
c 920 FORMAT(/,5X,' REFERENCE ARRAY ( FORCES & PRESCRIBED DISP. )',/)
  920 FORMAT(/,5X,' REFERENCE ARRAY ( FORCES )',/)
  925 FORMAT(5H NODE,6X,9H X-FORCES,3X,9H Y-FORCES)
  930 FORMAT(5H NODE,6X,9H R-FORCES,3X,9H Z-FORCES)
  935 FORMAT(5H NODE,6X,9H X-FORCES,3X,9H Y-FORCES,3X,9H Z-FORCES)
  940 FORMAT(I7,5X,10E12.4)
      END
