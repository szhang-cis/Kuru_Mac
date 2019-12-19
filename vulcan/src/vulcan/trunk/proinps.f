      SUBROUTINE PROINPS
C***********************************************************************
C
C**** THIS ROUTINE REDEFINES THE DEFAULT PROBLEM PARAMETERS 
C
C     PARAMETERS         CONTROLLING WORDS
C
C     KDYNA              'DYNAM','TRANS'
C     KPORE              'COUPL'
C     KPROB              'SEEPA', 'STRUC', 'WAVES', 'THERM'
C     LARGE              'LARGE'
C     others             'DIMEN'
C     COMMUNICATION      'COMMU'
C     ICONVT             'ADVEC'
C     IGALES             'WEAK_'
C     NPOROT             'POROS'
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      PARAMETER (MCOMMS=14)
      CHARACTER*5 COMMDS(MCOMMS)
C
      DATA COMMDS/'SEEPA','STRUC','WAVES','THERM','TRANS','DYNAM',
     .            'COUPL','LARGE','DIMEN','COMMU','ADVEC','WEAK_',
     .            'POROS','END_P'/
C
C**** READ PROBLEM PARAMETERS
C
      NPRINT=0
      ITAPET=LUDATS
      CALL LISTENS('PROINPS',NPRINT,ITAPET)
      IF(WORDSS(1).NE.'PROBL') GO TO 2000
 1000 CALL LISTENS('PROINPS',NPRINT,ITAPET)
C
C**** IDENTIFY COMMAND
C
      DO ICOMMS=1,MCOMMS
       IF(WORDSS(1).EQ.COMMDS(ICOMMS)) GO TO 100
      END DO
      GO TO 2000 
C
C**** EXECUTE APROPRIATE COMMAND
C
  100 CONTINUE
      GO TO (1,1,1,1,2,2,3,4,5,6,7,8,9,10), ICOMMS
C
    1 IF(WORDSS(1).EQ.'SEEPA')                          KPROBS=0
      IF(WORDSS(1).EQ.'STRUC')                          KPROBS=1
      IF(WORDSS(1).EQ.'STRUC'.AND.WORDSS(2).EQ.'SHELL') KPROBS=2
      IF(WORDSS(1).EQ.'WAVES')                          KPROBS=3
      IF(WORDSS(1).EQ.'THERM'.OR.WORDSS(1).EQ.'SPECI')  KPROBS=4
      GO TO 1000
C
    2 KDYNAS=1
      GO TO 1000
C
    3 KPORES=2
      IF(WORDSS(2).EQ.'DISCR') KSMUSS=-1
      IF(WORDSS(2).EQ.'LOCAL') KSMUSS=-2
      GO TO 1000
C
C**** LARGE DISPLACEMENTS/STRAINS ANALYSIS
C
C     Notes:
C
C     LARGES=1: Total Lagrangian Formulation
C     LARGES=2: Updated Lagrangian Formulation
C     LARGES=3: Eulerian Formulation
C
    4 LARGES=1                                   ! default
C
      IF(ITERME.LT.0) THEN                  ! uncoupled thermal analysis
       LARGES=0
       CALL RUNMENS('WARNING: LARGES NE 0 IS NOT CONSIDERED')
      ELSE
       IF(WORDSS(2).NE.'     ') THEN
        IF(WORDSS(2).NE.'TOTAL'.AND.WORDSS(2).NE.'UPDAT'.AND.
     .     WORDSS(2).NE.'EULER')
     .   CALL RUNENDS('PROINPT: WRONG LARGE STRAINS INPUT DATA')
c       IF(WORDSS(2).EQ.'TOTAL') LARGES=1
        IF(WORDSS(2).EQ.'TOTAL')
     .   CALL RUNENDS('ERROR: LARGES=1 NOT IMPLEMENTED')
c       IF(WORDSS(2).EQ.'UPDAT') LARGES=2
        IF(WORDSS(3).EQ.'TOTAL')
     .   CALL RUNENDS('ERROR: LARGES=2 NOT IMPLEMENTED')
        IF(WORDSS(2).EQ.'EULER') LARGES=3
       ENDIF                     ! wordss(2).ne....
      ENDIF                      ! iterme.lt.0
      GO TO 1000
C
    5 DWORDS(1)='NPOIN'
      DWORDS(2)='NELEM'
      DWORDS(3)='NDIME'
      DWORDS(4)='NNODE'
      DWORDS(5)='NGAUS'
      DWORDS(6)='NSETS'
      DWORDS(7)='NMATS'
      DWORDS(8)='NFUNC'
      DPARAS(1)=0.
      DPARAS(2)=0.
      DPARAS(3)=2.
      DPARAS(4)=4.
      DPARAS(5)=4.
      DPARAS(6)=1.
      DPARAS(7)=1.
      DPARAS(8)=1.
C
      CALL SORTERS(1)
C
      NPOINS=INT(DPARAS(1))
      NELEMS=INT(DPARAS(2))
      NDIMES=INT(DPARAS(3))
      NNODES=INT(DPARAS(4))
      NGAUSS=INT(DPARAS(5))
      NGRUPS=INT(DPARAS(6))
      NMATSS=INT(DPARAS(7))
      NFUNCS=INT(DPARAS(8))
C
      IF(NGRUPS.LT.NMATSS)
     . CALL RUNENDS('NGRUPS CAN NOT BE LESS THAN NMATSS')
      GO TO 1000
C
    6 CONTINUE
      GO TO 1000
C
    7 ICONVS=1                          ! advective effects
c     ICUBICS=1                         ! fpc cubical (see capcoft.f) ??
c     IF(PARAMS(1).NE.0.0) ICUBICS=0
      GO TO 1000
C
    8 IF(WORDSS(2).EQ.'GALER') IGALES=0 ! default value (see setdatt.f)
C
      IF(WORDSS(2).EQ.'PETRO') THEN     ! upwinding techniques
       IGALES=1
       IF(WORDSS(3).EQ.'UPWIN') THEN
        IF(PARAMS(1).NE.0.0) IUPWIS=INT(PARAMS(1))  ! upwind function
       ELSE
        IUPWIS=1                                    ! default
       ENDIF
       IF(WORDSS(4).EQ.'PERTU') THEN
        IF(PARAMS(2).NE.0.0) IPERTS=INT(PARAMS(2)) ! perturbat. function
       ELSE
        IPERTS=1                                    ! default
       ENDIF
       ICCONS=0
       IF(WORDSS(5).EQ.'PECLE') THEN
        ICCONS=1
        IF(PARAMS(3).NE.0.0) ISUPWS=INT(PARAMS(3))  ! Peclet definition
        IF(IUPWIS.LT.1.OR.IUPWIS.GT.5)
     .   CALL RUNENDS('ERROR IN UPWIND FUNCTION CHOICE')
        IF(IUPWIS.EQ.3)
     .   CALL RUNENDS('ERROR IN UPWIND FUNCTION CHOICE')
       ENDIF
       IF(WORDSS(5).EQ.'PARAM') THEN
        ICCONS=1
        IF(PARAMS(3).NE.0.0) ISUPWS=INT(PARAMS(3)) ! paramet. definition
        IF(IUPWIS.LT.3.OR.IUPWIS.GT.3)
     .   CALL RUNENDS('ERROR IN UPWIND FUNCTION CHOICE')
       ENDIF
       IF(ICCONS.EQ.0) ISUPW=1                      ! default
       IF(WORDSS(6).EQ.'FACTO') THEN
        IF(PARAMS(4).NE.0.0) EXTUPS=PARAMS(4)       ! factor
       ELSE
        EXTUPS=1.0                                  ! default
       ENDIF
C
C**** SOME CONTROLS
C
       IF(IUPWIS.LT.1.OR.IUPWIS.GT.5)
     .  CALL RUNENDS('ERROR IN UPWIND FUNCTION CHOICE')
       IF(IPERTS.LT.1.OR.IPERTS.GT.3)
     .  CALL RUNENDS('ERROR IN PERTURBATION FUNCTION CHOICE')
       IF(ISUPWS.LT.1.OR.ISUPWS.GT.6)
     .  CALL RUNENDS('ERROR IN PECLET DEFINITION CHOICE')
       IF(EXTUPS.LT.0.0)
     .  CALL RUNENDS('ERROR IN FACTOR CHOICE')
      ENDIF
      GO TO 1000
C
    9 NPOROS=4                         ! porosity criteria
      GO TO 1000
C
   10 CONTINUE
C
C**** ADDITIONAL CONTROL
C
      IF(ICONVS.EQ.1.AND.LARGES.NE.0) THEN
       CALL RUNMENS('WARNING : ICONVS=1 & LARGES NE 0 => INCOMPATIBLE;
     . LARGES IS SET TO 0')
       LARGES=0
      ENDIF
C
      RETURN 
 2000 CALL RUNENDS('PROINPS:ERROR IN PROBLEM DATA BLOCK')
      END
