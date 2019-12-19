      SUBROUTINE INTERV(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HTLOD,IFFIX,
     .                  PRESC,LNODS,LPNTN,MATNO,PROEL,PROPS,RLOAD,
     .                  RLOAH,FICTO,TFICT,INFRI,COFRI,COORD,HEADS,
     .                  TLOAD,REFOR,PWORK,TGAPS,PREAS,LACTI,NOPRF,
     .                  PREHF,VANIS,WORK1,
     .                  NEQNS)
C***********************************************************************
C
C**** THIS ROUTINE SETS UP THE TYPE OF ALGORITHM AND REFERENCE LOAD
C     TO BE USED IN THE CURRENT TIME INTERVAL
C     FOR THE FIRST INTERVAL, IT ALSO EVALUATES THE CONSTANT ELEMENT
C     MATRICES OF THE PROBLEM
C
C     Note: Dimensions of PREAS are inverted in order to properly
C           transfer it to output.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION MATNO(NELEM),         LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP),   PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),         ELPRE(NPREV),
     .          ELVAR(NSTAT),         ELMAT(NMATX)
      DIMENSION DISTO(*),
     .          HTLOD(NHLOD,NSUBF,*), IFFIX(NTOTV,*),
     .          LPNTN(*),             RLOAD(NTOTV)
      DIMENSION PRESC(NTOTV,2),       RLOAH(NTOTV,NFUNC),
     .          FICTO(NFUNC),         TFICT(NFUNC)
      DIMENSION INFRI(NPOIN),         COFRI(NSKEW,NDIME,*)
      DIMENSION COORD(NDIME,NPOIN)
      DIMENSION HEADS(NPOIN,4),       TLOAD(NTOTV,2),
     .          REFOR(NTOTV,2),       PWORK(NPOIN,2),
     .          TGAPS(NPOIN),         LACTI(NELEM),
     .          PREAS(NPOIN,NPREA)
      DIMENSION NOPRF(NNODE,NELEM),   PREHF(NNODE,NDIME,NELEM,NFUNC),
     .          VANIS(NANIV,NANIC,NELEM)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
      INCLUDE 'soladd.inc'
C
C**** LOOK FOR THE BEGINING OF THE INTERVAL DATA CARD  
C
      CALL INTRED(TFICT)
C
C**** OPEN RESTART FILE
C
      IF(NFURESM.EQ.2)
     . OPEN(UNIT=LURST,FILE=CK,STATUS='UNKNOWN',ACCESS='DIRECT',
     .      FORM='UNFORMATTED',RECL=LENRC)
C
C**** READ THE NEW LOAD, NEW BOUNDARY CONDITIONS AND SOLUTION STRATEGY
C     ( UNLESS THIS IS A "CONTINUE RESTART" )
C
      IF(IREST.NE.1.OR.(IREST.EQ.1.AND.ISKIP.NE.0))
     . CALL INTLOD(ELDAT,ELPRE,ELVAR,ELMAT,HTLOD,IFFIX,PRESC,
     .             LNODS,MATNO,PROEL,PROPS,RLOAD,RLOAH,INFRI,
     .             COFRI,COORD,LACTI,DISTO,NOPRF,PREHF,VANIS,
     .             WORK1)
C
C**** COMPUTE DISPLACEMENT AT FICTITIOUS TIMES
C
      IF(INITI.EQ.1.AND.KDYNA.EQ.1.AND.IPRCO.EQ.0)
     . CALL INIMDT(DISTO)
C
C**** PRINTS INITIAL VALUES (ONLY FOR UNCOUPLED PROBLEMS OR COUPLED
C     PROBLEMS WITH NO INITIAL CONDITIONS)
C
C     Note: the second DISTO should be DISIT but it is not used when
C           IPRCO=0 (DISIT is used to print contact forces for 
C           NPOIC > 0  &  IAUGM=2
C
      IF(ITERME.LT.0.OR.(ITERME.GE.0.AND.INICOU.EQ.1)) THEN
       IF(IPRCO.EQ.0)
     .  CALL OUTPUT(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HEADS,IFFIX,
     .              LNODS,MATNO,PROEL,PROPS,TLOAD,REFOR,COORD,
     .              PWORK(1,2),TGAPS,PREAS(1,1),DISTO,INFRI,LACTI,
     .              PREAS(1,NPRE1+1),PREAS(1,NPRE1+NPRE2+1),WORK1)
      ENDIF
C
C**** VIRTUAL MEMORY ALLOCATION FOR SOLVER
C                             ( when new boundary or when re-starting)
C
      IF(NEWBO.EQ.1.OR.IREST.NE.0)       
     . CALL SOLADD(IFFIX,LNODS,LPNTN,PRESC,
     .             WORK1(IFIXY( 1)),WORK1(IFIXY( 2)),WORK1(IFIXY( 3)),
     .             WORK1(IFIXY( 4)),WORK1(IFIXY( 5)),WORK1(IFIXY( 6)),
     .             WORK1(IFIXY( 7)),
     .             NEQNS,WORK1)
C
      IF(ITIME.EQ.1.OR.IREST.NE.0) THEN   ! 1st interval or when restart
C
C**** DEAL WITH VIRTUAL MEMORY ALLOCATION 
C
       CALL ADDDAT(DISTO,ELDAT,ELPRE,ELVAR)
       LTOTL=LPRIN+LWOR1+LSOLV+LDABA
       WRITE(LURES,900) LPRIN,LWOR1,LSOLV,LDABA,LTOTL
       WRITE(LUPRI,900) LPRIN,LWOR1,LSOLV,LDABA,LTOTL
C
C**** COMPUTE CONSTANT ELEMENT MATRICES
C
       CALL CONMTX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .             PROPS,WORK1)
C
      ENDIF
C
C**** CLOSE RESTART FILE
C
      IF(NFURESM.EQ.2) THEN
       IF(KSAVE.EQ.-1) THEN
        CLOSE(LURST,STATUS='DELETE')
       ELSE
        CLOSE(LURST,STATUS='KEEP')
       ENDIF
      ENDIF
C
      RETURN
  900 FORMAT(1H1,//5X,'SUMMARY OF MEMORY REQUIREMENTS :',/,
     .             5X,'==============================  ',/,
     .            15X,'REQUIRED PERMANENT MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'REQUIRED TEMPORARY MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'REQUIRED SOLVER    MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'ALLOCAT. DATA BASE MEMORY  =',I8,' WORDS (R*8)'/
     .            15X,'ALLOCAT. TOTAL     MEMORY  =',I8,' WORDS (R*8)'/)
      END
