      SUBROUTINE STARTS(DISTO,DISPR,ELDAT,ELPRE,ELVAR,ELMAT,
     .                  HEADS,HTLOD,IFFIX,LNODS,LPNTN,MATNO,
     .                  PROEL,PROPS,REFOR,RLOAD,TLOAD,COORD,
     .                  TEMPN,DTEMP,PWORK,PREAS,TGAPS,VNORM,
     .                  FPCHA,LACTI,PRESF,VANIS,WORK1,INDES)
C***********************************************************************
C
C**** THIS ROUTINE INITIALISES NODAL AND GAUSSIAN VARIABLES
C
C     Notes: NFURESM=2 (future restart) & IREST=1 (restart) have to be
C            improved
C            Dimensions of PREAS are inverted in order to properly
C            transfer it to prevos.f & setmtx.f
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
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          COORD(NDIME,NPOIN), WORK1(*)
      DIMENSION DISTO(NTOTV,*),     DISPR(NTOTV,*), 
     .          HEADS(NPOIN,*),     HTLOD(NHLOD,NSUBF,*), 
     .          IFFIX(NTOTV,*),     LPNTN(*),
     .          REFOR(*),           RLOAD(*), 
     .          TLOAD(NTOTV,*)
      DIMENSION TEMPN(NPOIN,2),     DTEMP(NPOIN),
     .          PWORK(NPOIN,2),     PREAS(NPOIN,NPREA),
     .          TGAPS(NPOIN),       VNORM(NTOTV),
     .          FPCHA(NFPCH,NPOIN), LACTI(NELEM)
      DIMENSION PRESF(NNODE,NDIME,NELEM),
     .          VANIS(NANIV,NANIC,NELEM)
C
      IF(IREST.EQ.0) THEN
C***********************************************************************
C
C**** THIS IS A FRESH RUN
C
C***********************************************************************
       IF(INDES.EQ.1) THEN
C
C**** READ IN MOST OF THE PROBLEM DATA
C
        CALL INPDAT(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .              PROPS,COORD,VANIS,WORK1)
C
C**** RENUMBER NODES
C
        CALL RENUMN(KRENU,LNODS,LPNTN,NELEM,NNODE,NPOIN,WORK1)  
C
C**** WRITE TO DATA BASE
C
        IF(NFURESM.EQ.2)
     .   CALL RSTAR2(LPNTN,LNODS,MATNO,PROEL,PROPS,    1)
C
C**** INITIALIZE GLOBAL AND ELEMENTAL VARIABLES
C
        CALL ZEROTE(ELPRE,ELVAR,DISTO,HEADS,TLOAD,RLOAD,REFOR,TEMPN,
     .              DTEMP,PWORK,PREAS,TGAPS,VNORM,FPCHA,LACTI,
     .              PRESF)
C
C**** READ IN INITIAL GLOBAL AND ELEMENTAL VARIABLES
C
        IF(INITI.EQ.1)
     .   CALL PREVOS(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HEADS,LNODS,
     .               MATNO,PROEL,PROPS,PREAS(1,NPRE1+1),
     .               PREAS(1,NPRE1+NPRE2+1),WORK1)
        IF(INITI.EQ.2)
     .   CALL RSTAR5(DISPR,DISTO,ELPRE,ELVAR,LNODS,TLOAD,LACTI)
       ELSE
C
C**** CALCULATE THE CONSTANT MATRICES AND ARRAYS OF THE PROBLEM
C
        CALL SETMTX(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,PROPS,
     .              COORD,VNORM,DISTO,PREAS(1,NPRE1+1),
     .              PREAS(1,NPRE1+NPRE2+1),WORK1,    1)
C
C**** COMPUTE INITIAL STRESSES DUE TO CONSTANT PORE WATER PRESSURE
C     DISTRIBUTION
C
        IF(INITI.EQ.1.AND.KPORE.EQ.1)
     .   CALL INISTR(ELDAT,ELPRE,HEADS,LNODS,MATNO,PROEL,PROPS)
       ENDIF                  ! indes.eq.1
C        
      ELSE
C
C***********************************************************************
C
C**** THIS IS A RESTART RUN
C
C***********************************************************************
       IF(INDES.EQ.1) THEN
C
C**** READ GENERAL DATA FROM DATA BASE
C
        CALL RSTAR2(LPNTN,LNODS,MATNO,PROEL,PROPS,    2)
C
C**** READ LAST CONVERGED VALUES FROM DATA BASE
C
        CALL RSTAR3(DISPR,DISTO,ELPRE,ELVAR,HEADS,HTLOD,IFFIX,
     .              REFOR,RLOAD,TLOAD,    2)
       ENDIF                  ! indes.eq.1
C
      ENDIF
C
C**** DUMP GPCOD TO TAPE FOR POSTPROCESS IF NECESSARY
C
      IF(INDES.EQ.1) RETURN
C
      IF(NWPOS.EQ.1) THEN
       IF(NMEMO1M.EQ.0) THEN
        CALL OUTPOS(ELDAT,LNODS,MATNO,PROEL,PROPS,
     .              WORK1(ISTAR(1)),WORK1(ISTAR(2)),
     .              ELPRE,ELVAR,ELMAT,WORK1)
       ELSE
        CALL OUTPOS(ELDAT,LNODS,MATNO,PROEL,PROPS,
     .              COORD,WORK1(ISTAR(1)),
     .              ELPRE,ELVAR,ELMAT,WORK1)
       ENDIF
      ENDIF
C
      RETURN
      END
