      SUBROUTINE FORCIN(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,ELOAD,HEADS,
     .                  LNODS,MATNO,PROEL,PROPS,WORK1,TEMPN,DTEMP,
     .                  DISPR,PWORK,PREAS,TGAPS,VNORM,COORD,TEMPI,
     .                  FPCHA,DISIT,LACTI,NOPRF,PRESF,VANIS,STDGA,
     .                  CTDGA)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE RESISTING FORCES
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
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION DISTO(*),           TEMPN(*),
     .          DTEMP(*),           ELOAD(*),
     .          HEADS(*),           DISPR(*)
      DIMENSION PWORK(NPOIN),       PREAS(*),
     .          TGAPS(*),           VNORM(*),
     .          COORD(NDIME,NPOIN), TEMPI(*),
     .          FPCHA(NFPCH,NPOIN), DISIT(*),
     .          LACTI(NELEM)
      DIMENSION NOPRF(NNODE,NELEM), PRESF(NNODE,NDIME,NELEM),
     .          VANIS(NANIV,NANIC,NELEM)
      DIMENSION STDGA(NPRE4,NPOIN), CTDGA(NPRE5,NPOIN)
C
      KUNLD=0
      KEROR=0
C
      NDOFT=1
C
C**** LOOP OVER ELEMENTS
C
      DO 1000 IELEM=1,NELEM
      IF(NACTI.EQ.1) THEN
       IF(LACTI(IELEM).EQ.0) GO TO 1000       ! skip non active elements
      ENDIF
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL(1,LGRUP))
      NNODL=INT(PROEL(2,LGRUP))
      LTYPE=INT(PROEL(5,LGRUP))
C
      NNODX=NNODL
      NNODY=NNODL
      NNODZ=NNODL
      IF(LTYPE.EQ.32) THEN
       NOCOL=INT(PROEL(22,LGRUP))
       IF(NOCOL.EQ.1) THEN                    ! non-coincident mesh
        NNODN=INT(PROEL(21,LGRUP))
        NNODX=NNODN
        NNODY=NNODN  ! with large displacements (ld) or with friction
        NNODZ=NNODN  ! with ld and linearized contact or with friction
       ENDIF
      ENDIF
C
      IEROR=0
C
C**** READ ELEM. DATA FROM DATA BASE
C
      IF(KPROB.EQ.3) THEN
       CALL DATBAS(ELMAT(IMATX(2)),    7,    2) ! ESTIF
      ELSE
       IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0)
     .  CALL DATBAS(ELDAT,    1,    2)  ! geometrical data
       IF(NMEMOM.EQ.1.OR.INITV.EQ.1)
     .  CALL DATBAS(ELPRE,    2,    2)  ! current values
       CALL DATBAS(ELVAR,    5,    2)   ! last converged
      END IF
C
      IF(NMEMO5M.EQ.0) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS(ELVAR) )
C
       CALL GATHER(DISTO,NDOFC,NPOIN,ELVAR,           NDOFC,NNODX,
     .             LNODS(1,IELEM))
      ELSE
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS )
C
       CALL GATHER(DISTO,NDOFC,NPOIN,WORK1(IFORC(37)),NDOFC,NNODX,
     .             LNODS(1,IELEM))
      ENDIF                   ! nmemo5m.eq.0
C
C**** SCALAR GATHER OPERATIONS ( DISPR ---> WORK1(IFORC(15))     )
C
      CALL GATHER(DISPR,NDOFC,NPOIN,WORK1(IFORC(15)),NDOFC,
     .            NNODZ,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( VNORM ---> WORK1(IFORC(18))     )
C
      IF(LTYPE.EQ.4) THEN
       CALL GATHER(VNORM,NDOFC,NPOIN,WORK1(IFORC(18)),NDIME,
     .             NNODL,LNODS(1,IELEM))
      ENDIF                    ! ltype.eq.4
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1M.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORD,NDIME,NPOIN,WORK1(IFORC(19)),NDIME,
     .             NNODY,LNODS(1,IELEM))
      ENDIF
C
      IF(ITERME.GE.0) THEN                ! uni or bidirectional coupled
C
C**** SCALAR GATHER OPERATIONS ( TEMPN ---> WORK1(IFORC(11)) )
C
       CALL GATHER(TEMPN,NDOFT,NPOIN,WORK1(IFORC(11)),NDOFT,
     .             NNODX,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( DTEMP ---> WORK1(IFORC(12)) )
C
       CALL GATHER(DTEMP,NDOFT,NPOIN,WORK1(IFORC(12)),NDOFT,
     .             NNODL,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( TEMPI ---> WORK1(IFORC(20)) )
C
       CALL GATHER(TEMPI,NDOFT,NPOIN,WORK1(IFORC(20)),NDOFT,
     .             NNODL,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( FPCHA ---> WORK1(IFORC(36)) )
C
       CALL GATHER(FPCHA,NFPCH,NPOIN,WORK1(IFORC(36)),NFPCH,
     .             NNODL,LNODS(1,IELEM))
      ENDIF                  ! iterme.ge.0
C
C**** SCALAR GATHER OPERATIONS ( DISIT ---> WORK1(IFORC(38)) )
C
      IIAUX=0
      IF(LTYPE.EQ.4) THEN
       IF(NPOIC.GT.0) THEN
        IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) IIAUX=1
       ENDIF
      ENDIF                  ! ltype.eq.4
      IF(LTYPE.EQ.32) THEN
       IIAUX=1               ! only needed for friction
      ENDIF
      IF(NLDSF.EQ.1) THEN
       IF(IIAUX.EQ.0) IIAUX=1
      ENDIF
      IF(IIAUX.EQ.1) THEN
       CALL GATHER(DISIT,NDOFC,NPOIN,WORK1(IFORC(38)),NDOFC,
     .             NNODY,LNODS(1,IELEM))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( STDGA,CTDGA ---> WORK1(IFORC(48,49)) )
C
      IF(LTYPE.EQ.33) THEN
       CALL GATHER(STDGA,NPRE4,NPOIN,WORK1(IFORC(48)),NPRE4,
     .             NNODL,LNODS(1,IELEM))           ! stress
       CALL GATHER(CTDGA,NPRE5,NPOIN,WORK1(IFORC(49)),NPRE5,
     .             NNODL,LNODS(1,IELEM))           ! constitutive tensor
      ENDIF                    ! ltype.eq.33
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .            WORK1,WORK1,NOPRF(1,IELEM),PRESF(1,1,IELEM),
     .            VANIS(1,1,IELEM),
     .            ELDAT,ELPRE,ELVAR,ELMAT,WORK1,    4)
C
      IF(IEROR.NE.0) GO TO 1000
C
C**** WRITE CURRENT STATE VARIABLES TO DATA BASE
C 
      CALL DATBAS(ELVAR,    3,    1)
C
C**** SCALAR SCATTER OPERATION ( WORK1 --> ELOAD )
C
      CALL SCATER(WORK1,NDOFN,NNODL,ELOAD,NDOFC,NPOIN,LNODS(1,IELEM))
C
      IF(ITERME.GT.0) THEN      ! bidirectional coupling
C
C**** SCALAR SCATTER OPERATION ( WORK1 --> PWORK )
C
       CALL SCATER(WORK1(IFORC(14)),NDOFT,NNODL,PWORK,NDOFT,NPOIN,
     .             LNODS(1,IELEM))
      ENDIF                     ! iterme.gt.0
C
      IF(KGAPC.EQ.1) THEN
C
C**** SCALAR SCATTER ASSIGN OPERATION ( WORK1 --> PREAS )
C
       IF(LTYPE.EQ.4)
     .  CALL SCATERA(WORK1(IFORC(16)),NDOFT,NNODL,PREAS,NDOFT,NPOIN,
     .               LNODS(1,IELEM))
       IF(LTYPE.EQ.32)
     .  CALL SCATER(WORK1(IFORC(16)),NDOFT,NNODL,PREAS,NDOFT,NPOIN,
     .              LNODS(1,IELEM))
C
C**** SCALAR SCATTER ASSIGN OPERATION ( WORK1 --> TGAPS )
C
       IF(LTYPE.EQ.4)
     .  CALL SCATERA(WORK1(IFORC(17)),NDOFT,NNODL,TGAPS,NDOFT,NPOIN,
     .               LNODS(1,IELEM))
       IF(LTYPE.EQ.32)
     .  CALL SCATER(WORK1(IFORC(17)),NDOFT,NNODL,TGAPS,NDOFT,NPOIN,
     .              LNODS(1,IELEM))
      ENDIF                     ! kgapc.eq.1
C
C**** SCALAR SCATTER ASSIGN OPERATION ( WORK1 --> DISIT )
C
      IF(LTYPE.EQ.4) THEN
       IF(NPOIC.GT.0) THEN
        IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) THEN
         CALL SCATERA(WORK1(IFORC(38)),NDOFC,NNODL,DISIT,NDOFC,NPOIN,
     .                LNODS(1,IELEM))
        ENDIF
       ENDIF
      ENDIF                    ! ltype.eq.4
C
C**** SCALAR SCATTER ASSIGN OPERATION ( WORK1 --> STDGA,CTDGA )
C
      IF(LTYPE.EQ.30.AND.IGALE.EQ.1) THEN
       CALL SCATERA(WORK1(IFORC(48)),NPRE4,NNODL,STDGA,NPRE4,NPOIN,
     .              LNODS(1,IELEM))                ! stress
       CALL SCATERA(WORK1(IFORC(49)),NPRE5,NNODL,CTDGA,NPRE5,NPOIN,
     .              LNODS(1,IELEM))                ! constitutive tensor
      ENDIF                    ! ltype.eq.30.and.igale.eq.1
C
 1000 CONTINUE
C
      IF(KEROR.NE.0)
     . CALL RUNEND('ERROR IN JACOBIAN MATRIX           ')
C
      RETURN
      END
