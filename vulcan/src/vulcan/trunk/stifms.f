      SUBROUTINE STIFMS(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,PROPS,
     .                  WORK1,TEMPN,DISPR,DTEMP,VNORM,COORD,DISTO,INFRI,
     .                  COFRI,LACTI,CSTIF)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE STIFFNESS & COUPLING MATRICES
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
      COMMON/JACOBSA/IEROR,KEROR
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          WORK1(*)
      DIMENSION TEMPN(*),           DISPR(*),
     .          DTEMP(*),           VNORM(*)
      DIMENSION DISTO(*),           COORD(NDIME,NPOIN)
      DIMENSION CSTIF(NEVAB,*)
      DIMENSION INFRI(*),           COFRI(*),
     .          LACTI(NELEM)
C
      CALL CPUTIM(TIME1)
C
C**** READS CONTACT CONTROL PARAMETERS (external file)
C
      IF(IITER.EQ.1) THEN
       IF(ICONC.EQ.1) THEN
        CALL INPCON
       ENDIF
      ENDIF
C
C**** LOOP ON ELEMENTS
C
      IF(NACTI.EQ.1) THEN
       IF(LACTI(IELEM).EQ.0) THEN       ! deals with non active elements
        N=NEVAB*NEVAB
        CALL VECZER(N,WORK1(ISTIF(41)))             ! C
        GO TO 1001
       ENDIF
      ENDIF
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL(1,LGRUP))
      NNODL=INT(PROEL(2,LGRUP))
      LTYPE=INT(PROEL(5,LGRUP))
C
      NNODX=NNODL
      NNODY=NNODL
      IF(NOCOI.GT.0) THEN               ! non-coincident mesh
       IF(LTYPE.EQ.32) THEN
        NOCOL=INT(PROEL(22,LGRUP))
        IF(NOCOL.EQ.1) THEN
         NNODN=INT(PROEL(21,LGRUP))
         NNODX=NNODN
         IF(LARGC.NE.0) NNODY=NNODN
        ENDIF
       ENDIF
      ENDIF
C
      IEROR=0                            ! KEROR not used (see stifmx.f)
C
C**** READ ELDAT & ELVAR FROM DATA BASE
C
      IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0)
     . CALL DATBAS(ELDAT,    1,    2)
      CALL DATBAS(ELVAR,    3,    2)      ! current
C
      IF(NMEMO5M.EQ.1) THEN
C
C**** SCALAR GATHER OPERATIONS ( DISTO ---> ELDIS )
C
       CALL GATHER(DISTO,NDOFC,NPOIN,WORK1(ISTIF(27)),NDOFC,NNODX,
     .             LNODS(1,IELEM))
      ENDIF
C
      IF(ITERME.GE.0) THEN                ! uni or bidirectional coupled
C
C**** SCALAR GATHER OPERATIONS ( TEMPN ---> WORK1(ISTIF(5))     )
C
       NDOFT=1
       CALL GATHER(TEMPN,NDOFT,NPOIN,WORK1(ISTIF(5)),NDOFT,
     .             NNODX,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( DTEMP ---> WORK1(ISTIF(6))     )
C
       CALL GATHER(DTEMP,NDOFT,NPOIN,WORK1(ISTIF(6)),NDOFT,
     .             NNODL,LNODS(1,IELEM))
      ENDIF
C
C**** SCALAR GATHER OPERATIONS ( DISPR ---> WORK1(ISTIF(7))     )
C
      CALL GATHER(DISPR,NDOFC,NPOIN,WORK1(ISTIF(7)),NDOFC,
     .            NNODX,LNODS(1,IELEM))
C
C**** SCALAR GATHER OPERATIONS ( VNORM ---> WORK1(ISTIF(8))     )
C
      IF(LTYPE.EQ.4) THEN
       CALL GATHER(VNORM,NDOFC,NPOIN,WORK1(ISTIF(8)),NDIME,
     .             NNODL,LNODS(1,IELEM))
      ENDIF
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1M.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORD,NDIME,NPOIN,WORK1(ISTIF(9)),NDIME,
     .             NNODY,LNODS(1,IELEM))
      ENDIF
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .            INFRI,COFRI,WORK1,WORK1,WORK1,
     .            ELDAT,ELPRE,ELVAR,ELMAT,WORK1,    3)
C
C**** ASSIGNS OPERATION ( WORK1(ISTIF(28)) --> CSTIF )
C
 1001 IKOVA=0
      DO IEVAB=1,NEVAB
       DO JEVAB=1,NEVAB
        IKOVA=IKOVA+1
        CSTIF(IEVAB,JEVAB)=WORK1(ISTIF(28)+IKOVA-1)
       ENDDO
      ENDDO
C
      IF(IEROR.NE.0)
     . CALL RUNEND('ERROR IN JACOBIAN MATRIX           ')
C
      CALL CPUTIM(TIME2)
      CPUSF=CPUSF+(TIME2-TIME1)
C
      RETURN
      END
