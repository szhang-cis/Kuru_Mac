      SUBROUTINE FORCDY(VELCT,ACCEL,LNODS,MATNO,PROEL,PROPS,ELOAD,
     .                  ELDAT,ELPRE,ELVAR,ELMAT,COORD,LACTI,WORK1)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE INERTIAL & DAMPING FORCES
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
C
      DIMENSION ACCEL(*),           VELCT(*),
     .          ELOAD(*),           LNODS(NNODE,*),
     .          MATNO(*),           PROEL(NPREL,*),
     .          PROPS(NPROP,*),     WORK1(*)
      DIMENSION COORD(NDIME,NPOIN)
      DIMENSION ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          LACTI(NELEM)
C
c     IF((ITIME*ISTEP.EQ.1).AND.(IITER.EQ.0)) RETURN     ! to be revised
C
C**** LOOP ON ELEMENTS
C
      DO 1000 IELEM=1,NELEM
      IF(NACTI.EQ.1) THEN
       IF(LACTI(IELEM).EQ.0) GO TO 1000       ! skip non active elements
      ENDIF
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL(1,LGRUP))
      NNODL=INT(PROEL(2,LGRUP))
      ITYPE=INT(PROEL(5,LGRUP))
C
C**** READ ELEM. DATA FROM DATA BASE
C
      IF(NMEMO1M.EQ.0.OR.NMEMO2M.EQ.0)
     . CALL DATBAS(ELDAT,    1,    2)  ! geometrical data
C
C**** GATHER NODAL COORDINATES
C
      IF(NMEMO1M.EQ.1) THEN    ! coordinates in a global array
       CALL GATHER(COORD,NDIME,NPOIN,WORK1(IFORD(20)),NDIME,NNODL,
     .             LNODS(1,IELEM))
      ENDIF
C
C**** GATHER NODAL ACCELERATIONS
C
      CALL GATHER(ACCEL,NDOFC,NPOIN,WORK1(IFORD(2)),NDOFN,NNODL,
     .            LNODS(1,IELEM))
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .            WORK1,WORK1,WORK1,WORK1,WORK1,
     .            ELDAT,ELPRE,ELVAR,ELMAT,WORK1,    6)
C
C**** SCATTER RESULTS INTO ELOAD
C
      CALL SCATER(WORK1(IFORD(1)),NDOFN,NNODL,ELOAD,NDOFC,NPOIN,
     .            LNODS(1,IELEM))
C
 1000 CONTINUE
C
      RETURN
      END    
