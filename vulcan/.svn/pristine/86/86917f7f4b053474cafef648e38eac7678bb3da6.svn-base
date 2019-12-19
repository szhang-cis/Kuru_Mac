      SUBROUTINE PREGAU(ELDAT,ELPRE,ELVAR,ELMAT,LNODS,MATNO,PROEL,
     .                  PROPS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE INITIAL STATE VARIABLES 
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/PRFILE/JTAPE,IFLAG
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION SUBTA(12)
C
      READ(JTAPE,905) SUBTA
      WRITE(LURES,906) SUBTA
C
      WRITE(LURES,915)
C
C**** LOOP ON ELEMENTS
C
      DO 1000 IELEM=1,NELEM
        LGRUP=MATNO(IELEM)
        LMATS=INT(PROEL(1,LGRUP))
C
C**** READ ELPRE, ELVAR FROM RESTART-FILE
C 
        IF(IFLAG.EQ.1) THEN
          CALL DATBAS(ELPRE,    2,    2)
        ELSE IF(IFLAG.EQ.2) THEN
          CALL DATBAS(ELPRE,    2,    2)
          CALL DATBAS(ELVAR,    3,    2)
        ELSE IF(IFLAG.EQ.3) THEN
          CALL DATBAS(ELVAR,    3,    2)
        ENDIF
C
C**** CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
        CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .              WORK1,WORK1,WORK1,WORK1,WORK1,
     .              ELDAT,ELPRE,ELVAR,ELMAT,WORK1,   11)
C
C**** WRITE ELPRE,ELVAR TO RESTART-FILE
C 
        IF(IFLAG.EQ.1) THEN
          CALL DATBAS(ELPRE,    2,    1)
        ELSE IF(IFLAG.EQ.2) THEN
          CALL DATBAS(ELPRE,    2,    1)
          CALL DATBAS(ELVAR,    3,    1)
        ELSE IF(IFLAG.EQ.3) THEN
          CALL DATBAS(ELVAR,    3,    1)
        ENDIF
C
 1000 CONTINUE
C
      RETURN
  905 FORMAT(12A6)
  906 FORMAT(5X,12A6)
  915 FORMAT(//,10X,'INITIAL STATE VARIABLES AT GAUSS POINTS',/)
      END
