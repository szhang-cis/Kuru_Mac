      SUBROUTINE INCSTN(ELDAT,ELPRE,ELVAR,ELMAT,HTLOD,LNODS,MATNO,
     .                  PROEL,PROPS)
C***********************************************************************
C
C****THIS ROUTINE INCREMENTS THE NON-TENSIONAL STRAINS
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
      COMMON/INCSTNA/IFUNC
      COMMON/INCSTNB/HFUNC(5)
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX)
      DIMENSION HTLOD(NHLOD,*)
C
      DO 100 IELEM=1,NELEM
C
      LGRUP=MATNO(IELEM)
      LMATS=INT(PROEL( 1,LGRUP))
      NGAUL=INT(PROEL( 4,LGRUP))
      IFUNC=INT(PROPS(37,LMATS))
C
ctm      IF((INITI.EQ.0).AND.       ! no initial data
ctm     .   (IFUNC.EQ.0)     )      ! no prescribed strains
      IF((INITI.EQ.0).OR.        ! no initial data
     .   (IFUNC.EQ.0)     )      ! no prescribed strains
     .  GO TO 100
C
      IF(NMEMOM.EQ.0)
     . CALL RUNEND('ERROR (INCSTN): MEMORY IS NOT DIMENSIONED FOR 
     .              NON TENSIONAL STRAINS')
C
C***READ ELPRE FROM DATA BASE ( AND ELVAR, IF NECESSARY ) 
C
                     CALL DATBAS(ELPRE,    4,    2) ! last converged
      IF(IFUNC.LT.0) CALL DATBAS(ELVAR,    5,    2) ! last converged
C
C***SELECT APPROPRIATE CURVE
C
      IF(IFUNC.NE.0) THEN
        DO INDEX=1,NHLOD
          HFUNC(INDEX)=HTLOD(INDEX,IABS(IFUNC))
        ENDDO
      ENDIF
C
C***CALL ELEMENT PROCESSOR TO PERFORM REQUIRED OPERATIONS
C
      CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS(1,LMATS),
     .            WORK1,WORK1,WORK1,WORK1,WORK1,
     .            ELDAT,ELPRE,ELVAR,ELMAT,WORK1,   10)
C
C***WRITE ELPRE TO DATA BASE 
C 
      CALL DATBAS(ELPRE,    2,    1)  ! current values
C
  100 CONTINUE
C
      RETURN
      END
