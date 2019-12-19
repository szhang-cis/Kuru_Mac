      SUBROUTINE FIXELM(LNODS,ELDAT,MATNO,PROEL,IFFIX,NWORK)
C***********************************************************************
C
C**** THIS ROUTINE CALLS THE ELEMENTAL ROUTINES FOR DOF ACTIVATION 
C     IF THESE OPERATIVE ROUTINES DO NOT EXIST ALL D.O.F. ARE ACTIVATED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION LNODS(NNODE,*),MATNO(NELEM),IFFIX(NTOTV)
      DIMENSION NWORK(*)
      REAL*8  ELDAT(*),PROEL(NPREL,*)
C
C**** LOOP OVER ELEMENTS
C
      DO ITOTV=1,NTOTV   ! Desactivate all degrees of freedom
       IFFIX(ITOTV)=1
      ENDDO
      DO IELEM=1,NELEM
       LGRUP=MATNO(IELEM)
       LMATS=INT(PROEL(1,LGRUP))
       NNODL=INT(PROEL(2,LGRUP))
       CALL DATBAS(ELDAT,    1,    2)
C
C**** The element must acrivate d.o.f.
C
       CALL ELMLIB(LNODS(1,IELEM),PROEL(1,LGRUP),PROPS,
     .             NWORK,NWORK,NWORK,NWORK,NWORK,
     .             ELDAT,ELPRE,ELVAR,ELMAT,NWORK,   17)
       CALL SCATIN(NWORK,NDOFN,NNODL,IFFIX,NDOFN,NPOIN,LNODS(1,IELEM))
      ENDDO
C
      RETURN
      END
