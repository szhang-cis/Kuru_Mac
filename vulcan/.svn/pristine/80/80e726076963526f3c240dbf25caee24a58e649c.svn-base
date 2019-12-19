      SUBROUTINE OUTMIC(FPCHAT)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS (OUTPUT & POSTPROCESSOR FILES) MACROSCOPICAL
C     PHASE-CHANGES
C
C***********************************************************************
C
C     KPRI6=0  DO NOT WRITE NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C           1  WRITE NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C     KFEMV=0  DO NOT WRITE TO POSTPROCESSOR FILE
C           1  WRITE TO POSTPROCESSOR FILE
C
C***********************************************************************
C
C             IPCFO=0 TOTAL PHASE-CHANGE FORMULATION
C                     f_pc(T) is bijective
C
C                     IPCMO=1 LINEAR
C                     IPCMO=2 IDEM 1 WITH f_s
C                     IPCMO=3 PARABOLIC
C                     IPCMO=4 CUBIC
C                     IPCMO=5 SCHEIL'S EQUATION
C                     IPCMO=6 LEVER'S EQUATION
C                     IPCMO=7 ....
C
C***********************************************************************
C
C     Index of variables:
C
C     NNUPT: total number of phase-changes (maximum)
C     NNUPC: number of macroscopic phase-changes (maximum)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
      INCLUDE 'nued_om.f'   ! thermal-microstructural
      INCLUDE 'nuef_om.f'   ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION FPCHAT(NFPCH,NPOINT)
C
C**** OUTPUT NODAL MACROSCOPIC PHASE-CHANGE FUNCTION
C
      IF(IPRCOT.GT.0) THEN
       IF(KPRI6T.GT.0) THEN
        IF(NNUPC.GT.0) THEN
         WRITE(LUREST,956)
C
         DO INUPC=1,NNUPC
          WRITE(LUREST,960) INUPC
          DO IPOINT=1,NPOINT
           WRITE(LUREST,905) IPOINT,FPCHAT(INUPC,IPOINT)
          ENDDO
         ENDDO
        ENDIF      ! nnupc.gt.0
       ENDIF       ! kpri6t.gt.0
      ENDIF        ! iprcot.gt.0
C
C**** WRITE TO POSTPROCESSOR FILE MACROSCOPIC PHASE-CHANGE FUNCTION &
C     MATERIAL INDEX (PSEUDOCONCENTRATION) FOR FILLING PROBLEMS WITH
C     THERMAL ANALYSIS (see outpost.f)
C
      IPSEU=0
      IF(ITERMEF.EQ.0) THEN
       IF(NFILL.EQ.1) THEN
        IF(IMICR.EQ.0) THEN
         IPSEU=2*NNUPT+1
        ELSE
         IPSEU=2*NNUPT+NNUPO+1
        ENDIF
       ENDIF
      ENDIF
C
      IF(NNUPC.GT.0.AND.IPSEU.EQ.0) THEN
       IF(KFEMVT.EQ.1)
     .  WRITE(LUPOST) ((SNGL(FPCHAT(INUPC,IPOINT)),INUPC=1,NNUPC),
     .                                             IPOINT=1,NPOINT)
      ENDIF
      IF(NNUPC.EQ.0.AND.IPSEU.GT.0) THEN
       IF(KFEMVT.EQ.1)
     .  WRITE(LUPOST) ((SNGL(FPCHAT(INUPC,IPOINT)),INUPC=IPSEU,IPSEU),
     .                                             IPOINT=1,NPOINT)
      ENDIF
      IF(NNUPC.GT.0.AND.IPSEU.GT.0) THEN
       IF(KFEMVT.EQ.1)
     .  WRITE(LUPOST) (((SNGL(FPCHAT(INUPC,IPOINT)),INUPC=1,NNUPC),
     .                   SNGL(FPCHAT(INUPX,IPOINT)),INUPX=IPSEU,IPSEU),
     .                                              IPOINT=1,NPOINT)
      ENDIF
C
      RETURN
  905 FORMAT(5X,I5,8E15.6)
  956 FORMAT(/,5X,38H PHASE-CHANGE FUNCTION AT NODAL POINTS)
  960 FORMAT(5X,5HIPOIN,2X,22H PHASE-CHANGE FUNCTION,I5)
      END
