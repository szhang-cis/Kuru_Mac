      SUBROUTINE OUTPRI(DISTO,ELDAT,ELPRE,ELVAR,ELMAT,HEADS,
     .                  LNODS,MATNO,PROEL,PROPS,SMSTS,SMSTN,SMSTP,
     .                  IFFIX,REFOR,COORD,TGAPS,PREAS,DISIT,INFRI,ACCPN,
     .                  WORK1)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS DISPLACEMENTS AND STRESSES
C     TO OUTPUT FILE
C
C      KPRI0= 0  DO NOT WRITE TO OUTPUT FILE
C             1  WRITE TO OUTPUT FILES
C      KPRI1= 0  DO NOT WRITE DISPLACEMENTS TO OUTPUT FILE
C             1  WRITE DISPLACEMENTS TO OUTPUT FILES
C      KPRI2= 0  DO NOT WRITE GAUSSIAN STRESSES TO OUTPUT FILE
C             1  WRITE GAUSSIAN STRESSES TO OUTPUT FILES
C      KPRI3= 0  DO NOT WRITE GAUSSIAN PRINCIPAL STRESSES TO OUTPUT FILE
C             1  WRITE GAUSSIAN PRINCIPAL STRESSES TO OUTPUT FILES
C      KPRI4= 0  DO NOT WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE GAUSSIAN INTERNAL VARIABLES TO OUTPUT FILES
C      KPRI5= 0  DO NOT WRITE NODAL STRESSES
C             1  WRITE NODAL STRESSES
C      KPRI6= 0  DO NOT WRITE NODAL PRINCIPAL STRESSES
C             1  WRITE NODAL PRINCIPAL STRESSES
C      KPRI7= 0  DO NOT WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILE
C             1  WRITE NODAL INTERNAL VARIABLES TO OUTPUT FILES
C      KPRI8= 0  DO NOT WRITE REACTIONS TO OUTPUT FILE
C             1  WRITE REACTIONS TO OUTPUT FILES
C      KPRI9 =0  DO NOT WRITE NODAL STRAINS
C             1  WRITE NODAL STRAINS
C      KPRI10=0  DO NOT WRITE NODAL PRINCIPAL STRAINS
C             1  WRITE NODAL PRINCIPAL STRAINS
C      KPRI11=0  DO NOT WRITE NORMAL GAP
C             1  WRITE NORMAL GAP
C      KFEMV= 0  DO NOT WRITE FOR POSTPROCESOR
C             1  WRITE FOR POSTPROCESOR
C      KPRIN= 1  WRITE TO OUTPUT FILES DISPLACEMENTS AND STRESSES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION MATNO(NELEM),       LNODS(NNODE,NELEM),
     .          PROEL(NPREL,NGRUP), PROPS(NPROP,NMATS),
     .          ELDAT(NDATA),       ELPRE(NPREV),
     .          ELVAR(NSTAT),       ELMAT(NMATX),
     .          IFFIX(*),           REFOR(NTOTV)
      DIMENSION DISTO(*),           HEADS(NPOIN,*),
     .          SMSTS(NSTR1,NPOIN), SMSTN(NSTR1,NPOIN),
     .          SMSTP(NNUIN,NPOIN)
      DIMENSION STRSP(3)
      DIMENSION COORD(NDIME,NPOIN), TGAPS(NPOIN),
     .          DISIT(*),           INFRI(NPOIN),
     .          ACCPN(NPOIN),       PREAS(NPOIN)
      DIMENSION WORK1(*)
C
      IF(KFEMV.NE.0)
     . OPEN(UNIT=LUPOS,FILE=CJ,STATUS='OLD',FORM='UNFORMATTED',
     .      ACCESS='APPEND')
C
C**** OUTPUT DISPLACEMENTS
C
      CALL OUTDIS(DISTO,HEADS,DISIT)
C
C**** OUTPUT REACTIONS
C
      IF((KPRI8+KFEMV).NE.0)
     . CALL OUTREA(IFFIX,REFOR,INFRI)
C
C**** OUTPUT NORMAL GAPS
C
      IF((KPRI11+KFEMV).NE.0)
     . CALL OUTGAP(MATNO,PROEL,TGAPS,PREAS,LNODS,ACCPN)
C
C**** OUTPUT GAUSSIAN VARIABLES
C
      IF((KPRI2+KPRI3+KPRI4+KFEMV).NE.0)
     . CALL OUTGAU(ELDAT,ELPRE,ELVAR,ELMAT,
     .             LNODS,MATNO,PROEL,PROPS,COORD,DISTO,WORK1)
C
C**** OUTPUT NODAL STRESSES, STRAINS & INTERNAL VARIABLES
C
      IF(((KPRI5+KPRI6+KPRI7+KPRI9+KPRI10+KFEMV).GT.0).AND.(KSGAU.NE.0))
     . CALL OUTNOD(SMSTS,SMSTN,SMSTP)
C
C**** CLOSE POSTPROCESS FILE
C
C     Warning: in CONVEX C-120, comment "close(lupos)" (february 1993)
C
      IF(KFEMV.NE.0) CLOSE(LUPOS)
C
      RETURN
      END