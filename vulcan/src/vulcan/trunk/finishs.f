      SUBROUTINE FINISHS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,
     .                   IFFIXS,REFORS,RLOADS,TLOADS,PWORKS,TEMPIS)
C***********************************************************************
C
C**** THIS ROUTINE SAVES THE CONVERGED VALUES AND CHECKS FOR ANALYSIS
C     TERMINATION
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION  DISPRS(NTOTVS,3),    DISTOS(NTOTVS,3),
     .           ELPRES(NPREVS),      ELVARS(NSTATS),          
     .           HEADSS(NPOINS,4),    HTLODS(NHLODS,NSUBFS,NFUNCS),    
     .           IFFIXS(NTOTVS,2),    REFORS(NTOTVS),          
     .           RLOADS(NTOTVS),      TLOADS(NTOTVS,2),
     .           PWORKS(NPOINS,3),    TEMPIS(NTOTVS,2)
C
C**** CHECK IF CPU TIME IS RUNNING OUT
C
      ILATES=0
      CALL CPUTIMS(CUCPUS)
      IF(CUCPUS.GE.0.95*TLIMTS) ILATES=1
C
C**** COMPUTE SOME CONVERGED VALUES 
C
      CALL STIPARS(DISPRS,IFFIXS,TLOADS,REFORS,PWORKS)
C
C**** SAVE CONVERGED VALUES AND WRITE TO RESTART FILE IF NECESSARY
C
      CALL STORESS(DISPRS,DISTOS,ELPRES,ELVARS,HEADSS,HTLODS,IFFIXS,
     .             REFORS,RLOADS,TLOADS,TEMPIS,ILATES)
C
C**** CHECK FOR ANALYSIS TERMINATION
C
      IF(CUCPUS.GE.TLIMTS)
     .     CALL RUNENDS(' CPU TIME EXPIRED                  ')
      IF(STIFIS.NE.0.0.AND.STICUS.LT.STIFIS) 
     .     CALL RUNENDS(' FINAL STIFFNESS REACHED           ')
C
      RETURN
      END