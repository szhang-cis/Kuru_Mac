      SUBROUTINE INCFAC(HTLOD,FICTO,TFICT)
C***********************************************************************
C
C**** THIS ROUTINE DETERMINES THE LOAD/DISPLACEMENT FACTOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
C
      DIMENSION HTLOD(NHLOD,NSUBF,*), FICTO(NFUNC),
     .          TFICT(NFUNC)
C
      DATA FACIN/0.2/
C
      IF(LAUTO.EQ.0) THEN                         ! LOAD HISTORY
       CALL FUNLOD(FICTO,HTLOD,NSUBF,NFUNC,NHLOD,TFICT,TTIME)
      ELSE                                        ! AUTOMATIC INCREMENT
       IF(ISTEP.EQ.1) FACTO=DTIME
       IF(ISTEP.EQ.2) FACTO=FACIN*FACTO 
       IF(ISTEP.GT.2) FACTO=FACTO*DSQRT(DITER/PITER)
       DTIME=FACTO
      ENDIF
C
C**** PRINT TOTAL AND INCREMENTAL LOAD FACTORS
C
      DO IFUNC=1,NFUNC
       TFICT(IFUNC)=TFICT(IFUNC)+FICTO(IFUNC)
       WRITE(LURES,900) IFUNC,TFICT(IFUNC),FICTO(IFUNC)
      ENDDO
      RETURN
C
  900 FORMAT(1X,'FUNCT.',I2,5X,'CURRENT TIME FUNCT. FACTOR =',E15.5,/
     .       14X,'INCREM. TIME FUNCT. FACTOR =',E15.5,/)
      END
