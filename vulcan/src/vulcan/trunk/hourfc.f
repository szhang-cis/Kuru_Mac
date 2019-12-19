
      SUBROUTINE HOURFC(BMSIG,ELDIS,EMASS,FMULT,NEVAB)
C*****************************************************************
C
C****THIS ROUTINE ADDS THE HOURGLASS STABILIZATION FORCES
C
C*****************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION BMSIG(*), ELDIS(*), EMASS(*)
C
      DO IEVAB=1,NEVAB                             ! upper triangle
        KPOSI=NEVAB*(IEVAB-1)-(IEVAB*IEVAB-IEVAB)/2
        DO JEVAB=IEVAB,NEVAB
          KPOSN=KPOSI+JEVAB
          BMSIG(IEVAB)=BMSIG(IEVAB)+EMASS(KPOSN)*ELDIS(JEVAB)
        ENDDO
      ENDDO
      DO JEVAB=1,NEVAB-1                           ! lower triangle
        KPOSJ=NEVAB*(JEVAB-1)-(JEVAB*JEVAB-JEVAB)/2
        DO IEVAB=JEVAB+1,NEVAB
          KPOSN=KPOSJ+IEVAB
          BMSIG(IEVAB)=BMSIG(IEVAB)+EMASS(KPOSN)*ELDIS(JEVAB)
        ENDDO
      ENDDO
      DO IEVAB=1,NEVAB
        BMSIG(IEVAB)=BMSIG(IEVAB)*FMULT
      ENDDO
C
      RETURN
      END
