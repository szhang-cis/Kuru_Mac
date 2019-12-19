      SUBROUTINE SCATIN(IELVEC,NROWS,NNODE,IGLVEC,NGLOB,NPOIN,LNODS)
C***********************************************************************
C
C**** THIS ROUTINE ACTIVATE GLOBAL DOF WITH SCATER OPERATION FROM LOCAL
C     DOF
C          IELVEC(NROWS,NNODE) ---> IGLVEC(NROWS,NPOIN) 
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IELVEC(NROWS*NNODE), IGLVEC(NGLOB,NPOIN), LNODS(NNODE)
C
      DO K=1,NROWS*NNODE
       IF(IELVEC(K).EQ.0) THEN  ! Activate global d.o.f.
        INODE=(K-1)/NROWS+1
        IROWS=K-(INODE-1)*NROWS
        LNODE=LNODS(INODE)
        IGLVEC(IROWS,LNODE)=0
       END IF
      ENDDO
C
      RETURN
      END
