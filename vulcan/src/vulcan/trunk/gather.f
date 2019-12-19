      SUBROUTINE GATHER(GLVEC,NGLOB,NPOIN,ELVEC,NROWS,NNODE,LNODS)
C***********************************************************************
C
C****THIS ROUTINE PERFORMS SCALAR GATHER OPERATIONS 
C          GLVEC(NGLOB,NPOIN) ---> ELVEC(NROWS,NNODE)
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ELVEC(NROWS*NNODE), GLVEC(NGLOB,NPOIN), LNODS(NNODE)
C
C$DIR NO_RECURRENCE
      DO K=1,NROWS*NNODE
        INODE=(K-1)/NROWS+1
        IROWS=K-(INODE-1)*NROWS
        LNODE=LNODS(INODE)
        ELVEC(K)=GLVEC(IROWS,LNODE)
      ENDDO
C
      RETURN
      END
