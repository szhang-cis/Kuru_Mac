      SUBROUTINE SCATER(ELVEC,NROWS,NNODE,GLVEC,NGLOB,NPOIN,LNODS)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS SCALAR SCALAR OPERATIONS 
C          ELVEC(NROWS,NNODE) ---> GLVEC(NROWS,NPOIN) 
C    
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ELVEC(NROWS*NNODE), GLVEC(NGLOB,NPOIN), LNODS(NNODE)
C
C$DIR NO_RECURRENCE
      DO K=1,NROWS*NNODE
        INODE=(K-1)/NROWS+1
        IROWS=K-(INODE-1)*NROWS
        LNODE=LNODS(INODE)
        GLVEC(IROWS,LNODE)=GLVEC(IROWS,LNODE)+ELVEC(K)
      ENDDO
C
      RETURN
      END
