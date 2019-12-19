      SUBROUTINE OUTASS(LNODS,SVECT,WVECT,FACTA,
     .                  NNODE,IELEM,NELEM,NPOIN,
     .                  NWWWW)
C***********************************************************************
C
C**** THIS ROUTINE ASSIGNS SOME SMOOTHING VARIABLES
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION LNODS(NNODE,NELEM)
      DIMENSION SVECT(NWWWW,NPOIN), WVECT(NWWWW,NNODE)
C
      DO IWWWW=1,NWWWW
       DO INODE=1,NNODE
        LPOIN=LNODS(INODE,IELEM)
        IF(LPOIN.NE.0) SVECT(IWWWW,LPOIN)=SVECT(IWWWW,LPOIN)+
     .                 WVECT(IWWWW,INODE)*FACTA       ! 23/2/94 ctm
       ENDDO
      ENDDO
C
      RETURN
      END
