      SUBROUTINE WMHEATT(DVOLU,NDOFN,NEVAB,NNODE,PROPS,SHAPE,WSTIF,
     .                   BASCC,BASMM,COUTT,WHAPE,KSYMM,EMATX)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE CONSISTENT HEAT CAPACITY ELEMENT
C     MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PROPS(*), SHAPE(*), WSTIF(*),
     .          WHAPE(*)
      DIMENSION EMATX(NEVAB,*)
C
      COEFC=BASCC
      COEFM=BASMM
      CONST=COEFC*COEFM
C
C**** INITIALISES E
C
      DO IEVAB=1,NEVAB
       DO JEVAB=1,NEVAB
        EMATX(IEVAB,JEVAB)=0.0D0
       ENDDO
      ENDDO
C
      DO INODE=1,NNODE
       INDEX=INODE               ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1    ! UNSYMMETRIC CASE
       DO JNODE=INDEX,NNODE
        EMATX(INODE,JNODE)=EMATX(INODE,JNODE)+WHAPE(INODE)*
     .                     CONST*SHAPE(JNODE)
       ENDDO
      ENDDO
C
C**** ADD THE CONTRIBUTION TO THE CAPACITY MATRIX
C
      CALL SQTOTA(WSTIF,EMATX,DVOLU,NEVAB,KSYMM)
C
      RETURN
      END
c-----------------------------------------------------------------
cctm
cctm   para agregar el termino de acoplamiento que no depende de la
cctm   deformacion termica (la matriz C_p en la tesis)
cctm
cctm      CONST=(COEFC*COEFM+COUTT)*WHAPE(INODE)*SHAPE(JNODE)
cctm
c      CONST=COEFC*COEFM*WHAPE(INODE)*SHAPE(JNODE)
cC
c
c       OJO: es necesario definir COUTT (en capcoft o en otra rutina)
c
