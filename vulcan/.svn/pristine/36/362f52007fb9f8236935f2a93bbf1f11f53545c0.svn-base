      SUBROUTINE WBOUNCT(DVOLU,NDOFN,NEVOB,NNODE,PROPS,SHAPE,WSTIF,
     .                   BASHH,ELDIS,NDOFC,KSYMM,EMATX,NEVAB,FMATX,
     .                   IAUXY,LNODS,IGABO)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE CONSISTENT BOUNDARY ELEMENT MATRIX
C     (FOR CONTACT ELEMENTS FOR NON-COINCIDENT MESHES)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PROPS(*), SHAPE(*), WSTIF(*), ELDIS(NDOFC,*)
      DIMENSION EMATX(NEVOB,*), FMATX(NEVAB,*)
      DIMENSION LNODS(*)
C
      COEFH=BASHH
      CONST=COEFH
C
C**** INITIALISES F
C
C     Note: this operation is necessary because the non-zero components
C           of the stiffness matrix are different, in general, for each
C           integration point
C
      DO IEVAB=1,NEVAB
       DO JEVAB=1,NEVAB
        FMATX(IEVAB,JEVAB)=0.0D0
       ENDDO
      ENDDO
C
C**** INITIALISES E
C
      DO IEVAB=1,NEVOB
       DO JEVAB=1,NEVOB
        EMATX(IEVAB,JEVAB)=0.0D0
       ENDDO
      ENDDO
C
C**** COMPUTES Nmaster^T h Nmaster
C
      DO INODE=1,NNODE
       INDEX=INODE               ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1    ! UNSYMMETRIC CASE
       DO JNODE=INDEX,NNODE
        EMATX(INODE,JNODE)=EMATX(INODE,JNODE)+SHAPE(INODE)*
     .                     CONST*SHAPE(JNODE)
       ENDDO
      ENDDO
C
C**** ADDS TO FMATX
C
      DO INODE=1,NNODE
       INDEX=INODE               ! SYMMETRIC CASE
       IF(KSYMM.EQ.0) INDEX=1    ! UNSYMMETRIC CASE
       DO JNODE=INDEX,NNODE
        FMATX(INODE,JNODE)=EMATX(INODE,JNODE)
       ENDDO
      ENDDO
C
      IF(IGABO.EQ.0) THEN        ! only for AIRGAP (skips BOUNDARY)
C
C**** INITIALISES E
C
       DO IEVAB=1,NEVOB
        DO JEVAB=1,NEVOB
         EMATX(IEVAB,JEVAB)=0.0D0
        ENDDO
       ENDDO
C
C**** COMPUTES Nmaster^T h Nslave
C
       DO INODE=1,NNODE
        INDEX=INODE              ! SYMMETRIC CASE
        IF(KSYMM.EQ.0) INDEX=1   ! UNSYMMETRIC CASE
        DO JNODE=INDEX,NNODE
         EMATX(INODE,JNODE)=EMATX(INODE,JNODE)+SHAPE(INODE)*
     .                      CONST*SHAPE(JNODE+NNODE)
        ENDDO
       ENDDO
C
C**** ADDS TO FMATX
C
       DO INODE=1,NNODE
        INDEX=INODE              ! SYMMETRIC CASE
        IF(KSYMM.EQ.0) INDEX=1   ! UNSYMMETRIC CASE
        DO JNODE=INDEX,NNODE
         JNODX=JNODE+NNODE+IAUXY
c        IF(IAUXY.GT.0) THEN     ! deals with possible repeated nodes
c         DO KNODE=NNODE+1,NNODE+IAUXY      ! it seems to be unnecessary
c          LNODE=NNODE+IAUXY+JNODE
c          IF(LNODS(KNODE).EQ.LNODS(LNODE)) JNODX=KNODE
c         ENDDO
c        ENDIF
         FMATX(INODE,JNODX)=-EMATX(INODE,JNODE)
        ENDDO
       ENDDO
      ENDIF                      ! igabo.eq.0
C
C**** ADD THE CONTRIBUTION TO THE INTERFACE MATRIX
C
      CALL SQTOTA(WSTIF,FMATX,DVOLU,NEVAB,KSYMM)
C
      RETURN
      END
