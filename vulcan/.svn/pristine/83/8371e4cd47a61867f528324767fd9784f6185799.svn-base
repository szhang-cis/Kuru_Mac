      SUBROUTINE LUMPMT(WSTIF,
     .                  NKOVA,NNODL,NDOFN,KSYMM,NEVAB)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A MATRIX WITH NDOFN
C     DEGREE OF FREEDOM PER NODE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION WSTIF(NKOVA)
C
      SUMDI=0.0
      SUMTO=0.0
C
      DO INODE=1,NNODL
       KEVAB=(INODE-1)*NDOFN
       INDEX=INODE                ! symmetric case
       IF(KSYMM.EQ.0) INDEX=1     ! unsymmetric case
       DO JNODE=INDEX,NNODL
        JEVAB=(JNODE-1)*NDOFN
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
C
         IF(IEVAB.EQ.JEVAB) THEN
          SUMDI=SUMDI+WSTIF(KLOCA)
          SUMTO=SUMTO+WSTIF(KLOCA)
         ELSE
          IF(KSYMM.EQ.0) THEN      ! unsymmetric case
           SUMTO=SUMTO+WSTIF(KLOCA)
          ELSE                     ! symmetric case
           SUMTO=SUMTO+2*WSTIF(KLOCA)
          ENDIF       
         ENDIF       
        ENDDO
       ENDDO
      ENDDO
C
      SUTOL=1.0D-10
      FACLU=0.0
      IF(SUMDI.GT.SUTOL) FACLU=SUMTO/SUMDI
C
      DO INODE=1,NNODL
       KEVAB=(INODE-1)*NDOFN
       INDEX=INODE                ! symmetric case
       IF(KSYMM.EQ.0) INDEX=1     ! unsymmetric case
       DO JNODE=INDEX,NNODL
        JEVAB=(JNODE-1)*NDOFN
        DO IDOFN=1,NDOFN
         IEVAB=KEVAB+IDOFN
         JEVAB=JEVAB+1
         KLOCA=(2*NEVAB-IEVAB)*(IEVAB-1)/2+JEVAB
C
         IF(IEVAB.EQ.JEVAB) THEN
          WSTIF(KLOCA)=WSTIF(KLOCA)*FACLU
         ELSE
          WSTIF(KLOCA)=0.0
         ENDIF       
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
