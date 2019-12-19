      SUBROUTINE LUMPMMT(WSTIFT,
     .                   NNODLT,KSYMMT)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A MATRIX WITH
C     DIMENSION NNODLT
C     (idem lumpmtt.f but WSTIFT is written in matrix form)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION WSTIFT(NNODLT,*)
C
      SUMDI=0.0D0
      SUMTO=0.0D0
C
      DO INODET=1,NNODLT
       INDEXT=INODET              ! symmetric case
       IF(KSYMMT.EQ.0) INDEXT=1   ! unsymmetric case
       DO JNODET=INDEXT,NNODLT
        IF(INODET.EQ.JNODET) THEN
         SUMDI=SUMDI+WSTIFT(INODET,JNODET)
         SUMTO=SUMTO+WSTIFT(INODET,JNODET)
        ELSE
         IF(KSYMMT.EQ.0) THEN     ! unsymmetric case
          SUMTO=SUMTO+WSTIFT(INODET,JNODET)
         ELSE                     ! symmetric case
          SUMTO=SUMTO+2*WSTIFT(INODET,JNODET)
         ENDIF       
        ENDIF       
       ENDDO
      ENDDO
C
      SUTOL=1.0D-10*SUMTO
      FACLU=0.0D0
      IF(SUMDI.GT.SUTOL) FACLU=SUMTO/SUMDI
C
      DO INODET=1,NNODLT
       INDEXT=INODET              ! symmetric case
       IF(KSYMMT.EQ.0) INDEXT=1   ! unsymmetric case
       DO JNODET=INDEXT,NNODLT
        IF(INODET.EQ.JNODET) THEN
         WSTIFT(INODET,JNODET)=WSTIFT(INODET,JNODET)*FACLU
        ELSE
         WSTIFT(INODET,JNODET)=0.0D0
        ENDIF       
       ENDDO
      ENDDO
C
      RETURN
      END
