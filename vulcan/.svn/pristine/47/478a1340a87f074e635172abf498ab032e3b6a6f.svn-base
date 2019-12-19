      SUBROUTINE LUMPMTT(WSTIFT,
     .                   NKOVAT,NNODLT,KSYMMT,NEVABT)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A MATRIX WITH ONE
C     DEGREE OF FREEDOM PER NODE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION WSTIFT(NKOVAT)
C
      SUMDI=0.0
      SUMTO=0.0
C
      DO INODET=1,NNODLT
       INDEXT=INODET              ! symmetric case
       IF(KSYMMT.EQ.0) INDEXT=1   ! unsymmetric case
       DO JNODET=INDEXT,NNODLT
        KLOCA=(2*NEVABT-INODET)*(INODET-1)/2+JNODET
        IF(INODET.EQ.JNODET) THEN
         SUMDI=SUMDI+WSTIFT(KLOCA)
         SUMTO=SUMTO+WSTIFT(KLOCA)
        ELSE
         IF(KSYMMT.EQ.0) THEN     ! unsymmetric case
          SUMTO=SUMTO+WSTIFT(KLOCA)
         ELSE                     ! symmetric case
          SUMTO=SUMTO+2*WSTIFT(KLOCA)
         ENDIF       
        ENDIF       
       ENDDO
      ENDDO
C
      SUTOL=1.0D-10
      FACLU=0.0
      IF(SUMDI.GT.SUTOL) FACLU=SUMTO/SUMDI
C
      DO INODET=1,NNODLT
       INDEXT=INODET              ! symmetric case
       IF(KSYMMT.EQ.0) INDEXT=1   ! unsymmetric case
       DO JNODET=INDEXT,NNODLT
        KLOCA=(2*NEVABT-INODET)*(INODET-1)/2+JNODET
        IF(INODET.EQ.JNODET) THEN
         WSTIFT(KLOCA)=WSTIFT(KLOCA)*FACLU
        ELSE
         WSTIFT(KLOCA)=0.0
        ENDIF       
       ENDDO
      ENDDO
C
      RETURN
      END
