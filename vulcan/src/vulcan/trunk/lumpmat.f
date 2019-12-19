      SUBROUTINE LUMPMAT(WSTIFT,WSTI2,WSTI1)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE DIAGONALIZATION OF THE PHASE-CHANGE
C     MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION WSTIFT(*),WSTI2(*),WSTI1(*)
C
      SUMDI=0.0
      SUMTO=0.0
C
      DO INODET=1,NNODLT
       DO JNODET=INODET,NNODLT
        KLOCA=(2*NNODLT-INODET)*(INODET-1)/2+JNODET
        WSTI2(KLOCA)=WSTI2(KLOCA)-WSTI1(KLOCA)
        IF(INODET.EQ.JNODET) THEN
         SUMDI=SUMDI+WSTI2(KLOCA)
         SUMTO=SUMTO+WSTI2(KLOCA)
        ELSE
         SUMTO=SUMTO+2*WSTI2(KLOCA)
        ENDIF       
       ENDDO
      ENDDO
C
      SUTOL=1.0D-10
      FACLU=0.0
      IF(SUMDI.GT.SUTOL) FACLU=SUMTO/SUMDI
C
      DO INODET=1,NNODLT
       DO JNODET=INODET,NNODLT
        KLOCA=(2*NNODLT-INODET)*(INODET-1)/2+JNODET
        IF(INODET.EQ.JNODET) THEN
         WSTI2(KLOCA)=WSTI2(KLOCA)*FACLU
        ELSE
         WSTI2(KLOCA)=0.0
        ENDIF       
        WSTIFT(KLOCA)=WSTIFT(KLOCA)+WSTI2(KLOCA)
       ENDDO
      ENDDO
C
      RETURN
      END
