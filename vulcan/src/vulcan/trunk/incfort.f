      SUBROUTINE INCFORT(HEADST,IFFIXT,RLOADT,RLOAHT,FICTOT,TLOADT)
C***********************************************************************
C
C****THIS ROUTINE INCREMENTS THE APPLIED FORCE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
C
      DIMENSION HEADST(NPOINT,*),      IFFIXT(*),      RLOADT(*),
     .          RLOAHT(NTOTVT,NFUNCT), FICTOT(NFUNCT), TLOADT(*)
C
      DO IDOFNT=1,NDOFNT
C$DIR NO_RECURRENCE
        DO IPOINT=1,NPOINT
          ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
          IF(IFFIXT(ITOTVT).EQ.0) THEN
           DO IFUNCT=1,NFUNCT
            RLOA1T=RLOAHT(ITOTVT,IFUNCT)*FICTOT(IFUNCT)
            TLOADT(ITOTVT)=TLOADT(ITOTVT)+RLOA1T
           ENDDO
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
