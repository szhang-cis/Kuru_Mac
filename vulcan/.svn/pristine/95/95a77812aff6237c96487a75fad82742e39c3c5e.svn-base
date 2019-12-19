      SUBROUTINE INCFORS(HEADSS,IFFIXS,RLOADS,RLOAHS,FICTOS,TLOADS)
C***********************************************************************
C
C****THIS ROUTINE INCREMENTS THE APPLIED FORCE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
C
      DIMENSION HEADSS(NPOINS,*),      IFFIXS(*),      RLOADS(*),
     .          RLOAHS(NTOTVS,NFUNCS), FICTOS(NFUNCS), TLOADS(*)
C
      DO IDOFNS=1,NDOFNS
C$DIR NO_RECURRENCE
        DO IPOINS=1,NPOINS
          ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
          IF(IFFIXS(ITOTVS).EQ.0) THEN
           DO IFUNCS=1,NFUNCS
            RLOA1S=RLOAHS(ITOTVS,IFUNCS)*FICTOS(IFUNCS)
            TLOADS(ITOTVS)=TLOADS(ITOTVS)+RLOA1S
           ENDDO
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
