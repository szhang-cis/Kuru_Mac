      SUBROUTINE SECANTS(DISITS,DSOLDS,REFORS,RFOLDS)
C***********************************************************************
C
C****THIS ROUTINE MODIFIES THE ITERATIVE DISPLACEMENTS IF A 
C    SECANT-NEWTON METHOD IS USED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION DISITS(*), DSOLDS(*), REFORS(*), RFOLDS(*)
      DATA RFAC1S,RFAC2S/4.00,0.25/
C
C***SET DEFAULT VALUES
C
      AACOES=0.0D+00
      BBCOES=0.0D+00
C
C***CHECK IF SECANT-NEWTON METHODS ARE TO BE USED
C
      IF((LACCES.NE.1.AND.LACCES.NE.2).OR.(IITERS.LT.2).OR.
     .   (KARCLS.NE.0)) RETURN
C
C.......SET COEFFICIENTS AND CONSIDER "CUT-OFFS"
C
        DIFFES=GZEROS-GCURNS
        IF(DIFFES.EQ.0.0)                                         RETURN
        AVALUS=GZEROS/DIFFES
        IF((AVALUS.GT.RFAC1S).OR.(AVALUS.LT.1./RFAC1S))           RETURN
        AACOES=AVALUS
C
        IF(LACCES.EQ.2) THEN
          VALUES=0.0
          DO IDOFNS=1,NDOFNS
            DO IPOINS=1,NPOINS
              ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
              VALUES=VALUES+DISITS(ITOTVS)*(REFORS(ITOTVS)-
     .               RFOLDS(ITOTVS))
            ENDDO
          ENDDO
C
          BBCOES=AACOES*(1.0-VALUES/(AALPHS*(GCURNS-GZEROS)))-1.0
          COMPAS=BBCOES/AACOES
          IF((COMPAS.GT.RFAC2S).OR.(COMPAS.LT.-0.50*RFAC2S))
     .        BBCOES=0.0D+00
C
        ENDIF
C
C.......MODIFY DISIT ACCORDING TO AACOE, BBCOE
C
        DO IDOFNS=1,NDOFCS
C$DIR NO_RECURRENCE
          DO IPOINS=1,NPOINS
            ITOTVS=(IPOINT-1)*NDOFCS+IDOFNS
            DISITS(ITOTVS)=AACOES*DISITS(ITOTVS)+BBCOES*AALPHS*
     .                     DSOLDS(ITOTVS)
          ENDDO
        ENDDO
C
      END
