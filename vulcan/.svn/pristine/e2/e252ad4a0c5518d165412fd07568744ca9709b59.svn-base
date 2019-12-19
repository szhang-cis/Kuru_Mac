      SUBROUTINE SECANT(DISIT,DSOLD,REFOR,RFOLD)
C***********************************************************************
C
C****THIS ROUTINE MODIFIES THE ITERATIVE DISPLACEMENTS IF A 
C    SECANT-NEWTON METHOD IS USED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISIT(*), DSOLD(*), REFOR(*), RFOLD(*)
      DATA RFAC1,RFAC2/4.00,0.25/
C
C***SET DEFAULT VALUES
C
      AACOE=0.0D+00
      BBCOE=0.0D+00
C
C***CHECK IF SECANT-NEWTON METHODS ARE TO BE USED
C
      IF((LACCE.NE.1.AND.LACCE.NE.2).OR.(IITER.LT.2).OR.(KARCL.NE.0))
     .                                                            RETURN
C
C.......SET COEFFICIENTS AND CONSIDER "CUT-OFFS"
C
        DIFFE=GZERO-GCURN
        IF(DIFFE.EQ.0.0)                                          RETURN
        AVALU=GZERO/DIFFE
        IF((AVALU.GT.RFAC1).OR.(AVALU.LT.1./RFAC1))               RETURN
        AACOE=AVALU
C
        IF(LACCE.EQ.2) THEN
          VALUE=0.0
          DO IDOFN=1,NDOFN
            DO IPOIN=1,NPOIN
              ITOTV=(IPOIN-1)*NDOFC+IDOFN
              VALUE=VALUE+DISIT(ITOTV)*(REFOR(ITOTV)-RFOLD(ITOTV))
            ENDDO
          ENDDO
C
          BBCOE=AACOE*(1.0-VALUE/(AALPH*(GCURN-GZERO)))-1.0
          COMPA=BBCOE/AACOE
          IF((COMPA.GT.RFAC2).OR.(COMPA.LT.-0.50*RFAC2)) BBCOE=0.0D+00
C
        ENDIF
C
C.......MODIFY DISIT ACCORDING TO AACOE, BBCOE
C
        DO IDOFN=1,NDOFC
C$DIR NO_RECURRENCE
          DO IPOIN=1,NPOIN
            ITOTV=(IPOIN-1)*NDOFC+IDOFN
            DISIT(ITOTV)=AACOE*DISIT(ITOTV)+BBCOE*AALPH*DSOLD(ITOTV)
          ENDDO
        ENDDO
C
      END
