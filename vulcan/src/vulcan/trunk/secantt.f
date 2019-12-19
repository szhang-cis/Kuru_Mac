      SUBROUTINE SECANTT(DISITT,DSOLDT,REFORT,RFOLDT)
C***********************************************************************
C
C****THIS ROUTINE MODIFIES THE ITERATIVE DISPLACEMENTS IF A 
C    SECANT-NEWTON METHOD IS USED
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISITT(*), DSOLDT(*), REFORT(*), RFOLDT(*)
      DATA RFAC1T,RFAC2T/4.00,0.25/
C
C***SET DEFAULT VALUES
C
      AACOET=0.0D+00
      BBCOET=0.0D+00
C
C***CHECK IF SECANT-NEWTON METHODS ARE TO BE USED
C
      IF((LACCET.NE.1.AND.LACCET.NE.2).OR.(IITERT.LT.2).OR.
     .   (KARCLT.NE.0)) RETURN
C
C.......SET COEFFICIENTS AND CONSIDER "CUT-OFFS"
C
        DIFFET=GZEROT-GCURNT
        IF(DIFFET.EQ.0.0)                                         RETURN
        AVALUT=GZEROT/DIFFET
        IF((AVALUT.GT.RFAC1T).OR.(AVALUT.LT.1./RFAC1T))           RETURN
        AACOET=AVALUT
C
        IF(LACCET.EQ.2) THEN
          VALUET=0.0
          DO IDOFNT=1,NDOFNT
            DO IPOINT=1,NPOINT
              ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
              VALUET=VALUET+DISITT(ITOTVT)*(REFORT(ITOTVT)-
     .               RFOLDT(ITOTVT))
            ENDDO
          ENDDO
C
          BBCOET=AACOET*(1.0-VALUET/(AALPHT*(GCURNT-GZEROT)))-1.0
          COMPAT=BBCOET/AACOET
          IF((COMPAT.GT.RFAC2T).OR.(COMPAT.LT.-0.50*RFAC2T))
     .        BBCOET=0.0D+00
C
        ENDIF
C
C.......MODIFY DISIT ACCORDING TO AACOE, BBCOE
C
        DO IDOFNT=1,NDOFCT
C$DIR NO_RECURRENCE
          DO IPOINT=1,NPOINT
            ITOTVT=(IPOINT-1)*NDOFCT+IDOFNT
            DISITT(ITOTVT)=AACOET*DISITT(ITOTVT)+BBCOET*AALPHT*
     .                     DSOLDT(ITOTVT)
          ENDDO
        ENDDO
C
      END
