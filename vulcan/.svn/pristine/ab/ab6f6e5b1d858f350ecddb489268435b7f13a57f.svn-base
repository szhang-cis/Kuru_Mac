      SUBROUTINE STIPAR(DISPR,IFFIX,TLOAD,REFOR)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE CURRENT STIFFNESS PARAMETER
C     & SAVES SOME CONVERGED VALUES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DISPR(*), IFFIX(*), TLOAD(NTOTV,*), REFOR(NTOTV)
C
      STIOL=STICU
      FACPR=FACTO
C
      PROD1=0.0
      PROD2=0.0
      PROD3=0.0
C
C**** STORE PREVIOUS VALUE OF TOTAL LOAD & ...
C
      DO 10 ITOTV=1,NTOTV
C
      DISIN=DISPR(ITOTV)           ! Increm. displ.
      TOTLO=TLOAD(ITOTV,1)         ! Total load
      IF(IFFIX(ITOTV).EQ.1) TOTLO=TOTLO-REFOR(ITOTV) ! Add reaction
      DEFOR=TOTLO-TLOAD(ITOTV,2)   ! Increm. load
      TLOAD(ITOTV,2)=TOTLO         ! Save total load
      PROD1=PROD1+DEFOR*DISIN
      PROD2=PROD2+DEFOR*DEFOR
      PROD3=PROD3+DISIN*DISIN
C
   10 CONTINUE
C
      IF(KDYNA.EQ.1.OR.KPORE.EQ.2) RETURN
C
C**** COMPUTE CURRENT STIFFNESS PARAMETER AND ARC-LENGTH
C
      IF(KARCL.EQ.4) THEN                      ! FOR DISPL. CONTROL
        ARCLN=DISPR(NCDIS)
      ELSE                                     ! FOR ARC-LENGTH
        ARCLN=DSQRT(PROD3)
      ENDIF
C
      IF(PROD1.NE.0.0) THEN
                              STICU=PROD2/PROD1
      ELSE
                              STICU=1.0D00
      ENDIF
      IF(ISTEP.EQ.1)          STIIN=STICU
      IF(STIIN.NE.0.0)        STICU=STICU/STIIN
      IF(ISTEP.EQ.1)          STIOL=STICU
      STICH=DABS(STIOL-STICU)
      IF(STIOL.NE.0.0)        RATIO=100.0*STICH/STIOL
      IF(ABS(RATIO).GE.100.0) RATIO=999.9999*RATIO/ABS(RATIO)
C
      WRITE(LURES,900) STICU,STICH,RATIO
C
      IF(ISTEP.EQ.1) THEN
        PITER=DITER
      ELSE
        IF(LAUTO.EQ.1) PITER=FLOAT(IITER)
        IF(LAUTO.EQ.2) PITER=STICH
        IF(LAUTO.EQ.3) PITER=STICH/STIOL
        IF(LAUTO.EQ.4) PITER=STICU/STIOL
      ENDIF
C
      RETURN
  900 FORMAT(/,132('='),/
     .         5X,'CURRENT STIFF. PARAMETER =',F10.5,
     .    2X,'***',2X,'STIFF. CHANGE =',F10.5,2X,'( ',F9.4,' % )',/,
     .        132('='),/)
      END
