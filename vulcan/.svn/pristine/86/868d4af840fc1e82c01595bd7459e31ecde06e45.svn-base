      SUBROUTINE CONVER(DISIT,DISPR,DISTO,REFOR,TLOAD,IFFIX)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS FOR CONVERGENCE OF THE ITERATION PROCESS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION DISIT(*),       DISPR(*),
     .          DISTO(NTOTV,*), IFFIX(*),
     .          REFOR(NTOTV),   TLOAD(NTOTV,*)
C
      DIMENSION PVALS(4), RATIS(4), SELEC(4), SIGNS(4)
      DIMENSION RATISC(2)
C
      CHARACTER*12 BLANK,CONVE,DIVER,SELEC,SIGNA,SIGNS,CHA12
C
      SAVE PVALS,GINIT
      SAVE PVALUM
C
      DATA CONVE,DIVER,BLANK,SIGNA
     .    /'  D -->> C  ','  D <<-- C  ','            ','  SELECTED  '/ 
C
      ZEROM=1.0E-08
      IF(ZEROM.GT.TOLER) ZEROM=TOLER
C
      CALL CPUTIM(CUCPU)
C
C**** INITIALISE FLAGS TO CONVERGED
C
      NCHEK=0
      IJUMP=0
C
      IF(IITER.EQ.1) THEN
       CHA12=BLANK
      ELSE
       CHA12=CONVE
      END IF
      DO I=1,4
       RATIS(I)=0.0D0
       SELEC(I)=BLANK
       SIGNS(I)=CHA12
      END DO
      SELEC(KCONV)=SIGNA
C
C**** COMPUTE RELEVANT NORMS
C
      RESID=0.0D0
      RETO1=0.0D0
      RETO2=0.0D0
      DIST1=0.0D0
      DIST2=0.0D0
      GITER=0.0D0
      GTOTL=0.0D0
C
      IF(NPOIC.GT.0) THEN
c      DO I=1,2
c       RATISC(I)=0.0D0
c      END DO
C
       RESIDC=0.0D0
       RETO1C=0.0D0
       DIST1C=0.0D0
       DIST2C=0.0D0
      ENDIF
C
      DO IDOFN=1,NDOFN
       DO IPOIN=1,NPOIN
        ITOTV=(IPOIN-1)*NDOFC+IDOFN
        RESLD=REFOR(ITOTV)
        TOTL1=TLOAD(ITOTV,1)
        IF(IFFIX(ITOTV).EQ.1) THEN
         TOTL1=TOTL1-RESLD               ! Add reaction
         RESLD=0.0D0
        ENDIF
        TOTL2=TOTL1-TLOAD(ITOTV,2)
C
        DIITE=DISIT(ITOTV)
        DIINC=DISPR(ITOTV)
        DTOTL=DISTO(ITOTV,1)
C
        IF(NPOIC.GT.0) THEN
         NPOI1=NPOIN-NPOIC
         IF(IPOIN.GT.NPOI1) THEN
C
C**** CORRECTION OF THE CONTACT FORCE FOR THE ALM (IAUGM=2)
C
          IF(IAUGM.EQ.2.OR.IAUGM.EQ.3)
     .     DISTO(ITOTV,1)=DISTO(ITOTV,1)+DISIT(ITOTV)
C
          DIITE=0.0D0
          DIINC=0.0D0
          DTOTL=0.0D0
          RESLD=0.0D0
          TOTL1=0.0D0
          TOTL2=0.0D0
C
          DIITEC=DISIT(ITOTV)
          DIINCC=DISPR(ITOTV)
          IF(IAUGM.EQ.2.OR.IAUGM.EQ.3) DIINCC=DIINCC+DIITEC
          DTOTLC=DISTO(ITOTV,1)
          RESLDC=REFOR(ITOTV)
C
          RESIDC=RESIDC+RESLDC*RESLDC
          RETO1C=RETO1C+DTOTLC*DTOTLC
          DIST1C=DIST1C+DIITEC*DIITEC
          DIST2C=DIST2C+DIINCC*DIINCC
         ENDIF
        ENDIF
C
        RESID=RESID+RESLD*RESLD
        RETO1=RETO1+TOTL1*TOTL1
        RETO2=RETO2+TOTL2*TOTL2
C
        DIST1=DIST1+DIITE*DIITE
        DIST2=DIST2+DIINC*DIINC
C
        GITER=GITER+RESLD*DIITE
        GTOTL=GTOTL+TOTL1*DTOTL
C
       END DO
      END DO
C
C**** FOR PROBLEMS WITH ZERO EXTERNAL & REACTION FORCES
C
      IF(RETO1.LE.ZEROM.AND.KCONR.EQ.1) RETO1=ZEROM
      IF(RETO2.LE.ZEROM.AND.KCONR.EQ.1) RETO2=ZEROM
      IF(RESID.LE.ZEROM.AND.RESID.LE.RETO1.AND.
     .                      RESID.LE.RETO2) RESID=0.0D0
C
C**** COMPUTE NORM RATIOS
C
      IF(RETO1.GE.ZEROM)       RATIS(1)=100.0D0*DSQRT(RESID/RETO1)
      IF(RETO2.GE.ZEROM)       RATIS(2)=100.0D0*DSQRT(RESID/RETO2)
      IF(DIST2.GT.ZEROM)       RATIS(3)=100.0D0*DSQRT(DIST1/DIST2)
      IF(DABS(GTOTL).GT.ZEROM) RATIS(4)=100.0D0*DABS(GITER/GTOTL)
      IF(NPOIC.GT.0) THEN
       if(iaugm.eq.3) then
        if(iaug3.eq.0) then
         IF(RETO1C.GE.ZEROM)   RATISC(1)=100.0D0*DSQRT(RESIDC/RETO1C)
         IF(DIST2C.GT.ZEROM)   RATISC(2)=100.0D0*DSQRT(DIST1C/DIST2C)
        endif
       else
        IF(RETO1C.GE.ZEROM)    RATISC(1)=100.0D0*DSQRT(RESIDC/RETO1C)
        IF(DIST2C.GT.ZEROM)    RATISC(2)=100.0D0*DSQRT(DIST1C/DIST2C)
       endif
      ENDIF
C
C**** CONTACT CONTROL (ALM)
C
      IF(NPOIC.GT.0) THEN
       if(iiter.eq.1) then
        call solcon(ratis(1))
c
c**** itruc < 0 ==> only tsoli varies (not nsol_i)
c
        if(itrucm.gt.0) then
         isolim=nsol1m-itrucm
         nsol1m=nsol1m-isolim
         nsol2m=nsol2m-isolim
         nsol3m=nsol3m-isolim
        endif
       endif
c
c**** options (better as input): a) tsoli=1.0  b) tsoli=truch
c
       tsolim=1.0d0
       if(iiter.ge.nsol1m) then
        if(iiter.lt.nsol2m) tsolim=tsol1m+(tsol2m-tsol1m)/
     .                      (nsol2m-nsol1m)*(iiter-nsol1m)
        if(iiter.ge.nsol2m.and.iiter.lt.nsol3m) tsolim=tsol2m+
     .                    (tsol3m-tsol2m)/(nsol3m-nsol2m)*(iiter-nsol2m)
        if(iiter.ge.nsol3m) tsolim=tsol3m
        tsolim=tsolim*truchm
       endif
c
       DO IDOFN=1,NDOFN
        DO IPOIN=1,NPOIN
         ITOTV=(IPOIN-1)*NDOFC+IDOFN
         refor(itotv)=refor(itotv)*tsolim
        ENDDO
       ENDDO
      ENDIF
C
      DO I=1,4
       IF(RATIS(I).GT.999.9999D0) RATIS(I)=999.9999D0
      END DO
C
C**** CHECK FOR CONVERGENCE WITH SPECIAL REFERENCE TO THE SELECTED ONE
C
      RATIO=RATIS(KCONV)
      IF(RATIO.GT.TOLER) NCHEK=1
C
      IF(IITER.GT.1) THEN
       DO 6 I=1,4
    6  IF(RATIS(I).GT.PVALS(I)) SIGNS(I)=DIVER
       IF(NCHEK.EQ.1.AND.RATIO.GT.PVALUM) NCHEK=999
      ENDIF
C
      DO 7 I=1,4
    7 PVALS(I)=RATIS(I)
      PVALUM=RATIO
C
C>>>> JUMPING TO NEXT STEP
C
      IF((NCHEK.NE.0).AND.(IITER.GT.MITER/2).AND.
     .                     (RATIO.LT.2.0D0*TOLER)) THEN
       NCHEK=0
       IJUMP=1
c      WRITE(LUPRI,905)
      ENDIF
C
C**** CHECKS CONVERGENCE OF CONTACT FORCES FOR ALM
C
      IF(NPOIC.GT.0) THEN
       IF(NCHEK.EQ.0.AND.IITER.EQ.1) THEN
        IF(RATISC(2).GT.TOLER) THEN
         NCHEK=1
         IJUMP=0
        ENDIF
       ENDIF
      ENDIF
C
C**** CONTROLS IMPENETRABILITY FOR ALM (IAUGM=3)
C
      IF(NPOIC.GT.0) THEN
       IF(IAUGM.EQ.3) THEN
        IAUG3=0
        IF(NCHEK.EQ.0) THEN
         IF(RATISC(2).GT.TOLERC) THEN
          NCHEK=1
          IJUMP=0
          IAUG3=0
         ENDIF
        ELSE
         IAUG3=1
        ENDIF
       ENDIF
      ENDIF
C
      IF(IJUMP.EQ.1) WRITE(LUPRI,905)
C
C**** PRINT RESULTS
C
      WRITE(LURES,900) (RATIS(I),SIGNS(I),SELEC(I),I=1,4)
      IF(NPOIC.GT.0) WRITE(LURES,1900) (RATISC(I),I=1,2)
C
      IF(NCHEK.EQ.0)   WRITE(LURES,901) CUCPU-CPUIN
      IF(NCHEK.EQ.1)   WRITE(LURES,902) CUCPU-CPUIN
      IF(NCHEK.EQ.999) WRITE(LURES,903) CUCPU-CPUIN
      IF(IJUMP.EQ.1)   WRITE(LURES,904) CUCPU-CPUIN
      WRITE(LUPRI,910) ITIME,ISTEP,IITER,CUCPU-CPUIN,RATIO/TOLER
      WRITE(LUINF,910) ITIME,ISTEP,IITER,CUCPU-CPUIN,RATIO/TOLER
C
      RETURN
  900 FORMAT(/10X,'CONVERGENCE CHECKS :',//,
     .    15X,'RESIDUAL FORCES/TOTAL   FORCES =',F10.4,' ( % )',2A12/,
     .    15X,'RESIDUAL FORCES/INCRMT. FORCES =',F10.4,' ( % )',2A12/,
     .    15X,'ITERATV. DISPL./INCRMT. DISPL. =',F10.4,' ( % )',2A12/,
     .    15X,'ITERATV. ENERGY/TOTAL   ENERGY =',F10.4,' ( % )',2A12/)
  901 FORMAT(10X,'>>>>>>>>> CONVERGED ',10X,'CPU TIME =',F10.3)
  902 FORMAT(10X,'>>>>>>>>> CONVERGING',10X,'CPU TIME =',F10.3)
  903 FORMAT(10X,'>>>>>>>>> DIVERGING ',10X,'CPU TIME =',F10.3)
  904 FORMAT(10X,'>>>>>>>>> SKIPPING  ',10X,'CPU TIME =',F10.3)
  905 FORMAT(10X,'>>>>>>>>> JUMPING TO NEXT STEP')
  910 FORMAT(5X,'ITIME =',I5,2X,'ISTEP =',I5,2X,'IITER =',I5,2X,/,
     .      10X,'CUCPU =',F10.3,2X,'RATIO/TOLER =',E15.6/)
C
 1900 FORMAT(
     .    15X,'RESIDUAL C.CON./TOTAL C.FORCES =',F10.4,' ( % )',/,
     .    15X,'ITERATV. C.FOR./INCRMT. C.FOR. =',F10.4,' ( % )',/)
      END
