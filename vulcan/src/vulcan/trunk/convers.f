      SUBROUTINE CONVERS(DISITS,DISPRS,DISTOS,REFORS,TLOADS,IFFIXS)
C***********************************************************************
C
C**** THIS ROUTINE CHECKS FOR CONVERGENCE OF THE ITERATION PROCESS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'                  ! variable strategy (Oct/97)
C
      DIMENSION DISITS(*),        DISPRS(*),        DISTOS(NTOTVS,*), 
     .          REFORS(NTOTVS,2), TLOADS(NTOTVS,*), IFFIXS(*)
C
      DIMENSION PVALS(4),RATIS(4),SELEC(4),SIGNS(4)
C
      CHARACTER*12 BLANK,CONVE,DIVER,SELEC,SIGNA,SIGNS
C
      SAVE PVALS,GINIT
      SAVE RETO4,DDDDD
      SAVE PVALU
      save mitertx,tolertx,itrucx,truchx
C
      DATA CONVE,DIVER,BLANK,SIGNA
     .    /'  D -->> C  ','  D <<-- C  ','            ','  SELECTED  '/ 
C
      ZEROM=1.0E-08
C
      CALL CPUTIMS(CUCPUS)
C
C**** INITIALISE FLAGS
C
      NCHEKS=0
      IJUMPS=0
C
      DO 1 I=1,4
      RATIS(I)=0.0D0
   1  SELEC(I)=BLANK
      SELEC(KCONVS)=SIGNA
C
      IF(IITERS.EQ.1) THEN
       DO 2 I=1,4
   2   SIGNS(I)=BLANK
      ELSE
       DO 3 I=1,4
   3   SIGNS(I)=CONVE
      ENDIF
C
C**** COMPUTE RELEVANT NORMS
C
      RESID=0.0D0
      RETO1=0.0D0
      RETO11=0.0D0
      RETO2=0.0D0
      RETO3=0.0D0
      DIST1=0.0D0
      DIST2=0.0D0
      GITER=0.0D0
      GTOTL=0.0D0
C
      DO 10 IDOFNS=1,NDOFNS
      DO 11 IPOINS=1,NPOINS
      ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
C
      RESLD=REFORS(ITOTVS,1)
      TOTL1=TLOADS(ITOTVS,1)              ! not used in thermal problems
      TOTL11=REFORS(ITOTVS,2)
C
      IF(IFFIXS(ITOTVS).EQ.1) THEN
       TOTL1=TOTL1-RESLD                  ! Add reaction; not used
       RESLD=0.0D0
      ENDIF
      TOTL2=TOTL1-TLOADS(ITOTVS,2)        ! not used
      TOTL3=REFORS(ITOTVS,1)+REFORS(ITOTVS,2)
C
      DIITE=DISITS(ITOTVS)
      DIINC=DISPRS(ITOTVS)
      DTOTL=DISTOS(ITOTVS,1)
C
      RESID=RESID+RESLD*RESLD
      RETO11=RETO11+TOTL11*TOTL11
      RETO1=RETO1+TOTL1*TOTL1
      RETO2=RETO2+TOTL2*TOTL2
      RETO3=RETO3+TOTL3*TOTL3
C
      DIST1=DIST1+DIITE*DIITE
      DIST2=DIST2+DIINC*DIINC
C
      IF(IITERS.EQ.1) THEN
       RETO4=RESID
       DDDDD=DIST1
      ENDIF
C
      GITER=GITER+RESLD*DIITE
      GTOTL=GTOTL+TOTL1*DTOTL
C
   11 CONTINUE
   10 CONTINUE 
C
C**** FOR PROBLEMS WITH ZERO CONDUCTIVITY
C
      IF(RETO11.LE.ZEROM) THEN
       IF(RESID.LE.RETO11) RESID=0.0D0
       RETO11=ZEROM
      ENDIF
C
C**** COMPUTE NORM RATIOS
C
      IF(RETO11.GE.ZEROM)      RATIS(1)=100.0D0*DSQRT(RESID/RETO11)
      IF(RETO4.GT.ZEROM)       RATIS(2)=100.0D0*DSQRT(RESID/RETO4)
      IF(DIST2.GT.ZEROM)       RATIS(3)=100.0D0*DSQRT(DIST1/DIST2)
      IF(DABS(GTOTL).GT.ZEROM) RATIS(4)=100.0D0*DABS(GITER/GTOTL)
      IF(DDDDD.GT.ZEROM)       CONNU=   100.0D0*DSQRT(DIST1/DDDDD)
C
C**** FRONT VELOCITY CONTROL ALGORITHM
C
c     if(iiters.eq.1) then
c      call solcons(ratis(1))
c
c**** read itruc & truch for every step
c
c      if(itrucv.eq.1) then                            ! Oct/97
c       call listent('convert',    1,  83)
c
c       isteptx=paramt(1)
c       mitertx=paramt(2)
c       tolertx=paramt(3)
c       itrucx=paramt(4)
c       truchx=paramt(5)
c
c       if(istept.ne.isteptx)
c    .   call runendt('istept ne isetptx in fort.83   ')
c      endif
c     endif
c
c**** transfers tolerance values
c
c     miterto=mitert                                   ! Oct/97
c     tolerto=tolert
c     itruco=itruc
c     trucho=truch
c     if(itrucv.eq.1) then
c      if(mitertx.gt.0.and.tolertx.gt.0.0) then
c       mitert=mitertx
c       tolert=tolertx
c       itruc=itrucx
c       truch=truchx
c       endif
c      endif
c
c**** itruc < 0 ==> only tsoli varies (not nsol_i)
c
c     if(iitert.eq.1) then
c      if(itruc.gt.0) then
c       isoli=nsol1-itruc
c       nsol1=nsol1-isoli
c       nsol2=nsol2-isoli
c       nsol3=nsol3-isoli
c      endif
c     endif
c
c**** options (better as input): a) tsoli=1.0  b) tsoli=truch
c
      tsoli=1.0d0
c     if(iitert.ge.nsol1) then
c      if(iitert.lt.nsol2) tsoli=tsol1+(tsol2-tsol1)/
c    .                     (nsol2-nsol1)*(iitert-nsol1)
c      if(iitert.ge.nsol2.and.iitert.lt.nsol3) tsoli=tsol2+
c    .                     (tsol3-tsol2)/(nsol3-nsol2)*(iitert-nsol2)
c      if(iitert.ge.nsol3) tsoli=tsol3
c      tsoli=tsoli*truch
c     endif
c
      DO IDOFNS=1,NDOFNS
       DO IPOINS=1,NPOINS
        ITOTVS=(IPOINS-1)*NDOFCS+IDOFNS
        refors(itotvs,1)=refors(itotvs,1)*tsoli
       ENDDO
      ENDDO
C
      DO 5 I=1,4
    5 IF(RATIS(I).GT.999.9999D0) RATIS(I)=999.9999D0
C
C**** CHECK FOR CONVERGENCE WITH SPECIAL REFERENCE TO THE SELECTED ONE
C
      RATIO=RATIS(KCONVS)
      IF(RATIO.GT.TOLERS) NCHEKS=1
C
      IF(IITERS.GT.1) THEN
       DO 6 I=1,4
    6  IF(RATIS(I).GT.PVALS(I)) SIGNS(I)=DIVER
       IF(NCHEKS.EQ.1.AND.RATIO.GT.PVALU) NCHEKT=999
      ENDIF
C
      DO 7 I=1,4
    7 PVALS(I)=RATIS(I)
      PVALU=RATIO
C
C>>>> JUMPING TO NEXT STEP
C
      IF((NCHEKS.NE.0).AND.(IITERS.GT.MITERS/2).AND.
     .                     (RATIO.LT.2.0D0*TOLERS)) THEN
       NCHEKS=0
       IJUMPS=1
       WRITE(LUPRIS,905)
      ENDIF
C
C**** PRINT RESULTS
C
      WRITE(LURESS,900) (RATIS(I),SIGNS(I),SELEC(I),I=1,4)
c     WRITE(LURESS,800) (RATIS(I),SIGNS(I),SELEC(I),I=1,4),CONNU
      IF(NCHEKS.EQ.0)   WRITE(LURESS,901) CUCPUS-CPUINS
      IF(NCHEKS.EQ.1)   WRITE(LURESS,902) CUCPUS-CPUINS
      IF(NCHEKS.EQ.999) WRITE(LURESS,903) CUCPUS-CPUINS
      IF(IJUMPS.EQ.1)   WRITE(LURESS,904) CUCPUS-CPUINS
      WRITE(LUPRIS,910) ITIMES,ISTEPS,IITERS,CUCPUS-CPUINS,RATIO/TOLERS
c
c**** transfers back tolerance values 
c
c     mitert=miterto                        ! Oct/97
c     tolert=tolerto
c     itruc=itruco
c     truch=trucho
C
      RETURN
C
  900 FORMAT(/10X,'CONVERGENCE CHECKS :',//,
     .    15X,'RESIDUAL FLUXES/TOTAL   FLUXES =',F10.4,' ( % )',2A12/,
     .    15X,'RESIDUAL FLUXES/INCRMT. FLUXES =',F10.4,' ( % )',2A12/,
     .    15X,'ITERATV. TEMPE./INCRMT. TEMPE. =',F10.4,' ( % )',2A12/,
     .    15X,'ITERATV. ENERGY/TOTAL   ENERGY =',F10.4,' ( % )',2A12/)
C
c 800 FORMAT(/10X,'CONVERGENCE CHECKS :',//,
c    .    15X,'RESIDUAL FLUXES/TOTAL   FLUXES =',F10.4,' ( % )',2A12/,
c    .    15X,'RESIDUAL FLUXES/INCRMT. FLUXES =',F10.4,' ( % )',2A12/,
c    .    15X,'ITERATV. TEMPE./INCRMT. TEMPE. =',F10.4,' ( % )',2A12/,
c    .    15X,'ITERATV. ENERGY/TOTAL   ENERGY =',F10.4,' ( % )',2A12/,
c    .    15X,'ITERATV. TEMPE./INC111. TEMPE. =',F10.4,' ( % )',2A12/)
c
  901 FORMAT(10X,'>>>>>>>>> CONVERGED ',10X,'CPU TIME =',F10.3)
  902 FORMAT(10X,'>>>>>>>>> CONVERGING',10X,'CPU TIME =',F10.3)
  903 FORMAT(10X,'>>>>>>>>> DIVERGING ',10X,'CPU TIME =',F10.3)
  904 FORMAT(10X,'>>>>>>>>> SKIPPING  ',10X,'CPU TIME =',F10.3)
  905 FORMAT(10X,'>>>>>>>>> JUMPING TO NEXT STEP')
  910 FORMAT(5X,'ITIME =',I5,2X,'ISTEP =',I5,2X,'IITER =',I5,2X,/,
     .      10X,'CUCPU =',F10.3,2X,'RATIO/TOLER =',E15.6/)
      END
