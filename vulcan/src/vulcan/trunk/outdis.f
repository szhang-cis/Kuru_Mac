      SUBROUTINE OUTDIS(DISTO,HEADS,DISIT)
C***********************************************************************
C
C**** THIS ROUTINE WRITES DISPLACEMENTS TO OUTPUT AND POSTPROCESSOR
C     FILES
C
C       KPRI1=0  DO NOT WRITE DISPLACEMENTS TO OUTPUT FILE
C             1  WRITE DISPLACEMENTS TO OUTPUT FILES
C       KFEMV=0  DO NOT WRITE DISPLACEMENTS TO POSTPROCESSOR FILE
C             1  WRITE DISPLACEMENTS TO POSTPROCESSOR FILE
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION DISTO(NTOTV,*), HEADS(NPOIN,*), DISIT(*)
      DIMENSION DIAUX(3)
C
C**** COMPUTES NUMBER OF POINTS WITHOUT CONTACT (AL METHOD)
C
      NPOI1=NPOIN-NPOIC
      NTOT1=NPOI1*NDOFC
C
C**** OUTPUT DISPLACEMENTS
C
      IF(IPRCO.GT.0) THEN
       IF(KPRI1.EQ.1) THEN
        WRITE(LURES,915)
        IF(NDIME.EQ.2.AND.KPORE.EQ.0) WRITE(LURES,920)
        IF(NDIME.EQ.2.AND.KPORE.NE.0) WRITE(LURES,925)
        IF(NDIME.EQ.3.AND.KPORE.EQ.0) WRITE(LURES,930)
        IF(NDIME.EQ.3.AND.KPORE.NE.0) WRITE(LURES,935)
C
        GAMAW=9999999.*GRAVY         !???????????????????
        IF(GAMAW.EQ.0.0) GAMAW=1.0
C
        DO IPOIN=1,NPOI1
         ITOT0=(IPOIN-1)*NDOFC+1
         ITOT1=(IPOIN-1)*NDOFC+NDOFN
         ITOTF=IPOIN*NDOFC
C
         IF(KPORE.EQ.0) TOPWP=0.
         IF(KPORE.EQ.1) TOPWP=HEADS(IPOIN,1)
         IF(KPORE.EQ.2) TOPWP=DISTO(ITOTF,1)
C
         IF(KPORE.NE.0) THEN
          EXPWP=TOPWP-HEADS(IPOIN,1)
          THEAD=HEADS(IPOIN,2)+EXPWP/GAMAW
         ENDIF
C
         IF(KPORE.EQ.0)
     .    WRITE(LURES,905) IPOIN,(DISTO(ITOTV,1),ITOTV=ITOT0,ITOT1)
         IF(KPORE.NE.0)
     .    WRITE(LURES,905) IPOIN,(DISTO(ITOTV,1),ITOTV=ITOT0,ITOT1),
     .                     EXPWP,TOPWP,THEAD
        ENDDO
C
C**** OUTPUT VELOCITIES & ACCELERATIONS
C
        IF(KDYNA.EQ.1) THEN
         WRITE(LURES,1915)
         IF(NDIME.EQ.2.AND.KPORE.EQ.0) WRITE(LURES,1920)
         IF(NDIME.EQ.2.AND.KPORE.NE.0) WRITE(LURES,1925)
         IF(NDIME.EQ.3.AND.KPORE.EQ.0) WRITE(LURES,1930)
         IF(NDIME.EQ.3.AND.KPORE.NE.0) WRITE(LURES,1935)
C
         GAMAW=9999999.*GRAVY         !???????????????????
         IF(GAMAW.EQ.0.0) GAMAW=1.0
C
         DO IPOIN=1,NPOI1
          ITOT0=(IPOIN-1)*NDOFC+1
          ITOT1=(IPOIN-1)*NDOFC+NDOFN
          ITOTF=IPOIN*NDOFC

          IF(KPORE.EQ.0) TOPWP=0.
          IF(KPORE.EQ.1) TOPWP=HEADS(IPOIN,1)
          IF(KPORE.EQ.2) TOPWP=DISTO(ITOTF,1)
C
          IF(KPORE.NE.0) THEN
           EXPWP=TOPWP-HEADS(IPOIN,1)
           THEAD=HEADS(IPOIN,2)+EXPWP/GAMAW
          ENDIF
C
          IF(KPORE.EQ.0)
     .     WRITE(LURES,905) IPOIN,(DISTO(ITOTV,2),ITOTV=ITOT0,ITOT1)
          IF(KPORE.NE.0)
     .     WRITE(LURES,905) IPOIN,(DISTO(ITOTV,2),ITOTV=ITOT0,ITOT1),
     .                      EXPWP,TOPWP,THEAD
         ENDDO
C
         WRITE(LURES,2915)
         IF(NDIME.EQ.2.AND.KPORE.EQ.0) WRITE(LURES,2920)
         IF(NDIME.EQ.2.AND.KPORE.NE.0) WRITE(LURES,2925)
         IF(NDIME.EQ.3.AND.KPORE.EQ.0) WRITE(LURES,2930)
         IF(NDIME.EQ.3.AND.KPORE.NE.0) WRITE(LURES,2935)
C
         GAMAW=9999999.*GRAVY         !???????????????????
         IF(GAMAW.EQ.0.0) GAMAW=1.0
C
         DO IPOIN=1,NPOI1
          ITOT0=(IPOIN-1)*NDOFC+1
          ITOT1=(IPOIN-1)*NDOFC+NDOFN
          ITOTF=IPOIN*NDOFC

          IF(KPORE.EQ.0) TOPWP=0.
          IF(KPORE.EQ.1) TOPWP=HEADS(IPOIN,1)
          IF(KPORE.EQ.2) TOPWP=DISTO(ITOTF,1)
C
          IF(KPORE.NE.0) THEN
           EXPWP=TOPWP-HEADS(IPOIN,1)
           THEAD=HEADS(IPOIN,2)+EXPWP/GAMAW
          ENDIF
C
          IF(KPORE.EQ.0)
     .     WRITE(LURES,905) IPOIN,(DISTO(ITOTV,3),ITOTV=ITOT0,ITOT1)
          IF(KPORE.NE.0)
     .     WRITE(LURES,905) IPOIN,(DISTO(ITOTV,3),ITOTV=ITOT0,ITOT1),
     .                      EXPWP,TOPWP,THEAD
         ENDDO
        ENDIF                 ! kdyna.eq.1
C
C**** PRINTS CONTACT FORCES FOR THE AUGMENTED LAGRANGIAN METHOD
C
        IF(NPOIC.NE.0) THEN
         WRITE(LURES,916)
         IF(NDIME.EQ.2.AND.KPORE.EQ.0) WRITE(LURES,921)
         IF(NDIME.EQ.3.AND.KPORE.EQ.0) WRITE(LURES,931)
C
         DO IPOIN=NPOI1+1,NPOIN
          ITOT0=(IPOIN-1)*NDOFC+1
          ITOT1=(IPOIN-1)*NDOFC+NDOFN
          ITOTF=IPOIN*NDOFC
          IPOIC=IPOIN-NPOI1
          IAUXX=(IPOIN-1)*NDOFC
C
          DO ITOTV=ITOT0,ITOT1                      ! no traction force
           DIAUX(ITOTV-IAUXX)=DISTO(ITOTV,1)
           IF(DIAUX(ITOTV-IAUXX).GT.0.0) DIAUX(ITOTV-IAUXX)=0.0
          ENDDO
C
          WRITE(LURES,905) IPOIC,(DIAUX(ITOTV-IAUXX),ITOTV=ITOT0,ITOT1)
         ENDDO
        ENDIF             ! npoic.ne.0
       ENDIF              ! kpri1.eq.1
      ENDIF               ! iprco.gt.0
C
      IF(KFEMV.NE.0)
     . WRITE(LUPOS) TITLE,SUBTI,ITIME,ISTEP,IITER,SNGL(TTIME),
     .             (SNGL(DISTO(ITOTV,1)),ITOTV=1,NTOT1)
C
      IF(KFEMV.NE.0) THEN
       IF(KDYNA.EQ.1) THEN
        WRITE(LUPOS) (SNGL(DISTO(ITOTV,2)),ITOTV=1,NTOT1)
        WRITE(LUPOS) (SNGL(DISTO(ITOTV,3)),ITOTV=1,NTOT1)
       ENDIF
      ENDIF
C
      RETURN
  905 FORMAT(5X,I5,8E15.6)
  915 FORMAT(/,5X,14H DISPLACEMENTS)
  916 FORMAT(/,5X,15H CONTACT FORCES)
  920 FORMAT(5X,5H NODE,7X,8H X-DISP.,7X,8H Y-DISP.)
  921 FORMAT(5X,5H NODE,7X,8H X-FORC.,7X,8H Y-FORC.)
  925 FORMAT(5X,5H NODE,7X,8H X-DISP.,7X,8H Y-DISP.,7X,8H EXP.PWP,
     .                  7X,8H TOT.PWP,6X,9H TOT.HEAD)
  930 FORMAT(5X,5H NODE,7X,8H X-DISP.,7X,8H Y-DISP.,7X,8H Z-DISP.)
  931 FORMAT(5X,5H NODE,7X,8H X-FORC.,7X,8H Y-FORC.,7X,8H Z-FORC.)
  935 FORMAT(5X,5H NODE,7X,8H X-DISP.,7X,8H Y-DISP.,7X,8H Z-DISP.,
     .                  7X,8H EXP.PWP,7X,8H TOT.PWP,6X,9H TOT.HEAD)
C
 1915 FORMAT(/,5X,11H VELOCITIES)
 1920 FORMAT(5X,5H NODE,7X,7H X-VEL.,7X,7H Y-VEL.)
 1925 FORMAT(5X,5H NODE,7X,7H X-VEL.,7X,7H Y-VEL.,7X,8H EXP.PWP,
     .                  7X,8H TOT.PWP,6X,9H TOT.HEAD)
 1930 FORMAT(5X,5H NODE,7X,7H X-VEL.,7X,7H Y-VEL.,7X,7H Z-VEL.)
 1935 FORMAT(5X,5H NODE,7X,7H X-VEL.,7X,7H Y-VEL.,7X,7H Z-VEL.,
     .                  7X,8H EXP.PWP,7X,8H TOT.PWP,6X,9H TOT.HEAD)
C
 2915 FORMAT(/,5X,14H ACCELERATIONS)
 2920 FORMAT(5X,5H NODE,7X,7H X-ACC.,7X,7H Y-ACC.)
 2925 FORMAT(5X,5H NODE,7X,7H X-ACC.,7X,7H Y-ACC.,7X,8H EXP.PWP,
     .                  7X,8H TOT.PWP,6X,9H TOT.HEAD)
 2930 FORMAT(5X,5H NODE,7X,7H X-ACC.,7X,7H Y-ACC.,7X,7H Z-ACC.)
 2935 FORMAT(5X,5H NODE,7X,7H X-ACC.,7X,7H Y-ACC.,7X,7H Z-ACC.,
     .                  7X,8H EXP.PWP,7X,8H TOT.PWP,6X,9H TOT.HEAD)
C
      END
