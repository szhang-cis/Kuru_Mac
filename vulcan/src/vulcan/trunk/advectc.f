      SUBROUTINE ADVECTC(ADVELT,FPCHAT,IPRINX)
C***********************************************************************
C
C**** THIS ROUTINE READS AND SETS DATA FOR ADVECTIVE VELOCITY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'             ! thermal-mechanical
      INCLUDE 'nued_om.f'             ! thermal-microstructural
      INCLUDE 'nuef_om.f'             ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
      INCLUDE 'prob_omt.f'
C
      COMMON/LDFILET/ITAPET
C
      DIMENSION ADVELT(NDIMET,NTOTVT), FPCHAT(NFPCH,NPOINT)
C
      IERR1T=1
      IF(ITAPET.NE.LUDATT.AND.ITAPET.NE.LUADVT) GOTO 1000
      IF(ITAPET.EQ.LUADVT) THEN       ! external .adv file
       IF(IOADVT.EQ.0) OPEN(UNIT=LUADVT,FILE=CG1T,STATUS='OLD',ERR=1000)
       IOADVT=1
       IERR1T=0
C
 1000  IF(IERR1T.NE.0) THEN
        IF(IERR1T.EQ.1) WRITE(LUREST,901)
        CALL RUNENDT('ERROR IN OPENING FILE ')
       ENDIF
      ENDIF      ! itapet.eq.luadvt
C
      WRITE(LUREST,900)
C
C**** READ THE ADVECTIVE VELOCITY
C
      NPRINT=1
      CALL LISTENT('ADVECTC',NPRINT,ITAPET)
      NPOI1T=INT(PARAMT(1))
      WRITE(LUREST,905) NPOI1T
      IF(NPOINT.LT.NPOI1T) THEN
       WRITE(LUREST,910)
       CALL RUNENDT('ADVECTC: ERROR WHEN READING THE NUMBER OF POINTS
     .                        WITH ADVECTIVE VELOCITY')
      ENDIF
C
C**** READ NODAL VARIABLES
C
      IAUXX=0
      IF(ITERMEF.EQ.0) THEN
       IF(NFILL.EQ.1) THEN
        IAUXX=1
        IF(IMICR.EQ.0) THEN
         IPSEU=2*NNUPT+1
        ELSE
         IPSEU=2*NNUPT+NNUPO+1
        ENDIF
       ENDIF
      ENDIF
C
      DO IPOINT=1,NPOI1T
       CALL LISTENT('ADVECTC',NPRINT,ITAPET)
       JPOINT=INT(PARAMT(1))
       DO IDIMET=1,NDIMET
        ADVELT(IDIMET,JPOINT)=PARAMT(1+IDIMET)
       ENDDO
       IF(IAUXX.EQ.1) FPCHAT(IPSEU,IPOINT)=PARAMT(1+NDIMET+1)   ! pseudo
      ENDDO                ! ipoint=1,npoint
C
C**** ECHO NODAL VALUES FOR THE WHOLE MESH
C
      IF(IPRINX.EQ.0) THEN
       IF(IAUXX.EQ.0) THEN
        WRITE(LUREST,921)
        DO IPOINT=1,NPOINT
         WRITE(LUREST,931)
     .                   IPOINT,(ADVELT(IDIMET,IPOINT),IDIMET=1,NDIMET)
        ENDDO
       ENDIF
       IF(IAUXX.EQ.1) THEN
        WRITE(LUREST,922)
        DO IPOINT=1,NPOINT
         WRITE(LUREST,932)
     .                   IPOINT,(ADVELT(IDIMET,IPOINT),IDIMET=1,NDIMET),
     .                   FPCHAT(IPSEU,IPOINT)
        ENDDO
       ENDIF
      ENDIF                ! iprinx.eq.0
C
      WRITE(LUREST,940)
C
      RETURN
C
  900 FORMAT('ADVECTIVE VELOCITY DATA',//)
  901 FORMAT(' ERROR IN OPENING ADVECTIVE VELOCITY INPUT FILE 246',
     . '(56 LINUX)')
  905 FORMAT(/,5X,'NUMBER OF NODES WITH ADVECTIVE VELOCITY DATA',I5,/)
  910 FORMAT(//,'ERROR: NUMBER OF INPUT POINTS GT TOTAL NO. POINTS')
  921 FORMAT(1X,'NODE',4X,'ADVECTIVE VELOCITY',/)
  922 FORMAT(1X,'NODE',4X,'ADVECTIVE VELOCITY',5X,'MATERIAL INDEX',/)
  931 FORMAT(I5,10X,3(E15.6,5X))
  932 FORMAT(I5,10X,4(E15.6,5X))
  940 FORMAT(/)
      END
