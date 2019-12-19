      SUBROUTINE PLOOUTS(DISTOS,ELVARS,TLOADS)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS RESULTS FOR X-Y CURVES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C
      DIMENSION ELVARS(NSTATS), DISTOS(NTOTVS,3), TLOADS(NTOTVS)
      DIMENSION BUFFE(2), PRINC(3)
C
      CHARACTER FORM1*13,NLINE*1
      CHARACTER GNUPLOT*1
C
      DATA GNUPLOT/'#'/
C
      NLINE=CHAR(10)
      FORM1='('//FORMAS//',A1)'     ! Format with return
      NRECHS=  2                    ! Head record with NPONT
C
      DO ICURVS=1,NCURVS            ! Loop on curves of points in curve
       NPONTS(ICURVS+NCOLDS)=NPONTS(ICURVS+NCOLDS)+1
       NRECLS=5+NPONTS(ICURVS+NCOLDS)                 ! Define record
       DO IAXISS=1,2                                  ! Loop on axes
        ICOMDS=MPLOTS(ICURVS,IAXISS,1)         ! Component & command
        ICOMPS=ICOMDS/100                      ! Define component
        ICOMDS=ICOMDS-ICOMPS*100               ! Define command
        IPOI1S=MPLOTS(ICURVS,IAXISS,2)         ! Define node1
        IPOI2S=MPLOTS(ICURVS,IAXISS,3)         ! Define node2
        ITOT1S=0
        ITOT2S=0
        IF(IPOI1S.NE.0) ITOT1S=(IPOI1S-1)*NDOFCS+ICOMPS    ! Define dof1
        IF(IPOI2S.NE.0) ITOT2S=(IPOI2S-1)*NDOFCS+ICOMPS    ! Define dof2
         IELEMS=IPOI1S                               ! Define element
         IGAUSS=IPOI2S                               ! Define Gauss P.
         BVALUS=0.0
C
         GOTO(1,2,3,4,5,6,7,8,9,10,11), ICOMDS          ! Do command
C
    1    AVALUS=TTIMES                                   ! Time
         GO TO 100
    2    AVALUS=ISTEPS                                   ! Nstep
         GO TO 100
    3    AVALUS=TFACTS                                   ! Lambda
         GO TO 100
    4    AVALUS=DISTOS(ITOT1S,1)                         ! Temperature
         IF(ITOT2S.NE.0) BVALUS=DISTOS(ITOT2S,1)
         GO TO 100
    5    AVALUS=TLOADS(ITOT1S)                           ! Force
         IF(ITOT2S.NE.0) BVALUS=TLOADS(ITOT2S)
         GO TO 100
    6    AVALUS=DISTOS(ITOT1S,1)                         ! Displacement
         IF(ITOT2S.NE.0) BVALUS=DISTOS(ITOT2S,1)
         GO TO 100
    7    AVALUS=DISTOS(ITOT1S,2)                         ! Velocity
         IF(ITOT2S.NE.0) BVALUS=DISTOS(ITOT2S,2)
         GO TO 100
    8    AVALUS=DISTOS(ITOT1S,3)                         ! Acceleration
         IF(ITOT2S.NE.0) BVALUS=DISTOS(ITOT2S,3)
         GO TO 100
    9    CONTINUE
         IF(NMEMO4S.EQ.0)
     .    CALL RUNENDS('PLOOUTS: ERROR IN WRITING HEAT FLUX')
         CALL DATBASS(ELVARS,    3,    2)                ! Stress
         ISTASS=ISTATS(4)+(IGAUSS-1)*NSTR1S
         IF(ICOMPS.LE.6) THEN                            !  - cartesian
          AVALUS=ELVARS(ISTASS-1+ICOMPS)
         ELSE                                            !  - principal
          ICOMPS=ICOMPS-6
c         CALL PRIVAL(NSTR1T,ELVART(ISTAST),PRINC)
          AVALUS=PRINC(ICOMPS)
         ENDIF
         GO TO 100
   10    CONTINUE
         IF(NMEMO4S.EQ.0)
     .    CALL RUNENDS('PLOOUTS: ERROR IN WRITING TEMP. GRAD.')
         CALL DATBASS(ELVARS,    3,    2)                ! Strain
         ISTASS=ISTATS(3)+(IGAUSS-1)*NSTR1S
         IF(ICOMPS.LE.6) THEN                            !  - cartesian
         AVALUS=ELVARS(ISTASS-1+ICOMPS)
        ELSE                                            !  - principal
         ICOMPS=ICOMPS-6
         CALL PRISTN(NSTR1S,ELVARS(ISTASS),PRINC)
         AVALUS=PRINC(ICOMPS)
        ENDIF
        GO TO 100
   11   CONTINUE
         IF(NMEMO3S.EQ.0)
     .    CALL RUNENDS('PLOOUTS: ERROR IN WRITING INTERNAL VAR.')
        CALL DATBASS(ELVARS,    3,    2)                ! Inter. Var.
        AVALUS=ELVARS(ISTATS(2)+(IGAUSS-1)*NHISTS-1+ICOMPS)
C
  100   BUFFE(IAXISS)=AVALUS-BVALUS
       ENDDO ! IAXIS
C
       NUNITS=LUCU1S-1+ICURVS+NCOLDS                 ! Define file
       WRITE(UNIT=NUNITS,REC=NRECLS,FMT=FORM1)       ! Write record
     .       BUFFE(1),BUFFE(2),NLINE                 ! Write head record
C
       IGNUPLOT=1                                    ! better as input
       IF(IGNUPLOT.EQ.0) THEN
        WRITE(UNIT=NUNITS,REC=NRECHS,FMT='(2I5,30X,A1)')
     .                        NPONTS(ICURVS+NCOLDS),0,NLINE
       ELSE
        WRITE(UNIT=NUNITS,REC=NRECHS,FMT='(A1,2I5,29X,A1)')
     .                GNUPLOT,NPONTS(ICURVS+NCOLDS),0,NLINE
       ENDIF
      ENDDO ! ICURV
C
      RETURN
      END
