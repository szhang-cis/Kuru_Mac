      SUBROUTINE PLOOUTT(DISTOT,ELVART,TLOADT)
C***********************************************************************
C
C**** THIS ROUTINE OUTPUTS RESULTS FOR X-Y CURVES
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inpo_omt.f'
C
      DIMENSION ELVART(NSTATT), DISTOT(NTOTVT,3), TLOADT(NTOTVT)
      DIMENSION BUFFE(2), PRINC(3)
C
      CHARACTER FORM1*13,NLINE*1
      CHARACTER GNUPLOT*1
C
      DATA GNUPLOT/'#'/
C
      NLINE=CHAR(10)
      FORM1='('//FORMAT//',A1)'     ! Format with return
      NRECHT=  2                    ! Head record with NPONT
C
      DO ICURVT=1,NCURVT            ! Loop on curves of points in curve
       NPONTT(ICURVT+NCOLDT)=NPONTT(ICURVT+NCOLDT)+1
       NRECLT=5+NPONTT(ICURVT+NCOLDT)                 ! Define record
       DO IAXIST=1,2                                  ! Loop on axes
        ICOMDT=MPLOTT(ICURVT,IAXIST,1)         ! Component & command
        ICOMPT=ICOMDT/100                      ! Define component
        ICOMDT=ICOMDT-ICOMPT*100               ! Define command
        IPOI1T=MPLOTT(ICURVT,IAXIST,2)         ! Define node1
        IPOI2T=MPLOTT(ICURVT,IAXIST,3)         ! Define node2
        ITOT1T=0
        ITOT2T=0
        IF(IPOI1T.NE.0) ITOT1T=(IPOI1T-1)*NDOFCT+ICOMPT    ! Define dof1
        IF(IPOI2T.NE.0) ITOT2T=(IPOI2T-1)*NDOFCT+ICOMPT    ! Define dof2
         IELEMT=IPOI1T                               ! Define element
         IGAUST=IPOI2T                               ! Define Gauss P.
         BVALUT=0.0
C
         GOTO(1,2,3,4,5,6,7,8,9,10,11), ICOMDT          ! Do command
C
    1    AVALUT=TTIMET                                   ! Time
         GO TO 100
    2    AVALUT=ISTEPT                                   ! Nstep
         GO TO 100
    3    AVALUT=TFACTT                                   ! Lambda
         GO TO 100
    4    AVALUT=DISTOT(ITOT1T,1)                         ! Temperature
         IF(ITOT2T.NE.0) BVALUT=DISTOT(ITOT2T,1)
         GO TO 100
    5    AVALUT=TLOADT(ITOT1T)                           ! Force
         IF(ITOT2T.NE.0) BVALUT=TLOADT(ITOT2T)
         GO TO 100
    6    AVALUT=DISTOT(ITOT1T,1)                         ! Displacement
         IF(ITOT2T.NE.0) BVALUT=DISTOT(ITOT2T,1)
         GO TO 100
    7    AVALUT=DISTOT(ITOT1T,2)                         ! Velocity
         IF(ITOT2T.NE.0) BVALUT=DISTOT(ITOT2T,2)
         GO TO 100
    8    AVALUT=DISTOT(ITOT1T,3)                         ! Acceleration
         IF(ITOT2T.NE.0) BVALUT=DISTOT(ITOT2T,3)
         GO TO 100
    9    CONTINUE
         IF(NMEMO4.EQ.0)
     .    CALL RUNENDT('PLOOUTT: ERROR IN WRITING HEAT FLUX')
         CALL DATBAST(ELVART,    3,    2)                ! Stress
         ISTAST=ISTATT(4)+(IGAUST-1)*NSTR1T
         IF(ICOMPT.LE.6) THEN                            !  - cartesian
          AVALUT=ELVART(ISTAST-1+ICOMPT)
         ELSE                                            !  - principal
          ICOMPT=ICOMPT-6
c         CALL PRIVAL(NSTR1T,ELVART(ISTAST),PRINC)
          AVALUT=PRINC(ICOMPT)
         ENDIF
         GO TO 100
   10    CONTINUE
         IF(NMEMO4.EQ.0)
     .    CALL RUNENDT('PLOOUTT: ERROR IN WRITING TEMP. GRAD.')
         CALL DATBAST(ELVART,    3,    2)                ! Strain
         ISTAST=ISTATT(3)+(IGAUST-1)*NSTR1T
         IF(ICOMPT.LE.6) THEN                            !  - cartesian
         AVALUT=ELVART(ISTAST-1+ICOMPT)
        ELSE                                            !  - principal
         ICOMPT=ICOMPT-6
         CALL PRISTN(NSTR1T,ELVART(ISTAST),PRINC)
         AVALUT=PRINC(ICOMPT)
        ENDIF
        GO TO 100
   11   CONTINUE
         IF(NMEMO3.EQ.0)
     .    CALL RUNENDT('PLOOUTT: ERROR IN WRITING INTERNAL VAR.')
        CALL DATBAST(ELVART,    3,    2)                ! Inter. Var.
        AVALUT=ELVART(ISTATT(2)+(IGAUST-1)*NHISTT-1+ICOMPT)
C
  100   BUFFE(IAXIST)=AVALUT-BVALUT
       ENDDO ! IAXIS
C
       NUNITT=LUCU1T-1+ICURVT+NCOLDT                 ! Define file
       WRITE(UNIT=NUNITT,REC=NRECLT,FMT=FORM1)       ! Write record
     .       BUFFE(1),BUFFE(2),NLINE                 ! Write head record
C
       IGNUPLOT=1                                    ! better as input
       IF(IGNUPLOT.EQ.0) THEN
        WRITE(UNIT=NUNITT,REC=NRECHT,FMT='(2I5,30X,A1)')
     .                        NPONTT(ICURVT+NCOLDT),0,NLINE
       ELSE
        WRITE(UNIT=NUNITT,REC=NRECHT,FMT='(A1,2I5,29X,A1)')
     .                GNUPLOT,NPONTT(ICURVT+NCOLDT),0,NLINE
       ENDIF
      ENDDO ! ICURV
C
      RETURN
      END
