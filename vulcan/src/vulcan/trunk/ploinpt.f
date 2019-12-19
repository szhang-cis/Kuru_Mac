      SUBROUTINE PLOINPT(NFILET,IFLAGT)
C***********************************************************************
C
C**** THIS ROUTINE READS AND INTERPRETS THE "PLOTTER" COMMANDS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** VAR
C
      INCLUDE 'inpo_omt.f'
      INCLUDE 'prob_omt.f'
C
      PARAMETER (MCOMOT=13,MCOMMT=11,IRECLET=41)
C
      CHARACTER CNAMET*80,NLINET*1,COMMDT(MCOMMT)*5,COMOST(MCOMOT)*2
C
      DATA FORMAT/'(2E20.8)'/
      DATA COMMDT/'TIME', 'NSTEP','LAMBD','TEMPE',
     .            'FORCE','DISPL','VELOC','ACCEL',
     .            'STRES','STRAI','INTER'/
      DATA COMOST/'X','Y','Z','P',
     .            'XX','YY','XY','ZZ','XZ','YZ',
     .            'P1','P2','P3'/
C
C**** BEGIN
C
      NLINET=CHAR(10)
   10 READ(NFILET,'(A80)',ERR=100,END=101) CNAMET
C                                             ! Begin loop reading TITLE
      IFIRSTT=1                               ! Drop initial blanks
      DO WHILE (CNAMET(IFIRSTT:IFIRSTT).EQ.' ')
       IFIRSTT=IFIRSTT+1
      ENDDO
      CNAMET=CNAMET(IFIRSTT:80)
      IF((CNAMET(1:1).EQ.'$').OR.           ! Jump comment & blank lines
     .   (CNAMET(1:1).EQ.'!').OR.
     .   (CNAMET(1:1).EQ.' ')    ) GO TO 10
      IF(CNAMET(1:8).EQ.'END_PLOT') GO TO 20     ! Go to end of loop
      IF(IFLAGT.EQ.1) THEN                       ! Save old curves
       NCOLDT=NCOLDT+NCURVT
       NCURVT=0
      ENDIF
      NCURVT=NCURVT+1                            ! Start with new curve
      IF(NCURVT.GT.MMCURT) GOTO 103              ! Too many curves
      NUNITT=NCURVT+NCOLDT+LUCU1T-1
C
      npichu=nunitt-(lucu1t-1)
      go to (31,32,33,34,35,36,37,38,39,40), npichu
   31 continue
      OPEN(UNIT=NUNITT,FILE=C1T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   32 continue
      OPEN(UNIT=NUNITT,FILE=C2T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   33 continue
      OPEN(UNIT=NUNITT,FILE=C3T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   34 continue
      OPEN(UNIT=NUNITT,FILE=C4T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   35 continue
      OPEN(UNIT=NUNITT,FILE=C5T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   36 continue
      OPEN(UNIT=NUNITT,FILE=C6T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   37 continue
      OPEN(UNIT=NUNITT,FILE=C7T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   38 continue
      OPEN(UNIT=NUNITT,FILE=C8T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   39 continue
      OPEN(UNIT=NUNITT,FILE=C9T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   40 continue
      OPEN(UNIT=NUNITT,FILE=C10T,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLET,ACCESS='DIRECT')
      go to 11
   11 continue
C
      IGNUPLOT=1                                 ! better as input
      IF(IGNUPLOT.EQ.0) THEN
       WRITE(UNIT=NUNITT,REC=1,FMT=901) FORMAT,NLINET
       WRITE(UNIT=NUNITT,REC=2,FMT=902) NLINET   !Save space for npoints
       WRITE(UNIT=NUNITT,REC=3,FMT=903) CNAMET,NLINET
      ELSE
       WRITE(UNIT=NUNITT,REC=1,FMT=1901) FORMAT,NLINET
       WRITE(UNIT=NUNITT,REC=2,FMT=1902) NLINET  !Save space for npoints
       WRITE(UNIT=NUNITT,REC=3,FMT=1903) CNAMET,NLINET
      ENDIF
C
      DO IAXIST=1,2
       NNRECT=IAXIST+3
       NPRINT=0
       ITAPET=LUDATT
       CALL LISTENT('PLOINPT',NPRINT,ITAPET)
       IF((IAXIST.EQ.1.AND.WORDST(1).NE.'X').OR.
     .    (IAXIST.EQ.2.AND.WORDST(1).NE.'Y')    ) GOTO 102 ! Error
       NPAR1T=0
       NPAR2T=0
       NPAR3T=0
       IF(NNPART.GE.1) NPAR1T=INT(PARAMT(1)) ! The paramet. are integers
       IF(NNPART.GE.2) NPAR2T=INT(PARAMT(2))
       IF(NNPART.GE.3) NPAR3T=INT(PARAMT(3))
       NCOMOT=0
       IF(WORDST(3).NE.'' ) THEN             ! Identify component
        DO ICOMOT=1,MCOMOT
         IF (WORDST(3).EQ.COMOST(ICOMOT)) NCOMOT=ICOMOT
        ENDDO
        IF(NCOMOT.EQ.0) GOTO 102 ! Error
       ENDIF
       NCOMMT=0                              ! Identify command
       ILISTT=0
       DO WHILE ((ILISTT.LT.MCOMMT).AND.(NCOMMT.EQ.0))
        ILISTT=ILISTT+1
        IF (WORDST(2).EQ.COMMDT(ILISTT)) NCOMMT=ILISTT
       ENDDO
       IF(NCOMMT.EQ.0) GOTO 102 ! Error
C
       IF(NCOMMT.LE.3) THEN                 ! TIME,NSTEP or LAMBDA
        IF(NNPART.NE.0) GOTO 102 ! Error
        MPLOTT(NCURVT,IAXIST,1)=NCOMMT      ! COMMAND
        IF(IGNUPLOT.EQ.0) THEN
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=910) WORDST(2),NLINET
        ELSE
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1910) WORDST(2),NLINET
        ENDIF
       ENDIF
C
       IF(NCOMMT.EQ.4) THEN                 ! TEMPE
        IF((NNPART.LT.1).OR.(NNPART.GT.2)) GOTO 102 ! Error
        NCOMOT=1
        MPLOTT(NCURVT,IAXIST,1)=NCOMOT*100+NCOMMT  ! COMPONENT & COMMAND
        MPLOTT(NCURVT,IAXIST,2)=NPAR1T             ! IPOIN1
        MPLOTT(NCURVT,IAXIST,3)=0                  ! IPOIN2
        IF(NNPART.EQ.1) THEN
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=920)
     .     WORDST(2),WORDST(3),NPAR1T,NLINET
         ELSE
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1920)
     .     WORDST(2),WORDST(3),NPAR1T,NLINET
         ENDIF
        ELSE ! NNPART=2  =>  Relative force,displ,veloc or accel
         MPLOTT(NCURVT,IAXIST,3)=NPAR2T            ! IPOIN2
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=930)
     .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
         ELSE
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1930)
     .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
         ENDIF
        ENDIF
       ENDIF
C
       IF((NCOMMT.GE.5).AND.(NCOMMT.LE.8)) THEN
C                                           ! FORCE,DISPL,VELOC or ACCEL
        IF((NNPART.LT.1).OR.(NNPART.GT.2)) GOTO 102 ! Error
        MPLOTT(NCURVT,IAXIST,1)=NCOMOT*100+NCOMMT ! COMPONENT & COMMAND
        MPLOTT(NCURVT,IAXIST,2)=NPAR1T            ! IPOIN1
        MPLOTT(NCURVT,IAXIST,3)=0                 ! IPOIN2
        IF(NNPART.EQ.1) THEN
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=920)
     .     WORDST(2),WORDST(3),NPAR1T,NLINET
         ELSE
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1920)
     .     WORDST(2),WORDST(3),NPAR1T,NLINET
         ENDIF
        ELSE ! NNPART=2  =>  Relative force,displ,veloc or accel
         MPLOTT(NCURVT,IAXIST,3)=NPAR2T           ! IPOIN2
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=930)
     .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
         ELSE
          WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1930)
     .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
         ENDIF
        ENDIF
       ENDIF
C
       IF((NCOMMT.EQ.9).OR.(NCOMMT.EQ.10)) THEN   ! STRESS or STRAIN
        NCOMOT=NCOMOT-4                           ! Jump X,Y,Z & P
        IF((NNPART.NE.2).OR.(NCOMOT.LT.1)) GOTO 102 ! Error
        MPLOTT(NCURVT,IAXIST,1)=NCOMOT*100+NCOMMT ! COMPONENT & COMMAND
        MPLOTT(NCURVT,IAXIST,2)=NPAR1T            ! IELEM
        MPLOTT(NCURVT,IAXIST,3)=NPAR2T            ! IGAUS
        IF(IGNUPLOT.EQ.0) THEN
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=940)
     .    WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
        ELSE
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1940)
     .    WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
        ENDIF
       ENDIF
C
       IF(NCOMMT.EQ.11) THEN                      ! INTERNAL VARIABLE
        IF(NNPART.NE.3) GOTO 102 ! Error
        MPLOTT(NCURVT,IAXIST,1)=NPAR1T*100+NCOMMT ! #VARIABLE & COMMAND
        MPLOTT(NCURVT,IAXIST,2)=NPAR2T            ! IELEM
        MPLOTT(NCURVT,IAXIST,3)=NPAR3T            ! IGAUS
        IF(IGNUPLOT.EQ.0) THEN
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=950)
     .    WORDST(2),NPAR1T,NPAR2T,NPAR3T,NLINET
        ELSE
         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1950)
     .    WORDST(2),NPAR1T,NPAR2T,NPAR3T,NLINET
        ENDIF
       ENDIF
      ENDDO
      GO TO 10 ! Loop until find 'END_PLOT'
   20 RETURN
C
  100 CALL RUNENDT('PLOINPT:ERROR WHEN READING TITLE   ')
  101 CALL RUNENDT('PLOINPT:END OF FILE DETECTED       ')
  102 CALL RUNENDT('PLOINPT:ERRORWHEN READING AXIS DATA')
  103 CALL RUNENDT('PLOINPT:MORE THAN 40 CURVES INPUT  ')
C
C**** THE INCLUSION OF '#' IS FOR "gnuplot" POSPROCESSING PROGRAM
C
  901 FORMAT(4X,'1 ',A8,26X,A1)                  !    1 (2E20.8)
  902 FORMAT(4X,'0',4X,'0',30X,A1)               !    0    0
  903 FORMAT(A40,A1)                             !title in 40 characters
  910 FORMAT(A5,35X,A1)                               !NSTEP
  920 FORMAT(A5,'[',A1,'],N',I4,26X,A1)               !FORCE[X],N  32
  930 FORMAT(A5,'[',A1,'],N',I4,' -',I4,20X,A1)       !TEMPE,N  32 -  45
  940 FORMAT(A5,'[',A2,'],E',I4,', P',I2,20X,A1)     !STRES[XX],E 45,P 6
  950 FORMAT(A5,'[',I2,'],E',I4,', P',I2,20X,A1)     !INTER[21],E 89,P 8
C
 1901 FORMAT('#',4X,'1 ',A8,25X,A1)              !    1 (2E20.8)
 1902 FORMAT('#',4X,'0',4X,'0',29X,A1)           !    0    0
 1903 FORMAT('#',A39,A1)                         !title in 40 characters
 1910 FORMAT('#',A5,34X,A1)                           !NSTEP
c 1920 FORMAT('#',A5,'[',A1,'],N',I4,25X,A1)           !FORCE[X],N  32
 1920 FORMAT('#',A5,'[',A1,'],N',I6,23X,A1)           !FORCE[X],N  32
 1930 FORMAT('#',A5,'[',A1,'],N',I4,' -',I4,19X,A1)   !TEMPE,N  32 -  45
 1940 FORMAT('#',A5,'[',A2,'],E',I4,', P',I2,19X,A1) !STRES[XX],E 45,P 6
 1950 FORMAT('#',A5,'[',I2,'],E',I4,', P',I2,19X,A1) !INTER[21],E 89,P 8
      END
