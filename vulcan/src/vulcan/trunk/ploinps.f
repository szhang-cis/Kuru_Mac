      SUBROUTINE PLOINPS(NFILES,IFLAGS)
C***********************************************************************
C
C**** THIS ROUTINE READS AND INTERPRETS THE "PLOTTER" COMMANDS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** VAR
C
      INCLUDE 'inpo_oms.f'
      INCLUDE 'prob_oms.f'
C
      PARAMETER (MCOMOS=13,MCOMMS=11,IRECLES=41)
C
      CHARACTER CNAMES*80,NLINES*1,COMMDS(MCOMMS)*5,COMOSS(MCOMOS)*2
C
      DATA FORMAS/'(2E20.8)'/
      DATA COMMDS/'TIME', 'NSTEP','LAMBD','TEMPE',
     .            'FORCE','DISPL','VELOC','ACCEL',
     .            'STRES','STRAI','INTER'/
      DATA COMOSS/'X','Y','Z','P',
     .            'XX','YY','XY','ZZ','XZ','YZ',
     .            'P1','P2','P3'/
C
C**** BEGIN
C
      NLINES=CHAR(10)
   10 READ(NFILES,'(A80)',ERR=100,END=101) CNAMES
C                                             ! Begin loop reading TITLE
      IFIRSTS=1                               ! Drop initial blanks
      DO WHILE (CNAMES(IFIRSTS:IFIRSTS).EQ.' ')
       IFIRSTS=IFIRSTS+1
      ENDDO
      CNAMES=CNAMES(IFIRSTS:80)
      IF((CNAMES(1:1).EQ.'$').OR.           ! Jump comment & blank lines
     .   (CNAMES(1:1).EQ.'!').OR.
     .   (CNAMES(1:1).EQ.' ')    ) GO TO 10
      IF(CNAMES(1:8).EQ.'END_PLOT') GO TO 20     ! Go to end of loop
      IF(IFLAGS.EQ.1) THEN                       ! Save old curves
       NCOLDS=NCOLDS+NCURVS
       NCURVS=0
      ENDIF
      NCURVS=NCURVS+1                            ! Start with new curve
      IF(NCURVS.GT.MMCURS) GOTO 103              ! Too many curves
      NUNITS=NCURVS+NCOLDS+LUCU1S-1
C
      npichu=nunits-(lucu1s-1)
      go to (31,32,33,34,35,36,37,38,39,40), npichu
   31 continue
      OPEN(UNIT=NUNITS,FILE=C1S,STATUS='NEW',FORM='FORMATTED',
     .     RECL=IRECLES,ACCESS='DIRECT')
      go to 11
   32 call runends('unit not available')
c  32 continue
c     OPEN(UNIT=NUNITT,FILE=C2T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   33 call runends('unit not available')
c  33 continue
c     OPEN(UNIT=NUNITT,FILE=C3T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   34 call runends('unit not available')
c  34 continue
c     OPEN(UNIT=NUNITT,FILE=C4T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   35 call runends('unit not available')
c  35 continue
c     OPEN(UNIT=NUNITT,FILE=C5T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   36 call runends('unit not available')
c  36 continue
c     OPEN(UNIT=NUNITT,FILE=C6T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   37 call runends('unit not available')
c  37 continue
c     OPEN(UNIT=NUNITT,FILE=C7T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   38 call runends('unit not available')
c  38 continue
c     OPEN(UNIT=NUNITT,FILE=C8T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   39 call runends('unit not available')
c  39 continue
c     OPEN(UNIT=NUNITT,FILE=C9T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   40 call runends('unit not available')
c  40 continue
c     OPEN(UNIT=NUNITT,FILE=C10T,STATUS='NEW',FORM='FORMATTED',
c    .     RECL=IRECLET,ACCESS='DIRECT')
c     go to 11
   11 continue
C
      IGNUPLOT=1                                 ! better as input
      IF(IGNUPLOT.EQ.0) THEN
       WRITE(UNIT=NUNITS,REC=1,FMT=901) FORMAS,NLINES
       WRITE(UNIT=NUNITS,REC=2,FMT=902) NLINES   !Save space for npoints
       WRITE(UNIT=NUNITS,REC=3,FMT=903) CNAMES,NLINES
      ELSE
       WRITE(UNIT=NUNITS,REC=1,FMT=1901) FORMAS,NLINES
       WRITE(UNIT=NUNITS,REC=2,FMT=1902) NLINES  !Save space for npoints
       WRITE(UNIT=NUNITS,REC=3,FMT=1903) CNAMES,NLINES
      ENDIF
C
      DO IAXISS=1,2
       NNRECS=IAXISS+3
       NPRINS=0
       ITAPES=LUDATS
       CALL LISTENS('PLOINPS',NPRINS,ITAPES)
       IF((IAXISS.EQ.1.AND.WORDSS(1).NE.'X').OR.
     .    (IAXISS.EQ.2.AND.WORDSS(1).NE.'Y')    ) GOTO 102 ! Error
       NPAR1S=0
       NPAR2S=0
       NPAR3S=0
       IF(NNPARS.GE.1) NPAR1S=INT(PARAMS(1)) ! The paramet. are integers
       IF(NNPARS.GE.2) NPAR2S=INT(PARAMS(2))
       IF(NNPARS.GE.3) NPAR3S=INT(PARAMS(3))
       NCOMOS=0
       IF(WORDSS(3).NE.'' ) THEN             ! Identify component
        DO ICOMOS=1,MCOMOS
         IF (WORDSS(3).EQ.COMOSS(ICOMOS)) NCOMOS=ICOMOS
        ENDDO
        IF(NCOMOS.EQ.0) GOTO 102 ! Error
       ENDIF
       NCOMMS=0                              ! Identify command
       ILISTS=0
       DO WHILE ((ILISTS.LT.MCOMMS).AND.(NCOMMS.EQ.0))
        ILISTS=ILISTS+1
        IF (WORDSS(2).EQ.COMMDS(ILISTS)) NCOMMS=ILISTS
       ENDDO
       IF(NCOMMS.EQ.0) GOTO 102 ! Error
C
       IF(NCOMMS.LE.3) THEN                 ! TIME,NSTEP or LAMBDA
        IF(NNPARS.NE.0) GOTO 102 ! Error
        MPLOTS(NCURVS,IAXISS,1)=NCOMMS      ! COMMAND
        IF(IGNUPLOT.EQ.0) THEN
         WRITE(UNIT=NUNITS,REC=NNRECS,FMT=910) WORDSS(2),NLINES
        ELSE
         WRITE(UNIT=NUNITS,REC=NNRECS,FMT=1910) WORDSS(2),NLINES
        ENDIF
       ENDIF
C
       IF(NCOMMS.EQ.4) THEN                 ! TEMPE
        IF((NNPARS.LT.1).OR.(NNPARS.GT.2)) GOTO 102 ! Error
        NCOMOS=1
        MPLOTS(NCURVS,IAXISS,1)=NCOMOS*100+NCOMMS  ! COMPONENT & COMMAND
        MPLOTS(NCURVS,IAXISS,2)=NPAR1S             ! IPOIN1
        MPLOTS(NCURVS,IAXISS,3)=0                  ! IPOIN2
        IF(NNPARS.EQ.1) THEN
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITS,REC=NNRECS,FMT=920)
     .     WORDSS(2),WORDSS(3),NPAR1S,NLINES
         ELSE
          WRITE(UNIT=NUNITS,REC=NNRECS,FMT=1920)
     .     WORDSS(2),WORDSS(3),NPAR1S,NLINES
         ENDIF
        ELSE ! NNPARS=2  =>  Relative force,displ,veloc or accel
         MPLOTS(NCURVS,IAXISS,3)=NPAR2S            ! IPOIN2
         IF(IGNUPLOT.EQ.0) THEN
          WRITE(UNIT=NUNITS,REC=NNRECS,FMT=930)
     .     WORDSS(2),WORDSS(3),NPAR1S,NPAR2S,NLINES
         ELSE
          WRITE(UNIT=NUNITS,REC=NNRECS,FMT=1930)
     .     WORDSS(2),WORDSS(3),NPAR1S,NPAR2S,NLINES
         ENDIF
        ENDIF
       ENDIF
C
c      IF((NCOMMT.GE.5).AND.(NCOMMT.LE.8)) THEN
C                                           ! FORCE,DISPL,VELOC or ACCEL
c       IF((NNPART.LT.1).OR.(NNPART.GT.2)) GOTO 102 ! Error
c       MPLOTT(NCURVT,IAXIST,1)=NCOMOT*100+NCOMMT ! COMPONENT & COMMAND
c       MPLOTT(NCURVT,IAXIST,2)=NPAR1T            ! IPOIN1
c       MPLOTT(NCURVT,IAXIST,3)=0                 ! IPOIN2
c       IF(NNPART.EQ.1) THEN
c        IF(IGNUPLOT.EQ.0) THEN
c         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=920)
c    .     WORDST(2),WORDST(3),NPAR1T,NLINET
c        ELSE
c         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1920)
c    .     WORDST(2),WORDST(3),NPAR1T,NLINET
c        ENDIF
c       ELSE ! NNPART=2  =>  Relative force,displ,veloc or accel
c        MPLOTT(NCURVT,IAXIST,3)=NPAR2T           ! IPOIN2
c        IF(IGNUPLOT.EQ.0) THEN
c         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=930)
c    .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
c        ELSE
c         WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1930)
c    .     WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
c        ENDIF
c       ENDIF
c      ENDIF
C
c      IF((NCOMMT.EQ.9).OR.(NCOMMT.EQ.10)) THEN   ! STRESS or STRAIN
c       NCOMOT=NCOMOT-4                           ! Jump X,Y,Z & P
c       IF((NNPART.NE.2).OR.(NCOMOT.LT.1)) GOTO 102 ! Error
c       MPLOTT(NCURVT,IAXIST,1)=NCOMOT*100+NCOMMT ! COMPONENT & COMMAND
c       MPLOTT(NCURVT,IAXIST,2)=NPAR1T            ! IELEM
c       MPLOTT(NCURVT,IAXIST,3)=NPAR2T            ! IGAUS
c       IF(IGNUPLOT.EQ.0) THEN
c        WRITE(UNIT=NUNITT,REC=NNRECT,FMT=940)
c    .    WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
c       ELSE
c        WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1940)
c    .    WORDST(2),WORDST(3),NPAR1T,NPAR2T,NLINET
c       ENDIF
c      ENDIF
C
c      IF(NCOMMT.EQ.11) THEN                      ! INTERNAL VARIABLE
c       IF(NNPART.NE.3) GOTO 102 ! Error
c       MPLOTT(NCURVT,IAXIST,1)=NPAR1T*100+NCOMMT ! #VARIABLE & COMMAND
c       MPLOTT(NCURVT,IAXIST,2)=NPAR2T            ! IELEM
c       MPLOTT(NCURVT,IAXIST,3)=NPAR3T            ! IGAUS
c       IF(IGNUPLOT.EQ.0) THEN
c        WRITE(UNIT=NUNITT,REC=NNRECT,FMT=950)
c    .    WORDST(2),NPAR1T,NPAR2T,NPAR3T,NLINET
c       ELSE
c        WRITE(UNIT=NUNITT,REC=NNRECT,FMT=1950)
c    .    WORDST(2),NPAR1T,NPAR2T,NPAR3T,NLINET
c       ENDIF
c      ENDIF
      ENDDO
      GO TO 10 ! Loop until find 'END_PLOT'
   20 RETURN
C
  100 CALL RUNENDS('PLOINPT:ERROR WHEN READING TITLE   ')
  101 CALL RUNENDS('PLOINPT:END OF FILE DETECTED       ')
  102 CALL RUNENDS('PLOINPT:ERRORWHEN READING AXIS DATA')
  103 CALL RUNENDS('PLOINPT:MORE THAN 40 CURVES INPUT  ')
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
c1920 FORMAT('#',A5,'[',A1,'],N',I4,25X,A1)           !FORCE[X],N  32
 1920 FORMAT('#',A5,'[',A1,'],N',I6,23X,A1)           !FORCE[X],N  32
 1930 FORMAT('#',A5,'[',A1,'],N',I4,' -',I4,19X,A1)   !TEMPE,N  32 -  45
 1940 FORMAT('#',A5,'[',A2,'],E',I4,', P',I2,19X,A1) !STRES[XX],E 45,P 6
 1950 FORMAT('#',A5,'[',I2,'],E',I4,', P',I2,19X,A1) !INTER[21],E 89,P 8
      END
