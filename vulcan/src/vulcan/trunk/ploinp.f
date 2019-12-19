      SUBROUTINE PLOINP(NFILE,IFLAG)
C***********************************************************************
C
C**** THIS ROUTINE READS AND INTERPRETS THE "PLOTTER" COMMANDS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'inpo_om.f'
      INCLUDE 'prob_om.f'
C
C**** VAR
C
      PARAMETER (MCOMO=13,MCOMM=10,iRECLE=41)
C
      CHARACTER CNAME*80,NLINE*1,COMMD(MCOMM)*5,COMOS(MCOMO)*2
C
      DATA FORMA/'(2E20.8)'/
      DATA COMMD/'TIME','NSTEP','LAMBD',
     .           'FORCE','DISPL','VELOC','ACCEL',
     .           'STRES','STRAI','INTER'/
      DATA COMOS/'X','Y','Z','P',
     .           'XX','YY','XY','ZZ','XZ','YZ',
     .           'P1','P2','P3'/
C
C**** BEGIN
C
      NLINE=CHAR(10)
   10 READ(NFILE,'(A80)',ERR=100,END=101) CNAME
C                                           ! Begin loop reading TITLE
      IFIRST=1                              ! Drop initial blanks
      DO WHILE (CNAME(IFIRST:IFIRST).EQ.' ')
       IFIRST=IFIRST+1
      ENDDO
      CNAME=CNAME(IFIRST:80)
      IF((CNAME(1:1).EQ.'$').OR.            ! Jump comment & blank lines
     .   (CNAME(1:1).EQ.'!').OR.
     .   (CNAME(1:1).EQ.' ')    ) GO TO 10
      IF(CNAME(1:8).EQ.'END_PLOT') GO TO 20     ! Go to end of loop
      IF(IFLAG.EQ.1) THEN                       ! Save old curves
       NCOLD=NCOLD+NCURV
       NCURV=0
      ENDIF
      NCURV=NCURV+1                             ! Start with new curve
      IF(NCURV.GT.MMCUR) GOTO 103               ! Too many curves
      NUNIT=NCURV+NCOLD+LUCU1-1
C
      npichu=nunit-(lucu1-1)
      go to (31,32,33,34,35,36,37,38,39,40), npichu
   31 continue
      OPEN(UNIT=NUNIT,FILE=C1M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   32 continue
      OPEN(UNIT=NUNIT,FILE=C2M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   33 continue
      OPEN(UNIT=NUNIT,FILE=C3M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   34 continue
      OPEN(UNIT=NUNIT,FILE=C4M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   35 continue
      OPEN(UNIT=NUNIT,FILE=C5M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   36 continue
      OPEN(UNIT=NUNIT,FILE=C6M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   37 continue
      OPEN(UNIT=NUNIT,FILE=C7M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   38 continue
      OPEN(UNIT=NUNIT,FILE=C8M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   39 continue
      OPEN(UNIT=NUNIT,FILE=C9M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   40 continue
      OPEN(UNIT=NUNIT,FILE=C10M,STATUS='NEW',FORM='FORMATTED',
     .     RECL=iRECLE,ACCESS='DIRECT')
      go to 11
   11 continue
C
      IGNUPLO=1                                 ! better as input
      IF(IGNUPLO.EQ.0) THEN
       WRITE(UNIT=NUNIT,REC=1,FMT=901) FORMA,NLINE
       WRITE(UNIT=NUNIT,REC=2,FMT=902) NLINE    ! Save space for npoints
       WRITE(UNIT=NUNIT,REC=3,FMT=903) CNAME,NLINE
      ELSE
       WRITE(UNIT=NUNIT,REC=1,FMT=1901) FORMA,NLINE
       WRITE(UNIT=NUNIT,REC=2,FMT=1902) NLINE   ! Save space for npoints
       WRITE(UNIT=NUNIT,REC=3,FMT=1903) CNAME,NLINE
      ENDIF
      DO IAXIS=1,2
       NNREC=IAXIS+3
       NPRIN=0
       ITAPE=LUDAT
       CALL LISTEN('PLOINP',NPRIN,ITAPE)
       IF((IAXIS.EQ.1.AND.WORDS(1).NE.'X').OR.
     .    (IAXIS.EQ.2.AND.WORDS(1).NE.'Y')    ) GOTO 102 ! Error
       NPAR1=0
       NPAR2=0
       NPAR3=0
       IF(NNPAR.GE.1) NPAR1=INT(PARAM(1))    ! The paramet. are integers
       IF(NNPAR.GE.2) NPAR2=INT(PARAM(2))
       IF(NNPAR.GE.3) NPAR3=INT(PARAM(3))
       NCOMO=0
       IF(WORDS(3).NE.'' ) THEN              ! Identify component
        DO ICOMO=1,MCOMO
         IF (WORDS(3).EQ.COMOS(ICOMO)) NCOMO=ICOMO
        ENDDO
        IF(NCOMO.EQ.0) GOTO 102 ! Error
       ENDIF
       NCOMM=0                               ! Identify command
       ILIST=0
       DO WHILE ((ILIST.LT.MCOMM).AND.(NCOMM.EQ.0))
        ILIST=ILIST+1
        IF (WORDS(2).EQ.COMMD(ILIST)) NCOMM=ILIST
       ENDDO
       IF(NCOMM.EQ.0) GOTO 102  ! Error
       IF(NCOMM.LE.3) THEN                   ! TIME,NSTEP or LAMBDA
        IF(NNPAR.NE.0) GOTO 102 ! Error
        MPLOT(NCURV,IAXIS,1)=NCOMM           ! COMMAND
        IF(IGNUPLO.EQ.0) THEN
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=910) WORDS(2),NLINE
        ELSE
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=1910) WORDS(2),NLINE
        ENDIF
       ENDIF
       IF((NCOMM.GE.4).AND.(NCOMM.LE.7)) THEN
C                                            ! FORCE,DISPL,VELOC,ACCEL
        IF((NNPAR.LT.1).OR.(NNPAR.GT.2)) GOTO 102 ! Error
        MPLOT(NCURV,IAXIS,1)=NCOMO*100+NCOMM ! COMPONENT & COMMAND
        MPLOT(NCURV,IAXIS,2)=NPAR1           ! IPOIN1
        MPLOT(NCURV,IAXIS,3)=0               ! IPOIN2
        IF(NNPAR.EQ.1) THEN
         IF(IGNUPLO.EQ.0) THEN
          WRITE(UNIT=NUNIT,REC=NNREC,FMT=920)
     .     WORDS(2),WORDS(3),NPAR1,NLINE
         ELSE
          WRITE(UNIT=NUNIT,REC=NNREC,FMT=1920)
     .     WORDS(2),WORDS(3),NPAR1,NLINE
         ENDIF
        ELSE ! NNPAR=2  =>  Relative force,displ,veloc or accel
         MPLOT(NCURV,IAXIS,3)=NPAR2          ! IPOIN2
         IF(IGNUPLO.EQ.0) THEN
          WRITE(UNIT=NUNIT,REC=NNREC,FMT=930)
     .     WORDS(2),WORDS(3),NPAR1,NPAR2,NLINE
         ELSE
          WRITE(UNIT=NUNIT,REC=NNREC,FMT=1930)
     .     WORDS(2),WORDS(3),NPAR1,NPAR2,NLINE
         ENDIF
        ENDIF
       ENDIF
       IF((NCOMM.EQ.8).OR.(NCOMM.EQ.9)) THEN ! STRESS or STRAIN
        NCOMO=NCOMO-4                        ! Jump X,Y,Z & P
        IF((NNPAR.NE.2).OR.(NCOMO.LT.1)) GOTO 102 ! Error
        MPLOT(NCURV,IAXIS,1)=NCOMO*100+NCOMM ! COMPONENT & COMMAND
        MPLOT(NCURV,IAXIS,2)=NPAR1           ! IELEM
        MPLOT(NCURV,IAXIS,3)=NPAR2           ! IGAUS
        IF(IGNUPLO.EQ.0) THEN
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=940)
     .    WORDS(2),WORDS(3),NPAR1,NPAR2,NLINE
        ELSE
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=1940)
     .    WORDS(2),WORDS(3),NPAR1,NPAR2,NLINE
        ENDIF
       ENDIF
       IF(NCOMM.EQ.10) THEN                  ! INTERNAL VARIABLE
        IF(NNPAR.NE.3) GOTO 102 ! Error
        MPLOT(NCURV,IAXIS,1)=NPAR1*100+NCOMM ! #VARIABLE & COMMAND
        MPLOT(NCURV,IAXIS,2)=NPAR2           ! IELEM
        MPLOT(NCURV,IAXIS,3)=NPAR3           ! IGAUS
        IF(IGNUPLO.EQ.0) THEN
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=950)
     .    WORDS(2),NPAR1,NPAR2,NPAR3,NLINE
        ELSE
         WRITE(UNIT=NUNIT,REC=NNREC,FMT=1950)
     .    WORDS(2),NPAR1,NPAR2,NPAR3,NLINE
        ENDIF
       ENDIF
      ENDDO
      GO TO 10 ! Loop until find 'END_PLOT'
   20 RETURN
C
  100 CALL RUNEND('PLOINP: ERROR WHEN READING TITLE   ')
  101 CALL RUNEND('PLOINP: END OF FILE DETECTED       ')
  102 CALL RUNEND('PLOINP:ERROR WHEN READING AXIS DATA')
  103 CALL RUNEND('PLOINP: MORE THAN 40 CURVES INPUT  ')
C
C**** THE INCLUSION OF '#' IS FOR "gnuplot" POSPROCESSING PROGRAM
C
  901 FORMAT(4X,'1 ',A8,26X,A1)                   !    1 (2E20.8)
  902 FORMAT(4X,'0',4X,'0',30X,A1)                !    0    0
  903 FORMAT(A40,A1)                              ! title in 40 charac.
  910 FORMAT(A5,35X,A1)                           !NSTEP
  920 FORMAT(A5,'[',A1,'],N',I4,26X,A1)           !FORCE[X],N 32
  930 FORMAT(A5,'[',A1,'],N',I4,' -',I4,20X,A1)   !DISPL[Y],N 32 - 45
  940 FORMAT(A5,'[',A2,'],E',I4,', P',I2,20X,A1)  !STRES[XX],E 45, P 6
  950 FORMAT(A5,'[',I2,'],E',I4,', P',I2,20X,A1)  !INTER[21],E 89, P 8
C
 1901 FORMAT('#',4X,'1 ',A8,25X,A1)                  !    1 (2E20.8)
 1902 FORMAT('#',4X,'0',4X,'0',29X,A1)               !    0    0
 1903 FORMAT('#',A39,A1)                             !title in 40 char.
 1910 FORMAT('#',A5,34X,A1)                          !NSTEP
 1920 FORMAT('#',A5,'[',A1,'],N',I6,23X,A1)          !FORCE[X],N 41
 1930 FORMAT('#',A5,'[',A1,'],N',I6,' -',I4,17X,A1)  !DISPL[Y],N 41 - 55
 1940 FORMAT('#',A5,'[',A2,'],E',I6,', P',I2,17X,A1) !STRES[XX],E 45,P 6
 1950 FORMAT('#',A5,'[',I2,'],E',I6,', P',I2,17X,A1) !INTER[21],E 89,P 8
      END
