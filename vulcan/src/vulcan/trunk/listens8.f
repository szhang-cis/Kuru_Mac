      SUBROUTINE LISTENS8(SUBNAT,NPRINT,ITAPET)
C***********************************************************************
C
C**** READS A STRING AND INTERPRETS IT AS WORDS AND PARAMETERS. (LINUX)
C
C   - MAXIMUM NUMBER OF WORDS AND PARAMETERS = MAXWP
C
C   - ONLY THE FIRST FIVE CHARACTERS OF EACH WORD ARE DECODED.
C
C   - EACH WORD OR PARAMETER MUST BE SEPARETED BY ' ', '=', ':' OR ','
C
C   - A "COMMENT" BEGIN WITH '$', '!', '/'  OR '\' IN ANY PLACE.
C
C   - A LINE THAN HAVE A COMMENT BEGINING WITH '/' OR '\'
C     CONTINUATE IN THE NEXT LINE
C
C     Note: to change the length of CARD, change in CHARACTER, READ and
C           DO.
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inpo_oms.f'
      INCLUDE 'prob_oms.f'
C
      CHARACTER SUBNAT*7, CARD*121 ! Must be 120+1
C
      LOGICAL*1 NEWLINE
C
      SAVE KPRIN    ! Save print status
C
C**** Default print status=ON. To change it use: ECHO OFF
C
      DATA KPRIN/1/
C
C**** Blank card as default
C
      DATA CARD/' '/
C
C**** FORMATS:
C
    1 FORMAT(1X,A7,' <-- ',A)
C
C**** BEGIN
C
      NNWORS=0         ! Initialize.
      NNPARS=0
      DO I=1,MAXWPS
       WORDSS(I)=' '
       PARAMS(I)=0
      END DO
      NEWLINE=.FALSE.
C
C**** Don't return without answer And continue reading if / or \
C
      DO WHILE((NNWORS.EQ.0).AND.(NNPARS.EQ.0)
     .                            .OR.NEWLINE)
C
       NEWLINE=.FALSE.            ! Initialize.
       LAST=0
       LASTP=0
C
       READ(ITAPET,'(A120)',END=2) CARD ! Read a card from input file.
c      LENG=LNBLNK(CARD)                ! Calculate the length.
       DO I=120,1,-1                    ! Calculate the length.
        IF(CARD(I:I).NE.' ') THEN
         LENG=I
         GO TO 9999
        ENDIF
       ENDDO
 9999  CONTINUE
C
       DO WHILE(LAST.LT.LENG)     ! Decode all the card.
        IFIRST=LAST+1
C
C**** Jump separators ( =:,)
C
        DO WHILE(IFIRST.LT.LENG           .AND.
     .           CARD(IFIRST:IFIRST).EQ.' '.OR.
     .           CARD(IFIRST:IFIRST).EQ.'='.OR.
     .           CARD(IFIRST:IFIRST).EQ.':'.OR.
     .           CARD(IFIRST:IFIRST).EQ.','    )
         IFIRST=IFIRST+1
        END DO
        IF(LAST.EQ.0) IFIRSP=IFIRST  ! Save IFIRST to print card
        LAST=IFIRST
C
C**** Look for end (LAST=LENG), Look for separator ( =:,),
C     Look for coment ($!/\).
C
        DO WHILE((LAST.LE.LENG)        .AND.
     .            CARD(LAST:LAST).NE.' '.AND.
     .            CARD(LAST:LAST).NE.'='.AND.
     .            CARD(LAST:LAST).NE.':'.AND.
     .            CARD(LAST:LAST).NE.','.AND.
     .            CARD(LAST:LAST).NE.'$'.AND.
     .            CARD(LAST:LAST).NE.'!'.AND.
     .            CARD(LAST:LAST).NE.'/'.AND.
     .            CARD(LAST:LAST).NE.'\'     )
         LAST=LAST+1
        END DO
        IF(CARD(LAST:LAST).EQ.'$'.OR.CARD(LAST:LAST).EQ.'!')
     .   LENG=LAST-1
        IF(CARD(LAST:LAST).EQ.'/'.OR.CARD(LAST:LAST).EQ.'\') THEN
         LENG=LAST-1           ! Deal with continuation characters.
         NEWLINE=.TRUE.        ! Set new line flag.
        END IF
        LAST=LAST-1
        IF(LAST.GE.IFIRST) THEN ! Is it a word or a parameter?
         LASTP=LAST            ! Save LAST to print card
         READ(CARD(IFIRST:LAST),'(F20.0)',IOSTAT=IFLAG) DIGIT
         IF(IFLAG.EQ.0) THEN   ! It is a parameter.
          NNPARS=NNPARS+1
          IF(NNPARS.GT.MAXWPS) GO TO 3 ! Error
          PARAMS(NNPARS)=DIGIT
         ELSE                  ! It is a word.
          NNWORS=NNWORS+1
          IF(NNWORS.GT.MAXWPS) GO TO 4 ! Error
          WORDSS(NNWORS)=CARD(IFIRST:LAST)
         END IF
        END IF ! (LAST.GE.IFIRST)
       END DO ! WHILE (LAST.LT.LENG)
C
       IF(WORDSS(1).EQ.'ECHO') THEN ! Deal with echo
        IF(WORDSS(2).EQ.'OFF') THEN
         KPRIN=0
         WRITE(LURESS,1) 'LISTEN','ECHO OFF'
        ELSE
         KPRIN=1
         WRITE(LURESS,1) 'LISTEN','ECHO ON'
        ENDIF
        NNWORS=0      ! Forget all
        NNPARS=0
        DO I=1,MAXWPS
         WORDSS(I)=' '
         PARAMS(I)=0
        END DO
       ELSE ! Print card
        IF((KPRIN.EQ.1).AND.(IFIRSP.LE.LASTP)) THEN
         IF(NEWLINE) THEN
          LASTP=LASTP+2
          CARD(LASTP-1:LASTP)=' \'
         END IF
         IF(NPRINT.EQ.0) WRITE(LURESS,1) SUBNAT,CARD(IFIRSP:LASTP)
        END IF
       ENDIF
C
      END DO ! WHILE ((NNWORS.EQ.0).AND.(NNPARS.EQ.0).OR.NEWLINE)
C
      RETURN
C
C**** ERRORS:
C
    2 IF(SUBNAT.NE.'LOADPST'.AND.ITAPET.NE.244)      ! to be implemented
     . CALL RUNENDS('LISTENT:END OF FILE DETECTED       ')
    3 CALL RUNENDS('LISTENT:TOO MANY PARAM IN COMMAND  ')
    4 CALL RUNENDS('LISTENT:TOO MANY WORDS IN COMMAND  ')
C
      END
