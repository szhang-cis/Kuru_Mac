      SUBROUTINE SORTERS(NSKIPS)
C***********************************************************************
C     
C     THIS ROUTINE SORT WORDS AND PARAMETERS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** VAR
C
      INCLUDE 'inpo_oms.f'
C
C**** BEGIN
C
      NFLAGS=NNWORS-NNPARS                        ! Number of flag words
      IF(NFLAGS.LE.0) GOTO 100                    ! Error detected
      DO IDWORS=1,NSKIPS-1                        ! Number of flag words
       DO IWORDS=2,NFLAGS                         ! Skip first word
        IF(WORDSS(IWORDS).EQ.DWORDS(IDWORS))
     .   DPARAS(IDWORS)=-DPARAS(IWORDS)+1.0       ! Toggle the flag
       ENDDO
      ENDDO
      DO IDWORS=NSKIPS,MAXWPS                     ! Number of data words
       DO IWORDS=NFLAGS+1,NNWORS                  ! Skip NSKIP words
        IF(WORDSS(IWORDS).EQ.DWORDS(IDWORS))
     .   DPARAS(IDWORS)=PARAMS(IWORDS-NFLAGS)   ! Move parameter to data
       ENDDO
      ENDDO
      RETURN
  100 CALL RUNENDS('SORTERS: ILLEGAL NUMBER OF PARAMETERS')
      END
