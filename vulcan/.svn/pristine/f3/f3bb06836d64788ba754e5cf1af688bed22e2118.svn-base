      SUBROUTINE SORTER(NSKIP)
C***********************************************************************
C     
C**** THIS ROUTINE SORT WORDS AND PARAMETERS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** VAR
C
      INCLUDE 'inpo_om.f'
C
C**** BEGIN
C
      NFLAG=NNWOR-NNPAR                           ! Number of flag words
      IF(NFLAG.LE.0) GOTO 100                     ! Error detected
      DO IDWOR=1,NSKIP-1                          ! Number of flag words
       DO IWORD=2,NFLAG                           ! Skip first word
        IF(WORDS(IWORD).EQ.DWORD(IDWOR))
     .   DPARA(IDWOR)=-DPARA(IDWOR)+1.0           ! Toggle the flag
       ENDDO
      ENDDO
      DO IDWOR=NSKIP,MAXWP                        ! Number of data words
       DO IWORD=NFLAG+1,NNWOR                     ! Skip NSKIP words
        IF(WORDS(IWORD).EQ.DWORD(IDWOR))
     .   DPARA(IDWOR)=PARAM(IWORD-NFLAG)        ! Move parameter to data
       ENDDO
      ENDDO
      RETURN
  100 CALL RUNEND('SORTER:ILLEGAL NUMBER OF PARAMETERS')
      END
