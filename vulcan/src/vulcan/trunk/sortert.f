      SUBROUTINE SORTERT(NSKIPT)
C***********************************************************************
C     
C     THIS ROUTINE SORT WORDS AND PARAMETERS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** VAR
C
      INCLUDE 'inpo_omt.f'
C
C**** BEGIN
C
      NFLAGT=NNWORT-NNPART                        ! Number of flag words
      IF(NFLAGT.LE.0) GOTO 100                    ! Error detected
      DO IDWORT=1,NSKIPT-1                        ! Number of flag words
       DO IWORDT=2,NFLAGT                         ! Skip first word
        IF(WORDST(IWORDT).EQ.DWORDT(IDWORT))
     .   DPARAT(IDWORT)=-DPARAT(IWORDT)+1.0       ! Toggle the flag
       ENDDO
      ENDDO
      DO IDWORT=NSKIPT,MAXWPT                     ! Number of data words
       DO IWORDT=NFLAGT+1,NNWORT                  ! Skip NSKIP words
        IF(WORDST(IWORDT).EQ.DWORDT(IDWORT))
     .   DPARAT(IDWORT)=PARAMT(IWORDT-NFLAGT)   ! Move parameter to data
       ENDDO
      ENDDO
      RETURN
  100 CALL RUNENDT('SORTERT: ILLEGAL NUMBER OF PARAMETERS')
      END
