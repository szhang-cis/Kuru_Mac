      SUBROUTINE NUMCRK(MCRAK,NCRAK,NEWCK,NEWSY,NTYPE,PROPS,SIGMA,
     .                  SGPRI,IFLAG)
C******************************************************************
C
C****THIS ROUTINES UPDATES THE NUMBER OF ACTIVE CRACKS
C
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PROPS(*), SIGMA(*), SGPRI(*)
C
      FTULT=PROPS(32)
C
      IF(IFLAG.EQ.0) THEN
C
        NCRAK=0
C
CC$DIR SCALAR
        DO ICRAK=1,MCRAK
          IF(SGPRI(ICRAK).GT.FTULT) NCRAK=NCRAK+1
        ENDDO
C
      ELSE
C
C       Find representative stresses in material system
C
        IROTA=0
C       IF(NTYPE.EQ.4.AND.NCRAK.EQ.1) IROTA=1
C
        CALL REPVAL(IROTA,SGPRI,SIGMA)     
C
C       Check for new cracks
C
C       NEWSY=0
        NEWCK=0
        NCRK0=NCRAK
        IF(NCRAK.EQ.MCRAK) RETURN
C
CC$DIR SCALAR
        DO ICRAK=NCRAK+1,MCRAK
          IF(SGPRI(ICRAK).GT.FTULT) NCRAK=NCRAK+1
        ENDDO
C
C       Establish if the material system has rotated 
C              and if new cracks have formed
C
C       IF(IROTA.EQ.1.AND.NCRAK.GT.1) NEWSY=1
        IF(NCRK0.NE.NCRAK)            NEWCK=1
C
      ENDIF
C
      RETURN
      END
