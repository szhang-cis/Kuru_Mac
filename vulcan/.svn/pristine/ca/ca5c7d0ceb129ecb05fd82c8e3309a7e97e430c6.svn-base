      SUBROUTINE LINSCH(AALPH,GCURN,GZERO,ITRIL,LSRCH)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES AALPH FOR NEXT LINE SEARCH TRIAL, IF
C     NECESSARY
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (MTRIL=5, MTRL1=MTRIL-1, ZEROM=1.0D-15)
      PARAMETER (STOLR=0.5D0, SMINM=0.10D0, SMAXM=10.0D0)
      REAL*8 VALFA(0:MTRL1), ENERG(0:MTRL1)
      SAVE IBRKT,VALFA,ENERG
C
C**** CHECK IF LINE SEARCH IS REQUIRED
C
      IF( (ITRIL.EQ.MTRIL)                 .OR.   ! Max. number exceeded
     .    (ABS(GCURN).LE.STOLR*ABS(GZERO)) .OR.   ! Condition satisfied
     .    (ABS(GCURN).LT.ZEROM)           ) RETURN
C
      IF(ITRIL.EQ.1) THEN                 ! Initialize start values
       IBRKT=0
       VALFA(0)=0.0D0
       ENERG(0)=GZERO
      END IF
      VALFA(ITRIL)=AALPH                  ! Save current values in array
      ENERG(ITRIL)=GCURN
C
C**** DEFINE VALUES OF GA, GB, SA, SB
C
      SAVAL=AALPH                         ! Retrieve current values
      GAVAL=GCURN
      ITRL1=ITRIL-1
      SBVAL=VALFA(ITRL1)                  ! Retrieve previous values
      GBVAL=ENERG(ITRL1)
C
      IF((GAVAL*GBVAL.GE.0.0D0).AND.
     .                        (IBRKT.EQ.1)) THEN ! Retrieve older values
       ITRL2=ITRIL-2
       SBVAL=VALFA(ITRL2)
       GBVAL=ENERG(ITRL2)
       VALFA(ITRL1)=SBVAL           ! Save older value as previous value
       ENERG(ITRL1)=GBVAL
      END IF
C
      IF(GAVAL*GZERO.LT.0.0D0) IBRKT=1            ! Find bracket on zero
C
C**** ILLINOIS ALGORITHM TO FIND ZERO
C
      AALPH=SAVAL-GAVAL*(SAVAL-SBVAL)/(GAVAL-GBVAL)
C
C**** CHECK THE CALCULATED VALUE OF AALPH
C
      IF((AALPH.GT.SMAXM).OR.(AALPH.LT.SMINM))
     . AALPH=SAVAL*GZERO/(GZERO-GCURN)
      IF(AALPH.LT.SMINM) THEN
       AALPH=SMINM
      ELSE IF(AALPH.GT.SMAXM) THEN
       AALPH=SMAXM
      END IF
C
      IF(DABS(AALPH-SAVAL).LT.0.25D0*STOLR*DABS(AALPH+SAVAL)) THEN
       AALPH=VALFA(ITRIL)
       RETURN                                 ! Reject very close values
      END IF
C
      IF(IBRKT.EQ.0.AND.
     .               ITRIL.EQ.(MTRIL-1)) THEN ! Choose AALPH for min. G.
       GMINM=ENERG(1)
       NOTRL=1
       DO JTRIL = 2,ITRIL
        IF(DABS(ENERG(JTRIL)).LT.DABS(GMINM)) GMINM=ENERG(JTRIL)
        IF(DABS(ENERG(JTRIL)).EQ.DABS(GMINM)) NOTRL=JTRIL
       END DO
       AALPH=VALFA(NOTRL)
       IF(NOTRL.EQ.ITRIL) RETURN
      END IF
C
      LSRCH=1                                      ! Set mark to iterate
C
      END
