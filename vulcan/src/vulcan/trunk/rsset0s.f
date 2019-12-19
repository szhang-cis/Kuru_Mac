      SUBROUTINE RSSET0S
C***********************************************************************
C
C**** THIS ROUTINE:
C
C     - SETS UP ELEMENTAL POINTERS:             IDATA, IPREV, ISTAT,
C                                               IMATX
C     - SETS UP PROCESSING AREA POINTERS:       IDATP, NLENP, NRECP
C     - SETS UP CONVERGED STEPS AREA POINTERS:  IDATC, NLENC, NRECC,
C                                               NRECG 
C
C***********************************************************************
C
      CALL ADDELMS
C
      CALL RSSETPS
C
      CALL RSSETCS
C
      RETURN
      END
