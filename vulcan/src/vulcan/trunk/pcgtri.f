      SUBROUTINE PCGTRI(GSTDI,IFFIX,NDOFN,NPOIN,KWIDT,TOLCX)
C***********************************************************************
C
C****THIS ROUTINE INVERTS THE GLOBAL PRECONDITIONING MATRIX
C
C....INPUT PARAMETERS
C
C       GSTDI(*)            - PRECONDITIONING MATRIX
C       KWITD               - FLAG FOR PRECONDITIONING
C                                =0 DIAGONAL
C                                <0 BLOCK-DIAGONAL
C
C.... OUTPUT PARAMETERS
C
C       GSTDI(*)           - INVERSE OF PRECONDITIONING MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION GSTDI(*), IFFIX(*)
C
      ZEROS=0.0D0
      UNITY=1.0D0
C
      IF(KWIDT.EQ.0) THEN                    ! diagonal preconditioning
        NTOTV=NPOIN*NDOFN
        DO ITOTV=1,NTOTV
          IF(DABS(GSTDI(ITOTV)).GT.TOLCX) THEN
            GSTDI(ITOTV)=UNITY/GSTDI(ITOTV)
          ELSE
            GSTDI(ITOTV)=ZEROS
          ENDIF
        ENDDO
      ELSE                                   ! block-diagonal precond.
        DO IPOIN=1,NPOIN
          NSTAR=(IPOIN-1)*(NDOFN*NDOFN)
          ITOTV=(IPOIN-1)*NDOFN
          DO IDOFN=1,NDOFN
            ITOTV=ITOTV+1
            IF(IFFIX(ITOTV).NE.0) THEN       ! delete row and column
              DO JDOFN=1,NDOFN
                GSTDI(NSTAR+(JDOFN-1)*NDOFN+IDOFN)=ZEROS ! row
                GSTDI(NSTAR+(IDOFN-1)*NDOFN+JDOFN)=ZEROS ! column
              ENDDO
              GSTDI(NSTAR+(IDOFN-1)*NDOFN+IDOFN)=UNITY ! one in diagonal
            ENDIF
          ENDDO
          CALL INVERT(GSTDI(NSTAR+1),NDOFN,NDOFN)
        ENDDO
      ENDIF
C
      RETURN
      END
