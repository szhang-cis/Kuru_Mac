      SUBROUTINE PCGPRE(ALOAD,GSTDI,WORKS,KWITD,NDOFN,NPOIN)
C***********************************************************************
C
C****THIS ROUTINE APPLIES PRECONDITIONING
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ALOAD(*), GSTDI(*), WORKS(*)
C
      IF(KWITD.EQ.0) THEN
        DO ITOTV=1,NPOIN*NDOFN
          ALOAD(ITOTV)=GSTDI(ITOTV)*ALOAD(ITOTV)
        ENDDO
      ELSE
        DO IPOIN=1,NPOIN
          NSTA1=(IPOIN-1)*(NDOFN*NDOFN)+1
          NSTA2=(IPOIN-1)*NDOFN+1
          CALL PCGMUL(GSTDI(NSTA1),ALOAD(NSTA2),WORKS,NDOFN)
        ENDDO
      ENDIF
C
      RETURN
      END   
