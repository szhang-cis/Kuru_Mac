      SUBROUTINE PLBMTX(DMATX,DMTEP,NSTRS,NSTR1,KSYMM)
C***********************************************************************
C
C**** THIS SUBROUTINE TRANSFORM THE ELASTO-PLASTIC VECTOR TO MATERIAL
C     TENSOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'inte_om.f'
C
      DIMENSION DMATX(NSTRS,*), DMTEP(*)
C
      IF(KUNLD.EQ.1) RETURN
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMATX(ISTRS,JSTRS)=DMTEP(IKONT)
        IF(KSYMM.EQ.1) DMATX(JSTRS,ISTRS)=DMATX(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
      RETURN
      END
