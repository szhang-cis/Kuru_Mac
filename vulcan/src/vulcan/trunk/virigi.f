      SUBROUTINE VIRIGI(HTRIX,DTIME,NSTRS,DMATX,DMTEP,NSTR1,DMATE,KSYMM)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATE THE ELASTO-PLASTIC MATERIAL TENSOR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DMATX(NSTRS,*), DMATE(6,6), DMTEP(*)
      DIMENSION HTRIX(6,6),     TMATX(6,6)
C
C**** INVERTS ELASTIC CONSTITUTIVE TENSOR
C
      CALL INVERT(DMATE,NSTRS,    6)
C
C**** ADDS THIS TENSOR TO THE CONTRIBUTION GIVEN BY "HTRIX"
C
      DO ISTRS=1,NSTRS
       DO JSTRS=1,NSTRS
        TMATX(ISTRS,JSTRS)=DMATE(ISTRS,JSTRS)+DTIME*HTRIX(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
C**** INVERTS THIS TENSOR => PSEUDO VISCOPLASTIC CONSTITUTIVE TENSOR
C
      CALL INVERT(TMATX,NSTRS,    6)
C
C**** LOAD IT IN "DMTEP"
C
      IKONT=0
      DO ISTRS=1,NSTRS
       INDEX=ISTRS
       IF(KSYMM.EQ.0) INDEX=1
       DO JSTRS=INDEX,NSTRS
        IKONT=IKONT+1
        DMTEP(IKONT)=TMATX(ISTRS,JSTRS)
       ENDDO
      ENDDO
C
      RETURN
      END