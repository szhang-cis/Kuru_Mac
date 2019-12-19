      SUBROUTINE WMATRR(WSTIR,SHAPR,DVOLU,NGAUR,NGAUL,NNODE)
C**********************************************************************
C
C**** THIS ROUTINE EVALUATES THE "REDUCED ELEMENT MASS" MATRIX
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION WSTIR(NNODE,*), SHAPR(NNODE,*), DVOLU(*),
     .          AUXXX(8,8)
C
      DO IGAUR=1,NGAUR
       DO JGAUR=1,NGAUR
        WSTIR(IGAUR,JGAUR)=0.0
        DO IGAUS=1,NGAUL
         WSTIR(IGAUR,JGAUR)=WSTIR(IGAUR,JGAUR)+SHAPR(IGAUR,IGAUS)*
     .                      SHAPR(JGAUR,IGAUS)*DVOLU(IGAUS)
        ENDDO
       ENDDO
      ENDDO
C
      CALL INVMTR(WSTIR,AUXXX,DETER,NGAUR,NNODE,    8)
C
      DO IGAUR=1,NGAUR
       DO JGAUR=1,NGAUR
        WSTIR(IGAUR,JGAUR)=AUXXX(IGAUR,JGAUR)
       ENDDO
      ENDDO
C
      RETURN
      END
