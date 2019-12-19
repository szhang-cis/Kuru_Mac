      SUBROUTINE HOURMX(CARTD,EPMTX,DVOLU,ELCOD,GPCOD,PROPS,STIFH,
     .                  THICK)
C*****************************************************************
C
C**** THIS ROUTINE COMPUTES THE HOURGLASS STABILIZATION MATRIX
C
C*****************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION CARTD(NDIME,NNODE,*), ELCOD(NDIME,*), EPMTX(*),
     .          GPCOD(NDIME,*),       PROPS(*),       STIFH(*)
      DIMENSION DMATX(6,6)
C
      IF((NDIME.EQ.2).AND.((NNODL.NE.4).OR.(NGAUL.NE.1))) RETURN
      IF((NDIME.EQ.3).AND.((NNODL.NE.8).OR.(NGAUL.NE.1))) RETURN
C
      DO IKOVA=1,NKOVA
       STIFH(IKOVA)=0.0
      ENDDO
C
      IF(NDIME.EQ.2) THEN
        CALL HOURK4(CARTD,EPMTX,DVOLU,ELCOD,GPCOD,HPARA,PROPS,
     .              STIFH,DMATX,NDIME,NEVAB,NNODE,NSTRS,NTYPE,
     .              THICK)
      ELSE ! (NDIME.EQ.3)
        CALL HOURK8(CARTD,EPMTX,DVOLU,ELCOD,GPCOD,HPARA,PROPS,
     .              STIFH,DMATX,NDIME,NEVAB,NNODE,NSTRS,NTYPE)
      ENDIF
C
      RETURN
      END
