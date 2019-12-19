      SUBROUTINE SMOMTX(ELCOD,DVOLU,SHAPE,EMASS)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE ELEMENT SMOOTHING MATRIX AND INVERTS IT
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION ELCOD(NDIME,*), EMASS(NNODS,*),
     .          SHAPE(NNODS,*), DVOLU(*)
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAUL.EQ.1) RETURN
C
C**** "MASS" MATRIX
C
      DO INODE=1,NNODS
       DO JNODE=INODE,NNODS
        EMASS(JNODE,INODE)=0.0
        DO IGAUL=1,NGAUL
         EMASS(JNODE,INODE)=EMASS(JNODE,INODE)+SHAPE(INODE,IGAUL)*
     .                      SHAPE(JNODE,IGAUL)*DVOLU(IGAUL)
         ENDDO
       ENDDO
      ENDDO
C
      DO INODE=2,NNODS
       DO JNODE=1,INODE-1
        EMASS(JNODE,INODE)=EMASS(INODE,JNODE)
       ENDDO
      ENDDO
C
      IF(ISMO2.EQ.1) CALL LUMPMMT(EMASS,NNODS,     0)
C
      CALL INVERT(EMASS,NNODS,NNODS)
C
      RETURN
      END
