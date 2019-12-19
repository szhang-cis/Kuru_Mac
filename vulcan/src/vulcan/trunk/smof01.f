      SUBROUTINE SMOF01(DVOLU,EMASS,SHAPE,EHIST,
     .                  SFISB,SFISA)
C***********************************************************************
C
C****THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL STRAINS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION DVOLU(*), SHAPE(NNODE,*), EHIST(NHIST,NGAUL),
     .          EMASS(NNODE,NNODE)
      DIMENSION SFISB(NDIME,NNODE), SFISA(NNODE,NDIME)
C
C***FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAUL.EQ.1) THEN
C$DIR SCALAR
        DO IDIME=1,NDIME
          DO INODE=1,NNODL
            SFISB(IDIME,INODE)=EHIST(9+IDIME,1)
          ENDDO
        ENDDO
        RETURN
      ENDIF
C
C***INTEGRAL PERFORMING
C
C$DIR SCALAR
      DO IDIME=1,NDIME
        DO INODE=1,NNODL
          DO IGAUS=1,NGAUL
            SFISA(INODE,IDIME)=SFISA(  INODE,IDIME)+
     .                         SHAPE(  INODE,IGAUS)*
     .                         EHIST(9+IDIME,IGAUS)*
     .                         DVOLU(IGAUS)
          ENDDO
        ENDDO
      ENDDO
C
C***PERFORMS MATRIX PRODUC
C
C$DIR SCALAR
      DO IDIME=1,NDIME
        DO INODE=1,NNODL
          DO JNODE=1,NNODL
            SFISB(IDIME,INODE)=SFISB(IDIME,INODE)+
     .                         SFISA(JNODE,IDIME)*EMASS(INODE,JNODE)
          ENDDO
        ENDDO
      ENDDO
      RETURN
C
      END
