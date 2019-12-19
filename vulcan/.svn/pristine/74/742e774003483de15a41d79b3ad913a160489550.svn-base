      SUBROUTINE SMOF30(DVOLU,EMASS,SHAPE,STRAN,
     .                  SFISB,SFISA)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL STRAINS
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      COMMON/SMOGAU/FACTA
C
      DIMENSION DVOLU(*),           SHAPE(NNODS,*),
     .          STRAN(NSTR1,*),     EMASS(NNODS,*)
      DIMENSION SFISB(NSTR1,*),     SFISA(NNODS,*)
C
C**** COMPUTE TOTAL VOLUME IF NECESSARY
C
      IF(IABS(KSGAU).EQ.2) THEN       ! local smoothing
       FACTA=0.D+00
       DO IGAUS=1,NGAUL
        FACTA=FACTA+DVOLU(IGAUS)
       ENDDO
      ENDIF
C
C**** FOR 'ONE-POINT' INTEGRATION RULES
C
      IF(NGAUL.EQ.1) THEN
C$DIR SCALAR
       DO ISTRE=1,NSTR1
        DO INODE=1,NNODS
         SFISB(ISTRE,INODE)=STRAN(ISTRE,1)
        ENDDO
       ENDDO
       RETURN
      ENDIF
C
C**** INTEGRAL PERFORMING
C
C$DIR SCALAR
      DO ISTRE=1,NSTR1
       DO INODE=1,NNODS
        DO IGAUS=1,NGAUL
         SFISA(INODE,ISTRE)=SFISA(INODE,ISTRE)+
     .                      SHAPE(INODE,IGAUS)*
     .                      STRAN(ISTRE,IGAUS)*
     .                      DVOLU(IGAUS)
        ENDDO
       ENDDO
      ENDDO
C
C**** PERFORMS MATRIX PRODUC
C
C$DIR SCALAR
      DO ISTRE=1,NSTR1
       DO INODE=1,NNODS
        DO JNODE=1,NNODS
         SFISB(ISTRE,INODE)=SFISB(ISTRE,INODE)+
     .                      SFISA(JNODE,ISTRE)*EMASS(INODE,JNODE)
        ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
