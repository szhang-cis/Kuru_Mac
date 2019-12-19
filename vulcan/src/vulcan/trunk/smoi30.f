      SUBROUTINE SMOI30(DVOLU,EMASS,SHAPE,EHIST,PROPS,
     .                  SFIPB,SFIPA)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ELEMENTAL CONTRIBUTION TO NODAL INTERNAL
C     VARIABLES
C
C     Note: the smoothing operations may depend on NCRIT (it is not
C           considered now)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
      INCLUDE 'inpo_om.f'
C
      COMMON/SMOGAU/FACTA
C
      DIMENSION DVOLU(*),           SHAPE(NNODS,*),
     .          EHIST(NHIST,*),     EMASS(NNODS,*),
     .          PROPS(*)
      DIMENSION SFIPB(NNUIN,*),     SFIPA(NNODS,*)
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
       DO INUIN=1,NNUIN
        IPUNS=IPLAO(INUIN)
C$DIR SCALAR
        DO INODE=1,NNODS
         SFIPB(INUIN,INODE)=EHIST(IPUNS,1)
        ENDDO
C
       ENDDO       ! inuin=1,nnuin
       RETURN
      ENDIF
C
C**** INTEGRAL PERFORMING       ! to be improved !!!!!!!
C
      DO INUIN=1,NNUIN
       IPUNS=IPLAO(INUIN)
C$DIR SCALAR
       DO INODE=1,NNODS
        DO IGAUS=1,NGAUL
         SFIPA(INODE,INUIN)=SFIPA(INODE,INUIN)+SHAPE(INODE,IGAUS)*
     .                      EHIST(IPUNS,IGAUS)*DVOLU(IGAUS)
        ENDDO
       ENDDO
C
      ENDDO        ! inuin=1,nnuin
C
C**** PERFORMS MATRIX PRODUCT
C
      DO INUIN=1,NNUIN
C$DIR SCALAR
       DO INODE=1,NNODS
        DO JNODE=1,NNODS
         SFIPB(INUIN,INODE)=SFIPB(INUIN,INODE)+
     .                      SFIPA(JNODE,INUIN)*EMASS(INODE,JNODE)
        ENDDO
       ENDDO
C
      ENDDO        ! inuin=1,nnuin
C
      RETURN
      END
