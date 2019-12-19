      SUBROUTINE BBMEAS(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
     .                  IGAUS)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE VOLUMETRIC-BAR STRAIN DISPLACEMENT
C     MATRIX CMEAN(3,NNODL)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION CARTD(NDIME,NNODL,*), CMEAN(3,*),
     .          GPCOD(NDIME,*),       SHAPE(NNODL,*),
     .          DVOLU(*),             SHAPR(NNODL,*),
     .          WSTIR(NNODL,*)
C
      DO INODE=1,NNODL
C
       CMEAN(1,INODE)=0.0D0
       CMEAN(2,INODE)=0.0D0
       CMEAN(3,INODE)=0.0D0
C
C**** 1
C
       DO IGAUR=1,NGAUR
        CMEAN(1,INODE)=CMEAN(1,INODE)+CARTD(1,INODE,IGAUR)*
     .                 SHAPR(IGAUR,IGAUS)
       ENDDO
C
C**** 2
C
       DO IGAUR=1,NGAUR
        CMEAN(2,INODE)=CMEAN(2,INODE)+CARTD(2,INODE,IGAUR)*
     .                 SHAPR(IGAUR,IGAUS)
       ENDDO
C
C**** AXISYMMETRIC
C
       IF(NTYPE.EQ.3) THEN
        DO IGAUR=1,NGAUR
         GPCOA=0.0D0                                    ! average radius
         DO IGAUX=1,NGAUR
          GPCOA=GPCOA+GPCOD(1,IGAUX)
         ENDDO
         GPCOA=GPCOA/NGAUR
         IF(GPCOD(1,IGAUR).GT.(1.0D-10*GPCOA)) THEN   ! see GPCOA in *.f
          CMEAN(3,INODE)=CMEAN(3,INODE)+
     .                   SHAPE(INODE,IGAUR)/GPCOD(1,IGAUR)*
     .                   SHAPR(IGAUR,IGAUS)
         ENDIF
        ENDDO
       ENDIF
C
C**** 3
C
       IF(NTYPE.EQ.4) THEN
        DO IGAUR=1,NGAUR
         CMEAN(3,INODE)=CMEAN(3,INODE)+CARTD(3,INODE,IGAUR)*
     .                  SHAPR(IGAUR,IGAUS)
        ENDDO
       ENDIF
C
      ENDDO
      RETURN
C
      END
