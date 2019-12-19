      SUBROUTINE BBMEAR(CARTD,CMEAN,GPCOD,SHAPE,DVOLU,SHAPR,WSTIR,NGAUR,
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
      DIMENSION AUXX1(27,8),          AUXX2(27,8)
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
        AUXX1(INODE,IGAUR)=0.0
        DO IGAUX=1,NGAUL
         AUXX1(INODE,IGAUR)=AUXX1(INODE,IGAUR)+SHAPR(IGAUR,IGAUX)*
     .                      CARTD(1,INODE,IGAUX)*DVOLU(IGAUX)
        ENDDO
       ENDDO
C
       DO IGAUR=1,NGAUR
        AUXX2(INODE,IGAUR)=0.0
        DO JGAUR=1,NGAUR
         AUXX2(INODE,IGAUR)=AUXX2(INODE,IGAUR)+AUXX1(INODE,JGAUR)*
     .                      WSTIR(JGAUR,IGAUR)
        ENDDO
       ENDDO
C
       DO IGAUR=1,NGAUR
        CMEAN(1,INODE)=CMEAN(1,INODE)+AUXX2(INODE,IGAUR)*
     .                 SHAPR(IGAUR,IGAUS)
       ENDDO
C
C**** 2
C
       DO IGAUR=1,NGAUR
        AUXX1(INODE,IGAUR)=0.0
        DO IGAUX=1,NGAUL
         AUXX1(INODE,IGAUR)=AUXX1(INODE,IGAUR)+SHAPR(IGAUR,IGAUX)*
     .                      CARTD(2,INODE,IGAUX)*DVOLU(IGAUX)
        ENDDO
       ENDDO
C
       DO IGAUR=1,NGAUR
        AUXX2(INODE,IGAUR)=0.0
        DO JGAUR=1,NGAUR
         AUXX2(INODE,IGAUR)=AUXX2(INODE,IGAUR)+AUXX1(INODE,JGAUR)*
     .                      WSTIR(JGAUR,IGAUR)
        ENDDO
       ENDDO
C
       DO IGAUR=1,NGAUR
        CMEAN(2,INODE)=CMEAN(2,INODE)+AUXX2(INODE,IGAUR)*
     .                 SHAPR(IGAUR,IGAUS)
       ENDDO
C
C**** AXISYMMETRIC
C
       IF(NTYPE.EQ.3) THEN
        DO IGAUR=1,NGAUR
         AUXX1(INODE,IGAUR)=0.0
         GPCOA=0.0D0                                    ! average radius
         DO IGAUX=1,NGAUL
          GPCOA=GPCOA+GPCOD(1,IGAUX)
         ENDDO
         GPCOA=GPCOA/NGAUL
         DO IGAUX=1,NGAUL
          IF(GPCOD(1,IGAUX).GT.(1.0D-10*GPCOA)) THEN  ! see GPCOA in *.f
           AUXX1(INODE,IGAUR)=AUXX1(INODE,IGAUR)+SHAPR(IGAUR,IGAUX)*
     .                        SHAPE(INODE,IGAUX)/GPCOD(1,IGAUX)*
     .                        DVOLU(IGAUX)
          ENDIF
         ENDDO
        ENDDO
C
        DO IGAUR=1,NGAUR
         AUXX2(INODE,IGAUR)=0.0
         DO JGAUR=1,NGAUR
          AUXX2(INODE,IGAUR)=AUXX2(INODE,IGAUR)+AUXX1(INODE,JGAUR)*
     .                       WSTIR(JGAUR,IGAUR)
         ENDDO
        ENDDO
C
        DO IGAUR=1,NGAUR
         CMEAN(3,INODE)=CMEAN(3,INODE)+AUXX2(INODE,IGAUR)*
     .                  SHAPR(IGAUR,IGAUS)
        ENDDO
       ENDIF
C
C**** 3
C
       IF(NTYPE.EQ.4) THEN
        DO IGAUR=1,NGAUR
         AUXX1(INODE,IGAUR)=0.0
         DO IGAUX=1,NGAUL
          AUXX1(INODE,IGAUR)=AUXX1(INODE,IGAUR)+SHAPR(IGAUR,IGAUX)*
     .                       CARTD(3,INODE,IGAUX)*DVOLU(IGAUX)
         ENDDO
        ENDDO
C
        DO IGAUR=1,NGAUR
         AUXX2(INODE,IGAUR)=0.0
         DO JGAUR=1,NGAUR
          AUXX2(INODE,IGAUR)=AUXX2(INODE,IGAUR)+AUXX1(INODE,JGAUR)*
     .                       WSTIR(JGAUR,IGAUR)
         ENDDO
        ENDDO
C
        DO IGAUR=1,NGAUR
         CMEAN(3,INODE)=CMEAN(3,INODE)+AUXX2(INODE,IGAUR)*
     .                  SHAPR(IGAUR,IGAUS)
        ENDDO
       ENDIF
C
      ENDDO
C
      RETURN
      END
