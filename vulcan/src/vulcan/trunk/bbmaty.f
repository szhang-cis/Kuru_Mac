      SUBROUTINE BBMATY(BSBAR,CMEAN,CARTD,GPCOD,SHAPE)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-SHEAR STRAIN-DISPLACEMENT MATRIX
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION BSBAR(NSTR1,*), CARTD(NDIME,*),
     .          CMEAN(3,*),     GPCOD(*),
     .          SHAPE(*)
C
      NGASH=0
C
      DO INODE=1,NNODL
       LGASH=NGASH+1
       NGASH=LGASH
       BSBAR(1,LGASH)=CARTD(1,INODE)
C
       IF(NTYPE.NE.5) THEN
        MGASH=LGASH+1
        NGASH=MGASH
        BSBAR(1,MGASH)=0.0D0
        BSBAR(2,LGASH)=0.0D0
        BSBAR(2,MGASH)=CARTD(2,INODE)
        BSBAR(3,LGASH)=CMEAN(2,INODE)
        BSBAR(3,MGASH)=CMEAN(1,INODE)
        BSBAR(4,LGASH)=0.0D0
        BSBAR(4,MGASH)=0.0D0
       ENDIF
C
       IF(NTYPE.EQ.3) BSBAR(4,LGASH)=BSBAR(4,LGASH)+
     .                                             SHAPE(INODE)/GPCOD(1)
C
       IF(NTYPE.EQ.4) THEN
        NGASH=MGASH+1
        BSBAR(1,NGASH)=0.0D0
        BSBAR(2,NGASH)=0.0D0
        BSBAR(3,NGASH)=0.0D0
        BSBAR(4,LGASH)=0.0D0
        BSBAR(4,MGASH)=0.0D0
        BSBAR(4,NGASH)=CARTD(3,INODE)
        BSBAR(5,LGASH)=CMEAN(3,INODE)
        BSBAR(5,MGASH)=0.0D0
        BSBAR(5,NGASH)=CMEAN(1,INODE)
        BSBAR(6,LGASH)=0.0D0
        BSBAR(6,MGASH)=CMEAN(3,INODE)
        BSBAR(6,NGASH)=CMEAN(2,INODE)
       ENDIF
C
      ENDDO
      RETURN
C
      END
