      SUBROUTINE LINEAB(BSBAR,DEPSV,DESIG,DMATX,DSTRA,ELDIS,GPCOD,
     .                  LARGE,NDIME,NDOFN,NNODE,NSTR1,PROPS,SHAPE,
     .                  STRAN,STRA0,TSTRA,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .                  ELCOD,CARTD,NEVAB,
     .                  XJACN,XJANI,XJA3N,XJ3NI,DETJN,DETJB,
     .                  ITASL)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES STRAINS AND INCREMENTAL ELASTIC STRESSES
C
C     Notes:
C
C     THE JACOBIAN MATRIX IS THE DEFORMATION GRADIENT FOR LARGE > 0
C
C     FOR LARGE=1 (XJACM = MATERIAL DISPLACEMENTS GRADIENT + IDENTITY)
C
C     FOR LARGE=2 (XJACM = UPDATED INCREMENTAL DISPLACEMENTS GRADIENT +
C                          IDENTITY)
C
C     IN BOTH CASES:
C     ELDIS = SPATIAL COORDINATES = MATERIAL COORDINATES + DISPLACEMENTS
C
C     THE SAME CONCEPTS ARE APPLIED TO THE MEAN JACOBIAN MATRIX (XJACN)
C
C     FOR AXISYMMETRIC CASES, THE IMPLEMENTED FORM TO COMPUTE TSTRA(4)
C     (SEE BELOW) IS ABSOLUTELY EQUIVALENT TO THE CLASSICAL ONE GIVEN
C     BY (see linear.f):
C 
C     TSTRA(4)=0.0D0
C     DO INODE=1,NNODE
C      ELDCO=ELDIS(1,INODE)-ELCOD(1,INODE)       ! total displacements
C      TSTRA(4)=TSTRA(4)+ELDCO*SHAPE(INODE)/GPCOD(1)
C     ENDDO
C     TSTRA(4)=TSTRA(4)+0.5D0*TSTRA(4)*TSTRA(4)
C
C     LARGE flag independent computation of TSTRA will be attempted for 
C     large strains
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION BSBAR(NSTR1,*), DMATX(NSTRS,*), DESIG(*), DSTRA(*),
     .          ELDIS(NDIME,*), PROPS(*),       GPCOD(*), SHAPE(*),
     .          STRAN(*),       STRA0(*),       TSTRA(*)
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION ELCOD(NDIME,*)
      DIMENSION CARTD(NDIME,*), XJACN(NDIME,*), XJANI(NDIME,*)
      DIMENSION XJACA(3,3)
C
C**** CALCULATE TOTAL STRAINS
C
      IF(LARGE.EQ.0) THEN
       DO ISTR1=1,NSTRE
        TSTRA(ISTR1)=0.0D0
        DO INODE=1,NNODE
         DO IDIME=1,NDIME
          IEVAB=(INODE-1)*NDIME+IDIME
          TSTRA(ISTR1)=TSTRA(ISTR1)+BSBAR(ISTR1,IEVAB)*
     .                 ELDIS(IDIME,INODE)
         ENDDO
        ENDDO
       ENDDO
       IF(NTYPE.EQ.1.OR.NTYPE.EQ.2)
     .  TSTRA(4)=0.0D0 ! for ntype=1, this comp. is computed in planes.f
C
      ELSE             ! large strains
C
       IF(LARGE.EQ.1) THEN                         ! large strains (TLF)
C
C**** CALCULATE THE JACOBIAN MATRIX
C
        CALL PROMA2(XJACM,ELDIS,CARTD,NDIME,NNODE)
C
C**** CALCULATE THE DETERMINANT AND INVERS OF JACOBIAN MATRIX
C
        CALL INVMTX(XJACM,XJACI,DETJM,NDIME)
C
        XJA3M=1.0D0
        XJA3I=1.0D0
        IF(NTYPE.EQ.3) THEN
         GPCOI=0.0D0
         GPCOA=0.0D0                 ! average radial spatial coordinate
         DO INODE=1,NNODE
          GPCOI=GPCOI+ELDIS(1,INODE)*SHAPE(INODE)
          GPCOA=GPCOA+ELDIS(1,INODE)
         ENDDO
         GPCOA=GPCOA/NNODE
         IF(GPCOD(1).GT.(1.0D-10*GPCOA).AND.
     .         GPCOI.GT.(1.0D-10*GPCOA)) THEN         ! see GPCOA in *.f
          DETJM=DETJM*GPCOI/GPCOD(1)
          XJA3M=GPCOI/GPCOD(1)
          XJA3I=1.0D0/XJA3M
         ENDIF
        ENDIF
        DETJB=DETJM
C
C**** CALCULATE THE MEAN JACOBIAN MATRIX
C
        CALL PROMA2X(XJACN,ELDIS,BSBAR,NDIME,NNODE,NSTR1,NEVAB)
C
C**** CALCULATE THE DETERMINANT AND INVERS OF THE MEAN JACOBIAN MATRIX
C
        CALL INVMTX(XJACN,XJANI,DETJN,NDIME)
C
        XJA3N=1.0D0
        XJ3NI=1.0D0
        IF(NTYPE.EQ.3) THEN
         GPCOI=0.0D0
         DO INODE=1,NNODE
          GPCOI=GPCOI+ELDIS(1,INODE)*BSBAR(3,INODE)
         ENDDO
         DETJN=DETJN*GPCOI
         XJA3N=GPCOI
         XJ3NI=1.0D0/XJA3N
        ENDIF
C
        KPARB=INT(PPARB)
        IF(KPARB.LT.0) KPARB=-KPARB
C
        GO TO (21,22,23,21,22,23,21,22,23,
     .         24,24,24,25,25,25,26,26,26,
     .         27,27,27,28,28,28), KPARB
C
C**** B-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   21   DEFRA=(DETJN/DETJM)**(1.0D0/3.0D0)
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)*DEFRA
          XJACI(IDIME,JDIME)=XJACI(IDIME,JDIME)/DEFRA
c         XJACM(IDIME,JDIME)=XJACN(IDIME,JDIME)               ! cylinder
c         XJACI(IDIME,JDIME)=XJANI(IDIME,JDIME)
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
        IF(NTYPE.EQ.3) THEN
         XJA3M=XJA3M*DEFRA
         XJA3I=XJA3I/DEFRA
c        XJA3M=XJA3N                                          ! cylinder
c        XJA3I=1.0D0/XJA3M
        ENDIF
        DETJB=DETJN
        DIFVO=0.0D0
        DSHEA=0.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** B-SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   22   DEFRA=1.0D0
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACA(IDIME,JDIME)=XJACN(IDIME,JDIME)
         ENDDO
        ENDDO
        DIFVO=0.0D0
        DSHEA=0.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** B-BAR-SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   23   DEFRA=(DETJN/DETJM)**(1.0D0/3.0D0)
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)*DEFRA
          XJACI(IDIME,JDIME)=XJACI(IDIME,JDIME)/DEFRA
          XJACA(IDIME,JDIME)=XJACN(IDIME,JDIME)
         ENDDO
        ENDDO
        IF(NTYPE.EQ.3) THEN
         XJA3M=XJA3M*DEFRA
         XJA3I=XJA3I/DEFRA
        ENDIF
        DETJB=DETJN
        DIFVO=0.0D0
        DSHEA=0.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** OTHER B-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   24   DEFRA=1.0D0
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
C
        TSTM1=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-1.0D0)
        TSTM2=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-1.0D0)
        TSTM4=0.0D0
        TSTN1=0.5D0*(XJACN(1,1)*XJACN(1,1)+XJACN(2,1)*XJACN(2,1)-1.0D0)
        TSTN2=0.5D0*(XJACN(1,2)*XJACN(1,2)+XJACN(2,2)*XJACN(2,2)-1.0D0)
        TSTN4=0.0D0
        IF(NTYPE.EQ.3) THEN
         TSTM4=0.5D0*(XJA3M*XJA3M-1.0D0)
         TSTN4=0.5D0*(XJA3N*XJA3N-1.0D0)
        ENDIF
        IF(NTYPE.EQ.4) THEN
         TSTM1=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)
     .             +XJACM(3,1)*XJACM(3,1)-1.0D0)
         TSTM2=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)
     .             +XJACM(3,2)*XJACM(3,2)-1.0D0)
         TSTM4=0.5D0*(XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)
     .             +XJACM(3,3)*XJACM(3,3)-1.0D0)
         TSTN1=0.5D0*(XJACN(1,1)*XJACN(1,1)+XJACN(2,1)*XJACN(2,1)
     .             +XJACN(3,1)*XJACN(3,1)-1.0D0)
         TSTN2=0.5D0*(XJACN(1,2)*XJACN(1,2)+XJACN(2,2)*XJACN(2,2)
     .             +XJACN(3,2)*XJACN(3,2)-1.0D0)
         TSTN4=0.5D0*(XJACN(1,3)*XJACN(1,3)+XJACN(2,3)*XJACN(2,3)
     .             +XJACN(3,3)*XJACN(3,3)-1.0D0)
        ENDIF
        DIFVO=1.0D0/3.0D0*(TSTN1+TSTN2+TSTN4-TSTM1-TSTM2-TSTM4)
        DSHEA=0.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** B-e_SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   25   DEFRA=1.0D0
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
        DIFVO=0.0D0
        DSHEA=1.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** B-(e_SHEAR-BAR) (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEG.)
C
   26   DEFRA=(DETJN/DETJM)**(1.0D0/3.0D0)
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)*DEFRA
          XJACI(IDIME,JDIME)=XJACI(IDIME,JDIME)/DEFRA
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
        IF(NTYPE.EQ.3) THEN
         XJA3M=XJA3M*DEFRA
         XJA3I=XJA3I/DEFRA
        ENDIF
        DETJB=DETJN
        DIFVO=0.0D0
        DSHEA=1.0D0
        DCO33=0.0D0
        GO TO 34
C
C**** B-BAR_F33 (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   27   DEFRA=(DETJN/DETJM)**(1.0D0/3.0D0)
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)*DEFRA
          XJACI(IDIME,JDIME)=XJACI(IDIME,JDIME)/DEFRA
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
        IF(NTYPE.EQ.3) THEN
         XJA3M=XJA3M*DEFRA
         XJA3I=XJA3I/DEFRA
        ENDIF
        DETJB=DETJN
        DIFVO=0.0D0
        DSHEA=0.0D0
        DCO33=1.0D0
        GO TO 34
C
C**** J-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   28   DEFRA=1.0D0
        DO IDIME=1,NDIME
         DO JDIME=1,NDIME
          XJACA(IDIME,JDIME)=XJACM(IDIME,JDIME)
         ENDDO
        ENDDO
        DETJB=DETJN
        DIFVO=0.0D0
        DSHEA=0.0D0
        DCO33=0.0D0
        GO TO 34
C
   34   CONTINUE
C
        IF(ITASL.EQ.1) RETURN
C
        GO TO (11,12,13,14,15), NTYPE
C
   11   TSTRA(1)=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,1)*XJACM(1,1)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,1)*XJACM(2,1)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,1)*XJACM(2,1)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(2)=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,2)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,2)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,2)*XJACM(2,2)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(3)=     XJACA(1,1)*XJACA(1,2)+XJACA(2,1)*XJACA(2,2)-
     .           DSHEA*(XJACM(1,1)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))+
     .                  XJACM(2,1)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))+
     .                 (XJACM(1,1)*XJACM(2,2)+XJACM(1,2)*XJACM(2,1))*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(4)=0.0D0          ! this component is computed in planes.f
        GO TO 7
   12   TSTRA(1)=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,1)*XJACM(1,1)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,1)*XJACM(2,1)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,1)*XJACM(2,1)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(2)=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,2)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,2)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,2)*XJACM(2,2)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(3)=     XJACA(1,1)*XJACA(1,2)+XJACA(2,1)*XJACA(2,2)-
     .           DSHEA*(XJACM(1,1)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))+
     .                  XJACM(2,1)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))+
     .                 (XJACM(1,1)*XJACM(2,2)+XJACM(1,2)*XJACM(2,1))*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(4)=0.0D0
        GO TO 7
   13   TSTRA(1)=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,1)*XJACM(1,1)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,1)*XJACM(2,1)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,1)*XJACM(2,1)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(2)=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)-
     .           1.0D0)+
     .           DIFVO+
     .           DSHEA*(0.5D0-
     .                  0.5D0*XJACM(1,2)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))-
     .                  0.5D0*XJACM(2,2)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))-
     .                  XJACM(1,2)*XJACM(2,2)*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(3)=     XJACA(1,1)*XJACA(1,2)+XJACA(2,1)*XJACA(2,2)-
     .           DSHEA*(XJACM(1,1)*XJACM(1,2)*
     .                    (XJACI(1,1)*XJACI(1,1)+XJACI(2,1)*XJACI(2,1))+
     .                  XJACM(2,1)*XJACM(2,2)*
     .                    (XJACI(1,2)*XJACI(1,2)+XJACI(2,2)*XJACI(2,2))+
     .                 (XJACM(1,1)*XJACM(2,2)+XJACM(1,2)*XJACM(2,1))*
     .                    (XJANI(1,1)*XJANI(1,2)+XJANI(2,1)*XJANI(2,2)))
        TSTRA(4)=(1.0D0-DCO33)*(0.5D0*(XJA3M*XJA3M-1.0D0))+DIFVO+
     .                   DCO33*(0.5D0*(XJA3N*XJA3N-1.0D0))
        GO TO 7
   14   TSTRA(1)=0.5D0*(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)
     .                 +XJACM(3,1)*XJACM(3,1)-1.0D0)+DIFVO
        TSTRA(2)=0.5D0*(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)
     .                 +XJACM(3,2)*XJACM(3,2)-1.0D0)+DIFVO
        TSTRA(3)=     XJACA(1,1)*XJACA(1,2)+XJACA(2,1)*XJACA(2,2)
     .               +XJACA(3,1)*XJACA(3,2)
        TSTRA(4)=0.5D0*(XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)
     .                 +XJACM(3,3)*XJACM(3,3)-1.0D0)+DIFVO
        TSTRA(5)=     XJACA(1,1)*XJACA(1,3)+XJACA(2,1)*XJACA(2,3)
     .               +XJACA(3,1)*XJACA(3,3)
        TSTRA(6)=     XJACA(1,2)*XJACA(1,3)+XJACA(2,2)*XJACA(2,3)
     .               +XJACA(3,2)*XJACA(3,3)
        GO TO 7
   15   TSTRA(1)=0.5D0*(XJACM(1,1)*XJACM(1,1)-1.0D0)
    7   CONTINUE
       ENDIF                  ! large.eq.1
C
       IF(LARGE.EQ.2) THEN                         ! large strains (ULF)
        CALL RUNEND('ERROR IN LINEAB - LARGE=2          ')
       ENDIF                  ! large.eq.2
C
      ENDIF                   ! large.eq.0
C
C**** AND THE CORRESPONDING INCREMENTAL STRAINS, NEW TOTAL STRAINS
C
      IF(NMEMOM.EQ.0) THEN
       DO ISTRE=1,NSTR1
        DSTRA(ISTRE)=TSTRA(ISTRE)-STRAN(ISTRE)    ! INCR. STRAINS
        STRAN(ISTRE)=TSTRA(ISTRE)                 ! UPDATE TOTAL STRAINS
       ENDDO
      ELSE
       DO ISTRE=1,NSTR1
        DSTRA(ISTRE)=TSTRA(ISTRE)-STRAN(ISTRE)
     .                           -STRA0(ISTRE)    ! INCR. STRAINS
        STRAN(ISTRE)=TSTRA(ISTRE)-STRA0(ISTRE)    ! UPDATE TOTAL STRAINS
        TSTRA(ISTRE)=STRAN(ISTRE)               
       ENDDO
      ENDIF
C
      DEPSV=DSTRA(1)+DSTRA(2)+DSTRA(4)
C
C**** AND THE INCREMENTAL ELASTIC STRESSES 
C
      DO ISTRE=1,NSTR1
       DESIG(ISTRE)=0.0D0
      ENDDO
      IF(NMEMOM.EQ.1) THEN
       DO ISTRE=1,NSTRS
        DO JSTRE=1,NSTRS
         DESIG(ISTRE)=DESIG(ISTRE)+DMATX(ISTRE,JSTRE)*DSTRA(JSTRE)
        ENDDO
       ENDDO
      ENDIF
C
      RETURN
      END
