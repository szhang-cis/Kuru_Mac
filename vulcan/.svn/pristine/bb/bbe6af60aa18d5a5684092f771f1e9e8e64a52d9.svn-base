      SUBROUTINE BBMATYLS(CARTD,ELDIS,GPCOD,
     .                   NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .                   NTYPE,SHAPE,NSTR1,NSTRS,BMATX,CMEAN,
     .                   XJACM,XJA3M,DETJM,
     .                   XJACN,XJA3N,DETJN,
     .                   XJACI,XJA3I,
     .                   XJANI,XJ3NI)
C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-e_SHEAR STRAIN-DISPLACEMENT MATRIX
C     FOR LARGE STRAINS (LARGE=1)
C
C     Notes:
C
C     This routine is based on blarge.f
C     CMEAN in this routine has the dimensions of BSBAR
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CARTD(NDIME,*), 
     .          ELDIS(NDOFC,*), SHAPE(*),
     .          XJACM(NDIME,*), XJACN(NDIME,*),
     .          XJACI(NDIME,*), XJANI(NDIME,*)
      DIMENSION BMATX(NSTRE,*), CMEAN(NSTR1,*)


      dimension rcgtt(6)         ! borrar

C***********************************************************************
C
C**** THIS SUBROUTINE EVALUATES THE B-(e_SHEAR-BAR) STRAIN-DISPLACEMENT
C     MATRIX FOR LARGE STRAINS (LARGE=1)
C
C     Notes:
C
C     This routine is based on blarge.f
C     CMEAN in this routine has the dimensions of BSBAR
C
C     E=Green-Lagrange deformation tensor
C     E=E_1+E_2, where E_1 & E_2 are the parts of E affected and not
C     affected by DEFRA, respectively.
C
C     dE/dP=dE_1/dP+dE_2/dP              P=U,V,W
C
C     dE_1i/dP=AEiP1*dN/dX+AEiP2*(dN/dX)-bar+
C              AEiP3*dN/dY+AEiP4*(dN/dY)-bar+
C              AEiP5*dN/dZ+AEiP6*(dN/dZ)-bar+
C              AEiP7*N/X+AEiP8*(N/X)-bar
C                                        i=1,2,3,4,5,6 (components of E)
C
C     dE_2i/dP=BEiP1*dN/dX+BEiP2*(dN/dX)-bar+
C              BEiP3*dN/dY+BEiP4*(dN/dY)-bar+
C              BEiP5*dN/dZ+BEiP6*(dN/dZ)-bar
C              BEiP7*N/X+BEiP8*(N/X)-bar
C                                        i=1,2,3,4,5,6 (components of E)
C
C     BEiP2, BEiP4, BEiP6 & BEiP8 are always zero
C
C***********************************************************************
C
C**** PRELIMINARY COMPUTATIONS
C
      DEFRA=1.0D0
      DEFR2=1.0D0
      B23=2.0D0/3.0D0
C
      F11=  XJACM(1,1)/DEFRA
      F11I= XJACI(1,1)*DEFRA
      F11B= XJACN(1,1)
      F11BI=XJANI(1,1)
C
      E11B1=0.5D0*F11*F11                              ! part 1 of E-bar
      AE1U1=F11
      AE1U2=0.0D0
      BE1U1=-F11*F11I*F11I+
     .       F11*F11*F11I*F11I*F11I
      BE1U2=0.0D0
      DXNU1=1.0D0
      DXMU1=1.0D0
      IF(NTYPE.NE.5) THEN
       F12=  XJACM(1,2)/DEFRA
       F12I= XJACI(1,2)*DEFRA
       F12B= XJACN(1,2)
       F12BI=XJANI(1,2)
       F21=  XJACM(2,1)/DEFRA
       F21I= XJACI(2,1)*DEFRA
       F21B= XJACN(2,1)
       F21BI=XJANI(2,1)
       F22=  XJACM(2,2)/DEFRA
       F22I= XJACI(2,2)*DEFRA
       F22B= XJACN(2,2)
       F22BI=XJANI(2,2)
       AUX21=F11BI*F12BI+F21BI*F22BI
       E11B1=0.5D0*F11*F11+0.5D0*F21*F21-F11*F21*AUX21
       E22B1=0.5D0*F12*F12+0.5D0*F22*F22-F12*F22*AUX21
       E12B1=F11*F12+F21*F22-(F11*F22+F12*F21)*AUX21
       A3=XJA3M/(DETJM/DEFRA)
       A3B=XJA3N/DETJN
       IF(NTYPE.EQ.3) A3=XJA3M/DETJM                      ! axisymmetric
       AE1U1=F11-F21*AUX21
       AE1U2=F11*F21*(2.0D0*F11BI*F11BI*F12BI+
     .                2.0D0*F11BI*F21BI*F22BI-
     .                      F21BI*A3B)
       AE1U3=0.0D0
       AE1U4=F11*F21*(2.0D0*F11BI*F12BI*F21BI+
     .                      F11BI*A3B+
     .                      F21BI*F21BI*F22BI+
     .                      F12BI*F21BI*F22BI)
       BE1U1=-F11*(F11I*F11I+F21I*F21I)+
     .        F11*F11*(F11I*F11I*F11I+F11I*F21I*F21I)+
     .        F21*F21*(F11I*F12I*F12I+F11I*F22I*F22I-F22I*A3)
       BE1U2=0.0D0
       BE1U3=F11*F11*(F11I*F11I*F21I+F21I*F21I*F21I)+
     .       F21*F21*(F12I*F12I*F21I+F12I*A3+F12I*F22I*F22I)
       BE1U4=0.0D0
       AE2U1=0.0D0
       AE2U2=F12*F22*(2.0D0*F11BI*F11BI*F12BI+
     .                2.0D0*F11BI*F21BI*F22BI-
     .                      F21BI*A3B)
       AE2U3=F12-F22*AUX21
       AE2U4=F12*F22*(2.0D0*F11BI*F12BI*F21BI+
     .                      F11BI*A3B+
     .                      F21BI*F21BI*F22BI+
     .                      F12BI*F21BI*F22BI)
       BE2U1=F12*F12*(F11I*F11I*F11I+F11I*F21I*F21I)+
     .       F22*F22*(F11I*F12I*F12I+F11I*F22I*F22I-F22I*A3)
       BE2U2=0.0D0
       BE2U3=-F12*(F11I*F11I+F21I*F21I)+
     .        F12*F12*(F11I*F11I*F21I+F21I*F21I*F21I)+
     .        F22*F22*(F12I*F12I*F21I+F12I*A3+F12I*F22I*F22I)
       BE2U4=0.0D0
       AE3U1=F12-F22*AUX21
       AE3U2=(F11*F22+F12*F21)*(2.0D0*F11BI*F11BI*F12BI+
     .                          2.0D0*F11BI*F21BI*F22BI-
     .                                F21BI*A3B)
       AE3U3=F11-F21*AUX21
       AE3U4=(F11*F22+F12*F21)*(2.0D0*F11BI*F12BI*F21BI+
     .                                F11BI*A3B+
     .                                F21BI*F21BI*F22BI+
     .                                F12BI*F21BI*F22BI)
       BE3U1=      -F12*(F11I*F11I+F21I*F21I)+
     .        2.0D0*F11*F12*(F11I*F11I*F11I+F11I*F21I*F21I)+
     .        2.0D0*F21*F22*(F11I*F12I*F12I+F11I*F22I*F22I-F22I*A3)
       BE3U2=0.0D0
       BE3U3=      -F11*(F11I*F11I+F21I*F21I)+
     .        2.0D0*F11*F12*(F11I*F11I*F21I+F21I*F21I*F21I)+
     .        2.0D0*F21*F22*(F12I*F12I*F21I+F12I*A3+F12I*F22I*F22I)
       BE3U4=0.0D0
       AE1V1=F21-F11*AUX21
       AE1V2=F11*F21*(2.0D0*F11BI*F12BI*F12BI+
     .                2.0D0*F12BI*F21BI*F22BI+
     .                      F22BI*A3B)
       AE1V3=0.0D0
       AE1V4=F11*F21*(2.0D0*F11BI*F12BI*F22BI-
     .                      F12BI*A3B+
     .                2.0D0*F21BI*F22BI*F22BI)
       BE1V1=F11*F11*(F11I*F11I*F12I+F12I*F21I*F21I+F21I*A3)-
     .       F21*(F12I*F12I+F22I*F22I)+
     .       F21*F21*(F12I*F12I*F12I+F12I*F22I*F22I)
       BE1V2=0.0D0
       BE1V3=F11*F11*(F11I*F11I*F22I-F11I*A3+F21I*F21I*F22I)+
     .       F21*F21*(F12I*F12I*F22I+F22I*F22I*F22I)
       BE1V4=0.0D0
       AE2V1=0.0D0
       AE2V2=F12*F22*(2.0D0*F11BI*F12BI*F12BI+
     .                2.0D0*F12BI*F21BI*F22BI+
     .                      F22BI*A3B)
       AE2V3=F22-F12*AUX21
       AE2V4=F12*F22*(2.0D0*F11BI*F12BI*F22BI-
     .                      F12BI*A3B+
     .                2.0D0*F21BI*F22BI*F22BI)
       BE2V1=F12*F12*(F11I*F11I*F12I+F12I*F21I*F21I+F21I*A3)+
     .       F22*F22*(F12I*F12I*F12I+F12I*F22I*F22I)
       BE2V2=0.0D0
       BE2V3=F12*F12*(F11I*F11I*F22I-F11I*A3+F21I*F21I*F22I)-
     .       F22*(F12I*F12I+F22I*F22I)+
     .       F22*F22*(F12I*F12I*F22I+F22I*F22I*F22I)
       BE2V4=0.0D0
       AE3V1=F22-F12*AUX21
       AE3V2=(F11*F22+F12*F21)*(2.0D0*F11BI*F12BI*F12BI+
     .                          2.0D0*F12BI*F21BI*F22BI+
     .                                F22BI*A3B)
       AE3V3=F21-F11*AUX21
       AE3V4=(F11*F22+F12*F21)*(2.0D0*F11BI*F12BI*F22BI-
     .                                F12BI*A3B+
     .                          2.0D0*F21BI*F22BI*F22BI)
       BE3V1=     -F22*(F12I*F12I+F22I*F22I)+
     .       2.0D0*F11*F12*(F11I*F11I*F12I+F12I*F21I*F21I+F21I*A3)+
     .       2.0D0*F21*F22*(F12I*F12I*F12I+F12I*F22I*F22I)
       BE3V2=0.0D0
       BE3V3=     -F21*(F12I*F12I+F22I*F22I)+
     .       2.0D0*F11*F12*(F11I*F11I*F22I-F11I*A3+F21I*F21I*F22I)+
     .       2.0D0*F21*F22*(F12I*F12I*F22I+F22I*F22I*F22I)
       BE3V4=0.0D0



c      write(7,*) 'ae1u1,ae1u2,ae1u3,ae1u4',ae1u1,ae1u2,ae1u3,ae1u4
c      write(7,*) 'ae2u1,ae2u2,ae2u3,ae2u4',ae2u1,ae2u2,ae2u3,ae2u4
c      write(7,*) 'ae3u1,ae3u2,ae3u3,ae3u4',ae3u1,ae3u2,ae3u3,ae3u4

c      write(7,*) 'be1u1,be1u2,be1u3,be1u4',be1u1,be1u2,be1u3,be1u4
c      write(7,*) 'be2u1,be2u2,be2u3,be2u4',be2u1,be2u2,be2u3,be2u4
c      write(7,*) 'be3u1,be3u2,be3u3,be3u4',be3u1,be3u2,be3u3,be3u4

c      write(7,*) 'ae1v1,ae1v2,ae1v3,ae1v4',ae1v1,ae1v2,ae1v3,ae1v4
c      write(7,*) 'ae2v1,ae2v2,ae2v3,ae2v4',ae2v1,ae2v2,ae2v3,ae2v4
c      write(7,*) 'ae3v1,ae3v2,ae3v3,ae3v4',ae3v1,ae3v2,ae3v3,ae3v4

c      write(7,*) 'be1v1,be1v2,be1v3,be1v4',be1v1,be1v2,be1v3,be1v4
c      write(7,*) 'be2v1,be2v2,be2v3,be2v4',be2v1,be2v2,be2v3,be2v4
c      write(7,*) 'be3v1,be3v2,be3v3,be3v4',be3v1,be3v2,be3v3,be3v4

c      write(7,*) '-----------'




c      DXNU1= XJACN(2,2)
c      DXNU2=-XJACN(2,1)
c      DXMU1= XJACM(2,2)/DEFRA
c      DXMU2=-XJACM(2,1)/DEFRA
c      DXNV1=-XJACN(1,2)
c      DXNV2= XJACN(1,1)
c      DXMV1=-XJACM(1,2)/DEFRA
c      DXMV2= XJACM(1,1)/DEFRA
       IF(NTYPE.EQ.3) THEN
        F33=  XJA3M/DEFRA
        F33I= XJA3I*DEFRA
        F33B= XJA3N
        F33BI=XJ3NI
        E33B1=0.5D0*F33*F33

        AE1U7=0.                           ! axisymmetric coefficients
        AE1U8=0.
        BE1U7=0.
        BE1U8=0.

        AE2U7=0.
        AE2U8=0.
        BE2U7=0.
        BE2U8=0.

        AE3U7=0.
        AE3U8=0.
        BE3U7=0.
        BE3U8=0.

        AE4U7=XJA3M/DEFRA
        AE4U8=0.0D0
        BE4U7=0.0D0
        BE4U8=0.0D0

        AE1V7=0.
        AE1V8=0.
        BE1V7=0.
        BE1V8=0.

        AE2V7=0.
        AE2V8=0.
        BE2V7=0.
        BE2V8=0.

        AE3V7=0.
        AE3V8=0.
        BE3V7=0.
        BE3V8=0.

        AE4V7=0.0D0
        AE4V8=0.0D0
        BE4V7=0.0D0
        BE4V8=0.0D0

c       DXNU1=DXNU1*XJA3N
c       DXNU2=DXNU2*XJA3N
c       DXMU1=DXMU1*XJA3M/DEFRA
c       DXMU2=DXMU2*XJA3M/DEFRA
c       DXNV1=DXNV1*XJA3N
c       DXNV2=DXNV2*XJA3N
c       DXMV1=DXMV1*XJA3M/DEFRA
c       DXMV2=DXMV2*XJA3M/DEFRA
c       DE2DN=DETJN/XJA3N
c       DE2DM=(DETJM/DEFRA)/(XJA3M/DEFRA)
       ENDIF
       IF(NTYPE.EQ.4) THEN
        call runend('error: ndime=3 not implemented in bbmatzls.f')
c       E11B1=E11B1+
c       E22B1=E22B1+
c       E12B1=E12B1+
c       E33B1=E33B1+
c       E13B1=
c       E23B1=

c       AE1V1=
c       AE1V2=
c       AE1V3=
c       AE1V4=
c       BE1V1=
c       BE1V2=0.0D0
c       BE1V3=
c       BE1V4=0.0D0

        DXNU1= (XJACN(2,2)*XJACN(3,3)-XJACN(3,2)*XJACN(2,3))
        DXNU2=-(XJACN(2,1)*XJACN(3,3)-XJACN(3,1)*XJACN(2,3))
        DXNU3= (XJACN(2,1)*XJACN(3,2)-XJACN(3,1)*XJACN(2,2))
        DXMU1= (XJACM(2,2)*XJACM(3,3)-XJACM(3,2)*XJACM(2,3))/DEFR2
        DXMU2=-(XJACM(2,1)*XJACM(3,3)-XJACM(3,1)*XJACM(2,3))/DEFR2
        DXMU3= (XJACM(2,1)*XJACM(3,2)-XJACM(3,1)*XJACM(2,2))/DEFR2
        DXNV1=-(XJACN(1,2)*XJACN(3,3)-XJACN(3,2)*XJACN(1,3))
        DXNV2= (XJACN(1,1)*XJACN(3,3)-XJACN(1,3)*XJACN(3,1))
        DXNV3=-(XJACN(1,1)*XJACN(3,2)-XJACN(3,1)*XJACN(1,2))
        DXMV1=-(XJACM(1,2)*XJACM(3,3)-XJACM(3,2)*XJACM(1,3))/DEFR2
        DXMV2= (XJACM(1,1)*XJACM(3,3)-XJACM(1,3)*XJACM(3,1))/DEFR2
        DXMV3=-(XJACM(1,1)*XJACM(3,2)-XJACM(3,1)*XJACM(1,2))/DEFR2
        DXNW1= (XJACN(1,2)*XJACN(3,3)-XJACN(2,2)*XJACN(1,3))
        DXNW2=-(XJACN(1,1)*XJACN(2,3)-XJACN(2,1)*XJACN(1,3))
        DXNW3= (XJACN(1,1)*XJACN(2,2)-XJACN(2,1)*XJACN(1,2))
        DXMW1= (XJACM(1,2)*XJACM(3,3)-XJACM(2,2)*XJACM(1,3))/DEFR2
        DXMW2=-(XJACM(1,1)*XJACM(2,3)-XJACM(2,1)*XJACM(1,3))/DEFR2
        DXMW3= (XJACM(1,1)*XJACM(2,2)-XJACM(2,1)*XJACM(1,2))/DEFR2
       ENDIF           ! ntype.eq.4
      ENDIF            ! ntype.ne.5
C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        AU=(DXNU1*CMEAN(1,INODE))/DETJN-
     .     (DXMU1*CARTD(1,INODE))/(DETJM/DEFRA)

        au=0.0d0

        BMATX(1,INODE)=DEFR2*
     .                (AE1U1*CARTD(1,INODE)+AE1U2*CMEAN(1,INODE)+
     .                 B23*AU*E11B1)+
     .                 BE1U1*CARTD(1,INODE)+BE1U2*CMEAN(1,INODE)
       ENDDO
      ENDIF
C
C**** BUILD UP THE BMATX
C
      LGASH=0
      DO 30 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
C
      AU=(DXNU1*CMEAN(1,INODE)+DXNU2*CMEAN(2,INODE))/DETJN-
     .   (DXMU1*CARTD(1,INODE)+DXMU2*CARTD(2,INODE))/(DETJM/DEFRA)
      AV=(DXNV1*CMEAN(1,INODE)+DXNV2*CMEAN(2,INODE))/DETJN-
     .   (DXMV1*CARTD(1,INODE)+DXMV2*CARTD(2,INODE))/(DETJM/DEFRA)

        au=0.0d0
        av=0.0d0

C
      BMATX(1,MGASH)=DEFR2*
     .              (AE1U1*CARTD(1,INODE)+AE1U2*CMEAN(1,INODE)+
     .               AE1U3*CARTD(2,INODE)+AE1U4*CMEAN(2,INODE)+
     .               B23*AU*E11B1)+
     .               BE1U1*CARTD(1,INODE)+BE1U2*CMEAN(1,INODE)+
     .               BE1U3*CARTD(2,INODE)+BE1U4*CMEAN(2,INODE)
      BMATX(1,NGASH)=DEFR2*
     .              (AE1V1*CARTD(1,INODE)+AE1V2*CMEAN(1,INODE)+
     .               AE1V3*CARTD(2,INODE)+AE1V4*CMEAN(2,INODE)+
     .               B23*AV*E11B1)+
     .               BE1V1*CARTD(1,INODE)+BE1V2*CMEAN(1,INODE)+
     .               BE1V3*CARTD(2,INODE)+BE1V4*CMEAN(2,INODE)
      BMATX(2,MGASH)=DEFR2*
     .              (AE2U1*CARTD(1,INODE)+AE2U2*CMEAN(1,INODE)+
     .               AE2U3*CARTD(2,INODE)+AE2U4*CMEAN(2,INODE)+
     .               B23*AU*E22B1)+
     .               BE2U1*CARTD(1,INODE)+BE2U2*CMEAN(1,INODE)+
     .               BE2U3*CARTD(2,INODE)+BE2U4*CMEAN(2,INODE)
      BMATX(2,NGASH)=DEFR2*
     .              (AE2V1*CARTD(1,INODE)+AE2V2*CMEAN(1,INODE)+
     .               AE2V3*CARTD(2,INODE)+AE2V4*CMEAN(2,INODE)+
     .               B23*AV*E22B1)+
     .               BE2V1*CARTD(1,INODE)+BE2V2*CMEAN(1,INODE)+
     .               BE2V3*CARTD(2,INODE)+BE2V4*CMEAN(2,INODE)
      BMATX(3,MGASH)=DEFR2*
     .              (AE3U1*CARTD(1,INODE)+AE3U2*CMEAN(1,INODE)+
     .               AE3U3*CARTD(2,INODE)+AE3U4*CMEAN(2,INODE)+
     .               B23*AU*E12B1)+
     .               BE3U1*CARTD(1,INODE)+BE3U2*CMEAN(1,INODE)+
     .               BE3U3*CARTD(2,INODE)+BE3U4*CMEAN(2,INODE)
      BMATX(3,NGASH)=DEFR2*
     .              (AE3V1*CARTD(1,INODE)+AE3V2*CMEAN(1,INODE)+
     .               AE3V3*CARTD(2,INODE)+AE3V4*CMEAN(2,INODE)+
     .               B23*AV*E12B1)+
     .               BE3V1*CARTD(1,INODE)+BE3V2*CMEAN(1,INODE)+
     .               BE3V3*CARTD(2,INODE)+BE3V4*CMEAN(2,INODE)
C
C**** COMPLETE FOR THE AXI-SYMMETRIC CASE
C
      IF(NTYPE.NE.3) GOTO 40
C
      AU33=DE2DN*CMEAN(3,INODE)/DETJN-
     .    (DE2DM*SHAPE(INODE)/GPCOD)/(DETJM/DEFRA)

      au33=0.0d0

C
      BMATX(1,MGASH)=BMATX(1,MGASH)+
     .               DEFR2*
     .              (AE1U7*SHAPE(INODE)/GPCOD+AE1U8*CMEAN(3,INODE)+
     .               B23*AU33*E11B1)+
     .               BE1U7*SHAPE(INODE)/GPCOD+BE1U8*CMEAN(3,INODE)
      BMATX(2,MGASH)=BMATX(2,MGASH)+
     .               DEFR2*
     .              (AE2U7*SHAPE(INODE)/GPCOD+AE2U8*CMEAN(3,INODE)+
     .               B23*AU33*E22B1)+
     .               BE2U7*SHAPE(INODE)/GPCOD+BE2U8*CMEAN(3,INODE)
      BMATX(3,MGASH)=BMATX(3,MGASH)+
     .               DEFR2*
     .              (AE3U7*SHAPE(INODE)/GPCOD+AE3U8*CMEAN(3,INODE)+
     .               B23*AU33*E12B1)+
     .               BE3U7*SHAPE(INODE)/GPCOD+BE3U8*CMEAN(3,INODE)
C
      BMATX(4,MGASH)=DEFR2*
     .              (AE4U7*SHAPE(INODE)/GPCOD+AE4U8*CMEAN(3,INODE)+
     .               B23*AU*E33B1+
     .               B23*AU33*E33B1)+
     .               BE4U7*SHAPE(INODE)/GPCOD+BE4U8*CMEAN(3,INODE)
      BMATX(4,NGASH)=DEFR2*
     .              (AE4V7*SHAPE(INODE)/GPCOD+AE4V8*CMEAN(3,INODE)+
     .               B23*AV*E33B1)+
     .               BE4V7*SHAPE(INODE)/GPCOD+BE4V8*CMEAN(3,INODE)
C
C**** COMPLETE FOR THE 3-D CASE
C
   40 IF(NTYPE.NE.4) GOTO 30
C
      LGASH=NGASH+1
C
      AUX=(DXNU3*CMEAN(3,INODE))/DETJN-
     .    (DXMU3*CARTD(3,INODE))/(DETJM/DEFRA)
      AVX=(DXNV3*CMEAN(3,INODE))/DETJN-
     .    (DXMV3*CARTD(3,INODE))/(DETJM/DEFRA)
      AU=AU+AUX
      AV=AV+AVX
      AW=(DXNW1*CMEAN(1,INODE)+DXNW2*CMEAN(2,INODE)+
     .    DXNW3*CMEAN(3,INODE))/DETJN-
     .   (DXMW1*CARTD(1,INODE)+DXMW2*CARTD(2,INODE)+
     .    DXMW3*CARTD(3,INODE))/(DETJM/DEFRA)
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+1.0/3.0*RCGTT(1)*AUX*DEFR2
      BMATX(1,NGASH)=BMATX(1,NGASH)+1.0/3.0*RCGTT(1)*AVX*DEFR2
      BMATX(2,MGASH)=BMATX(2,MGASH)+1.0/3.0*RCGTT(2)*AUX*DEFR2
      BMATX(2,NGASH)=BMATX(2,NGASH)+1.0/3.0*RCGTT(2)*AVX*DEFR2
      BMATX(3,MGASH)=BMATX(3,MGASH)+2.0/3.0*RCGTT(3)*AUX*DEFR2
      BMATX(3,NGASH)=BMATX(3,NGASH)+2.0/3.0*RCGTT(3)*AVX*DEFR2
C
      BMATX(1,LGASH)=(CARTD(1,INODE)*XJACM(3,1)/DEFRA+
     .                                        1.0/3.0*RCGTT(1)*AW)*DEFR2
      BMATX(2,LGASH)=(CARTD(2,INODE)*XJACM(3,2)/DEFRA+
     .                                        1.0/3.0*RCGTT(2)*AW)*DEFR2
      BMATX(3,LGASH)=(CARTD(2,INODE)*XJACM(3,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(3,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(3)*AW)*DEFR2
      BMATX(4,MGASH)=(CARTD(3,INODE)*XJACM(1,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AU)*DEFR2
      BMATX(4,NGASH)=(CARTD(3,INODE)*XJACM(2,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AV)*DEFR2
      BMATX(4,LGASH)=(CARTD(3,INODE)*XJACM(3,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AW)*DEFR2
      BMATX(5,MGASH)=(CARTD(1,INODE)*XJACM(1,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(1,1)/DEFRA+
     .                                        2.0/3.0*RCGTT(5)*AU)*DEFR2
      BMATX(5,NGASH)=(CARTD(1,INODE)*XJACM(2,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(2,1)/DEFRA+
     .                                        2.0/3.0*RCGTT(5)*AV)*DEFR2
      BMATX(5,LGASH)=(CARTD(1,INODE)*XJACM(3,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(3,1)/DEFRA+
     .                                        2.0/3.0*RCGTT(5)*AW)*DEFR2
      BMATX(6,MGASH)=(CARTD(2,INODE)*XJACM(1,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(1,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(6)*AU)*DEFR2
      BMATX(6,NGASH)=(CARTD(2,INODE)*XJACM(2,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(2,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(6)*AV)*DEFR2
      BMATX(6,LGASH)=(CARTD(2,INODE)*XJACM(3,3)/DEFRA+
     .                CARTD(3,INODE)*XJACM(3,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(6)*AW)*DEFR2
C
   30 CONTINUE
C


c     return



C
C**** PRELIMINARY COMPUTATIONS
C
      DEFR2=DEFRA*DEFRA
C
      RCGTT(1)=XJACM(1,1)*XJACM(1,1)/DEFR2   ! right Cauchy-Green tensor
      DXNU1=1.0
      DXMU1=1.0
      IF(NTYPE.NE.5) THEN
       RCGTT(1)=(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1))/DEFR2
       RCGTT(2)=(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2))/DEFR2
       RCGTT(3)=(XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2))/DEFR2
       RCGTT(4)=1.0
       DXNU1= XJACN(2,2)
       DXNU2=-XJACN(2,1)
       DXMU1= XJACM(2,2)/DEFRA
       DXMU2=-XJACM(2,1)/DEFRA
       DXNV1=-XJACN(1,2)
       DXNV2= XJACN(1,1)
       DXMV1=-XJACM(1,2)/DEFRA
       DXMV2= XJACM(1,1)/DEFRA
       IF(NTYPE.EQ.3) THEN
        RCGTT(4)=XJA3M*XJA3M/DEFR2
        DXNU1=DXNU1*XJA3N
        DXNU2=DXNU2*XJA3N
        DXMU1=DXMU1*XJA3M/DEFRA
        DXMU2=DXMU2*XJA3M/DEFRA
        DXNV1=DXNV1*XJA3N
        DXNV2=DXNV2*XJA3N
        DXMV1=DXMV1*XJA3M/DEFRA
        DXMV2=DXMV2*XJA3M/DEFRA
        DE2DN=DETJN/XJA3N
        DE2DM=(DETJM/DEFRA)/(XJA3M/DEFRA)
       ENDIF
       IF(NTYPE.EQ.4) THEN
        RCGTT(1)=(XJACM(1,1)*XJACM(1,1)+XJACM(2,1)*XJACM(2,1)+
     .            XJACM(3,1)*XJACM(3,1))/DEFR2
        RCGTT(2)=(XJACM(1,2)*XJACM(1,2)+XJACM(2,2)*XJACM(2,2)+
     .            XJACM(3,2)*XJACM(3,2))/DEFR2
        RCGTT(3)=(XJACM(1,1)*XJACM(1,2)+XJACM(2,1)*XJACM(2,2)+
     .            XJACM(3,1)*XJACM(3,2))/DEFR2
        RCGTT(4)=(XJACM(1,3)*XJACM(1,3)+XJACM(2,3)*XJACM(2,3)+
     .            XJACM(3,3)*XJACM(3,3))/DEFR2
        RCGTT(5)=(XJACM(1,1)*XJACM(1,3)+XJACM(2,1)*XJACM(2,3)+
     .            XJACM(3,1)*XJACM(3,3))/DEFR2
        RCGTT(6)=(XJACM(1,2)*XJACM(1,3)+XJACM(2,2)*XJACM(2,3)+
     .            XJACM(3,2)*XJACM(3,3))/DEFR2
        DXNU1= (XJACN(2,2)*XJACN(3,3)-XJACN(3,2)*XJACN(2,3))
        DXNU2=-(XJACN(2,1)*XJACN(3,3)-XJACN(3,1)*XJACN(2,3))
        DXNU3= (XJACN(2,1)*XJACN(3,2)-XJACN(3,1)*XJACN(2,2))
        DXMU1= (XJACM(2,2)*XJACM(3,3)-XJACM(3,2)*XJACM(2,3))/DEFR2
        DXMU2=-(XJACM(2,1)*XJACM(3,3)-XJACM(3,1)*XJACM(2,3))/DEFR2
        DXMU3= (XJACM(2,1)*XJACM(3,2)-XJACM(3,1)*XJACM(2,2))/DEFR2
        DXNV1=-(XJACN(1,2)*XJACN(3,3)-XJACN(3,2)*XJACN(1,3))
        DXNV2= (XJACN(1,1)*XJACN(3,3)-XJACN(1,3)*XJACN(3,1))
        DXNV3=-(XJACN(1,1)*XJACN(3,2)-XJACN(3,1)*XJACN(1,2))
        DXMV1=-(XJACM(1,2)*XJACM(3,3)-XJACM(3,2)*XJACM(1,3))/DEFR2
        DXMV2= (XJACM(1,1)*XJACM(3,3)-XJACM(1,3)*XJACM(3,1))/DEFR2
        DXMV3=-(XJACM(1,1)*XJACM(3,2)-XJACM(3,1)*XJACM(1,2))/DEFR2
        DXNW1= (XJACN(1,2)*XJACN(3,3)-XJACN(2,2)*XJACN(1,3))
        DXNW2=-(XJACN(1,1)*XJACN(2,3)-XJACN(2,1)*XJACN(1,3))
        DXNW3= (XJACN(1,1)*XJACN(2,2)-XJACN(2,1)*XJACN(1,2))
        DXMW1= (XJACM(1,2)*XJACM(3,3)-XJACM(2,2)*XJACM(1,3))/DEFR2
        DXMW2=-(XJACM(1,1)*XJACM(2,3)-XJACM(2,1)*XJACM(1,3))/DEFR2
        DXMW3= (XJACM(1,1)*XJACM(2,2)-XJACM(2,1)*XJACM(1,2))/DEFR2
       ENDIF           ! ntype.eq.4
      ENDIF            ! ntype.ne.5
C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        AU=(DXNU1*CMEAN(1,INODE))/DETJN-
     .     (DXMU1*CARTD(1,INODE))/(DETJM/DEFRA)
        BMATX(1,INODE)=(CARTD(1,INODE)*XJACM(1,1)/DEFRA+
     .                                       1.0/3.0*RCGTT(1)*AU)*DEFRA2
       ENDDO
      ENDIF
C
C**** CALCULATE CONSTANT TERM FOR AXI-SYMMETRIC CASE IN FINITE
C     DEFORMATIONS
C
      FMULT=1.0
      IF(NTYPE.NE.3) GOTO 110
      FMULT=0.0
      DO 120 JNODE=1,NNODE
  120 FMULT=FMULT+ELDIS(1,JNODE)*SHAPE(JNODE)
      FMULT=FMULT/GPCOD                              ! = XJA3M
  110 CONTINUE
C
C**** BUILD UP THE BMATX
C
      LGASH=0
      DO 130 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
C
      AU=(DXNU1*CMEAN(1,INODE)+DXNU2*CMEAN(2,INODE))/DETJN-
     .   (DXMU1*CARTD(1,INODE)+DXMU2*CARTD(2,INODE))/(DETJM/DEFRA)
      AV=(DXNV1*CMEAN(1,INODE)+DXNV2*CMEAN(2,INODE))/DETJN-
     .   (DXMV1*CARTD(1,INODE)+DXMV2*CARTD(2,INODE))/(DETJM/DEFRA)
C
      BMATX(1,MGASH)=(CARTD(1,INODE)*XJACM(1,1)/DEFRA+
     .                                        1.0/3.0*RCGTT(1)*AU)*DEFR2
      BMATX(1,NGASH)=(CARTD(1,INODE)*XJACM(2,1)/DEFRA+
     .                                        1.0/3.0*RCGTT(1)*AV)*DEFR2
      BMATX(2,MGASH)=(CARTD(2,INODE)*XJACM(1,2)/DEFRA+
     .                                        1.0/3.0*RCGTT(2)*AU)*DEFR2
      BMATX(2,NGASH)=(CARTD(2,INODE)*XJACM(2,2)/DEFRA+
     .                                        1.0/3.0*RCGTT(2)*AV)*DEFR2
c     BMATX(3,MGASH)=CMEAN(2,INODE)*XJACN(1,1)+
c    .               CMEAN(1,INODE)*XJACN(1,2)
c     BMATX(3,NGASH)=CMEAN(2,INODE)*XJACN(2,1)+
c    .               CMEAN(1,INODE)*XJACN(2,2)
      BMATX(3,MGASH)=(CARTD(2,INODE)*XJACM(1,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(1,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(3)*AU)*DEFR2
      BMATX(3,NGASH)=(CARTD(2,INODE)*XJACM(2,1)/DEFRA+
     .                CARTD(1,INODE)*XJACM(2,2)/DEFRA+
     .                                        2.0/3.0*RCGTT(3)*AV)*DEFR2
C
C**** COMPLETE FOR THE AXI-SYMMETRIC CASE
C
      IF(NTYPE.NE.3) GOTO 140
C
      AU33=DE2DN*CMEAN(3,INODE)/DETJN-
     .    (DE2DM*SHAPE(INODE)/GPCOD)/(DETJM/DEFRA)
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+1.0/3.0*RCGTT(1)*AU33*DEFR2
      BMATX(2,MGASH)=BMATX(2,MGASH)+1.0/3.0*RCGTT(2)*AU33*DEFR2
c
      BMATX(3,MGASH)=BMATX(3,MGASH)+2.0/3.0*RCGTT(3)*AU33*DEFR2
C
      BMATX(4,MGASH)=(SHAPE(INODE)*FMULT/GPCOD+
     .                                      1.0/3.0*RCGTT(4)*AU+
     .                                      1.0/3.0*RCGTT(4)*AU33)*DEFR2
      BMATX(4,NGASH)=1.0/3.0*RCGTT(4)*AV*DEFR2
C
C**** COMPLETE FOR THE 3-D CASE
C
  140 IF(NTYPE.NE.4) GOTO 130
C
      LGASH=NGASH+1
C
      AUX=(DXNU3*CMEAN(3,INODE))/DETJN-
     .    (DXMU3*CARTD(3,INODE))/(DETJM/DEFRA)
      AVX=(DXNV3*CMEAN(3,INODE))/DETJN-
     .    (DXMV3*CARTD(3,INODE))/(DETJM/DEFRA)
      AU=AU+AUX
      AV=AV+AVX
      AW=(DXNW1*CMEAN(1,INODE)+DXNW2*CMEAN(2,INODE)+
     .    DXNW3*CMEAN(3,INODE))/DETJN-
     .   (DXMW1*CARTD(1,INODE)+DXMW2*CARTD(2,INODE)+
     .    DXMW3*CARTD(3,INODE))/(DETJM/DEFRA)
C
      BMATX(1,MGASH)=BMATX(1,MGASH)+1.0/3.0*RCGTT(1)*AUX*DEFR2
      BMATX(1,NGASH)=BMATX(1,NGASH)+1.0/3.0*RCGTT(1)*AVX*DEFR2
      BMATX(2,MGASH)=BMATX(2,MGASH)+1.0/3.0*RCGTT(2)*AUX*DEFR2
      BMATX(2,NGASH)=BMATX(2,NGASH)+1.0/3.0*RCGTT(2)*AVX*DEFR2
C
      BMATX(1,LGASH)=(CARTD(1,INODE)*XJACM(3,1)/DEFRA+
     .                                        1.0/3.0*RCGTT(1)*AW)*DEFR2
      BMATX(2,LGASH)=(CARTD(2,INODE)*XJACM(3,2)/DEFRA+
     .                                        1.0/3.0*RCGTT(2)*AW)*DEFR2
      BMATX(3,LGASH)=CMEAN(2,INODE)*XJACN(3,1)+
     .               CMEAN(1,INODE)*XJACN(3,2)
      BMATX(4,MGASH)=(CARTD(3,INODE)*XJACM(1,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AU)*DEFR2
      BMATX(4,NGASH)=(CARTD(3,INODE)*XJACM(2,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AV)*DEFR2
      BMATX(4,LGASH)=(CARTD(3,INODE)*XJACM(3,3)/DEFRA+
     .                                        1.0/3.0*RCGTT(4)*AW)*DEFR2
      BMATX(5,MGASH)=CMEAN(1,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,1)
      BMATX(5,NGASH)=CMEAN(1,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,1)
      BMATX(5,LGASH)=CMEAN(1,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,1)
      BMATX(6,MGASH)=CMEAN(2,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,2)
      BMATX(6,NGASH)=CMEAN(2,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,2)
      BMATX(6,LGASH)=CMEAN(2,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,2)
C
  130 CONTINUE
C
































C
C**** 1-D CASE
C
      IF(NTYPE.EQ.5) THEN
       DO INODE=1,NNODE
        BMATX(1,INODE)=CARTD(1,INODE)*XJACM(1,1)
       ENDDO
      ENDIF
C
C**** CALCULATE CONSTANT TERM FOR AXI-SYMMETRIC CASE IN FINITE
C     DEFORMATIONS
C
      FMULT=1.0D0
      IF(NTYPE.NE.3) GOTO 10
      FMULT=0.0D0
      DO 20 JNODE=1,NNODE
   20 FMULT=FMULT+ELDIS(1,JNODE)*SHAPE(JNODE)
      FMULT=FMULT/GPCOD
   10 CONTINUE
C
C**** BUILD UP THE BMATX
C
      LGASH=0
      DO 230 INODE=1,NNODE
      MGASH=LGASH+1
      NGASH=MGASH+1
      LGASH=NGASH
      BMATX(1,MGASH)=CARTD(1,INODE)*XJACM(1,1)
      BMATX(1,NGASH)=CARTD(1,INODE)*XJACM(2,1)
      BMATX(2,MGASH)=CARTD(2,INODE)*XJACM(1,2)
      BMATX(2,NGASH)=CARTD(2,INODE)*XJACM(2,2)

c     BMATX(3,MGASH)=CMEAN(2,INODE)*XJACN(1,1)+
c    .               CMEAN(1,INODE)*XJACN(1,2)
c     BMATX(3,NGASH)=CMEAN(2,INODE)*XJACN(2,1)+
c    .               CMEAN(1,INODE)*XJACN(2,2)

      BMATX(3,MGASH)=CARTD(2,INODE)*XJACM(1,1)+
     .               CARTD(1,INODE)*XJACM(1,2)
      BMATX(3,NGASH)=CARTD(2,INODE)*XJACM(2,1)+
     .               CARTD(1,INODE)*XJACM(2,2)

C
C**** COMPLETE FOR THE AXI-SYMMETRIC CASE
C
      IF(NTYPE.NE.3) GOTO 240
C
      BMATX(4,MGASH)=SHAPE(INODE)*FMULT/GPCOD
      BMATX(4,NGASH)=0.0D0
C
C**** COMPLETE FOR THE 3-D CASE
C
  240 IF(NTYPE.NE.4) GOTO 230
C
      LGASH=NGASH+1
      BMATX(1,LGASH)=CARTD(1,INODE)*XJACM(3,1)
      BMATX(2,LGASH)=CARTD(2,INODE)*XJACM(3,2)
      BMATX(3,LGASH)=CMEAN(2,INODE)*XJACN(3,1)+
     .               CMEAN(1,INODE)*XJACN(3,2)
      BMATX(4,MGASH)=CARTD(3,INODE)*XJACM(1,3)
      BMATX(4,NGASH)=CARTD(3,INODE)*XJACM(2,3)
      BMATX(4,LGASH)=CARTD(3,INODE)*XJACM(3,3)
      BMATX(5,MGASH)=CMEAN(1,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,1)
      BMATX(5,NGASH)=CMEAN(1,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,1)
      BMATX(5,LGASH)=CMEAN(1,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,1)
      BMATX(6,MGASH)=CMEAN(2,INODE)*XJACN(1,3)+
     .               CMEAN(3,INODE)*XJACN(1,2)
      BMATX(6,NGASH)=CMEAN(2,INODE)*XJACN(2,3)+
     .               CMEAN(3,INODE)*XJACN(2,2)
      BMATX(6,LGASH)=CMEAN(2,INODE)*XJACN(3,3)+
     .               CMEAN(3,INODE)*XJACN(3,2)
C
  230 CONTINUE
C
      RETURN
      END
