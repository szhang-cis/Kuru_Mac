      SUBROUTINE FOURTH(DMATX,XJACM,XJACI,XJA3M,XJA3I,DETJM,
     .                  POISC1,POISC2,POISC3,D1PLS,D2PLS,D3PLS,
     .                  BULKM,DISTM,CPOL1,CPOL2,CPOL3,
     .                  EPOL1,EPOL2,EPOL3,
     .                  CPOL4,CPOL5,CPOL6,CPOL7,CPOL8,CPENI,
     .                  CPOL9,CPOL10,
     .                  CFIB1,CFIB2,CFIB3,CFIB4,CFIB5,
     .                  CEMIN,CEMAX,CEETA,
     .                  ALAMB,ILAEQ,FVOL1,
     .                  RCGTT,RCGTI,TRACC,DETJC,SEINC,
     .                  BULKA,NAUXI,UNOMA,FACTJ,DCOFA,UNORA,
     .                  STILD,ATILD,VANIS,
     .                  STRAP,RCGPI,FINET,TRABE,DCINV,EBASE)
#ifndef restricted
C***********************************************************************
C
C**** THIS ROUTINE TRANSFORMS THE CONSTITUTIVE TENSOR FOR STANDARD
C     MODELS OR COMPUTES THE CONSTITUTIVE TENSOR FOR NON-STANDARD ONES
C
C     Notes:
C
C     The transformation is valid for general (isotropic or
C     anisotropic) symmetric constitutive tensors
C
C     The transformation is only valid for LARGE=1,2 and IFREN=2,3,5,6
C
C     The computation is only valid for LARGE=1,2 and IFREN=7,8,
C                                                51,52,53,54,55,56,57,58
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION XJACM(NDIME,*), XJACI(NDIME,*)
      DIMENSION DMATX(NAUXI,*), RCGTI(*),       RCGTT(*), ALAMB(*),
     .          UNOMA(*),       VANIS(NANIV,*)
      DIMENSION DMATT(3,3,3,3), AUXIL(3,3,3,3)
      DIMENSION DCOFA(6,*),     UNORA(6,*),     STILD(*), BTILD(6),
     .          AMAA1(6),       AMAA2(6),       AMAA3(6),
     .          ADET1(6),       ADET2(6),       ADET3(6),
     .          ATHE1(6,6),     ATHE2(6,6),     ATHE3(6,6),
     .          ANIS1(6),       ANIS2(6),       ANIS3(6)
      DIMENSION STRAP(*),       RCGPI(*),       FINET(*),
     .          DCINV(6,*),     EBASE(*)
      DIMENSION CPOL9(*), CPOL10(*)
      DIMENSION VANIX(6), S0XS0(6,6)
C
      IF(LARGE.EQ.3) RETURN
      IF(IFREN.EQ.1.OR.IFREN.EQ.4.OR.IFREN.EQ.9.OR.IFREN.EQ.10) RETURN
C
      ISTAN=2
      IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .   IFREN.EQ.5.OR.IFREN.EQ.6) ISTAN=1
C
      IF(ISTAN.EQ.1) THEN          ! standard models
       IF(IFREN.EQ.2.OR.IFREN.EQ.5) DETJX=1.0D0
       IF(IFREN.EQ.3.OR.IFREN.EQ.6) DETJX=DETJM
C
C**** BUILDS CONSTITUTIVE TENSOR
C
C     Notes:
C
C     DMATT for I,J,K or L equals 3 (coorresponding to DMATX for I or J
C     equals 4) must be transformed for planes stress, planes strain &
C     axisymmetric problems.
C
C     For plane stress, the transformation of DMATX for I or J equals 3
C     is useful to compute the out-of-plane strain (constants POISCi).
C
C     For plane strain, the transformation of DMATX for I or J equals 3
C     is useful to compute the out-of-plane stress.
C
C     For axisymmetry, DETJM=DVOLI/DVOLU (this form of DETJM includes
C     the radius variation and it is used for all the components of
C     DMATX, where DVOLI & DVOLU are the differential volumes at the
C     spatial and material configurations, respectively. Other
C     equivalent form is DETJM=DETJM_(from PROMA2 & INVMTX) * r / r_o,
C     where r & r_o are the radii at the spatial and material
C     configurations, respectively. 
C
       DMATT(1,1,1,1)=DMATX(1,1)
       IF(NTYPE.NE.5) THEN
        DMATT(1,1,1,2)=DMATX(1,3)
        DMATT(1,1,2,1)=DMATX(1,3)
        DMATT(1,1,2,2)=DMATX(1,2)
        DMATT(1,2,1,1)=DMATX(3,1)
        DMATT(1,2,1,2)=DMATX(3,3)
        DMATT(1,2,2,1)=DMATX(3,3)
        DMATT(1,2,2,2)=DMATX(3,2)
        DMATT(2,1,1,1)=DMATX(3,1)
        DMATT(2,1,1,2)=DMATX(3,3)
        DMATT(2,1,2,1)=DMATX(3,3)
        DMATT(2,1,2,2)=DMATX(3,2)
        DMATT(2,2,1,1)=DMATX(2,1)
        DMATT(2,2,1,2)=DMATX(2,3)
        DMATT(2,2,2,1)=DMATX(2,3)
        DMATT(2,2,2,2)=DMATX(2,2)
        IF(NTYPE.NE.1) THEN
         DMATT(1,1,3,3)=DMATX(1,4)
         DMATT(1,2,3,3)=DMATX(3,4)
         DMATT(2,1,3,3)=DMATX(3,4)
         DMATT(2,2,3,3)=DMATX(2,4)
         DMATT(3,3,1,1)=DMATX(4,1)
         DMATT(3,3,1,2)=DMATX(4,3)
         DMATT(3,3,2,1)=DMATX(4,3)
         DMATT(3,3,2,2)=DMATX(4,2)
         DMATT(3,3,3,3)=DMATX(4,4)
         IF(NTYPE.EQ.4) THEN
          DMATT(1,1,1,3)=DMATX(1,5)
          DMATT(1,1,2,3)=DMATX(1,6)
          DMATT(1,1,3,1)=DMATX(1,5)
          DMATT(1,1,3,2)=DMATX(1,6)
          DMATT(1,2,1,3)=DMATX(3,5)
          DMATT(1,2,2,3)=DMATX(3,6)
          DMATT(1,2,3,1)=DMATX(3,5)
          DMATT(1,2,3,2)=DMATX(3,6)
          DMATT(1,3,1,1)=DMATX(5,1)
          DMATT(1,3,1,2)=DMATX(5,3)
          DMATT(1,3,1,3)=DMATX(5,5)
          DMATT(1,3,2,1)=DMATX(5,3)
          DMATT(1,3,2,2)=DMATX(5,2)
          DMATT(1,3,2,3)=DMATX(5,6)
          DMATT(1,3,3,1)=DMATX(5,5)
          DMATT(1,3,3,2)=DMATX(5,6)
          DMATT(1,3,3,3)=DMATX(5,4)
          DMATT(2,1,1,3)=DMATX(3,5)
          DMATT(2,1,2,3)=DMATX(3,6)
          DMATT(2,1,3,1)=DMATX(3,5)
          DMATT(2,1,3,2)=DMATX(3,6)
          DMATT(2,2,1,3)=DMATX(2,5)
          DMATT(2,2,2,3)=DMATX(2,6)
          DMATT(2,2,3,1)=DMATX(2,5)
          DMATT(2,2,3,2)=DMATX(2,6)
          DMATT(2,3,1,1)=DMATX(6,1)
          DMATT(2,3,1,2)=DMATX(6,3)
          DMATT(2,3,1,3)=DMATX(6,5)
          DMATT(2,3,2,1)=DMATX(6,3)
          DMATT(2,3,2,2)=DMATX(6,2)
          DMATT(2,3,2,3)=DMATX(6,6)
          DMATT(2,3,3,1)=DMATX(6,5)
          DMATT(2,3,3,2)=DMATX(6,6)
          DMATT(2,3,3,3)=DMATX(6,4)
          DMATT(3,1,1,1)=DMATX(5,1)
          DMATT(3,1,1,2)=DMATX(5,3)
          DMATT(3,1,1,3)=DMATX(5,5)
          DMATT(3,1,2,1)=DMATX(5,3)
          DMATT(3,1,2,2)=DMATX(5,2)
          DMATT(3,1,2,3)=DMATX(5,6)
          DMATT(3,1,3,1)=DMATX(5,5)
          DMATT(3,1,3,2)=DMATX(5,6)
          DMATT(3,1,3,3)=DMATX(5,4)
          DMATT(3,2,1,1)=DMATX(6,1)
          DMATT(3,2,1,2)=DMATX(6,3)
          DMATT(3,2,1,3)=DMATX(6,5)
          DMATT(3,2,2,1)=DMATX(6,3)
          DMATT(3,2,2,2)=DMATX(6,2)
          DMATT(3,2,2,3)=DMATX(6,6)
          DMATT(3,2,3,1)=DMATX(6,5)
          DMATT(3,2,3,2)=DMATX(6,6)
          DMATT(3,2,3,3)=DMATX(6,4)
          DMATT(3,3,1,3)=DMATX(4,5)
          DMATT(3,3,2,3)=DMATX(4,6)
          DMATT(3,3,3,1)=DMATX(4,5)
          DMATT(3,3,3,2)=DMATX(4,6)
         ENDIF                     ! ntype.eq.4
        ENDIF                      ! ntype.ne.1
       ENDIF                       ! ntype.ne.5
C
       IF(NTYPE.EQ.1) THEN
        DMATT(1,1,3,3)=D2PLS
        DMATT(1,2,3,3)=0.0D0
        DMATT(2,1,3,3)=0.0D0
        DMATT(2,2,3,3)=D3PLS
        DMATT(3,3,1,1)=D2PLS
        DMATT(3,3,1,2)=0.0D0
        DMATT(3,3,2,1)=0.0D0
        DMATT(3,3,2,2)=D3PLS
        DMATT(3,3,3,3)=D1PLS
       ENDIF
C                                  t        -1  -1 t   -T  -T
C**** COMPUTE CONSTITUTIVE TENSOR:  C = J (F   F    C F   F  )
C                                  0               t
       DO 70 I=1,NDIME
       DO 70 J=1,NDIME
       DO 70 K=1,NDIME
       DO 70 L=1,NDIME
       AUXIL(I,J,K,L)=0.D0
       DO 70 I1=1,NDIME
       DO 70 J1=1,NDIME
       DO 70 K1=1,NDIME
       DO 70 L1=1,NDIME
   70  AUXIL(I,J,K,L)=AUXIL(I,J,K,L)+
     .                XJACI(I,I1)*XJACI(J,J1)*DMATT(I1,J1,K1,L1)*
     .                XJACI(K,K1)*XJACI(L,L1)
C
       IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3) THEN
        DO 71 K=1,NDIME
        DO 71 L=1,NDIME
        AUXIL(3,3,K,L)=0.D0
        DO 71 K1=1,NDIME
        DO 71 L1=1,NDIME
   71   AUXIL(3,3,K,L)=AUXIL(3,3,K,L)+
     .                 XJA3I*XJA3I*DMATT(3,3,K1,L1)*
     .                 XJACI(K,K1)*XJACI(L,L1)
        DO 72 I=1,NDIME
        DO 72 J=1,NDIME
        AUXIL(I,J,3,3)=0.D0
        DO 72 I1=1,NDIME
        DO 72 J1=1,NDIME
   72   AUXIL(I,J,3,3)=AUXIL(I,J,3,3)+
     .                 XJACI(I,I1)*XJACI(J,J1)*DMATT(I1,J1,3,3)*
     .                 XJA3I*XJA3I
        AUXIL(3,3,3,3)=
     .                 XJA3I*XJA3I*DMATT(3,3,3,3)*
     .                 XJA3I*XJA3I
       ENDIF
C
       DO 80 I=1,NDIME
       DO 80 J=1,NDIME
       DO 80 K=1,NDIME
       DO 80 L=1,NDIME
   80  DMATT(I,J,K,L)=AUXIL(I,J,K,L)*DETJX
C
       IF(NTYPE.EQ.1.OR.NTYPE.EQ.2.OR.NTYPE.EQ.3) THEN
        DO 81 K=1,NDIME
        DO 81 L=1,NDIME
   81   DMATT(3,3,K,L)=AUXIL(3,3,K,L)*DETJX
        DO 82 I=1,NDIME
        DO 82 J=1,NDIME
   82   DMATT(I,J,3,3)=AUXIL(I,J,3,3)*DETJX
        DMATT(3,3,3,3)=AUXIL(3,3,3,3)*DETJX
       ENDIF
C
       IF(NTYPE.EQ.1) THEN
        POISC1=-DMATT(3,3,1,1)/DMATT(3,3,3,3)
        POISC2=-DMATT(3,3,2,2)/DMATT(3,3,3,3)
        POISC3=-DMATT(3,3,1,2)/DMATT(3,3,3,3)      ! (3,3,1,2)=(3,3,2,1)
       ENDIF
C
       DMATX(1,1)=DMATT(1,1,1,1)
       IF(NTYPE.NE.5) THEN
        DMATX(1,3)=DMATT(1,1,1,2)
        DMATX(1,3)=DMATT(1,1,2,1)
        DMATX(1,2)=DMATT(1,1,2,2)
        DMATX(3,1)=DMATT(1,2,1,1)
        DMATX(3,3)=DMATT(1,2,1,2)
        DMATX(3,3)=DMATT(1,2,2,1)
        DMATX(3,2)=DMATT(1,2,2,2)
        DMATX(3,1)=DMATT(2,1,1,1)
        DMATX(3,3)=DMATT(2,1,1,2)
        DMATX(3,3)=DMATT(2,1,2,1)
        DMATX(3,2)=DMATT(2,1,2,2)
        DMATX(2,1)=DMATT(2,2,1,1)
        DMATX(2,3)=DMATT(2,2,1,2)
        DMATX(2,3)=DMATT(2,2,2,1)
        DMATX(2,2)=DMATT(2,2,2,2)
        IF(NTYPE.NE.1) THEN
         DMATX(1,4)=DMATT(1,1,3,3)
         DMATX(3,4)=DMATT(1,2,3,3)
         DMATX(3,4)=DMATT(2,1,3,3)
         DMATX(2,4)=DMATT(2,2,3,3)
         DMATX(4,1)=DMATT(3,3,1,1)
         DMATX(4,3)=DMATT(3,3,1,2)
         DMATX(4,3)=DMATT(3,3,2,1)
         DMATX(4,2)=DMATT(3,3,2,2)
         DMATX(4,4)=DMATT(3,3,3,3)
         IF(NTYPE.EQ.4) THEN
          DMATX(1,5)=DMATT(1,1,1,3)
          DMATX(1,6)=DMATT(1,1,2,3)
          DMATX(3,5)=DMATT(1,2,1,3)
          DMATX(3,6)=DMATT(1,2,2,3)
          DMATX(1,5)=DMATT(1,1,3,1)
          DMATX(1,6)=DMATT(1,1,3,2)
          DMATX(3,5)=DMATT(1,2,3,1)
          DMATX(3,6)=DMATT(1,2,3,2)
          DMATX(5,1)=DMATT(1,3,1,1)
          DMATX(5,3)=DMATT(1,3,1,2)
          DMATX(5,5)=DMATT(1,3,1,3)
          DMATX(5,3)=DMATT(1,3,2,1)
          DMATX(5,2)=DMATT(1,3,2,2)
          DMATX(5,6)=DMATT(1,3,2,3)
          DMATX(5,5)=DMATT(1,3,3,1)
          DMATX(5,6)=DMATT(1,3,3,2)
          DMATX(5,4)=DMATT(1,3,3,3)
          DMATX(3,5)=DMATT(2,1,1,3)
          DMATX(3,6)=DMATT(2,1,2,3)
          DMATX(3,5)=DMATT(2,1,3,1)
          DMATX(3,6)=DMATT(2,1,3,2)
          DMATX(2,5)=DMATT(2,2,1,3)
          DMATX(2,6)=DMATT(2,2,2,3)
          DMATX(2,5)=DMATT(2,2,3,1)
          DMATX(2,6)=DMATT(2,2,3,2)
          DMATX(6,1)=DMATT(2,3,1,1)
          DMATX(6,3)=DMATT(2,3,1,2)
          DMATX(6,5)=DMATT(2,3,1,3)
          DMATX(6,3)=DMATT(2,3,2,1)
          DMATX(6,2)=DMATT(2,3,2,2)
          DMATX(6,6)=DMATT(2,3,2,3)
          DMATX(6,5)=DMATT(2,3,3,1)
          DMATX(6,6)=DMATT(2,3,3,2)
          DMATX(6,4)=DMATT(2,3,3,3)
          DMATX(5,1)=DMATT(3,1,1,1)
          DMATX(5,3)=DMATT(3,1,1,2)
          DMATX(5,5)=DMATT(3,1,1,3)
          DMATX(5,3)=DMATT(3,1,2,1)
          DMATX(5,2)=DMATT(3,1,2,2)
          DMATX(5,6)=DMATT(3,1,2,3)
          DMATX(5,5)=DMATT(3,1,3,1)
          DMATX(5,6)=DMATT(3,1,3,2)
          DMATX(5,4)=DMATT(3,1,3,3)
          DMATX(6,1)=DMATT(3,2,1,1)
          DMATX(6,3)=DMATT(3,2,1,2)
          DMATX(6,5)=DMATT(3,2,1,3)
          DMATX(6,3)=DMATT(3,2,2,1)
          DMATX(6,2)=DMATT(3,2,2,2)
          DMATX(6,6)=DMATT(3,2,2,3)
          DMATX(6,5)=DMATT(3,2,3,1)
          DMATX(6,6)=DMATT(3,2,3,2)
          DMATX(6,4)=DMATT(3,2,3,3)
          DMATX(4,5)=DMATT(3,3,1,3)
          DMATX(4,6)=DMATT(3,3,2,3)
          DMATX(4,5)=DMATT(3,3,3,1)
          DMATX(4,6)=DMATT(3,3,3,2)
         ENDIF                     ! ntype.ne.5
        ENDIF                      ! ntype.ne.1
       ENDIF                       ! ntype.eq.4
C
      ELSE                         ! istan=2 (non-standard models)
C
C**** SIMO'S MODEL
C
       IF(IFREN.EQ.7) THEN
        CALL DESIMO(XJACM,XJA3M,DETJM,STRAP,
     .              RCGPI,FINET,TRABE)
C
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress - Simo model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=
     .      BULKM*(RCGTI(ISTRS)*RCGTI(JSTRS)-
     .             DLOG(DETJM)*DCINV(ISTRS,JSTRS))+
     .      DISTM*2.0D0/3.0D0*FACTJ*(TRABE*(
     .                           1.0D0/2.0D0*DCINV(ISTRS,JSTRS)+
     .                           1.0D0/3.0D0*RCGTI(ISTRS)*RCGTI(JSTRS))-
     .                           RCGPI(ISTRS)*RCGTI(JSTRS)-
     .                           RCGTI(ISTRS)*UNOMA(JSTRS))
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
       ENDIF                       ! ifren.eq.7
C
C**** HYPERELASTIC LINEAR SPATIAL MODEL (Hooke's law)
C
       IF(IFREN.EQ.8) THEN         ! 1D
        D1=BULKM+(4.0D0/3.0D0)*DISTM
        D2=BULKM-(2.0D0/3.0D0)*DISTM
        D3=BULKM+(1.0D0/3.0D0)*DISTM
        DMATX(1,1)=D1*RCGTI(1)*RCGTI(1)
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR IN FOURTH')
         DMATX(1,1)=0.0D0          ! to be implemented
         DMATX(1,2)=0.0D0
         DMATX(1,3)=0.0D0
         DMATX(2,2)=0.0D0
         DMATX(2,3)=0.0D0
C
         POISC1=0.0D0
         POISC2=0.0D0
         POISC3=0.0D0
        ELSE                       ! plane strain, axis. & 3D
         DMATX(1,2)=D2*RCGTI(1)*RCGTI(2)+2.0D0*DISTM*RCGTI(3)*RCGTI(3)
         DMATX(1,3)=D1*RCGTI(1)*RCGTI(3)
         DMATX(1,4)=D2*RCGTI(1)*RCGTI(4)
         DMATX(2,2)=D1*RCGTI(2)*RCGTI(2)
         DMATX(2,3)=D1*RCGTI(2)*RCGTI(3)
         DMATX(2,4)=D2*RCGTI(2)*RCGTI(4)
         DMATX(3,3)=DISTM*RCGTI(1)*RCGTI(2)+D3*RCGTI(3)*RCGTI(3)
         DMATX(3,4)=D2*RCGTI(3)*RCGTI(4)
         DMATX(4,4)=D1*RCGTI(4)*RCGTI(4)
         IF(NTYPE.EQ.4) THEN       ! 3D
          DMATX(1,4)=DMATX(1,4)+2.0D0*DISTM*RCGTI(5)*RCGTI(5)
          DMATX(1,5)=D1*RCGTI(1)*RCGTI(5)
          DMATX(1,6)=D2*RCGTI(1)*RCGTI(6)+2.0D0*DISTM*RCGTI(3)*RCGTI(5)
          DMATX(2,4)=DMATX(2,4)+2.0D0*DISTM*RCGTI(6)*RCGTI(6)
          DMATX(2,5)=D2*RCGTI(2)*RCGTI(5)+2.0D0*DISTM*RCGTI(3)*RCGTI(6)
          DMATX(2,6)=D1*RCGTI(2)*RCGTI(6)
          DMATX(3,4)=DMATX(3,4)+2.0D0*DISTM*RCGTI(5)*RCGTI(6)
          DMATX(3,5)=DISTM*RCGTI(1)*RCGTI(6)+D3*RCGTI(3)*RCGTI(5)
          DMATX(3,6)=DISTM*RCGTI(2)*RCGTI(5)+D3*RCGTI(3)*RCGTI(6)
          DMATX(4,5)=D1*RCGTI(4)*RCGTI(5)
          DMATX(4,6)=D1*RCGTI(4)*RCGTI(6)
          DMATX(5,5)=DISTM*RCGTI(1)*RCGTI(4)+D3*RCGTI(5)*RCGTI(5)
          DMATX(5,6)=DISTM*RCGTI(3)*RCGTI(4)+D3*RCGTI(5)*RCGTI(6)
          DMATX(6,6)=DISTM*RCGTI(2)*RCGTI(4)+D3*RCGTI(6)*RCGTI(6)
         ENDIF                     ! ntype.eq.4
        ENDIF                      ! ntype.eq.1
       ENDIF                       ! ifren.eq.8
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.
     .    IFREN.EQ.58) THEN        ! M-R, Y, O, D, H, M, G & K
        IF(NMODI.EQ.1) FACTJ=1.0D0 ! dev-vol decompos. is not considered
       ENDIF                       ! ifren.eq.51....
C
C**** MOONEY-RIVLIN MODEL
C
       IF(IFREN.EQ.51) THEN
        D1=4.0D0*CPOL2
        D2=4.0D0*CPOL3
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress - MR-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=
     .      D1*(UNOMA(ISTRS)*UNOMA(JSTRS)-UNORA(ISTRS,JSTRS))+
     .      D2*(2.0D0*FACTJ*TRACC*UNOMA(ISTRS)*UNOMA(JSTRS)-
     .  UNOMA(ISTRS)*FACTJ*RCGTT(JSTRS)-FACTJ*RCGTT(ISTRS)*UNOMA(JSTRS)+
     .      (FACTJ*TRACC-3.0D0)*(UNOMA(ISTRS)*UNOMA(JSTRS)-
     .                                              UNORA(ISTRS,JSTRS)))
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=2.0D0*CPOL1*UNOMA(ISTRS)+
     .              2.0D0*CPOL2*FACTJ*(TRACC*UNOMA(ISTRS)-RCGTT(ISTRS))+
     .              2.0D0*CPOL3*((FACTJ*FACTJ*SEINC-3.0D0)*UNOMA(ISTRS)+
     .      (FACTJ*TRACC-3.0D0)*FACTJ*(TRACC*UNOMA(ISTRS)-RCGTT(ISTRS)))
        ENDDO
C
        PSI0=CPOL1*(FACTJ*TRACC-3.D0)+CPOL2*(FACTJ*FACTJ*SEINC-3.D0)+
     .       CPOL3*(FACTJ*TRACC-3.D0)*(FACTJ*FACTJ*SEINC-3.D0)
       ENDIF                       ! ifren.eq.51
C
C**** YEOH MODEL
C
       IF(IFREN.EQ.52) THEN
        D1=8.0D0*CPOL2+24.0D0*CPOL3*(FACTJ*TRACC-3.0D0)
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress-Yeoh-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=D1*UNOMA(ISTRS)*UNOMA(JSTRS)
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=(2.0D0*CPOL1+4.0D0*CPOL2*(FACTJ*TRACC-3.0D0)+
     . 6.0D0*CPOL3*(FACTJ*TRACC-3.0D0)*(FACTJ*TRACC-3.0D0))*UNOMA(ISTRS)
        ENDDO
C
        PSI0=CPOL1*(FACTJ*TRACC-3.D0)+CPOL2*(FACTJ*TRACC-3.D0)**2+
     .       CPOL3*(FACTJ*TRACC-3.D0)**3
       ENDIF                       ! ifren.eq.52
C
C**** OGDEN'S MODEL
C
       IF(IFREN.EQ.53) THEN
        CALL PRIRCG(NTYPE,RCGTT,ALAMB,ILAEQ)
        FACT3=FACTJ*FACTJ*FACTJ
C
        IF(ILAEQ.EQ.1) THEN        ! three different principal stretches
         AL1A1=1.0D0
         AL2A1=1.0D0
         AL3A1=1.0D0
         IF(EPOL1.NE.0.0D0) THEN
          AL1A1=(DSQRT(FACTJ)*ALAMB(1))**EPOL1
          AL2A1=(DSQRT(FACTJ)*ALAMB(2))**EPOL1
          AL3A1=(DSQRT(FACTJ)*ALAMB(3))**EPOL1
         ENDIF
         AL1A2=1.0D0
         AL2A2=1.0D0
         AL3A2=1.0D0
         IF(EPOL2.NE.0.0D0) THEN
          AL1A2=(DSQRT(FACTJ)*ALAMB(1))**EPOL2
          AL2A2=(DSQRT(FACTJ)*ALAMB(2))**EPOL2
          AL3A2=(DSQRT(FACTJ)*ALAMB(3))**EPOL2
         ENDIF
         AL1A3=1.0D0
         AL2A3=1.0D0
         AL3A3=1.0D0
         IF(EPOL3.NE.0.0D0) THEN
          AL1A3=(DSQRT(FACTJ)*ALAMB(1))**EPOL3
          AL2A3=(DSQRT(FACTJ)*ALAMB(2))**EPOL3
          AL3A3=(DSQRT(FACTJ)*ALAMB(3))**EPOL3
         ENDIF
C
         BETA1=CPOL1*AL1A1+CPOL2*AL1A2+CPOL3*AL1A3
         BETA2=CPOL1*AL2A1+CPOL2*AL2A2+CPOL3*AL2A3
         BETA3=CPOL1*AL3A1+CPOL2*AL3A2+CPOL3*AL3A3
C
         AL12=ALAMB(1)*ALAMB(1)
         AL22=ALAMB(2)*ALAMB(2)
         AL32=ALAMB(3)*ALAMB(3)
C
         DETA1=2.0D0*FACTJ*AL12*FACTJ*AL12-FACTJ*TRACC*FACTJ*AL12+
     .                                          FACT3*DETJC/(FACTJ*AL12)
         DETA2=2.0D0*FACTJ*AL22*FACTJ*AL22-FACTJ*TRACC*FACTJ*AL22+
     .                                          FACT3*DETJC/(FACTJ*AL22)
         DETA3=2.0D0*FACTJ*AL32*FACTJ*AL32-FACTJ*TRACC*FACTJ*AL32+
     .                                          FACT3*DETJC/(FACTJ*AL32)
C
         DO ISTRS=1,NSTRS
          AMAA1(ISTRS)=(FACTJ*RCGTT(ISTRS)-
     .                                  UNOMA(ISTRS)*FACTJ*(TRACC-AL12)+
     .                FACT3*DETJC*RCGTI(ISTRS)/FACTJ/(FACTJ*AL12))/DETA1
          AMAA2(ISTRS)=(FACTJ*RCGTT(ISTRS)-
     .                                  UNOMA(ISTRS)*FACTJ*(TRACC-AL22)+
     .                FACT3*DETJC*RCGTI(ISTRS)/FACTJ/(FACTJ*AL22))/DETA2
          AMAA3(ISTRS)=(FACTJ*RCGTT(ISTRS)-
     .                                  UNOMA(ISTRS)*FACTJ*(TRACC-AL32)+
     .                FACT3*DETJC*RCGTI(ISTRS)/FACTJ/(FACTJ*AL32))/DETA3
         ENDDO
C
         BETE1=EPOL1*CPOL1*AL1A1+EPOL2*CPOL2*AL1A2+EPOL3*CPOL3*AL1A3
         BETE2=EPOL1*CPOL1*AL2A1+EPOL2*CPOL2*AL2A2+EPOL3*CPOL3*AL2A3
         BETE3=EPOL1*CPOL1*AL3A1+EPOL2*CPOL2*AL3A2+EPOL3*CPOL3*AL3A3
C
         DO ISTRS=1,NSTRS
          ADET1(ISTRS)=(4.0D0*FACTJ*AL12*FACTJ*AL12-
     .    FACTJ*TRACC*FACTJ*AL12-FACT3*DETJC/(FACTJ*AL12))*AMAA1(ISTRS)-
     .                FACTJ*AL12*UNOMA(ISTRS)+
     .                FACT3*DETJC/(FACTJ*AL12)*RCGTI(ISTRS)/FACTJ
          ADET2(ISTRS)=(4.0D0*FACTJ*AL22*FACTJ*AL22-
     .    FACTJ*TRACC*FACTJ*AL22-FACT3*DETJC/(FACTJ*AL22))*AMAA2(ISTRS)-
     .                FACTJ*AL22*UNOMA(ISTRS)+
     .                FACT3*DETJC/(FACTJ*AL22)*RCGTI(ISTRS)/FACTJ
          ADET3(ISTRS)=(4.0D0*FACTJ*AL32*FACTJ*AL32-
     .    FACTJ*TRACC*FACTJ*AL32-FACT3*DETJC/(FACTJ*AL32))*AMAA3(ISTRS)-
     .                FACTJ*AL32*UNOMA(ISTRS)+
     .                FACT3*DETJC/(FACTJ*AL32)*RCGTI(ISTRS)/FACTJ
         ENDDO
C
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           ATHE1(ISTRS,JSTRS)=UNORA(ISTRS,JSTRS)-
     .                                        UNOMA(ISTRS)*UNOMA(JSTRS)+
     .                             FACTJ*AL12*UNOMA(ISTRS)*AMAA1(JSTRS)-
     .         FACT3*DETJC/(FACTJ*AL12)*RCGTI(ISTRS)/FACTJ*AMAA1(JSTRS)+
     .                                     1.0D0/AL12*DCOFA(ISTRS,JSTRS)
           ATHE2(ISTRS,JSTRS)=UNORA(ISTRS,JSTRS)-
     .                                        UNOMA(ISTRS)*UNOMA(JSTRS)+
     .                             FACTJ*AL22*UNOMA(ISTRS)*AMAA2(JSTRS)-
     .         FACT3*DETJC/(FACTJ*AL22)*RCGTI(ISTRS)/FACTJ*AMAA2(JSTRS)+
     .                                     1.0D0/AL22*DCOFA(ISTRS,JSTRS)
           ATHE3(ISTRS,JSTRS)=UNORA(ISTRS,JSTRS)-
     .                                        UNOMA(ISTRS)*UNOMA(JSTRS)+
     .                             FACTJ*AL32*UNOMA(ISTRS)*AMAA3(JSTRS)-
     .         FACT3*DETJC/(FACTJ*AL32)*RCGTI(ISTRS)/FACTJ*AMAA3(JSTRS)+
     .                                     1.0D0/AL32*DCOFA(ISTRS,JSTRS)
          ENDDO
         ENDDO
C
         IF(NTYPE.EQ.1) THEN       ! plane stress
          CALL RUNEND('ERROR: plane stress-Ogden model not implemented')
         ELSE                      ! 1D, plane strain, axis. & 3D
          DO ISTRS=1,NSTRS
           DO JSTRS=ISTRS,NSTRS
            DMATX(ISTRS,JSTRS)=AMAA1(ISTRS)*BETE1*AMAA1(JSTRS)+
     . BETA1*2.0D0/DETA1*(ATHE1(ISTRS,JSTRS)-AMAA1(ISTRS)*ADET1(JSTRS))+
     .                         AMAA2(ISTRS)*BETE2*AMAA2(JSTRS)+
     . BETA2*2.0D0/DETA2*(ATHE2(ISTRS,JSTRS)-AMAA2(ISTRS)*ADET2(JSTRS))+
     .                         AMAA3(ISTRS)*BETE3*AMAA3(JSTRS)+
     . BETA3*2.0D0/DETA3*(ATHE3(ISTRS,JSTRS)-AMAA3(ISTRS)*ADET3(JSTRS))
           ENDDO
          ENDDO
         ENDIF                     ! ntype.eq.1
C
         DO ISTRS=1,NSTRS          ! standard deviatoric term (S-tilde)
          STILD(ISTRS)=BETA1*AMAA1(ISTRS)+BETA2*AMAA2(ISTRS)+
     .                 BETA3*AMAA3(ISTRS)
         ENDDO
        ENDIF                      ! ilaeq.eq.1
C
        IF(ILAEQ.EQ.2) THEN        ! two equal principal stretches
         AL1A1=1.0D0
         AL2A1=1.0D0
         IF(EPOL1.NE.0.0D0) THEN
          AL1A1=(DSQRT(FACTJ)*ALAMB(1))**EPOL1
          AL2A1=(DSQRT(FACTJ)*ALAMB(2))**EPOL1
         ENDIF
         AL1A2=1.0D0
         AL2A2=1.0D0
         IF(EPOL2.NE.0.0D0) THEN
          AL1A2=(DSQRT(FACTJ)*ALAMB(1))**EPOL2
          AL2A2=(DSQRT(FACTJ)*ALAMB(2))**EPOL2
         ENDIF
         AL1A3=1.0D0
         AL2A3=1.0D0
         IF(EPOL3.NE.0.0D0) THEN
          AL1A3=(DSQRT(FACTJ)*ALAMB(1))**EPOL3
          AL2A3=(DSQRT(FACTJ)*ALAMB(2))**EPOL3
         ENDIF
C
         AL12=ALAMB(1)*ALAMB(1)
         AL22=ALAMB(2)*ALAMB(2)
         ALDIF=AL12-AL22
C
         AASTE=CPOL1*AL1A1+CPOL2*AL1A2+CPOL3*AL1A3
         BASTE=CPOL1*AL2A1+CPOL2*AL2A2+CPOL3*AL2A3
C
         DPAL1=-2.0D0*DSQRT(FACTJ)*ALAMB(1)/(FACTJ*ALDIF*FACTJ*ALDIF)
         DPAL2= 2.0D0*DSQRT(FACTJ)*ALAMB(2)/(FACTJ*ALDIF*FACTJ*ALDIF)
C
         AASTD=(CPOL1*EPOL1*AL1A1+CPOL2*EPOL2*AL1A2+CPOL3*EPOL3*AL1A3)/
     .                                           (DSQRT(FACTJ)*ALAMB(1))
         BASTD=(CPOL1*EPOL1*AL2A1+CPOL2*EPOL2*AL2A2+CPOL3*EPOL3*AL2A3)/
     .                                           (DSQRT(FACTJ)*ALAMB(2))
C
         DO ISTRS=1,NSTRS
          AMAA1(ISTRS)= DSQRT(FACTJ)*ALAMB(1)/(2.0D0*FACTJ*ALDIF)*
     .                                  (UNOMA(ISTRS)-AL22*RCGTI(ISTRS))
          AMAA2(ISTRS)=-DSQRT(FACTJ)*ALAMB(2)/(4.0D0*FACTJ*ALDIF)*
     .                                  (UNOMA(ISTRS)-AL12*RCGTI(ISTRS))
         ENDDO
C
         IF(NTYPE.EQ.1) THEN       ! plane stress
          CALL RUNEND('ERROR: plane stress-Ogden model not implemented')
         ELSE                      ! 1D, plane strain, axis. & 3D
          DO ISTRS=1,NSTRS
           DO JSTRS=ISTRS,NSTRS
            DMATX(ISTRS,JSTRS)=2.0D0*(
     .                          AASTE*(UNOMA(ISTRS)-AL22*RCGTI(ISTRS))-
     .                          BASTE*(UNOMA(ISTRS)-AL12*RCGTI(ISTRS)))*
     .                          (DPAL1*AMAA1(JSTRS)+DPAL2*AMAA2(JSTRS))+
     .                         2.0D0/(FACTJ*ALDIF)*(
     .              (UNOMA(ISTRS)-AL22*RCGTI(ISTRS))*AASTD*AMAA1(JSTRS)-
     .      2.0D0*AASTE*ALAMB(2)*RCGTI(ISTRS)/DSQRT(FACTJ)*AMAA2(JSTRS)-
     .         AASTE*FACTJ*AL22*(FACTJ*DCOFA(ISTRS,JSTRS)/(FACT3*DETJC)-
     .                           RCGTI(ISTRS)/FACTJ*RCGTI(JSTRS)/FACTJ)-
     .              (UNOMA(ISTRS)-AL12*RCGTI(ISTRS))*BASTD*AMAA2(JSTRS)+
     .      2.0D0*BASTE*ALAMB(1)*RCGTI(ISTRS)/DSQRT(FACTJ)*AMAA1(JSTRS)+
     .         BASTE*FACTJ*AL12*(FACTJ*DCOFA(ISTRS,JSTRS)/(FACT3*DETJC)-
     .                           RCGTI(ISTRS)/FACTJ*RCGTI(JSTRS)/FACTJ))
           ENDDO
          ENDDO
         ENDIF                     ! ntype.eq.1
C
         DO ISTRS=1,NSTRS          ! standard deviatoric term (S-tilde)
          STILD(ISTRS)=1.0D0/(FACTJ*ALDIF)*(
     .                           AASTE*(UNOMA(ISTRS)-AL22*RCGTI(ISTRS))-
     .                           BASTE*(UNOMA(ISTRS)-AL12*RCGTI(ISTRS)))
         ENDDO
        ENDIF                      ! ilaeq.eq.2
C
        IF(ILAEQ.EQ.3) THEN        ! three equal principal stretches
         AL1A1=1.0D0
         IF(EPOL1.NE.0.0D0) THEN
          AL1A1=(DSQRT(FACTJ)*ALAMB(1))**EPOL1
         ENDIF
         AL1A2=1.0D0
         IF(EPOL2.NE.0.0D0) THEN
          AL1A2=(DSQRT(FACTJ)*ALAMB(1))**EPOL2
         ENDIF
         AL1A3=1.0D0
         IF(EPOL3.NE.0.0D0) THEN
          AL1A3=(DSQRT(FACTJ)*ALAMB(1))**EPOL3
         ENDIF
C
         AL12=ALAMB(1)*ALAMB(1)
C
         ALX=(CPOL1*AL1A1*(EPOL1-2.0D0)+
     .        CPOL2*AL1A2*(EPOL2-2.0D0)+
     .        CPOL3*AL1A3*(EPOL3-2.0D0))/(3.0D0*FACTJ*AL12*FACTJ*AL12)
C
         IF(NTYPE.EQ.1) THEN       ! plane stress
          CALL RUNEND('ERROR: plane stress-Ogden model not implemented')
         ELSE                      ! 1D, plane strain, axis. & 3D
          DO ISTRS=1,NSTRS
           DO JSTRS=ISTRS,NSTRS
            DMATX(ISTRS,JSTRS)=ALX*UNOMA(ISTRS)*UNOMA(JSTRS)
           ENDDO
          ENDDO
         ENDIF                     ! ntype.eq.1
C
         DO ISTRS=1,NSTRS          ! standard deviatoric term (S-tilde)
          STILD(ISTRS)=(CPOL1*AL1A1+CPOL2*AL1A2+CPOL3*AL1A3)/
     .                                         (FACTJ*AL12)*UNOMA(ISTRS)
         ENDDO
        ENDIF                      ! ilaeq.eq.3
       ENDIF                       ! ifren.eq.53
C
C**** DELFINO MODEL
C
       IF(IFREN.EQ.54) THEN
        D1=DEXP(CPOL2/2.0D0*(FACTJ*TRACC-3.0D0))
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress-Delf-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=CPOL1*CPOL2*D1*UNOMA(ISTRS)*UNOMA(JSTRS)
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=CPOL1*D1*UNOMA(ISTRS)
        ENDDO
       ENDIF                       ! ifren.eq.54
C
C**** HOLZAPFEL MODEL
C
       IF(IFREN.EQ.55) THEN
        ANIS1(1)=VANIS(1,1)*VANIS(1,1)       ! a0 o a0
        ANIS1(2)=VANIS(1,2)*VANIS(1,2)
        ANIS1(3)=VANIS(1,1)*VANIS(1,2)
        ANIS1(4)=VANIS(1,3)*VANIS(1,3)
        ANIS1(5)=VANIS(1,1)*VANIS(1,3)
        ANIS1(6)=VANIS(1,2)*VANIS(1,3)
        ANIS2(1)=VANIS(2,1)*VANIS(2,1)       ! b0 o b0
        ANIS2(2)=VANIS(2,2)*VANIS(2,2)
        ANIS2(3)=VANIS(2,1)*VANIS(2,2)
        ANIS2(4)=VANIS(2,3)*VANIS(2,3)
        ANIS2(5)=VANIS(2,1)*VANIS(2,3)
        ANIS2(6)=VANIS(2,2)*VANIS(2,3)
        AINV4=RCGTT(1)*ANIS1(1)+RCGTT(2)*ANIS1(2)+RCGTT(4)*ANIS1(4)+
     .                    2.0D0*RCGTT(3)*ANIS1(3)
        AINV6=RCGTT(1)*ANIS2(1)+RCGTT(2)*ANIS2(2)+RCGTT(4)*ANIS2(4)+
     .                    2.0D0*RCGTT(3)*ANIS2(3)
        IF(NTYPE.EQ.4) THEN        ! 3D
         AINV4=AINV4+     2.0D0*RCGTT(5)*ANIS1(5)+
     .                    2.0D0*RCGTT(6)*ANIS1(6)
         AINV6=AINV6+     2.0D0*RCGTT(5)*ANIS2(5)+
     .                    2.0D0*RCGTT(6)*ANIS2(6)
        ENDIF
        AINV4=FACTJ*AINV4
        AINV6=FACTJ*AINV6
C
        D42=(AINV4-1.0D0)*(AINV4-1.0D0)
        D62=(AINV6-1.0D0)*(AINV6-1.0D0)
        D4=DEXP(CPOL3*D42)
        D6=DEXP(CPOL3*D62)
        IF(IFREM.EQ.1) THEN        ! fibers resist only tension
         IF(AINV4.LT.1.0D0) D4=0.0D0
         IF(AINV6.LT.1.0D0) D6=0.0D0
        ENDIF
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress-Holz-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=4.0D0*CPOL2*(
     .             D4*(2.0D0*CPOL3*D42+1.0D0)*ANIS1(ISTRS)*ANIS1(JSTRS)+
     .             D6*(2.0D0*CPOL3*D62+1.0D0)*ANIS2(ISTRS)*ANIS2(JSTRS))
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=CPOL1*UNOMA(ISTRS)+
     .          2.0D0*CPOL2*D4*(AINV4-1.0D0)*ANIS1(ISTRS)+
     .          2.0D0*CPOL2*D6*(AINV6-1.0D0)*ANIS2(ISTRS)
        ENDDO
       ENDIF                       ! ifren.eq.55
C
C**** MYOCARDIUM MODEL
C
       IF(IFREN.EQ.56) THEN
        ANIS1(1)=VANIS(1,1)*VANIS(1,1)       ! f0 o f0
        ANIS1(2)=VANIS(1,2)*VANIS(1,2)
        ANIS1(3)=VANIS(1,1)*VANIS(1,2)
        ANIS1(4)=VANIS(1,3)*VANIS(1,3)
        ANIS1(5)=VANIS(1,1)*VANIS(1,3)
        ANIS1(6)=VANIS(1,2)*VANIS(1,3)
        ANIS2(1)=VANIS(2,1)*VANIS(2,1)       ! s0 o s0
        ANIS2(2)=VANIS(2,2)*VANIS(2,2)
        ANIS2(3)=VANIS(2,1)*VANIS(2,2)
        ANIS2(4)=VANIS(2,3)*VANIS(2,3)
        ANIS2(5)=VANIS(2,1)*VANIS(2,3)
        ANIS2(6)=VANIS(2,2)*VANIS(2,3)
c       ANIS3(1)=VANIS(1,1)*VANIS(2,1)      ! f0 o s0 (Holzapfel; wrong)
c       ANIS3(2)=VANIS(1,2)*VANIS(2,2)
c       ANIS3(3)=VANIS(1,1)*VANIS(2,2)
c       ANIS3(4)=VANIS(1,3)*VANIS(2,3)
c       ANIS3(5)=VANIS(1,1)*VANIS(2,3)
c       ANIS3(6)=VANIS(1,2)*VANIS(2,3)
        ANIS3(1)=VANIS(1,1)*VANIS(2,1)      ! sym (f0 o s0) (Goktepe)
        ANIS3(2)=VANIS(1,2)*VANIS(2,2)
        ANIS3(3)=(VANIS(1,1)*VANIS(2,2)+VANIS(1,2)*VANIS(2,1))*0.5D0
        ANIS3(4)=VANIS(1,3)*VANIS(2,3)
        ANIS3(5)=(VANIS(1,1)*VANIS(2,3)+VANIS(1,3)*VANIS(2,1))*0.5D0
        ANIS3(6)=(VANIS(1,2)*VANIS(2,3)+VANIS(1,3)*VANIS(2,2))*0.5D0
        AIN4F=RCGTT(1)*ANIS1(1)+RCGTT(2)*ANIS1(2)+RCGTT(4)*ANIS1(4)+
     .                    2.0D0*RCGTT(3)*ANIS1(3)+
     .                    2.0D0*RCGTT(5)*ANIS1(5)+
     .                    2.0D0*RCGTT(6)*ANIS1(6)
        AIN4S=RCGTT(1)*ANIS2(1)+RCGTT(2)*ANIS2(2)+RCGTT(4)*ANIS2(4)+
     .                    2.0D0*RCGTT(3)*ANIS2(3)+
     .                    2.0D0*RCGTT(5)*ANIS2(5)+
     .                    2.0D0*RCGTT(6)*ANIS2(6)
        AI8FS=RCGTT(1)*ANIS3(1)+RCGTT(2)*ANIS3(2)+RCGTT(4)*ANIS3(4)+
     .                    2.0D0*RCGTT(3)*ANIS3(3)+
     .                    2.0D0*RCGTT(5)*ANIS3(5)+
     .                    2.0D0*RCGTT(6)*ANIS3(6)
        AIN4F=FACTJ*AIN4F
        AIN4S=FACTJ*AIN4S
        AI8FS=FACTJ*AI8FS
C
        D4F=(AIN4F-1.0D0)*(AIN4F-1.0D0)
        D4S=(AIN4S-1.0D0)*(AIN4S-1.0D0)
        D8FS=AI8FS*AI8FS
        D1 =DEXP(CPOL2*(TRACC-3.0D0))
        DF =DEXP(CPOL4*D4F)
        DS =DEXP(CPOL6*D4S)
        DFS=DEXP(CPOL8*D8FS)
        IF(IFREM.EQ.1) THEN        ! fibers resist only tension
         IF(AIN4F.LT.1.0D0) DF=0.0D0
         IF(AIN4S.LT.1.0D0) DS=0.0D0
         IF(AI8FS.LT.1.0D0) DFS=0.0D0
        ENDIF
        DO ISTRS=1,NSTRS
         DO JSTRS=ISTRS,NSTRS
          DMATX(ISTRS,JSTRS)=
     . 2.0D0*CPOL1*CPOL2*D1*UNOMA(ISTRS)*UNOMA(JSTRS)+
     . 4.0D0*CPOL3*DF*(2.0D0*CPOL4*D4F+1.0D0)*ANIS1(ISTRS)*ANIS1(JSTRS)+
     . 4.0D0*CPOL5*DS*(2.0D0*CPOL6*D4S+1.0D0)*ANIS2(ISTRS)*ANIS2(JSTRS)+
     . 4.0D0*CPOL7*DFS*(2.0D0*CPOL8*D8FS+1.0D0)*
     .                                        ANIS3(ISTRS)*ANIS3(JSTRS)
         ENDDO
        ENDDO
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=CPOL1*D1*UNOMA(ISTRS)+
     .          2.0D0*CPOL3*DF*(AIN4F-1.0D0)*ANIS1(ISTRS)+
     .          2.0D0*CPOL5*DS*(AIN4S-1.0D0)*ANIS2(ISTRS)+
     .          2.0D0*CPOL7*DFS*AI8FS*ANIS3(ISTRS)
        ENDDO
       ENDIF                       ! ifren.eq.56
C
C**** GASSER MODEL
C
       IF(IFREN.EQ.57) THEN
        ANIS1(1)=VANIS(1,1)*VANIS(1,1)*(1.0D0-3.0D0*CPOL4)+    ! a0 o a0
     .           UNOMA(1)*CPOL4
        ANIS1(2)=VANIS(1,2)*VANIS(1,2)*(1.0D0-3.0D0*CPOL4)+
     .           UNOMA(2)*CPOL4
        ANIS1(3)=VANIS(1,1)*VANIS(1,2)*(1.0D0-3.0D0*CPOL4)
        ANIS1(4)=VANIS(1,3)*VANIS(1,3)*(1.0D0-3.0D0*CPOL4)+
     .           UNOMA(4)*CPOL4
        ANIS1(5)=VANIS(1,1)*VANIS(1,3)*(1.0D0-3.0D0*CPOL4)
        ANIS1(6)=VANIS(1,2)*VANIS(1,3)*(1.0D0-3.0D0*CPOL4)
        ANIS2(1)=VANIS(2,1)*VANIS(2,1)*(1.0D0-3.0D0*CPOL4)+    ! b0 o b0
     .           UNOMA(1)*CPOL4
        ANIS2(2)=VANIS(2,2)*VANIS(2,2)*(1.0D0-3.0D0*CPOL4)+
     .           UNOMA(2)*CPOL4
        ANIS2(3)=VANIS(2,1)*VANIS(2,2)*(1.0D0-3.0D0*CPOL4)
        ANIS2(4)=VANIS(2,3)*VANIS(2,3)*(1.0D0-3.0D0*CPOL4)+
     .           UNOMA(4)*CPOL4
        ANIS2(5)=VANIS(2,1)*VANIS(2,3)*(1.0D0-3.0D0*CPOL4)
        ANIS2(6)=VANIS(2,2)*VANIS(2,3)*(1.0D0-3.0D0*CPOL4)
        AINV4=(RCGTT(1)*ANIS1(1)+RCGTT(2)*ANIS1(2)+RCGTT(4)*ANIS1(4)+
     .                     2.0D0*RCGTT(3)*ANIS1(3))*(1.0D0-3.0D0*CPOL4)+
     .                          (RCGTT(1)+RCGTT(2)+RCGTT(4))*CPOL4
        AINV6=(RCGTT(1)*ANIS2(1)+RCGTT(2)*ANIS2(2)+RCGTT(4)*ANIS2(4)+
     .                     2.0D0*RCGTT(3)*ANIS2(3))*(1.0D0-3.0D0*CPOL4)+
     .                          (RCGTT(1)+RCGTT(2)+RCGTT(4))*CPOL4
        IF(NTYPE.EQ.4) THEN        ! 3D
         AINV4=AINV4+     (2.0D0*RCGTT(5)*ANIS1(5)+
     .                     2.0D0*RCGTT(6)*ANIS1(6))*(1.0D0-3.0D0*CPOL4)
         AINV6=AINV6+     (2.0D0*RCGTT(5)*ANIS2(5)+
     .                     2.0D0*RCGTT(6)*ANIS2(6))*(1.0D0-3.0D0*CPOL4)
        ENDIF
        AINV4=FACTJ*AINV4
        AINV6=FACTJ*AINV6
C
        D42=(AINV4-1.0D0)*(AINV4-1.0D0)
        D62=(AINV6-1.0D0)*(AINV6-1.0D0)
        D4=DEXP(CPOL3*D42)
        D6=DEXP(CPOL3*D62)
        IF(IFREM.EQ.1) THEN        ! fibers resist only tension
         IF(AINV4.LT.1.0D0) D4=0.0D0
         IF(AINV6.LT.1.0D0) D6=0.0D0
        ENDIF
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress-Holz-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=4.0D0*CPOL2*(
     .             D4*(2.0D0*CPOL3*D42+1.0D0)*ANIS1(ISTRS)*ANIS1(JSTRS)+
     .             D6*(2.0D0*CPOL3*D62+1.0D0)*ANIS2(ISTRS)*ANIS2(JSTRS))
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=CPOL1*UNOMA(ISTRS)+
     .          2.0D0*CPOL2*D4*(AINV4-1.0D0)*ANIS1(ISTRS)+
     .          2.0D0*CPOL2*D6*(AINV6-1.0D0)*ANIS2(ISTRS)
        ENDDO
       ENDIF                       ! ifren.eq.57
C
C**** KAN MODEL
C
       IF(IFREN.EQ.58) THEN
        D1=4.0D0*CPOL1*(CPOL2*DEXP(CPOL2*(FACTJ*TRACC-3.0D0))-
     .            CPOL3/(FACTJ*TRACC-2.0D0))     ! 2nd derivative
        IF(NTYPE.EQ.1) THEN        ! plane stress
         CALL RUNEND('ERROR: plane stress-KAN-tilde model not implem.')
        ELSE                       ! 1D, plane strain, axis. & 3D
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=D1*UNOMA(ISTRS)*UNOMA(JSTRS)
          ENDDO
         ENDDO
        ENDIF                      ! ntype.eq.1
C
        DO ISTRS=1,NSTRS           ! standard deviatoric term (S-tilde)
         STILD(ISTRS)=2.0D0*(CPOL1*(DEXP(CPOL2*(FACTJ*TRACC-3.0D0))-
     .                CPOL3*DLOG(FACTJ*TRACC-2.0D0))*UNOMA(ISTRS))
        ENDDO
C
        PSI0=CPOL1*(1.0D0/CPOL2*DEXP(CPOL2*(FACTJ*TRACC-3.0D0))+
     .       CPOL3*(FACTJ*TRACC-2.0D0)*(1.0D0-DLOG(FACTJ*TRACC-2.0D0))-
     .       1.0D0/CPOL2-CPOL3)
       ENDIF                       ! ifren.eq.58
C
       IF(IFREN.EQ.51.OR.IFREN.EQ.52.OR.IFREN.EQ.53.OR.IFREN.EQ.54.OR.
     .    IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57.OR.
     .    IFREN.EQ.58) THEN         ! M-R, Y, O, D, H, M, G & K
        IDV=0                       ! dev-vol decomposition index
C
C**** DAMAGE MODELS
C
        IF(IDAMG.EQ.1) THEN         ! HFUNC independent of DAMAG
         HFUNC=0.D0                 ! initialization
         DO I=1,NSTRS
          DO J=I,NSTRS
           S0XS0(I,J)=0.D0          ! symmetric tensor
          END DO
         END DO
         DAMAG=EBASE(1)   ! damage
         THRES=EBASE(2)   ! maximum value of EQUIV over the past history
         IF(PSI0.GT.0.D0) EQUIV=DSQRT(2.D0*PSI0)     ! equivalent strain
         IF(EQUIV.GT.THRES) THEN
          IF(EQUIV.GE.CEMIN.AND.EQUIV.LE.CEMAX) THEN
           XN=EQUIV-CEMIN
           XX=CEMAX-CEMIN
           HFUNC=2.D0*XN/(XX*XX)*(1.D0+CEETA-2.D0*CEETA*XN/XX*XN/XX)
          ENDIF
          DAMAG=DAMAG+(EQUIV-THRES)*HFUNC
          DO I=1,NSTRS
           DO J=I,NSTRS
            S0XS0(I,J)=HFUNC*STILD(I)*STILD(J)
           END DO
          END DO
          THRES=EQUIV               ! updates THRES
         END IF
C
         DO I=1,NSTRS
          STILD(I)=(1-DAMAG)*STILD(I)
         END DO
C
         DO I=1,NSTRS
          DO J=I,NSTRS
           DMATX(I,J)=(1-DAMAG)*DMATX(I,J)-S0XS0(I,J)
          END DO
         END DO
        ENDIF
        IF(IDAMG.EQ.2) THEN

        ENDIF
C
C**** VISCOELASTIC MODELS
C
        IF(NCHAI.GT.0) THEN         ! viscoelastic response
         IF(NMODI.NE.1) THEN        ! dev-vol decompos. is considered
          IDV=1
          ATILD=0.0D0
          DO ISTRS=1,NSTRS          ! auxiliar term (A-tilde)
           ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
           IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .      ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
          ENDDO
          DO ISTRS=1,NSTRS          ! auxiliar term (B-tilde)
           BTILD(ISTRS)=0.0D0
           DO JSTRS=1,NSTRS
            I1=ISTRS
            I2=JSTRS
            IF(JSTRS.LT.ISTRS) THEN ! only upper triangle of D is used
             I1=JSTRS
             I2=ISTRS
            ENDIF
            BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
            IF(I2.EQ.3.OR.I2.EQ.5.OR.I2.EQ.6)
     .       BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
           ENDDO
          ENDDO
          AHATT=0.0D0
          DO ISTRS=1,NSTRS          ! auxiliar term (A-hat)
           AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
           IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .      AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
          ENDDO
C
          DO ISTRS=1,NSTRS          ! deviatoric correction terms
           DO JSTRS=ISTRS,NSTRS
            DMATX(ISTRS,JSTRS)=FACTJ*(FACTJ*DMATX(ISTRS,JSTRS)-
     . 1.0D0/3.0D0*RCGTI(ISTRS)*(FACTJ*BTILD(JSTRS)+2.0D0*STILD(JSTRS))-
     . 1.0D0/3.0D0*(FACTJ*BTILD(ISTRS)+2.0D0*STILD(ISTRS))*RCGTI(JSTRS)-
     . 2.0D0/3.0D0*ATILD/DETJC*DCOFA(ISTRS,JSTRS)+
     . 1.0D0/9.0D0*(8.0D0*ATILD+FACTJ*AHATT)*RCGTI(ISTRS)*RCGTI(JSTRS))
           ENDDO
          ENDDO
C
          DO ISTRS=1,NSTRS
           STILD(ISTRS)=            ! dev. correction
     .      FACTJ*(STILD(ISTRS)-1.0D0/3.0D0*ATILD*RCGTI(ISTRS))
          ENDDO
         ENDIF                      ! nmodi.ne.1
C
         DO IN=1,NCHAI
          DO ISTRS=1,NSTRS
           IX=IN*NSTRS+ISTRS
           EBASE(IX)=DEXP(-DTIME/CPOL9(IN))*EBASE(IX)+       ! Q_a
     .               DEXP(-0.5D0*DTIME/CPOL9(IN))*CPOL10(IN)*
     .                                       (STILD(ISTRS)-EBASE(ISTRS))
          ENDDO
         ENDDO
         DO ISTRS=1,NSTRS
          EBASE(ISTRS)=STILD(ISTRS)                          ! S^inf_iso
         ENDDO
         DO ISTRS=1,NSTRS
          DO IN=1,NCHAI
           IX=IN*NSTRS+ISTRS
           STILD(ISTRS)=STILD(ISTRS)+EBASE(IX)               ! S_iso
          ENDDO
         ENDDO
C
         DTABE=0.0D0
         DO IN=1,NCHAI
          DTABE=DTABE+CPOL10(IN)*DEXP(-0.5D0*DTIME/CPOL9(IN))
         ENDDO
         DO ISTRS=1,NSTRS
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=(1.0D0+DTABE)*DMATX(ISTRS,JSTRS)
          ENDDO
         ENDDO
        ENDIF                       ! nchai.gt.0
C
C**** FIBER MODELS
C
        IF(NFIBN.GT.0) THEN         ! fiber-reinforced response
         IF(IDV.EQ.0.AND.
     .      NMODI.NE.1) THEN        ! dev-vol decompos. is considered
          IDV=1
          ATILD=0.0D0
          DO ISTRS=1,NSTRS          ! auxiliar term (A-tilde)
           ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
           IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .      ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
          ENDDO
          DO ISTRS=1,NSTRS          ! auxiliar term (B-tilde)
           BTILD(ISTRS)=0.0D0
           DO JSTRS=1,NSTRS
            I1=ISTRS
            I2=JSTRS
            IF(JSTRS.LT.ISTRS) THEN ! only upper triangle of D is used
             I1=JSTRS
             I2=ISTRS
            ENDIF
            BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
            IF(I2.EQ.3.OR.I2.EQ.5.OR.I2.EQ.6)
     .       BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
           ENDDO
          ENDDO
          AHATT=0.0D0
          DO ISTRS=1,NSTRS          ! auxiliar term (A-hat)
           AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
           IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .      AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
          ENDDO
C
          DO ISTRS=1,NSTRS          ! deviatoric correction terms
           DO JSTRS=ISTRS,NSTRS
            DMATX(ISTRS,JSTRS)=FACTJ*(FACTJ*DMATX(ISTRS,JSTRS)-
     . 1.0D0/3.0D0*RCGTI(ISTRS)*(FACTJ*BTILD(JSTRS)+2.0D0*STILD(JSTRS))-
     . 1.0D0/3.0D0*(FACTJ*BTILD(ISTRS)+2.0D0*STILD(ISTRS))*RCGTI(JSTRS)-
     . 2.0D0/3.0D0*ATILD/DETJC*DCOFA(ISTRS,JSTRS)+
     . 1.0D0/9.0D0*(8.0D0*ATILD+FACTJ*AHATT)*RCGTI(ISTRS)*RCGTI(JSTRS))
           ENDDO
          ENDDO
C
          DO ISTRS=1,NSTRS
           STILD(ISTRS)=            ! dev. correction
     .      FACTJ*(STILD(ISTRS)-1.0D0/3.0D0*ATILD*RCGTI(ISTRS))
          ENDDO
         ENDIF                      ! idv.eq.0.and.nmodi.ne.1
C
         IF(NFIBM.EQ.1) THEN
          DO IN=1,NFIBN
           VANIX(1)=VANIS(IN,1)*VANIS(IN,1)          ! N_theta x N_theta
           VANIX(2)=VANIS(IN,2)*VANIS(IN,2)
           VANIX(3)=VANIS(IN,1)*VANIS(IN,2)
           VANIX(4)=VANIS(IN,3)*VANIS(IN,3)
           VANIX(5)=VANIS(IN,1)*VANIS(IN,3)
           VANIX(6)=VANIS(IN,2)*VANIS(IN,3)
C
           ALAMF=0.0D0                               ! lambda_f
           DO ISTRS=1,NSTRS
            XX=1.0D0
            IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) XX=2.0D0
            ALAMF=ALAMF+XX*VANIX(ISTRS)*RCGTT(ISTRS)
           ENDDO
           ALAMF=DSQRT(ALAMF)
C
           SFIBR=0.0D0                      ! stress (non-linear spring)
           DFIBR=0.0D0                      ! stress derivative
           IF(ALAMF.GT.1.0D0) THEN
            SFIBR=CFIB1*(ALAMF-1.0D0)**CFIB2
            DFIBR=CFIB1*CFIB2*(ALAMF-1.0D0)**(CFIB2-1.0D0)
           ENDIF
           IF((ALAMF-1.0D0).GT.0.0D0) THEN  ! tension
            SFIBS=CFIB3                     ! stress (skidding block)
            DFIBS=0.0D0                     ! stress derivative
            IF((ALAMF-1.0D0).LT.(CFIB3/CFIB4)) THEN
             SFIBS=CFIB4*(ALAMF-1.0D0)
             DFIBS=CFIB4
            ENDIF
           ENDIF
           IF((ALAMF-1.0D0).LE.0.0D0) THEN  ! compression
            SFIBS=-CFIB3                    ! stress (skidding block)
            DFIBS=0.0D0                     ! stress derivative
            IF(DABS(ALAMF-1.0D0).LT.(CFIB3/CFIB4)) THEN
             SFIBS=CFIB4*(ALAMF-1.0D0)
             DFIBS=CFIB4
            ENDIF
           ENDIF
           SFIBF=SFIBR+SFIBS
           DFIBF=DFIBR+DFIBS
C
           DO ISTRS=1,NSTRS
            STILD(ISTRS)=STILD(ISTRS)+               ! S_(matrix+fibers)
     .       VANIS(IN,4)*SFIBF/ALAMF*VANIX(ISTRS)
           ENDDO
C
           DO ISTRS=1,NSTRS
            DO JSTRS=ISTRS,NSTRS
             DMATX(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)+  ! C_(matrix+fibers)
     .        2.0D0*VANIS(IN,4)*(DFIBF-SFIBF/ALAMF)/(ALAMF*ALAMF)*
     .        VANIX(ISTRS)*VANIX(JSTRS)
            ENDDO
           ENDDO
          ENDDO
         ENDIF                      ! nfibm.eq.1
C
         IF(NFIBM.EQ.2) THEN
          DO IN=1,NFIBN
           VANIX(1)=VANIS(IN,1)*VANIS(IN,1)          ! N_theta x N_theta
           VANIX(2)=VANIS(IN,2)*VANIS(IN,2)
           VANIX(3)=VANIS(IN,1)*VANIS(IN,2)
           VANIX(4)=VANIS(IN,3)*VANIS(IN,3)
           VANIX(5)=VANIS(IN,1)*VANIS(IN,3)
           VANIX(6)=VANIS(IN,2)*VANIS(IN,3)
C
           ALAMF=0.0D0                               ! lambda_f
           DO ISTRS=1,NSTRS
            XX=1.0D0
            IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6) XX=2.0D0
            ALAMF=ALAMF+XX*VANIX(ISTRS)*RCGTT(ISTRS)
           ENDDO
           ALAMF=DSQRT(ALAMF)
C
           SFIBR=0.0D0                      ! stress (non-linear spring)
           DFIBR=0.0D0                      ! stress derivative
           IF(ALAMF.GT.1.0D0) THEN
            SFIBR=CFIB1*(ALAMF-1.0D0)+
     .            CFIB2*(ALAMF-1.0D0)*(ALAMF-1.0D0)+
     .            CFIB3*(ALAMF-1.0D0)*(ALAMF-1.0D0)*(ALAMF-1.0D0)+
     .            CFIB4*(ALAMF-1.0D0)*(ALAMF-1.0D0)*(ALAMF-1.0D0)*
     .                  (ALAMF-1.0D0)+
     .            CFIB5*(ALAMF-1.0D0)*(ALAMF-1.0D0)*(ALAMF-1.0D0)*
     .                  (ALAMF-1.0D0)*(ALAMF-1.0D0)
            DFIBR=CFIB1+2.0D0*CFIB2*(ALAMF-1.0D0)+
     .            3.0D0*CFIB3*(ALAMF-1.0D0)*(ALAMF-1.0D0)+
     .            4.0D0*CFIB4*(ALAMF-1.0D0)*(ALAMF-1.0D0)*(ALAMF-1.0D0)+
     .            5.0D0*CFIB5*(ALAMF-1.0D0)*(ALAMF-1.0D0)*(ALAMF-1.0D0)*
     .                        (ALAMF-1.0D0)
           ENDIF
           SFIBF=SFIBR
           DFIBF=DFIBR
C
           DO ISTRS=1,NSTRS
            STILD(ISTRS)=STILD(ISTRS)+               ! S_(matrix+fibers)
     .       VANIS(IN,4)*SFIBF/ALAMF*VANIX(ISTRS)
           ENDDO
C
           DO ISTRS=1,NSTRS
            DO JSTRS=ISTRS,NSTRS
             DMATX(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)+  ! C_(matrix+fibers)
     .        2.0D0*VANIS(IN,4)*(DFIBF-SFIBF/ALAMF)/(ALAMF*ALAMF)*
     .        VANIX(ISTRS)*VANIX(JSTRS)
            ENDDO
           ENDDO
          ENDDO
         ENDIF                      ! nfibm.eq.2
        ENDIF                       ! nfibn.gt.0
C
        IF(IDV.EQ.0.AND.
     .     NMODI.NE.1) THEN         ! dev-vol decompos. is considered
         IDV=1
         ATILD=0.0D0
         DO ISTRS=1,NSTRS           ! auxiliar term (A-tilde)
          ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
          IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .     ATILD=ATILD+STILD(ISTRS)*RCGTT(ISTRS)
         ENDDO
         DO ISTRS=1,NSTRS           ! auxiliar term (B-tilde)
          BTILD(ISTRS)=0.0D0
          DO JSTRS=1,NSTRS
           I1=ISTRS
           I2=JSTRS
           IF(JSTRS.LT.ISTRS) THEN  ! only upper triangle of D is used
            I1=JSTRS
            I2=ISTRS
           ENDIF
           BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
           IF(I2.EQ.3.OR.I2.EQ.5.OR.I2.EQ.6)
     .      BTILD(ISTRS)=BTILD(ISTRS)+DMATX(I1,I2)*RCGTT(I2)
          ENDDO
         ENDDO
         AHATT=0.0D0
         DO ISTRS=1,NSTRS           ! auxiliar term (A-hat)
          AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
          IF(ISTRS.EQ.3.OR.ISTRS.EQ.5.OR.ISTRS.EQ.6)
     .     AHATT=AHATT+BTILD(ISTRS)*RCGTT(ISTRS)
         ENDDO
C
         DO ISTRS=1,NSTRS           ! deviatoric correction terms
          DO JSTRS=ISTRS,NSTRS
           DMATX(ISTRS,JSTRS)=FACTJ*(FACTJ*DMATX(ISTRS,JSTRS)-
     . 1.0D0/3.0D0*RCGTI(ISTRS)*(FACTJ*BTILD(JSTRS)+2.0D0*STILD(JSTRS))-
     . 1.0D0/3.0D0*(FACTJ*BTILD(ISTRS)+2.0D0*STILD(ISTRS))*RCGTI(JSTRS)-
     . 2.0D0/3.0D0*ATILD/DETJC*DCOFA(ISTRS,JSTRS)+
     . 1.0D0/9.0D0*(8.0D0*ATILD+FACTJ*AHATT)*RCGTI(ISTRS)*RCGTI(JSTRS))
          ENDDO
         ENDDO
C
         DO ISTRS=1,NSTRS
          STILD(ISTRS)=             ! dev. correction
     .     FACTJ*(STILD(ISTRS)-1.0D0/3.0D0*ATILD*RCGTI(ISTRS))
         ENDDO
        ENDIF                       ! idv.eq.0.and.nmodi.ne.1
C
        FVOL1=DETJC-1.0D0           ! incompressibility model: defaults
        FVOL2=1.0D0
        IF(NMODI.EQ.2) THEN
         FVOL1=0.5D0*(DETJM-1.0D0)/DETJM
         FVOL2=0.25D0/(DETJC*DETJM)
        ENDIF
        IF(NMODI.EQ.3) THEN
         FVOL1=0.5D0*DLOG(DETJM)/DETJC
         FVOL2=0.25D0*(1.0D0-2.0D0*DLOG(DETJM))/(DETJC*DETJC)
        ENDIF
        IF(NMODI.EQ.4) THEN
         FVOL1=0.25D0*(1.0D0-1.0D0/DETJC)
         FVOL2=0.25D0/(DETJC*DETJC)
        ENDIF
        DO ISTRS=1,NSTRS            ! volumetric terms
         DO JSTRS=ISTRS,NSTRS
          DMATX(ISTRS,JSTRS)=DMATX(ISTRS,JSTRS)+4.0D0*CPENI*
     .                     (DCOFA(ISTRS,JSTRS)*FVOL1+
     .                      DETJC*RCGTI(ISTRS)*DETJC*RCGTI(JSTRS)*FVOL2)
         ENDDO
        ENDDO                       ! SIGMA corrected in calcst.f
       ENDIF                        ! ifren.eq.51....
C
       DO ISTRS=1,NSTRS
        DO JSTRS=ISTRS,NSTRS
         DMATX(JSTRS,ISTRS)=DMATX(ISTRS,JSTRS)
        ENDDO
       ENDDO
      ENDIF                         ! istan.eq.1
C
      RETURN
#endif
      END

c
c tension plana para Mooney-Rivilin cuando no se considera la
c descomposicion F-tilde=J^(-1/3)F (lo que se hizo para la tesis
c de Cristian Nunez)
c
c OJO: implementar tension plana considerando la descomposicion
c      de F-tilde es complicado!
c
c       IF(NTYPE.EQ.1) THEN                   ! plane stress
c        call runend('ERROR: plane stress - MR-tilde model not implem.')
c        IF(NMODI.NE.1)
c    .    CALL RUNEND('ERROR: NMODI NE 1 IN PLANE STRESS - IFREN=51')
c        XXX=RCGTT(1)*RCGTT(2)-RCGTT(3)*RCGTT(3)   ! only for NO F-tilde
c        YYY=RCGTT(1)*RCGTT(2)-RCGTT(3)*RCGTT(3)-3.0D0+
c    .       (RCGTT(1)+RCGTT(2))*(RCGTT(1)+RCGTT(2))-
c    .       3.0D0*(RCGTT(1)+RCGTT(2))
c        AAA=-2.0D0*CPOL1-2.0D0*CPOL2*(RCGTT(1)+RCGTT(2))+
c    .        2.0D0*CPENI*XXX-
c    .        2.0D0*CPOL3*YYY
c        BBB=2.0D0*CPENI*XXX*XXX+    ! AAA/BBB=C_33 obtained from S_33=0
c    .       4.0D0*CPOL3*(RCGTT(1)+RCGTT(2))
C
c        DETJC=DETJC*AAA/BBB                  ! redefines det(C)
c        TRACC=TRACC+AAA/BBB                  ! redefines tr(C)
c        SEINC=0.5D0*(TRACC*TRACC-            ! redefines 2nd invariant
c    .                RCGTT(1)*RCGTT(1)-RCGTT(2)*RCGTT(2)-
c    .                2.0D0*RCGTT(3)*RCGTT(3)-AAA/BBB*AAA/BBB)
C
c        PS1=D1+4.0D0*CPENI*((DETJC-1.0D0)*RCGTT(2)+      ! 2dS_11/dC_33
c    .                        DETJC*RCGTI(1)*DETJC*AAA/BBB)+
c    .       4.0D0*CPOL3*(2.0D0*RCGTT(1)+3.0D0*RCGTT(2)+2.0D0*AAA/BBB-
c    .                    3.0D0)
c        PS2=D1+4.0D0*CPENI*((DETJC-1.0D0)*RCGTT(1)+      ! 2dS_22/dC_33
c    .                        DETJC*RCGTI(2)*DETJC*AAA/BBB)+
c    .       4.0D0*CPOL3*(3.0D0*RCGTT(1)+2.0D0*RCGTT(2)+2.0D0*AAA/BBB-
c    .                    3.0D0)
c        PS3=-4.0D0*CPENI*((DETJC-1.0D0)*RCGTT(3)-        ! 2dS_12/dC_33
c    .                      DETJC*RCGTI(3)*DETJC*AAA/BBB)-
c    .       4.0D0*CPOL3*RCGTT(3)
C
c        AAA11=-2.0D0*CPOL2+2.0D0*CPENI*RCGTT(2)-
c    .          2.0D0*CPOL3*(RCGTT(2)+2.0D0*(RCGTT(1)+RCGTT(2))-3.0D0)
c        AAA22=-2.0D0*CPOL2+2.0D0*CPENI*RCGTT(1)-
c    .          2.0D0*CPOL3*(RCGTT(1)+2.0D0*(RCGTT(1)+RCGTT(2))-3.0D0)
c        AAA12=            -2.0D0*CPENI*RCGTT(3)+
c    .                      4.0D0*CPOL3*RCGTT(3)
c        BBB11= 4.0D0*CPENI*XXX*RCGTT(2)+
c    .          4.0D0*CPOL3
c        BBB22= 4.0D0*CPENI*XXX*RCGTT(1)+
c    .          4.0D0*CPOL3
c        BBB12=-4.0D0*CPENI*XXX*RCGTT(3)
C
c        D11=(AAA11*BBB-AAA*BBB11)/(BBB*BBB)              ! dC_33/dC_11
c        D22=(AAA22*BBB-AAA*BBB22)/(BBB*BBB)              ! dC_33/dC_22
c        D12=(AAA12*BBB-AAA*BBB12)/(BBB*BBB)              ! dC_33/dC_12
C
c        DMATX(1,1)=4.0D0*CPENI*DETJC*RCGTI(1)*DETJC*RCGTI(1)+
c    .              PS1*D11
c        DMATX(1,2)=D1+4.0D0*CPENI*((DETJC-1.0D0)*AAA/BBB+
c    .                               DETJC*RCGTI(1)*DETJC*RCGTI(2))+
c    .              PS1*D22
c        DMATX(1,3)=4.0D0*CPENI*DETJC*RCGTI(1)*DETJC*RCGTI(3)+
c    .              PS1*D12
c        DMATX(2,2)=4.0D0*CPENI*DETJC*RCGTI(2)*DETJC*RCGTI(2)+
c    .              PS2*D22
c        DMATX(2,3)=4.0D0*CPENI*DETJC*RCGTI(2)*DETJC*RCGTI(3)+
c    .              PS2*D12
c        DMATX(3,3)=-D1-4.0D0*CPENI*((DETJC-1.0D0)*AAA/BBB+
c    .                                DETJC*RCGTI(3)*DETJC*RCGTI(3))+
c    .              PS3*D12
C
c        POISC1=0.0D0
c        POISC2=0.0D0
c        POISC3=0.0D0
c       ELSE                                  ! plane strain, axis. & 3D
