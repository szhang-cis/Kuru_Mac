      SUBROUTINE KMATRB(BMATX,CARTD,DMATX,DVOLU,ESTIF,GPCOD,LARGE,
     .                  KSYMM,NDIME,NDOFN,NEVAB,NKOVA,NNODE,NSTRE,
     .                  NSTRS,NTYPE,SHAPE,NSTR1,NDOFC,ELDIS,EMATX,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .                  BSBAR,
     .                  XJACN,XJANI,XJA3N,XJ3NI,DETJN)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE MATERIAL STIFFNESS MATRIX
C
C     Notes:
C
C     LARGE flag independent computation of TSTRA is attempted for large
C     strains
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION BMATX(NSTRE,*), CARTD(NDIME,*), DMATX(NSTRS,*),
     .          ESTIF(*),       GPCOD(*),       SHAPE(*),
     .          EMATX(NEVAB,*)
      DIMENSION ELDIS(NDOFC,*), XJACM(NDIME,*), XJACN(NDIME,*),
     .                          XJACI(NDIME,*), XJANI(NDIME,*)
      DIMENSION BSBAR(NSTR1,*)
C
      NEVAL=NNODE*NDOFN      ! local NEVAB
C
C                        T
C**** PERFORME THE (B-BAR   D  B-BAR)   PRODUCT
C
      IF(LARGE.EQ.0) THEN
       CALL BDBCO1(BSBAR,DMATX,EMATX,KSYMM,NEVAB,NSTR1,NSTRS,NEVAL,
     .             NSTRE)
      ELSE                   ! large strains
       KPARB=INT(PPARB)      ! kparb > 0 here
       IF(KPARB.LT.0) KPARB=-KPARB
C
       GO TO (21,22,23,21,22,23,21,22,23,
     .        24,24,24,25,25,25,26,26,26,
     .        21,21,21,27,27,27), KPARB
C
C**** B-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   21  CALL BBMATXL(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,DEFRA,
     .              XJACN,XJA3N,DETJN,
     .              KPARB)
       GO TO 34
C
C**** B-SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   22  CALL BBMATYL(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,
     .              XJACN,XJA3N,DETJN)
       GO TO 34
C
C**** B-BAR-SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   23  CALL BBMATZL(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,DEFRA,
     .              XJACN,XJA3N,DETJN)
       GO TO 34
C
C**** OTHER B-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   24  CALL BBMATWL(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,
     .              XJACN,XJA3N,DETJN)
       GO TO 34
C
C**** B-e_SHEAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEGRATION)
C
   25  CALL BBMATYLS(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,
     .              XJACN,XJA3N,DETJN,
     .              XJACI,XJA3I,
     .              XJANI,XJ3NI)
       GO TO 34
C
C**** B-(e_SHEAR-BAR) (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEG.)
C
   26  CALL BBMATZLS(CARTD,ELDIS,GPCOD,
     .              NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .              NTYPE,SHAPE,NSTR1,NSTRS,BMATX,BSBAR,
     .              XJACM,XJA3M,DETJM,DEFRA,
     .              XJACN,XJA3N,DETJN,
     .              XJACI,XJA3I,
     .              XJANI,XJ3NI)
       GO TO 34
C
C**** J-BAR (SMOOTHING TECH., AVERAGE TECH. & SELECT. INTEG.)
C
   27  CALL BLARGE(BMATX,CARTD,ELDIS,GPCOD,
     .             LARGE,NDIME,NDOFC,NNODE,NSTRE,NTYPE,SHAPE,
     .             XJACM)
       GO TO 34
C
   34  CONTINUE
C                        T
C**** PERFORME THE (B-BAR   D  B-BAR)   PRODUCT
C
       CALL BDBCO1(BMATX,DMATX,EMATX,KSYMM,NEVAB,NSTRE,NSTRS,NEVAL,
     .             NSTRE)
      ENDIF             ! large.eq.0
C
C**** ADD THE CONTRIBUTION TO THE STIFFNESS MATRIX
C
      CALL SQTOTA(ESTIF,EMATX,DVOLU,NEVAB,KSYMM)
C
      RETURN
      END
