      SUBROUTINE EQLOAB(BSBAR,BMSIG,CARTD,DVOLU,ELDIS,GPCOD,
     .                  LARGE,NDIME,NDOFC,NDOFN,NEVAB,NNODE,NSTRE,
     .                  NTYPE,SHAPE,SGTOT,
     .                  XJACM,XJACI,XJA3M,XJA3I,DETJM,DEFRA,
     .                  NSTR1,NSTRS,BMATX,
     .                  XJACN,XJANI,XJA3N,XJ3NI,DETJN)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE EQUIVALENT NODAL FORCES DUE TO
C     THE EQUILIBRATED STRESS
C
C     Notes:
C
C     LARGE flag independent computation of BMSIG is attempted for large
C     strains
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
C
      DIMENSION BSBAR(NSTR1,*), BMSIG(*),      CARTD(NDIME,*), 
     .          ELDIS(NDOFC,*), SHAPE(*),      GPCOD(*),       SGTOT(*),
     .          XJACM(NDIME,*), XJACN(NDIME,*),
     .          XJACI(NDIME,*), XJANI(NDIME,*)
      DIMENSION BMATX(NSTRE,*)
C
      NEVAL=NNODE*NDOFN      ! local NEVAB
C
      IF(LARGE.EQ.0) THEN
       DO ISTRS=1,NSTRE
        DO IEVAB=1,NEVAL
         BMSIG(IEVAB)=BMSIG(IEVAB)+BSBAR(ISTRS,IEVAB)*SGTOT(ISTRS)*DVOLU
        ENDDO
       ENDDO
C
      ELSE                    ! large strains
C
       KPARB=INT(PPARB)
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
   27  CALL BLARGE(BMATX,CARTD,ELDIS,GPCOD,LARGE,NDIME,NDOFC,
     .             NNODE,NSTRE,NTYPE,SHAPE,XJACM)
       GO TO 34
C
   34  CONTINUE
C
       DO ISTRS=1,NSTRE
        DO IEVAB=1,NEVAL
         BMSIG(IEVAB)=BMSIG(IEVAB)+BMATX(ISTRS,IEVAB)*SGTOT(ISTRS)*DVOLU
        ENDDO
       ENDDO
C
      ENDIF                   ! large.eq.0
C
      RETURN
      END
