      SUBROUTINE ESEPAST(DISITT,DISPRT,IFFIXT,LNODST,REFORT,RLOADT,
     .                   TLOADT,WORK1T,
     .                   MATNOT,PROELT,PROPST,ELDATT,ELPRET,
     .                   ELVART,ELMATT,
     .                   DISTOT,COORDT,
     .                   ADVELT,TEMPIT,PREAST,TGAPST,DISPLT,FPCHAT,
     .                   LACTIT)
C***********************************************************************
C
C**** THIS ROUTINE PERFOMS THE SOLUTION FOR UNCOUPLED PROBLEM
C     (CONSTRAINED METHODS)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION DISITT(NTOTVT,*), DISPRT(NTOTVT,*), IFFIXT(NTOTVT,*),
     .          LNODST(NNODET,*), REFORT(NTOTVT,*), RLOADT(*),
     .          TLOADT(NTOTVT,*)
      DIMENSION WORK1T(*)
C
      DIMENSION MATNOT(NELEMT),
     .          PROELT(NPRELT,NGRUPT), PROPST(NPROPT,NMATST),
     .          ELDATT(NDATAT),        ELPRET(NPREVT),
     .          ELVART(NSTATT),        ELMATT(NMATXT)
      DIMENSION DISTOT(NTOTVT,3),      COORDT(NDIMET,NPOINT)
      DIMENSION ADVELT(NTOTVT*NDIMET), TEMPIT(NPOINT,2)
      DIMENSION PREAST(NPOINT),        TGAPST(NPOINT)
      DIMENSION DISPLT(NTOTVM),        FPCHAT(NFPCH,NPOINT),
     .          LACTIT(NELEMT)
C
C**** SOLVE EQUATIONS  ( FIRST PASS )
C
      CALL ANTESOT(REFORT(1,1),REFORT(1,2),DISITT(1,2),
     .             WORK1T(IQUAST(1)),WORK1T(IQUAST(2)))
C
      CALL SOLVERT(DISITT(1,1),REFORT(1,1),IFFIXT(1,1),
     .             IMPOST,KPORET,KRESLT,
     .             LNODST,NDOFNT,NELEMT,NEVABT,NNODET,NPOINT,
     .             NTOTVT,     1,IITERT,WORK1T,mitert,
     .             NPRELT,NGRUPT,NPROPT,NMATST,
     .             NDATAT,NPREVT,NSTATT,NMATXT,NDIMET,NTOTVM,NFPCH,
     .             MATNOT,PROELT,PROPST,ELDATT,ELPRET,
     .             ELVART,ELMATT,
     .             DISTOT,COORDT,
     .             ADVELT,TEMPIT,PREAST,TGAPST,DISPLT,FPCHAT,
     .             NACTIT,LACTIT)
C
      CALL POSTSOT(DISITT(1,1),DISITT(1,2),REFORT(1,1),REFORT(1,2),
     .             WORK1T(IQUAST(1)),WORK1T(IQUAST(2)))
C
      RETURN
      END
