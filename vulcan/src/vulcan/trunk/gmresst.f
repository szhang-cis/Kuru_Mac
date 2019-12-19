C=======================================================================
      subroutine gmresst(CSTIF,NEQNS,BGMRE,XGMRE,WGMRE,
     .                   UNKNW,RHSID,DIAGO,LNUEQ,LNODS,NDOFN,
     .                   NELEM,NEVAB,NNODE,IPASS,NTOTV,NPOIN,
     .                   NLAST,IDSK2,LPONT,
     .                   MITGM,MKRYL,IPGMR,TOLGM,
     .                   LURES)
C-----------------------------------------------------------------------
*     Description : Solves non-symmetric linear system of equations by
*                   preconditioned GMRES method
*
*     Parameters :  CSTIF   IN    elemental matrix
*                   NEQNS   IN    number of equations
*                   BGMRE   IN    constant vector
*                   XGMRE   OUT   solution  of  a*x = b
*                   TOLGM   IN    solver tolerance
*
*     Version : gmres_so  2.00  (Oct - 93)
C-----------------------------------------------------------------------
C
      REAL*8  CSTIF(*),BGMRE(*),XGMRE(*),WGMRE(*),UNKNW(*),
     .        RHSID(*),DIAGO(*)
C
      INTEGER*4 LNUEQ(NTOTV),LNODS(*),NDOFN,NELEM,NEVAB,NNODE,IPASS,
     .          NTOTV,
     .          NEQNS,NPOIN,
     .          NLAST,IDSK2,LPONT(NTOTV),LURES
C
      INTEGER*4 MITGM,MKRYL,IPGMR
      REAL*8 TOLGM
C-----------------------------------------------------------------------
C
C**** IDENTIFY THE DESTINATION OF THE EQUATIONS INTO PROFILE
C     NEQNS,NLAST,LNUEQ,LPONT (as in skygrat.f)
C
C     Note: NLAST & LPONT are not needed but they are read here because
C           skydest.f was used in soladdt.f
C
      REWIND IDSK2
      READ(IDSK2) NEQNS,NLAST,LNUEQ,LPONT
C-----------------------------------------------------------------------
C
      call gmrespt(CSTIF,DIAGO,NEQNS,WGMRE,XGMRE,MKRYL,MITGM,
     .                   TOLGM,BGMRE,IPGMR,UNKNW,RHSID,LNUEQ,LNODS,
     .                   NDOFN,NELEM,NEVAB,NNODE,IPASS,NTOTV,NPOIN,
     .                   LURES)                         ! added variable
C
      do itotv=1,ntotv
       idest=lnueq(itotv)
       if(idest.gt.0) unknw(itotv)=xgmre(idest)
      end do
C
      end
