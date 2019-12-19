      SUBROUTINE FRMF1DT(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,CARDDT,
     .                   GPCDDT,LNODST,ELEL1T,ELEL2T,INPCCT,HACHET,
     .                   WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE (FOR ELEMENT 5, WHEN IT HAS TWO OR MORE PHASES) FOR 
C     NDIME=1
C    
C***********************************************************************
C     EHIST(1) = Density
C     EHIST(2) = Temperature Derivative of Density
C     EHIST(3) = Specific Heat coefficient)
C     EHIST(4) = Temperature Derivative of the Specific Heat coefficient
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'   ! thermal-mechanical
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*),        VELCMT(*),
     .          ELELMT(*),        DVOLUT(*),
     .          SHAPET(NNODLT,*), EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), DISIMT(*),
     .          ELELTT(*)
      DIMENSION DVOLDT(*),        SHAPDT(NNODLT,*),
     .          POSGPT(NDIMET,*)
      DIMENSION GPCODT(NDIMET,*), ELCODT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION LNODST(*)
      DIMENSION ELEL1T(*),        ELEL2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*),        ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
      IF(NNODLT.EQ.2) THEN
       nr=nsubd(1,1)
       dr=2.00D0/nr
       rir=-1.00D0
       rirm=-1.00D0
       do ir=1,nr
        rir=rirm
        rirm=rir+dr
        r1=rir
        r2=rirm
        frrrc=0.50D0*(r1+r2)
        frrr1=0.50D0*(r2-r1)
        fdvol=frrr1
        CALL PHCH01T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,SHAPET,
     .               ELELMT,VELCMT,ELEL1T,ELEL2T,INPCCT,
     .                frrrc, frrr1, fdvol,
     .               WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .               DISIMT,ELELTT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .               FPCHLT,DVOLIT)
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.3) THEN
       CALL RUNENDT('FRMF1DT:ERROR DETECTED WITH NDIME=1')
       RETURN
      ENDIF
C
      END
