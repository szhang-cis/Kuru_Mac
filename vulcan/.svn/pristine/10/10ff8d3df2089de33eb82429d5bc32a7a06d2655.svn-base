      SUBROUTINE FRMC1DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                   GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                   STRSGT,STRA0T,STRS0T,
     .                   BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                   SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   ELEL1T,ELEL2T,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     ADVECTIVE PHASE-CHANGE FOR NDIME=1
C     (FOR ELEMENT 5, WHEN IT HAS TWO OR MORE PHASES)
C    
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
      DIMENSION CARTDT(NDIMET,NNODLT,*), DVOLUT(*),
     .          EHISTT(NHISTT,*),        ELCODT(NDIMET,*),
     .          ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*), 
     .          LNODST(*),               PROPST(*),
     .          RMAT1T(NDIMET,*),        SHAPET(NNODLT,*),
     .          STRANT(NSTR1T,*),        STRSGT(NSTR1T,*),
     .          STRA0T(NSTR1T,*),        STRS0T(NSTR1T,*)
      DIMENSION BMATXT(NSTR1T,*),        BMSIGT(*),
     .          DESIGT(*),               DMATXT(NSTR1T,*),
     .          DSTRAT(*),               PRESGT(*),
     .          SGTOTT(*),               SIGMAT(*),
     .          TSTRAT(*),               XJACMT(NDIMET,*),
     .          ELELTT(*),               VELCMT(*)
C
      DIMENSION DVOLDT(*),        SHAPDT(NNODLT,*),
     .          POSGPT(NDIMET,*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .                            CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION ELEL1T(*),        ELEL2T(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*)
      DIMENSION CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*)
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
        CALL PHCC01T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .               GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .               STRSGT,STRA0T,STRS0T,
     .               BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .               SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .               SHAPDT,DVOLDT,
     .               POSGPT,WEIGPT,
     .               DERIDT,
     .               CARDDT,GPCDDT,
     .               ELEL1T,ELEL2T,
     .                frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,
     .               ADVEMT,TEINIT,FPCHLT)
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.3) THEN
       CALL RUNENDT('FRMC1DT:ERROR DETECTED WITH NDIME=1')
       RETURN
      ENDIF
C
      END
