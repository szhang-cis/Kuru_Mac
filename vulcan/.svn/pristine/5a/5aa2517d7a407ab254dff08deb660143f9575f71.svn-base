      SUBROUTINE STMC1DT(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                   PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                   DMATXT,SIGMAT,XJACMT,
     .                   SHAPDT,DVOLDT,
     .                   POSGPT,WEIGPT,
     .                   DERIDT,
     .                   CARDDT,GPCDDT,
     .                   VELCMT,WSTI1T,
     .                   WSTI2T,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .                   LNODST,ELCODT,FPCHLT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES THE ADVECTIVE PHASE-CHANGE MATRIX 
C     (CONVECTIVE LATENT HEAT CONTRIBUTION INTO THE JACOBIAN MATRIX FOR
C     ELEMENT NO. 5 ) FOR NDIME=1
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
     .          EHISTT(NHISTT,*),        ELDIST(NDOFCT,*),
     .          EPMTXT(*),               GPCODT(NDIMET,*),
     .          PROPST(*),               SHAPET(NNODLT,*),
     .          STRSGT(NSTR1T,*),        ESTIFT(*),
     .          HSTIFT(NEVABT,NNODET)
      DIMENSION BMATXT(NSTR1T,*),        DMATXT(NSTR1T,*),
     .          SIGMAT(*),               XJACMT(NDIMET,*)
C
      DIMENSION DVOLDT(*),        POSGPT(NDIMET,*),
     .          SHAPDT(NNODLT,*)
C
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .                            CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*), VELCMT(*)
      DIMENSION WSTI1T(*),        WSTI2T(*)
      DIMENSION LNODST(*)
      DIMENSION ELCODT(NDIMET,*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*)
      DIMENSION CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*)
C
      IF(NNODLT.EQ.2) THEN
       nr=nsubd(1,2)
       dr=2.00D0/nr
       rir=-1.00D0
       rirm=-1.00D0
       do ir=1,nr
        rir=rirm
        rirm=rir+dr
        r1=rir
        r2=rirm
        frrrc=0.5D0*(r1+r2)
        frrr1=0.5D0*(r2-r1)
        fdvol=frrr1
        CALL LAHC01T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .               PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .               DMATXT,SIGMAT,XJACMT,
     .               SHAPDT,DVOLDT,
     .               POSGPT,WEIGPT,
     .               DERIDT,
     .               CARDDT,GPCDDT,
     .               VELCMT,WSTI1T,
     .               WSTI2T,
     .               LNODST,ELCODT,
     .                frrrc, frrr1, fdvol,HACHET,WHADEL,CENTRO,
     .               ADVEMT,TEINIT,FPCHLT)
       enddo
       RETURN    
      ENDIF 
C
      IF(NNODLT.EQ.3) THEN
       CALL RUNENDT('MAMF1DT:ERROR DETECTED WITH NDIME=1')
       RETURN    
      ENDIF
C
      END
