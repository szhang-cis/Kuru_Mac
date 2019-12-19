      SUBROUTINE MAMF1DT(DVOLUT,PROPST,SHAPET,WSTIFT,EHISTT,SHAPDT,
     .                   DVOLDT,POSGPT,ELDIST,ELCODT,WEIGPT,DERIDT,
     .                   XJACMT,CARDDT,GPCDDT,LNODST,VELCMT,WSTI1T,
     .                   WSTI2T,DISIMT,HACHET,WHADEL,CENTRO,ADVEMT,
     .                   TEINIT,FPCHLT,DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE COMPUTES PHASE CHANGE MATRIX (LATENT HEAT 
C     CONTRIBUTION INTO THE JACOBIAN MATRIX FOR ELEMENT NO. 5 ) FOR 
C     NDIME=1
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
      DIMENSION DVOLUT(*),        PROPST(*),
     .          SHAPET(NNODLT,*)
      DIMENSION DVOLDT(*),        POSGPT(NDIMET,*),
     .          SHAPDT(NNODLT,*)
      DIMENSION WSTIFT(*),        EHISTT(NHISTT,*),
     .          ELDIST(NDOFCT,*), ELCODT(NDIMET,*)
      DIMENSION VELCMT(*),        WSTI1T(*),
     .          WSTI2T(*)
      DIMENSION WEIGPT(*),        DERIDT(NDIMET,NNODLT,*),
     .          XJACMT(NDIMET,*), CARDDT(NDIMET,NNODLT,*),
     .          GPCDDT(NDIMET,*)
      DIMENSION LNODST(*),        DISIMT(*)
      DIMENSION HACHET(NNODLT),   WHADEL(NNODLT,*),
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
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
        CALL LAHE01T(POSGPT,SHAPDT,DVOLDT,PROPST,EHISTT,ELDIST,SHAPET,
     .               WSTIFT,VELCMT,WSTI1T,WSTI2T,
     .                frrrc, frrr1, fdvol,
     .               WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .               DISIMT,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,
     .               DVOLIT)
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
