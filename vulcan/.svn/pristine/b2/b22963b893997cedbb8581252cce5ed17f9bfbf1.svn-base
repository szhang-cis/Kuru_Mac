      SUBROUTINE FRMX2DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
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
C     ADVECTIVE PHASE-CHANGE FOR NDIME=2
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
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*)
C
      IF(NNODLT.EQ.3) THEN
       nr=nsubd(1,1)
       ns=nr
       dr=1.00D0/nr
       ds=1.00D0/ns
       do ic=1,2
        if(ic.eq.1) then
         nrr=nr
        else
         nrr=nr-1
        endif
        rir=0.00D0
        rirm=0.00D0
        do ir=1,nrr
         rir=rirm
         rirm=rir+dr
         if(ic.eq.1) then
          rx=rir
          rxr=rirm
         else
          rx=rirm
          rxr=rir
         endif
         r1=rx
         r2=rxr
         r3=rx
         sis=0.00D0
         sism=0.00D0
         do is=1,nrr-ir+1
          sis=sism
          sism=sis+ds
          if(ic.eq.1) then
           sx=sis
           sxs=sism
          else
           sx=sism
           sxs=sis
          endif
          s1=sx
          s2=sx
          s3=sxs
          frrrc=r1
          frrr1=r2-r1
          frrr2=r3-r1
          fsssc=s1
          fsss1=s2-s1
          fsss2=s3-s1
          fdvol=frrr1*fsss2-fsss1*frrr2
          CALL PHCX10T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                 GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                 STRSGT,STRA0T,STRS0T,
     .                 BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                 SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                 SHAPDT,DVOLDT,
     .                 POSGPT,WEIGPT,
     .                 DERIDT,
     .                 CARDDT,GPCDDT,
     .                 ELEL1T,ELEL2T,
     .                  frrrc, frrr1, frrr2,
     .                  fsssc, fsss1, fsss2,
     .                  fdvol,HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,
     .                 FPCHLT)
         enddo
        enddo
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.4) THEN
       nr=nsubd(1,1)
       ns=nsubd(2,1)
       dr=2.00D0/nr
       ds=2.00D0/ns
       rir=-1.00D0
       rirm=-1.00D0
       do ir=1,nr
        rir=rirm
        rirm=rir+dr
        r1=rir
        r2=rirm
        r3=rirm
        r4=rir
        sis=-1.00D0
        sism=-1.00D0
        do is=1,ns
         sis=sism
         sism=sis+ds
         s1=sis
         s2=sis
         s3=sism
         s4=sism
         frrrc=0.25D0*(r1+r2+r3+r4)
         frrr1=0.25D0*(-r1+r2+r3-r4)
         frrr2=0.25D0*(-r1-r2+r3+r4)
         frr12=0.25D0*(r1-r2+r3-r4)
         fsssc=0.25D0*(s1+s2+s3+s4)
         fsss1=0.25D0*(-s1+s2+s3-s4)
         fsss2=0.25D0*(-s1-s2+s3+s4)
         fss12=0.25D0*(s1-s2+s3-s4)
         fdvoc=frrr1*fsss2-frrr2*fsss1
         fdvo1=frrr1*fss12-fsss1*frr12
         fdvo2=fsss2*frr12-frrr2*fss12
         CALL PHCX05T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                STRSGT,STRA0T,STRS0T,
     .                BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                SHAPDT,DVOLDT,
     .                POSGPT,WEIGPT,
     .                DERIDT,
     .                CARDDT,GPCDDT,
     .                ELEL1T,ELEL2T,
     .                 frrrc, frrr1, frrr2, frr12,
     .                 fsssc, fsss1, fsss2, fss12,
     .                 fdvoc, fdvo1, fdvo2,HACHET,WHADEL,CENTRO,
     .                ADVEMT,TEINIT,FPCHLT)
        enddo
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.6) THEN
       CALL RUNENDT('FRMC2DT:ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.8) THEN
       CALL RUNENDT('FRMC2DT:ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.9) THEN
       CALL RUNENDT('FRMC2DT:ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
C
      END
