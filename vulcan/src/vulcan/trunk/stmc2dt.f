      SUBROUTINE STMC2DT(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
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
C     ELEMENT NO. 5 ) FOR NDIME=2
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
      IF(NNODLT.EQ.3) THEN
       nr=nsubd(1,2)
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
          CALL LAHC07T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                 PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                 DMATXT,SIGMAT,XJACMT,
     .                 SHAPDT,DVOLDT,
     .                 POSGPT,WEIGPT,
     .                 DERIDT,
     .                 CARDDT,GPCDDT,
     .                 VELCMT,WSTI1T,
     .                 WSTI2T,
     .                 LNODST,ELCODT,
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
       nr=nsubd(1,2)
       ns=nsubd(2,2)
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
C
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
         CALL LAHC05T(CARTDT,DVOLUT,EHISTT,ELDIST,EPMTXT,GPCODT,
     .                PROPST,SHAPET,STRSGT,ESTIFT,HSTIFT,BMATXT,
     .                DMATXT,SIGMAT,XJACMT,
     .                SHAPDT,DVOLDT,
     .                POSGPT,WEIGPT,
     .                DERIDT,
     .                CARDDT,GPCDDT,
     .                VELCMT,WSTI1T,
     .                WSTI2T,
     .                LNODST,ELCODT,
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
       CALL RUNENDT('STMC2DT: ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.8) THEN
       CALL RUNENDT('STMC2DT: ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.9) THEN
       CALL RUNENDT('STMC2DT: ERROR DETECTED WITH NDIME=2')
       RETURN
      ENDIF
      RETURN
C
      END
