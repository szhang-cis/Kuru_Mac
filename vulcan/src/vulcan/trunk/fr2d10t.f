      SUBROUTINE FR2D10T(PROPST,VELCMT,ELELMT,DVOLUT,SHAPET,EHISTT,
     .                   ELDIST,DISIMT,ELELTT,SHAPDT,DVOLDT,GPCODT,
     .                   POSGPT,ELCODT,WEIGPT,DERIDT,XJACMT,CARDDT,
     .                   GPCDDT,LNODST,ELEL1T,ELEL2T,INPCCT,TEMPLT,
     .                   HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT,
     .                   DVOLIT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE "THERMAL DYNAMIC FORCES" DUE TO 
C     PHASE-CHANGE WHEN IT HAS TWO PHASES ONLY FOR 2D ELEMENTS
C
C     FOR THE CASE: NPHCHT=10
C    
C***********************************************************************
C***********************************************************************
C    CASE
C
C NDIMET   NNODET   NISOTT                   CASE              NPHCHT
C***********************************************************************
C                                   4                    3
C   2        4        1             *---------R----------*       10
C                                   !          B         !
C                                   !                    !
C                                   !                    !
C                                   S                    S
C                                   !A                   !B
C                                   !                    !
C                                   *---------R----------*
C                                   1          A         2
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
     .          CENTRO(NDIMET,*), ADVEMT(NDIMET,*),
     .          TEINIT(NDOFCT,*), FPCHLT(NFPCH,*),
     .          DVOLIT(*)
C
C**** CALCULATE THE INTERFACE COORDINATES AND THE SHAPE FUNCTIONS FOR 
C     BIPHASE ELEMENTS
C
      IF(INPCCT.EQ.1) THEN
       ELDI1=ELDIST(1,1)-VELCMT(1)*DTIMET
       ELDI2=ELDIST(1,2)-VELCMT(2)*DTIMET
       ELDI3=ELDIST(1,3)-VELCMT(3)*DTIMET
       ELDI4=ELDIST(1,4)-VELCMT(4)*DTIMET
      ELSE
       ELDI1=ELDIST(1,1)
       ELDI2=ELDIST(1,2)
       ELDI3=ELDIST(1,3)
       ELDI4=ELDIST(1,4)
      ENDIF
      CALL COIN01T(ELDI2,ELDI1,TEMPLT,RCOOA,ISOLI)
      CALL COIN01T(ELDI3,ELDI4,TEMPLT,RCOOB,ISOLI)
      CALL COIN01T(ELDI4,ELDI1,TEMPLT,SCOOA,ISOLI)
      CALL COIN01T(ELDI3,ELDI2,TEMPLT,SCOOB,ISOLI)
C
      r1=-1.0D0
      r2=rcooa
      r3=-1.0D0
      s1=-1.0D0
      s2=-1.0D0
      s3=scooa
      frrrc=r1
      frrr1=r2-r1
      frrr2=r3-r1
      fsssc=s1
      fsss1=s2-s1
      fsss2=s3-s1
      fdvol=frrr1*fsss2-fsss1*frrr2
      CALL PHCH04T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,
     .              frrrc, frrr1, frrr2, fsssc, fsss1, fsss2, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             POSGPT,ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
C
      r1=1.0D0
      r2=1.0D0
      r3=rcooa
      s1=-1.0D0
      s2=scoob
      s3=-1.0D0
      frrrc=r1
      frrr1=r2-r1
      frrr2=r3-r1
      fsssc=s1
      fsss1=s2-s1
      fsss2=s3-s1
      fdvol=frrr1*fsss2-fsss1*frrr2
      CALL PHCH04T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,
     .              frrrc, frrr1, frrr2, fsssc, fsss1, fsss2, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             POSGPT,ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
C
      r1=1.0D0
      r2=rcoob
      r3=1.0D0
      s1=1.0D0
      s2=1.0D0
      s3=scoob
      frrrc=r1
      frrr1=r2-r1
      frrr2=r3-r1
      fsssc=s1
      fsss1=s2-s1
      fsss2=s3-s1
      fdvol=frrr1*fsss2-fsss1*frrr2
      CALL PHCH04T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,
     .              frrrc, frrr1, frrr2, fsssc, fsss1, fsss2, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             POSGPT,ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
C
      r1=-1.0D0
      r2=-1.0D0
      r3=rcoob
      s1=1.0D0
      s2=scooa
      s3=1.0D0
      frrrc=r1
      frrr1=r2-r1
      frrr2=r3-r1
      fsssc=s1
      fsss1=s2-s1
      fsss2=s3-s1
      fdvol=frrr1*fsss2-fsss1*frrr2
      CALL PHCH04T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,
     .              frrrc, frrr1, frrr2, fsssc, fsss1, fsss2, fdvol,
     .             WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,ELCODT,LNODST,
     .             POSGPT,ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
C
      r1=rcooa
      r2=1.0D0
      r3=rcoob
      r4=-1.0D0
      s1=-1.0D0
      s2=scoob
      s3=1.0D0
      s4=scooa
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
      CALL PHCH05T(PROPST,VELCMT,ELELMT,SHAPDT,DVOLDT,EHISTT,ELDIST,
     .             DISIMT,ELELTT,WEIGPT,DERIDT,XJACMT,CARDDT,GPCDDT,
     .             ELCODT,LNODST,POSGPT,
     .              frrrc, frrr1, frrr2, frr12, fsssc, fsss1, fsss2,
     .              fss12, fdvoc, fdvo1, fdvo2,
     .             ELEL1T,ELEL2T,INPCCT,HACHET,WHADEL,CENTRO,
     .             ADVEMT,TEINIT,FPCHLT,DVOLIT)
C
      RETURN
C
      END
