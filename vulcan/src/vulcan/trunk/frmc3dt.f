      SUBROUTINE FRMC3DT(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
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
C     ADVECTIVE PHASE-CHANGE (FOR ELEMENT 5, WHEN IT HAS TWO OR MORE
C     PHASES) FOR NDIME=3
C    
C***********************************************************************
C
C     Index of variables:
C
C     EHIST(1) = Density
C     EHIST(2) = Specific Heat coefficient
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
C**** TETRAHEDRON
C
      IF(NNODLT.EQ.4) THEN
C
       nr=nsubd(1,1)
       if(nr.eq.1) ng=1
       if(nr.eq.2) ng=7
       if(nr.gt.2) call runendt('error in frmc3dt with nnodlt=4')
c
       do ig=1,ng
        if(nr.eq.1) then
         x1=0.0
         x3=1.0
         y1=0.0
         y3=1.0
         z1=0.0
         z3=1.0
        endif
        if(nr.eq.2) then
         go to(10,20,30,40,50,60,70) ig
   10    continue
         x1=0.5
         x3=1.0
         y1=0.0
         y3=0.5
         z1=0.0
         z3=0.5
         go to 80
c
   20    continue
         x1=0.0
         x3=0.5
         y1=0.5
         y3=1.0
         z1=0.0
         z3=0.5
         go to 80
c
   30    continue
         x1=0.0
         x3=0.5
         y1=0.0
         y3=0.5
         z1=0.5
         z3=1.0
         go to 80
c
   40    continue
         x1=0.5
         x3=0.0
         y1=0.0
         y3=0.5
         z1=0.0
         z3=0.5
         go to 80
c
   50    continue
         x1=0.0
         x3=0.5
         y1=0.5
         y3=0.0
         z1=0.0
         z3=0.5
         go to 80
c
   60    continue
         x1=0.0
         x3=0.5
         y1=0.0
         y3=0.5
         z1=0.5
         z3=0.0
         go to 80
c
   70    continue
         x1=0.5
         x3=0.0
         y1=0.0
         y3=0.5
         z1=0.5
         z3=0.0
c
   80    continue
        endif
c
        x2=0.5*(x1+x3)
        y2=0.5*(y1+y3)
        z2=0.5*(z1+z3)
c
        do ic=1,7
         go to(1,2,3,4,5,6,7) ic
    1    continue
         r1=x2
         r2=x3
         r3=x2
         r4=x2
         s1=y1
         s2=y1
         s3=y2
         s4=y1
         t1=z1
         t2=z1
         t3=z1
         t4=z2
         go to 8
c
    2    continue
         r1=x1
         r2=x2
         r3=x1
         r4=x1
         s1=y2
         s2=y2
         s3=y3
         s4=y2
         t1=z1
         t2=z1
         t3=z1
         t4=z2
         go to 8
c
    3    continue
         r1=x1
         r2=x2
         r3=x1
         r4=x1
         s1=y1
         s2=y1
         s3=y2
         s4=y1
         t1=z2
         t2=z2
         t3=z2
         t4=z3
         go to 8
c
    4    continue
         r1=x2
         r2=x2
         r3=x1
         r4=x2
         s1=y1
         s2=y2
         s3=y1
         s4=y1
         t1=z1
         t2=z1
         t3=z1
         t4=z2
         go to 8
c
    5    continue
         r1=x1
         r2=x1
         r3=x1
         r4=x2
         s1=y2
         s2=y2
         s3=y1
         s4=y2
         t1=z1
         t2=z2
         t3=z1
         t4=z1
         go to 8
c
    6    continue
         r1=x1
         r2=x2
         r3=x1
         r4=x1
         s1=y1
         s2=y1
         s3=y1
         s4=y2
         t1=z2
         t2=z2
         t3=z1
         t4=z2
         go to 8
c
    7    continue
         r1=x2
         r2=x1
         r3=x2
         r4=x1
         s1=y1
         s2=y2
         s3=y2
         s4=y1
         t1=z2
         t2=z2
         t3=z1
         t4=z1
c
    8    continue
c
         frrrc=r1
         frrr1=r2-r1
         frrr2=r3-r1
         frrr3=r4-r1
         fsssc=s1
         fsss1=s2-s1
         fsss2=s3-s1
         fsss3=s4-s1
         ftttc=t1
         fttt1=t2-t1
         fttt2=t3-t1
         fttt3=t4-t1
         fdvol=frrr1*(fsss2*fttt3-fsss3*fttt2)-
     .         fsss1*(frrr2*fttt3-frrr3*fttt2)+
     .         fttt1*(frrr2*fsss3-frrr3*fsss2)
C
         CALL PHCC12T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                STRSGT,STRA0T,STRS0T,
     .                BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                SHAPDT,DVOLDT,
     .                POSGPT,WEIGPT,
     .                DERIDT,
     .                CARDDT,GPCDDT,
     .                ELEL1T,ELEL2T,
     .                 frrrc, frrr1, frrr2, frrr3,
     .                 fsssc, fsss1, fsss2, fsss3, 
     .                 ftttc, fttt1, fttt2, fttt3, 
     .                 fdvol,
     .                HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
        enddo
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.8) THEN
       nr=nsubd(1,1)
       ns=nsubd(2,1)
       nt=nsubd(3,1)
       dr=2.00D0/nr
       ds=2.00D0/ns
       dt=2.00D0/nt
       rir=-1.00D0
       rirm=-1.00D0
       do ir=1,nr
        rir=rirm
        rirm=rir+dr
        r1=rir
        r2=rirm
        r3=rirm
        r4=rir
        r5=rir
        r6=rirm
        r7=rirm
        r8=rir
        sis=-1.00D0
        sism=-1.00D0
        do is=1,ns
         sis=sism
         sism=sis+ds
         s1=sis
         s2=sis
         s3=sism
         s4=sism
         s5=sis
         s6=sis
         s7=sism
         s8=sism
         tit=-1.00D0
         titm=-1.00D0
         do it=1,nt
          tit=titm
          titm=tit+dt
          t1=tit
          t2=tit
          t3=tit
          t4=tit
          t5=titm
          t6=titm
          t7=titm
          t8=titm
C
          frrrc=0.125D0*( r1+r2+r3+r4+r5+r6+r7+r8)
          frrr1=0.125D0*(-r1+r2+r3-r4-r5+r6+r7-r8)
          frrr2=0.125D0*(-r1-r2+r3+r4-r5-r6+r7+r8)
          frrr3=0.125D0*(-r1-r2-r3-r4+r5+r6+r7+r8)
          frr12=0.125D0*( r1-r2+r3-r4+r5-r6+r7-r8)
          frr13=0.125D0*( r1-r2-r3+r4-r5+r6+r7-r8)
          frr23=0.125D0*( r1+r2-r3-r4-r5-r6+r7+r8)
          fr123=0.125D0*(-r1+r2-r3+r4+r5-r6+r7-r8)
c
          fsssc=0.125D0*( s1+s2+s3+s4+s5+s6+s7+s8)
          fsss1=0.125D0*(-s1+s2+s3-s4-s5+s6+s7-s8)
          fsss2=0.125D0*(-s1-s2+s3+s4-s5-s6+s7+s8)
          fsss3=0.125D0*(-s1-s2-s3-s4+s5+s6+s7+s8)
          fss12=0.125D0*( s1-s2+s3-s4+s5-s6+s7-s8)
          fss13=0.125D0*( s1-s2-s3+s4-s5+s6+s7-s8)
          fss23=0.125D0*( s1+s2-s3-s4-s5-s6+s7+s8)
          fs123=0.125D0*(-s1+s2-s3+s4+s5-s6+s7-s8)
c
          ftttc=0.125D0*( t1+t2+t3+t4+t5+t6+t7+t8)
          fttt1=0.125D0*(-t1+t2+t3-t4-t5+t6+t7-t8)
          fttt2=0.125D0*(-t1-t2+t3+t4-t5-t6+t7+t8)
          fttt3=0.125D0*(-t1-t2-t3-t4+t5+t6+t7+t8)
          ftt12=0.125D0*( t1-t2+t3-t4+t5-t6+t7-t8)
          ftt13=0.125D0*( t1-t2-t3+t4-t5+t6+t7-t8)
          ftt23=0.125D0*( t1+t2-t3-t4-t5-t6+t7+t8)
          ft123=0.125D0*(-t1+t2-t3+t4+t5-t6+t7-t8)
C
          fdvoc=   FRRR1*FSSS2*FTTT3 - FRRR1*FSSS3*FTTT2 - 
     .             FSSS1*FRRR2*FTTT3 + FSSS1*FRRR3*FTTT2 + 
     .             FTTT1*FRRR2*FSSS3 - FTTT1*FRRR3*FSSS2

          fdvo1=   FRRR1*FSS12*FTTT3 - FRRR1*FSS13*FTTT2 - 
     .             FRRR1*FTT12*FSSS3 + FRRR1*FTT13*FSSS2 - 
     .             FRR12*FSSS1*FTTT3 + FRR12*FTTT1*FSSS3 + 
     .             FRR13*FSSS1*FTTT2 - FRR13*FTTT1*FSSS2 + 
     .             FSSS1*FTT12*FRRR3 - FSSS1*FTT13*FRRR2 - 
     .             FSS12*FTTT1*FRRR3 + FSS13*FTTT1*FRRR2

          fdvo2= - FRRR1*FSS23*FTTT2 + FRRR1*FTT23*FSSS2 + 
     .             FRR12*FSSS2*FTTT3 - FRR12*FSSS3*FTTT2 + 
     .             FRR23*FSSS1*FTTT2 - FRR23*FTTT1*FSSS2 - 
     .             FSSS1*FTT23*FRRR2 - FSS12*FRRR2*FTTT3 + 
     .             FSS12*FRRR3*FTTT2 + FSS23*FTTT1*FRRR2 + 
     .             FTT12*FRRR2*FSSS3 - FTT12*FRRR3*FSSS2

          fdvo3=   FRRR1*FSS23*FTTT3 - FRRR1*FTT23*FSSS3 + 
     .             FRR13*FSSS2*FTTT3 - FRR13*FSSS3*FTTT2 - 
     .             FRR23*FSSS1*FTTT3 + FRR23*FTTT1*FSSS3 + 
     .             FSSS1*FTT23*FRRR3 - FSS13*FRRR2*FTTT3 + 
     .             FSS13*FRRR3*FTTT2 - FSS23*FTTT1*FRRR3 + 
     .             FTT13*FRRR2*FSSS3 - FTT13*FRRR3*FSSS2

          fdv12=   FRRR1*FSS12*FTT23 - FRRR1*FSS23*FTT12 - 
     .             FRRR1*FS123*FTTT2 + FRRR1*FT123*FSSS2 - 
     .             FRR12*FSSS1*FTT23 - FRR12*FSS13*FTTT2 + 
     .             FRR12*FSS23*FTTT1 + FRR12*FTT13*FSSS2 + 
     .             FRR13*FSS12*FTTT2 - FRR13*FTT12*FSSS2 + 
     .             FRR23*FSSS1*FTT12 - FRR23*FSS12*FTTT1 + 
     .             FR123*FSSS1*FTTT2 - FR123*FTTT1*FSSS2 - 
     .             FSSS1*FT123*FRRR2 - FSS12*FTT13*FRRR2 + 
     .             FSS13*FTT12*FRRR2 + FS123*FTTT1*FRRR2

          fdv13= - FRRR1*FSS13*FTT23 + FRRR1*FSS23*FTT13 + 
     .             FRRR1*FS123*FTTT3 - FRRR1*FT123*FSSS3 - 
     .             FRR12*FSS13*FTTT3 + FRR12*FTT13*FSSS3 + 
     .             FRR13*FSSS1*FTT23 + FRR13*FSS12*FTTT3 - 
     .             FRR13*FSS23*FTTT1 - FRR13*FTT12*FSSS3 - 
     .             FRR23*FSSS1*FTT13 + FRR23*FSS13*FTTT1 - 
     .             FR123*FSSS1*FTTT3 + FR123*FTTT1*FSSS3 +
     .             FSSS1*FT123*FRRR3 - FSS12*FTT13*FRRR3 + 
     .             FSS13*FTT12*FRRR3 - FS123*FTTT1*FRRR3

          fdv23=   FRR12*FSS23*FTTT3 - FRR12*FTT23*FSSS3 - 
     .             FRR13*FSS23*FTTT2 + FRR13*FTT23*FSSS2 - 
     .             FRR23*FSS12*FTTT3 + FRR23*FSS13*FTTT2 + 
     .             FRR23*FTT12*FSSS3 - FRR23*FTT13*FSSS2 +
     .             FR123*FSSS2*FTTT3 - FR123*FSSS3*FTTT2 + 
     .             FSS12*FTT23*FRRR3 - FSS13*FTT23*FRRR2 - 
     .             FSS23*FTT12*FRRR3 + FSS23*FTT13*FRRR2 - 
     .             FS123*FRRR2*FTTT3 + FS123*FRRR3*FTTT2 + 
     .             FT123*FRRR2*FSSS3 - FT123*FRRR3*FSSS2

          fdv11=   FRRR1*FSS12*FTT13 - FRRR1*FSS13*FTT12 - 
     .             FRR12*FSSS1*FTT13 + FRR12*FSS13*FTTT1 + 
     .             FRR13*FSSS1*FTT12 - FRR13*FSS12*FTTT1

          fdv22= - FRR12*FSS23*FTTT2 + FRR12*FTT23*FSSS2 + 
     .             FRR23*FSS12*FTTT2 - FRR23*FTT12*FSSS2 - 
     .             FSS12*FTT23*FRRR2 + FSS23*FTT12*FRRR2

          fdv33=   FRR13*FSS23*FTTT3 - FRR13*FTT23*FSSS3 - 
     .             FRR23*FSS13*FTTT3 + FRR23*FTT13*FSSS3 + 
     .             FSS13*FTT23*FRRR3 - FSS23*FTT13*FRRR3

          fd123=   2.0D0*(
     .           - FRR12*FSS13*FTT23 + FRR12*FSS23*FTT13 + 
     .             FRR13*FSS12*FTT23 - FRR13*FSS23*FTT12 - 
     .             FRR23*FSS12*FTT13 + FRR23*FSS13*FTT12
     .             )

          fd112=   FRRR1*FSS12*FT123 - FRRR1*FS123*FTT12 - 
     .             FRR12*FSSS1*FT123 + FRR12*FS123*FTTT1 + 
     .             FR123*FSSS1*FTT12 - FR123*FSS12*FTTT1

          fd113= - FRRR1*FSS13*FT123 + FRRR1*FS123*FTT13 + 
     .             FRR13*FSSS1*FT123 - FRR13*FS123*FTTT1 - 
     .             FR123*FSSS1*FTT13 + FR123*FSS13*FTTT1

          fd122= - FRR12*FS123*FTTT2 + FRR12*FT123*FSSS2 + 
     .             FR123*FSS12*FTTT2 - FR123*FTT12*FSSS2 -
     .             FSS12*FT123*FRRR2 + FS123*FTT12*FRRR2

          fd133=   FRR13*FS123*FTTT3 - FRR13*FT123*FSSS3 - 
     .             FR123*FSS13*FTTT3 + FR123*FTT13*FSSS3 + 
     .             FSS13*FT123*FRRR3 - FS123*FTT13*FRRR3

          fd223=   FRR23*FS123*FTTT2 - FRR23*FT123*FSSS2 - 
     .             FR123*FSS23*FTTT2 + FR123*FTT23*FSSS2 + 
     .             FSS23*FT123*FRRR2 - FS123*FTT23*FRRR2

          fd233= - FRR23*FS123*FTTT3 + FRR23*FT123*FSSS3 + 
     .             FR123*FSS23*FTTT3 - FR123*FTT23*FSSS3 - 
     .             FSS23*FT123*FRRR3 + FS123*FTT23*FRRR3

          f1123= - FRR12*FSS13*FT123 + FRR12*FS123*FTT13 + 
     .             FRR13*FSS12*FT123 - FRR13*FS123*FTT12 - 
     .             FR123*FSS12*FTT13 + FR123*FSS13*FTT12

          f1223=   FRR12*FSS23*FT123 - FRR12*FS123*FTT23 - 
     .             FRR23*FSS12*FT123 + FRR23*FS123*FTT12 + 
     .             FR123*FSS12*FTT23 - FR123*FSS23*FTT12

          f1233= - FRR13*FSS23*FT123 + FRR13*FS123*FTT23 + 
     .             FRR23*FSS13*FT123 - FRR23*FS123*FTT13 - 
     .             FR123*FSS13*FTT23 + FR123*FSS23*FTT13
C
          CALL PHCC08T(CARTDT,DVOLUT,EHISTT,ELCODT,ELDIST,EPMTXT,
     .                 GPCODT,LNODST,PROPST,RMAT1T,SHAPET,STRANT,
     .                 STRSGT,STRA0T,STRS0T,
     .                 BMSIGT,BMATXT,DESIGT,DMATXT,DSTRAT,PRESGT,
     .                 SGTOTT,SIGMAT,TSTRAT,XJACMT,ELELTT,VELCMT,
     .                 SHAPDT,DVOLDT,
     .                 POSGPT,WEIGPT,
     .                 DERIDT,
     .                 CARDDT,GPCDDT,
     .                 ELEL1T,ELEL2T,
     .                  frrrc, frrr1, frrr2, frrr3, frr12, frr13, frr23,
     .                  fr123,
     .                  fsssc, fsss1, fsss2, fsss3, fss12, fss13, fss23,
     .                  fs123,
     .                  ftttc, fttt1, fttt2, fttt3, ftt12, ftt13, ftt23,
     .                  ft123,
     .                  fdvoc, fdvo1, fdvo2, fdvo3, fdv12, fdv13, fdv23,
     .                  fdv11, fdv22, fdv33, fd123, fd112, fd113, fd122,
     .                  fd133, fd223, fd233, f1123, f1223, f1233,
     .                 HACHET,WHADEL,CENTRO,ADVEMT,TEINIT,FPCHLT)
         enddo
        enddo
       enddo
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.10) THEN
       CALL RUNENDT('FRMC3DT:ERROR DETECTED WITH NDIME=3')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.20) THEN
       CALL RUNENDT('FRMC3DT:ERROR DETECTED WITH NDIME=3')
       RETURN
      ENDIF
C
      IF(NNODLT.EQ.27) THEN
       CALL RUNENDT('FRMC3DT:ERROR DETECTED WITH NDIME=3')
       RETURN
      ENDIF
C
      END
