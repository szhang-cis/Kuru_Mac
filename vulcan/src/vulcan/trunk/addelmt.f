      SUBROUTINE ADDELMT
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE POINTERS FOR ELEMENTAL ARRAYS
C
C     Notes: VNORI is the normal vector stored at element level. This
C            array is only useful when NMEMO2=0. 
C
C***********************************************************************
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_omt.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_omt.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'         ! thermal-mechanical
C
C**** THERMAL VARIABLES
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      ITERME1=1                              ! unidirectional coupled
      ITERME2=1                              ! bidirectional coupled
      IF(ITERME.LT.0) THEN                   ! uncoupled problems
       ITERME1=0                             ! not used
       ITERME2=0
      ENDIF
      IF(ITERME.EQ.0) THEN                   ! unidirectional coupled
       ITERME1=1                             ! not used
       ITERME2=0
      ENDIF
C
C**** GENERAL CHECKS OF VARIABLES DEFINED IN setdatt.f, setdats.f &
C     consett.f
C
      IF(MPROPT.LT.NPROPT)
     . CALL RUNENDT('ERROR: MPROPT LT NPROPT')
C
      IF(MHLODT.LT.NHLODT)
     . CALL RUNENDT('ERROR: MHLODT LT NHLODT')
C
      IF(MSUBFT.LT.NSUBFT)
     . CALL RUNENDT('ERROR: MSUBFT LT NSUBFT')
C
      IF(MPRELT.LT.NPRELT)
     . CALL RUNENDT('ERROR: MPRELT LT NPRELT')
C
      IF(MHISTT.LT.NHISTT)
     . CALL RUNENDT('ERROR: MHISTT LT NHISTT')
C
      IF(MNUINT.LT.NNUINT)
     . CALL RUNENDT('ERROR: MNUINT LT NNUINT')
C
C**** GENERAL CHECKS OF VARIABLES DEFINED IN INPUT DATA
C     (now in setdatt.f)
C
      IF(MDISKD.NE.NDISKD)
     . CALL RUNENDT('ERROR: MDISKD =/ NDISKD')
C
      IF(MFURES.LT.NFURES)               ! this control is not necessary
     . CALL RUNENDT('ERROR: MFURES =/ NFURES')
C
      IF(MMACHI.NE.NMACHI)
     . CALL RUNENDT('ERROR: MMACHI =/ NMACHI')
C
      IF(MMEMO.NE.NMEMO)
     . CALL RUNENDT('ERROR: MMEMO =/ NMEMO')
C
      IF(MMEMO1.NE.NMEMO1)
     . CALL RUNENDT('ERROR: MMEMO1 =/ NMEMO1')
C
      IF(MMEMO2.NE.NMEMO2)
     . CALL RUNENDT('ERROR: MMEMO2 =/ NMEMO2')
C
      IF(MMEMO3.LT.NMEMO3)
     . CALL RUNENDT('ERROR: MMEMO3 LT NMEMO3')
C
      IF(MMEMO4.LT.NMEMO4)
     . CALL RUNENDT('ERROR: MMEMO4 =/ NMEMO4')
C
      IF(MMEMO5.NE.NMEMO5)
     . CALL RUNENDT('ERROR: MMEMO5 =/ NMEMO5')
C
      IF(MMEMO6.NE.NMEMO6)
     . CALL RUNENDT('ERROR: MMEMO6 =/ NMEMO6')
C
      IF(MMEMO7.NE.NMEMO7)
     . CALL RUNENDT('ERROR: MMEMO7 =/ NMEMO7')
C
      IF(MMEMO8.NE.NMEMO8)
     . CALL RUNENDT('ERROR: MMEMO8 =/ NMEMO8')
C
      IF(MMEMO9.NE.NMEMO9)
     . CALL RUNENDT('ERROR: MMEMO9 =/ NMEMO9')
C
      IF(MMEMO10.LT.NMEMO10)
     . CALL RUNENDT('ERROR: MMEMO10 LT NMEMO10')
C
      IF(MMEMO11.NE.NMEMO11)
     . CALL RUNENDT('ERROR: MMEMO11 =/ NMEMO11')
C
C**** INITIAL DATA
C
      IPOONT=0
      IDYONT=0
      ISMOMT=0
C
C**** ELEMENT DATA ARRAY ( ELDATT )
C
      NKOSTT=(NSTR1T+1)*NSTR1T/2
      IF(KSGAUT.NE.0) ISMOMT=1
C
      IF (KPROBT.EQ.2) THEN       ! STRUCTURAL SHELLS AND PLATES
       IDATAT(1)=1                                  ! ELCOD(NDIME,NNODE)
       IDATAT(2)=IDATAT(1)+NDIMET*NNODET          ! DERIV(2,NNODE,NGAUS)
       IDATAT(3)=IDATAT(2)+2*NNODET*NGAUST          ! DVOLU(NGAUS)
       IDATAT(4)=IDATAT(3)+NGAUST                   ! GPCOD(NDIME,NGAUS)
       IDATAT(5)=IDATAT(4)+NDIMET*NGAUST            ! SHAPE(NNODE,NGAUS)
       IDATAT(6)=IDATAT(5)+NNODET*NGAUST            ! EPMTX(NKOST,NGAUS)
       IDATAT(7)=IDATAT(6)+NKOSTT*NGAUST            ! TRNOD(NDIME,NNODE)
       IDATAT(8)=IDATAT(7)+NDIMET*NNODET*(NDIMET-2) ! EMASS(NNODE,NNODE)
       NDATAT   =IDATAT(8)+NNODET*NNODET*ISMOMT
      ELSE                        ! OTHERS PROBLEMS
       NMEMA1=1-NMEMO1
       NMEMA2=1-NMEMO2
       IDATAA=1
       IF(NMEMA1.EQ.0.AND.NMEMA2.EQ.0.AND.NMEMO.EQ.0.AND.
     .    (NMEMA2*ISMOMT).EQ.0.AND.NHOURT.EQ.0) IDATAA=0+IWINDT
       IDATAT(1)=IDATAA                             ! ELCOD(NDIME,NNODE)
       IDATAT(2)=IDATAT(1)+
     .           NDIMET*NNODET*NMEMA1         ! CARTD(NDIME,NNODE,NGAUS)
       IDATAT(3)=IDATAT(2)+NDIMET*NNODET*NGAUST*NMEMA2    ! DVOLU(NGAUS)
       IDATAT(4)=IDATAT(3)+NGAUST*NMEMA2            ! GPCOD(NDIME,NGAUS)
       IDATAT(5)=IDATAT(4)+NDIMET*NGAUST*NMEMA2     ! SHAPE(NNODE,NGAUS)
       IDATAT(6)=IDATAT(5)+NNODET*NGAUST*NMEMA2     ! EPMTX(NKOST,NGAUS)
       IDATAT(7)=IDATAT(6)+NKOSTT*NGAUST*NMEMO      ! RMAT1(NDIME,NDIME)
       IDATAT(8)=IDATAT(7)+NDIMET*NDIMET*NMEMO      ! EMASS(NNODE,NNODE)
       IDATAT(9)=IDATAT(8)+NNODET*NNODET*ISMOMT*NMEMA2
C                                                   ! STIFH(NEVAB,NEVAB)
       IDATAT(10)=IDATAT( 9)+NEVABT*NEVABT*NHOURT   ! VNORI(NDIME,NGAUS)
       NDATAT    =IDATAT(10)+NDIMET*NGAUST*NMEMA2*ITERME2
      ENDIF
C
C**** ELEMENT PRESCRIBED VARIABLES ARRAY ( ELPRE )
C
      IPREVTA=0+IWINDT
      IF(NMEMO.EQ.1) IPREVTA=1
      IPREVT(1)=IPREVTA                             ! STRA0(NSTR1,NGAUS)
      IPREVT(2)=IPREVT(1)+NSTR1T*NGAUST*NMEMO       ! STRS0(NSTR1,NGAUS)
      IPREVT(3)=IPREVT(2)+NSTR1T*NGAUST*NMEMO       ! TEMPC(4)
      NPREVT   =IPREVT(3)+4*NMEMO
C
C**** ELEMENT STATE VARIABLES ARRAY ( ELVAR )
C
      NMEMA5=1-NMEMO5
      ISTATA=0+IWINDT
      IF(NMEMO3.EQ.1.OR.NMEMO4.EQ.1.OR.NMEMO5.EQ.0)
     . ISTATA=1
      ISTATT(1)=ISTATA                              ! ELDIS(NDOFC,NNODE)
      ISTATT(2)=ISTATT(1)+NDOFCT*NNODET*NMEMA5      ! EHIST(NHIST,NGAUS)
      ISTATT(3)=ISTATT(2)+NHISTT*NGAUST*NMEMO3      ! STRAN(NSTR1,NGAUS)
      ISTATT(4)=ISTATT(3)+NSTR1T*NGAUST*NMEMO4      ! STRSG(NSTR1,NGAUS)
      NSTATT   =ISTATT(4)+NSTR1T*NGAUST*NMEMO4
C                           
C**** ELEMENT MATRICES ARRAY ( ELMAT )
C
      IF(KDYNAT.EQ.1) IDYONT=1
      IF(KPORET.EQ.2) IPOONT=1
C
      NMEMA6=1-NMEMO6
      NMEMA7=1-NMEMO7
      IMATXTA=1
      IF(NMEMA7.EQ.0.AND.IPOONT.EQ.0) IMATXTA=0+IWINDT
      IMATXT(1)=IMATXTA                             ! CSTIF(NEVAC,NEVAC)
      IMATXT(2)=IMATXT(1)+NEVACT*NEVACT*NMEMA7      ! KSTIF(NKOVA)
      IMATXT(3)=IMATXT(2)+NKOVAT*NMEMA6             ! MSTIF(NKOVA)
      IMATXT(4)=IMATXT(3)+NKOVAT*IDYONT*NMEMA6      ! PSTIF(NKOND)
      IMATXT(5)=IMATXT(4)+NKONDT*IPOONT             ! QSTIF(NKOND)
      IMATXT(6)=IMATXT(5)+NKONDT*IPOONT             ! HSTIF(NEVAB*NNODE)
      NMATXT   =IMATXT(6)+NEVABT*NNODET*IPOONT
C
C**** MEMORY CONTROL
C
      IF(NDATAT.GT.MDATAT)
     . CALL RUNENDT('ERROR IN ADDELMT: NDATAT          ')
      IF(NPREVT.GT.MPREVT)
     . CALL RUNENDT('ERROR IN ADDELMT: NPREVT          ')
      IF(NSTATT.GT.MSTATT)
     . CALL RUNENDT('ERROR IN ADDELMT: NSTATT          ')
      IF(NMATXT.GT.MMATXT)
     . CALL RUNENDT('ERROR IN ADDELMT: NMATXT          ')
C
      RETURN
      END
