C***********************************************************************
C
C**** THERMAL PARAMETERS
C
C     !m mechanical variable
C
C***********************************************************************
C
      INCLUDE 'gend_omt.f'
C
      PARAMETER(
c    . MNODGT=max(MNODET,MGAUST) )                ! sg
     . MNODGT=(MNODET+MGAUST) )                   ! PC & linux
C
C**** ADDELMT
C
      PARAMETER(
c    . MDATA1A=max(MMEMA1,MMEMA2,MMEMO,(MSMOMT*MMEMA2),MHOURT) )  ! sg
     . MDATA1A=0+MWINDT )                           ! PC & linux
C
      PARAMETER(
     . MDATA1T=MDATA1A,                             ! ELCOD(NDIME,NNODE)
     . MDATA2T=MDATA1T+MDIMET*MNODET*MMEMA1,  ! CARTD(NDIME,NNODE,NGAUS)
     . MDATA3T=MDATA2T+MDIMET*MNODET*MGAUST*MMEMA2,       ! DVOLU(NGAUS)
     . MDATA4T=MDATA3T+MGAUST*MMEMA2,               ! GPCOD(NDIME,NGAUS)
     . MDATA5T=MDATA4T+MDIMET*MGAUST*MMEMA2,        ! SHAPE(NNODE,NGAUS)
     . MDATA6T=MDATA5T+MNODET*MGAUST*MMEMA2,        ! EPMTX(NKOST,NGAUS)
     . MDATA7T=MDATA6T+MKOSTT*MGAUST*MMEMO,         ! RMAT1(NDIME,NDIME)
     . MDATA8T=MDATA7T+MDIMET*MDIMET*MMEMO,         ! EMASS(NNODE,NNODE)
     . MDATA9T=MDATA8T+MNODET*MNODET*MSMOMT*MMEMA2, ! STIFH(NEVAB,NEVAB)
     . MDATA10T=MDATA9T+MEVABT*MEVABT*MHOURT,       ! VNORI(NDIME,NGAUS)
     . MDATAT  =MDATA10T+MDIMET*MGAUST*MMEMA2 )
C
      PARAMETER(
     . MPREV1T=1*MMEMO+MWINDT,                      ! STRA0(NSTR1,NGAUS)
     . MPREV2T=MPREV1T+MSTR1T*MGAUST*MMEMO,         ! STRS0(NSTR1,NGAUS)
     . MPREV3T=MPREV2T+MSTR1T*MGAUST*MMEMO,         ! TEMPC(4)
     . MPREVT =MPREV3T+4*MMEMO )
C
      PARAMETER(
c    . MSTAT1A=max(MMEMO3,MMEMO4,MMEMA5) )   ! sg
     . MSTAT1A=1+MWINDT )                    ! PC & linux
C
      PARAMETER(
     . MSTAT1T=MSTAT1A,                             ! ELDIS(NDOFC,NNODE)
     . MSTAT2T=MSTAT1T+MDOFCT*MNODET*MMEMA5,        ! EHIST(NHIST,NGAUS)
     . MSTAT3T=MSTAT2T+MHISTT*MGAUST*MMEMO3,        ! STRAN(NSTR1,NGAUS)
     . MSTAT4T=MSTAT3T+MSTR1T*MGAUST*MMEMO4,        ! STRSG(NSTR1,NGAUS)
     . MSTATT =MSTAT4T+MSTR1T*MGAUST*MMEMO4 )
C
      PARAMETER(
c    . MPORETA=max(MPORET-1,0) )             ! sg
     . MPORETA=0 )                           ! PC & linux
C
      PARAMETER(
c    . MMATXTA=max(MMEMA7,MPORETA) )         ! sg
     . MMATXTA=1+MWINDT )                    ! PC & linux
C
      PARAMETER(
     . MMATX1T=MMATXTA,                             ! CSTIF(NEVAC,NEVAC)
     . MMATX2T=MMATX1T+MEVACT*MEVACT*MMEMA7,        ! KSTIF(NKOVA)
     . MMATX3T=MMATX2T+MKOVAT*MMEMA6,               ! MSTIF(NKOVA)
     . MMATX4T=MMATX3T+MKOVAT*MDYNAT*MMEMA6,        ! PSTIF(NKOND)
     . MMATX5T=MMATX4T+MKONDT*MPORETA,              ! QSTIF(NKOND)
     . MMATX6T=MMATX5T+MKONDT*MPORETA,              ! HSTIF(NEVAB*NNODE)
     . MMATXT =MMATX6T+MEVABT*MNODET*MPORETA )
C
C**** ADDPRIT
C
      PARAMETER(
c    . MTER12=max((MTERME+1),0) )            ! sg
     . MTER12=2 )                            ! PC & linux
      PARAMETER(
c    . MTER2=max(MTERME,0) )                 ! sg
     . MTER2=1 )                             ! PC & linux
      PARAMETER(
     . MTER1=MTER12/(1+MTER2) )
C
      PARAMETER(              ! simplification
     . MTER3=MTER2,
c    . MTER113=max(MTER3,MMEMO11) )          ! sg
     . MTER113=1 )                           ! PC & linux
C
      PARAMETER(              ! simplification
c    . MTER2X=max(MTER113,MOCOIT) )          ! sg
     . MTER2X=1 )                            ! PC & linux
C
      PARAMETER(
     . MPREAS=MTER2*MMEMA11, 
     . MPREAT=MPREAS+MPOROT )
C
      PARAMETER(                             ! simplified
c    . MITERCT=max(MTERME,0) )               ! sg
     . MITERCT=1 )                           ! PC & linux
C
      PARAMETER(
c    . MPORETB=max(MPORET,1) )               ! sg
     . MPORETB=1 )                           ! PC & linux
      PARAMETER(
     . MPORETC=MPORET/MPORETB )
C
      PARAMETER(
c    . MEGX=max(MCONVT*MMEMO3,0) )           ! sg
     . MEGX=1 )                              ! PC & linux
C
      PARAMETER(
c    . MEGX1=max(MCONVT*MMEMO3,0) )          ! sg
     . MEGX1=1 )                             ! PC & linux
C
      PARAMETER(
     . mprin1t=1,                             ! LNODS(NNODE,NELEM)
     . mprin2t =mprin1t +(MNODET*MELEMT*MCHALT+4)/8, ! MATNO(NELEM)
     . mprin3t =mprin2t +(MELEMT*MCHALT+4)/8, ! PROEL(NPREL,NGRUP)
     . mprin4t =mprin3t +MPRELT*MGRUPT,       ! PROPS(NPROP,NMATS)
     . mprin5t =mprin4t +MPROPT*MMATST,       ! COORD(NDIME,NPOIN)#
     . mprin6t =mprin5t +MDIMET*MPOINT*MMEMO1,! HTLOD(NHLOD,NSUBF,NFUNC)
     . mprin7t =mprin6t +MHLODT*MSUBFT*MFUNCT,! IFFIX(NTOTV,2)
     . mprin8t =mprin7t +(MTOTVT*MCHALT+4)/8*
     .                      (1+MPORETB/2)*(1+MACTIT),   ! PRESC(NTOTV,2)
     . mprin9t =mprin8t +MTOTVT*2,            ! RLOAD(NTOTV)
     . mprin10t=mprin9t +MTOTVT,              ! RLOAH(NTOTV,NFUNC)
     . mprin11t=mprin10t+MTOTVT*MFUNCT,       ! FICTO(NFUNC)
     . mprin12t=mprin11t+MFUNCT )             ! TFICT(NFUNC)
      PARAMETER(
     . mprin13t=mprin12t+MFUNCT,              ! DISIT(NTOTV,2)
     . mprin14t=mprin13t+MTOTVT*(1+MMEMO9),   ! DISPR(NTOTV,3)
     . mprin15t=mprin14t+MTOTVT*(1+MDYNAT+MMEMO8),      ! DISTO(NTOTV,3)
     . mprin16t=mprin15t+
     .          MTOTVT*(1+MDYNAT+MMEMO8+MEGX1),         ! HEADS(NPOIN,4)
     . mprin17t=mprin16t+MPOINT*4*MPORETC,    ! REFOR(NTOTV,2)
     . mprin18t=mprin17t+MTOTVT*2,            ! TLOAD(NTOTV,2)
     . mprin19t=mprin18t+MTOTVT*2,            ! LPNTN(NPOIN)
     . mprin20t=mprin19t+(MPOINT*MCHALT+4)/8*MRENUT,      ! ELDAT(NDATA)
     . mprin21t=mprin20t+MDATAT,              ! ELPRE(NPREV)
     . mprin22t=mprin21t+MPREVT,              ! ELVAR(NSTAT)
     . mprin23t=mprin22t+MSTATT,              ! ELMAT(NMATX)
     . mprin24t=mprin23t+MMATXT )             ! DISPL(NTOTV) !m
      PARAMETER(
     . mprin25t=mprin24t+MTOTVM*MTER2,                  ! PWORK(NPOIN,3)
     . mprin26t=mprin25t+MPOINT*(2+MITERCT)*MTER2,      ! PREAS(NPREA)
     . mprin27t=mprin26t+MPOINT*MPREAT,                 ! TGAPS(NPOIN)
     . mprin28t=mprin27t+MPOINT*MTER2*MMEMA11,          ! TEMPI(NPOIN,2)
     . mprin29t=mprin28t+MPOINT*(1+MMEMO10),        ! ADVEL(NTOTV) !m
     . mprin30t=mprin29t+MTOTVM*MCONVT,             ! FPCHA(NFPCH,NPOIN)
     . mprin31t=mprin30t+MFPCH*MPOINT,              ! LACTI(NELEM)
     . MPRINT  =mprin31t+MACTIT*MELEMT )
C
C**** RSSETPT
C
      PARAMETER(
     .       i11t =MDATAT,                    ! ELDAT
     .       i21t =MPREVT,                    ! ELPRE
     .       i31t =MSTATT,                    ! ELVAR
     .       i41t =MPREVT,                    ! ELPRE
     .       i51t =MSTATT,                    ! ELVAR
     .       i61t =MEVACT*MEVACT*MMEMA7,      ! CSTIF
     .       i71t =MKOVAT*MMEMA6,             ! ESTIF
     .       i81t =MKOVAT*MDYNAT*MMEMA6,      ! WSTIF
     .       i91t =MKONDT*MPORETA,            ! PSTIF
     .       i101t=MKONDT*MPORETA,            ! QSTIF
     .       i111t=MEVABT*MNODET*MPORETA,     ! HSTIF
     .       i121t=MTOTVT)                    ! DISTO
C
      PARAMETER(
     .       i23t =MDATAT,
     .       i33t =(i23t+MPREVT),
     .       i43t =(i33t+MSTATT),
     .       i53t =(i43t+MPREVT),
     .       i63t =(i53t+MSTATT),
     .       i73t =(i63t+MEVACT*MEVACT*MMEMA7),
     .       i83t =(i73t+MDYNAT*i71t),
     .       i93t =(i83t+i81t),
     .       i103t=(i93t+i91t),
     .       i113t=(i103t+i101t),
     .       i123t=(i113t+i111t) )
C
C**** DATBAST (BUFFERT)
C
C     DISTOT is always IN-CORE for MBUFFERT
C
      PARAMETER(
c    . MWORPT =max(i23t,i33t,i43t,i53t,i63t,i73t,            ! sg
c    .             i83t,i93t,i103t,i113t,i123t),
     . MWORPT =i123t,                                        !PC & linux
c    . index1t=max(MDATAT,MPREVT,MSTATT,(MEVACT*MEVACT)*MMEMA7,    ! sg
c    .             MKOVAT*MMEMA6,MKONDT*MPORETA,(MEVABT*MNODET),MTOTVT),
     . index1t=(MDATAT+MPREVT+MSTATT+(MEVACT*MEVACT)*MMEMA7+ !PC & linux
     .          MKOVAT*MMEMA6+MKONDT*MPORETA+(MEVABT*MNODET)+MTOTVT),
     . index2t=   (MDATAT+MPREVT+MSTATT+MEVACT*MEVACT*MMEMA7+
     .             MKOVAT*MMEMA6+MKONDT*MPORETA+MEVABT*MNODET+MTOTVT),
     . MBUFFERT=MTOTVT+(MELEMT*MWORPT+INDEX2T+1+INDEX1T)*MDISKD )
C
C**** ADDWORT
C
      PARAMETER(
     . ISTA1T =1,                                   ! COORD(NDIME,NPOIN)
     . ISTA2T =ISTA1T +MDIMET*MPOINT*MMEMA1,        ! GPCOD(NDIME,NGAUS)
     . ISTA3T =ISTA2T +MDIMET*MGAUST,         ! CARTD(NDIME,NNODE,NGAUS)
     . ISTA4T =ISTA3T +MDIMET*MNODET*MGAUST*MMEMO2, ! DVOLU(NGAUS)
     . ISTA5T =ISTA4T +MGAUST*MMEMO2,               ! GPCOD(NDIME,NGAUS)
     . ISTA6T =ISTA5T +MDIMET*MGAUST*MMEMO2,        ! SHAPE(NNODE,NGAUS)
     . ISTA7T =ISTA6T +MNODET*MGAUST*MMEMO2,  ! DERIV(NDIME,NNODE,NGAUS)
     . ISTA8T =ISTA7T +MDIMET*MNODET*MGAUST*MMEMO2, ! POSGP(NDIME,NGAUS)
     . ISTA9T =ISTA8T +MDIMET*MGAUST*MMEMO2,        ! WEIGP(NGAUS)
     . ISTA10T=ISTA9T +MGAUST*MMEMO2,               ! XJACM(NDIME,NDIME)
     . ISTA11T=ISTA10T+MDIMET*MDIMET*MMEMO2,        ! ELCO1(NDIME,NNODE)
     . ISTA12T=ISTA11T+MDIMET*MNODET*MMEMO2,        ! EMASS(NNODE,NNODE)
     . ISTA13T=ISTA12T+MNODET*MNODET*MMEMO2,        ! ELDI1(NDOFN,NNODE)
     . ISTA14T=ISTA13T+MDOFCT*MNODET*MMEMO2,        ! VNORL(NDIME,NGAUS)
     . ISTA15T=ISTA14T+MDIMET*MGAUST*MTER2*MMEMO11, ! DVOLI(NGAUS)
     . ISTA16T=ISTA15T+MGAUST*MTER3,                ! ELCOI(NDIME,NNODE)
     . ISTA17T=ISTA16T+MDIMET*MNODET*MTER3,      ! DISPT(NDOFC,NNODE) !m
     . LINDT  =ISTA17T+MDOFCM*MNODET*MTER3 )
C
      PARAMETER(
     . LRENT=MRENUT*MPOINT*(10*MNODET+3)/2 )
C
      PARAMETER(
     . IQUA1T=1,                              ! VVECT(NTOTV)
     . IQUA2T=IQUA1T+MTOTVT,                  ! WVECT(NTOTV)
     . LQUAT =IQUA2T+MTOTVT )
C
      PARAMETER(
     . IFIX1T=1,                              ! LEQNS(NEVAC)
     . IFIX2T=IFIX1T+(MEVACT*MCHALT+4)/8,     ! LNUEQ(NTOTV)
     . IFIX3T=IFIX2T+(MTOTVT*MCHALT+4)/8,     ! LOCEL(NEVAC)
     . IFIX4T=IFIX3T+(MEVACT*MCHALT+4)/8,     ! LPONT(NTOTV)
     . IFIX5T=IFIX4T+(MTOTVT*MCHALT+4)/8,     ! NACVA(NTOTV)
     . IFIX6T=IFIX5T+(MTOTVT*MCHALT+4)/8,     ! NDEST(NEVAC)
     . IFIX7T=IFIX6T+(MEVACT*MCHALT+4)/8,     ! NDFRO(NELEM)
     . IFIX8T=IFIX7T+(MELEMT*MCHALT+4)/8,     ! PRESC(NDOFC)
     . LFIXT =IFIX8T+MDOFCT )
C
      PARAMETER(
     . IFOT1T=1,                              ! BMSIG(NEVAB)
     . IFOT2T=IFOT1T+MEVABT,                  ! BMATX(NSTR1,NEVAB)
     . IFOT3T=IFOT2T+MSTR1T*MEVABT,           ! DESIG(NSTR1) 
     . IFOT4T=IFOT3T+MSTR1T,                  ! DMATX(NSTR1,NSTR1)
     . IFOT5T=IFOT4T+MSTR1T*MSTR1T,           ! DSTRA(NSTR1)
     . IFOT6T=IFOT5T+MSTR1T,                  ! PRESG(NSTR1)
     . IFOT7T=IFOT6T+MSTR1T,                  ! SGTOT(NSTR1)
     . IFOT8T=IFOT7T+MSTR1T,                  ! SIGMA(NSTR1)
     . IFOT9T=IFOT8T+MSTR1T,                  ! TSTRA(NSTR1)
     . IFOT10T=IFOT9T+MSTR1T,                 ! XJACM(NDIME,NDIME)
     . IFOT11T=IFOT10T+MDIMET*MDIMET,         ! ELELT(NEVAB)
     . IFOT12T=IFOT11T+MEVABT,                ! VELCM(NEVAB)
     . IFOT13T=IFOT12T+MEVABT,                ! DISPT(NDOFC,NNODET) !m
     . IFOT14T=IFOT13T+MDOFCM*MNODET*MTER2*MMEMO11, ! PREASL(NNODET)
     . IFOT15T=IFOT14T+MNODET*MTER2*MMEMA11,  ! TGAPSL(NNODET)
     . IFOT16T=IFOT15T+MNODET*MTER2*MMEMA11 ) ! ELCO1(NDIME,NNODE)
      PARAMETER(
c    . MAUX1=max(MCONVT,MGALET) )             ! sg
     . MAUX1=1 )                              ! PC & linux
      PARAMETER(
     . IFOT17T=IFOT16T+MDIMET*MNODET,         ! SHAPD(NNODE,NGAUS)
     . IFOT18T=IFOT17T+MNODET*MGAUST*MAUX1,   ! DVOLD(NGAUS)
     . IFOT19T=IFOT18T+MGAUST*MCONVT,         ! POSGP(NDIME,NGAUS)
     . IFOT20T=IFOT19T+MDIMET*MGAUST*MAUX1,   ! WEIGP(NGAUS)
     . IFOT21T=IFOT20T+MGAUST*MAUX1,          ! DERID(NDIME,NNODE,NGAUS)
     . IFOT22T=IFOT21T+MDIMET*MNODET*MGAUST*MAUX1,  ! XJACM(NDIME,NDIME)
     . IFOT23T=IFOT22T+MDIMET*MDIMET*MAUX1,   ! CARDD(NDIME,NNODE,NGAUS)
     . IFOT24T=IFOT23T+MDIMET*MNODET*MGAUST*MAUX1,  ! GPCDD(NDIME,NGAUS)
     . IFOT25T=IFOT24T+MDIMET*MGAUST*MAUX1,   ! ELEL1(NEVAB)
     . IFOT26T=IFOT25T+MEVABT*MCONVT,         ! ELEL2(NEVAB)
     . IFOT27T=IFOT26T+MEVABT*MCONVT,         ! WHAPE(NNODE,NGAUS)
     . IFOT28T=IFOT27T+MNODET*MGAUST,         ! HACHE(NNODE)
     . IFOT29T=IFOT28T+MNODET*MGALET,         ! WARTD(NDIME,NNODE,NGAUS)
     . IFOT30T=IFOT29T+MDIMET*MNODET*MGAUST,  ! WERIV(NDIME,NNODE,NGAUS)
     . IFOT31T=IFOT30T+MDIMET*MNODET*MGAUST*MGALET, ! WHADE(NNODE,NGAUS)
     . IFOT32T=IFOT31T+MNODET*MGAUST*MCONVT ) ! CENTR(NDIME,NNODE)
      PARAMETER(
     . IFOT33T=IFOT32T+MDIMET*MNODET*MGALET,  ! VEL1M(NEVAB)
     . IFOT34T=IFOT33T+MEVABT*MMICR,          ! ADVEM(NDIME*NNODE)
     . IFOT35T=IFOT34T+MDIMET*MNODET*MCONVT,  ! CARTD(NDIME,NNODE,NGAUS)
     . IFOT36T=IFOT35T+MDIMET*MNODET*MGAUST,  ! DVOLU(NGAUS)
     . IFOT37T=IFOT36T+MGAUST,                ! GPCOD(NDIME,NGAUS)
     . IFOT38T=IFOT37T+MDIMET*MGAUST,         ! SHAPE(NNODE,NGAUS)
     . IFOT39T=IFOT38T+MNODET*MGAUST,         ! DERIV(NDIME,NNODE,NGAUS)
     . IFOT40T=IFOT39T+MDIMET*MNODET*MGAUST,  ! POSGP(NDIME,NGAUS)
     . IFOT41T=IFOT40T+MDIMET*MGAUST,         ! WEIGP(NGAUS)
     . IFOT42T=IFOT41T+MGAUST,                ! XJACM(NDIME,NDIME)
     . IFOT43T=IFOT42T+MDIMET*MDIMET,         ! EMASS(NNODE,NNODE)
     . IFOT44T=IFOT43T+MNODET*MNODET,         ! ELDI1(NDOFN,NNODE)
     . IFOT45T=IFOT44T+MDOFCT*MNODET,         ! TEINI(NDOFN,NNODE)
     . IFOT46T=IFOT45T+MDOFCT*MNODET,         ! BOUCH(NDOFN,NNODE)
     . IFOT47T=IFOT46T+MDOFCT*MNODET*MMEMO10, ! VNORL(NDIME,NGAUS)
     . IFOT48T=IFOT47T+MDIMET*MGAUST*MTER2*MMEMO11 ) !FPCHL(NFPCH,NNODE)
      PARAMETER(
     . IFOT49T=IFOT48T+MFPCH*MNODET,          ! DVOLI(NGAUS)
     . IFOT50T=IFOT49T+MGAUST*MTER3,          ! ELCOI(NDIME,NNODE)
     . IFOT51T=IFOT50T+MDIMET*MNODET*MTER3,   ! ACCPN(NPOIN)
     . IFOT52T=IFOT51T+MPOINT*MEGX,           ! SMSTP(2,NPOIN)
     . IFOT53T=IFOT52T+2*MPOINT*MEGX,         ! SFIPA(NNODE,2)
     . IFOT54T=IFOT53T+MNODET*2*MEGX,         ! SFIPA(2,NNODE)
     . IFOT55T=IFOT54T+2*MNODET*MEGX,         ! ELDIQ(NDOFN,NNODE)
     . LFOTT  =IFOT55T+MDOFCT*MNODET*MEGX1 )
C
      PARAMETER(
     . IFOR1T=1,                              ! ELELM(NEVAB)
     . IFOR2T=IFOR1T+MEVABT,                  ! VELCM(NEVAB)
     . IFOR3T=IFOR2T+MEVABT,                  ! DISIM(NEVAB)
     . IFOR4T=IFOR3T+MEVABT,                  ! ESTIF(NKOVA)
     . IFOR5T=IFOR4T+MKOVAT,                  ! DSTIF(NKOVA)
     . IFOR6T=IFOR5T+MKOVAT,                  ! ELELT(NEVAB)
     . IFOR7T=IFOR6T+MEVABT,                  ! SHAPD(NNODE,NGAUS)
     . IFOR8T=IFOR7T+MNODET*MGAUST,           ! DVOLD(NGAUS)
     . IFOR9T=IFOR8T+MGAUST,                  ! POSGP(NDIME,NGAUS)
     . IFOR10T=IFOR9T+MDIMET*MGAUST,          ! WEIGP(NGAUS)
     . IFOR11T=IFOR10T+MGAUST,                ! DERID(NDIME,NNODE,NGAUS)
     . IFOR12T=IFOR11T+MDIMET*MNODET*MGAUST,  ! XJACM(NDIME,NDIME)
     . IFOR13T=IFOR12T+MDIMET*MDIMET,         ! CARDD(NDIME,NNODE,NGAUS)
     . IFOR14T=IFOR13T+MDIMET*MNODET*MGAUST,  ! GPCDD(NDIME,NGAUS)
     . IFOR15T=IFOR14T+MDIMET*MGAUST,         ! ELEL1(NEVAB)
     . IFOR16T=IFOR15T+MEVABT )               ! ELEL2(NEVAB)
      PARAMETER(
     . IFOR17T=IFOR16T+MEVABT,                ! ELCO1(NDIME,NNODE)
     . IFOR18T=IFOR17T+MDIMET*MNODET,         ! WHAPE(NNODE,NGAUS)
     . IFOR19T=IFOR18T+MNODET*MGAUST,         ! HACHE(NNODE)
     . IFOR20T=IFOR19T+MNODET*MGALET,         ! WHADE(NNODE,NGAUS)
     . IFOR21T=IFOR20T+MNODET*MGAUST,         ! CENTR(NDIME,NNODE)
     . IFOR22T=IFOR21T+MDIMET*MNODET*MGALET,  ! ADVEM(NDIME,NNODE)
     . IFOR23T=IFOR22T+MDIMET*MNODET*MCONVT,  ! CARTD(NDIME,NNODE,NGAUS)
     . IFOR24T=IFOR23T+MDIMET*MNODET*MGAUST,  ! DVOLU(NGAUS)
     . IFOR25T=IFOR24T+MGAUST,                ! GPCOD(NDIME,NGAUS)
     . IFOR26T=IFOR25T+MDIMET*MGAUST,         ! SHAPE(NNODE,NGAUS)
     . IFOR27T=IFOR26T+MNODET*MGAUST,         ! DERIV(NDIME,NNODE,NGAUS)
     . IFOR28T=IFOR27T+MDIMET*MNODET*MGAUST,  ! POSGP(NDIME,NGAUS)
     . IFOR29T=IFOR28T+MDIMET*MGAUST,         ! WEIGP(NGAUS)
     . IFOR30T=IFOR29T+MGAUST,                ! XJACM(NDIME,NDIME)
     . IFOR31T=IFOR30T+MDIMET*MDIMET,         ! EMASS(NNODE,NNODE)
     . IFOR32T=IFOR31T+MNODET*MNODET )        ! ELDI1(NDOFN,NNODE)
      PARAMETER(
     . IFOR33T=IFOR32T+MDOFCT*MNODET,         ! TEINI(NDOFN,NNODE)
     . IFOR34T=IFOR33T+MDOFCT*MNODET,         ! FPCHL(NFPCH,NNODE)
     . IFOR35T=IFOR34T+MFPCH*MNODET,                ! DVOLI(NGAUS)
     . IFOR36T=IFOR35T+MGAUST*MTER3,                ! ELCOI(NDIME,NNODE)
     . IFOR37T=IFOR36T+MDIMET*MNODET*MTER3,      ! DISPT(NDOFC,NNODE) !m
     . LFORT  =(IFOR37T+MDOFCM*MNODET*MTER3)*MDYNAT )
C
      PARAMETER(
     . ILOS1T=1,                              ! ALOAD(NEVAB)
     . ILOS2T=ILOS1T+MEVABT,                  ! DERIV(NDIME,NNODE,NGAUS)
     . ILOS3T=ILOS2T+MDIMET*MNODET*MGAUST,    ! ELEDG(NDIME,NNODE)
     . ILOS4T=ILOS3T+MDIMET*MNODET,           ! GPCOD(NDIME)
     . ILOS5T=ILOS4T+MDIMET,                  ! GVECT(NDIME)
     . ILOS6T=ILOS5T+MDIMET,                  ! NOPRS(NNODE)
     . ILOS7T=ILOS6T+MNODET,                  ! POSGP(NDIME,NGAUS)
     . ILOS8T=ILOS7T+MDIMET*MGAUST,           ! PRESS(NNODE,NDIME)
     . ILOS9T=ILOS8T+MNODET*MDIMET,           ! PVECT(NDIME)
     . ILOS10T=ILOS9T+MDIMET,                 ! SHAPE(NNODE)
     . ILOS11T=ILOS10T+MNODET,                ! VAREA(NDIME)
     . ILOS12T=ILOS11T+MDIMET,                ! WEIGP(NGAUS)
     . ILOS13T=ILOS12T+MGAUST,                ! XJACM(NDIME*NDIME)
     . ILOS14T=ILOS13T+MDIMET*MDIMET,         ! CARTD(NDIME,NNODE,NGAUS)
     . ILOS15T=ILOS14T+MDIMET*MNODET*MGAUST,  ! PRESC(NDOFC)
     . ILOS16T=ILOS15T+MDOFCT )               ! VELCM(NEVAB)
      PARAMETER(
     . ILOS17T=ILOS16T+MEVABT,                ! ELCO1(NDIME,NNODE)
     . ILOS18T=ILOS17T+MDIMET*MNODET,         ! WHAPE(NNODE,NGAUS)
     . ILOS19T=ILOS18T+MNODET*MGAUST,         ! HACHE(NNODE)
     . ILOS20T=ILOS19T+MNODET*MGALET,         ! CENTR(NDIME,NNODE)
     . ILOS21T=ILOS20T+MDIMET*MNODET*MGALET,  ! CARTD(NDIME,NNODE,NGAUS)
     . ILOS22T=ILOS21T+MDIMET*MNODET*MGAUST,  ! DVOLU(NGAUS)
     . ILOS23T=ILOS22T+MGAUST,                ! GPCOD(NDIME,NGAUS)
     . ILOS24T=ILOS23T+MDIMET*MGAUST,         ! SHAPE(NNODE,NGAUS)
     . ILOS25T=ILOS24T+MNODET*MGAUST,         ! DERIV(NDIME,NNODE,NGAUS)
     . ILOS26T=ILOS25T+MDIMET*MNODET*MGAUST,  ! POSGP(NDIME,NGAUS)
     . ILOS27T=ILOS26T+MDIMET*MGAUST,         ! WEIGP(NGAUS)
     . ILOS28T=ILOS27T+MGAUST,                ! XJACM(NDIME,NDIME)
     . ILOS29T=ILOS28T+MDIMET*MDIMET,         ! EMASS(NNODE,NNODE)
     . ILOS30T=ILOS29T+MNODET*MNODET,         ! ELDI1(NDOFN,NNODE)
     . ILOS31T=ILOS30T+MDOFCT*MNODET,         ! TEINI(NDOFN,NNODE)
     . ILOS32T=ILOS31T+MDOFCT*MNODET )        ! ADVEM(NDIME,NNODE)
      PARAMETER(
     . ILOS33T=ILOS32T+MDIMET*MNODET*MCONVT,  ! BOUCH(NDOFN,NNODE)
     . ILOS34T=ILOS33T+MDOFCT*MNODET*MMEMO10, ! FPCHL(NFPCH,NNODE)
     . ILOS35T=ILOS34T+MFPCH*MNODET,          ! DVOLI(NGAUS)
     . ILOS36T=ILOS35T+MGAUST*MTER3,          ! ELCOI(NDIME,NNODE)
     . ILOS37T=ILOS36T+MDIMET*MNODET*MTER3,   ! DISPT(NDOFC,NNODE) !m
     . LLOST  =ILOS37T+MDOFCM*MNODET*MTER3 )
C
      PARAMETER(
     . IGSM1T=1,                                    ! SMSTS(NSTR1,NPOIN)
     . IGSM2T=IGSM1T+MSTR1T*MPOINT*MSMOMT*MMEMO4,   ! SMSTN(NDIME,NPOIN)
     . IGSM3T=IGSM2T+MDIMET*MPOINT*MSMOMT*MMEMO4,   ! ACCPN(NPOIN)
     . IGSM4T=IGSM3T+MPOINT*MSMOMT,                 ! STREB(NSTR1,NNODE)
     . IGSM5T=IGSM4T+MSTR1T*MNODET*MSMOMT*MMEMO4,   ! STREA(NNODE,NSTR1)
     . IGSM6T=IGSM5T+MNODET*MSTR1T*MSMOMT*MMEMO4,   ! SFISB(NSTR1,NNODE)
     . IGSM7T=IGSM6T+MSTR1T*MNODET*MSMOMT*MMEMO4,   ! SFISA(NNODE,NSTR1)
     . IGSM8T=IGSM7T+MNODET*MSTR1T*MSMOMT*MMEMO4,   ! SMSTP(NNUIN,NPOIN)
     . IGSM9T=IGSM8T+MNUINT*MPOINT*MSMOMT,          ! SFIPA(NNODE,NNUIN)
     . IGSM10T=IGSM9T+MNODET*MNUINT*MSMOMT,         ! SFIPB(NNUIN,NNODE)
     . IGSM11T=IGSM10T+MNUINT*MNODET*MSMOMT,  ! CARTD(NDIME,NNODE,NGAUS)
     . IGSM12T=IGSM11T+MDIMET*MNODET*MGAUST,  ! DVOLU(NGAUS)
     . IGSM13T=IGSM12T+MGAUST,                ! GPCOD(NDIME,NGAUS)
     . IGSM14T=IGSM13T+MDIMET*MGAUST,         ! SHAPE(NNODE,NGAUS)
     . IGSM15T=IGSM14T+MNODET*MGAUST,         ! DERIV(NDIME,NNODE,NGAUS)
     . IGSM16T=IGSM15T+MDIMET*MNODET*MGAUST ) ! POSGP(NDIME,NGAUS)
      PARAMETER(
     . IGSM17T=IGSM16T+MDIMET*MGAUST,         ! WEIGP(NGAUS)
     . IGSM18T=IGSM17T+MGAUST,                ! XJACM(NDIME,NDIME)
     . IGSM19T=IGSM18T+MDIMET*MDIMET,         ! ELCO1(NDIME,NNODE)
     . IGSM20T=IGSM19T+MDIMET*MNODET,         ! EMASS(NNODE,NNODE)
     . IGSM21T=IGSM20T+MNODET*MNODET,         ! ELDI1(NDOFN,NNODE)
     . IGSM22T=IGSM21T+MDOFCT*MNODET,         ! SMSPP(NPORO,NPOIN)
     . IGSM23T=IGSM22T+MPOROT*MPOINT,         ! SPIPA(NNODE,NPORO)
     . IGSM24T=IGSM23T+MNODET*MPOROT,         ! SPIPB(NPORO,NPORO)
     . IGSM25T=IGSM24T+MPOROT*MNODGT,         ! SPIPC(NPORO,NNODG)
     . IGSM26T=IGSM25T+MPOROT*MNODGT,         ! FPCHL(NFPCH,NNODE)
     . IGSM27T=IGSM26T+MFPCH*MNODET,          ! VELCM(NDOFN,NNODE)
     . IGSM28T=IGSM27T+MDOFCT*MNODET,         ! DVOLI(NGAUS)
     . IGSM29T=IGSM28T+MGAUST*MTER3,          ! ELCOI(NDIME,NNODE)
     . IGSM30T=IGSM29T+MDIMET*MNODET*MTER3,   ! DISPT(NDOFC,NNODE) !m
     . LGSMT  =IGSM30T+MDOFCM*MNODET*MTER3 )
C
      PARAMETER(
     . ISET1T=1,                              ! DERIV(NDIME,NNODE,NGAUS)
     . ISET2T=ISET1T+MDIMET*MNODET*MGAUST,    ! POSGP(NDIME,NGAUS)
     . ISET3T=ISET2T+MGAUST*MDIMET,           ! WEIGP(NGAUS)
     . ISET4T=ISET3T+MGAUST,                  ! XJACM(NDIME,NDIME)
     . ISET5T=ISET4T+MDIMET*MDIMET,           ! RMAT2(NSTR1,NSTR1)
     . ISET6T=ISET5T+MSTR1T*MSTR1T,           ! ELCO1(NDIME,NNODE)
     . ISET7T=ISET6T+MDIMET*MNODET,           ! CARTD(NDIME,NNODE,NGAUS)
     . ISET8T=ISET7T+MDIMET*MNODET*MGAUST,    ! DVOLU(NGAUS)
     . ISET9T=ISET8T+MGAUST,                  ! GPCOD(NDIME,NGAUS)
     . ISET10T=ISET9T+MDIMET*MGAUST,          ! SHAPE(NNODE,NGAUS)
     . ISET11T=ISET10T+MNODET*MGAUST,         ! EMASS(NNODE,NNODE)
     . ISET12T=ISET11T+MNODET*MNODET,         ! ELDI1(NDOFN,NNODE)
     . ISET13T=ISET12T+MDOFCT*MNODET,         ! TEINI(NDOFN,NNODE)
     . ISET14T=ISET13T+MDOFCT*MNODET,         ! FPCHL(NFPCH,NNODE)
     . ISET15T=ISET14T+MFPCH*MNODET,          ! VNORL(NDIME,NGAUS)
     . ISET16T=ISET15T+MDIMET*MGAUST*MTER2X ) ! DVOLI(NGAUS)
      PARAMETER(
     . ISET17T=ISET16T+MGAUST*MTER3,          ! ELCOI(NDIME,NNODE)
     . ISET18T=ISET17T+MDIMET*MNODET*MTER3,   ! DISPT(NDOFC,NNODE) !m
     . ISET19T=ISET18T+MDOFCM*MNODET*MTER3,   ! COOSM(NDIME,NGAUS)
     . ISET20T=ISET19T+MDIMET*MGAUST*MOCOIT,  ! COOSS(NDIME,NGAUS)
     . ISET21T=ISET20T+MDIMET*MGAUST*MOCOIT,  ! VNOSM(NDIME,NGAUS)
     . ISET22T=ISET21T+MDIMET*MGAUST*MOCOIT,  ! GAPMI(NGAUS)
     . LSETT  =ISET22T+MGAUST*MOCOIT )
C
      PARAMETER(
     . ISTI1T=1,                              ! BMATX(NSTR1,NEVAB)
     . ISTI2T=ISTI1T+MSTR1T*MEVABT,           ! DMATX(NSTR1,NSTR1)
     . ISTI3T=ISTI2T+MSTR1T*MSTR1T,           ! SIGMA(NSTR1)
     . ISTI4T=ISTI3T+MSTR1T,                  ! XJACM(NDIME,NDIME)
     . ISTI5T=ISTI4T+MDIMET*MDIMET,           ! SHAPD(NNODE,NGAUS)
     . ISTI6T=ISTI5T+MNODET*MGAUST,           ! DVOLD(NGAUS)
     . ISTI7T=ISTI6T+MGAUST,                  ! POSGP(NDIME,NGAUS)
     . ISTI8T=ISTI7T+MDIMET*MGAUST,           ! WEIGP(NGAUS)
     . ISTI9T=ISTI8T+MGAUST,                  ! DERID(NDIME,NNODE,NGAUS)
     . ISTI10T=ISTI9T+MDIMET*MNODET*MGAUST,   ! XJACM(NDIME,NDIME)
     . ISTI11T=ISTI10T+MDIMET*MDIMET,         ! CARDD(NDIME,NNODE,NGAUS)
     . ISTI12T=ISTI11T+MDIMET*MNODET*MGAUST,  ! GPCDD(NDIME,NGAUS)
     . ISTI13T=ISTI12T+MDIMET*MGAUST,         ! VELCM(NEVAB)
     . ISTI14T=ISTI13T+MEVABT,                ! WSTI1(NKOVA)
     . ISTI15T=ISTI14T+MKOVAT,                ! WSTI2(NKOVA)
     . ISTI16T=ISTI15T+MKOVAT )               ! VELAN(NEVAB)
      PARAMETER(
     . ISTI17T=ISTI16T+MEVABT,                ! DISIT(NEVAB)
     . ISTI18T=ISTI17T+MEVABT,                ! ELCO1(NDIME,NNODE)
     . ISTI19T=ISTI18T+MDIMET*MNODET,         ! WHAPE(NNODE,NGAUS)
     . ISTI20T=ISTI19T+MNODET*MGAUST,         ! HACHE(NNODE)
     . ISTI21T=ISTI20T+MNODET*MGALET,         ! WARTD(NDIME,NNODE,NGAUS)
     . ISTI22T=ISTI21T+MDIMET*MNODET*MGAUST,  ! WERIV(NDIME,NNODE,NGAUS)
     . ISTI23T=ISTI22T+MDIMET*MNODET*MGAUST*MGALET, ! WHADE(NNODE,NGAUS)
     . ISTI24T=ISTI23T+MNODET*MGAUST,               ! CENTR(NDIME,NNODE)
     . ISTI25T=ISTI24T+MDIMET*MNODET*MGALET,        ! ADVEM(NDIME,NNODE)
     . ISTI26T=ISTI25T+MDIMET*MNODET*MCONVT,  ! CARTD(NDIME,NNODE,NGAUS)
     . ISTI27T=ISTI26T+MDIMET*MNODET*MGAUST,  ! DVOLU(NGAUS)
     . ISTI28T=ISTI27T+MGAUST,                ! GPCOD(NDIME,NGAUS)
     . ISTI29T=ISTI28T+MDIMET*MGAUST,         ! SHAPE(NNODE,NGAUS)
     . ISTI30T=ISTI29T+MNODET*MGAUST,         ! DERIV(NDIME,NNODE,NGAUS)
     . ISTI31T=ISTI30T+MDIMET*MNODET*MGAUST,  ! POSGP(NDIME,NGAUS)
     . ISTI32T=ISTI31T+MDIMET*MGAUST )        ! WEIGP(NGAUS)
      PARAMETER(
     . ISTI33T=ISTI32T+MGAUST,                ! XJACM(NDIME,NDIME)
     . ISTI34T=ISTI33T+MDIMET*MDIMET,         ! EMASS(NNODE,NNODE)
     . ISTI35T=ISTI34T+MNODET*MNODET,         ! PREASL(NNODE)
     . ISTI36T=ISTI35T+MNODET*MTER2*MMEMA11,  ! TGAPSL(NNODE)
     . ISTI37T=ISTI36T+MNODET*MTER2*MMEMA11,  ! ESTII(NKOVA)
     . ISTI38T=ISTI37T+MKOVAT,                ! WSTII(NKOVA)
     . ISTI39T=ISTI38T+MKOVAT*MDYNAT,         ! ELDI1(NDOFN,NNODE)
     . ISTI40T=ISTI39T+MDOFCT*MNODET,         ! ELDI1(NDOFN,NNODE)
     . ISTI41T=ISTI40T+MDOFCT*MNODET,         ! CSTI1(NEVAB,NEVAB)
     . ISTI42T=ISTI41T+MEVABT*MEVABT*MMEMO7,  ! BOUCH(NDOFN,NNODE)
     . ISTI43T=ISTI42T+MDOFCT*MNODET*MMEMO10, ! VNORL(NDIME,NGAUS)
     . ISTI44T=ISTI43T+MDIMET*MGAUST*MTER2*MMEMO11,
C                                             ! DISPT(NDOFC,NNODE) !m
     . ISTI45T=ISTI44T+MDOFCM*MNODET*MTER113, ! FPCHL(NFPCH,NNODE)
     . ISTI46T=ISTI45T+MFPCH*MNODET,          ! DVOLI(NGAUS)
     . ISTI47T=ISTI46T+MGAUST*MTER3,          ! ELCOI(NDIME,NNODE)
     . ISTI48T=ISTI47T+MDIMET*MNODET*MTER3,   ! EMATX(NEVAB,NEVAB)
     . ISTI49T=ISTI48T+MEVABT*MEVABT,         ! FMATX(NEVAB,NEVAB)
     . LSTIT  =ISTI49T+MEVABT*MEVABT )
C
C**** ADDSOLT
C
C**** 1) SKYLINE SOLVER
C
      PARAMETER(
     . ISOL1TS=1+LSTIT*MMEMO7,                 ! GSTDI(NEQNS)
     . ISOL2TS=ISOL1TS+MEQNST,                 ! GSTLO(NLAST*IUNSY)
     . ISOL3TS=ISOL2TS+MLASTT*MUNSY2,          ! GSTUP(NLAST)
     . ISOL4TS=ISOL3TS+MLASTT,                 ! CSTIF(NEVAC,NEVAC)
     . ISOL5TS=ISOL4TS+MEVACT*MEVACT,          ! ELOAD(NEQNS)
     . ISOL6TS=ISOL5TS+MEQNST,                 ! LNUEQ(NTOTV)
     . ISOL7TS=ISOL6TS+(MTOTVT*MCHALT+4)/8,    ! LPONT(NTOTV)
     . ISOL8TS=ISOL7TS+(MTOTVT*MCHALT+4)/8,    ! DISIM(NEVAC)
     . ISOL9TS=ISOL8TS+MEVACT,                 ! FOREL(NEVAC)
     . ISOL10TS=ISOL9TS+MEVACT,                ! ALOAD(NEQNS*IFITE)
     . ISOL11TS=ISOL10TS+MEQNST*MFITET,        ! DELTA(NEQNS*IFITE)
     . ISOL12TS=ISOL11TS+MEQNST*MFITET,        ! DISIC(NEQNS*IFITE)
     . ISOL13TS=ISOL12TS+MEQNST*MFITET,        ! LOCAL(NEVAC*IFITE)
     . LSOLTS  =ISOL13TS+(MEVACT*MFITET*MCHALT+4)/8 )
C
      PARAMETER(
     . LSOLT1=MSOLV1*LSOLTS )
C
C**** 2) FRONTAL SOLVER
C
      PARAMETER(
     . ISOL1TF=1+LSTIT*MMEMO7,                 ! EQRHS(NBUFA)
     . ISOL2TF=ISOL1TF+MBUFAT,                 ! EQUAT(NFRON,NBUFA)
     . ISOL3TF=ISOL2TF+MFRONT*MBUFAT,          ! GLOAD(NFRON)
     . ISOL4TF=ISOL3TF+MFRONT,                 ! GSTIF(NSTIF)
     . ISOL5TF=ISOL4TF+MSTIFT,                 ! LOCEL(NEVAC)
     . ISOL6TF=ISOL5TF+(MEVACT*MCHALT+4)/8,    ! NACVA(NFRON)
     . ISOL7TF=ISOL6TF+(MFRONT*MCHALT+4)/8,    ! NAMEV(NBUFA)
     . ISOL8TF=ISOL7TF+(MBUFAT*MCHALT+4)/8,    ! NDEST(NEVAC)
     . ISOL9TF=ISOL8TF+(MEVACT*MCHALT+4)/8,    ! NPIVO(NBUFA)
     . ISOL10TF=ISOL9TF+(MBUFAT*MCHALT+4)/8,   ! VECRV(NFRON)
     . ISOL11TF=ISOL10TF+MFRONT,               ! CSTIF(NEVAC,NEVAC)
     . ISOL12TF=ISOL11TF+MEVACT*MEVACT,       ! EQCOL(IUNSY*NFRON*NBUFA)
     . LSOLTF  =ISOL12TF+MUNSY2*MFRONT*MBUFAT  )
C
      PARAMETER(
     . LSOLT2=MSOLV2*LSOLTF )
C
C**** 3) PCG SOLVER
C
      PARAMETER(
     . ISOL1TP=1+LSTIT*MMEMO7,                 ! GSTDI(NPOIN*NSIZE)
     . ISOL2TP=ISOL1TP+MPOINT*MSIZET,          ! CSTIF(NEVAC,NEVAC)
     . ISOL3TP=ISOL2TP+MEVACT*MEVACT,          ! DISIM(NEVAC)
     . ISOL4TP=ISOL3TP+MEVACT,                 ! FOREL(NEVAC)
     . ISOL5TP=ISOL4TP+MEVACT,                 ! ALOAD(NTOTV)
     . ISOL6TP=ISOL5TP+MTOTVT,                 ! DELTA(NTOTV)
     . LSOLTP =ISOL6TP+MTOTVT )
C
      PARAMETER(
     . LSOLT3=MSOLV3*LSOLTP )
C
C**** 4) GMRES SOLVER
C
      PARAMETER(
     . ISOL1TG=1+LSTIT*MMEMO7,                 ! GSTDI(NEQNS)
     . ISOL2TG=ISOL1TG+MEQNST,                 ! GSTLO(NLAST*IUNSY)
     . ISOL3TG=ISOL2TG+MLASTT*MUNSY2*MCEROT,   ! GSTUP(NLAST)
     . ISOL4TG=ISOL3TG+MLASTT*MCEROT,          ! CSTIF(NEVAC,NEVAC)
     . ISOL5TG=ISOL4TG+MEVACT*MEVACT,          ! ELOAD(NEQNS)
     . ISOL6TG=ISOL5TG+MEQNST*MCEROT,          ! LNUEQ(NTOTV)
     . ISOL7TG=ISOL6TG+(MTOTVT*MCHALT+4)/8,    ! LPONT(NTOTV)
     . ISOL8TG=ISOL7TG+(MTOTVT*MCHALT+4)/8,    ! DISIM(NEVAC)
     . ISOL9TG=ISOL8TG+MEVACT*MCEROT,          ! FOREL(NEVAC)
     . ISOL10TG=ISOL9TG+MEVACT*MCEROT,         ! ALOAD(NEQNS*IFITE)
     . ISOL11TG=ISOL10TG+MEQNST*MFITET*MCEROT, ! DELTA(NEQNS*IFITE)
     . ISOL12TG=ISOL11TG+MEQNST*MFITET*MCEROT, ! DISIC(NEQNS*IFITE)
     . ISOL13TG=ISOL12TG+MEQNST*MFITET*MCEROT, ! LOCAL(NEVAC*IFITE)
     . ISOL14TG=ISOL13TG+(MEVACT*MFITET*MCHALT+4)/8*MCEROT,
     . ISOL15TG=ISOL14TG,                      ! BGMRE(NEQNS)
     . ISOL16TG=ISOL15TG+MEQNST,               ! XGMRE(NEQNS)
     . ISOL17TG=ISOL16TG+MEQNST,               ! WGMRE(NEQNS)
     . LSOLTG  =ISOL17TG+MWORKGT )
C
      PARAMETER(
     . LSOLT4=MSOLV4*LSOLTG )
C
C**** 5) EXPLICIT SOLVER
C
      PARAMETER(
     . ISOL1TE=1+LSTIT*MMEMO7,                 ! CSTIF(NEVAC,NEVAC)
     . ISOL2TE=ISOL1TE+MEVACT*MEVACT,          ! DISIM(NEVAC)
     . ISOL3TE=ISOL2TE+MEVACT,                 ! FOREL(NEVAC)
     . ISOL4TE=ISOL3TE+MEVACT,                 ! CSTII(NEVAC)
     . ISOL5TE=ISOL4TE+MEVACT,                 ! CSTIT(NTOTV)
     . LSOLTE =ISOL5TE+MTOTVT )
C
      PARAMETER(
     . LSOLT5=MSOLV5*LSOLTE )
C
      PARAMETER(
c    . LSOLT=max(LSOLT1,LSOLT2,LSOLT3,LSOLT4,LSOLT5) )      ! sg
     . LSOLT=(LSOLT1+LSOLT2+LSOLT3+LSOLT4+LSOLT5) )         ! PC & linux
C
C**** SCRATCH VECTOR (WORK1T)
C
      PARAMETER(
c    . MWORK1T=max(LINDT,LRENT,LQUAT,LFIXT,LFOTT,LFORT,     ! sg
c    .             LLOST,LGSMT,LSETT,LSTIT,LSOLT) )
     . MWORK1T=(LINDT+LRENT+LQUAT+LFIXT+LFOTT+LFORT+        ! PC & linux
     .          LLOST+LGSMT+LSETT+LSTIT+LSOLT) )
C