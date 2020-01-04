      SUBROUTINE ADDWORT
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE POINTERS FOR TEMPORARY ARRAYS
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
      INCLUDE 'nuec_om.f'     ! thermal-mechanical
      INCLUDE 'nued_om.f'     ! thermal-microstructural
      INCLUDE 'nuef_om.f'     ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      NNODGT=NNODET
      IF(NGAUST.GT.NNODET) NNODGT=NGAUST
C
      ITERME1=1                              ! unidirectional coupled
      ITERME2=1                              ! bidirectional coupled
      ITERME3=1                              ! deformed shape
      IF(ITERME.LT.0) THEN                   ! uncoupled problems
       ITERME1=0                             ! not used
       ITERME2=0
       ITERME3=0
      ENDIF
      IF(ITERME.EQ.0) THEN                   ! unidirectional coupled
       ITERME1=1                             ! not used
       ITERME2=0
       ITERME3=0
      ENDIF
      IF(ITERME.GT.0) THEN
       IF(ITERMD.EQ.0) ITERME3=0
      ENDIF
      NMEMA11=1-NMEMO11
C
      ITERME2X=ITERME2*NMEMO11
      NOCOIX=0
      IF(NOCOIT.GT.0) THEN
       NOCOIX=1
       IF(ITERME2X.EQ.0) ITERME2X=1
      ENDIF
C
      MMTER113=0
      IF(ITERME2.GT.0.AND.(NMEMO11.GT.0.OR.ITERME3.GT.0)) MMTER113=1
C
      IEGX=0
      IF(ICONVT.EQ.1.AND.NMEMO3.EQ.1) THEN
       IF(IEGFPC.EQ.1) IEGX=1
      ENDIF
C
      IEGX1=0
      IF(ITERMEF.GT.0) IEGX1=1   ! ktem1f is not available at this level
C
      LWOR1T=0
C
C**** ROUTINE STARTS
C
      IF(NMEMO1.EQ.0) THEN           ! coordinates in an elemental array
       ISTART(1)=1                               ! COORD(NDIME,NPOIN)
       ISTART(2)=ISTART(1)+NDIMET*NPOINT         ! GPCOD(NDIME,NGAUS)
       LINDAT   =ISTART(2)+NDIMET*NGAUST
      ELSE                           ! coordinates in a global array
       ISTART(1)=1                               ! GPCOD(NDIME,NGAUS)
       LINDAT   =ISTART(1)+NDIMET*NGAUST
      ENDIF
      IF(NMEMO2.NE.0) THEN
       ISTART(3) =LINDAT                      ! CARTD(NDIME,NNODE,NGAUS)
       ISTART(4) =ISTART( 3)+NDIMET*NNODET*NGAUST   ! DVOLU(NGAUS)
       ISTART(5) =ISTART( 4)+NGAUST                 ! GPCOD(NDIME,NGAUS)
       ISTART(6) =ISTART( 5)+NDIMET*NGAUST          ! SHAPE(NNODE,NGAUS)
       ISTART(7) =ISTART( 6)+NNODET*NGAUST    ! DERIV(NDIME,NNODE,NGAUS)
       ISTART(8) =ISTART( 7)+NDIMET*NNODET*NGAUST   ! POSGP(NDIME,NGAUS)
       ISTART(9) =ISTART( 8)+NDIMET*NGAUST          ! WEIGP(NGAUS)
       ISTART(10)=ISTART( 9)+NGAUST                 ! XJACM(NDIME,NDIME)
       ISTART(11)=ISTART(10)+NDIMET*NDIMET          ! ELCO1(NDIME,NNODE)
       ISTART(12)=ISTART(11)+NDIMET*NNODET          ! EMASS(NNODE,NNODE)
       ISTART(13)=ISTART(12)+NNODET*NNODET          ! ELDI1(NDOFN,NNODE)
       ISTART(14)=ISTART(13)+NDOFCT*NNODET          ! VNORL(NDIME,NGAUS)
       ISTART(15)=ISTART(14)+NDIMET*NGAUST*ITERME2*NMEMO11 !DVOLI(NGAUS)
       ISTART(16)=ISTART(15)+NGAUST*ITERME3         ! ELCOI(NDIME,NNODE)
       ISTART(17)=ISTART(16)+NDIMET*NNODET*ITERME3  ! DISPT(NDOFC,NNODE)
       LINDAT    =ISTART(17)+NDOFCM*NNODET*ITERME3
      ENDIF
C
      IF(LINDAT.GT.LWOR1T) LWOR1T=LINDAT
C
      IF(KRENUT.NE.0) THEN                      ! renumbering
       LRENUT=NPOINT*(10*NNODET+3)/2
       IF(LRENUT.GT.LWOR1T) LWOR1T=LRENUT
      ENDIF
C
C**** ROUTINES ANTESO & POSTSO
C
      IQUAST(1)=1                               ! VVECT(NTOTV)
      IQUAST(2)=IQUAST(1)+NTOTVT                ! WVECT(NTOTV)
      LQUAST   =IQUAST(5)+NTOTVT
C
      IF(LQUAST.GT.LWOR1T) LWOR1T=LQUAST
C
C**** ROUTINE AQUAMX
C
      IF(KPORET.NE.2)GOTO 10
C
      IAQUAT(1)=1                               ! PMATX
      IAQUAT(2)=IAQUAT(1)+NDIMET*(NDIMET+1)/2   ! PMAT1
      IAQUAT(3)=IAQUAT(2)+NDIMET*(NDIMET+1)/2   ! BPBMX
      IAQUAT(4)=IAQUAT(3)+NNODET*(NNODET+1)/2   ! BMATX
      LNAQUT   =IAQUAT(4)+NEVABT*NSTR1T
C
      IF(LNAQUT.GT.LWOR1T) LWOR1T=LNAQUT
C
   10 CONTINUE
C
C**** ROUTINE SOLADDT
C
      IFIXYT(1)=1                               !  LEQNS(NEVAC)
      IFIXYT(2)=IFIXYT(1)+(NEVACT*ICHALT+4)/8   !  LNUEQ(NTOTV)
      IFIXYT(3)=IFIXYT(2)+(NTOTVT*ICHALT+4)/8   !  LOCEL(NEVAC)
      IFIXYT(4)=IFIXYT(3)+(NEVACT*ICHALT+4)/8   !  LPONT(NTOTV)
      IFIXYT(5)=IFIXYT(4)+(NTOTVT*ICHALT+4)/8   !  NACVA(NTOTV)
      IFIXYT(6)=IFIXYT(5)+(NTOTVT*ICHALT+4)/8   !  NDEST(NEVAC)
      IFIXYT(7)=IFIXYT(6)+(NEVACT*ICHALT+4)/8   !  NDFRO(NELEM)
      IFIXYT(8)=IFIXYT(7)+(NELEMT*ICHALT+4)/8   !  PRESC(NDOFC)
      LFIXYT   =IFIXYT(8)+NDOFCT
C
      IF(LFIXYT.GT.LWOR1T) LWOR1T=LFIXYT
C
C**** ROUTINE FORCIN
C
      IFORCT( 1)=1                              ! BMSIG(NEVAB)
      IFORCT( 2)=IFORCT( 1)+NEVABT              ! BMATX(NSTR1,NEVAB)
      IFORCT( 3)=IFORCT( 2)+NSTR1T*NEVABT       ! DESIG(NSTR1) 
      IFORCT( 4)=IFORCT( 3)+NSTR1T              ! DMATX(NSTR1,NSTR1)
      IFORCT( 5)=IFORCT( 4)+NSTR1T*NSTR1T       ! DSTRA(NSTR1)
      IFORCT( 6)=IFORCT( 5)+NSTR1T              ! PRESG(NSTR1)
      IFORCT( 7)=IFORCT( 6)+NSTR1T              ! SGTOT(NSTR1)
      IFORCT( 8)=IFORCT( 7)+NSTR1T              ! SIGMA(NSTR1)
      IFORCT( 9)=IFORCT( 8)+NSTR1T              ! TSTRA(NSTR1)
      IFORCT(10)=IFORCT( 9)+NSTR1T              ! XJACM(NDIME,NDIME)
      IFORCT(11)=IFORCT(10)+NDIMET*NDIMET       ! ELELT(NEVAB)
      IFORCT(12)=IFORCT(11)+NEVABT              ! VELCM(NEVAB)
      IFORCT(13)=IFORCT(12)+NEVABT              ! DISPT(NDOFC,NNODET)
      IFORCT(14)=IFORCT(13)+NDOFCM*NNODET*MMTER113  ! PREASL(NNODET)
      IFORCT(15)=IFORCT(14)+NNODET*ITERME2*NMEMA11  ! TGAPSL(NNODET)
      IFORCT(16)=IFORCT(15)+NNODET*ITERME2*NMEMA11  ! ELCO1(NDIME,NNODE)
C
      IAUXX=0
      IF(ICONVT.EQ.1.OR.IGALET.EQ.1) IAUXX=1
      IFORCT(17)=IFORCT(16)+NDIMET*NNODET           ! SHAPD(NNODE,NGAUS)
      IFORCT(18)=IFORCT(17)+NNODET*NGAUST*IAUXX     ! DVOLD(NGAUS)
      IFORCT(19)=IFORCT(18)+NGAUST*ICONVT           ! POSGP(NDIME,NGAUS)
      IFORCT(20)=IFORCT(19)+NDIMET*NGAUST*IAUXX     ! WEIGP(NGAUS)
      IFORCT(21)=IFORCT(20)+NGAUST*IAUXX      ! DERID(NDIME,NNODE,NGAUS)
      IFORCT(22)=IFORCT(21)+NDIMET*NNODET*NGAUST*IAUXX
C                                                   ! XJACM(NDIME,NDIME)
      IFORCT(23)=IFORCT(22)+NDIMET*NDIMET*IAUXX
C                                             ! CARDD(NDIME,NNODE,NGAUS)
      IFORCT(24)=IFORCT(23)+NDIMET*NNODET*NGAUST*IAUXX
C                                                   ! GPCDD(NDIME,NGAUS)
      IFORCT(25)=IFORCT(24)+NDIMET*NGAUST*IAUXX     ! ELEL1(NEVAB)
      IFORCT(26)=IFORCT(25)+NEVABT*ICONVT           ! ELEL2(NEVAB)
C
      IFORCT(27)=IFORCT(26)+NEVABT*ICONVT           ! WHAPE(NNODE,NGAUS)
      IFORCT(28)=IFORCT(27)+NNODET*NGAUST           ! HACHE(NNODE)
      IFORCT(29)=IFORCT(28)+NNODET*IGALET     ! WARTD(NDIME,NNODE,NGAUS)
      IFORCT(30)=IFORCT(29)+NDIMET*NNODET*NGAUST
C                                             ! WERIV(NDIME,NNODE,NGAUS)
      IFORCT(31)=IFORCT(30)+NDIMET*NNODET*NGAUST*IGALET
C                                                   ! WHADE(NNODE,NGAUS)
      IFORCT(32)=IFORCT(31)+NNODET*NGAUST*ICONVT    ! CENTR(NDIME,NNODE)
C
      IAUX1=0
      IF(IMICR.EQ.1) IAUX1=1
      IFORCT(33)=IFORCT(32)+NDIMET*NNODET*IGALET    ! VEL1M(NEVAB)
      IFORCT(34)=IFORCT(33)+NEVABT*IAUX1            ! ADVEM(NDIME*NNODE)
C
      IFORCT(35)=IFORCT(34)+NDIMET*NNODET*ICONVT
C                                             ! CARTD(NDIME,NNODE,NGAUS)
      IFORCT(36)=IFORCT(35)+NDIMET*NNODET*NGAUST    ! DVOLU(NGAUS)
      IFORCT(37)=IFORCT(36)+NGAUST                  ! GPCOD(NDIME,NGAUS)
      IFORCT(38)=IFORCT(37)+NDIMET*NGAUST           ! SHAPE(NNODE,NGAUS)
      IFORCT(39)=IFORCT(38)+NNODET*NGAUST     ! DERIV(NDIME,NNODE,NGAUS)
      IFORCT(40)=IFORCT(39)+NDIMET*NNODET*NGAUST    ! POSGP(NDIME,NGAUS)
      IFORCT(41)=IFORCT(40)+NDIMET*NGAUST           ! WEIGP(NGAUS)
      IFORCT(42)=IFORCT(41)+NGAUST                  ! XJACM(NDIME,NDIME)
      IFORCT(43)=IFORCT(42)+NDIMET*NDIMET           ! EMASS(NNODE,NNODE)
      IFORCT(44)=IFORCT(43)+NNODET*NNODET           ! ELDI1(NDOFN,NNODE)
      IFORCT(45)=IFORCT(44)+NDOFCT*NNODET           ! TEINI(NDOFN,NNODE)
      IFORCT(46)=IFORCT(45)+NDOFCT*NNODET           ! BOUCH(NDOFN,NNODE)
      IFORCT(47)=IFORCT(46)+NDOFCT*NNODET*NMEMO10   ! VNORL(NDIME,NGAUS)
      IFORCT(48)=IFORCT(47)+NDIMET*NGAUST*ITERME2*NMEMO11
C                                                   ! FPCHL(NFPCH,NNODE)
      IFORCT(49)=IFORCT(48)+NFPCH*NNODET            ! DVOLI(NGAUS)
      IFORCT(50)=IFORCT(49)+NGAUST*ITERME3          ! ELCOI(NDIME,NNODE)
      IFORCT(51)=IFORCT(50)+NDIMET*NNODET*ITERME3   ! ACCPN(NPOIN)
      IFORCT(52)=IFORCT(51)+NPOINT*IEGX             ! SMSTP(2,NPOIN)
      IFORCT(53)=IFORCT(52)+2*NPOINT*IEGX           ! SFIPA(NNODE,2)
      IFORCT(54)=IFORCT(53)+NNODET*2*IEGX           ! SFIPB(2,NNODE)
      IFORCT(55)=IFORCT(54)+2*NNODET*IEGX           ! ELDIQ(NDOFN,NNODE)
      LFORCT    =IFORCT(55)+NDOFCT*NNODET*IEGX1
C
      IF(LFORCT.GT.LWOR1T) LWOR1T=LFORCT
C
C**** ROUTINE FORCDY AND THERDY
C
      IF(KDYNAT.EQ.1)THEN
       IFORDT(1)=1                              ! ELELM(NEVAB)
       IFORDT(2)=IFORDT(1)+NEVABT               ! VELCM(NEVAB)
       IFORDT(3)=IFORDT(2)+NEVABT               ! DISIM(NEVAB)
       IFORDT(4)=IFORDT(3)+NEVABT               ! ESTIF(NKOVA)
       IFORDT(5)=IFORDT(4)+NKOVAT               ! DSTIF(NKOVA)
       IFORDT(6)=IFORDT(5)+NKOVAT               ! ELELT(NEVAB)
       IFORDT(7)=IFORDT(6)+NEVABT               ! SHAPD(NNODE,NGAUS)
       IFORDT(8)=IFORDT(7)+NNODET*NGAUST        ! DVOLD(NGAUS)
       IFORDT(9)=IFORDT(8)+NGAUST               ! POSGP(NDIME,NGAUS)
       IFORDT(10)=IFORDT(9)+NDIMET*NGAUST       ! WEIGP(NGAUS)
       IFORDT(11)=IFORDT(10)+NGAUST           ! DERID(NDIME,NNODE,NGAUS)
       IFORDT(12)=IFORDT(11)+NDIMET*NNODET*NGAUST ! XJACM(NDIME,NDIME)
       IFORDT(13)=IFORDT(12)+NDIMET*NDIMET    ! CARDD(NDIME,NNODE,NGAUS)
       IFORDT(14)=IFORDT(13)+NDIMET*NNODET*NGAUST ! GPCDD(NDIME,NGAUS)
       IFORDT(15)=IFORDT(14)+NDIMET*NGAUST      ! ELEL1(NEVAB)
       IFORDT(16)=IFORDT(15)+NEVABT             ! ELEL2(NEVAB)
       IFORDT(17)=IFORDT(16)+NEVABT             ! ELCO1(NDIME,NNODE)
       IFORDT(18)=IFORDT(17)+NDIMET*NNODET          ! WHAPE(NNODE,NGAUS)
       IFORDT(19)=IFORDT(18)+NNODET*NGAUST          ! HACHE(NNODE)
       IFORDT(20)=IFORDT(19)+NNODET*IGALET          ! WHADE(NNODE,NGAUS)
       IFORDT(21)=IFORDT(20)+NNODET*NGAUST          ! CENTR(NDIME,NNODE)
       IFORDT(22)=IFORDT(21)+NDIMET*NNODET*IGALET   ! ADVEM(NDIME,NNODE)
       IFORDT(23)=IFORDT(22)+NDIMET*NNODET*ICONVT
C                                             ! CARTD(NDIME,NNODE,NGAUS)
       IFORDT(24)=IFORDT(23)+NDIMET*NNODET*NGAUST   ! DVOLU(NGAUS)
       IFORDT(25)=IFORDT(24)+NGAUST                 ! GPCOD(NDIME,NGAUS)
       IFORDT(26)=IFORDT(25)+NDIMET*NGAUST          ! SHAPE(NNODE,NGAUS)
       IFORDT(27)=IFORDT(26)+NNODET*NGAUST    ! DERIV(NDIME,NNODE,NGAUS)
       IFORDT(28)=IFORDT(27)+NDIMET*NNODET*NGAUST   ! POSGP(NDIME,NGAUS)
       IFORDT(29)=IFORDT(28)+NDIMET*NGAUST          ! WEIGP(NGAUS)
       IFORDT(30)=IFORDT(29)+NGAUST                 ! XJACM(NDIME,NDIME)
       IFORDT(31)=IFORDT(30)+NDIMET*NDIMET          ! EMASS(NNODE,NNODE)
       IFORDT(32)=IFORDT(31)+NNODET*NNODET          ! ELDI1(NDOFN,NNODE)
       IFORDT(33)=IFORDT(32)+NDOFCT*NNODET          ! TEINI(NDOFN,NNODE)
       IFORDT(34)=IFORDT(33)+NDOFCT*NNODET          ! FPCHL(NFPCH,NNODE)
       IFORDT(35)=IFORDT(34)+NFPCH*NNODET           ! DVOLI(NGAUS)
       IFORDT(36)=IFORDT(35)+NGAUST*ITERME3         ! ELCOI(NDIME,NNODE)
       IFORDT(37)=IFORDT(36)+NDIMET*NNODET*ITERME3  ! DISPT(NDOFC,NNODE)
       LFORDT    =IFORDT(37)+NDOFCM*NNODET*ITERME3
C
       IF(LFORDT.GT.LWOR1T) LWOR1T=LFORDT
      ENDIF                                  ! kdynat.eq.1
C
C**** ROUTINE FORCWP
C
      IF(KPORET.EQ.2) THEN
       IFORWT( 1)=1                            ! HSTIF(NEVAB,NNODE)
       IFORWT( 2)=IFORWT( 1)+NEVABT*NNODET     ! PSTIF(NKOND)
       IFORWT( 3)=IFORWT( 2)+NKONDT            ! QSTIF(NKOND)
       LFORWT    =IFORWT( 3)+NKONDT
C
       IF(LFORWT.GT.LWOR1T) LWOR1T=LFORWT
      ENDIF
C
C**** ROUTINE LOADPST & FIXITYT
C
      NGAULT=4          
      ILOSUT( 1)=1                              ! ALOAD(NEVAB)
      ILOSUT( 2)=ILOSUT( 1)+NEVABT            ! DERIV(NDIME,NNODE,NGAUS)
      ILOSUT( 3)=ILOSUT( 2)+NDIMET*NNODET*NGAUST    ! ELEDG(NDIME,NNODE)
      ILOSUT( 4)=ILOSUT( 3)+NDIMET*NNODET       ! GPCOD(NDIME)
      ILOSUT( 5)=ILOSUT( 4)+NDIMET              ! GVECT(NDIME)
      ILOSUT( 6)=ILOSUT( 5)+NDIMET              ! NOPRS(NNODE)
      ILOSUT( 7)=ILOSUT( 6)+NNODET              ! POSGP(NDIME,NGAUS)
      ILOSUT( 8)=ILOSUT( 7)+NDIMET*NGAUST       ! PRESS(NNODE,NDIME)
      ILOSUT( 9)=ILOSUT( 8)+NNODET*NDIMET       ! PVECT(NDIME)
      ILOSUT(10)=ILOSUT( 9)+NDIMET              ! SHAPE(NNODE)
      ILOSUT(11)=ILOSUT(10)+NNODET              ! VAREA(NDIME)
      ILOSUT(12)=ILOSUT(11)+NDIMET              ! WEIGP(NGAUS)
      ILOSUT(13)=ILOSUT(12)+NGAUST              ! XJACM(NDIME*NDIME)
      ILOSUT(14)=ILOSUT(13)+NDIMET*NDIMET     ! CARTD(NDIME,NNODE,NGAUS)
      ILOSUT(15)=ILOSUT(14)+NDIMET*NNODET*NGAUST ! PRESC(NDOFC)
      ILOSUT(16)=ILOSUT(15)+NDOFCT              ! VELCM(NEVAB)
      ILOSUT(17)=ILOSUT(16)+NEVABT              ! ELCO1(NDIME,NNODE)
      ILOSUT(18)=ILOSUT(17)+NDIMET*NNODET       ! WHAPE(NNODE,NGAUS)
      ILOSUT(19)=ILOSUT(18)+NNODET*NGAUST       ! HACHE(NNODE)
      ILOSUT(20)=ILOSUT(19)+NNODET*IGALET       ! CENTR(NDIME,NNODE)
      ILOSUT(21)=ILOSUT(20)+NDIMET*NNODET*IGALET
C                                             ! CARTD(NDIME,NNODE,NGAUS)
      ILOSUT(22)=ILOSUT(21)+NDIMET*NNODET*NGAUST    ! DVOLU(NGAUS)
      ILOSUT(23)=ILOSUT(22)+NGAUST                  ! GPCOD(NDIME,NGAUS)
      ILOSUT(24)=ILOSUT(23)+NDIMET*NGAUST           ! SHAPE(NNODE,NGAUS)
      ILOSUT(25)=ILOSUT(24)+NNODET*NGAUST     ! DERIV(NDIME,NNODE,NGAUS)
      ILOSUT(26)=ILOSUT(25)+NDIMET*NNODET*NGAUST    ! POSGP(NDIME,NGAUS)
      ILOSUT(27)=ILOSUT(26)+NDIMET*NGAUST           ! WEIGP(NGAUS)
      ILOSUT(28)=ILOSUT(27)+NGAUST                  ! XJACM(NDIME,NDIME)
      ILOSUT(29)=ILOSUT(28)+NDIMET*NDIMET           ! EMASS(NNODE,NNODE)
      ILOSUT(30)=ILOSUT(29)+NNODET*NNODET           ! ELDI1(NDOFN,NNODE)
      ILOSUT(31)=ILOSUT(30)+NDOFCT*NNODET           ! TEINI(NDOFN,NNODE)
      ILOSUT(32)=ILOSUT(31)+NDOFCT*NNODET           ! ADVEM(NDIME,NNODE)
      ILOSUT(33)=ILOSUT(32)+NDIMET*NNODET*ICONVT    ! BOUCH(NDOFN,NNODE)
      ILOSUT(34)=ILOSUT(33)+NDOFCT*NNODET*NMEMO10   ! FPCHL(NFPCH,NNODE)
      ILOSUT(35)=ILOSUT(34)+NFPCH*NNODET            ! DVOLI(NGAUS)
      ILOSUT(36)=ILOSUT(35)+NGAUST*ITERME3          ! ELCOI(NDIME,NNODE)
      ILOSUT(37)=ILOSUT(36)+NDIMET*NNODET*ITERME3   ! DISPT(NDOFC,NNODE)
      LLOSUT    =ILOSUT(37)+NDOFCM*NNODET*ITERME3
C
      IF(LLOSUT.GT.LWOR1T) LWOR1T=LLOSUT
C
C**** ROUTINE OUTSMOT & OUTGAUT
C
      ISMOMT=0
      IF(KSGAUT.NE.0) ISMOMT=1
C
      IGSMOT(1)=1                                   ! SMSTS(NSTR1,NPOIN)
      IGSMOT(2)=IGSMOT(1)+NSTR1T*NPOINT*ISMOMT*NMEMO4
C                                                   ! SMSTN(NSTR1,NPOIN)
      IGSMOT(3)=IGSMOT(2)+NSTR1T*NPOINT*ISMOMT*NMEMO4
C                                                   ! ACCPN(NPOIN)
C
C**** NODAL FLUXES:      STREB(NSTR1,NNODE) & STREA(NNODE,NSTR1)
C     NODAL TEMP. GRAD.: SFISB(NSTR1,NNODE) & SFISA(NNODE,NSTR1)
C
      IGSMOT(4)=IGSMOT(3)+NPOINT*ISMOMT             ! STREB(NSTR1,NNODE)
      IGSMOT(5)=IGSMOT(4)+NSTR1T*NNODET*ISMOMT*NMEMO4
C                                                   ! STREA(NNODE,NSTR1)
      IGSMOT(6)=IGSMOT(5)+NNODET*NSTR1T*ISMOMT*NMEMO4
C                                                   ! SFISB(NSTR1,NNODE)
      IGSMOT(7)=IGSMOT(6)+NSTR1T*NNODET*ISMOMT*NMEMO4
C                                                   ! SFISA(NNODE,NSTR1)
C
C**** NODAL INT. VARIABLES: SFIPB(NNUIN,NNODE) & SFIPA(NNODE,NNUIN)
C
      IGSMOT( 8)=IGSMOT( 7)+NNODET*NSTR1T*ISMOMT*NMEMO4
C                                                   ! SMSTP(NNUIN,NPOIN)
      IGSMOT( 9)=IGSMOT( 8)+NPOINT*NNUINT*ISMOMT    ! SFIPA(NNODE,NNUIN)
      IGSMOT(10)=IGSMOT( 9)+NNODET*NNUINT*ISMOMT    ! SFIPB(NNUIN,NNODE)
C
      IGSMOT(11)=IGSMOT(10)+NNUINT*NNODET*ISMOMT
C                                             ! CARTD(NDIME,NNODE,NGAUS)
      IGSMOT(12)=IGSMOT(11)+NDIMET*NNODET*NGAUST    ! DVOLU(NGAUS)
      IGSMOT(13)=IGSMOT(12)+NGAUST                  ! GPCOD(NDIME,NGAUS)
      IGSMOT(14)=IGSMOT(13)+NDIMET*NGAUST           ! SHAPE(NNODE,NGAUS)
      IGSMOT(15)=IGSMOT(14)+NNODET*NGAUST     ! DERIV(NDIME,NNODE,NGAUS)
      IGSMOT(16)=IGSMOT(15)+NDIMET*NNODET*NGAUST    ! POSGP(NDIME,NGAUS)
      IGSMOT(17)=IGSMOT(16)+NDIMET*NGAUST           ! WEIGP(NGAUS)
      IGSMOT(18)=IGSMOT(17)+NGAUST                  ! XJACM(NDIME,NDIME)
      IGSMOT(19)=IGSMOT(18)+NDIMET*NDIMET           ! ELCO1(NDIME,NNODE)
      IGSMOT(20)=IGSMOT(19)+NDIMET*NNODET           ! EMASS(NNODE,NNODE)
      IGSMOT(21)=IGSMOT(20)+NNODET*NNODET           ! ELDI1(NDOFN,NNODE)
      IGSMOT(22)=IGSMOT(21)+NDOFCT*NNODET           ! SMSPP(NPORO,NPOIN)
      IGSMOT(23)=IGSMOT(22)+NPOROT*NPOINT           ! SPIPA(NNODE,NPORO)
      IGSMOT(24)=IGSMOT(23)+NNODET*NPOROT           ! SPIPB(NPORO,NNODG)
      IGSMOT(25)=IGSMOT(24)+NPOROT*NNODGT           ! SPIBC(NPORO,NNODG)
      IGSMOT(26)=IGSMOT(25)+NPOROT*NNODGT           ! FPCHL(NFPCH,NNODE)
      IGSMOT(27)=IGSMOT(26)+NFPCH*NNODET            ! VELCM(NDOFN,NNODE)
      IGSMOT(28)=IGSMOT(27)+NDOFCT*NNODET           ! DVOLI(NGAUS)
      IGSMOT(29)=IGSMOT(28)+NGAUST*ITERME3          ! ELCOI(NDIME,NNODE)
      IGSMOT(30)=IGSMOT(29)+NDIMET*NNODET*ITERME3   ! DISPT(NDOFC,NNODE)
      LGSMOT    =IGSMOT(30)+NDOFCM*NNODET*ITERME3
C
      IF(LGSMOT.GT.LWOR1T) LWOR1T=LGSMOT
C
C**** ROUTINE SETMTX
C
      ISETMT( 1)=1                            ! DERIV(NDIME,NNODE,NGAUS)
      ISETMT( 2)=ISETMT( 1)+NDIMET*NNODET*NGAUST    ! POSGP(NDIME,NGAUS)
      ISETMT( 3)=ISETMT( 2)+NGAUST*NDIMET           ! WEIGP(NGAUS)
      ISETMT( 4)=ISETMT( 3)+NGAUST                  ! XJACM(NDIME,NDIME)
      ISETMT( 5)=ISETMT( 4)+NDIMET*NDIMET           ! RMAT2(NSTR1,NSTR1)
      ISETMT( 6)=ISETMT( 5)+NSTR1T*NSTR1T           ! ELCO1(NDIME,NNODE)
      ISETMT( 7)=ISETMT( 6)+NDIMET*NNODET     ! CARTD(NDIME,NNODE,NGAUS)
      ISETMT( 8)=ISETMT( 7)+NDIMET*NNODET*NGAUST    ! DVOLU(NGAUS)
      ISETMT( 9)=ISETMT( 8)+NGAUST                  ! GPCOD(NDIME,NGAUS)
      ISETMT(10)=ISETMT( 9)+NDIMET*NGAUST           ! SHAPE(NNODE,NGAUS)
      ISETMT(11)=ISETMT(10)+NNODET*NGAUST           ! EMASS(NNODE,NNODE)
      ISETMT(12)=ISETMT(11)+NNODET*NNODET           ! ELDI1(NDOFN,NNODE)
      ISETMT(13)=ISETMT(12)+NDOFCT*NNODET           ! TEINI(NDOFN,NNODE)
      ISETMT(14)=ISETMT(13)+NDOFCT*NNODET           ! FPCHL(NFPCH,NNODE)
      ISETMT(15)=ISETMT(14)+NFPCH*NNODET            ! VNORL(NDIME,NGAUS)
      ISETMT(16)=ISETMT(15)+NDIMET*NGAUST*ITERME2X  ! DVOLI(NGAUS)
      ISETMT(17)=ISETMT(16)+NGAUST*ITERME3          ! ELCOI(NDIME,NNODE)
      ISETMT(18)=ISETMT(17)+NDIMET*NNODET*ITERME3   ! DISPT(NDOFC,NNODE)
      ISETMT(19)=ISETMT(18)+NDOFCM*NNODET*ITERME3   ! COOSM(NDIME,NGAUS)
      ISETMT(20)=ISETMT(19)+NDIMET*NGAUST*NOCOIX    ! COOSS(NDIME,NGAUS)
      ISETMT(21)=ISETMT(20)+NDIMET*NGAUST*NOCOIX    ! VNOSM(NDIME,NGAUS)
      ISETMT(22)=ISETMT(21)+NDIMET*NGAUST*NOCOIX    ! GAPMI(NGAUS)
      LSETMT    =ISETMT(22)+NGAUST*NOCOIX
C
      IF(LSETMT.GT.LWOR1T) LWOR1T=LSETMT
C
C**** ROUTINE STIFMX
C
      ISTIFT(1)=1                               ! BMATX(NSTR1,NEVAB)
      ISTIFT(2)=ISTIFT(1)+NSTR1T*NEVABT         ! DMATX(NSTR1,NSTR1)
      ISTIFT(3)=ISTIFT(2)+NSTR1T*NSTR1T         ! SIGMA(NSTR1)
      ISTIFT(4)=ISTIFT(3)+NSTR1T                ! XJACM(NDIME,NDIME)
      ISTIFT(5)=ISTIFT(4)+NDIMET*NDIMET         ! SHAPD(NNODE,NGAUS)
      ISTIFT(6)=ISTIFT(5)+NNODET*NGAUST         ! DVOLD(NGAUS)
      ISTIFT(7)=ISTIFT(6)+NGAUST                ! POSGP(NDIME,NGAUS)
      ISTIFT(8)=ISTIFT(7)+NDIMET*NGAUST         ! WEIGP(NGAUS)
      ISTIFT(9)=ISTIFT(8)+NGAUST              ! DERID(NDIME,NNODE,NGAUS)
      ISTIFT(10)=ISTIFT(9)+NDIMET*NNODET*NGAUST ! XJACM(NDIME,NDIME)
      ISTIFT(11)=ISTIFT(10)+NDIMET*NDIMET     ! CARDD(NDIME,NNODE,NGAUS)
      ISTIFT(12)=ISTIFT(11)+NDIMET*NNODET*NGAUST ! GPCDD(NDIME,NGAUS)
      ISTIFT(13)=ISTIFT(12)+NDIMET*NGAUST       ! VELCM(NEVAB)
      ISTIFT(14)=ISTIFT(13)+NEVABT              ! WSTI1(NKOVA)
      ISTIFT(15)=ISTIFT(14)+NKOVAT              ! WSTI2(NKOVA)
      ISTIFT(16)=ISTIFT(15)+NKOVAT              ! VELAN(NEVAB)
      ISTIFT(17)=ISTIFT(16)+NEVABT              ! DISIT(NEVAB)
      ISTIFT(18)=ISTIFT(17)+NEVABT              ! ELCO1(NDIME,NNODE)
      ISTIFT(19)=ISTIFT(18)+NDIMET*NNODET       ! WHAPE(NNODE,NGAUS)
      ISTIFT(20)=ISTIFT(19)+NNODET*NGAUST       ! HACHE(NNODE)
      ISTIFT(21)=ISTIFT(20)+NNODET*IGALET     ! WARTD(NDIME,NNODE,NGAUS)
      ISTIFT(22)=ISTIFT(21)+NDIMET*NNODET*NGAUST
C                                             ! WERIV(NDIME,NNODE,NGAUS)
      ISTIFT(23)=ISTIFT(22)+NDIMET*NNODET*NGAUST*IGALET
C                                                   ! WHADE(NNODE,NGAUS)
      ISTIFT(24)=ISTIFT(23)+NNODET*NGAUST           ! CENTR(NDIME,NNODE)
      ISTIFT(25)=ISTIFT(24)+NDIMET*NNODET*IGALET    ! ADVEM(NDIME*NNODE)
      ISTIFT(26)=ISTIFT(25)+NDIMET*NNODET*ICONVT
C                                             ! CARTD(NDIME,NNODE,NGAUS)
      ISTIFT(27)=ISTIFT(26)+NDIMET*NNODET*NGAUST    ! DVOLU(NGAUS)
      ISTIFT(28)=ISTIFT(27)+NGAUST                  ! GPCOD(NDIME,NGAUS)
      ISTIFT(29)=ISTIFT(28)+NDIMET*NGAUST           ! SHAPE(NNODE,NGAUS)
      ISTIFT(30)=ISTIFT(29)+NNODET*NGAUST     ! DERIV(NDIME,NNODE,NGAUS)
      ISTIFT(31)=ISTIFT(30)+NDIMET*NNODET*NGAUST    ! POSGP(NDIME,NGAUS)
      ISTIFT(32)=ISTIFT(31)+NDIMET*NGAUST           ! WEIGP(NGAUS)
      ISTIFT(33)=ISTIFT(32)+NGAUST                  ! XJACM(NDIME,NDIME)
      ISTIFT(34)=ISTIFT(33)+NDIMET*NDIMET           ! EMASS(NNODE,NNODE)
      ISTIFT(35)=ISTIFT(34)+NNODET*NNODET           ! PREASL(NNODE)
      ISTIFT(36)=ISTIFT(35)+NNODET*ITERME2*NMEMA11  ! TGAPSL(NNODE)
      ISTIFT(37)=ISTIFT(36)+NNODET*ITERME2*NMEMA11  ! ESTII(NKOVA)
      ISTIFT(38)=ISTIFT(37)+NKOVAT                  ! WSTII(NKOVA)
      ISTIFT(39)=ISTIFT(38)+NKOVAT*KDYNAT           ! ELDI1(NDOFN,NNODE)
      ISTIFT(40)=ISTIFT(39)+NDOFCT*NNODET           ! TEINI(NDOFN,NNODE)
      ISTIFT(41)=ISTIFT(40)+NDOFCT*NNODET           ! CSTI1(NEVAB,NEVAB)
      ISTIFT(42)=ISTIFT(41)+NEVABT*NEVABT*NMEMO7    ! BOUCH(NDOFN,NNODE)
      ISTIFT(43)=ISTIFT(42)+NDOFCT*NNODET*NMEMO10   ! VNORL(NDIME,NGAUS)
      ISTIFT(44)=ISTIFT(43)+NDIMET*NGAUST*ITERME2*NMEMO11
C                                                   ! DISPT(NDOFC,NNODE)
      ISTIFT(45)=ISTIFT(44)+NDOFCM*NNODET*MMTER113  ! FPCHL(NFPCH,NNODE)
      ISTIFT(46)=ISTIFT(45)+NFPCH*NNODET            ! DVOLI(NGAUS)
      ISTIFT(47)=ISTIFT(46)+NGAUST*ITERME3          ! ELCOI(NDIME,NNODE)
      ISTIFT(48)=ISTIFT(47)+NDIMET*NNODET*ITERME3   ! EMATX(NEVAB,NEVAB)
      ISTIFT(49)=ISTIFT(48)+NEVABT*NEVABT           ! FMATX(NEVAB,NEVAB)
      LSTIFT    =ISTIFT(49)+NEVABT*NEVABT
C
      IF(LSTIFT.GT.LWOR1T) LWOR1T=LSTIFT
C
C**** MEMORY CONTROL
C
      IF(LWOR1T.GT.MWORK1T)
     . CALL RUNENDT('ERROR IN ADDWORT: LWOR1T          ')
C
C**** OTHER CONTROLS
C
      IF(KDYNAT.GT.MDYNAT)
     . CALL RUNENDT('ERROR IN ADDWORT: KDYNAT          ')
      IF(ISMOMT.GT.MSMOMT)
     . CALL RUNENDT('ERROR IN ADDWORT: ISMOMT          ')
      IF(KRENUT.GT.MRENUT)
     . CALL RUNENDT('ERROR IN ADDWORT: KRENUT          ')
      IF(KPORET.GT.MPORET)
     . CALL RUNENDT('ERROR IN ADDWORT: KPORET          ')
      IF(NHOURT.GT.MHOURT)
     . CALL RUNENDT('ERROR IN ADDWORT: NHOURT          ')
C
      ITERMA=1
      IF(ITERME.LE.0) ITERMA=0                   ! uncoupled prob.
C
      IF(ITERMA.GT.MTERME)
     . CALL RUNENDT('ERROR IN ADDWORT: ITERME          ')
      IF(IMICR.GT.MMICR)
     . CALL RUNENDT('ERROR IN ADDWORT: IMICR           ')
      IF(IMICR.EQ.1) THEN
       IF(IMICO.LT.MMICO)
     .  CALL RUNENDT('ERROR IN ADDWORT: IMICO          ')
      ENDIF
      IF(ICONVT.GT.MCONVT)
     . CALL RUNENDT('ERROR IN ADDWORT: ICONVT          ')
      IF(IGALET.GT.MGALET)
     . CALL RUNENDT('ERROR IN ADDWORT: IGALET          ')
C
      IF(NMEMO3.EQ.0.AND.IMICR.NE.0)
     . CALL RUNENDT('ERROR: NMEMO3=0 WITH IMICR=1')
      IF(NMEMO3.EQ.0.AND.NITERC.EQ.1)
     . CALL RUNENDT('ERROR: NMEMO3=0 WITH NITERC=1')
C
      IF(NPOROT.GT.MPOROT)
     . CALL RUNENDT('ERROR: NPOROT GT MPOROT      ')
      IF(NFPCH.GT.MFPCH)
     . CALL RUNENDT('ERROR: NFPCH GT MFPCH        ')
C
      IF(IMICR.EQ.1) THEN
       IF(NNUM4T.GT.MNUM4T)
     .  CALL RUNENDT('ERROR: NNUM4T GT MNUM4T      ')
       IF(NNUM4T.GT.MNU4TX)
     .  CALL RUNENDT('ERROR: NNUM4T GT MNU4TX      ')
       IF(NNUM4TM.GT.MNUM4TM)
     .  CALL RUNENDT('ERROR: NNUM4TM GT MNUM4TM    ')
       IF(NNUM4TM.GT.MNU4TMX)
     .  CALL RUNENDT('ERROR: NNUM4TM GT MNU4TMX    ')
      ENDIF
C
C**** DEFINE STARTING POSITION OF MEMORY
C
      LBYW1T=LWOR1T*8
C
      RETURN
      END