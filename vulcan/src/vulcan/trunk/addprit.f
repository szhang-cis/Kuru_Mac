      SUBROUTINE ADDPRIT(IPRINT)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE POINTERS FOR PERMANENT ARRAYS
C
C
C***********************************************************************
C
C     Notes:
C
C     # variables subject to special memory conditions
C     !m mechanical variable
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
      INCLUDE 'nuef_om.f'         ! thermal-flow
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION IPRINT(50)
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
      NMEMA11=1-NMEMO11
C
      NITERCM=0
      IF(ITERME.GT.0) THEN
       IF(NITERC.EQ.2.OR.NITERC.EQ.3) NITERCM=1
      ENDIF
C
C**** COMPUTES NPREAT (PREAST(NPREAT): pressure + porosity criteria)
C
      NPRE1T=ITERME2*NMEMA11
      NPREAT=NPRE1T+NPOROT
C
      KPORETB=1
      IF(KPORET-1.GT.0) KPORETB=2          ! smoothing of water pressure
      KPORETC=0
      IF(KPORET.NE.0) KPORETC=1
C
      IEGX1=0
      IF(ITERMEF.GT.0) IEGX1=1   ! ktem1f is not available at this level
C
      NCOMPT=1+KDYNAT+NMEMO8+IEGX1   ! DISTOT(NTOTVT,NCOMPT) must be
C                           considered in the DIMENSION of the routines;
C                                NCOMPT must be included in a COMMON !!!
C
      IPRINT( 1)=1                                  ! LNODS(NNODE,NELEM)
      IPRINT( 2)=IPRINT( 1)+(NNODET*NELEMT*ICHALT+4)/8  ! MATNO(NELEM)
      IPRINT( 3)=IPRINT( 2)+(NELEMT*ICHALT+4)/8     ! PROEL(NPREL,NGRUP)
      IPRINT( 4)=IPRINT( 3)+NPRELT*NGRUPT           ! PROPS(NPROP,NMATS)
      IPRINT( 5)=IPRINT( 4)+NPROPT*NMATST          ! COORD(NDIME,NPOIN)#
      IPRINT( 6)=IPRINT( 5)+NDIMET*NPOINT*NMEMO1
C                                             ! HTLOD(NHLOD,NSUBF,NFUNC)
      IPRINT( 7)=IPRINT( 6)+NHLODT*NSUBFT*NFUNCT ! IFFIX(NTOTV,2)
      IPRINT( 8)=IPRINT( 7)+(NTOTVT*ICHALT+4)/8*
     .                         (1+KPORETB/2)*(1+NACTIT) ! PRESC(NTOTV,2)
      IPRINT( 9)=IPRINT( 8)+NTOTVT*2             ! RLOAD(NTOTV)
      IPRINT(10)=IPRINT( 9)+NTOTVT               ! RLOAH(NTOTV,NFUNC)
      IPRINT(11)=IPRINT(10)+NTOTVT*NFUNCT        ! FICTO(NFUNC)
      IPRINT(12)=IPRINT(11)+NFUNCT               ! TFICT(NFUNC)
C
      IPRINT(13)=IPRINT(12)+NFUNCT               ! DISIT(NTOTV,2)
      IPRINT(14)=IPRINT(13)+NTOTVT*(1+NMEMO9)    ! DISPR(NTOTV,3)
      IPRINT(15)=IPRINT(14)+NTOTVT*(1+KDYNAT+NMEMO8)    ! DISTO(NTOTV,3)
      IPRINT(16)=IPRINT(15)+NTOTVT*NCOMPT        ! HEADS(NPOIN,4)
      IPRINT(17)=IPRINT(16)+NPOINT*4*KPORETC     ! REFOR(NTOTV,2)
      IPRINT(18)=IPRINT(17)+NTOTVT*2             ! TLOAD(NTOTV,2)
C
      IPRINT(19)=IPRINT(18)+NTOTVT*2             ! LPNTN(NPOIN)
      IPRINT(20)=IPRINT(19)+(NPOINT*ICHALT+4)/8*KRENUT    ! ELDAT(NDATA)
      IPRINT(21)=IPRINT(20)+NDATAT               ! ELPRE(NPREV)
      IPRINT(22)=IPRINT(21)+NPREVT               ! ELVAR(NSTAT)
      IPRINT(23)=IPRINT(22)+NSTATT               ! ELMAT(NMATX)
      IPRINT(24)=IPRINT(23)+NMATXT               ! DISPLT(NTOTV) !m
C
      IPRINT(25)=IPRINT(24)+NTOTVM*ITERME2       ! PWORKT(NPOIN,3)
      IPRINT(26)=IPRINT(25)+NPOINT*(2+NITERCM)*ITERME2 ! PREAST(NPREA)
      IPRINT(27)=IPRINT(26)+NPOINT*NPREAT              ! TGAPST(NPOIN)
      IPRINT(28)=IPRINT(27)+NPOINT*ITERME2*NMEMA11     ! TEMPIT(NPOIN,2)
      IPRINT(29)=IPRINT(28)+NPOINT*(1+NMEMO10)   ! ADVELT(NTOTV) !m
      IPRINT(30)=IPRINT(29)+NTOTVM*ICONVT        ! FPCHAT(NFPCH,NPOIN)
      IPRINT(31)=IPRINT(30)+NFPCH*NPOINT         ! LACTIT(NELEM)
C
      LPRINT=IPRINT(31)+NACTIT*NELEMT
C
C**** MEMORY CONTROL
C
      IF(LPRINT.GT.MPRINT)
     . CALL RUNENDT('ERROR IN ADDPRIT: LPRINT          ')
C
C**** DEFINE STARTING POSITION OF MEMORY
C
      LBYPRT=LPRINT*8
C              
      RETURN
      END
