      SUBROUTINE ADDPRIS(IPRINS)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE POINTERS FOR PERMANENT ARRAYS
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
      INCLUDE 'para_oms.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_oms.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'         ! thermal-mechanical
      INCLUDE 'nued_om.f'         ! thermal-microstructural
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION IPRINS(50)
C
      ITERME1=1                              ! unidirectional coupled
      ITERME2=1                              ! bidirectional coupled
      IF(ITERME.LT.0) THEN                   ! uncoupled problems
       ITERME1=0
       ITERME2=0
      ENDIF
      IF(ITERME.EQ.0) THEN                   ! unidirectional coupled
       ITERME1=1
       ITERME2=0
      ENDIF
      NMEMA11S=1-NMEMO11S
C
      NITERCM=0
      IF(ITERME.GT.0) THEN
       IF(NITERCS.EQ.2.OR.NITERCS.EQ.3) NITERCM=1
      ENDIF
C
C**** COMPUTES NPREAS (PREASS(NPREAS): temperature gradient)
C
      NPRE1S=NDIMES
      NPREAS=NPRE1S
C
      KPORESB=1
      IF(KPORES-1.GT.0) KPORESB=2          ! smoothing of water pressure
      KPORESC=0
      IF(KPORES.NE.0) KPORESC=1
C
      IF(IEVFI.EQ.0) THEN
C
       IPRINS( 1)=1                                ! LNODS(NNODE,NELEM)
       IPRINS( 2)=IPRINS( 1)                             ! MATNO(NELEM)
       IPRINS( 3)=IPRINS( 2)                       ! PROEL(NPREL,NGRUP)
       IPRINS( 4)=IPRINS( 3)                       ! PROPS(NPROP,NMATS)
       IPRINS( 5)=IPRINS( 4)+NPROPS*NMATSS         ! COORD(NDIME,NPOIN)#
       IPRINS( 6)=IPRINS( 5)                      
C                                            ! HTLOD(NHLOD,NSUBF,NFUNC)
       IPRINS( 7)=IPRINS( 6)                       ! IFFIX(NTOTV,2)
       IPRINS( 8)=IPRINS( 7)                       ! PRESC(NTOTV,2)
       IPRINS( 9)=IPRINS( 8)                       ! RLOAD(NTOTV)
       IPRINS(10)=IPRINS( 9)                       ! RLOAH(NTOTV,NFUNC)
       IPRINS(11)=IPRINS(10)                       ! FICTO(NFUNC)
       IPRINS(12)=IPRINS(11)                       ! TFICT(NFUNC)
C
       IPRINS(13)=IPRINS(12)                       ! DISIT(NTOTV,2)
       IPRINS(14)=IPRINS(13)                       ! DISPR(NTOTV,3)
       IPRINS(15)=IPRINS(14)                       ! DISTO(NTOTV,3)
       IPRINS(16)=IPRINS(15)                       ! HEADS(NPOIN,4)
       IPRINS(17)=IPRINS(16)                       ! REFOR(NTOTV,2)
C
       IPRINS(18)=IPRINS(17)                       ! TLOAD(NTOTV,2)
       IPRINS(19)=IPRINS(18)                       ! LPNTN(NPOIN)
       IPRINS(20)=IPRINS(19)                       ! ELDAT(NDATA)
       IPRINS(21)=IPRINS(20)                       ! ELPRE(NPREV)
       IPRINS(22)=IPRINS(21)                       ! ELVAR(NSTAT)
       IPRINS(23)=IPRINS(22)                       ! ELMAT(NMATX)
       IPRINS(24)=IPRINS(23)                       ! DISPLT(NTOTV) !m
C
       IPRINS(25)=IPRINS(24)                       ! PWORKT(NPOIN,3)
       IPRINS(26)=IPRINS(25)                       ! PREAST(NPOIN)
       IPRINS(27)=IPRINS(26)                       ! TGAPST(NPOIN)
       IPRINS(28)=IPRINS(27)                       ! TEMPIT(NPOIN,2)
       IPRINS(29)=IPRINS(28)                       ! ADVELT(NTOTV) !m
       IPRINS(30)=IPRINS(29)                       ! FPCHAT(NFPCH,NPOIN)
       IPRINS(31)=IPRINS(30)                       ! LACTI(NELEM)
C
       LPRINS=IPRINS(31)
C
      ELSE              ! ievfi=1
C
       IPRINS( 1)=1                                ! LNODS(NNODE,NELEM)
       IPRINS( 2)=IPRINS( 1)+(NNODES*NELEMS*ICHALS+4)/8  ! MATNO(NELEM)
       IPRINS( 3)=IPRINS( 2)+(NELEMS*ICHALS+4)/8   ! PROEL(NPREL,NGRUP)
       IPRINS( 4)=IPRINS( 3)+NPRELS*NGRUPS         ! PROPS(NPROP,NMATS)
       IPRINS( 5)=IPRINS( 4)+NPROPS*NMATSS         ! COORD(NDIME,NPOIN)#
       IPRINS( 6)=IPRINS( 5)+NDIMES*NPOINS*NMEMO1S
C                                            ! HTLOD(NHLOD,NSUBF,NFUNC)
       IPRINS( 7)=IPRINS( 6)+NHLODS*NSUBFS*NFUNCS  ! IFFIX(NTOTV,2)
       IPRINS( 8)=IPRINS( 7)+(NTOTVS*ICHALS+4)/8*
     .                         (1+KPORESB/2)*(1+NACTIS) ! PRESC(NTOTV,2)
       IPRINS( 9)=IPRINS( 8)+NTOTVS*2              ! RLOAD(NTOTV)
       IPRINS(10)=IPRINS( 9)+NTOTVS                ! RLOAH(NTOTV,NFUNC)
       IPRINS(11)=IPRINS(10)+NTOTVS*NFUNCS         ! FICTO(NFUNC)
       IPRINS(12)=IPRINS(11)+NFUNCS                ! TFICT(NFUNC)
C
       IPRINS(13)=IPRINS(12)+NFUNCS                ! DISIT(NTOTV,2)
       IPRINS(14)=IPRINS(13)+NTOTVS*(1+NMEMO9S)    ! DISPR(NTOTV,3)
       IPRINS(15)=IPRINS(14)+NTOTVS*(1+KDYNAS+NMEMO8S)  ! DISTO(NTOTV,3)
       IPRINS(16)=IPRINS(15)+NTOTVS*(1+KDYNAS+NMEMO8S)  ! HEADS(NPOIN,4)
       IPRINS(17)=IPRINS(16)+NPOINS*4*KPORESC           ! REFOR(NTOTV,2)
       IPRINS(18)=IPRINS(17)+NTOTVS*2              ! TLOAD(NTOTV,2)
C
       IPRINS(19)=IPRINS(18)+NTOTVS*2              ! LPNTN(NPOIN)
       IPRINS(20)=IPRINS(19)+(NPOINS*ICHALS+4)/8*KRENUS   ! ELDAT(NDATA)
       IPRINS(21)=IPRINS(20)+NDATAS                ! ELPRE(NPREV)
       IPRINS(22)=IPRINS(21)+NPREVS                ! ELVAR(NSTAT)
       IPRINS(23)=IPRINS(22)+NSTATS                ! ELMAT(NMATX)
       IPRINS(24)=IPRINS(23)+NMATXS                ! DISPLT(NTOTV) !m
C
       IPRINS(25)=IPRINS(24)+NTOTVMS*ITERME2       ! PWORKT(NPOIN,3)
       IPRINS(26)=IPRINS(25)+NPOINS*(2+NITERCM)*ITERME2  ! PREAST(NPOIN)
       IPRINS(27)=IPRINS(26)+NPOINS*NPREAS               ! TGAPST(NPOIN)
       IPRINS(28)=IPRINS(27)+NPOINS*ITERME2*NMEMA11S   ! TEMPIT(NPOIN,2)
       IPRINS(29)=IPRINS(28)+NPOINS*(1+NMEMO10S)   ! ADVELT(NTOTV) !m
       IPRINS(30)=IPRINS(29)+NTOTVMS*ICONVS        ! FPCHAT(NFPCH,NPOIN)
       IPRINS(31)=IPRINS(30)+NFPCH*NPOINS          ! LACTIT(NELEM)
C
       LPRINS=IPRINS(31)+NACTIS*NELEMS
C
      ENDIF               ! ievfi.eq.0
C
C**** MEMORY CONTROL
C
      IF(LPRINS.GT.MPRINS)
     . CALL RUNENDS('ERROR IN ADDPRIS: LPRINS          ')
C
C**** DEFINE STARTING POSITION OF MEMORY
C
      LBYPRS=LPRINS*8
C              
      RETURN
      END
