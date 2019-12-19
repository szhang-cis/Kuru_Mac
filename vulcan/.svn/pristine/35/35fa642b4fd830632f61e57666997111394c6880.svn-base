      SUBROUTINE ADDPRI(IPRIN)
C***********************************************************************
C
C**** THIS ROUTINE CALCULATES THE POINTERS FOR PERMANENT ARRAYS
C
C     Notes: in order to print PREAS & TGAPS, the dimensions of them
C            should be multipled by:
C            ITERME2:  for bidirectional coupled problems
C            ITERME1:  for unidirectional coupled problems
C            ITERME1A: for uncoupled problems
C
C            VNORM is only necessary for node to node contact elements
C            (ITYPE=4)
C
C***********************************************************************
C
C**** GENERAL DIMENSIONS
C
      INCLUDE 'para_om.f'
C
C**** ADDITIONAL PARAMETERS
C
      INCLUDE 'addi_om.f'
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'         ! thermal-mechanical
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'auxl_om.f'
C
      DIMENSION IPRIN(50)                    ! see note in vul-*m*.f
C
      ITERME1A=1                             ! uncoupled
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
C
      NITERCM=0
      IF(ITERME.GT.0) THEN
       IF(NITERC.EQ.2.OR.NITERC.EQ.3) NITERCM=1
      ENDIF
C
      KPOREB=1
      IF(KPORE-1.GT.0) KPOREB=2            ! smoothing of water pressure
      KPOREC=0
      IF(KPORE.NE.0) KPOREC=1
C
      NSKE1=0
      IF(NSKEW.GT.0) NSKE1=1
C
      IF(ITERME1.EQ.0) NFPCH=0
C
      IPRIN( 1)= 1                            ! LNODS(NNODE,NELEM)
      IPRIN( 2)=IPRIN( 1)+                    ! MATNO(NELEM)
     .          (NNODE*NELEM*ICHAL+4)/8
      IPRIN( 3)=IPRIN( 2)+(NELEM*ICHAL+4)/8   ! PROEL(NPREL,NGRUP)
      IPRIN( 4)=IPRIN( 3)+NPREL*NGRUP         ! PROPS(NPROP,NMATS)
      IPRIN( 5)=IPRIN( 4)+NPROP*NMATS         ! COORD(NDIME,NPOIN) 
      IPRIN( 6)=IPRIN( 5)+NDIME*NPOIN*NMEMO1M ! HTLOD(NHLOD,NSUBF,NFUNC)
C
      IPRIN( 7)=IPRIN( 6)+NHLOD*NSUBF*NFUNC   ! IFFIX(NTOTV,2)
      IPRIN( 8)=IPRIN( 7)+(NTOTV*ICHAL+4)/8*  ! PRESC(NTOTV,2)
     .          (1+KPOREB/2)*(1+NACTI)
      IPRIN( 9)=IPRIN( 8)+NTOTV*2             ! RLOAD(NTOTV)
      IPRIN(10)=IPRIN( 9)+NTOTV               ! RLOAH(NTOTV,NFUNC)
      IPRIN(11)=IPRIN(10)+NTOTV*NFUNC         ! FICTO(NFUNC)
      IPRIN(12)=IPRIN(11)+NFUNC               ! TFICT(NFUNC)
C
      IPRIN(13)=IPRIN(12)+NFUNC               ! DISIT(NTOTV,2)
      IPRIN(14)=IPRIN(13)+NTOTV*(1+NMEMO9M)   ! DISPR(NTOTV,NDISR)
      IPRIN(15)=IPRIN(14)+NTOTV*NDISR         ! DISTO(NTOTV,NDISO)
      IPRIN(16)=IPRIN(15)+NTOTV*NDISO         ! HEADS(NPOIN,4)
      IPRIN(17)=IPRIN(16)+NPOIN*4*KPOREC      ! REFOR(NTOTV,2)
      IPRIN(18)=IPRIN(17)+NTOTV*2             ! TLOAD(NTOTV,2)
C
      IPRIN(19)=IPRIN(18)+NTOTV*2             ! LPNTN(NPOIN)
      IPRIN(20)=IPRIN(19)+                    ! ELDAT(NDATA)
     .          (NPOIN*ICHAL+4)/8*KRENU
      IPRIN(21)=IPRIN(20)+NDATA               ! ELPRE(NPREV)
      IPRIN(22)=IPRIN(21)+NPREV               ! ELVAR(NSTAT)
      IPRIN(23)=IPRIN(22)+NSTAT               ! ELMAT(NMATX)
      IPRIN(24)=IPRIN(23)+NMATX               ! TEMPN(NPOIN,2)
C
      IPRIN(25)=IPRIN(24)+NPOIN*2*ITERME1     ! DTEMP(NPOIN)
      IPRIN(26)=IPRIN(25)+NPOIN*ITERME1       ! INFRI(NPOIN)
      IPRIN(27)=IPRIN(26)+                    ! COFRI(NSKEW,NDIME,NDIME)
     .          (NPOIN*ICHAL+4)/8*NSKE1
      IPRIN(28)=IPRIN(27)+NSKEW*NDIME*NDIME   ! PWORK(NPOIN,2)
      IPRIN(29)=IPRIN(28)+                    ! PREAS(NPREA,NPOIN)
     .          NPOIN*(1+NITERCM)*ITERME2
      IPRIN(30)=IPRIN(29)+                    ! TGAPS(NPOIN)
     .          NPREA*NPOIN*ITERME1A
      IPRIN(31)=IPRIN(30)+NPOIN*ITERME1A      ! VNORM(NTOTV)
      IPRIN(32)=IPRIN(31)+NTOTV               ! FPCHA(NFPCH,NPOIN)
C
      IPRIN(33)=IPRIN(32)+NFPCH*NPOIN         ! LACTI(NELEM)
      IPRIN(34)=IPRIN(33)+                    ! NOPRF(NNODE,NELEM)
     .          (NELEM*ICHAL+4)/8*NACTI
      IPRIN(35)=IPRIN(34)+                    ! PRESF(NNODE,NDIME,NELEM)
     .          (NNODE*NELEM*ICHAL+4)/8*NLDSF
      IPRIN(36)=IPRIN(35)+              ! PREHF(NNODE,NDIME,NELEM,NFUNC)
     .          NNODE*NDIME*NELEM*NLDSF
      IPRIN(37)=IPRIN(36)+                    ! VANIS(NANIV,NANIC,NELEM)
     .          NNODE*NDIME*NELEM*NFUNC*NLDSF
C
      LPRIN=IPRIN(37)+NANIV*NANIC*NELEM*NANIS
C
C**** MEMORY CONTROL
C
c     IF(LPRIN.GT.MPRIN)
c    . CALL RUNEND('ERROR IN ADDPRI: LPRIN            ')
C
C**** OTHER CONTROLS
C
c     IF(NANIS.GT.MANIS)
c    . CALL RUNEND('ERROR IN ADDPRI: NANIS            ')
c     IF(NANIV.GT.MANIV)
c    . CALL RUNEND('ERROR IN ADDPRI: NANIV            ')
c     IF(NANIC.GT.MANIC)
c    . CALL RUNEND('ERROR IN ADDPRI: NANIC            ')
C
C**** DEFINE STARTING POSITION OF MEMORY
C
      LBYPR=LPRIN*8
C              
      RETURN
      END
