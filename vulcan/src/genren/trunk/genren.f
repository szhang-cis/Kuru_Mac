C***********************************************************************
C
C     MASTER GENTOP (GENREN)
C
C**** PROGRAM FOR THE AUTOMATIC GENERATION OF FINITE ELEMENT MESHES
C
C
C     WARNINGS: 1) CHECK MAXIMUM DIMENSIONS IN genren_i.f
C               2) ONLY QUADRILATERAL ELEMENTS ARE AVAILABLE IN 2D &
C                  BRICK ELEMENTS IN 3D
C               3) MPOIN is normally greater than the NPOIN generated
C                  (see genren_i.f)
C
C     LINUX VERSION:
C
C     compile: forv-t genren.f (O2, cl)
C     execute: genren.c
C
C
C     W95 VERSION:
C
C     1) mv genren.f genren.for
C     2) mv genren_i.f genren_i.for
C     3) comment "linux lines" & uncomment "w95 lines" in genren.for
C
C     compile: w95 facilities
C     execute: genren name      (genren means genren.bat)
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
c     INCLUDE 'FLIB.FD'                                    ! w95
C
      CHARACTER*80 C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12, ! linux & w95
     .             C13,C14,C15,C16
C
c     CHARACTER*6 VARIABLE                                 ! w95
c     INTEGER*4 LVAR                                       ! w95
C
C**** THIS SUBROUTINE ACCEPTS DATA DEFINING THE STRUCTURAL OUTLINE
C
C     DEFAULTS VALUES:
C
      NOPTI=0         ! vulcan-t, vulcan-f or vulcan-tf, vulcan-m,
C                     ! vulcan-tm
      NOPTJ=1         ! thermal elements 101 (1=no, 2=yes)
      NOPTK=1         ! velocity boundary conditions (1=no, 2=yes)
      NOPTKK=0        ! for NOPTK=2:
C                     ! =0 (11 in 2D; 111 in 3D)
C                     ! =1 (10 in 2D; 100 in 3D)
C                     ! =2 (01 in 2D; 010 in 3D)
C                     ! =3 (        ; 001 in 3D)
      NOPTL=0         ! for NOPTJ=2, convection conditions + prescribed
C                     ! normal heat flux (minimum material number)
      NOPTM=1         ! maximum material number for fluid
      NOPTN=0         ! for NOPTJ=2, prescribed normal heat flux
C                     ! (minimum material number)
      NOPTO=1         ! thermal elements 104 (1=no, 2=yes)
      NOPTP=2         ! for NOPTO=2, truncation in the coordinates
C                     ! NOPTP gives the digits considered in the coord.
      NOPTQ=1         ! for NOPTI=1, advective file (1=no, 2=yes)
C                     ! for NOPTI=3, NOPTQ=2 always
      NOPTR=1         ! for NOPTI=4,5, mechanical initial conditions
C                     ! (1=no, 2=yes)
      NOPTS=1         ! for NOPTO=2, indicator for coincident or non
C                     ! coincident (nc) nodes for thermomechanical
C                     ! contact (1=no, 2=yes)
      NOPTT=0         ! for NOPTS=2, thermomechanical contact (minimum
C                     ! material number)
      NOPTU=1         ! 1: no material front; 2: material front
C
      NSKEW=1         ! skew system, only for fluids (1=no, 2=yes)
C
      NACTI=0         ! active elements option
C
      NTRAN=1         ! transition elements (1: no; 2:yes)
C
      NFEMV=0         ! output for postprocess (0=no, 1=FEMVIEW,
C                     ! 2=FLAVIA, 3=PATRAN, 4=GiD)
C
C-----------------------------------------------------------------------
C
C     Notes:
C
C     Typical element numbering: first:  fluid materials
C                                second: additional thermal materials
C                                third:  thermal elements 104
C                                fourth: thermal elements 101 with
C                                        convection conditions
C                                fifth:  thermal elements 101 with
C                                        prescribed heat flux
C
C     NOPTL must be less or equal than NOPTN when NOPTN > 0
C
C     For non coincident nodes at the thermomechanical contact
C     (NOPTS=2), NOPTT must be less than NOPTL when NOPTL > 0
C
C-----------------------------------------------------------------------
C
      write(6,99)
   99 format(' Quiere una salida para VULCAN-T : 1',/,
     .       ' Quiere una salida para VULCAN-F : 2',/,
     .       ' Quiere una salida para VULCAN-TF: 3',/,
     .       ' Quiere una salida para VULCAN-M : 4',/,
     .       ' Quiere una salida para VULCAN-TM: 5'   )
      read (5,*) NOPTI
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       write(6,93)
   93  format(' No hay elementos de contacto termomec. entre mat.: 1',/,
     .        ' Hay elementos de contacto termomec. entre mat.:    2')
       read (5,*) NOPTO
C
       IF(NOPTO.EQ.2) THEN
        write(6,92)
   92   format(' Numero de digitos decimales en las coordenadas',
     . ' (disponibles: 2,3,4 o 6)')
        read (5,*) NOPTP
        if(NOPTP.LT.2.OR.NOPTP.GT.6) then
         WRITE(8,*) '**** NUMERO DE DIGITOS ERRONEO ****'
         STOP
        endif
C
        write(6,87)
   87   format(' Nodos coincidentes en el contacto:    1',/,
     .         ' Nodos NO coincidentes en el contacto: 2')
        read (5,*) NOPTS
C
        IF(NOPTS.EQ.2) THEN
         write(6,86)
   86    format(' Material/es de contacto termomecanico con nodos',
     . ' no coincidentes (menor numeracion)')
         read (5,*) NOPTT
        ENDIF
       ENDIF
      ENDIF
C
      IF(NOPTI.EQ.4) THEN
       write(6,89)
   89  format(' No hay bloques de contorno termico (a eliminar): 1',/,
     .        ' Hay bloques de contorno termico (a eliminar):    2')
       read (5,*) NOPTJ
C
       IF(NOPTJ.EQ.2) THEN
        write(6,88)
   88   format(' Material/es de contorno termico a eliminar',
     . ' (menor numeracion)')
        read (5,*) NOPTL
       ENDIF
      ENDIF
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       write(6,98)
   98  format(' No asigna flujo de calor normal en los contornos: 1',/,
     .        ' Asigna flujo de calor normal en los contornos:    2')
       read (5,*) NOPTJ
C
       IF(NOPTJ.EQ.2) THEN
        write(6,96)
   96   format(' Material/es de contorno con flujo de calor normal',
     . ' [conveccion y valor prescripto] (menor numeracion)')
        read (5,*) NOPTL
C
        write(6,94)
   94   format(' Material/es de contorno con flujo de calor normal',
     . ' prescripto (menor numeracion)')
        read (5,*) NOPTN
       ENDIF
      ENDIF
C
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) THEN
       write(6,97)
   97  format(' No asigna condiciones de contorno de velocidad: 1',/,
     .        ' Asigna condiciones de contorno de velocidad:    2')
       read (5,*) NOPTK
C
       IF(NOPTK.EQ.2) THEN
        write(6,101)
  101   format(' =0 (11 in 2D; 111 in 3D)',/,
     .         ' =1 (10 in 2D; 100 in 3D)',/,
     .         ' =2 (01 in 2D; 010 in 3D)',/,
     .         ' =3 (        ; 001 in 3D)')
        read (5,*) NOPTKK
       ENDIF
C
       write(6,100)
  100  format(' No se consideran "skew systems": 1',/,
     .        ' Se consideran "skew systems":    2')
       read (5,*) nskew
C
       write(6,102)
  102  format(' No hay frente material: 1',/,
     .        ' Hay frente material:    2')
       read (5,*) noptu
      ENDIF
C
      IF(NOPTI.EQ.3) THEN
       write(6,95)
   95  format(' Material/es del fluido (mayor numeracion): ')
       read (5,*) NOPTM
      ENDIF
C
      IF(NOPTI.EQ.1) THEN
       write(6,91)
   91  format(' No hay efectos advectivos: 1',/,
     .        ' Hay efectos advectivos:    2')
       read (5,*) NOPTQ
      ENDIF
C
      IF(NOPTI.EQ.3) THEN
       NOPTQ=2
      ENDIF
C
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       write(6,90)
   90  format(' No hay efectos inerciales: 1',/,
     .        ' Hay efectos inerciales:    2')
       read (5,*) NOPTR
      ENDIF
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       write(6,103)
  103  format(' No se considera la opcion ACTIVE_ELEMENTS: 1',/,
     .        ' Se considera la opcion ACTIVE_ELEMENTS:    2')
       read (5,*) NACTI
       IF(NACTI.EQ.2) THEN
        IF(NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
         write(6,*) 'Warning: only one file is created (*-t.act)'
        ENDIF
       ENDIF
      ENDIF
C
      write(6,104)
  104 format(' No hay elementos de transicion: 1',/,
     .       ' Hay elementos de transicion:    2')
      read (5,*) NTRAN
C
      write(6,85)
   85 format(' No quiere una salida para postproceso: 0',/,
     .       ' Quiere una salida para postproceso:',/,
     .       '    1=FEMVIEW, 2=FLAVIA, 3=PATRAN, 4=GiD')
      read (5,*) nfemv
C
C**** CONTROLS
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       IF(NOPTJ.EQ.2) THEN
        IF(NOPTN.GT.0) THEN
         IF(NOPTN.LT.NOPTL) THEN
          WRITE(8,*) '**** NOPTN.LT.NOPTL ****'
          STOP
         ENDIF
        ENDIF
       ENDIF
C
       IF(NOPTS.EQ.2) THEN
        IF(NOPTT.EQ.0) THEN
         WRITE(8,*) '**** NOPTT=0 FOR NOPTS=2 ****'
         STOP
        ENDIF
C
        IF(NOPTL.GT.0) THEN
         IF(NOPTL.LE.NOPTT) THEN
          WRITE(8,*) '**** NOPTL.LT.NOPTT ****'
          STOP
         ENDIF
        ENDIF
       ENDIF
      ENDIF
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       IF(NOPTJ.EQ.2) THEN
        IF(NOPTL.EQ.0) THEN
         WRITE(8,*) '**** NOPTL=0 FOR NOPTJ=2 ****'
         STOP
        ENDIF
       ENDIF
      ENDIF
C
C**** FILE ASSIGMENT
C
      call getenv ('FOR007',c1)                ! sg, sun, convex & linux
      call getenv ('FOR008',c2)                ! sg, sun, convex & linux
C
c     VARIABLE='FOR007'                        ! w95
c     LVAR=GETENVQQ (VARIABLE,C1)              ! w95
c     VARIABLE='FOR008'                        ! w95
c     LVAR=GETENVQQ (VARIABLE,C2)              ! w95
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       call getenv ('FOR003',c3)               ! sg, sun, convex & linux
       call getenv ('FOR004',c4)               ! sg, sun, convex & linux
       call getenv ('FOR009',c5)               ! sg, sun, convex & linux
       IF(NOPTQ.EQ.2) THEN
        call getenv ('FOR010',c6)              ! sg, sun, convex & linux
       ENDIF
C
c      VARIABLE='FOR003'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C3)             ! w95
c      VARIABLE='FOR004'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C4)             ! w95
c      VARIABLE='FOR009'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C5)             ! w95
c      VARIABLE='FOR010'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C6)             ! w95
C
       IF(NOPTJ.EQ.2) THEN
        IF(NOPTN.GT.0) THEN
         call getenv ('FOR013',c9)             ! sg, sun, convex & linux
C
c        VARIABLE='FOR013'                     ! w95
c        LVAR=GETENVQQ (VARIABLE,C9)           ! w95
        ENDIF
       ENDIF
      ENDIF
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) THEN
       call getenv ('FOR011',c7)               ! sg, sun, convex & linux
       call getenv ('FOR012',c8)               ! sg, sun, convex & linux
C
c      VARIABLE='FOR011'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C7)             ! w95
c      VARIABLE='FOR012'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C8)             ! w95
      ENDIF
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       call getenv ('FOR014',c10)              ! sg, sun, convex & linux
       call getenv ('FOR015',c11)              ! sg, sun, convex & linux
       call getenv ('FOR016',c12)              ! sg, sun, convex & linux
C
c      VARIABLE='FOR014'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C10)            ! w95
c      VARIABLE='FOR015'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C11)            ! w95
c      VARIABLE='FOR016'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C12)            ! w95
      ENDIF
      IF(NFEMV.EQ.1) THEN                      ! FEMVIEW
       WRITE(8,*) '**** POSTPROCESS FOR FEMVIEW NOT IMPLEMENTED ****'
       STOP
       call getenv ('FOR017',c13)              ! sg, sun, convex & linux
c      VARIABLE='FOR017'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C13)            ! w95
      ENDIF
      IF(NFEMV.EQ.2) THEN                      ! FLAVIA
       WRITE(8,*) '**** POSTPROCESS FOR FLAVIA NOT IMPLEMENTED ****'
       STOP
      ENDIF
      IF(NFEMV.EQ.3) THEN                      ! PATRAN
       WRITE(8,*) '**** POSTPROCESS FOR PATRAN NOT IMPLEMENTED ****'
       STOP
      ENDIF
      IF(NFEMV.EQ.4) THEN                      ! GiD
       call getenv ('FOR018',c14)              ! sg, sun, convex & linux
c      VARIABLE='FOR018'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C14)            ! w95
       call getenv ('FOR019',c15)              ! sg, sun, convex & linux
c      VARIABLE='FOR019'                       ! w95
c      LVAR=GETENVQQ (VARIABLE,C15)            ! w95
      ENDIF
      IF(NACTI.EQ.2) THEN
       call getenv ('FOR020',c16)              ! sg, sun, convex & linux
      ENDIF
C
      OPEN(7,FILE=C1,STATUS='OLD')
      OPEN(8,FILE=C2,STATUS='UNKNOWN')
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       OPEN(3,FILE=C3,STATUS='UNKNOWN')
       OPEN(4,FILE=C4,STATUS='UNKNOWN')
       OPEN(9,FILE=C5,STATUS='UNKNOWN')
       IF(NOPTQ.EQ.2) OPEN(10,FILE=C6,STATUS='UNKNOWN')
       IF(NOPTJ.EQ.2) THEN
        IF(NOPTN.GT.0) THEN
         OPEN(13,FILE=C9,STATUS='UNKNOWN')
        ENDIF
       ENDIF
      ENDIF
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) THEN
       OPEN(11,FILE=C7,STATUS='UNKNOWN')
       OPEN(12,FILE=C8,STATUS='UNKNOWN')
      ENDIF
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       OPEN(14,FILE=C10,STATUS='UNKNOWN')
       OPEN(15,FILE=C11,STATUS='UNKNOWN')
       IF(NOPTR.EQ.2) OPEN(16,FILE=C12,STATUS='UNKNOWN')
      ENDIF
      IF(NFEMV.EQ.2) OPEN(17,FILE=C13,STATUS='UNKNOWN')
      IF(NFEMV.EQ.4) THEN
                     OPEN(18,FILE=C14,STATUS='UNKNOWN')
                     OPEN(19,FILE=C15,STATUS='UNKNOWN')
      ENDIF
      IF(NACTI.EQ.2) THEN
                     OPEN(20,FILE=C16,STATUS='UNKNOWN')
      ENDIF
C
      LNODE=9
      NITER=1
C
      NP=7                      ! input: .dat
      NT=8                      ! output: .sal
C
C**** GEOMETRY OUTPUT FILE (.geo): CONNECTIVITIES & COORDINATES
C
      NR=3                      ! thermal
      NRF=11                    ! flow
      NRM=14                    ! mechanical
C
C**** BOUNDARY CONDITIONS OUTPUT FILE (.fix)
C
      NS=4                      ! thermal
      NSF=12                    ! flow
      NSM=15                    ! mechanical
C
C**** INITIAL CONDITIONS OUTPUT FILE (.ini .adv): TEMPERATURE, VELOCITY
C     OR DISPLACEMENT/VELOCITY
C
      NI=9                      ! thermal
      NA=10                     ! flow
      NIM=16                    ! mechanical
C
C**** NORMAL HEAT FLUXES OUTPUT FILE (.loa)
C
      NG=13                     ! prescribed normal heat flux
C
C**** FEMVIEW FILE
C
      NFEM=17
C
C**** ACTIVE ELEMENTS FILE
C
      NACT=20
C
C**** THIS ROUTINE ACCEPTS MOST OF THE INPUT DATA
C
      CALL INBLOC
C
C**** THIS SUBROUTINE UNDERTAKES THE MESH SUBDIVISION
C
      CALL GENER
C
C**** THIS SUBROUTINE REDUCES THE FRONTWIDTH
C
C     IF(NOPTI.NE.1) CALL FRONTWIDTH
C
C**** THIS SUBROUTINE OUTPUTS THE GENERATED MESH
C
      CALL GMESH(NOPTI,NOPTJ,NOPTK,NOPTKK,NOPTL,NOPTM,NOPTN,NOPTO,
     .           NOPTP,NOPTQ,NOPTR,NOPTS, NOPTT,NOPTU,NSKEW,NACTI,
     .           NFEMV)
C
      STOP
      END
C_______________________________________________________________________
      SUBROUTINE INBLOC
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
      CHARACTER TITLE*80
C
C**** THIS SUBROUTINE ACCEPTS THE INPUT DATA DEFINING THE STRUCTURAL
C     OUTLINE
C
c     WRITE(NR,960)                  ! not necessary
c 960 FORMAT(' START')
C
      READ(NP,900) TITLE
c     WRITE(NR,900) TITLE            ! not necessary
  900 FORMAT(1A80)
C
      READ(NP,*) NPOIN,NELEM,NDIME,NNODB,NTYPE,TOLCO
      WRITE(NT,910) NPOIN,NELEM
  905 FORMAT(15I5)
  910 FORMAT(1H0,5X,30HAUTOMATIC DATA GENERATION FROM,I4,21H  LOCATION P
     .OINTS AND,I3,13H  DATA BLOCKS)
      WRITE(NT,945) NNODB
  945 FORMAT(1H0,5X,'NUMBER OF NODES PER BLOCK   =',I3)
      WRITE(NT,920)
  920 FORMAT(1H0,5X,11HDATA BLOCKS)
      WRITE(NT,925)
  925 FORMAT(1H0,9HBLOCK NO.,7X,17HDEFINITION POINTS,11X,8HMATERIAL)
C
C**** CONTROLS DIMENSIONS (INPUT)
C
      IF(NPOIN.GT.MPOIN) THEN
       WRITE(NT,*) '*** NPOIN.GT.MPOIN ***',NPOIN,MPOIN
       STOP
      ENDIF
      IF(NDIME.GT.MDIME) THEN
       WRITE(NT,*) '*** NDIME.GT.MDIME ***',NDIME,MDIME
       STOP
      ENDIF
      IF(NELEM.GT.MBLOC) THEN
       WRITE(NT,*) '*** NELEM.GT.MBLOC ***',NELEM,MBLOC
       STOP
      ENDIF
      IF(NNODB.GT.MNODM) THEN
       WRITE(NT,*) '*** NNODB.GT.MNODM ***',NNODB,MNODM
       STOP
      ENDIF
C
      DO 10 IELEM=1,NELEM
      READ(NP,*) NUMEL,NNODE(NUMEL),NINTE(NUMEL),NGAUS(NUMEL),
     .    MATNO(NUMEL),(LNODS(NUMEL,INODE),INODE=1,NNODB)
   10 WRITE(NT,905) NUMEL,NNODE(NUMEL),NINTE(NUMEL),NGAUS(NUMEL),
     .    MATNO(NUMEL),(LNODS(NUMEL,INODE),INODE=1,NNODB)
      WRITE(NT,930)
      WRITE(NT,935)
  930 FORMAT(1H0,5X,15HLOCATION POINTS)
  935 FORMAT(1H0,9HPOINT NO.,5X,10H  X-COORD.,5X,10H  Y-COORD.,
     .5X,10H  Z-COORD.)
        DO 55 IPOIN=1,NPOIN
        DO 55 IDIME=1,NDIME
   55 COORD(IPOIN,IDIME)=55555555555.555555
      I=0
   33 READ(NP,*) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
      I=I+1
      IF(JPOIN.LT.NPOIN) GO TO 33
      IF(I.EQ.NPOIN) GO TO 40
C
C**** SERENDIPITY & LAGRANGEAN ELEMENTS (ONLY 2D)
C
      IF(NDIME.EQ.2) THEN
       NNODG=NNODB
       IF(NNODB.EQ.9) NNODB=8
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NNODB.GT.MNODM) THEN
        WRITE(NT,*) '*** NNODB.GT.MNODM ***',NNODB,MNODM
        STOP
       ENDIF
C
       DO 50 IELEM=1,NELEM
       DO 50 INODE=1,NNODB
       IPOIN=LNODS(IELEM,INODE)
       IF(COORD(IPOIN,1).NE.55555555555.555555) GO TO 50
       JNODE=INODE+1
       KNODE=INODE-1
       IF(INODE.EQ.8) JNODE=1
       JPOIN=LNODS(IELEM,JNODE)
       KPOIN=LNODS(IELEM,KNODE)
       DO 60 IDIME=1,NDIME
   60  COORD(IPOIN,IDIME)=(COORD(JPOIN,IDIME)+COORD(KPOIN,IDIME))/2.
   50  CONTINUE
       NNODB=NNODG
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NNODB.GT.MNODM) THEN
        WRITE(NT,*) '*** NNODB.GT.MNODM ***',NNODB,MNODM
        STOP
       ENDIF
C
       IF(NNODB.EQ.9) THEN
       DO 70 IELEM=1,NELEM
       LPOIN=LNODS(IELEM,9)
       IF(COORD(LPOIN,1).EQ.5555.55) THEN
       EXISP=0.0
       ETASP=0.0
       DO 5 IDIME=1,NDIME
    5  COORD(LPOIN,IDIME)=0.0
       NNODG=NNODB
       NNODB=8
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NNODB.GT.MNODM) THEN
        WRITE(NT,*) '*** NNODB.GT.MNODM ***',NNODB,MNODM
        STOP
       ENDIF
C
       CALL SFRQ(EXISP,ETASP)
       NNODB=NNODG
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NNODB.GT.MNODM) THEN
        WRITE(NT,*) '*** NNODB.GT.MNODM ***',NNODB,MNODM
        STOP
       ENDIF
C
       DO 80 INODE=1,8
       JTEMP=LNODS(IELEM,INODE)
       DO 80 IDIME=1,NDIME
   80  COORD(LPOIN,IDIME)=COORD(LPOIN,IDIME)+SHAPE(INODE)*
     .                    COORD(JTEMP,IDIME)
       END IF
   70  CONTINUE
       END IF
      ENDIF             ! ndime.eq.2
C
      IF(NDIME.EQ.3) THEN
       WRITE(NT,*) '*** 20 & 27 NODED ELEMENTS NOT IMPLEMENTED ***'
       STOP
      ENDIF             ! ndime.eq.3
C
   40 CONTINUE
      DO 20 IPOIN=1,NPOIN
   20 WRITE(NT,940) IPOIN,(COORD(IPOIN,IDIME),IDIME=1,NDIME)
 9999 FORMAT(I5,3F10.3)
  940 FORMAT(I10,3F15.5)
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE GENER
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
C**** THIS SUBROUTINE PERFORMS THE MESH SUBDIVISION
C
C**** INITIALIZATION SECTION
C
      DO 10 IREPN=1,MREPN
   10 LREPN(IREPN)=0
      NPONT=NPOIN
      NBLOC=NELEM
      NPOIN=0
      NELEM=0
C
      DO 20 IPONT=1,NPONT
      DO 20 IDIME=1,NDIME
   20 TCORD(IPONT,IDIME)=COORD(IPONT,IDIME)
      DO 30 IPOIN=1,MPOIN
      DO 30 IDIME=1,NDIME
   30 COORD(IPOIN,IDIME)=0.0
      DO 40 IBLOC=1,NBLOC
      MMGAU(IBLOC)=NGAUS(IBLOC)
      MMINT(IBLOC)=NINTE(IBLOC)
      MMATO(IBLOC)=MATNO(IBLOC)
      NBNOD(IBLOC)=NNODE(IBLOC)
      DO 40 INODE=1,NNODB
      MLNOD(IBLOC,INODE)=LNODS(IBLOC,INODE)
   40 LNODS(IBLOC,INODE)=0
C
C***  READ AND WRITE BLOCK SUBDIVISION DATA
C
      DO 170 IBLOC=1,NBLOC
      DO 25 IJK=1,MDIVX
       WEITX(IJK)=0.0
   25 CONTINUE
      DO 26 IJK=1,MDIVY
       WEITY(IJK)=0.0
   26 CONTINUE
      IF(NDIME.EQ.3) THEN
       DO 27 IJK=1,MDIVZ
        WEITZ(IJK)=0.0
   27  CONTINUE
      ENDIF
C
      IF(NDIME.EQ.2) THEN
       MNODE=4
       IF(NBNOD(IBLOC).EQ.3) MNODE=4
       IF(NBNOD(IBLOC).EQ.6.OR.NBNOD(IBLOC).EQ.8.OR.
     .    NBNOD(IBLOC).EQ.9) MNODE=8
       KNODE=MNODE/4
       FNODE=KNODE
C
       IF(NBNOD(IBLOC).EQ.3.OR.NBNOD(IBLOC).EQ.4) 
     .    LNODE=4
       IF(NBNOD(IBLOC).EQ.8) LNODE=8
      ENDIF
C
      IF(NDIME.EQ.3) THEN
       MNODE=8
       KNODE=MNODE/8
       FNODE=KNODE
C
       LNODE=8
      ENDIF
C
      IF(NTRAN.EQ.1) THEN                  ! wihtout transition elements
       IF(NDIME.EQ.2) READ(NP,*) KBLOC,NDIVX,NDIVY
       IF(NDIME.EQ.3) READ(NP,*) KBLOC,NDIVX,NDIVY,NDIVZ
      ELSE                                 ! with transition elements
       IF(NDIME.EQ.2) READ(NP,*) KBLOC,NDIVX,NDIVY,KTRAN
       IF(NDIME.EQ.3) READ(NP,*) KBLOC,NDIVX,NDIVY,NDIVZ,KTRAN
C
       IF(KTRAN.EQ.1) THEN                 ! 3 to 1 transition
        IF(((NDIVX/3)*3-NDIVX).NE.0) THEN
         WRITE(8,*) '**** NRO. DE DIV. EN X NO ES MULTIPLO DE 3 ****'
         STOP
        ENDIF
        IF(NDIVY.NE.2) THEN
         WRITE(8,*) '**** NRO. DE DIV. EN Y NO ES IGUAL A 2 ****'
         STOP
        ENDIF
        KTRAT=3
       ENDIF
       IF(KTRAN.EQ.2) THEN                 ! 4 to 2 transition
        IF(((NDIVX/4)*4-NDIVX).NE.0) THEN
         WRITE(8,*) '**** NRO. DE DIV. EN X NO ES MULTIPLO DE 4 ****'
         STOP
        ENDIF
        IF(NDIVY.NE.2) THEN
         WRITE(8,*) '**** NRO. DE DIV. EN Y NO ES IGUAL A 2 ****'
         STOP
        ENDIF
        KTRAT=4
       ENDIF
C
       MMDIV(KBLOC)=NDIVX
       NELEM3=0
       NELEM3A=0
      ENDIF
C
C**** CONTROLS DIMENSIONS (INPUT)
C
      IF(NDIVX.GT.MDIVX) THEN
       WRITE(NT,*) '*** NDIVX.GT.MDIVX ***',NDIVX,MDIVX
       STOP
      ENDIF
      IF(NDIVY.GT.MDIVY) THEN
       WRITE(NT,*) '*** NDIVY.GT.MDIVY ***',NDIVY,MDIVY
       STOP
      ENDIF
      IF(NDIME.EQ.3) THEN
       IF(NDIVZ.GT.MDIVZ) THEN
        WRITE(NT,*) '*** NDIVZ.GT.MDIVZ ***',NDIVZ,MDIVZ
        STOP
       ENDIF
      ENDIF
C
      READ(NP,*) (WEITX(IDIVX),IDIVX=1,NDIVX)
      READ(NP,*) (WEITY(IDIVY),IDIVY=1,NDIVY)
      IF(NDIME.EQ.3) READ(NP,*) (WEITZ(IDIVZ),IDIVZ=1,NDIVZ)
      WRITE(NT,910) KBLOC
      WRITE(NT,915)
      IF(NDIME.EQ.2) WRITE(NT,920) NDIVX,NDIVY
      IF(NDIME.EQ.3) WRITE(NT,921) NDIVX,NDIVY,NDIVZ
      WRITE(NT,925)
      WRITE(NT,905) (WEITX(IDIVX),IDIVX=1,NDIVX)
      WRITE(NT,930)
      WRITE(NT,905) (WEITY(IDIVY),IDIVY=1,NDIVY)
      IF(NDIME.EQ.3) THEN
       WRITE(NT,931) 
       WRITE(NT,905) (WEITZ(IDIVZ),IDIVZ=1,NDIVZ)
      ENDIF
C
  900 FORMAT(3I5)
  905 FORMAT(8F10.3)
 9998 FORMAT(20F5.2)
  910 FORMAT(1H0,5X,14HDATA BLOCK NO.,I3)
  915 FORMAT(6X,17H-----------------)
  920 FORMAT(1H0,5X,32HNO. OF DIV. IN FIRST DIRECTION =,I3,2X,
     .33HNO. OF DIV. IN SECOND DIRECTION =,I3)
  921 FORMAT(1H0,5X,32HNO. OF DIV. IN FIRST DIRECTION =,I3,2X,
     .33HNO. OF DIV. IN SECOND DIRECTION =,I3,2X,
     .32HNO. OF DIV. IN THIRD DIRECTION =,I3)
  925 FORMAT(1HO,5X,44HLIST OF WEIGHTING FACTORS IN FIRST DIRECTION)
  930 FORMAT(1H0,5X,45HLIST OF WEIGHTING FACTORS IN SECOND DIRECTION)
  931 FORMAT(1H0,5X,44HLIST OF WEIGHTING FACTORS IN THIRD DIRECTION)
C
C**** DIVIDE EACH BLOCK INTO ELEMENTS
C
      TOTAL=0.0
      DO 50 IDIVX=1,NDIVX
      IF(WEITX(IDIVX).EQ.0.0) WEITX(IDIVX)=1.0
   50 TOTAL=TOTAL+WEITX(IDIVX)
      XNORM=2.0/TOTAL
      TOTAL=0.0
      DO 60 IDIVY=1,NDIVY
      IF(WEITY(IDIVY).EQ.0.0) WEITY(IDIVY)=1.0
   60 TOTAL=TOTAL+WEITY(IDIVY)
      YNORM=2.0/TOTAL
      IF(NDIME.EQ.3) THEN
       TOTAL=0.0
       DO 61 IDIVZ=1,NDIVZ
       IF(WEITZ(IDIVZ).EQ.0.0) WEITZ(IDIVZ)=1.0
   61  TOTAL=TOTAL+WEITZ(IDIVZ)
       ZNORM=2.0/TOTAL
      ENDIF
C
      NXTWO=NDIVX*KNODE+1
      NYTWO=NDIVY*KNODE+1
      IF(NDIME.EQ.3) NZTWO=NDIVZ*KNODE+1
C
      IF(NDIME.EQ.2) THEN
       IASEY=0
       ETASP=-1.0
       KWETY=0
       KONTY=-1
       DO 160 IYTWO=1,NYTWO
       IASEY=IASEY+1
       IF(MNODE.NE.8.AND.IASEY.EQ.3) IASEY=2
       IF(MNODE.EQ.8.AND.IASEY.EQ.4) IASEY=2
       IASEX=0
       EXISP=-1.0
       KWETX=0
       KONTX=-1
       DO 130 IXTWO=1,NXTWO
       IASEX=IASEX+1
       IF(MNODE.NE.8.AND.IASEX.EQ.3) IASEX=2
       IF(MNODE.EQ.8.AND.IASEX.EQ.4) IASEX=2
       IF(IASEX.EQ.2.AND.IASEY.EQ.2.AND.NBNOD(IBLOC).EQ.8) GO TO 100
       NPOIN=NPOIN+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NPOIN.GT.MPOIN) THEN
        WRITE(NT,*) '*** NPOIN.GT.MPOIN ***',NPOIN,MPOIN
        STOP
       ENDIF
C
       CALL SFRQ(EXISP,ETASP)
       DO 70 INODE=1,NNODB
       JTEMP=MLNOD(IBLOC,INODE)
       DO 70 IDIME=1,NDIME
   70  COORD(NPOIN,IDIME)=COORD(NPOIN,IDIME)+
     .                    SHAPE(INODE)*TCORD(JTEMP,IDIME)
       GO TO(75,75,75,80,80,85,85,90,95)NBNOD(IBLOC)
   75  IF(IASEX.NE.2.OR.IASEY.NE.2) GO TO 100
       JPOIN=NPOIN-NXTWO
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-1
       LNODS(NELEM,2)=JPOIN
       LNODS(NELEM,3)=NPOIN
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-1
       LNODS(NELEM,2)=NPOIN
       LNODS(NELEM,3)=NPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 100
   80  IF(IASEX.NE.2.OR.IASEY.NE.2) GO TO 100
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       JPOIN=NPOIN-NXTWO
       LNODS(NELEM,1)=JPOIN-1
       LNODS(NELEM,2)=JPOIN
       LNODS(NELEM,3)=NPOIN
       LNODS(NELEM,4)=NPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       IF(NTRAN.EQ.2) THEN                       ! deals with transition
        LTRAN(NELEM,1)=NELEM
        LTRAN(NELEM,2)=0
        LTRAN(NELEM,3)=0
        LTRAN(NELEM,4)=0
        IF(KTRAN.NE.0) THEN
         NELEM3=NELEM3+1
         NELEM3B=(NELEM3-1)/MMDIV(IBLOC)+1
         LTRAN(NELEM,2)=NELEM3B                  ! 1=1st row; 2=2nd row
C
         NELEM3A=NELEM3A+1
         I3A=((NELEM3A-1)/KTRAT)*KTRAT+1
         IF((NELEM3A-I3A).EQ.0) NELEM3A=1
         LTRAN(NELEM,3)=NELEM3A                  ! (1,2,3) or (1,2,3,4)
C
         LTRAN(NELEM,4)=MMDIV(IBLOC)
         LTRAN(NELEM,5)=KTRAN                    ! both transitions
        ENDIF
       ENDIF
       GO TO 100
   85  IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 100
       IPOIN=NPOIN-NXTWO
       JPOIN=IPOIN-NXTWO
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=IPOIN-1
       LNODS(NELEM,3)=NPOIN
       LNODS(NELEM,4)=NPOIN-1
       LNODS(NELEM,5)=NPOIN-2
       LNODS(NELEM,6)=IPOIN-2
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 100
   90  IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 100
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       IPOIN=NPOIN-IXTWO-NDIVX+(IXTWO-1)/2
       JPOIN=NPOIN-NXTWO-NDIVX-1
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=NPOIN-1
       LNODS(NELEM,7)=NPOIN-2
       LNODS(NELEM,8)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 100
   95  IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 100
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       IPOIN=NPOIN-NXTWO
       JPOIN=IPOIN-NXTWO
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=NPOIN-1
       LNODS(NELEM,7)=NPOIN-2
       LNODS(NELEM,8)=IPOIN-2
       LNODS(NELEM,9)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
  100  CONTINUE
       GO TO (110,110,110,110,110,120,120,120,120)NBNOD(IBLOC)
  110  KWETX=KWETX+1
       GO TO 130
  120  IF(KONTX.LT.0) KWETX=KWETX+1
       KONTX=KONTX*(-1)
  130  EXISP=EXISP+XNORM*WEITX(KWETX)/FNODE
       GO TO (140,140,140,140,140,150,150,150,150)NBNOD(IBLOC)
  140  KWETY=KWETY+1
       GO TO 160
  150  IF(KONTY.LT.0) KWETY=KWETY+1
       KONTY=KONTY*(-1)
  160  ETASP=ETASP+YNORM*WEITY(KWETY)/FNODE
      ENDIF                ! ndime.eq.2
C
      IF(NDIME.EQ.3) THEN
       IASEZ=0
       EZETA=-1.0
       KWETZ=0
       KONTZ=-1
       DO 560 IZTWO=1,NZTWO
       IASEZ=IASEZ+1
       IF(MNODE.NE.20.AND.IASEZ.EQ.3) IASEZ=2
       IF(MNODE.EQ.20.AND.IASEZ.EQ.4) IASEZ=2
       IASEY=0
       ETASP=-1.0
       KWETY=0
       KONTY=-1
       DO 460 IYTWO=1,NYTWO
       IASEY=IASEY+1
       IF(MNODE.NE.20.AND.IASEY.EQ.3) IASEY=2
       IF(MNODE.EQ.20.AND.IASEY.EQ.4) IASEY=2
       IASEX=0
       EXISP=-1.0
       KWETX=0
       KONTX=-1
       DO 430 IXTWO=1,NXTWO
       IASEX=IASEX+1
       IF(MNODE.NE.20.AND.IASEX.EQ.3) IASEX=2
       IF(MNODE.EQ.20.AND.IASEX.EQ.4) IASEX=2
       IF(IASEX.EQ.2.AND.IASEY.EQ.2.AND.IASEZ.EQ.2.AND.
     .                                     NBNOD(IBLOC).EQ.20) GO TO 500
       NPOIN=NPOIN+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NPOIN.GT.MPOIN) THEN
        WRITE(NT,*) '*** NPOIN.GT.MPOIN ***',NPOIN,MPOIN
        STOP
       ENDIF
C
       CALL SFRQ3D(EXISP,ETASP,EZETA)
       DO 470 INODE=1,NNODB
       JTEMP=MLNOD(IBLOC,INODE)
       DO 470 IDIME=1,NDIME
  470  COORD(NPOIN,IDIME)=COORD(NPOIN,IDIME)+
     .                    SHAPE(INODE)*TCORD(JTEMP,IDIME)
       GO TO(475,475,475,475,475,475,475,
     .       480,480,485,485,490,495) NBNOD(IBLOC)
  475  WRITE(NT,*) '*** 475 NOT IMPLEMENTED ***'
       STOP
       IF(IASEX.NE.2.OR.IASEY.NE.2) GO TO 500      ! not implemented yet
       JPOIN=NPOIN-NXTWO
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-1
       LNODS(NELEM,2)=JPOIN
       LNODS(NELEM,3)=NPOIN
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-1
       LNODS(NELEM,2)=NPOIN
       LNODS(NELEM,3)=NPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 500
  480  IF(IASEX.NE.2.OR.IASEY.NE.2.OR.IASEZ.NE.2) GO TO 500
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       JPOIN=NPOIN-NXTWO*NYTWO
       KPOIN=JPOIN-NXTWO
       LNODS(NELEM,1)=KPOIN-1
       LNODS(NELEM,2)=KPOIN
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=JPOIN-1
       JPOIN=NPOIN
       KPOIN=JPOIN-NXTWO
       LNODS(NELEM,5)=KPOIN-1
       LNODS(NELEM,6)=KPOIN
       LNODS(NELEM,7)=JPOIN
       LNODS(NELEM,8)=JPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       IF(NTRAN.EQ.2) THEN                       ! deals with transition
        LTRAN(NELEM,1)=NELEM
        LTRAN(NELEM,2)=0
        LTRAN(NELEM,3)=0
        LTRAN(NELEM,4)=0
        IF(KTRAN.NE.0) THEN
         NELEM3=NELEM3+1
         NLAYZ=(NELEM3-1)/(2*MMDIV(IBLOC))+1     ! layers in z
         NELEM3B=((NELEM3-1)/MMDIV(IBLOC)+1)/NLAYZ
         LTRAN(NELEM,2)=NELEM3B                  ! 1=1st row; 2=2nd row
C
         NELEM3A=NELEM3A+1
         I3A=((NELEM3A-1)/KTRAT)*KTRAT+1
         IF((NELEM3A-I3A).EQ.0) NELEM3A=1
         LTRAN(NELEM,3)=NELEM3A                  ! (1,2,3) or (1,2,3,4)
C
         LTRAN(NELEM,4)=MMDIV(IBLOC)
         LTRAN(NELEM,5)=KTRAN                    ! both transitions
        ENDIF
       ENDIF
       GO TO 500
  485  WRITE(NT,*) '*** 485 NOT IMPLEMENTED ***'
       STOP
       IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 500      ! not implemented yet
       IPOIN=NPOIN-NXTWO
       JPOIN=IPOIN-NXTWO
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=IPOIN-1
       LNODS(NELEM,3)=NPOIN
       LNODS(NELEM,4)=NPOIN-1
       LNODS(NELEM,5)=NPOIN-2
       LNODS(NELEM,6)=IPOIN-2
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 500
  490  WRITE(NT,*) '*** 490 NOT IMPLEMENTED ***'
       STOP
       IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 500      ! not implemented yet
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       IPOIN=NPOIN-IXTWO-NDIVX+(IXTWO-1)/2
       JPOIN=NPOIN-NXTWO-NDIVX-1
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=NPOIN-1
       LNODS(NELEM,7)=NPOIN-2
       LNODS(NELEM,8)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
       GO TO 500
  495  WRITE(NT,*) '*** 495 NOT IMPLEMENTED ***'
       STOP
       IF(IASEX.NE.3.OR.IASEY.NE.3) GO TO 500      ! not implemented yet
       NELEM=NELEM+1
C
C**** CONTROLS DIMENSIONS (GENERATED)
C
       IF(NELEM.GT.MELEM) THEN
        WRITE(NT,*) '*** NELEM.GT.MELEM ***',NELEM,MELEM
        STOP
       ENDIF
C
       IPOIN=NPOIN-NXTWO
       JPOIN=IPOIN-NXTWO
       LNODS(NELEM,1)=JPOIN-2
       LNODS(NELEM,2)=JPOIN-1
       LNODS(NELEM,3)=JPOIN
       LNODS(NELEM,4)=IPOIN
       LNODS(NELEM,5)=NPOIN
       LNODS(NELEM,6)=NPOIN-1
       LNODS(NELEM,7)=NPOIN-2
       LNODS(NELEM,8)=IPOIN-2
       LNODS(NELEM,9)=IPOIN-1
       NGAUS(NELEM)=MMGAU(IBLOC)
       NINTE(NELEM)=MMINT(IBLOC)
       MATNO(NELEM)=MMATO(IBLOC)
       NNODE(NELEM)=NBNOD(IBLOC)
  500  CONTINUE
       GO TO (410,410,410,410,410,410,410,410,
     .        420,420,420,420) NBNOD(IBLOC)
  410  KWETX=KWETX+1
       GO TO 430
  420  IF(KONTX.LT.0) KWETX=KWETX+1
       KONTX=KONTX*(-1)
  430  EXISP=EXISP+XNORM*WEITX(KWETX)/FNODE
       GO TO (440,440,440,440,440,440,440,440,
     .        450,450,450,450) NBNOD(IBLOC)
  440  KWETY=KWETY+1
       GO TO 460
  450  IF(KONTY.LT.0) KWETY=KWETY+1
       KONTY=KONTY*(-1)
  460  ETASP=ETASP+YNORM*WEITY(KWETY)/FNODE
       GO TO (540,540,540,540,540,540,540,540,
     .        550,550,550,550) NBNOD(IBLOC)
  540  KWETZ=KWETZ+1
       GO TO 560
  550  IF(KONTZ.LT.0) KWETZ=KWETZ+1
       KONTZ=KONTZ*(-1)
  560  EZETA=EZETA+ZNORM*WEITZ(KWETZ)/FNODE
      ENDIF                ! ndime.eq.3
  170 CONTINUE
C
C**** ELIMINATE REPEATED NODES AT BLOCK INTERFACES
C
      NREPN=0
      DO 210 IPOIN=1,NPOIN
      IF(NREPN.EQ.0) GO TO 190
      DO 180 IREPN=1,NREPN
      IF(IPOIN.EQ.LREPN(IREPN)) GO TO 210
  180 CONTINUE
  190 CONTINUE
      LPOIN=IPOIN+1
      DO 200 JPONT=LPOIN,NPOIN
      IF(NDIME.EQ.2) THEN
       TOTAL=ABS(COORD(IPOIN,1)-COORD(JPONT,1))+
     .       ABS(COORD(IPOIN,2)-COORD(JPONT,2))
      ENDIF
      IF(NDIME.EQ.3) THEN
       TOTAL=ABS(COORD(IPOIN,1)-COORD(JPONT,1))+
     .       ABS(COORD(IPOIN,2)-COORD(JPONT,2))+
     .       ABS(COORD(IPOIN,3)-COORD(JPONT,3))
      ENDIF
      IF(TOTAL.GT.TOLCO) GO TO 200
      NREPN=NREPN+1
      LREPN(NREPN)=JPONT
      LASOC(NREPN)=IPOIN
  200 CONTINUE
  210 CONTINUE
      IF(NREPN.EQ.0) GO TO 360
      INDEX=0
      DO 240 IPOIN=1,NPOIN
      DO 220 IREPN=1,NREPN
      IF(LREPN(IREPN).EQ.IPOIN) GO TO 230
  220 CONTINUE
      GO TO 240
  230 INDEX=INDEX+1
      LFINN(INDEX)=LREPN(IREPN)
      LFASC(INDEX)=LASOC(IREPN)
  240 CONTINUE
      DO 250 IREPN=1,NREPN
      LREPN(IREPN)=LFINN(IREPN)
  250 LASOC(IREPN)=LFASC(IREPN)
      DO 260 IREPN=1,NREPN
      DO 260 IELEM=1,NELEM
      DO 260 INODE=1,NNODE(IELEM)
  260 IF(LNODS(IELEM,INODE).EQ.LREPN(IREPN))
     . LNODS(IELEM,INODE)=LASOC(IREPN)
      DO 310 IPOIN=1,NPOIN
      DO 270 IREPN=1,NREPN
      IF(IPOIN.EQ.LREPN(IREPN)) GO TO 310
  270 CONTINUE
      IF(IPOIN.LT.LREPN(1)) GO TO 310
      IDIFF=IPOIN-NREPN
      IF(IPOIN.GT.LREPN(NREPN)) GO TO 290
      DO 280 IREPN=1,NREPN
      KREPN=NREPN-IREPN+1
  280 IF(IPOIN.LT.LREPN(KREPN)) IDIFF=IPOIN-KREPN+1
  290 DO 300 IDIME=1,NDIME
  300 COORD(IDIFF,IDIME)=COORD(IPOIN,IDIME)
  310 CONTINUE
      DO 350 IELEM=1,NELEM
      DO 350 INODE=1,NNODE(IELEM)
      NPOSI=LNODS(IELEM,INODE)
      DO 320 IREPN=1,NREPN
      IF(NPOSI.EQ.LREPN(IREPN)) GO TO 350
  320 CONTINUE
      IF(NPOSI.LT.LREPN(1)) GO TO 350
      IDIFF=NPOSI-NREPN
      IF(NPOSI.GT.LREPN(NREPN)) GO TO 340
      DO 330 IREPN=1,NREPN
      KREPN=NREPN-IREPN+1
  330 IF(NPOSI.LT.LREPN(KREPN)) IDIFF=NPOSI-KREPN+1
  340 LNODS(IELEM,INODE)=IDIFF
  350 CONTINUE
  360 CONTINUE
      NPOIN=NPOIN-NREPN
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE GMESH(NOPTI,NOPTJ,NOPTK,NOPTKK,NOPTL,NOPTM,NOPTN,NOPTO,
     .                 NOPTP,NOPTQ,NOPTR,NOPTS, NOPTT,NOPTU,NSKEW,NACTI,
     .                 NFEMV)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
C**** THIS SUBROUTINE OUTPUTS THE COORDINATES AND ELEMENT TOPOLOGIES OF
C     THE GENERATED MESH & OTHER THINGS !
C
c     WRITE(NT,900)
c     WRITE(NT,905) NPOIN                 ! generated points by genren
c     WRITE(NT,910) NELEM                 ! generated elements by genren
  900 FORMAT(////10X,29H**** GENERATED MESH DATA ****//)
  905 FORMAT(1H0,10X,24HNUMBER OF NODAL POINTS =,I7)
  910 FORMAT(1H0,14X,20HNUMBER OF ELEMENTS =,I7)
  930 FORMAT(1H0,10X,25HELEMENT NODAL CONNECTIONS)
  935 FORMAT(1H0,9H  ELEMENT,14X,7HMAT NO.,11X,12HNODAL POINTS)
  940 FORMAT(13I5)
  915 FORMAT(1H0,4X,23HNODAL POINT COORDINATES)
  920 FORMAT(1H0,9H NODE NO.,10H  X-COORD.,10H  Y-COORD.,
     .10H  Z-COORD.)
  925 FORMAT(I5,3F10.3)
C
C**** DEALS WITH TRANSITION ELEMENTS
C
C     2D & 3D: 3 to 1 (MMINT(IBLOC)=3)
C              4 to 2 (MMINT(IBLOC)=4)
C
      IF(NTRAN.EQ.2) THEN                        ! deals with transition
       DO IPOIN=1,NPOIN
        LPUNT(IPOIN)=IPOIN                       ! initialization
        DO IDIME=1,NDIME
         COORT(IPOIN,IDIME)=COORD(IPOIN,IDIME)
        ENDDO
       ENDDO
C
       DO IELEM=1,NELEM
        DO INODE=1,NNODE(IELEM)
         LNODT(IELEM,INODE)=LNODS(IELEM,INODE)
        ENDDO
       ENDDO
C
       DO IELEM=1,NELEM
        IF(LTRAN(IELEM,5).EQ.1) THEN        ! 3 to 1 transition
         IF(LTRAN(IELEM,2).EQ.1) THEN       ! 1st row transition element
          IF(LTRAN(IELEM,3).EQ.1) THEN
           LPUNT(LNODS(IELEM,4))=0
           LNODS(IELEM,4)=LNODS(IELEM+LTRAN(IELEM,4),4)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,8))=0
            LNODS(IELEM,8)=LNODS(IELEM+LTRAN(IELEM,4),8)
           ENDIF
          ENDIF
          IF(LTRAN(IELEM,3).EQ.3) THEN
           LPUNT(LNODS(IELEM,3))=0
           LNODS(IELEM,3)=LNODS(IELEM+LTRAN(IELEM,4),3)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,7))=0
            LNODS(IELEM,7)=LNODS(IELEM+LTRAN(IELEM,4),7)
           ENDIF
          ENDIF
         ENDIF
         IF(LTRAN(IELEM,2).EQ.2) THEN       ! 2nd row transition element
          IF(LTRAN(IELEM,3).EQ.1.OR.LTRAN(IELEM,3).EQ.3) THEN
           LTRAN(IELEM,1)=0
          ENDIF
          IF(LTRAN(IELEM,3).EQ.2) THEN
           LPUNT(LNODS(IELEM,3))=0
           LPUNT(LNODS(IELEM,4))=0
           LNODS(IELEM,3)=LNODS(IELEM-LTRAN(IELEM,4)+1,3)
           LNODS(IELEM,4)=LNODS(IELEM-LTRAN(IELEM,4)-1,4)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,7))=0
            LPUNT(LNODS(IELEM,8))=0
            LNODS(IELEM,7)=LNODS(IELEM-LTRAN(IELEM,4)+1,7)
            LNODS(IELEM,8)=LNODS(IELEM-LTRAN(IELEM,4)-1,8)
           ENDIF
          ENDIF
         ENDIF
        ENDIF
        IF(LTRAN(IELEM,5).EQ.2) THEN        ! 4 to 2 transition
         IF(LTRAN(IELEM,2).EQ.1) THEN       ! 1st row transition element
          IF(LTRAN(IELEM,3).EQ.1) THEN
           LPUNT(LNODS(IELEM,4))=0
           LNODS(IELEM,4)=LNODS(IELEM+LTRAN(IELEM,4),4)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,8))=0
            LNODS(IELEM,8)=LNODS(IELEM+LTRAN(IELEM,4),8)
           ENDIF
          ENDIF
          IF(LTRAN(IELEM,3).EQ.4) THEN
           LPUNT(LNODS(IELEM,3))=0
           LNODS(IELEM,3)=LNODS(IELEM+LTRAN(IELEM,4),3)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,7))=0
            LNODS(IELEM,7)=LNODS(IELEM+LTRAN(IELEM,4),7)
           ENDIF
          ENDIF
         ENDIF
         IF(LTRAN(IELEM,2).EQ.2) THEN       ! 2nd row transition element
          IF(LTRAN(IELEM,3).EQ.1.OR.LTRAN(IELEM,3).EQ.4) THEN
           LTRAN(IELEM,1)=0
          ENDIF
          IF(LTRAN(IELEM,3).EQ.2) THEN
           LPUNT(LNODS(IELEM,4))=0
           LNODS(IELEM,4)=LNODS(IELEM-LTRAN(IELEM,4)-1,4)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,8))=0
            LNODS(IELEM,8)=LNODS(IELEM-LTRAN(IELEM,4)-1,8)
           ENDIF
          ENDIF
          IF(LTRAN(IELEM,3).EQ.3) THEN
           LPUNT(LNODS(IELEM,3))=0
           LNODS(IELEM,3)=LNODS(IELEM-LTRAN(IELEM,4)+1,3)
           IF(NDIME.EQ.3) THEN
            LPUNT(LNODS(IELEM,7))=0
            LNODS(IELEM,7)=LNODS(IELEM-LTRAN(IELEM,4)+1,7)
           ENDIF
          ENDIF
         ENDIF
        ENDIF
       ENDDO                     ! ielem=1,nelem
C
       NPOIN3=0                  ! remove fictitious transition nodes
       DO IPOIN=1,NPOIN
        IF(LPUNT(IPOIN).NE.0) THEN
         NPOIN3=NPOIN3+1
         DO IDIME=1,NDIME
          COORD(NPOIN3,IDIME)=COORT(IPOIN,IDIME)
         ENDDO
        ENDIF
       ENDDO
       NPOIN=NPOIN3
C
       NELEM3C=0                 ! remove fictitious transition elements
       DO IELEM=1,NELEM
        IF(LTRAN(IELEM,1).NE.0) THEN
         NELEM3C=NELEM3C+1
         DO INODE=1,NNODE(IELEM)
          LNODT(NELEM3C,INODE)=LNODS(IELEM,INODE)
         ENDDO
         MATNT(NELEM3C)=MATNO(IELEM)
         NINTT(NELEM3C)=NINTE(IELEM)
         NGAUT(NELEM3C)=NGAUS(IELEM)
        ENDIF
       ENDDO
       NELEM=NELEM3C
       DO IELEM=1,NELEM
        DO INODE=1,NNODE(IELEM)
         LNODS(IELEM,INODE)=LNODT(IELEM,INODE)
        ENDDO
        MATNO(IELEM)=MATNT(IELEM)
        NINTE(IELEM)=NINTT(IELEM)
        NGAUS(IELEM)=NGAUT(IELEM)
       ENDDO
       DO IELEM=1,NELEM
        DO INODE=1,NNODE(IELEM)
         IX1=LNODS(IELEM,INODE)
         IX3=0
         DO IX2=1,IX1
          IF(LPUNT(IX2).EQ.0) IX3=IX3+1
         ENDDO
         LNODS(IELEM,INODE)=LNODS(IELEM,INODE)-IX3
        ENDDO
       ENDDO
      ENDIF               ! ntran.eq.2
C
C**** INITIALIZATION
C
      DO IPOIN=1,NPOIN
       JPUNT(IPOIN)=IPOIN   ! useful for thermal problems (101 elements)
       KPUNT(IPOIN)=IPOIN   ! useful for flow problems (101 elements)
      ENDDO
      NODMAXFL=1            ! maximum node numbering for flow problem
      NELMAXFL=1            ! maximum element numbering for flow prob.
C
      NUMATMAX=1            ! maximum material numbering
      DO IELEM=1,NELEM
       IF(MATNO(IELEM).GT.NUMATMAX) NUMATMAX=MATNO(IELEM)
      ENDDO 
C
C**** GEOMETRY, BOUNDARY & INITIAL CONDITIONS
C
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) WRITE(NRF,*) 'ELEMENTS'
C
      DO 40 IELEM=1,NELEM
C
      KNODE=1           ! reorders connectivities for quadratic elements
      JNODE=2
      KKY=1
      KDIME=4
      IF(NDIME.EQ.3) KDIME=8
      IF(NNODE(IELEM).GT.KDIME) KKY=2
      DO INODE=1,NNODE(IELEM)
       IF(INODE.LE.KDIME) THEN
        NUEVN(INODE)=LNODS(IELEM,KNODE)
        KNODE=KNODE+KKY
       ELSE
        NUEVN(INODE)=LNODS(IELEM,JNODE)
        KNODE=KNODE+2
       ENDIF
      ENDDO
C
      NOPTJA=NOPTJ
      NOPTLA=NOPTL
      NOPTLB=0
      IF(NOPTS.EQ.2) THEN          ! reassigment for nc nodes at contact
       NOPTJA=2
       NOPTLA=NOPTT
       NOPTLB=NUMATMAX+1
       IF(NOPTL.GT.0) NOPTLB=NOPTL
      ENDIF
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       IF(NOPTJA.EQ.1) THEN               ! no thermal elements 101
        WRITE(NR,950) IELEM,MATNO(IELEM),
     .               (NUEVN(INODE),INODE=1,NNODE(IELEM))
       ELSE                               ! thermal elements 101 &/or nc
        IF(NOPTLB.EQ.0) THEN                   ! coincident mesh
         IF(MATNO(IELEM).LT.NOPTLA) THEN       ! solid + gap elements
          WRITE(NR,950) IELEM,MATNO(IELEM),
     .                 (NUEVN(INODE),INODE=1,NNODE(IELEM))
         ELSE                                  ! boundary elements
          WRITE(NR,950) IELEM,MATNO(IELEM),
     .                 (NUEVN(INODE),INODE=1,NNODE(IELEM)/2)
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).lt.noptla
        ELSE                                   ! non-coincident mesh
         IF(MATNO(IELEM).LT.NOPTLA) THEN       ! solid elements
          WRITE(NR,950) IELEM,MATNO(IELEM),
     .                 (NUEVN(INODE),INODE=1,NNODE(IELEM))
         ENDIF       ! matno(ielem).lt.noptla
         IF(MATNO(IELEM).GE.NOPTLA.AND.
     .      MATNO(IELEM).LT.NOPTLB) THEN       ! gap elements
          WRITE(NR,950) IELEM,MATNO(IELEM),
     .                 (NUEVN(INODE),INODE=1,NNODE(IELEM)/2),
     .                  NINTE(IELEM),NGAUS(IELEM)
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).ge.noptla.and.matno(ielem).lt.noptlb
         IF(MATNO(IELEM).GE.NOPTLB) THEN       ! boundary elements
          WRITE(NR,950) IELEM,MATNO(IELEM),
     .                 (NUEVN(INODE),INODE=1,NNODE(IELEM)/2)
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).ge.noptlb
        ENDIF        ! noptlb.gt.0
C
        IF(NOPTN.GT.0) THEN                ! prescribed normal heat flux
         IF(MATNO(IELEM).GE.NOPTN) THEN
          WRITE(NG,957) IELEM,1,0.0
         ENDIF       ! matno(ielem).ge.noptn
        ENDIF        ! noptn.gt.0
       ENDIF         ! noptja.eq.1
      ENDIF          ! nopti.eq.1.or.nopti.eq.3.or.nopti.eq.5
C
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) THEN
       IF(MATNO(IELEM).LE.NOPTM) THEN
        WRITE(NRF,941) IELEM,
     .                (NUEVN(INODE),INODE=1,NNODE(IELEM))
        DO INODE=1,NNODE(IELEM)          ! compute nodes of flow problem
         KPOIN=NUEVN(INODE)
         IF(KPOIN.GT.NODMAXFL) NODMAXFL=KPOIN
        ENDDO
        NELMAXFL=IELEM
       ELSE
        DO INODE=1,NNODE(IELEM)
         KPOIN=NUEVN(INODE)
         IF(KPOIN.GT.NODMAXFL) KPUNT(KPOIN)=0
        ENDDO
       ENDIF         ! matno(ielem).eq.noptm
      ENDIF          ! nopti.eq.2.or.nopti.eq.3
C
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       IF(NOPTJA.EQ.1) THEN               ! no thermal elements 101
        WRITE(NRM,950) IELEM,MATNO(IELEM),
     .                (nuevn(INODE),INODE=1,NNODE(IELEM))
       ELSE                               ! thermal elements 101 &/or nc
        IF(NOPTLB.EQ.0) THEN                   ! coincident mesh
         IF(MATNO(IELEM).LT.NOPTLA) THEN       ! solid + gap elements
          WRITE(NRM,950) IELEM,MATNO(IELEM),
     .                  (NUEVN(INODE),INODE=1,NNODE(IELEM))
         ELSE                                  ! boundary elements
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).lt.noptla
        ELSE                                   ! non-coincident mesh
         IF(MATNO(IELEM).LT.NOPTLA) THEN       ! solid elements
          WRITE(NRM,950) IELEM,MATNO(IELEM),
     .                  (NUEVN(INODE),INODE=1,NNODE(IELEM))
         ENDIF       ! matno(ielem).lt.noptla
         IF(MATNO(IELEM).GE.NOPTLA.AND.
     .      MATNO(IELEM).LT.NOPTLB) THEN       ! gap elements
          IF(NDIME.EQ.2) THEN
           WRITE(NRM,950) IELEM,MATNO(IELEM),
     .                   (NUEVN(INODE),INODE=1,NNODE(IELEM)/2),
     .                    NINTE(IELEM),NGAUS(IELEM)
          ENDIF
          IF(NDIME.EQ.3) THEN                  ! useful for auxiliar
           IDUMM=0                             ! program reorder-mesh.f
           WRITE(NRM,950) IELEM,MATNO(IELEM),
     .                   (NUEVN(INODE),INODE=1,NNODE(IELEM)/2),
     .                    NINTE(IELEM),NGAUS(IELEM),IDUMM,IDUMM
          ENDIF
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).ge.noptla.and.matno(ielem).lt.noptlb
         IF(MATNO(IELEM).GE.NOPTLB) THEN       ! boundary elements
          DO INODE=NNODE(IELEM)/2+1,NNODE(IELEM)
           JPOIN=NUEVN(INODE)
           JPUNT(JPOIN)=0
          ENDDO
         ENDIF       ! matno(ielem).ge.noptlb
        ENDIF        ! noptlb.gt.0
       ENDIF         ! noptja.eq.1
      ENDIF          ! nopti.eq.4.or.nopti.eq.5
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       IF(NACTI.EQ.2) THEN
        IACTI=1      ! all active elements are assumed
        WRITE(NACT,941) IELEM,IACTI
       ENDIF
      ENDIF          ! nopti.eq.1,3,4 or 5
C
   40 CONTINUE
C
      NPOIX=0        ! computes generated nodes for the analysis
      DO IPOIN=1,NPOIN
       IF(JPUNT(IPOIN).GT.0) NPOIX=NPOIX+1
      ENDDO
C
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) WRITE(NRF,*) 'END_ELEMENTS'
      IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) WRITE(NRF,*) 'COORDINATES'
C
      IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
       WRITE(NI,*) NPOIX
       IF(NOPTQ.EQ.2) WRITE(NA,*) NPOIX
      ENDIF
C
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
       IF(NOPTR.EQ.2) THEN
        WRITE(NIM,*) NPOIX
       ENDIF
      ENDIF
C
      DO 30 IPOIN=1,NPOIN
      JPOIN=JPUNT(IPOIN)
      KPOIN=KPUNT(IPOIN)
C
      IF(JPOIN.GT.0) THEN
       IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
        WRITE(NS,949) JPOIN,0,1,0.0
       ENDIF             ! nopti.eq.2.or.nopti.eq.3
      ENDIF              ! jpoin.gt.0
C
      IF(KPOIN.GT.0) THEN
       IF(NOPTI.EQ.2) THEN                 ! flow boundary conditions
        IF(NOPTK.EQ.1) THEN
         IF(NOPTU.EQ.1) THEN
          IF(NDIME.EQ.2) WRITE(NSF,956) KPOIN,0,0,0,0.0,0.0,0.0
          IF(NDIME.EQ.3) WRITE(NSF,958) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
         ELSE
          IF(NDIME.EQ.2) WRITE(NSF,958) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
          IF(NDIME.EQ.3) WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,
     .                                                              0.0
         ENDIF
        ELSE
         ICOUNT=0
         DO IELEM=1,NELMAXFL
          DO INODE=1,NNODE(IELEM)
           IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
          ENDDO
         ENDDO
         IF(NDIME.EQ.2) THEN
          IF(ICOUNT.LT.4) THEN             ! quadrilateral elements only
           IF(NOPTU.EQ.1) THEN
            IF(NOPTKK.EQ.0) WRITE(NSF,956) KPOIN,1,1,0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,956) KPOIN,1,0,0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,956) KPOIN,0,1,0,0.0,0.0,0.0
           ELSE
            IF(NOPTKK.EQ.0) WRITE(NSF,958) KPOIN,1,1,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,958) KPOIN,1,0,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,958) KPOIN,0,1,0,0,0.0,0.0,0.0,0.0
           ENDIF
          ELSE
           IF(NOPTU.EQ.1) THEN
            WRITE(NSF,956) KPOIN,0,0,0,0.0,0.0,0.0
           ELSE
            WRITE(NSF,958) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
           ENDIF
          ENDIF          ! icount.lt.4
         ENDIF           ! ndime.eq.2
         IF(NDIME.EQ.3) THEN
          IF(ICOUNT.LT.8) THEN             ! brick elements only 
           IF(NOPTU.EQ.1) THEN
            IF(NOPTKK.EQ.0) WRITE(NSF,958) KPOIN,1,1,1,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,958) KPOIN,1,0,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,958) KPOIN,0,1,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.3) WRITE(NSF,958) KPOIN,0,0,1,0,0.0,0.0,0.0,0.0
           ELSE
            IF(NOPTKK.EQ.0) WRITE(NSF,959) KPOIN,1,1,1,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,959) KPOIN,1,0,0,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,959) KPOIN,0,1,0,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
            IF(NOPTKK.EQ.3) WRITE(NSF,959) KPOIN,0,0,1,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
           ENDIF
          ELSE
           IF(NOPTU.EQ.1) THEN
            WRITE(NSF,958) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
           ELSE
            WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0
           ENDIF
          ENDIF          ! icount.lt.8
         ENDIF           ! ndime.eq.3
        ENDIF            ! noptk.eq.1
       ENDIF             ! nopti.eq.2
C
       IF(NOPTI.EQ.3) THEN                 ! thermal-flow boundary cond.
        IF(NOPTK.EQ.1) THEN
         IF(NOPTU.EQ.1) THEN
          IF(NDIME.EQ.2) WRITE(NSF,946) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
          IF(NDIME.EQ.3) WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,
     .                                                              0.0
         ELSE
          IF(NDIME.EQ.2) WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,
     .                                                              0.0
          IF(NDIME.EQ.3) WRITE(NSF,969) KPOIN,0,0,0,0,0,0,0.0,0.0,0.0,
     .                                                    0.0,0.0,0.0
         ENDIF
        ELSE
         ICOUNT=0
         DO IELEM=1,NELMAXFL
          DO INODE=1,NNODE(IELEM)
           IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
          ENDDO
         ENDDO
         IF(NDIME.EQ.2) THEN
          IF(ICOUNT.LT.4) THEN             ! quadrilateral elements only
           IF(NOPTU.EQ.1) THEN
            IF(NOPTKK.EQ.0) WRITE(NSF,946) KPOIN,1,1,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,946) KPOIN,1,0,0,0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,946) KPOIN,0,1,0,0,0.0,0.0,0.0,0.0
           ELSE
            IF(NOPTKK.EQ.0) WRITE(NSF,959) KPOIN,1,1,0,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,959) KPOIN,1,0,0,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,959) KPOIN,0,1,0,0,0,0.0,0.0,0.0,
     .                                                         0.0,0.0
           ENDIF
          ELSE
           IF(NOPTU.EQ.1) THEN
            WRITE(NSF,946) KPOIN,0,0,0,0,0.0,0.0,0.0,0.0
           ELSE
            WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0
           ENDIF
          ENDIF          ! icount.lt.4
         ENDIF           ! ndime.eq.2
         IF(NDIME.EQ.3) THEN
          IF(ICOUNT.LT.8) THEN             ! brick elements only
           IF(NOPTU.EQ.1) THEN
            IF(NOPTKK.EQ.0) WRITE(NSF,959)
     .                               KPOIN,1,1,1,0,0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,959)
     .                               KPOIN,1,0,0,0,0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,959)
     .                               KPOIN,0,1,0,0,0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.3) WRITE(NSF,959)
     .                               KPOIN,0,0,1,0,0,0.0,0.0,0.0,0.0,0.0
           ELSE
            IF(NOPTKK.EQ.0) WRITE(NSF,969)
     .                         KPOIN,1,1,1,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.1) WRITE(NSF,969)
     .                         KPOIN,1,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.2) WRITE(NSF,969)
     .                         KPOIN,0,1,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0
            IF(NOPTKK.EQ.3) WRITE(NSF,969)
     .                         KPOIN,0,0,1,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0
           ENDIF
          ELSE
           IF(NOPTU.EQ.1) THEN
            WRITE(NSF,959) KPOIN,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0
           ELSE
            WRITE(NSF,969) KPOIN,0,0,0,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0
           ENDIF
          ENDIF          ! icount.lt.8
         ENDIF           ! ndime.eq.2
        ENDIF            ! noptk.eq.1
       ENDIF             ! nopti.eq.3
      ENDIF              ! kpoin.gt.0
C
      IF(JPOIN.GT.0) THEN                  ! mech. & thermomech. bc
       IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
        IF(NDIME.EQ.2) WRITE(NSM,951) JPOIN,0,0,1,0.0,0.0
        IF(NDIME.EQ.3) WRITE(NSM,952) JPOIN,0,0,0,1,0.0,0.0,0.0
       ENDIF             ! nopti.eq.4.or.nopti.eq.5
      ENDIF              ! jpoin.gt.0
C
      IF(JPOIN.GT.0) THEN                  ! initial temp. & veloc.
       IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
        WRITE(NI,947) JPOIN,0.0
        IF(NOPTQ.EQ.2) THEN
         IF(NDIME.EQ.2) WRITE(NA,948) JPOIN,0.0,0.0
         IF(NDIME.EQ.3) WRITE(NA,960) JPOIN,0.0,0.0,0.0
        ENDIF
       ENDIF
C
       IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN   ! initial disp.
        IF(NOPTR.EQ.2) THEN
         IF(NDIME.EQ.2) WRITE(NIM,948) JPOIN,0.0,0.0
         IF(NDIME.EQ.3) WRITE(NIM,960) JPOIN,0.0,0.0,0.0
        ENDIF
       ENDIF
C
       IF(NOPTI.EQ.1.OR.NOPTI.EQ.3.OR.NOPTI.EQ.5) THEN
        IF(NOPTO.EQ.1) THEN
         WRITE(NR,945) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
        ELSE
         IF(NOPTP.EQ.2)
     .    WRITE(NR,942) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.3)
     .    WRITE(NR,943) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.4)
     .    WRITE(NR,944) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.6)
     .    WRITE(NR,945) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
        ENDIF            ! noptj.eq.1
       ENDIF             ! nopti.eq.1.or.nopti.eq.3.or.nopti.eq.5
C
       IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
        IF(NOPTO.EQ.1) THEN
         WRITE(NRM,945) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
        ELSE
         IF(NOPTP.EQ.2)
     .    WRITE(NRM,942) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.3)
     .    WRITE(NRM,943) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.4)
     .    WRITE(NRM,944) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
         IF(NOPTP.EQ.6)
     .    WRITE(NRM,945) JPOIN,(COORD(JPOIN,IDIME),IDIME=1,NDIME)
        ENDIF            ! noptj.eq.1
       ENDIF             ! nopti.eq.4.or.nopti.eq.5
      ENDIF              ! jpoin.gt.0
C
      IF(KPOIN.GT.0) THEN
       IF(NOPTI.EQ.2.OR.NOPTI.EQ.3) THEN
        IF(NOPTO.EQ.1) THEN
         IF(NSKEW.EQ.1) THEN
          WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
         ELSE
          ICOUNT=0
          DO IELEM=1,NELMAXFL
           DO INODE=1,NNODE(IELEM)
            IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
           ENDDO
          ENDDO
          IF(NDIME.EQ.2) THEN
           IF(ICOUNT.LT.4) THEN            ! quadrilateral elements only
            WRITE(NRF,961) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
           ELSE
            WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
           ENDIF         ! icount.lt.4
          ENDIF          ! ndime.eq.2
          IF(NDIME.EQ.3) THEN
           IF(ICOUNT.LT.8) THEN            ! brick elements only
            WRITE(NRF,962) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
           ELSE
            WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
           ENDIF         ! icount.lt.8
          ENDIF          ! ndime.eq.3
         ENDIF           ! nskew.eq.0
        ELSE
         IF(NOPTP.EQ.2) THEN
          IF(NSKEW.EQ.1) THEN
           WRITE(NRF,942) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
          ELSE
           ICOUNT=0
           DO IELEM=1,NELMAXFL
            DO INODE=1,NNODE(IELEM)
             IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
            ENDDO
           ENDDO
           IF(NDIME.EQ.2) THEN
            IF(ICOUNT.LT.4) THEN           ! quadrilateral elements only
             WRITE(NRF,963) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,942) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.4
           ENDIF         ! ndime.eq.2
           IF(NDIME.EQ.3) THEN
            IF(ICOUNT.LT.8) THEN           ! brick elements only
             WRITE(NRF,964) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,942) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.8
           ENDIF         ! ndime.eq.3
          ENDIF          ! nskew.eq.0
         ENDIF           ! noptp.eq.2
         IF(NOPTP.EQ.3) THEN
          IF(NSKEW.EQ.1) THEN
           WRITE(NRF,943) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
          ELSE
           ICOUNT=0
           DO IELEM=1,NELMAXFL
            DO INODE=1,NNODE(IELEM)
             IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
            ENDDO
           ENDDO
           IF(NDIME.EQ.2) THEN
            IF(ICOUNT.LT.4) THEN           ! quadrilateral elements only
             WRITE(NRF,965) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,943) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.4
           ENDIF         ! ndime.eq.2
           IF(NDIME.EQ.3) THEN
            IF(ICOUNT.LT.8) THEN           ! brick elements only
             WRITE(NRF,966) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,943) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.8
           ENDIF         ! ndime.eq.3
          ENDIF          ! nskew.eq.0
         ENDIF           ! noptp.eq.3
         IF(NOPTP.EQ.4) THEN
          IF(NSKEW.EQ.1) THEN
           WRITE(NRF,944) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
          ELSE
           ICOUNT=0
           DO IELEM=1,NELMAXFL
            DO INODE=1,NNODE(IELEM)
             IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
            ENDDO
           ENDDO
           IF(NDIME.EQ.2) THEN
            IF(ICOUNT.LT.4) THEN           ! quadrilateral elements only
             WRITE(NRF,967) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,944) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.4
           ENDIF         ! ndime.eq.2
           IF(NDIME.EQ.3) THEN
            IF(ICOUNT.LT.8) THEN           ! brick elements only
             WRITE(NRF,968) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,944) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.8
           ENDIF         ! ndime.eq.3
          ENDIF          ! nskew.eq.0
         ENDIF           ! noptp.eq.4
         IF(NOPTP.EQ.6) THEN
          IF(NSKEW.EQ.1) THEN
           WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
          ELSE
           ICOUNT=0
           DO IELEM=1,NELMAXFL
            DO INODE=1,NNODE(IELEM)
             IF(LNODS(IELEM,INODE).EQ.KPOIN) ICOUNT=ICOUNT+1
            ENDDO
           ENDDO
           IF(NDIME.EQ.2) THEN
            IF(ICOUNT.LT.4) THEN           ! quadrilateral elements only
             WRITE(NRF,961) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.4
           ENDIF         ! ndime.eq.2
           IF(NDIME.EQ.3) THEN
            IF(ICOUNT.LT.8) THEN           ! brick elements only
             WRITE(NRF,962) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME),0
            ELSE
             WRITE(NRF,945) KPOIN,(COORD(KPOIN,IDIME),IDIME=1,NDIME)
            ENDIF        ! icount.lt.8
           ENDIF         ! ndime.eq.3
          ENDIF          ! nskew.eq.0
         ENDIF           ! noptp.eq.6
        ENDIF            ! noptj.eq.1
       ENDIF             ! nopti.eq.2.or.nopti.eq.3
      ENDIF              ! kpoin.gt.0
   30 CONTINUE
      IF (NOPTI.EQ.2.OR.NOPTI.EQ.3) WRITE(NRF,*) 'END_COORDINATES'
C
      IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN    ! initial veloc.
       IF(NOPTR.EQ.2) THEN
        WRITE(NIM,*) NPOIX
       ENDIF
      ENDIF
C
      DO 50 IPOIN=1,NPOIN
      JPOIN=JPUNT(IPOIN)
      KPOIN=KPUNT(IPOIN)
C
      IF(JPOIN.GT.0) THEN
       IF(NOPTI.EQ.4.OR.NOPTI.EQ.5) THEN
        IF(NOPTR.EQ.2) THEN
         IF(NDIME.EQ.2) WRITE(NIM,948) JPOIN,0.0,0.0
         IF(NDIME.EQ.3) WRITE(NIM,960) JPOIN,0.0,0.0,0.0
        ENDIF
       ENDIF
      ENDIF
   50 CONTINUE
C
      WRITE(NT,900)
      WRITE(NT,905) NPOIX          ! generated points for the analysis
      WRITE(NT,910) NELEM          ! generated elements for the analysis
C
c     IF(NFEMV.EQ.2) THEN
c      call writegeom(mem(plnods),mem(pmatno),mem(pproel),
c    .                mem(pprops),mem(pcoord),
c    .                ndime,nnode,nelem,ngrup,nmats,npoin,npoic,
c    .                nprel,nprop,name,jpoin,npoin1)
c
c      write(NFEM,*)'9999'
c     ENDIF
      IF(NFEMV.EQ.4) THEN
       write(18,*) 'titulo: problema'
       write(18,*) 'subtitulo: problema'
       write(18,*) 'a'
       write(18,*) 'a'
       write(18,*) 'a'
c
       write(18,*) 'c'
c
       nnode1=nnode(1)     ! chapuza: to be improved !!!
c
c----------------------------------------------------------- Write model
c
c  nnodx=1 => hexahedra
c        1 => hexahedra 20-noded
c        3 => tetrahedra
c        3 => tetrahedra 10-noded
c        3 => toblerone  >>  transformed to 3 tetrahedra
c        3 => triangle
c        3 => triangle 6-noded
c        4 => quadrilateral
c        4 => quadrilateral 8-noded
c       11 => line
c
       nnodx=0
       if(nnode1.eq. 8.and.ndime.eq.3) nnodx=1   ! hexahedra
       if(nnode1.eq.20.and.ndime.eq.3) nnodx=1   ! hexahedra 20-noded
       if(nnode1.eq. 4.and.ndime.eq.3) nnodx=3   ! tetrahedra
       if(nnode1.eq.10.and.ndime.eq.3) nnodx=3   ! tetrahedra 10-noded
       if(nnode1.eq. 6.and.ndime.eq.3) nnodx=3   ! toblerone
       if(nnode1.eq. 3.and.ndime.eq.2) nnodx=3   ! triangle
       if(nnode1.eq. 6.and.ndime.eq.2) nnodx=3   ! triangle 6-noded
       if(nnode1.eq. 4.and.ndime.eq.2) nnodx=4   ! quadrilateral
       if(nnode1.eq. 8.and.ndime.eq.2) nnodx=4   ! quadrilateral 8-noded
       if(nnode1.eq. 2) nnodx=11                 ! line
c
c**** deals with toblerone
c
       ntoble=0
       if(nnode1.eq.6.and.ndime.eq.3) ntoble=1
c
       if(ntoble.eq.0) then
        write(18,*) nelem, npoin, nnodx
       else
        nelem3=nelem*3
        write(18,*) nelem3, npoin, nnodx
       endif
c
       write(18,*) 'c'
       do ipoin=1,npoin
          x=coord(ipoin,1)
          y=coord(ipoin,2)
          if (ndime.eq.3) then
            z=coord(ipoin,3)
            write(18,901) ipoin,x,y,z      ! formatted
          else
            write(18,901) ipoin,x,y        ! formatted
          endif
       end do
  901  format(i8,3e15.6)
       write(18,*) 'c'
c
       if(ntoble.eq.0) then
        do ielem=1,nelem
c        write(18,*) ielem,(lnods(ielem,i),i=1,nnode1),matno(ielem)
         write(18,902) ielem,(lnods(ielem,i),i=1,nnode1),matno(ielem)
        enddo
  902   format(i8,27i8,i5)
       else
        do ielem=1,nelem
         do j=1,3
          ielem3=(ielem-1)*3+j
          if(j.eq.1) write(18,900) ielem3,lnods(ielem,1),lnods(ielem,2),
     .                                    lnods(ielem,3),lnods(ielem,4),
     .                                    matno(ielem)
          if(j.eq.2) write(18,900) ielem3,lnods(ielem,3),lnods(ielem,4),
     .                                    lnods(ielem,6),lnods(ielem,2),
     .                                    matno(ielem)
          if(j.eq.3) write(18,900) ielem3,lnods(ielem,4),lnods(ielem,2),
     .                                    lnods(ielem,5),lnods(ielem,6),
     .                                    matno(ielem)
         enddo
        enddo
       endif           ! ntoble.eq.0
      ENDIF            ! NFEMV.EQ.4
C
c 941 FORMAT(15I6)                          ! original
  941 FORMAT(I7,I7,13I7)                    ! large meshes
  942 FORMAT(I6,3F15.2)
  943 FORMAT(I6,3F15.3)
  944 FORMAT(I6,3F15.4)
c 945 FORMAT(I6,3F15.6)                     ! original
  945 FORMAT(I7,3F15.6)                     ! large meshes
  946 FORMAT(I6,2x,4I1,4(2x,F8.6))
  947 FORMAT(I6,2x,F10.6)
  948 FORMAT(I6,2(2x,F8.6))
  949 FORMAT(I6,2x,I1,2x,I2,(2x,F8.6))
c 950 FORMAT(15I6)                          ! original
  950 FORMAT(I7,I5,13I7)                    ! large meshes
  951 FORMAT(I6,2x,2I1,2x,I2,2(2x,F8.6))
c 952 FORMAT(I6,2x,3I1,2x,I2,3(2x,F8.6))    ! original
  952 FORMAT(I7,2x,3I1,2x,I2,3(2x,F8.6))    ! large meshes
  956 FORMAT(I6,2x,3I1,3(2x,F8.6))
  957 FORMAT(I6,2X,I2,2X,F10.6)
  958 FORMAT(I6,2x,4I1,4(2x,F8.6))
  959 FORMAT(I6,2x,5I1,5(2x,F8.6))
  960 FORMAT(I6,3(2x,F8.6))
  961 FORMAT(I6,2F15.6,I5)
  962 FORMAT(I6,3F15.6,I5)
  963 FORMAT(I6,2F15.2,I5)
  964 FORMAT(I6,3F15.2,I5)
  965 FORMAT(I6,2F15.3,I5)
  966 FORMAT(I6,3F15.3,I5)
  967 FORMAT(I6,2F15.4,I5)
  968 FORMAT(I6,3F15.4,I5)
  969 FORMAT(I6,2x,6I1,6(2x,F8.6))
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE SFRQ(S,T)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
C***  CALCULATES SHAPE FUNCTIONS
C
      SS=S*S
      TT=T*T
      ST=S*T
      SST=S*S*T
      STT=S*T*T
C
C*** SHAPE FUNCTIONS
      IF(NNODB.GT.4) GO TO 12
      SHAPE(1)=.25*(1.-S)*(1.-T)
      SHAPE(2)=.25*(1.+S)*(1.-T)
      SHAPE(3)=.25*(1.+S)*(1.+T)
      SHAPE(4)=.25*(1.-S)*(1.+T)
C
      RETURN
 12   IF(NNODB.GT.8) GO TO 13
      SHAPE(1)=(-1.0+ST+SS+TT-SST-STT)/4.0
      SHAPE(2)=(1.0-T-SS+SST)/2.0
      SHAPE(3)=(-1.0-ST+SS+TT-SST+STT)/4.0
      SHAPE(4)=(1.0+S-TT-STT)/2.0
      SHAPE(5)=(-1.0+ST+SS+TT+SST+STT)/4.0
      SHAPE(6)=(1.0+T-SS-SST)/2.0
      SHAPE(7)=(-1.0-ST+SS+TT+SST-STT)/4.0
      SHAPE(8)=(1.0-S-TT+STT)/2.0
      RETURN
   13 CONTINUE
      S1=S+1.0
      T1=T+1.0
      S2=S*2.0
      T2=T*2.0
      S9=S-1.0
      T9=T-1.0
      SHAPE(1)=0.25*S9*ST*T9
      SHAPE(2)=0.5*(1.0-SS)*T*T9
      SHAPE(3)=0.25*S1*ST*T9
      SHAPE(4)=0.5*S*S1*(1.0-TT)
      SHAPE(5)=0.25*S1*ST*T1
      SHAPE(6)=0.5*(1.0-SS)*T*T1
      SHAPE(7)=0.25*S9*ST*T1
      SHAPE(8)=0.5*S*S9*(1.0-TT)
      SHAPE(9)=(1.0-SS)*(1.0-TT)
      RETURN
      END
C_______________________________________________________________________
      SUBROUTINE SFRQ3D(S,T,Z)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
C***  CALCULATES SHAPE FUNCTIONS
C
c     SS=S*S
c     TT=T*T
c     ST=S*T
c     SST=S*S*T
c     STT=S*T*T
C
C*** SHAPE FUNCTIONS
C
      IF(NNODB.GT.8) GO TO 12
      SHAPE(1)=.125*(1.-S)*(1.-T)*(1.-Z)
      SHAPE(2)=.125*(1.+S)*(1.-T)*(1.-Z)
      SHAPE(3)=.125*(1.+S)*(1.+T)*(1.-Z)
      SHAPE(4)=.125*(1.-S)*(1.+T)*(1.-Z)
      SHAPE(5)=.125*(1.-S)*(1.-T)*(1.+Z)
      SHAPE(6)=.125*(1.+S)*(1.-T)*(1.+Z)
      SHAPE(7)=.125*(1.+S)*(1.+T)*(1.+Z)
      SHAPE(8)=.125*(1.-S)*(1.+T)*(1.+Z)
C
      RETURN
 12   WRITE(NT,*) '*** ERROR WITH NDIME=3 ***'
      STOP

      END
c-------------------------------------------------------- Write Geometry
c     subroutine writegeom(lnods,matno,proel,props,coord,ndime,
c    .                     nnode,nelem,ngrup,nmats,npoin,npoic,
c    .                     nprel,nprop,name,jpoin,npoin1)
c     implicit none
c     integer   ndime,nnode,nelem,ngrup,nmats,npoin,nprel,nprop,
c    .          numstp,ipoin,isystm,ielem,key,
c    .          lgrup,nnodp,imats,ltype,ntype,i,jpoin,npoin1,
c    .          lnods(nnode,nelem),matno(nelem),key1,key2,npoic,itype,
c    .          iskip
c     real      proel(nprel,ngrup),props(nprop,nmats),
c    .          coord(ndime,jpoin),x,y,z
c     character name*6
c     character opcode*1
c-------------------------------------- Write model elements  (Topology)
c     key=2
c     opcode='C'
c     write(11,3) key,opcode,name
c   3 format(1x,i4,A1,A6)
c
c     npoin1=npoin-npoic
c---------------------------------------------------------- Write model`
c        key=-1
c     do ipoin=1,npoin1
c        x=coord(1,ipoin)
c        y=coord(2,ipoin)
c        z=0.0
c        if (ndime.eq.3) z=coord(3,ipoin)
c        write(11,5) key,ipoin,x,y,z
c   5 format(1x,i2,i5,3E12.5)
c     end do
c      key=-3
c      write(11,6) key
c      key=3
c     write(11,3) key,opcode,name
c      key=-1
c     do ielem=1,nelem

c        iskip=0

c        lgrup=matno(ielem)
c        imats=int(proel(1,lgrup))
c        ITYPE=INT(PROEL(5,LGRUP))
c        if(itype.eq.4) proel(2,lgrup)=2.0   ! 2 nodes for contact elem.

c        if(itype.eq.32) then
c         iskip=1                            ! skip contact elements
c         nnodp=int(proel(2,lgrup))
c         if(ndime.eq.1) proel(2,lgrup)=2.0
c         if(ndime.eq.2) then
c          if(nnodp.eq.6) proel(2,lgrup)=4.0
c          if(nnodp.eq.9) proel(2,lgrup)=6.0
c         endif
c         if(ndime.eq.3) then
c          if(nnodp.eq.9) proel(2,lgrup)=6.0
c          if(nnodp.eq.12) proel(2,lgrup)=8.0
c         endif
c        endif

c        nnodp=int(proel(2,lgrup))
c        call Eltype_fgv(ndime,proel(1,lgrup),ntype)
c        key1=-1
c        key2=-2
c        imats=1                                      !para que funcione
c        if(iskip.eq.0) then
c         write(11,8) key1,ielem,ntype,lgrup,imats
c        endif
c   8 format(1x,i2,4i5)

c        if(iskip.eq.0) then
c        if(nnodp.gt.15) then
c         write(11,9) key2,(lnods(i,ielem),i=1,15)
c         write(11,9) key2,(lnods(i,ielem),i=16,nnodp)
c        else
c         write(11,9) key2,(lnods(i,ielem),i=1,nnodp)
c        endif
c        endif

c   9 format(1x,i2,15i5)
c     end do
c      key=-3
c      write(11,6) key
c   6 format(1x,i2)
c     end ! Write geometry
C_______________________________________________________________________
c     SUBROUTINE FRONTWIDTH
C
c     IMPLICIT REAL*8 (A-H,O-Z)
C
c     INCLUDE 'genren_i.f'                                 ! linux
c     INCLUDE 'genren_i.for'                               ! w95
C
C *** INITIALIZE
C
c     NFRON=1
c     DO 20 IELEM=1,NELEM
c     MNUME(IELEM)=0
c  20 NFWID(IELEM)=IELEM
C
C *** CONSIDER EACH ELEMENT IN TURN
C
c     DO 80 IELEM=1,NELEM
C
C *** AND CONSIDER EACH NODE OF THAT ELEMENT IN TURN
C
c     DO 70 INODE=1,LNODE
C
C *** IDENTIFY THE NODE
C
c     IPOIN=LNODS(IELEM,INODE)
c     DO 60 JELEM=1,NELEM
C
C *** ALL THE NODES OF JELEM ARE RELATED TO IPOIN
C
c     IF(JELEM.EQ.IELEM) GO TO 60
c     DO 50 JNODE=1,LNODE
C
C *** IDENTIFY A RELATED NODE
C
c     JPOIN=LNODS(JELEM,JNODE)
c     IF(IPOIN.EQ.JPOIN) GO TO 25
c  50 CONTINUE
c     GO TO 60
c  25 NNUME=MNUME(IELEM)
c     IF(NNUME.EQ.0) GO TO 40
c     DO 30 INUME=1,NNUME
c     IF(LENLS(IELEM,INUME).EQ.JELEM) GO TO 60
c  30 CONTINUE
C
C *** IF SUCH A RELATIONSHIP DOES NOT EXIST, PROCEED
C     WITH FORMING LENLS
C
c  40 MNUME(IELEM)=MNUME(IELEM)+1
c     LENLS(IELEM,MNUME(IELEM))=JELEM
C
C *** FIND THE WIDTH OF THE ORIGINAL MATRIX BAND
C
c     NWORK=IABS(IELEM-JELEM)
c     IF(NFRON.LT.NWORK) NFRON=NWORK
c  60 CONTINUE
c  70 CONTINUE
c  80 CONTINUE
c     MFRON=NFRON
c     NFROM=NFRON
c     DO 140 IELEM=1,NELEM
c     DO 90 JELEM=1,NELEM
c     NENEW(JELEM)=0
c  90 NEOLD(JELEM)=0
C
C *** INITIALIZE FOR NEW NODE NUMBER NITER
C
c     MAXIM=0
c     KOUNT=NITER
c     JELEM=NITER
c     NEOLD(JELEM)=IELEM
c     NENEW(IELEM)=JELEM
C
C *** EACH NODE RELATED TO NEW NODE NUMBER JPOIN IS THEN
C     ASSIGNED A NEW NUMBER
C
c 100 NNUME=MNUME(NEOLD(JELEM))
c     DO 110 INUME=1,NNUME
c     KELEM=LENLS(NEOLD(JELEM),INUME)
C
C *** IF THE RELATED NODE HAS ALREADY BEEN RENUMBERED,
C     SKIP TO THE NEXT RELATED NODE
C
c     IF(NENEW(KELEM).GT.0) GO TO 110
C
C *** ASSIGN A NEW NUMBER TO OLD NODE NUMBER KPOIN
C
c     KOUNT=KOUNT+1
c     NEOLD(KOUNT)=KELEM
c     NENEW(KELEM)=KOUNT
C
C *** THE DIFFERENCE BETWEEN THE NEW NUMBERS OF
C     THE RELATED NODES IS CHECKED
C
c     NFRON=IABS(JELEM-KOUNT)
C
C *** IF IT IS GREATER THAN THE BANDWIDTH PRODUCED BY A
C     PREVIOUS NUMBERING SCHEME (INCLUDING THE ORIGINAL
C     SCHEME), THE CURRENT SCHEME IS ABANDONED AND A NEW
C     ONE STARTED
C
c     IF(NFRON.GE.MFRON) GO TO 140
C
C *** IF THE DIFFERENCE IS THE BIGGEST VALUE YET ENCOUNTERED
C     FOR THE CURRENT NUMBERING SCHEME, ITS VALUE IS KEPT
C
c     IF(NFRON.GT.MAXIM) MAXIM=NFRON
c 110 CONTINUE
c     IF(KOUNT.EQ.NELEM) GO TO 120
C
C *** JPOIN IS NOW INCREMENTED SO THAT THE NODES RELATED TO THE
C     NEXT NEW JOINT NUMBER ARE ASSIGNED NEW NUMBERS
C
c     JELEM=JELEM+1
c     GO TO 100
C
C *** AT THIS POINT A COMPLETELY NEW NUMBERING SCHEME HAS BEEN
C     COMPLETED WHICH IS SUPERIOR TO PREVIOS NUMBERING SCHEMES
C     THE BANDWIDTH IS KEPT FOR COMPARISON WITH SUBSEQUENT
C     NUMBERING SCHEMES
C
c 120 MFRON=MAXIM
C
C *** THE POINTERS TO THE NEW NUMBERS ARE ALSO KEPT
C
c     DO 130 KELEM=1,NELEM
c     NFWID(KELEM)=NENEW(KELEM)
c 130 CONTINUE
c 140 CONTINUE
C
C *** WRITE THE SOLUTION OF THE PROBLEM
C
c     WRITE(NT,960) NFROM,NFRON
C
c     DO 190 IELEM=1,NELEM
c     NNODE1(NFWID(IELEM))=NNODE(IELEM)
c 190 MATNC(NFWID(IELEM))=MATNO(IELEM)
c     DO 200 IELEM=1,NELEM
c     NNODE(IELEM)=NNODE1(IELEM)
c 200 MATNO(IELEM)=MATNC(IELEM)
C
c     DO 150 IELEM=1,NELEM
c     DO 150 INODE=1,NNODE(IELEM)
c 150 LNOTS(NFWID(IELEM),INODE)=LNODS(IELEM,INODE)
c     DO 170 IELEM=1,NELEM
c     DO 180 INODE=1,NNODE(IELEM)
c 180 LNODS(IELEM,INODE)=LNOTS(IELEM,INODE)
c 170 CONTINUE
C
C *** ALL THE FORMATS USED IN THE PROGRAM
C
c 900 FORMAT()
c 955 FORMAT(10I5)
c 930 FORMAT(/,3X,18HNUMERO DE PUNTOS =,I4,
c    .       5X,21HNUMERO DE ELEMENTOS =,I4)
c 940 FORMAT(/,3X,8HELEMENTO,8(I3,X,'NODO'),12(I3),/)
c 950 FORMAT(2X,I6,3X,8(2X,I6),12(I3))
c 960 FORMAT(//,5X,22HEL SEMIANCHO DE FRENTE,I4,X,10HPASA A SER,I3)
C    .      //,5X,48HCON LA RENUMERACION DE LOS ELEMENTOS SIGUIENTE :,/)
c 970 FORMAT(/,500(2X,4(I4,2X,6HPASA A,I4,4X),/))
c 980 FORMAT(/,5X,35HLA RENUMERACION DEJA A LNODS COMO :)
c     RETURN
c     END
