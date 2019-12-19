C=============================================================== INCLUDE
C
      INTEGER*4 MPOIN,MDIME,MELEM,MNODM,MBLOC,MDIVX,MDIVY,MDIVZ,MREPN
C
      PARAMETER(MPOIN=8000000,          ! Change dimensions if necessary
     .          MDIME=3,
     .          MELEM=8000000,
     .          MNODM=9,                ! only this fixed
     .          MBLOC=800,
     .          MDIVX=800,
     .          MDIVY=800,
     .          MDIVZ=800,
     .          MREPN=MPOIN )
C
      INTEGER*4    LNODS(MELEM,MNODM), MATNO(MELEM), NINTE(MELEM),
     .             NGAUS(MELEM),       NNODE(MELEM), LTRAN(MELEM,5),
     .             JPUNT(MPOIN),       KPUNT(MPOIN), LPUNT(MPOIN),
     .             LNODT(MELEM,MNODM), MATNT(MELEM), NINTT(MELEM),
     .             NGAUT(MELEM)
      COMMON/MESHA/LNODS,              MATNO,        NINTE,
     .             NGAUS,              NNODE,        LTRAN,
     .             JPUNT,              KPUNT,        LPUNT,
     .             LNODT,              MATNT,        NINTT,
     .             NGAUT
C
      REAL*8       COORD(MPOIN,MDIME), SHAPE(MNODM), COORT(MPOIN,MDIME)
      COMMON/MESHB/COORD,              SHAPE,        COORT
C
      INTEGER*4    NPOIN,NELEM,NDIME,MNODE,NP,NT,NR,NNODB,NTYPE,NS,NI,
     .             NA,NRF,NSF,NG,NRM,NSM,NIM,NFEM,NACT,
     .             LNODE,NITER,NTRAN
      COMMON/MESH1/NPOIN,NELEM,NDIME,MNODE,NP,NT,NR,NNODB,NTYPE,NS,NI,
     .             NA,NRF,NSF,NG,NRM,NSM,NIM,NFEM,NACT,
     .             LNODE,NITER,NTRAN
C
      INTEGER*4    MLNOD(MELEM,MNODM), LREPN(MPOIN), LASOC(MPOIN),
     .             LFINN(MPOIN),       LFASC(MPOIN), MMATO(MBLOC),
     .             MMINT(MBLOC),       MMGAU(MBLOC), NBNOD(MBLOC),
     .             MMDIV(MBLOC),       NUEVN(MNODM)
      COMMON/MESH2/MLNOD,              LREPN,        LASOC,
     .             LFINN,              LFASC,        MMATO,
     .             MMINT,              MMGAU,        NBNOD,
     .             MMDIV,              NUEVN
C
      REAL*8       WEITX(MDIVX), WEITY(MDIVY), WEITZ(MDIVZ),
     .             TCORD(MPOIN,MDIME), TOLCO
      COMMON/MESH3/WEITX,        WEITY,        WEITZ,
     .             TCORD,              TOLCO
C
c     INTEGER*4    MNUME(MELEM), LENLS(MELEM,MELEM), LNOTS(MELEM,MNODM),
c    .             NEOLD(MELEM), NENEW(MELEM),       NFWID(MELEM),
c    .             MATNC(MELEM), NNODE1(MELEM)
c     COMMON/MESH4/MNUME,        LENLS,              LNOTS,
c    .             NEOLD,        NENEW,              NFIWD,
c    .             MATNC,        NNODE1
C
C=============================================================== INCLUDE
