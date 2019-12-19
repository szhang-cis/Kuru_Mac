      SUBROUTINE MICROS2(TGAUAT,TGAUST,TGAUIT,DTEMPT,
     .                    BASMM, BASCC, BASKK,
     .                    TSOE2, TSOE1, TSOC1,
     .                    IPLAT,
     .                   ALPHAM, INUPC,INDEX3)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES THE PHASE-CHANGE FUNCTION, DENSITY,
C     CAPACITY & CONDUCTIVITY ACCORDING WITH THE MICROSTRUCTURAL
C     MODEL NUMBER 2 (IPCMO=2)
C
C     GRAY CAST IRON MICROSTRUCTURAL MODEL: FOR HYPOEUTECTIC, EUTECTIC &
C     HYPEREUTECTIC COMPOSITIONS
C
C     Solidification of the primary austenite (explicit f_pc-T function)
C     Eutectic solidification (nucleation and growth models)
C
C***********************************************************************
C
C     Index of variables
C
C     Input:
C     TGAUAT= Temperature at time t
C     TGAUST= Temperature at time t+dt
C     TGAUIT= Initial temperature
C     DTEMPT= Temperature rate (or material derivative)
C
C     Input/output:
C     BASMM = Density
C     BASCC = Capacity coefficient
C     BASKK = Conductivity
C
C     Output:
C     TSOE2 = each L*phase-change function at time t
C     TSOE1 = each L*phase-change function at time t+dt
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C     ALPHAM(IN+1)=liquid fraction
C
C     ALPHAM(IN+2)=primary austenite or graphite fraction (solid)
C                  (for hypoeutectic or hypereutectic compositions)
C
C     ALPHAM(IN+3)=graphite (gray) eutectic fraction (solid)
C     ALPHAM(IN+4)=graphite eutectic grain density
C     ALPHAM(IN+5)=graphite eutectic grain radius
C
C     ALPHAM(IN+6)=cementite (white) eutectic fraction (solid)
C     ALPHAM(IN+7)=cementite eutectic grain density
C     ALPHAM(IN+8)=cementite eutectic grain radius
C
C     ALPHAM(IN+9)=austenite fraction
C     ALPHAM(IN+10)=graphite fraction
C     ALPHAM(IN+11)=average graphite lamellar spacing
C     ALPHAM(IN+12)=average graphite lamellae thickness
C
C     Auxiliar microstructural variables (not printed; see pointes.f):
C     ALPHAM(IN+13)=index for graphite nucleation
C     ALPHAM(IN+14)=index for cementite nucleation
C
C
C     Primary-austenite or primary-graphite solidification
C
C     Eutectic solidification (graphite and cementite)
C
C
C     Notes:
C
C     IAFLOJ=fraction correction form (=0, relaxed correction;
C                                      =1, full correction)
C
C     More global iterations are obtained when this micro computation
C     is only allowed for iitert>0
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES    ! thermal-microstructural
C
      INCLUDE 'nued_om.f'
C
C**** THERMAL VARIABLES
C
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION ALPHAM(NHISTM)
      DIMENSION BASKK(NDIMETO)
C
      DIMENSION TSOE1(5), TSOE2(5), TSOC1(5)
C
C**** CHECKS LIQUID FRACTION
C
      IN=INUPC
      FLIQD=ALPHAM(IN+1)
C
      IF(FLIQD.EQ.0.0D0) THEN      ! gets out if liquid fraction is zero
       INUPC=INUPC+14              ! less than NBASES; see pointes.f
       RETURN
      ENDIF
C
      TWOPIT=6.283185307179586D0
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      TEINF=VPLAT(IPLAT,1)
      TESUP=VPLAT(IPLAT,2)
      HENER=VPLAT(IPLAT,3)
      DELEE=TESUP-TEINF
C
      CA=VPLAT(IPLAT,6)                    ! carbon
      SI=VPLAT(IPLAT,7)                    ! silicon
      PH=VPLAT(IPLAT,8)                    ! phosporus
      CU=VPLAT(IPLAT,9)                    ! copper
      QM=VPLAT(IPLAT,10)                   ! manganese
      QG=VPLAT(IPLAT,11)                   ! magnesium
      CR=VPLAT(IPLAT,12)                   ! chromium
      AN=VPLAT(IPLAT,13)                   ! nickel
C
      DENS=VPLAT(IPLAT,14)
      ESPH=VPLAT(IPLAT,15)
      COND=VPLAT(IPLAT,16)
      HENER=VPLAT(IPLAT,17)
C
      IAUSM=INT(VPLAT(IPLAT,18))
      AUSTS=VPLAT(IPLAT,19)
      AUSTL=VPLAT(IPLAT,20)
      AUSLA=VPLAT(IPLAT,21)
C
      ANUCG=VPLAT(IPLAT,22)
      ANUCC=VPLAT(IPLAT,23)
      BGROG=VPLAT(IPLAT,24)
      BGROC=VPLAT(IPLAT,25)
C
      AK1LA=VPLAT(IPLAT,26)
      AK2LA=VPLAT(IPLAT,27)
      MAXNGL=INT(VPLAT(IPLAT,28))
C
      IFPCDT=INT(VPLAT(IPLAT,29))
      IAFLOJ=INT(VPLAT(IPLAT,30))
C
      IV=30                      ! number of VPLAT defined in input data
C
      TL=VPLAT(IPLAT,IV+1)
      TS=VPLAT(IPLAT,IV+2)
      TEX=VPLAT(IPLAT,IV+3)
      TEG=VPLAT(IPLAT,IV+4)
      TEC=VPLAT(IPLAT,IV+5)
      PK=VPLAT(IPLAT,IV+6)
      CEQ=VPLAT(IPLAT,IV+7)
      FGREUT=VPLAT(IPLAT,IV+8)
C
C**** TRANSFER (REST OF) "ALPHAM" TO MICROSCOPICAL VARIABLES
C
      FSOLA=ALPHAM(IN+2)
C
      FSOLG=ALPHAM(IN+3)
      XXNGN=ALPHAM(IN+4)
      RRRGM=ALPHAM(IN+5)
C
      FSOLC=ALPHAM(IN+6)
      XXNCN=ALPHAM(IN+7)
      RRRCM=ALPHAM(IN+8)
C
      FAUST=ALPHAM(IN+9)
      FGRAP=ALPHAM(IN+10)
      AGLAS=ALPHAM(IN+11)
      AGLAT=ALPHAM(IN+12)
C
      INDEXG= INT(ALPHAM(IN+13))
      INDEXC= INT(ALPHAM(IN+14))       ! less than NBASES; see pointes.f
C
C**** LOAD LAST CONVERGED VALUES
C
      FLIQDO=FLIQD
C
      FSOLAO=FSOLA
C
      FSOLGO=FSOLG
      XXNGNO=XXNGN
      RRRGMO=RRRGM
C
      FSOLCO=FSOLC
      XXNCNO=XXNCN
      RRRCMO=RRRCM
C
      FAUST0=FAUST
      FGRAP0=FGRAP
      AGLASO=AGLAS
      AGLATO=AGLAT
C
      INDEXGO=INDEXG
      INDEXCO=INDEXC
C
C**** LOAD CONVECTIVE TERM ASSOCIATED WITH MICROSTRUCTURAL VARIABLES
C
      RRRGMC=0.0D0
      RRRCMC=0.0D0
c     IF(ICONVT.EQ.1) THEN      ! to be revised
c      RRRGMC=APLUOT(INDEX3+1)
c      RRRCMC=APLUOT(INDEX3+2)
c     ENDIF
C
C**** ESTABLISH VARIABLES OF MODEL 2
C
      FSOLID=1.0D0-FLIQD        ! solid fraction
      FSOLIDO=1.0D0-FLIQDO
C
C**** SOLVES MODEL 2
C
      TSOCA =0.0D0
      TSOCGN=0.0D0
      TSOCGG=0.0D0
      TSOCCN=0.0D0
      TSOCCG=0.0D0
C
C**** PRIMARY-AUSTENITE SOLIDIFICATION
C
C     Note: it is assumed that the primary-austenite solidification ends
C           when the graphite eutectic temperature is reached.
C
      TLX=TL
      IF(IAUSM.EQ.1) TLX=AUSTL
C
      GUA=TLX-TGAUST                   ! austenite undercooling
      ICORREC=0
C
      IF(GUA.GT.0.0D0.AND.FSOLA.LT.1.0D0.AND.FSOLID.LT.1.0D0.AND.
     .   CEQ.LT.4.26D0) THEN
       IF(IAUSM.EQ.0) THEN
        IF(TGAUST.GT.TEG) THEN
         FSOLA=(1.0D0/(1.0D0-PK))*(TGAUST-TL)/(TGAUST-TEX)
         TSOCA=-(TL-TEX)/((1.0D0-PK)*(TGAUST-TEX)**2.0D0)
        ELSE
         FSOLA=(1.0D0/(1.0D0-PK))*(TEG-TL)/(TEG-TEX)
         TSOCA=0.0D0
        ENDIF
       ENDIF
       IF(IAUSM.EQ.1) THEN
        IF(TGAUST.GT.AUSTS) THEN
         FSOLA=1.0D0-(TGAUST-AUSTS)/(AUSTL-AUSTS)
         TSOCA=1.0D0/(AUSTL-AUSTS)
        ELSE
         FSOLA=1.0D0
         TSOCA=0.0D0
        ENDIF
        FSOLA=FSOLA*(1.0D0/(1.0D0-PK))*(TEG-TL)/(TEG-TEX)
        TSOCA=TSOCA*(1.0D0/(1.0D0-PK))*(TEG-TL)/(TEG-TEX)
       ENDIF
C
       IF(FSOLA.LT.FSOLAO) THEN        ! control
        FSOLA=FSOLAO                   ! no remelting
        TSOCA=0.0D0
       ENDIF
C
       FSOLIX=FSOLA+FSOLG+FSOLC        ! control
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        F1=FSOLA
        FSOLA=1.0D0-FSOLG-FSOLC
        TSOCA=TSOCA*FSOLA/F1
       ENDIF
C
       ICORREC=1
      ENDIF                            ! gua.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FSOLID=FSOLA+FSOLG+FSOLC
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
C**** PRIMARY-GRAPHITE SOLIDIFICATION
C
C     Note: it is assumed that the primary-graphite solidification ends
C           when the graphite eutectic temperature is reached.
C
      GUA=TL-TGAUST                    ! graphite undercooling
      ICORREC=0
C
      IF(GUA.GT.0.0D0.AND.FSOLA.LT.1.0D0.AND.FSOLID.LT.1.0D0.AND.
     .   CEQ.GE.4.26D0) THEN
       IF(IAUSM.EQ.0) THEN             ! not implemented yet
        IF(TGAUST.GT.TEG) THEN
         FSOLA=(1.0D0/(1.0D0-PK))*(TGAUST-TL)/(TGAUST-TEX)
         TSOCA=-(TL-TEX)/((1.0D0-PK)*(TGAUST-TEX)**2.0D0)
        ELSE
         FSOLA=(1.0D0/(1.0D0-PK))*(TEG-TL)/(TEG-TEX)
         TSOCA=0.0D0
        ENDIF
       ENDIF
       IF(IAUSM.EQ.1) THEN
        AUSLA=AUSTS                    ! pg latent heat
        IF(TGAUST.GT.TEG) THEN
         FSOLA=1.0D0-(TGAUST-TEG)/(TL-TEG)
         TSOCA=1.0D0/(TL-TEG)
        ELSE
         FSOLA=1.0D0
         TSOCA=0.0D0
        ENDIF
        FSOLA=FSOLA*AUSTL
        TSOCA=TSOCA*AUSTL
       ENDIF
C
       IF(FSOLA.LT.FSOLAO) THEN        ! control
        FSOLA=FSOLAO                   ! no remelting
        TSOCA=0.0D0
       ENDIF
C
       FSOLIX=FSOLA+FSOLG+FSOLC        ! control
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        F1=FSOLA
        FSOLA=1.0D0-FSOLG-FSOLC
        TSOCA=TSOCA*FSOLA/F1
       ENDIF
C
       ICORREC=1
      ENDIF                            ! gua.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FSOLID=FSOLA+FSOLG+FSOLC
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
C**** GRAPHITE EUTECTIC SOLIDIFICATION
C
      GUG=TEG-TGAUST                   ! graphite eutectic undercooling
      ICORREC=0
C
      IF(GUG.GT.0.0D0.AND.FSOLG.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
       IF(XXNGN.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEXG=1    ! end of nucl.
       IF(INDEXG.EQ.0) THEN
        XXNGN=ANUCG*GUG*GUG            ! nucleation
        TSOCGN=ANUCG*2.0D0*GUG
       ENDIF
C
       IF(XXNGN.LT.XXNGNO) THEN        ! control
        XXNGN=XXNGNO
        TSOCGN=0.0D0
       ENDIF
C
c      IF(DTEMPT.LT.0.0D0) THEN        ! incremental form-to be revised
c       if(fsolg.lt.0.05d0) then
c        XXNGN=XXNGN-                  ! nucleation
c    .         ANUCG*2.0D0*GUG*DTEMPT*DTIMET
c        TSOCGN=ANUCG*2.0D0*GUG
c       endif
c      ENDIF
C
       FSOLG=2.0D0/3.0D0*TWOPIT*XXNGN*RRRGM**3.0D0             ! control
       FSOLIX=FSOLA+FSOLG+FSOLC
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        FSOLG=1.0D0-FSOLA-FSOLC
        F1=XXNGN
        XXNGN=3.0D0*FSOLG/(2.0D0*TWOPIT*RRRGM**3.0D0)
        TSOCGN=TSOCGN*XXNGN/F1
        IF(XXNGN.LT.XXNGNO) THEN
         XXNGN=XXNGNO
         TSOCGN=0.0D0
        ENDIF
       ENDIF
C
       RRRGM=RRRGM+
     .  (BGROG*GUG*GUG-RRRGMC)*DTIMET  ! growth
       TSOCGG=BGROG*2.0D0*GUG*DTIMET
C
       IF(RRRGM.LT.RRRGMO) THEN        ! control
        RRRGM=RRRGMO
        TSOCGG=0.0D0
       ENDIF
C
       FSOLG=2.0D0/3.0D0*TWOPIT*XXNGN*RRRGM**3.0D0             ! control
       FSOLIX=FSOLA+FSOLG+FSOLC
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        FSOLG=1.0D0-FSOLA-FSOLC
        F1=RRRGM
        RRRGM=(3.0D0*FSOLG/(2.0D0*TWOPIT*XXNGN))**(1.0D0/3.0D0)
        TSOCGG=TSOCGG*RRRGM/F1
        IF(RRRGM.LT.RRRGMO) THEN
         RRRGM=RRRGMO
         TSOCGG=0.0D0
        ENDIF
       ENDIF
C
       FAUST=(1.0D0-FGREUT)*FSOLG      ! austenite in grap. eutec. frac.
       FGRAP=FGREUT*FSOLG              ! graphite in grap. eutec. frac.
C
       IF(AK1LA.GT.0.0D0.AND.AK2LA.GT.0.0D0) THEN     ! lamellar spacing
        IF(RRRGM.GT.0.D0) THEN
         DO IGL=1,MAXNGL
          ALAMGR=TWOPIT*RRRGM/IGL
          GUI=AK1LA*(RRRGM-RRRGMO)/DTIMET*ALAMGR+AK2LA/ALAMGR
          IF(GUG.GT.GUI) THEN
           IGX=IGL
           GO TO 11
          ENDIF
         ENDDO
         IGX=MAXNGL                    ! max. number of graph. lamellas
   11    ALAMGR=TWOPIT*RRRGM/IGX
         AGLAS=AGLAS*RRRGMO+ALAMGR*(RRRGM-RRRGMO)
         AGLAS=AGLAS/RRRGM             ! aver. graphite lamellar spacing
         AGLAT=FGREUT*AGLAS            ! aver. graphite lamellae thickn.
        ENDIF
       ENDIF
C
       ICORREC=1
      ENDIF                            ! gug.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FSOLID=FSOLA+FSOLG+FSOLC
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
      TSOCGN=TSOCGN*2.0D0/3.0D0*TWOPIT*RRRGM**3.0D0
      TSOCGG=TSOCGG*2.0D0*TWOPIT*XXNGN*RRRGM**2.0D0
C
C**** CEMENTITE EUTECTIC SOLIDIFICATION
C
      GUC=TEC-TGAUST                   ! cementite eutectic undercooling
      ICORREC=0
C
      IF(GUC.GT.0.0D0.AND.FSOLC.LT.1.0D0.AND.FSOLID.LT.1.0D0) THEN
       IF(XXNCN.GT.0.0D0.AND.DTEMPT.GT.0.0D0) INDEXC=1    ! end of nucl.
       IF(INDEXC.EQ.0) THEN
        XXNCN=ANUCC*GUC*GUC            ! nucleation
        TSOCCN=ANUCC*2.0D0*GUC
       ENDIF
C
       IF(XXNCN.LT.XXNCNO) THEN        ! control
        XXNCN=XXNCNO
        TSOCCN=0.0D0
       ENDIF
C
       FSOLC=2.0D0/3.0D0*TWOPIT*XXNCN*RRRCM**3.0D0             ! control
       FSOLIX=FSOLA+FSOLG+FSOLC
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        FSOLC=1.0D0-FSOLA-FSOLG
        F1=XXNCN
        XXNCN=3.0D0*FSOLC/(2.0D0*TWOPIT*RRRCM**3.0D0)
        TSOCCN=TSOCCN*XXNCN/F1
        IF(XXNCN.LT.XXNCNO) THEN
         XXNCN=XXNCNO
         TSOCCN=0.0D0
        ENDIF
       ENDIF
C
       RRRCM=RRRCM+
     .  (BGROC*GUC*GUC-RRRCMC)*DTIMET  ! growth
       TSOCCG=BGROC*2.0D0*GUC*DTIMET
C
       IF(RRRCM.LT.RRRCMO) THEN        ! control
        RRRCM=RRRCMO
        TSOCCG=0.0D0
       ENDIF
C
       FSOLC=2.0D0/3.0D0*TWOPIT*XXNCN*RRRCM**3.0D0             ! control
       FSOLIX=FSOLA+FSOLG+FSOLC
       IF(FSOLIX.GT.1.0D0) THEN        ! assumption valid for small dt
        FSOLC=1.0D0-FSOLA-FSOLG
        F1=RRRCM
        RRRCM=(3.0D0*FSOLC/(2.0D0*TWOPIT*XXNCN))**(1.0D0/3.0D0)
        TSOCCG=TSOCCG*RRRCM/F1
        IF(RRRCM.LT.RRRCMO) THEN
         RRRCM=RRRCMO
         TSOCCG=0.0D0
        ENDIF
       ENDIF
C
       ICORREC=1
      ENDIF                            ! guc.gt.0.0.and....
C
      IF(ICORREC.EQ.1) THEN            ! correction
       IF(IAFLOJ.EQ.1) THEN            ! strict conditions
        FSOLID=FSOLA+FSOLG+FSOLC
       ENDIF                           ! iafloj.eq.1
      ENDIF                            ! icorrec.eq.1
C
      TSOCCN=TSOCCN*2.0D0/3.0D0*TWOPIT*RRRCM**3.0D0
      TSOCCG=TSOCCG*2.0D0*TWOPIT*XXNCN*RRRCM**2.0D0
C
      FSOLID=FSOLA+FSOLG+FSOLC         ! always
C
C**** LIQUID FRACTION
C
      FLIQD=1.0D0-FSOLID                   ! liquid fraction
C
C**** DEFINES THE "MACROSCOPICAL" PROPERTIES (DENS., CAPACITY & CONDUC.)
C
c     BASMM=DENS                           ! no change is assumed
c     BASCC=ESPH
c     BASKK(1)=COND
C
C**** ESTABLISHES LATENT HEAT * PHASE-CHANGE FUNCTION (at time t & t+dt)
C
      IF(IAUSM.EQ.0) THEN                  ! classical
       TSOE1(IPLAT)=FLIQD                  ! f_pc at time t+dt
       TSOE2(IPLAT)=FLIQDO                 ! f_pc at time t
       TSOE2(IPLAT)=TSOE2(IPLAT)*HENER     ! L*f_pc at time t
       TSOE1(IPLAT)=TSOE1(IPLAT)*HENER     ! L*f_pc at time t+dt
      ENDIF
      IF(IAUSM.EQ.1) THEN
       TSOE2(IPLAT)=-HENER*(FSOLGO+FSOLCO)-AUSLA*FSOLAO ! L*f_pc at t
       TSOE1(IPLAT)=-HENER*(FSOLG+FSOLC)  -AUSLA*FSOLA  ! L*f_pc at t+dt
      ENDIF
C
C**** OPTIONS IN THE EVALUATION OF df_pc/dT
C
C     1) standard df_pc/dT = delta f_pc / delta T (Diego's thesis)
C     2) standard df_pc/dT = delta f_pc / delta T (always ge 0)
C     3) df_pc/dT "exact"
C     4) df_pc/dT "exact" (always ge 0)
C     5) df_pc/dT = 0
C
C     Notes:
C
C     L is constant
C     One possible option to obtain recalescense is to consider a
C     graphite/cementite-dependent latent heat value (L increases with
C     FG/FC)
C     Only one value of ICACOT can be defined, i.e., ICACOT does not
C     depend of any particular phase-change. This can not be useful for
C     more than one microstructural phase-changes.
C
      IF(ICONVT.EQ.1) THEN
       IF(IEGFPC.EQ.0) THEN
        IF(IFPCDT.EQ.2.OR.IFPCDT.EQ.4.OR.IFPCDT.EQ.5)
     .   CALL RUNENDT('ERROR: IFPCDT NE 1-3 WHEN ICONVT=1')
       ENDIF
      ENDIF
C
      GO TO (1,2,3,4,5) IFPCDT
C
    1 ICACOT=0
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
c      IF(IITERT.GT.0)
c    .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
       TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C
    2 ICACOT=1
      VELTOT=1.D-10
      VELA1T=DABS(DTEMPT)
      IF(VELA1T.GT.VELTOT) THEN
c      IF(IITERT.GT.0)                     ! to be revised
c    .  TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
        TSOC1(IPLAT)=(TSOE1(IPLAT)-TSOE2(IPLAT))/(DTEMPT*DTIMET)
      ENDIF
      GO TO 10
C
    3 ICACOT=0
      IF(IAUSM.EQ.0) THEN
       TSOC1(IPLAT)=TSOCA+TSOCGN+TSOCGG+TSOCCN+TSOCCG
       TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      ENDIF
      IF(IAUSM.EQ.1) THEN
       TSOC1(IPLAT)=TSOCA*AUSLA+(TSOCGN+TSOCGG+TSOCCN+TSOCCG)*HENER
      ENDIF
      GO TO 10
C
    4 ICACOT=1
      IF(IAUSM.EQ.0) THEN
       TSOC1(IPLAT)=TSOCA+TSOCGN+TSOCGG+TSOCCN+TSOCCG
       TSOC1(IPLAT)=TSOC1(IPLAT)*HENER
      ENDIF
      IF(IAUSM.EQ.1) THEN
       TSOC1(IPLAT)=TSOCA*AUSLA+(TSOCGN+TSOCGG+TSOCCN+TSOCCG)*HENER
      ENDIF
      GO TO 10
C
    5 ICACOT=1
      TSOC1(IPLAT)=0.0D0
      GO TO 10
C
   10 CONTINUE
C
C**** TRANSFER MICROSCOPICAL VARIABLES TO "ALPHAM"
C
      ALPHAM(IN+1)=FLIQD
C
      ALPHAM(IN+2)=FSOLA
C
      ALPHAM(IN+3)=FSOLG
      ALPHAM(IN+4)=XXNGN
      ALPHAM(IN+5)=RRRGM
C
      ALPHAM(IN+6)=FSOLC
      ALPHAM(IN+7)=XXNCN
      ALPHAM(IN+8)=RRRCM
C
      ALPHAM(IN+9)=FAUST
      ALPHAM(IN+10)=FGRAP
      ALPHAM(IN+11)=AGLAS
      ALPHAM(IN+12)=AGLAT
C
      ALPHAM(IN+13)=FLOAT(INDEXG)
      ALPHAM(IN+14)=FLOAT(INDEXC)
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+14                   ! less than NBASES; see pointes.f
C
C**** INCREMENTS "APLUOT" INDEX
C
      IF(ICONVT.EQ.1) INDEX3=INDEX3+2
C
      RETURN
      END
