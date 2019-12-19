      SUBROUTINE IDEPROSS(PROPSS)
C***********************************************************************
C
C**** THIS ROUTINE ORDER THE PROPERTIES FOR THERMAL SOLID ELEMENTS
C
C***********************************************************************
C
C     NDENS
C     VDENS(1)=Density values at       VDENS(2)=temperature
C
C     NCAPA
C     VCAPA(1)=Specific Heat values at VCAPA(2)=temperature
C
C     NCOND
C     VCOND(1)=Conductivity values at  VCOND(2)=temperature
C
C     NPLAT
C     VPLAT(1)=Solidus Temperature
C     VPLAT(2)=Liquidus Temperature
C     VPLAT(3)=Latent Heat
C
C     NSUBD(IDIME,1)=Number of divisions to integrate the forces
C     NSUBD(IDIME,2)=Number of divisions to integrate the phase-change
C                    matrix
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'
C
C**** MICROSTRUCTURAL VARIABLES
C
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
C
      DIMENSION PROPSS(*)
C
      I1=NPROPM+3      ! ISOTRS=PROPSS(NPROPM+1),IFILLS=PROPSS(NPROPM+2)
      NDENSS=INT(PROPSS(I1))
C
      DO IDENS=1,NDENSS
       IA=I1-1+2*IDENS
       IB=IA+1
       VDENSS(IDENS,1)=PROPSS(IA)
       VDENSS(IDENS,2)=PROPSS(IB)
      ENDDO
C
      I2=I1+2*NDENSS+1
      NCAPAS=INT(PROPSS(I2))
C
      DO ICAPA=1,NCAPAS
       IA=I2-1+2*ICAPA
       IB=IA+1
       VCAPAS(ICAPA,1)=PROPSS(IA)
       VCAPAS(ICAPA,2)=PROPSS(IB)
      ENDDO
C
      I3=I2+2*NCAPAS+1
      NCONDS=INT(PROPSS(I3))
C
      DO ICOND=1,NCONDS
       IA=I3-1+2*ICOND
       IB=IA+1
       VCONDS(ICOND,1,1)=PROPSS(IA)
       VCONDS(ICOND,2,1)=PROPSS(IB)
      ENDDO
C
      I4=I3+2*NCONDS+1
      NPLATS=INT(PROPSS(I4))
C
      IF(NPLATS.NE.0) THEN

       call runends('nplats > 0 not implemented in idepross')

c      IX=I4
c      DO IPLAT=1,NPLAT
c       I4A=IX+1
c       IPCFO=INT(PROPST(I4A))
c       VPLAT(IPLAT,4)=PROPST(I4A)
c       I4B=I4A+1
c       IPCMO=INT(PROPST(I4B))
c       VPLAT(IPLAT,5)=PROPST(I4B)
C
c       NPCPC=2                               ! 2=ipcfo,ipcmo
c       IF(IPCFO.EQ.0) THEN
C
c        IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1) THEN
c         IA=I4B+1
c         IB=IA+1
c         IC=IB+1
c         VPLAT(IPLAT,1)=PROPST(IA)
c         VPLAT(IPLAT,2)=PROPST(IB)
c         VPLAT(IPLAT,3)=PROPST(IC)
c         IX=IX+NPCPC+3                       ! 2=ipcmo,ipcfo, 3=Ts,Tl,L
c        ENDIF
C
c        IF(IPCMO.EQ.1.OR.IPCMO.EQ.2) THEN
c         I4C=I4B+1
c         NPCFU=INT(PROPST(I4C))
c         VPLAT(IPLAT,6)=PROPST(I4C)
c         I4D=I4C+1
c         NLATE=INT(PROPST(I4D))
c         VPLAT(IPLAT,7)=PROPST(I4D)
c         IF(NLATE.EQ.0) THEN                 ! constant L
c          DO IPCFU=1,NPCFU
c           IA=I4D-1+2*IPCFU
c           IB=IA+1
c           VPCFU(IPLAT,IPCFU,1)=PROPST(IA)
c           VPCFU(IPLAT,IPCFU,2)=PROPST(IB)
c          ENDDO
c          VPLAT(IPLAT,1)=VPCFU(IPLAT,1,2)      ! Ts
c          VPLAT(IPLAT,2)=VPCFU(IPLAT,NPCFU,2)  ! Tl
c          I4E=IB+1
c          VPLAT(IPLAT,3)=PROPST(I4E)
c          IX=IX+2*NPCFU+NPCPC+3              ! 3=npcfu,nlate,L
c         ELSE                                ! temp.-dep. L
c          I4D=I4D+1
c          NTYPC=INT(PROPST(I4D))             ! secant or tangent L
c          VPLAT(IPLAT,8)=PROPST(I4D)
c          DO IPCFU=1,NPCFU
c           IA=I4D-2+3*IPCFU
c           IB=IA+1
c           IC=IB+1
c           VPCFU(IPLAT,IPCFU,1)=PROPST(IA)
c           VPCFU(IPLAT,IPCFU,2)=PROPST(IB)
c           VPCFU(IPLAT,IPCFU,3)=PROPST(IC)
c          ENDDO
C
c          IF(NTYPC.EQ.1) THEN                ! tangent L
c           VPCAU(IPLAT,1)=VPCFU(IPLAT,1,3)   ! computes L^s
c           DO IPCFU=2,NPCFU
c            VAUXX=0.0D+00
c            DO JPCFU=2,IPCFU
c             I1=JPCFU-1
c             I2=JPCFU
c             IF(IPCMO.EQ.1) THEN
c              VAUXX=VAUXX+(VPCFU(IPLAT,I2,3)+VPCFU(IPLAT,I1,3))*
c    .               (VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
c    .               VPCFU(IPLAT,IPCFU,1)
c             ELSE                            ! considers f_pc
c              VAUXX=VAUXX+(VPCFU(IPLAT,I2,3)+VPCFU(IPLAT,I1,3))*
c    .               (-VPCFU(IPLAT,I2,1)+VPCFU(IPLAT,I1,1))/
c    .               (1.0-VPCFU(IPLAT,IPCFU,1))
c             ENDIF
c            ENDDO
c            VPCAU(IPLAT,IPCFU)=0.5*VAUXX
cc           VPCAU(IPLAT,IPCFU)=VAUXX-VPCAU(IPLAT,1)    ! old version
c           ENDDO
c          ENDIF                             ! ntypc.eq.1
C
c          VPLAT(IPLAT,1)=VPCFU(IPLAT,1,2)      ! Ts
c          VPLAT(IPLAT,2)=VPCFU(IPLAT,NPCFU,2)  ! Tl
c          VPLAT(IPLAT,3)=VPCFU(IPLAT,NPCFU,3)  ! last L
c          IX=IX+3*NPCFU+NPCPC+3                ! 3=npcfu,nlate,ntypc
c         ENDIF
c        ENDIF
C
c        IF(IPCMO.EQ.3) THEN
c         IA=I4B+1
c         IB=IA+1
c         IC=IB+1
c         VPLAT(IPLAT,1)=PROPST(IA)
c         VPLAT(IPLAT,2)=PROPST(IB)
c         VPLAT(IPLAT,3)=PROPST(IC)
c         IX=IX+NPCPC+3                       ! 3=Ts,Tl,L
c        ENDIF
C
c        IF(IPCMO.EQ.4) THEN
c         call runendt('error in ideprot')
c        ENDIF
C
c        IF(IPCMO.EQ.5) THEN
c         IA=I4B+1
c         IB=IA+1
c         IC=IB+1
c         ID=IC+1
c         IE=ID+1
c         VPLAT(IPLAT,1)=PROPST(IA)    ! Ts
c         VPLAT(IPLAT,7)=PROPST(IB)    ! Te
c         VPLAT(IPLAT,2)=PROPST(IC)    ! Tl
c         VPLAT(IPLAT,6)=PROPST(ID)    ! Tm
c         VPLAT(IPLAT,3)=PROPST(IE)    ! L
c         IX=IX+NPCPC+5                       ! 5=Ts,Te,Tl,Tm,L
c        ENDIF
C
c        IF(IPCMO.EQ.6) THEN
c         IA=I4B+1
c         IB=IA+1
c         IC=IB+1
c         ID=IC+1
c         VPLAT(IPLAT,1)=PROPST(IA)    ! Ts
c         VPLAT(IPLAT,2)=PROPST(IB)    ! Tl
c         VPLAT(IPLAT,6)=PROPST(IC)    ! Tm
c         VPLAT(IPLAT,3)=PROPST(ID)    ! L
c         IX=IX+NPCPC+4                       ! 4=Ts,Tl,Tm,L
c        ENDIF
C
c        IF(IPCMO.EQ.7) THEN
c         call runendt('error in ideprot')
c        ENDIF
C
c       ELSE            ! ipcfo=1
c        CALL RUNENDT('IDEPROT: ERROR WITH IPCFO=1 >> NOT POSSIBLE')
c       ENDIF           ! ipcfo.eq.0
c      ENDDO            ! do iplat=1,nplat
C
c      I5=IX
c      DO IDIMET=1,NDIMET
c       IA=I5-1+2*IDIMET
c       IB=IA+1
c       NSUBD(IDIMET,1)=INT(PROPST(IA))
c       NSUBD(IDIMET,2)=INT(PROPST(IB))
C
C**** CHECK NSUBD
C
c       IF(NDIMET.EQ.3.AND.NNODLT.EQ.4) THEN  ! linear tetrahedra
c        NSUBD(IDIMET,1)=1
c        NSUBD(IDIMET,2)=1
c       ENDIF
c      ENDDO
C
C**** POROSITY CRITERIA
C
c      I6=I5+2*NDIMET+1
c      KPOROT=INT(PROPST(I6))
      ELSE              ! nplat=0
       I6=I4
       KPOROS=0
      ENDIF             ! nplat.ne.0
C
C**** FREE ENERGY MODEL
C
      I11=I6+1
      IFRENTS=INT(PROPSS(I11))
C
C**** HEAT FLUX MODEL
C
      I12=I11+1
      IFREKTS=INT(PROPSS(I12))
C
      RETURN
      END
