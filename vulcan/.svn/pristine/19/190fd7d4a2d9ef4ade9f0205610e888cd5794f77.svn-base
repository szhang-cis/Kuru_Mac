      SUBROUTINE IDEPROTO(PROPST)
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
      INCLUDE 'prob_omt.f'
      INCLUDE 'inte_omt.f'
      INCLUDE 'auxl_omt.f'
C
      DIMENSION PROPST(*)
C
      I1=3                            ! ISOTRT=PROPST(1),IFILL=PROPST(2)
      NDENS=INT(PROPST(I1))
C
      DO IDENS=1,NDENS
       IA=I1-1+2*IDENS
       IB=IA+1
       VDENS(IDENS,1)=PROPST(IA)
       VDENS(IDENS,2)=PROPST(IB)
      ENDDO
C
      I2=I1+2*NDENS+1
      NCAPA=INT(PROPST(I2))
C
      DO ICAPA=1,NCAPA
       IA=I2-1+2*ICAPA
       IB=IA+1
       VCAPA(ICAPA,1)=PROPST(IA)
       VCAPA(ICAPA,2)=PROPST(IB)
      ENDDO
C
      I3=I2+2*NCAPA     ! note: no +1
      NCOND=0
      DO IDIMET=1,NDIMET
       I3=I3+2*NCOND+1
       NCOND=INT(PROPST(I3))
C
       DO ICOND=1,NCOND
        IA=I3-1+2*ICOND
        IB=IA+1
        VCOND(ICOND,1,IDIMET)=PROPST(IA)
        VCOND(ICOND,2,IDIMET)=PROPST(IB)
       ENDDO
      ENDDO             ! idimet=1,ndimet
C
      I4=I3+2*NCOND+1
      NPLAT=INT(PROPST(I4))
C
      IF(NPLAT.NE.0) THEN
       I4=I4+1
       LINPC=INT(PROPST(I4))
C
       IX=I4
       DO IPLAT=1,NPLAT
        I4A=IX+1
        IPCFO=INT(PROPST(I4A))
        VPLAT(IPLAT,4)=PROPST(I4A)
        I4B=I4A+1
        IPCMO=INT(PROPST(I4B))
        VPLAT(IPLAT,5)=PROPST(I4B)
C
        NPCPC=2                               ! 2=ipcfo,ipcmo
        IF(IPCFO.EQ.0) THEN
C
         IF(IPCMO.EQ.0.OR.IPCMO.EQ.-1.OR.IPCMO.EQ.-2) THEN
          IA=I4B+1
          IB=IA+1
          IC=IB+1
          VPLAT(IPLAT,1)=PROPST(IA)
          VPLAT(IPLAT,2)=PROPST(IB)
          VPLAT(IPLAT,3)=PROPST(IC)
          IX=IX+NPCPC+3                       ! 2=ipcmo,ipcfo, 3=Ts,Tl,L
         ENDIF
C
         IF(IPCMO.EQ.1.OR.IPCMO.EQ.2) THEN
          I4C=I4B+1
          NPCFU=INT(PROPST(I4C))
          VPLAT(IPLAT,6)=PROPST(I4C)
          I4D=I4C+1
          NLATE=INT(PROPST(I4D))
          VPLAT(IPLAT,7)=PROPST(I4D)
          IF(NLATE.EQ.0) THEN                 ! constant L
           DO IPCFU=1,NPCFU
            IA=I4D-1+2*IPCFU
            IB=IA+1
            VPCFU(IPLAT,IPCFU,1)=PROPST(IA)
            VPCFU(IPLAT,IPCFU,2)=PROPST(IB)
           ENDDO
           VPLAT(IPLAT,1)=VPCFU(IPLAT,1,2)      ! Ts
           VPLAT(IPLAT,2)=VPCFU(IPLAT,NPCFU,2)  ! Tl
           I4E=IB+1
           VPLAT(IPLAT,3)=PROPST(I4E)
           IX=IX+2*NPCFU+NPCPC+3              ! 3=npcfu,nlate,L
          ELSE                                ! temp.-dep. L
           I4D=I4D+1
           NTYPC=INT(PROPST(I4D))             ! secant or tangent L
           VPLAT(IPLAT,8)=PROPST(I4D)
           DO IPCFU=1,NPCFU
            IA=I4D-2+3*IPCFU
            IB=IA+1
            IC=IB+1
            VPCFU(IPLAT,IPCFU,1)=PROPST(IA)
            VPCFU(IPLAT,IPCFU,2)=PROPST(IB)
            VPCFU(IPLAT,IPCFU,3)=PROPST(IC)
           ENDDO
C
           IF(NTYPC.EQ.1) THEN                ! tangent L
            VPCAU(IPLAT,1)=VPCFU(IPLAT,1,3)   ! computes L^s
            DO IPCFU=2,NPCFU
             VAUXX=0.0D+00
             DO JPCFU=2,IPCFU
              I1=JPCFU-1
              I2=JPCFU
              IF(IPCMO.EQ.1) THEN
               VAUXX=VAUXX+(VPCFU(IPLAT,I2,3)+VPCFU(IPLAT,I1,3))*
     .               (VPCFU(IPLAT,I2,1)-VPCFU(IPLAT,I1,1))/
     .               VPCFU(IPLAT,IPCFU,1)
              ELSE                            ! considers f_pc
               VAUXX=VAUXX+(VPCFU(IPLAT,I2,3)+VPCFU(IPLAT,I1,3))*
     .               (-VPCFU(IPLAT,I2,1)+VPCFU(IPLAT,I1,1))/
     .               (1.0-VPCFU(IPLAT,IPCFU,1))
              ENDIF
             ENDDO
             VPCAU(IPLAT,IPCFU)=0.5*VAUXX
            ENDDO
           ENDIF                             ! ntypc.eq.1
C
           VPLAT(IPLAT,1)=VPCFU(IPLAT,1,2)      ! Ts
           VPLAT(IPLAT,2)=VPCFU(IPLAT,NPCFU,2)  ! Tl
           VPLAT(IPLAT,3)=VPCFU(IPLAT,NPCFU,3)  ! last L
           IX=IX+3*NPCFU+NPCPC+3                ! 3=npcfu,nlate,ntypc
          ENDIF
         ENDIF
C
         IF(IPCMO.EQ.3) THEN
          IA=I4B+1
          IB=IA+1
          IC=IB+1
          VPLAT(IPLAT,1)=PROPST(IA)
          VPLAT(IPLAT,2)=PROPST(IB)
          VPLAT(IPLAT,3)=PROPST(IC)
          IX=IX+NPCPC+3                       ! 3=Ts,Tl,L
         ENDIF
C
         IF(IPCMO.EQ.4) THEN
          call runendt('error in ideprot')
         ENDIF
C
         IF(IPCMO.EQ.5) THEN
          IA=I4B+1
          IB=IA+1
          IC=IB+1
          ID=IC+1
          IE=ID+1
          VPLAT(IPLAT,1)=PROPST(IA)    ! Ts
          VPLAT(IPLAT,7)=PROPST(IB)    ! Te
          VPLAT(IPLAT,2)=PROPST(IC)    ! Tl
          VPLAT(IPLAT,6)=PROPST(ID)    ! Tm
          VPLAT(IPLAT,3)=PROPST(IE)    ! L
          IX=IX+NPCPC+5                       ! 5=Ts,Te,Tl,Tm,L
         ENDIF
C
         IF(IPCMO.EQ.6) THEN
          IA=I4B+1
          IB=IA+1
          IC=IB+1
          ID=IC+1
          VPLAT(IPLAT,1)=PROPST(IA)    ! Ts
          VPLAT(IPLAT,2)=PROPST(IB)    ! Tl
          VPLAT(IPLAT,6)=PROPST(IC)    ! Tm
          VPLAT(IPLAT,3)=PROPST(ID)    ! L
          IX=IX+NPCPC+4                       ! 4=Ts,Tl,Tm,L
         ENDIF
C
         IF(IPCMO.EQ.7) THEN
          call runendt('error in ideprot')
         ENDIF
C
        ELSE            ! ipcfo=1
         CALL RUNENDT('IDEPROT: ERROR WITH IPCFO=1 >> NOT POSSIBLE')
        ENDIF           ! ipcfo.eq.0
       ENDDO            ! do iplat=1,nplat
C
       I5=IX
       DO IDIMET=1,NDIMET
        IA=I5-1+2*IDIMET
        IB=IA+1
        NSUBD(IDIMET,1)=INT(PROPST(IA))
        NSUBD(IDIMET,2)=INT(PROPST(IB))
C
C**** CHECK NSUBD
C
        IF(NDIMET.EQ.3.AND.NNODLT.EQ.4) THEN  ! linear tetrahedra
         NSUBD(IDIMET,1)=1
         NSUBD(IDIMET,2)=1
        ENDIF
       ENDDO
C
C**** POROSITY CRITERIA
C
       I6=I5+2*NDIMET+1
       KPOROT=INT(PROPST(I6))
      ELSE              ! nplat=0
       I6=I4
       KPOROT=0
      ENDIF             ! nplat.ne.0
C
C**** FREE ENERGY MODEL
C
      I11=I6+1
      IFRENT=INT(PROPST(I11))
C
C**** HEAT FLUX MODEL
C
      I12=I11+1
      IFREKT=INT(PROPST(I12))
C
      RETURN
      END
