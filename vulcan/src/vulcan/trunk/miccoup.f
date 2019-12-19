      SUBROUTINE MICCOUP
C***********************************************************************
C
C**** THIS ROUTINE DEALS WITH COUPLING BETWEEN MICROSTRUCTURAL MODELS
C
C     Notes:
C     Implemented coupling (may/07): models 4, 5 or 9 (solidification)
C     provide initial condition (ic) to model 12 (SG eutectoid transf.).
C     Implemented coupling (nov/12): model 2 (solidification)
C     provides initial condition to model 11 (gray eutectoid transf.).
C     Implemented coupling (may/14): model 14 (reaustenization)
C     provides initial condition to model 6 version 3 (austempered
C     transf.).
C     To increase the third dimension of IPLAC, increase NPLAC in
C     setdatd.f (see also nued_om.f)
C     For changes in INUPC, aldo see pointes.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nued_om.f'
C
C**** MICROSCOPICAL VARIABLES
C
      INCLUDE 'auxl_oms.f'
C
      DO INUPM=1,NNUPM
       NCOPC=IPLAC(INUPM,1,1)             ! total number of ph-ch for ic
       IF(NCOPC.NE.0) THEN
        DO ICOPC=1,NCOPC
         INUPX=IPLAC(INUPM,ICOPC+1,1)     ! ph-ch numb. for ic
         INUPC=0
C                                         ! all phase-change models
         IF(INUPX.GT.1) THEN
          DO JNUPX=1,INUPX-1
           IF(IPLLLS(JNUPX).EQ.2)
     .      INUPC=INUPC+14                ! see inmicro2.f & micros2.f
           IF(IPLLLS(JNUPX).EQ.4)
     .      INUPC=INUPC+ 4+2*NNUM4S+2     ! see inmicro4.f & micros4.f
           IF(IPLLLS(JNUPX).EQ.5)
     .      INUPC=INUPC+ 4+3*NNUM4S+2     ! see inmicro5.f & micros5.f
           IF(IPLLLS(JNUPX).EQ.6)
     .      INUPC=INUPC+6                 ! see inmicro6.f & micros6.f
           IF(IPLLLS(JNUPX).EQ.9)
     .      INUPC=INUPC+26+5*NNUM4S+5     ! see inmicro9.f & micros9.f
           IF(IPLLLS(JNUPX).EQ.11)
     .      INUPC=INUPC+7                 ! see inmicro11.f & micros11.f
           IF(IPLLLS(JNUPX).EQ.12)
     .      INUPC=INUPC+4+6*NNUM4S+11     ! see inmicro12.f & micros12.f
           IF(IPLLLS(JNUPX).EQ.14)
     .      INUPC=INUPC+14                ! see inmicro14.f & micros14.f
          ENDDO
         ENDIF
C                                         ! only ic models
         IF(IPLLLS(INUPX).EQ.2) THEN            ! IPCMO=2
          IPLAC(INUPM,ICOPC+1, 2)=2             ! Gray cast iron model
          IPLAC(INUPM,ICOPC+1, 3)=INUPC+3       ! graph. eutec. fraction
          IPLAC(INUPM,ICOPC+1, 4)=INUPC+9       ! austenite fraction
          IPLAC(INUPM,ICOPC+1, 5)=INUPC+10      ! graphite fraction
          IPLAC(INUPM,ICOPC+1, 6)=INUPC+11      ! lamellar spacing
          IPLAC(INUPM,ICOPC+1, 7)=INUPC+12      ! lamellae thickness
          IPLAC(INUPM,ICOPC+1, 8)=INUPC+1       ! liquid fraction
          IPLAC(INUPM,ICOPC+1, 9)=INUPC+6       ! white eutec. fraction
         ENDIF
         IF(IPLLLS(INUPX).EQ.4) THEN            ! IPCMO=4
          IPLAC(INUPM,ICOPC+1, 2)=4             ! Boeri's model
          IPLAC(INUPM,ICOPC+1, 3)=INUPC+2       ! austenite fraction
          IPLAC(INUPM,ICOPC+1, 4)=0             ! nothing
          IPLAC(INUPM,ICOPC+1, 5)=INUPC+3       ! graphite fraction
          IPLAC(INUPM,ICOPC+1, 6)=INUPC+5       ! first graphite gr. N
          IPLAC(INUPM,ICOPC+1, 7)=INUPC+6+      ! first graphite gr. R
     .                            (NNUM4S-1)
          IPLAC(INUPM,ICOPC+1, 8)=0             ! nothing
          IPLAC(INUPM,ICOPC+1, 9)=0             ! nothing
          IPLAC(INUPM,ICOPC+1,10)=0             ! nothing
          IPLAC(INUPM,ICOPC+1,11)=INUPC+7+      ! graphite nucl. index
     .                            2*(NNUM4S-1)
         ENDIF
         IF(IPLLLS(INUPX).EQ.5) THEN            ! IPCMO=5
          IPLAC(INUPM,ICOPC+1, 2)=5             ! Su's model
          IPLAC(INUPM,ICOPC+1, 3)=INUPC+2       ! austenite fraction
          IPLAC(INUPM,ICOPC+1, 4)=INUPC+7+      ! first austenite gr. R
     .                            2*(NNUM4S-1)
          IPLAC(INUPM,ICOPC+1, 5)=INUPC+3       ! graphite fraction
          IPLAC(INUPM,ICOPC+1, 6)=INUPC+5       ! first graphite gr. N
          IPLAC(INUPM,ICOPC+1, 7)=INUPC+6+      ! first graphite gr. R
     .                            (NNUM4S-1)
          IPLAC(INUPM,ICOPC+1, 8)=0             ! nothing
          IPLAC(INUPM,ICOPC+1, 9)=0             ! nothing
          IPLAC(INUPM,ICOPC+1,10)=0             ! nothing
          IPLAC(INUPM,ICOPC+1,11)=INUPC+8+      ! graphite nucl. index
     .                            3*(NNUM4S-1)
         ENDIF
         IF(IPLLLS(INUPX).EQ.9) THEN            ! IPCMO=9
          IPLAC(INUPM,ICOPC+1, 2)=9             ! Patricia's model
          IPLAC(INUPM,ICOPC+1, 3)=INUPC+2       ! austenite fraction
          IPLAC(INUPM,ICOPC+1, 4)=INUPC+5       ! austenite gr. R
          IPLAC(INUPM,ICOPC+1, 5)=INUPC+3       ! graphite fraction
          IPLAC(INUPM,ICOPC+1, 6)=INUPC+27      ! first graphite gr. Nz1
          IPLAC(INUPM,ICOPC+1, 7)=INUPC+28+     ! first graphite gr. Nz2
     .                            (NNUM4S-1)
          IPLAC(INUPM,ICOPC+1, 8)=INUPC+29+     ! first graphite gr. Nz3
     .                            2*(NNUM4S-1)
          IPLAC(INUPM,ICOPC+1, 9)=INUPC+30+     ! first graphite R RNU
     .                            3*(NNUM4S-1)
          IPLAC(INUPM,ICOPC+1,10)=INUPC+31+     ! first graphite R RNUZ1
     .                            4*(NNUM4S-1)
          IPLAC(INUPM,ICOPC+1,11)=INUPC+32+     ! graphite nucl. index
     .                            5*(NNUM4S-1)
          IPLAC(INUPM,ICOPC+1,12)=INUPC+11      ! eutectic cell R RADC
          IPLAC(INUPM,ICOPC+1,13)=INUPC+4       ! silicon content
          IPLAC(INUPM,ICOPC+1,14)=INUPC+14      ! primary & secondary
C                                               ! arm spacing PDAS,SDAS
          IPLAC(INUPM,ICOPC+1,15)=INUPC+17      ! copper content
          IPLAC(INUPM,ICOPC+1,16)=INUPC+18      ! manganese content
          IPLAC(INUPM,ICOPC+1,17)=INUPC+19      ! magnesium content
          IPLAC(INUPM,ICOPC+1,18)=INUPC+20      ! niobium content
          IPLAC(INUPM,ICOPC+1,19)=INUPC+21      ! tin content
          IPLAC(INUPM,ICOPC+1,20)=INUPC+22      ! chromium content
          IPLAC(INUPM,ICOPC+1,21)=INUPC+23      ! molybdenum content
          IPLAC(INUPM,ICOPC+1,22)=INUPC+24      ! nickel content
          IPLAC(INUPM,ICOPC+1,23)= INUPC+9      ! cpro
         ENDIF
         IF(IPLLLS(INUPX).EQ.14) THEN
          IPLAC(INUPM,ICOPC+1,2)=14             ! eutectoid inverse mod.
          IPLAC(INUPM,ICOPC+1,3)=INUPC+2        ! vfGr
          IPLAC(INUPM,ICOPC+1,4)=INUPC+3        ! vfFs
          IPLAC(INUPM,ICOPC+1,5)=INUPC+4        ! vfAs
          IPLAC(INUPM,ICOPC+1,6)=INUPC+5        ! vfP
          IPLAC(INUPM,ICOPC+1,7)=INUPC+6        ! vfCm
          IPLAC(INUPM,ICOPC+1,8)=INUPC+7        ! vfFm
          IPLAC(INUPM,ICOPC+1,9)=INUPC+8        ! vfAm
          IPLAC(INUPM,ICOPC+1,10)=INUPC+9       ! vfA
          IPLAC(INUPM,ICOPC+1,11)=INUPC+12      ! cA (micros14.f)
C
          IPLAC(INUPM,ICOPC+1,13)=6             ! c
          IPLAC(INUPM,ICOPC+1,14)=7             ! si
          IPLAC(INUPM,ICOPC+1,15)=8             ! mn
          IPLAC(INUPM,ICOPC+1,16)=9             ! cr
          IPLAC(INUPM,ICOPC+1,17)=10            ! ni
          IPLAC(INUPM,ICOPC+1,18)=11            ! cu
          IPLAC(INUPM,ICOPC+1,19)=12            ! mo (inmicro14.f)
         ENDIF
        ENDDO                ! icopc=1,ncopc
       ENDIF                 ! ncopc.ne.0
      ENDDO                  ! inupm=1,nnupm
C
      RETURN
      END
