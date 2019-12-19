      SUBROUTINE INMICRO11(ALPHAM,INUPC,IPLAT,INUPM)
C***********************************************************************
C
C**** THIS ROUTINE INITIALITES MICROSTRUCTURAL PARAMETERS FOR MODEL 11
C
C***********************************************************************
C
C     Index of variables
C
C     ALPHAM= array of microstructural (microscopical) variables
C
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C**** COUPLING VARIABLES (thermal-microstructural)
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
C
C**** TRANSFERS "VPLAT" TO MICROSCOPICAL VARIABLES
C
      IX=0
      NCOPC=IPLAC(INUPM,1,1)
      IF(NCOPC.EQ.0) THEN
       CA=VPLAT(IPLAT, 6)
       SI=VPLAT(IPLAT, 7)
       CU=VPLAT(IPLAT, 8)
       QM=VPLAT(IPLAT, 9)
       CR=VPLAT(IPLAT,10)
       AN=VPLAT(IPLAT,11)
       IX=6
      ENDIF
C
      TEFERI=VPLAT(IPLAT, 6+IX)          ! reference ferrite T
      RFER0 =VPLAT(IPLAT,11+IX)          ! initial ferrite grain radius
      TEPERI=VPLAT(IPLAT,12+IX)          ! reference pearlite T
C
      IFPCDT=INT(VPLAT(IPLAT,20+IX))     ! index for temp. derivative
      IAFLOJ=INT(VPLAT(IPLAT,21+IX))     ! index for fraction correction
C
      IV=21+IX                   ! number of VPLAT defined in input data
C
C**** DEALS WITH COUPLING BETWEEN MICROSTRUCTURAL MODELS
C
      IPLX2=0
      IF(NCOPC.NE.0) THEN       ! initial condition from other model
       IPLX1=IPLAC(INUPM,2,1)   ! ph-ch that provides initial condition
       IPLX2=IPLAC(INUPM,2,2)   ! model that provides initial condition
      ENDIF                     ! ncopc.ne.0
C
      IF(IPLX2.EQ.2) THEN
       CA=VPLAT(IPLX1,6)                 ! carbon
       SI=VPLAT(IPLX1,7)                 ! silicon
       CU=VPLAT(IPLX1,9)                 ! copper
       QM=VPLAT(IPLX1,10)                ! manganese
       CR=VPLAT(IPLX1,12)                ! chromium
       AN=VPLAT(IPLX1,13)                ! nickel
       IVX=30
       FGREUT=VPLAT(IPLX1,IVX+8)
      ENDIF
C
      TEFER=TEFERI+31.5D0*SI-7.7D0*CU-18.7D0*QM-10.7D0*CR-26.0D0*AN
C     TEPER=TEPERI+30.07D0*SI-1.98D0*SI*SI-10.7D0*CU-13.7*QM+24.3D0*CR-
C    .      12.0D0*AN
      TEPER=TEPERI+21.6D0*SI+0.023D0*SI*SI-21.0D0*CU-25.0*QM+13.0D0*CR-
     .      33.0D0*AN                                               ! AU
      TALFA=TEFERI+18.4D0*SI+2.0D0*SI*SI-14.0D0*CU-45.0D0*QM-24.0D0*CR-
     .      27.5D0*AN
      TCURIE=1043.D0-10.0D0*SI
C
C     CEUD=18.76D0-0.04112D0*TEFER+0.0000226*TEFER*TEFER+0.125D0*SI
      CEUD=18.76D0-0.04112D0*TALFA+0.0000226*TALFA*TALFA+0.125D0*SI ! AU
      CAMAX=0.067D0-5.0D-8*TALFA*TALFA-2.8D-5*TALFA+(1.2D-8*TALFA*TALFA-
     .                                                        4.5D-3)*SI
C
      VPLAT(IPLAT,IV+1)=CA
      VPLAT(IPLAT,IV+2)=SI
      VPLAT(IPLAT,IV+3)=QM
C
      VPLAT(IPLAT,IV+4)=FGREUT
C
      VPLAT(IPLAT,IV+5)=TEFER
      VPLAT(IPLAT,IV+6)=TEPER
      VPLAT(IPLAT,IV+7)=TALFA
      VPLAT(IPLAT,IV+8)=TCURIE
      VPLAT(IPLAT,IV+9)=CEUD
      VPLAT(IPLAT,IV+10)=CAMAX
C
C**** INITIALITES MICROSTRUCTURAL VARIABLES ("ALPHAM")
C
      IN=INUPC
      DO INU=1,7
       ALPHAM(IN+INU)=0.0D0
      ENDDO
      ALPHAM(IN+3)=RFER0
C
C**** INCREMENTS "ALPHAM" INDEX
C
      INUPC=INUPC+7                    ! less than NBASES; see pointes.f
C
      RETURN
      END
