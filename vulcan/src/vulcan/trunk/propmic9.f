      SUBROUTINE PROPMIC9(ITAPET,PROPSS,NPOI2T,INUPM)
C***********************************************************************
C     
C**** THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
C     MODEL 9
C     
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C     
C**** MICROSCOPICAL VARIABLES
C     
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
C     
      DIMENSION PROPSS(*)       ! inicializado en algunos de los include
C     
C**** COMPOSITION (INITIAL)
C     
C**   carbon
      NPRINT= 0
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'COMPO') THEN
         WRITE(LURESS,810)
         WRITE(LURESS,809)
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'CARBO') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de C (carbono)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO CARBON CARD')
         ENDIF
C     control on carbon weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG CARBON COMPOSITION')
C     control on carbon (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG CARBON PART. COEF.')
C     control on carbon (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG CARBON PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID CARBON 
     .        INDEX PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID CARBON 
     .        INDEX PART. COEF.')
         WRITE(LURESS,811) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2), INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   silicon         
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'SILIC') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Si (carbono)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Si (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Si (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO SILICON CARD')
         ENDIF
C     control on carbon weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG SILICON COMPOSITION')
C     control on carbon (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG SILICON PART. COEF.')
C     control on carbon (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG SILICON PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID SILICON 
     .        INDEX PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID SILICON 
     .        INDEX PART. COEF.')
         WRITE(LURESS,812) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   phosporus
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'PHOSP') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de P (phosphorus)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del P (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del P (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO PHOSPHORUS CARD')
         ENDIF
C     control on carbon weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG PHOSPHORUS COMPOSITION')
C     control on carbon (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG PHOSPHORUS PART. COEF.')
C     control on carbon (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG PHOSPHORUS PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID PHOSPHORUS 
     .        INDEX PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID PHOSPHORUS 
     .        INDEX PART. COEF.')
         WRITE(LURESS,813) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'COPPE') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Cu (copper)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Cu (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Cu (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO COPPER CARD')
         ENDIF
C     control on carbon weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG COPPER COMPOSITION')
C     control on carbon (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG COPPER PART. COEF.')
C     control on carbon (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP COPPER PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID COPPER 
     .        INDEX PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID COPPER 
     .        INDEX PART. COEF.')
         WRITE(LURESS,814) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   manganese
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'MANGA') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Mn (Manganese)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Mn (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Mn (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO MANGANESE CARD')
         ENDIF
C     control on carbon weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MANGANESE COMPOSITION')
C     control on carbon (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MANGANESE PART. COEF.')
C     control on carbon (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP MANGANESE PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID MANGANESE 
     .        INDEX PART. COEF.')
C     control on coeficiente de particion
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID MANGANESE 
     .        INDEX PART. COEF.')
         WRITE(LURESS,815) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   magnesium
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'MAGNE') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Mg (magnesium)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Mg (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Mg (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO MAGNESIUM CARD')
         ENDIF
C     control on magnesium weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MAGNESIUM COMPOSITION')
C     control on magnesium (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MAGNESIUM PART. COEF.')
C     control on magnesium (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP MAGNESIUM PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID MAGNESIUM 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID MAGNESIUM 
     .        INDEX PART. COEF.')
         WRITE(LURESS,816) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   niobium
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'NIOBI') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Ni (niobium)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Ni (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del MNi (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO NIOBIUM CARD')
         ENDIF
C     control on niobium weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG NIOBIUM COMPOSITION')
C     control on niobium (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG NIOBIUM PART. COEF.')
C     control on niobium (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP NIOBIUM PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID NIOBIUM 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID NIOBIUM 
     .        INDEX PART. COEF.')
         WRITE(LURESS,817) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   tin == estano
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'TIN_C') THEN ! pearlite promoter in SGI
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Tin (estano)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Tin (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Tin (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO TIN CARD')
         ENDIF
C     control on tin weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG TIN COMPOSITION')
C     control on tin (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG TIN PART. COEF.')
C     control on tin (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP TIN PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID TIN 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID TIN 
     .        INDEX PART. COEF.')
         WRITE(LURESS,818) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)),
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C**   chromium
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'CHROM') THEN ! carbide promoter in SGI
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Cr (cromo)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Cr (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Cr (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO CHROMIUM CARD')
         ENDIF
C     control on chromium weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG CHROMIUM COMPOSITION')
C     control on chromium (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG CHROMIUM PART. COEF.')
C     control on chromium (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP CHROMIUM PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID CHROMIUM 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID CHROMIUM 
     .        INDEX PART. COEF.')
         WRITE(LURESS,819) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C**   molybdenum
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'MOLYB') THEN ! carbide promoter in SGI
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Mo (molibdeno)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Mo (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Mo (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO MOLYBDENUM CARD')
         ENDIF
C     control on molybdenum weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MOLYBDENUM COMPOSITION')
C     control on molybdenum (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG MOLYBDENUM PART. COEF.')
C     control on molybdenum (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP MOLYBDENUM PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID MOLYBDENUM 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID MOLYBDENUM 
     .        INDEX PART. COEF.')
         WRITE(LURESS,820) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
C     
C**   nickel
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'NICKE') THEN ! carbide promoter in SGI
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! percentaje en peso de Ni (niquel)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! coeficiente de particion del Ni (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! coeficiente de particion del Ni (solid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(4) ! modelo de particion (liquid state)
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(5) ! modelo de particion (solid state)
         ELSE
            CALL RUNENDS('PROPMIC9: NO NICKEL CARD')
         ENDIF
C     control on nickel weight porcent
         IF(PARAMS(1).LT.0.0D0.OR.PARAMS(1).GT.6.3D0)
     .        CALL RUNENDS('PROPMIC9: WRONG NICKEL COMPOSITION')
C     control on nickel (liquid) partition coeficient
         IF(PARAMS(2).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONG NICKEL PART. COEF.')
C     control on nickel (solid) partition coeficient
         IF(PARAMS(3).LT.0.0D0.OR.PARAMS(2).GT.2.0D0)
     .        CALL RUNENDS('PROPMIC9: WRONGP NICKEL PART. COEF.')
C     control on partition coefficient in Liquid state
         IF(INT(PARAMS(4)).LT.1.OR.INT(PARAMS(4)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG LIQUID NICKEL 
     .        INDEX PART. COEF.')
C     control on partition coefficient in Solid state
         IF(INT(PARAMS(5)).LT.1.OR.INT(PARAMS(5)).GT.6)
     .        CALL RUNENDS('PROPMIC9: WRONG SOLID NICKEL 
     .        INDEX PART. COEF.')
         WRITE(LURESS,821) PROPSS(NPOI2T-4),PROPSS(NPOI2T-3),
     .        PROPSS(NPOI2T-2),INT(PROPSS(NPOI2T-1)), 
     .        INT(PROPSS(NPOI2T))
         WRITE(LURESS,809)
      ELSE
         CALL RUNENDS('PROPMIC9: NO COMPOSITION CARD')
      ENDIF
C     
C**** PRIMARY GRAPHITE SOLIDIFICATION (HIPOEUTECTIC COMPOSITION)
C     (not implemented yet, only eutectic composition)
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'HIPOE') THEN
         WRITE(LURESS,825)
         WRITE(LURESS,809)
      ENDIF
C     
C**** PRIMARY AUSTENITE SOLIDIFICATION  (HIPEREUTECTIC COMPOSITION)
C     (not implemented yet, only eutectic composition)
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'HIPER') THEN
         WRITE(LURESS,830)
         WRITE(LURESS,809)
      ENDIF
C     
C**** EUTECTIC SOLIDIFICATION
C     
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'EUTEC') THEN
         WRITE(LURESS,840)
         WRITE(LURESS,809)
C     
C**** NUCLEATION PROPERTIES
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'NUCLE'.AND.WORDSS(2).EQ.'MODEL') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! nucleation model index
            INUCMX= INT(PARAMS(1))
            WRITE(LURESS,841) INUCMX
            WRITE(LURESS,809)
            IF(INUCMX.EQ.1.OR.INUCMX.EQ.2.OR.INUCMX.EQ.3
     .           .OR.INUCMX.EQ.4) THEN ! Su & Boeri & Rappaz & Stefanescu
               CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(1) ! coefficient A
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(2) ! coefficient b
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(3) ! coefficient c
               WRITE(LURESS,842) PROPSS(NPOI2T-2),PROPSS(NPOI2T-1),
     .              PROPSS(NPOI2T)
               WRITE(LURESS,809)
            ELSE
               CALL RUNENDS('PROPMIC9: INCORRECT NUCLEATION MODEL 
     .              NUMBER')
            ENDIF
         ELSE
            CALL RUNENDS('PROPMIC9: NO NUCLEATION OR MODEL CARD')
         ENDIF
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'ARRES') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! nucleation arrest criterion index
            INUCAX= INT(PARAMS(1))
            WRITE(LURESS,843) INUCAX
            WRITE(LURESS,809)
            IF(INUCAX.NE.1.AND.INUCAX.NE.2) ! not Trecal & Tmin
     .           CALL RUNENDS('PROPMIC9: INCORRECT ARREST CRITERION 
     .           NUMBER')
         ELSE
            CALL RUNENDS('PROPMIC9: NO ARREST CRITERION CARD')
         ENDIF
C     
C**** GROWTH PROPERTIES
C     
C     Note: diffusion coef. of C in the liquid & 
C     auste. shell radius/graphite radius NOT USED NOW !!
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
C     GRAPHITE
         IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'GRAPH'.AND.
     .        WORDSS(3).EQ.'MODEL')THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! growth model index
            IGRGMX= INT(PARAMS(1))
            WRITE(LURESS,844) IGRGMX
            WRITE(LURESS,809)
            IF(IGRGMX.EQ.1.OR.IGRGMX.EQ.2) THEN ! Patri's or Zener
               CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(1) ! diffusion coef. of C in liq.
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(2) ! initial r of aust. shell
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(3) ! initial r of graphite
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(4) ! relation ra/rg
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(5) ! austenite density
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(6) ! graphite density
               WRITE(LURESS,845) PROPSS(NPOI2T-5),PROPSS(NPOI2T-4),
     .              PROPSS(NPOI2T-3), PROPSS(NPOI2T-2), 
     .              PROPSS(NPOI2T-1), PROPSS(NPOI2T)
               WRITE(LURESS,809)
            ELSE
               CALL RUNENDS('PROPMIC9: INCORRECT GRAPHITE GROWTH MODEL')
            ENDIF
         ELSE
            CALL RUNENDS('PROPMIC9: NOT GRAPHITE GROWTH CARD')
         ENDIF
C     AUSTENITE
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'AUSTE'.AND.
     .        WORDSS(3).EQ.'MODEL')THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! austenite growth model index
            IGRAMX= INT(PARAMS(1))
            WRITE(LURESS,846) IGRAMX
            WRITE(LURESS,809)
            IF(IGRAMX.EQ.1.OR.IGRAMX.EQ.2) THEN ! Rappaz or Stefanescu
               CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(1) ! diff. coef. of C in the aust.
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(2) ! liquidus pending line
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(3) ! part. coef. of C - PartCc
               NPOI2T= NPOI2T+1
               PROPSS(NPOI2T)= PARAMS(4) ! Gibbs-Thompson coefficient
               WRITE(LURESS,847) PROPSS(NPOI2T-3),PROPSS(NPOI2T-2),
     .              PROPSS(NPOI2T-1),PROPSS(NPOI2T)
               WRITE(LURESS,809)
            ELSE
               CALL RUNENDS('PROPMIC9: INCORRECT AUSTENITE GROWTH 
     .              MODEL')
            ENDIF
         ELSE
            CALL RUNENDS('PROPMIC9: NOT AUSTENITE GROWTH CARD')
         ENDIF
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'DELRN'.AND.
     .        WORDSS(3).EQ.'MODEL') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1)
            IDELRN= INT(PARAMS(1))
            WRITE(LURESS,848) IDELRN
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC9: INCORRECT DELRN GROWTH MODEL')
         ENDIF
C     sec. dendrite arm spacing
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'GROWT'.AND.WORDSS(2).EQ.'ISDAS'.AND.
     .        WORDSS(3).EQ.'MODEL') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1)
            WRITE(LURESS,849) INT(PARAMS(1))
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC9: INCORRECT ISDAS GROWTH MODEL')
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC9: NO EUTECTIC SOLIDIFICATION CARD')
      ENDIF
C     
C**** COUPLING MICRO-MACRO
C     
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET) 
      IF(WORDSS(1).EQ.'COUPL'.AND.WORDSS(2).EQ.'IMICM'.AND.
     .     WORDSS(3).EQ.'MODEL') THEN ! Boeri, old or Smooth
         NPOI2T= NPOI2T+1
         PROPSS(NPOI2T)= PARAMS(1)
         IMICOUP= INT(PARAMS(1))
         WRITE(LURESS,860) IMICOUP
      ELSE
         CALL RUNENDS('PROPMIC9: NO CUOPLED MICRO-MARCRO CARD')
      ENDIF
C     
C**** MICROSTRUCTURE-DEPENDENT THERMAL PROPERTIES
C     
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'CONDU'.AND.WORDSS(2).EQ.'MODEL') THEN
         NPOI2T= NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! micro-dependent conductivity index
         IKMICX=INT(PARAMS(1))
         IF(IKMICX.EQ.1.OR.IKMICX.EQ.2.OR.IKMICX.EQ.3) THEN
            CALL LISTENS('PROPMIC4',NPRINT,ITAPET)
            PROPSS(NPOI2T)= 1.0D0 ! conductivity depends on microstructure
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1) ! solid conductivity
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(2) ! mushy conductivity
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(3) ! liquid conductivity
            WRITE(LURESS,871) PROPSS(NPOI2T-2),PROPSS(NPOI2T-1),
     .           PROPSS(NPOI2T)
         ELSE
            CALL RUNENDS('PROPMIC9: INCORRECT CONDUCTIVITY MODEL 
     .           NUMBER')
         ENDIF
         WRITE(LURESS,809)
      ELSE
         NPOI2T= NPOI2T+1
         PROPSS(NPOI2T)= 0.0D0  ! conductivity does not depend on micro.
         GO TO 1
      ENDIF
C     
C**** NUMERICAL STRATEGY PARAMETERS
C     
      CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
 1    IF(WORDSS(1).EQ.'NUMER') THEN
         WRITE(LURESS,880)
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'TEMPE') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1)
            IFPCDT= INT(PARAMS(1))
            IF(IFPCDT.LT.1.OR.IFPCDT.GT.5)
     .           CALL RUNENDS('PROPMIC9: WRONG IFPCDT VALUE')
            WRITE(LURESS,881) IFPCDT
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC9: NO TEMPERATURE DERIVATIVE 
     .           FORM OF FL')
         ENDIF
C     
         CALL LISTENS('PROPMIC9',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'FRACT') THEN
            NPOI2T= NPOI2T+1
            PROPSS(NPOI2T)= PARAMS(1)
            IAFLOJ= INT(PARAMS(1))
            WRITE(LURESS,882) IAFLOJ
            WRITE(LURESS,809)
         ELSE
            CALL RUNENDS('PROPMIC9: NO FRACTION CORRECTION CARD')
         ENDIF
      ELSE
         CALL RUNENDS('PROPMIC9: NO NUMERICAL STRATEGY CARD')
      ENDIF
C     
C**** OUTPUT (see pointes.f)
C     
      KPLA9S= 1
      IPLLLS(INUPM)= 9
C     
      RETURN
 809  FORMAT(/)
 810  FORMAT(3X,'COMPOSITION (INITIAL):')
 811  FORMAT(3X,' CARBON ELEMENT                        :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 812  FORMAT(3X,' SILICON ELEMENT                       :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 813  FORMAT(3X,'PHOSPHORUS  ELEMENT                    :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 814  FORMAT(3X,'COPPER  ELEMENT                        :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 815  FORMAT(3X,'MANGANESE  ELEMENT                     :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 816  FORMAT(3X,'MAGNESIUM  ELEMENT                     :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 817  FORMAT(3X,'NIOBIUM  ELEMENT                       :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 818  FORMAT(3X,'TIN  ELEMENT                           :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
c$$$     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',E15.6,/,
c$$$     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',E15.6)
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 819  FORMAT(3X,'CHROMIUM  ELEMENT                      :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 820  FORMAT(3X,'MOLYBDENUM  ELEMENT                           :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
 821  FORMAT(3X,'NICKEL ELEMENT                         :',/,
     .     3X,'  WT%                                    = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (LIQUID)          = ',E15.6,/,
     .     3X,'  PARTITION COEFICIENT (SOLID)           = ',E15.6,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (LIQUID)   : ',I5,/,
     .     3X,'  INDEX FOR MODEL SEGREGATION (SOLID)    : ',I5)
C     
 825  FORMAT(3X,'HIPO EUTECTIC SOLIDIFI PARAMETERS      :')
 830  FORMAT(3X,'HIPER EUTECTIC SOLIDIFI PARAMETERS     :')
C     
 840  FORMAT(3X,'EUTECTIC SOLIDIFICATION PARAMETERS     :')
 841  FORMAT(3X,'GRAPHITE NUCLEATION MODEL              = ', I5)
 842  FORMAT(3X,'NUCLEATION PARAMETERS                  :',/,
     .     3X,'  COEFFICIENT A (austenite)              = ',E15.6,/,
     .     3X,'  COEFFICIENT b (graphite)               = ',E15.6,/,
     .     3X,'  COEFFICIENT c (graphite)               = ',E15.6)
 843  FORMAT(3X,'ARREST CRITERION MODEL                 : ', I5)
C     
 844  FORMAT(3X,'GRAPHITE GROWTH MODEL                  : ', I5)
 845  FORMAT(3X,'GRAPHITE GROWTH PARAMETERS             :',/,
     .     3X,'  DIFFUSION COEF. OF C IN THE LIQUID     = ',E15.6,/,
     .     3X,'  INITIAL RADIUS OF AUSTENITE SHELL      = ',E15.6,/,
     .     3X,'  AUSTENITE SHELL RADIUS/GRAPHITE RADIUS = ',E15.6,/,
     .     3X,'  INITIAL RADIUS OF NODULE               = ',E15.6,/,
     .     3X,'  AUSTENITE DENSITY                      = ',E15.6,/,
     .     3X,'  GRAPHITE DENSITY                       = ',E15.6)
 846  FORMAT(3X,'AUSTENITE GROWTH MODEL                 : ', I5)
 847  FORMAT(3X,'DIFFUSION COEF. OF C IN THE AUSTENITE  = ',E15.6,/,
     .     3X,'  PENDING LIQUIDUS LINE                  = ',E15.6,/,
     .     3X,'  CARBON PARTITION COEFFICIENT           = ',E15.6,/,
     .     3X,'  GIBBS-THOMPSON COEFFICIENT             = ',E15.6)
 848  FORMAT(3X,'RADIUS ZONE 2 GROWTH MODEL             : ',I5)
 849  FORMAT(3X,'SDAS INDEX GROWTH MODEL                : ',I5)
C     
 860  FORMAT(3X,'CLOUPING MICRO-MACRO MODEL             : ',I5)
C     
 871  FORMAT(3X,'MICROSTRUCTURE-DEPENDENT CONDUCTIVITY  : ',/,
     .     3X,' SOLID CONDUCTIVITY                      = ',E15.6,/,
     .     3X,' MUSHY CONDUCTIVITY                      = ',E15.6,/,
     .     3X,' LIQUID CONDUCTIVITY                     = ',E15.6)
C     
 880  FORMAT(3X,'NUMERICAL STRATEGY PARAMETERS          :')
 881  FORMAT(3X,'TEMPERATURE DERIVATIVE FORM OF FL      = ',I5)
 882  FORMAT(3X,'FRACTION CORRECTION FORM               = ',I5)
      END

