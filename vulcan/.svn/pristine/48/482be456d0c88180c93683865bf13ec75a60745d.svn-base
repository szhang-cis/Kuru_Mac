      SUBROUTINE PROPMIC14(ITAPET,PROPSS,NPOI2T,INUPM)
c*******************************************************************
c     THIS ROUTINE READS THE MICROSTRUCTURAL MATERIAL PROPERTIES FOR
c     MODEL 14
c*******************************************************************
c
c
c
c
c
c
c*******************************************************************
c     autor: adrian boccardo
c     date: last write 03-2014
c*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'prob_oms.f'
      INCLUDE 'inte_oms.f'
      INCLUDE 'auxl_oms.f'
      INCLUDE 'inpo_oms.f'
      DIMENSION PROPSS(*)
      NPRINT=0
c     
c--   CHEMICAL COMPOSITION
c
      CALL LISTENS('PROPM14',NPRINT,ITAPET) ! lee una linea completa
      IF(WORDSS(1).NE.'CHEMI') THEN
         CALL RUNENDS('PROPM14: CHEMICAL CARD NOT FOUND')
      ENDIF
c     CARBON
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'CARBO') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! c
         WRITE(LURESS,830) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: CARBON CARD NOT FOUND')
      ENDIF
c     SILICON
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'SILIC') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! si
         WRITE(LURESS,831) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: SILICON CARD NOT FOUND')
      ENDIF
c     MANGANESE
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'MANGA') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! mn
         WRITE(LURESS,832) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: MANGANESE CARD NOT FOUND')
      ENDIF
c     CHROMIUM
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'CHROM') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! cr
         WRITE(LURESS,833) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: CHROMIUM CARD NOT FOUND')
      ENDIF
c     NICKEL
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'NICKE') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! ni
         WRITE(LURESS,834) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: MICKEL CARD NOT FOUND')
      ENDIF
c     CUPPER
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'CUPPE') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! cu
         WRITE(LURESS,835) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: CUPPER CARD NOT FOUND')
      ENDIF
c     MOLYBDENUM
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'MOLYB') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! mo
         WRITE(LURESS,836) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: MOLYBDENUM CARD NOT FOUND')
      ENDIF

c     
c--   initial microstructure
c     
      CALL LISTENS('PROPM14',NPRINT,ITAPET) ! lee una linea completa
      IF(WORDSS(1).NE.'INITI') THEN
         CALL RUNENDS('PROPM14: INI_MICROS CARD NOT FOUND')
      ENDIF
C     volume fraction graphite
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'GRAPH') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! vfGr
         WRITE(LURESS,840) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: VOL_FRAC_GRAPHITE CARD NOT FOUND')
      ENDIF
C     volume fraction ferrite
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'FERRI') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! vfF
         WRITE(LURESS,841) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: VOL_FRAC_FERRITE CARD NOT FOUND')
      ENDIF
C     volume fraction pearlite
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'PEARL') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) ! vfP
         WRITE(LURESS,842) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: VOL_FRAC_PEARLITE CARD NOT FOUND')
      ENDIF

c     
c--   microstructure_feature
c     
      CALL LISTENS('PROPMIC14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'MICRO') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) !sip
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(2) !grnden
         WRITE(LURESS,845) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPMIC14: MICRO_FEATURE CARD NOT FOUND')
      ENDIF

c     
c--   latent heat
c     
      CALL LISTENS('PROPMIC14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'LATEN') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) !lhf
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(2) !lhp
         WRITE(LURESS,847) PROPSS(NPOI2T-1),PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPMIC14: LATEN HEAT CARD NOT FOUND')
      ENDIF
      
c     
c--   numerical strategy parameters
c     
      CALL LISTENS('PROPM14',NPRINT,ITAPET) ! lee una linea completa
      IF(WORDSS(1).NE.'NUMER') THEN
         CALL RUNENDS('PROPM14: NUMERICAL STRATEGY CARD NOT FOUND')
      ENDIF
C     subdivition step
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'SUBDI') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) !subdivition step
         WRITE(LURESS,849) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: SUBDIVITION STEP CARD NOT FOUND')
      ENDIF
C     temperatute derivative of fpc
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'TEMPE') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) !fpc
         WRITE(LURESS,850) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: TEMPERATURE_DERIVATIVE CARD NOT FOUND')
      ENDIF

c     
c--   mechanicals parameters
c     
      CALL LISTENS('PROPM14',NPRINT,ITAPET) ! lee una linea completa
      IF(WORDSS(1).NE.'MECHA') THEN
         CALL RUNENDS('PROPM14: MECHANICALS PARAM CARD NOT FOUND')
      ENDIF
C-    version mechanical model
      CALL LISTENS('PROPM14',NPRINT,ITAPET)
      IF(WORDSS(1).EQ.'MECHA') THEN
         NPOI2T=NPOI2T+1
         PROPSS(NPOI2T)=PARAMS(1) !mechanical model
         WRITE(LURESS,855) PROPSS(NPOI2T)
         WRITE(LURESS,899)
      ELSE
         CALL RUNENDS('PROPM14: MECHANICAL MODEL CARD NOT FOUND')
      ENDIF
C-    mechanicals parameter for model 1
      IF (PARAMS(1).EQ.1) THEN
C     secant volumetric change
         CALL LISTENS('PROPM14',NPRINT,ITAPET)
         IF(WORDSS(1).EQ.'SECAN') THEN
            NPOI2T=NPOI2T+1
            PROPSS(NPOI2T)=PARAMS(1) !svc
            WRITE(LURESS,856) PROPSS(NPOI2T)
            WRITE(LURESS,899)
         ELSE
            CALL RUNENDS('PROPM14: SECANT_VOLUMETRIC_C CARD NOT FOUND')
         ENDIF
         
      ELSE
         CALL RUNENDS('PROPM14: ONLY MM1 HAS BEEN IMPLEMENTED')
      ENDIF
C     
C--   OUTPUT (see pointes.f)
C     
c     se lee en la subrutina pointes, se lo hace igual
c     a 1 por que es una bandera que se usa para guardar
c     los datos en el pos.
      KPLA14S=1       
      
c     modelo          
      IPLLLS(INUPM)=14
      
      RETURN
 830  FORMAT(/,'CARBON=',E15.6)
 831  FORMAT(/,'SILICON=',E15.6)
 832  FORMAT(/,'MANGANESE=',E15.6)
 833  FORMAT(/,'CHROMIUM=',E15.6)
 834  FORMAT(/,'MICKEL=',E15.6)
 835  FORMAT(/,'CUPPER=',E15.6)
 836  FORMAT(/,'MOLYBDENUM=',E15.6)
c     
 840  FORMAT(/,'FVGRAPHITE=',E15.6)
 841  FORMAT(/,'FVFERRITE=',E15.6)
 842  FORMAT(/,'FVPEARLITE=',E15.6)
c     
 845  FORMAT(/,'SPACIN INTER PEARLITE=',E15.6,
     .     /,'GRAPHITE NODULE DENSITY=',E15.6)
c     
 847  FORMAT(/,'LATEN HEAT FERRITE=',E15.6,
     .     /,'LATEN HEAT PEARLITE=',E15.6)
c
 849  FORMAT(/,'SUBDIVITION STEP=',E15.6)
 850  FORMAT(/,'TEMP_DERIVATIVE_FPC=',E15.6)
c
 855  FORMAT(/,'MECHANICAL MODEL=',E15.6)
 856  FORMAT(/,'SECANT VOLUMETRIC CHANGE=',E15.6)
c     
 899  FORMAT(/)
      END
      
