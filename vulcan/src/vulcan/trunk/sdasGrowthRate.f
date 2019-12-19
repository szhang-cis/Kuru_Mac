      subroutine sdasGrowthRate(isdasgr,difcl,cogth,cla,ca,amlq,akca,tf,
     *     sdas,sdasini)
C     no admitimos variables declaradas implicitamente
      implicit none
C     
C***  subroutine: calcula el espaciamiento entre las ramas secundarias 
C                 de las dendritas de austenita al final de la
C                 solidificacion
C                 
C     modelos de segregacion implementados -- 03/11/2010 @ buenaluna
C     
C     isdasgr= 1 -- kattamis-flemings --
C     isdasgr= 2 --  --
C     isdasgr= 3 --  --
C     isdasgr= 4 --  --
C     isdasgr= 5 --  --
C     isdasgr= 6 --  --
C     
C***  NOTA IMPORTANTE: extraidas del libro de Stefanescu (page 181).
C     
C***  author: Fernando Diego Carazo Rodriguez (buenaluna)
C     
C***  date  : mié dic  8 17:26:04 ART 2010
C             Today is Boomtime, the 50th day of The Aftermath in the 
C             YOLD 3176 Celebrate Afflux
      
C     declaracion de las variables
      integer*4 isdasgr
      real*8  difcl, cogth, cla, ca, amlq, akca, tf, sdas, sdasini
C     seleccionamos el modelo de particion
      if(isdasgr.eq.1) then
C***  kattamis - flemings (1965)
C     
         if(cla.le.ca) then
            sdas= sdasini
         else
            sdas= 10.0*((difcl*-1.0D0*cogth* dlog(cla/ca))/
     .           (amlq*(1.0D0-akca)*(cla-ca))*tf)**(1.0D0/3.0D0) 
         endif
C***  
C     
      elseif(isdasgr.eq.2) then
C***  ardell (1972)
C     
      elseif(isdasgr.eq.3) then
C***  voorhees - glicksman (1984) 
C     
      elseif(isdasgr.eq.4) then
C***  mortensen (1991)
C     
      else
C***  no exsite crecimiento del espaciado de las ramas secundarias
         sdas= 0.0D0
      endif
C***  
C     fin sdasFrowthRate
C     
      return
      end subroutine
