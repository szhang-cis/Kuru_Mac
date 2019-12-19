      subroutine micrSegrL(chemVal, iMicrSegr, partitionCoef, 
     .     initChemComp, fsolid, tf, Ds, PDAs, SDAs)
C     no admitimos variables declaradas implicitamente
      implicit none
C     
C***  subroutine: calcula la segreagacion de los elementos de aleacion a 
C     escala micro
C     
C     modelos de segregacion implementados -- 03/11/2010 @ buenaluna
C     
C     imicrsegr= 1 -- lever rule --
C     imicrsegr= 2 -- scheil rule (1942) --
C     imicrsegr= 3 -- brody-flemings rule (1996) --
C     imicrsegr= 4 -- clyne - kurz  -- 1981
C     imicrsegr= 5 -- ohnaka rule (1986) --
C     imicrsegr= 6 -- without segregation --
C     
C***  NOTA IMPORTANTE: a todas les he sacado del segundo miembro el 
C     coeficiente de particion y lo he pasado al 
C     primer miembro para obtener la CL (liquido) y
C     no la CS (solido).
C     
C***  author    : Fernando Diego Carazo Rodriguez (buenaluna)
C     
C***  date      : mié nov  3 12:33:47 ART 2010
C     Today is Boomtime, the 15th day of The Aftermath 
C     in the YOLD 3176
      
C     declaracion de las variables
      integer*4 iMicrSegr
      real*8  partitionCoef, initChemComp, alpha, omega, 
     .     beta, gamma, fsolid, Ds, tf, SDAs, PDAs, 
     .     chemVal
C     
      chemval= 0.0D0
C     
C     seleccionamos el modelo de particion
      if(iMicrSegr.eq.1) then
C***  lever rule
C     
         chemVal=initChemComp/((1.0D0-fsolid)+
     .        partitionCoef*fsolid)
C     
      elseif(iMicrSegr.eq.2) then
C***  scheil equation -- 1942
C     
         chemVal= initChemComp*(1.0D0-fsolid)**(partitionCoef-1.0D0)
C     
      elseif(iMicrSegr.eq.3) then
C***  brody-flemings -- 1996
C     variable auxiliar
         alpha= (4.0D0*Ds*tf) / ((SDAs/2.0D0)**2.0D0)
C     
         chemVal= initChemComp*(1.0D0-
     .        (1.0D0-2.0D0*alpha*partitionCoef)*fsolid)
     .        **((partitionCoef-1.0D0)/(1.0D0-2.0D0*alpha*
     .        partitionCoef))
C     
      elseif(iMicrSegr.eq.4) then
C***  clyne - kurz  -- 1981
C     variable auxiliar
         alpha= 4.0D0*Ds*tf/PDAs**2.0D0
         omega= alpha*(1- dexp(-1.0D0/alpha))-5.0D-1*
     .        dexp(-1.0D0/(2.0D0*alpha))
C     
         chemVal= initChemComp*(1.0D0-
     .        (1.0D0-2.0D0*omega*partitionCoef)*fsolid)
     .        **((partitionCoef-1.0D0)/(1.0D0-2.0D0*omega*
     .        partitionCoef))
C     
      elseif(iMicrSegr.eq.5) then
C***  ohnaka  -- 1986
C     variable auxiliar
         gamma= 8.0D0*Ds*tf/PDAs**2.0D0
         beta= 2.0D0*gamma/(1.0D0+2.0D0*gamma)
C     
         chemVal= initChemComp*(1.0D0-
     .        (1.0D0-2.0D0*beta*partitionCoef)*fsolid)
     .        **((partitionCoef-1.0D0)/(1.0D0-2.0D0*beta*
     .        partitionCoef))
      else
C***  no hay segregacion
         chemVal= initChemComp
      endif
C     fin micrSegrL
      return
      end
