c------------------------------------------------------------------------
      subroutine trtv(ntype,nstrs,st,sv,c1,c2,ifl)
c*****************************************************************
c
c****transforma un tensor de tension
c                       vector  ---> tensor    PARA: ifl=1
c                       tensor  ---> vector    PARA: ifl=2
c
c*****************************************************************
      implicit real*8(a-h,o-z)
      dimension st(3,3), sv(nstrs),c1(6),c2(6)
C
      if(ifl.eq.1)then
c
c...construye el tensor
c
      do 20 i=1,3
      do 20 j=1,3
       st(i,j)=0.0
  20  continue
c
      do 50 istrs=1,nstrs
       sv(istrs)=sv(istrs)*c1(istrs)
  50  continue
c
       st(1,1)=sv(1)
       if(ntype.eq.4)then
         st(1,2)=sv(3)
         st(1,3)=sv(5)
         st(2,2)=sv(2)
         st(2,3)=sv(6)
         st(3,3)=sv(4)
       else if(ntype.ne.5) then
         st(1,2)=sv(3)
         st(2,2)=sv(2)
         if(ntype.ne.1) st(3,3)=sv(4)
       endif
       do 60  i=1,3
       do 60  j=i,3
   60  st(j,i)=st(i,j)
      else
c
c...construye el vector
c
       sv(1)=st(1,1)
       if(ntype.eq.4)then
         sv(3)=st(1,2)
         sv(5)=st(1,3)
         sv(2)=st(2,2)
         sv(6)=st(2,3)
         sv(4)=st(3,3)
       else if(ntype.ne.5) then
         sv(3)=st(1,2)
         sv(2)=st(2,2)
         if(ntype.ne.1) sv(4)=st(3,3)
       endif
c
      do 61 istrs=1,nstrs
       sv(istrs)=sv(istrs)*c2(istrs)
  61  continue
c
      endif
c
      return
      end
