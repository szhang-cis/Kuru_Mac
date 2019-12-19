c------------------------------------------------------------------------
      subroutine tenvec(ndime,nstrs,vect,rmat1,index,ifl)
c*********************************************************************
c
c**** esta rutina rota el tensor de tensiones y deformaciones
c
c             index  =1  tens TENSIONES
c             index  =2  tens DEFORMACIONES
c
c             ifl    =1  Glob >>>> Loc
c             ifl    =2  Loc  >>>> Glo
c
c*********************************************************************
      implicit real*8(a-h,o-z)
c
      dimension vect(*),tens(3,3),tran(3,3),tena(3,3),rmat1(ndime,*),
     .          c1(6),c2(6)
c
c*** cambio de las componentes de corte en el tensor de deformaciones
c
       do 5 istrs=1,nstrs
        c1(istrs)=1.0
        c2(istrs)=1.0
   5   continue
       if(index.eq.2)then
        c1(3)=0.5
        c1(5)=0.5
        c1(6)=0.5
        c2(3)=2.0
        c2(5)=2.0
        c2(6)=2.0
       endif
c
c*** construccion de la matriz de cambio de base
c
       do 10 idime=1,3
       do 10 jdime=1,3
  10   tran(idime,jdime)=0.0
c
       do 11 idime=1,ndime
       do 11 jdime=1,ndime
       if(ifl.eq.1)tran(idime,jdime)=rmat1(idime,jdime)
       if(ifl.eq.2)tran(jdime,idime)=rmat1(idime,jdime)
  11   continue
       if(ndime.eq.2)tran(3,3)=1.0
c
c*** transformacion de vector a tensor
c
       call trtv(ntype,nstrs,tens,vect,c1,c2,     1)    
c
c*** producto tensorial simplemente contraido
c
       do 38 i=1,3
       do 38 j=1,3
       tena(i,j)=0.0
       do 38 k=1,3
       do 38 l=1,3
       tena(i,j)=tena(i,j)+tran(k,i)*tens(k,l)*tran(l,j)
  38   continue
c
       call trtv(ntype,nstrs,tena,vect,c1,c2,     2)    
c
      return
      end
