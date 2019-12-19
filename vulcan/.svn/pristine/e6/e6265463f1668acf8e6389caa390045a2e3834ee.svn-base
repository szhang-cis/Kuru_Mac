c=======================================================================
      subroutine scalpr(x,y,idim,prod)
c-----------------------------------------------------------------------
c   computes the scalar product of two real vectors
c         
c    datas :
c           x    : the first vector
c           y    : the second vector
c           idim : dimension of the two vectors
c          
c    result :
c           prod : the scalar product of x and y
c       
      real*8 x,y,prod
      integer*4 idim,i
c         
      dimension x(*),y(*)
c
      prod=0.d0       
      do 10 i=1,idim
         prod=prod+x(i)*y(i)     
 10   continue
c
      return
      end
c=======================================================================
      subroutine sctvec(a,x,z,idim)
c-----------------------------------------------------------------------
c   computes a scalar times a vector : z=a*x
c
c   datas :
c          a    : the scalar
c          x    : the vector
c          idim : dimension of x 
c
c   result :
c          z : a times x
c
      implicit none
      real*8 a,x,z
      integer*4 idim,i
c
      dimension x(*),z(*)
c 
      do 10 i=1,idim
         z(i)=a*x(i)
 10   continue
c
      return
      end
c=======================================================================
      subroutine advect(a,x,b,y,z,idim)
c-----------------------------------------------------------------------
c this subroutine computes the linear combinaison of two vectors :
c      z=a*x+b*y
c
c   datas :
c          a    : coefficient of the vector x
c          x    : the first vector considered
c          b    : coefficient of the second vector
c          y    : the second vector considered
c          idim : the dimension of all the vectors
c   
c   results :
c          z : the linear combinaison of the vectors x and y
c  
      implicit none

      real*8 a,x,b,y,z
      integer*4 idim,i
c
      dimension x(*),y(*),z(*)     
c
      do 10 i=1,idim
         z(i)=a*x(i)+b*y(i)
 10   continue
c
      return
      end
c
c=======================================================================
      subroutine solve(h,y,e,idim)
c-----------------------------------------------------------------------
c     solves an upper linear system : hy=e
c
c     datas :
c             h : the triangular martrix stored column by column
c             e : the second member of the linear system
c             idim : the dimension of e
c
c     result :
c            y : hy=e
      implicit none
c    
      real*8 h, y, e
      integer*4 idim, i, j, ly, lhij, lhii
c
      dimension h(*), y(*), e(*)
c
      do 20 i=1,idim
         ly=idim-i+1
         y(ly)=e(ly)
         do 10 j=ly+1,idim
            lhij=((j-1)*j)/2+ly
            y(ly)=y(ly)-h(lhij)*y(j)
 10      continue
         lhii=((ly-1)*ly)/2+ly
         y(ly)=y(ly)/h(lhii)
 20   continue
c                
      return
      end         
c
c=======================================================================
