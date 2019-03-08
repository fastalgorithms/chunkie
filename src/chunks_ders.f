c     Copyright (C) 2019: Travis Askham
c     Contact: askham@njit.edu
c      
cc This software is being released under a modified FreeBSD license
cc (see LICENSE in home directory). 

      subroutine chunkderf(f,df,ndim,k,nch,chunks,ders,hs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the derivative with respect to arc-length of the real 
c     function f (which can be vector valued)
c
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ndim, k, nch
      real *8 f(ndim,k,*), df(ndim,k,*), chunks(2,k,*), ders(2,k,*)
      real *8 hs(*)
c     local variables
      real *8, allocatable :: u(:), v(:), x(:), whts(:)
      real *8, allocatable :: pol(:), dpol(:), ftemp(:,:)
      integer i, j, l
      integer korder, itype
      real *8 dsdt

      korder = k-1

      allocate(u(k*k),v(k*k),x(k),whts(k),pol(k),dpol(k),ftemp(k,ndim))

      itype = 2
      call legeexps(itype,k,x,u,v,whts)

      do i = 1,nch
         do j = 1,k
         do l = 1,ndim
            ftemp(j,l) = f(l,j,i)
         enddo
         enddo
         do l = 1,ndim
c     get coeffs, compute coeffs of derivative, get values of derivative
            call amatchunkders(u,ftemp(1,l),pol,k)
            call legediff(pol,korder,dpol)
            call amatchunkders(v,dpol,ftemp(1,l),k)
         enddo
         do j = 1,k
            dsdt = dsqrt(ders(1,j,i)**2+ders(2,j,i)**2)*hs(i)
            do l = 1,ndim
               df(l,j,i) = ftemp(j,l)/dsdt
            enddo
         enddo
      enddo

      return
      end

      subroutine chunkderf_trans(f,df,ndim,k,nch,chunks,ders,hs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the transpose of the differentiation matrix
c     with respect to arc-length applied to the function f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ndim, k, nch
      real *8 f(ndim,k,*), df(ndim,k,*), chunks(2,k,*), ders(2,k,*)
      real *8 hs(*)
c     local variables
      real *8, allocatable :: u(:), v(:), x(:), whts(:)
      real *8, allocatable :: pol(:), dpol(:), ftemp(:,:)
      integer i, j, l
      integer korder, itype
      real *8 dsdt

      korder = k-1

      allocate(u(k*k),v(k*k),x(k),whts(k),pol(k),dpol(k),ftemp(k,ndim))

      itype = 2
      call legeexps(itype,k,x,u,v,whts)

      do i = 1,nch
         do j = 1,k
         do l = 1,ndim
            ftemp(j,l) = f(l,j,i)
         enddo
         enddo
         do l = 1,ndim
c     get coeffs, compute coeffs of derivative, get values of derivative
            call amatchunkders(u,ftemp(1,l),pol,k)
            call legediff(pol,korder,dpol)
            call amatchunkders(v,dpol,ftemp(1,l),k)
         enddo
         do j = 1,k
            dsdt = dsqrt(ders(1,j,i)**2+ders(2,j,i)**2)*hs(i)
            do l = 1,ndim
               df(l,j,i) = ftemp(j,l)/dsdt
            enddo
         enddo
      enddo

      return
      end

      subroutine amatchunkders(a,x,y,n)
c     evaluates y = a*x
      implicit none
      real *8 a(n,n), x(*), y(*)
      integer n, i, j
      
      do i = 1,n
         y(i) = 0.0d0
      enddo

      do j = 1,n
         do i = 1,n
            y(i) = y(i)+a(i,j)*x(j)
         enddo
      enddo

      return
      end

