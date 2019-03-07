
      
subroutine zbuildmat_vec(ndim,norder, wgeo, fkernel, q1, q2, &
     fgreens, par0, pars1, pars2, ntot, cmat)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: ndim, norder
  real *8, intent(in) :: wgeo(*)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  integer, intent(out) :: ntot
  complex *16, intent(out) :: cmat(*)
  external fkernel, fgreens

  !
  ! build **an arbitrary** layer potential matrix using bremer-style
  ! chunk quadratures. this assumed that each chunk is *locally*
  ! oriented counter-clockwise. 
  !
  ! input:
  !   norder - the order of the quadrature routine to use even
  !       though right now it enforces that norder is equal to
  !       the number of points on each chunk
  !   wgeo - the geometry structure created by chunks routines
  !   fkernel - the kernel evaluation routine, must have its
  !       calling sequence of the form:
  !           fkernel(par0, src, targ, src_normal, targ_normal,
  !               q1, q2, fgreens, pars1, pars2, val)
  !   fgreens- subroutine generating the greens function, its
  !       calling sequence must be of the form:
  !           fgreens(par0, src, targ, pars1, pars2, val, grad,
  !               hess)
  !     NOTE: be sure to take note that fgreens is **passed** to
  !         to the subroutine fkernel - that is, we make a
  !         distinction between the greens function and the kernel
  !
  ! output:
  !   ntot - the dimension of the matrix, ntot = k*nch*ndim
  !   cmat - the matrix which applies the kernel fkernel
  !
  !     This is mainly a memory management and unpacking routine
  !

  !
  ! first get some chunk info
  !
  call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders, &
       iders2, ihs)

  ntot = k*nch*ndim
  if (norder .ne. k) then
     call prinf('norder = *', norder, 1)
     call prinf('k = *', k, 1)
     call prinf('bomb!!!!*', k, 0)
     stop
  endif

  call zbuildmat_vec0(ndim,norder, k, nch, wgeo, fkernel, &
       q1, q2, &
       fgreens, par0, pars1, pars2, ntot, cmat)

  return
end subroutine zbuildmat_vec





subroutine zbuildmat_vec0(ndim,norder, k, nch, wgeo, &
     fkernel, q1, q2, &
     fgreens, par0, pars1, pars2, ntot, zmat)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: ndim,norder, k, nch, ntot
  real *8, intent(in) :: wgeo(*)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  complex *16, intent(out) :: zmat(ntot,ntot)

  external fkernel, fgreens

  real *8, allocatable :: xs(:), whts(:)
  real *8, allocatable :: xs1(:), whts1(:)
  real *8, allocatable :: xs0(:,:), whts0(:,:)
  complex *16, allocatable :: tempmat(:,:)
  complex *16, allocatable :: tempnear(:,:,:)
  complex *16, allocatable :: tempself(:,:,:)
  integer ntmat

  !
  ! actually go ahead and build the matrix, block by block
  !
  done = 1
  ntmat = ndim*k
  allocate(tempmat(ntmat,ntmat))

  itype = 1
  allocate(xs(k), whts(k))
  call legerts(itype, k, xs, whts)

  !
  ! get the bremer near and self quadratures
  !

  do j = 1,ntot
     do i = 1,ntot
        zmat(i,j) = 0
     enddo
  enddo

  call getquadsinfo(k, nquad1, nquad0)
  allocate(xs1(nquad1), whts1(nquad1))
  allocate(xs0(k,k), whts0(k,k))
  call getquads(k, nquad1, xs1, whts1, nquad0, xs0, whts0)

  allocate(tempnear(nquad1,ndim,ndim))
  allocate(tempself(k,ndim,ndim))

  !$omp parallel do default(shared) &
  !$omp     private(itarget, imat, ibefore, iafter) &
  !$omp     private(jsource,jmat,tempmat,tempnear,tempself)

  do itarget = 1,nch
     imat = 1 + (itarget-1)*ntmat

     call chunk_neighbors(itarget, wgeo, ibefore, iafter)

     !
     ! for all target chunks, loop over source chunks
     !

     do jsource = 1,nch
        jmat = 1 + (jsource-1)*ntmat

        ifoff = 1

        if (jsource .eq. ibefore) then
           call znear_buildmat_vec(ndim,itarget, jsource, &
                k, wgeo, &
                fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                nquad1, xs1, whts1, tempmat,tempnear)
           call zinsertmat(ntmat,ntmat, tempmat, imat, jmat, &
                ntot, ntot, zmat)
           ifoff = -1
        endif

        if (jsource .eq. iafter) then
           call znear_buildmat_vec(ndim,itarget, jsource, k, wgeo, &
                fkernel, q1, q2, fgreens, par0, pars1, pars2, &
                nquad1, xs1, whts1, tempmat,tempnear)
           call zinsertmat(ntmat,ntmat, tempmat, imat, jmat, &
                ntot, ntot, zmat)
           ifoff = -1
        endif

        if (jsource .eq. itarget) then
           call zdiagonal_buildmat_vec(ndim,itarget, k, &
                wgeo, fkernel, &
                q1, q2, fgreens, par0, pars1, pars2, xs0, &
                whts0, tempmat, tempself)
           call zinsertmat(ntmat,ntmat, tempmat, imat, jmat, &
                ntot, ntot, zmat)
           ifoff = -1
        endif

        if (ifoff .eq. 1) then
           call zfar_buildmat_vec(ndim,itarget, jsource, &
                k, wgeo, &
                fkernel, q1, q2, fgreens, par0, pars1, pars2, xs, &
                whts, tempmat)
           call zinsertmat(ntmat,ntmat, tempmat, imat, &
                jmat, ntot, &
                ntot, zmat)
        endif

     enddo
  enddo

  !$omp end parallel do

end subroutine zbuildmat_vec0





subroutine zfar_buildmat_vec(ndim,ich_targ, ich_src, k, wgeo, &
     fkernel, q1, q2, fgreens, par0, pars1, pars2, xs, &
     whts, amat)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim
  real *8, intent(in) :: wgeo(*), xs(k), whts(k)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  external fkernel, fgreens

  real *8 :: chunk_src(2000), chunk_targ(2000)
  real *8 :: der_src(2000), der_targ(2000)

  !
  ! build a kth order smooth quadrature on for the far
  ! field interation between targ and src chunks
  !

  call chunkaddr(ich_src, wgeo, k_src, ichunk_src, &
       iadj_src, ider_src, ider2_src, h_src)

  call chunkaddr(ich_targ, wgeo, k_targ, ichunk_targ, &
       iadj_targ, ider_targ, ider2_targ, h_targ)

  len = 2*k
  call chunks_dcopy(len, wgeo(ichunk_src), chunk_src)
  call chunks_dcopy(len, wgeo(ider_src), der_src)

  call chunks_dcopy(len, wgeo(ichunk_targ), chunk_targ)
  call chunks_dcopy(len, wgeo(ider_targ), der_targ)

  call zfar_buildmat_vec0(ndim,k, chunk_targ, der_targ, &
       h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
       fgreens, par0, pars1, pars2, xs, whts, amat)

  return
end subroutine zfar_buildmat_vec


subroutine zfar_buildmat_vec0(ndim,k, chunk_targ, &
     der_targ, h_targ, chunk_src, der_src, h_src, fkernel, &
     q1, q2, fgreens, par0, pars1, pars2, xs, whts, amat)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim
  real *8, intent(in) :: chunk_targ(2,k), der_targ(2,k), h_targ
  real *8, intent(in) :: chunk_src(2,k), der_src(2,k), h_src
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  real *8, intent(in) :: xs(k), whts(k)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  external fgreens, fkernel

  real *8 :: src(10), targ(10)
  real *8 :: target_normal(10), source_normal(10)
  complex *16 :: val(ndim,ndim)

  !
  ! fill up the matrix with the smooth quadrature
  !

  done = 1

  do i = 1,k

     targ(1) = chunk_targ(1,i)
     targ(2) = chunk_targ(2,i)

     xp_targ = der_targ(1,i)
     yp_targ = der_targ(2,i)
     ds_targ = sqrt(xp_targ**2 + yp_targ**2)

     target_normal(1) = yp_targ/ds_targ
     target_normal(2) = -xp_targ/ds_targ

     do j = 1,k

        src(1) = chunk_src(1,j)
        src(2) = chunk_src(2,j)

        xp = der_src(1,j)
        yp = der_src(2,j)

        dd = sqrt(xp**2 + yp**2)
        source_normal(1) = yp/dd
        source_normal(2) = -xp/dd

        call fkernel(par0, src, targ, source_normal, &
             target_normal, q1, q2, fgreens, pars1, pars2, val)

        rnorm = sqrt(xp**2+yp**2)
        dsdt = rnorm*h_src
        do l1 = 1,ndim
           do l2 = 1,ndim
              amat((i-1)*ndim+l1,(j-1)*ndim+l2) = &
                   val(l1,l2)*whts(j)*dsdt
           enddo
        enddo
     enddo
  enddo

  return
end subroutine zfar_buildmat_vec0





subroutine znear_buildmat_vec(ndim,ich_targ, ich_src, k, wgeo, &
     fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
     xs, whts, amat, tempnear)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim, nquad
  real *8, intent(in) :: wgeo(*), xs(nquad), whts(nquad)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  complex *16, intent(inout) :: tempnear(nquad,ndim,ndim)
  external fkernel, fgreens

  !
  ! build the matrix performing a neighbor interaction using
  ! bremer style quadratures
  !

  real *8 :: chunk_src(1000)
  real *8 :: chunk_targ(1000), der_src(1000), der_targ(1000)

  call chunkaddr(ich_src, wgeo, k_src, ichunk_src, &
       iadj_src, ider_src, ider2_src, h_src)

  call chunkaddr(ich_targ, wgeo, k_targ, ichunk_targ, &
       iadj_targ, ider_targ, ider2_targ, h_targ)

  len = 2*k
  !!call chunks_dcopy(len, wgeo(ichunk_src), chunk_src)
  !!call chunks_dcopy(len, wgeo(ider_src), der_src)

  !!call chunks_dcopy(len, wgeo(ichunk_targ), chunk_targ)
  !!call chunks_dcopy(len, wgeo(ider_targ), der_targ)

  call znear_buildmat_vec0(ndim,k, wgeo(ichunk_targ), wgeo(ider_targ), &
       h_targ, wgeo(ichunk_src), wgeo(ider_src), h_src, &
       fkernel, q1, q2, fgreens, par0, pars1, pars2, nquad, &
       xs, whts, amat, tempnear)

  return
end subroutine znear_buildmat_vec


subroutine znear_buildmat_vec0(ndim,k, chunk_targ, der_targ, &
     h_targ, chunk_src, der_src, h_src, fkernel, q1, q2, &
     fgreens, par0, pars1, pars2, nquad, xs, whts, amat,temp)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim, nquad
  real *8, intent(in) :: chunk_targ(2,k), der_targ(2,k), h_targ
  real *8, intent(in) :: chunk_src(2,k), der_src(2,k), h_src
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  real *8, intent(in) :: xs(nquad), whts(nquad)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  external fgreens, fkernel
  complex *16, intent(inout) :: temp(nquad,ndim,ndim)

  real *8 :: src(10), targ(10), source_normal(10)
  real *8 :: ts(1000), ainterp(1000000), target_normal(10)
  real *8 :: xcoefs_src(1000), ycoefs_src(1000)
  real *8 :: xpcoefs_src(1000), ypcoefs_src(1000)
  complex *16 :: val(ndim,ndim), w(1000000)

  !
  ! build the k by k near interaction matrix.
  !

  done = 1
  call chunksexps(k, chunk_src, xcoefs_src, ycoefs_src)
  call chunksexps(k, der_src, xpcoefs_src, ypcoefs_src)

  call lematrin(k, nquad, xs, ainterp, ts, w)

  do i = 1,k

     targ(1) = chunk_targ(1,i)
     targ(2) = chunk_targ(2,i)

     xp_targ = der_targ(1,i)
     yp_targ = der_targ(2,i)
     ds_targ = sqrt(xp_targ**2 + yp_targ**2)

     target_normal(1) = yp_targ/ds_targ
     target_normal(2) = -xp_targ/ds_targ


     do j = 1,nquad

        t_src = xs(j)
        call legeexev(t_src, src(1), xcoefs_src, k-1)
        call legeexev(t_src, src(2), ycoefs_src, k-1)
        call legeexev(t_src, xp, xpcoefs_src, k-1)
        call legeexev(t_src, yp, ypcoefs_src, k-1)

        dd = sqrt(xp**2 + yp**2)
        source_normal(1) = yp/dd
        source_normal(2) = -xp/dd

        call fkernel(par0, src, targ, source_normal, &
             target_normal, q1, q2, fgreens, pars1, pars2, val)

        dsdt = sqrt(xp**2+yp**2)*h_src
        do l1 = 1,ndim
           do l2 = 1,ndim
              temp(j,l1,l2) = val(l1,l2)*whts(j)*dsdt
           enddo
        enddo
     enddo

!!!!call lematrin(k, nquad, xs, ainterp, ts, w)
     i1 = 1
     do l1 = 1,ndim
        do l2 = 1,ndim
           call zmat_by_dmat(i1, nquad, temp(1,l1,l2), k, &
                ainterp, w)
           do j = 1,k
              amat((i-1)*ndim+l1,(j-1)*ndim+l2) = w(j)
           enddo
        enddo
     enddo
  enddo
  return
end subroutine znear_buildmat_vec0



subroutine zdiagonal_buildmat_vec(ndim,ich, k, wgeo, fkernel, q1, q2, &
     fgreens, par0, pars1, pars2, xs0, whts0, amat,tempself)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim
  real *8, intent(in) :: wgeo(*), xs0(k,k), whts0(k,k)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  complex *16, intent(inout) :: tempself(k,ndim,ndim)
  external fkernel, fgreens

  !
  ! build the k by k self interaction matrix. this routine (will) builds
  ! matrices which are square-root weights by default, since it's
  ! the morally correct thing to do. this uses Bremer quadratures.
  !
  ! input:
  !   ich - the chunk number
  !   k - dimension of amat, i.e. the number of points on the chunk
  !   wgeo - geometry structure created by chunks.f
  !   fkernel - the arbitrary kernel routine, see above
  !       documention since this routine receives the greens function
  !   fgreens - the single layer greens function routine
  !       it must have the calling sequence
  !           fgreens(par0, src, targ, pars1, pars2, val)
  !       and val is assumed to be real
  !   par0, pars1, pars2 - inputs to fgreens
  !
  ! output:
  !   amat - the self-interaction matrix
  !

  call chunkaddr(ich, wgeo, k, ichunk, iadj, ider, ider2, h)

  call zdiagonal_buildmat_vec0(ndim,k, wgeo(ichunk), &
       wgeo(ider), h, fkernel, q1, q2, fgreens, par0, &
       pars1, pars2, xs0, whts0, amat,tempself)

  return
end subroutine zdiagonal_buildmat_vec


subroutine zdiagonal_buildmat_vec0(ndim,k, xys, dxys, h, &
     fkernel, q1, q2, fgreens, par0, pars1, pars2, &
     xs, whts, amat,temp)
  implicit real *8 (a-h,o-z)
  integer, intent(in) :: k,ndim 
  real *8, intent(in) :: xys(2,k), dxys(2,k), h
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(in) :: q1(*), q2(*)
  real *8, intent(in) :: xs(k,k), whts(k,k)
  complex *16, intent(out) :: amat(ndim*k,ndim*k)
  complex *16, intent(inout) :: temp(k,ndim,ndim)
  external fgreens, fkernel

  real *8 :: src(10), targ(10), source_normal(10)
  real *8 :: ts(1000), ainterp(1000000), target_normal(10)
  real *8 :: xcoefs(1000)
  real *8 :: ycoefs(1000), xpcoefs(1000), ypcoefs(1000)
  complex *16 :: val(ndim,ndim), w(1000000)

  !
  ! build the k by k self interaction matrix. no square-root
  ! weighting yet... that is handled in the outer routine,
  ! xxx_xxx_direct.
  ! Note: k should not be larger than 800
  !

  done = 1
  call chunksexps(k, xys, xcoefs, ycoefs)
  call chunksexps(k, dxys, xpcoefs, ypcoefs)

  do i = 1,k

     targ(1) = xys(1,i)
     targ(2) = xys(2,i)

     xp_targ = dxys(1,i)
     yp_targ = dxys(2,i)
     ds_targ = sqrt(xp_targ**2 + yp_targ**2)

     target_normal(1) = yp_targ/ds_targ
     target_normal(2) = -xp_targ/ds_targ

     do j = 1,k

        t_src = xs(j,i)
        call legeexev(t_src, src(1), xcoefs, k-1)
        call legeexev(t_src, src(2), ycoefs, k-1)
        call legeexev(t_src, xp, xpcoefs, k-1)
        call legeexev(t_src, yp, ypcoefs, k-1)

        dd = sqrt(xp**2 + yp**2)
        source_normal(1) = yp/dd
        source_normal(2) = -xp/dd

        call fkernel(par0, src, targ, source_normal, &
             target_normal, q1, q2, fgreens, pars1, pars2, val)

        dsdt = sqrt(xp**2+yp**2)*h
        do l1=1,ndim
           do l2=1,ndim
              temp(j,l1,l2) = val(l1,l2)*whts(j,i)*dsdt
           enddo
        enddo
     enddo

     call lematrin(k, k, xs(1,i), ainterp, ts, w)
     i1 = 1

     do l1 = 1,ndim
        do l2 = 1,ndim
           call zmat_by_dmat(i1, k, temp(1,l1,l2), k, &
                ainterp, w)
           do j = 1,k
              amat((i-1)*ndim+l1,(j-1)*ndim+l2) = w(j)
           enddo
        enddo
     enddo

  enddo

  return
end subroutine zdiagonal_buildmat_vec0

