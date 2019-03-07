c       Copyright (C) 2014: Michael O'Neil
c       Contact: oneil@cims.nyu.edu
c 
c       This program is free software; you can redistribute it and/or
c       modify it under the terms of the GNU General Public License as
c       published by the Free Software Foundation; either version 2 of
c       the License, or (at your option) any later version.  This
c       program is distributed in the hope that it will be useful, but
c       WITHOUT ANY WARRANTY; without even the implied warranty of
c       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c       GNU General Public License for more details. You should have
c       received a copy of the GNU General Public License along with
c       this program; if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       $Date$
c       $Revision$
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the code below contains several user-callable routines
c       for the construction of refind polygons chunked up smooth
c       curves. The most useful routines are:
c
c         chunk2trap - takes a chunked discretization and outputs a
c             global discretization which is sampled in arc length - 
c             note that this is most often NOT an optimal discretization
c     
c         chunks_dist - returns the distance from a point to a chunk, along
c             with the point on the chunk which is the closest (most likely
c             not one of the discretization points)
c
c         chunksplot - plots either the chunks or the chunks and normals
c
c         chunkfunc_pack - same as chunkfunc, except returns
c             a packaged work array
c
c         chunkfunc - passing in a subroutine which describes a 
c             curve analytically, returns an adpative chunked
c             version given by Legendre nodes
c
c         chunknormals - calculate the normal at each point
c             on the chunks
c
c         chunkdsdt - calculates the value of dsdt at each point on
c           the curve so that the usual gaussian weights integrate
c           functions correctly
c
c         chunkrls - calculates the cumulative archlength at each
c           of the discretization points, starting with chunk 1
c
c         chunkcurvatures - calculates the curvature
c           at each point on the chunk
c
c         chunkpoly - constructs a chunked version of a polygon or
c             an open curve
c
c         chunksplit_pack - splits a chunk that is already stored
c             inside a work array created by chunkpack
c
c         chunksplit1 - using only chunk points (not an underlying
c             parameterization), split the chunk into two pieces
c
c         chunksplit2 - using only chunk points (not an underlying
c             parameterization), split the chunk into two pieces
c             and resample densities as well
c
c         chunksort - arranges the chunks in memory so that they are
c             oriented end-to-end counterclockwise instead of
c             scattered based on the adaptive construction.
c
c         chunkpack - packs up the chunk into into one array
c
c         chunkunpack - unpacks the chunk from the array created
c             by chunkpack
c
c         chunkunpack1 - unpacks the chunk from the array created
c             by chunkpack, returns only the array indices
c
c         chunksize - takes one discretized chunk and computes its
c             arclength and midpoint
c
c         chunklength - receives the same subroutine that chunkfunc
c             does, calculates the arclength between
c             two specified points (in parameterization space)
c
c         chunkwhts - returns smooth quadrature weights for every
c             point on the curve
c
c         chunkinterp2 - interpolates 2d function from k nodes
c             on a chunk to kout nodes
c
c         chunkinterp2f - interpolates 2d function from k nodes
c                   on a chunk to kout nodes fast. Requires
c                   function to legendre coefficient for k nodes
c                   and kout legendre coefficent to function
c                   matrices as input
c
c         chunkres - computes accuracy of a real-valued interpolated 
c             functio on the chunk nodes by checking the tail of
c             the expansions
c
c         chunkzres - computes accuracy of a complex-valued interpolated 
c             function on the chunk nodes by checking the tail of
c             the expansions
c
c         NOTE: all these routines assume that the number of 
c             nodes that you put on the chunks is something like
c             k=16 instead of k=1000. things will break if you insist
c             on k being huge, namely probably over 100.
c
c       the follow routines are utility routines, and may or may not
c       be useful to the user:
c
c         chunk2to1 - this routine separates a 2 by n array into 2
c             separate length n arrays, i.e. grabs real and imaginary
c             parts
c
c         chunkreverse - reverses the orientation of the chunks
c            from counterclockwise to clockwise. you won't need this
c            routine unless you really know what you're doing.
c
c         chunkmatvec - matrix vector multiplication
c
c         chunkcopy - copies one array to another
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c        
        subroutine chunk2trap(wgeo, nout, rltot, xys, ders, ders2)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1), xys(2,1), ders(2,1), ders2(2,1)
c
c       using the geometry in wgeo, output a global discretization
c       in terms of arclength sampled at nout points
c
        done = 1
        pi = 4*atan(done)
        


c
c       if the curve is not closed, add in the last point
c


c
        return
        end
c
c
c
c
c
        subroutine chunk_reorder(wgeo, iorder, wgeosort)
        implicit real *8 (a-h,o-z)
        integer iorder(*)
        real *8 wgeo(*), wgeosort(*)
c
c       this routine reorders the chunks in wgeo so that they
c       are ordered *in memory* according to iorder
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     $      iders2, ihs)

        do i = 1,ichunks-1
          wgeosort(i) = wgeo(i)
        enddo
c
        call chunk_reorder0(k, nch, iorder, wgeo(ichunks), wgeo(iadjs),
     $      wgeo(iders), wgeo(iders2), wgeo(ihs), wgeosort(ichunks),
     $      wgeosort(iadjs), wgeosort(iders), wgeosort(iders2),
     $      wgeosort(ihs))
c
        return
        end
c
c
c
        subroutine chunk_reorder0(k, nch, iorder, chunks, adjs, ders, 
     $      ders2, hs, chunksout, adjsout, dersout, ders2out, hsout)
        implicit real *8 (a-h,o-z)
        integer :: iorder(nch), adjs(2,nch), adjsout(2,nch)
        real *8 :: chunks(2,k,nch), ders(2,k,nch), ders2(2,k,nch)
        real *8 :: hs(nch)
        real *8 :: chunksout(2,k,nch), dersout(2,k,nch)
        real *8 :: ders2out(2,k,nch), hsout(nch)
c
c       reorder the chunks, and the adjacency information, according
c       to the ordering scheme iorder
c
        do i = 1,nch
          ind = iorder(i)
          do j = 1,k
            chunksout(1,j,i) = chunks(1,j,ind)
            chunksout(2,j,i) = chunks(2,j,ind)
            dersout(1,j,i) = ders(1,j,ind)
            dersout(2,j,i) = ders(2,j,ind)
            ders2out(1,j,i) = ders2(1,j,ind)
            ders2out(2,j,i) = ders2(2,j,ind)
            hsout(i) = hs(ind)
          enddo
          ibefore = adjs(1,i)
          iafter = adjs(2,i)
          adjsout(1,i) = iorder(ibefore)
          adjsout(2,i) = iorder(iafter)
        enddo
c        
        return
        end
c
c
c
c
c
cccc        call chunk_reorder_zvals(k, nch, wlists(iisource), charge,
cccc     $      chargesort)

        subroutine chunk_reorder_zvals(k, nch, iorder, sigma, sigmaout)
        implicit real *8 (a-h,o-z)
        integer :: iorder(nch)
        complex *16 :: sigma(k,nch), sigmaout(k,nch)
c
c       this routine reorders the values contained in sigma according
c       to the list in iorder
c
cccc        call prinf('iorder = *', iorder, nch)
        
        do i = 1,nch
          ind = iorder(i)
          do j = 1,k
            sigmaout(j,i) = sigma(j,ind)
          enddo
        enddo
c
        return
        end
c       
c        
c
c
c        
        subroutine chunk_dist(ich, wgeo, targ, dist, xy, xy_norm)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1), targ(2), xy(2), xy_norm(2)
        real *8 xnodes(1000), whts(1000)
        real *8, allocatable :: u(:,:), v(:,:)
c
c       calculate the distance from chunk ich to targ - this routine 
c       is just a wrapper for chunk_dist1
c
c       input:
c         ich - the number number
c         wgeo - geometry structure created by chunkfunc or some
c             similar routine
c         targ - the distance from chunk ich to targ is desired
c
c       output:
c         dist - the distance from the chunk to targ
c         xy - the point lying on the chunk that is the closest, this
c             point is determined by newton
c
c       first get the address of the chunk inside wgeo
c
        call chunkaddr(ich, wgeo, k, ichunk, iadj, ider,
     1      ider2, h)
c
cccc        ifwhts = 0
cccc        call legewhts(k, xnodes, whts, ifwhts)

        itype = 2
        allocate(u(k,k), v(k,k))
        call legeexps(itype, k, xnodes, u, v, whts)
c
        call chunk_dist1(k, xnodes, u, wgeo(ichunk), wgeo(ider), 
     1      wgeo(ider2), h, targ, dist, xy, xy_norm)
c
        return
        end
c
c
c
c
c
        subroutine chunkaddr(ich, wgeo, k, ichunk, iadj, ider,
     1      ider2, h)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1)
c
c       get the address of all the arrays stored in wgeo
c
c       input:
c         ich - the chunk number
c         wgeo - the storage array
c
c       all others are output values
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1      iders2, ihs)
c
        ichunk = ichunks + (ich-1)*2*k
        ider = iders + (ich-1)*2*k
        ider2 = iders2 + (ich-1)*2*k
        h = wgeo(ihs+(ich-1))
c
        return
        end
c
c
c
c
c
        subroutine chunk_dist1(k, xnodes, u, chunk, der, der2, h, 
     1      targ, dist, xy, xy_norm)
        implicit real *8 (a-h,o-z)
        real *8 chunk(2,k), targ(2), xy(2), xy_norm(2)
        real *8 xcoefs(1000), ycoefs(1000), xnodes(*), u(k,k)
        real *8 xpcoefs(1000), ypcoefs(1000)
        real *8 xppcoefs(1000), yppcoefs(1000)
c
c       unpackaged routine for computing the distance from chunk
c       to targ
c
c       input:
c         k - number of points per chunk
c         xnodes - the k Legendre points
c         u - the matrix created by legewhts mapping points to coefs
c         chunk - points on the chunk at legendre nodes
c         der - derivatices w.r.t. "some" parameterization
c         der2 - 2nd derivatives w.r.t. the same parameterization
c         h - normalizing factor for the parameterization
c         targ - the distance from chunk ich to targ is desired
c
c       output:
c         dist - the distance from the chunk to targ
c         xy - the point lying on the chunk that is the closest, this
c             point is determined by newton
c         xy_norm - the outward normal, i.e. into the unbounded domain,
c             assuming that the chunk is parameterized counter-clockwise
c
c
        done = 1
c
c       find the closest legendre node
c
cccc        ifwhts = 0
cccc        call legewhts(k, xnodes, whts, ifwhts)

c
        dist = 1.0d15
c
        do i = 1,k
          x = chunk(1,i)
          y = chunk(2,i)
          d = sqrt((x-targ(1))**2 + (y-targ(2))**2)
          if (d .lt. dist) then
            tt = xnodes(i)
            dist = d
            xy(1) = x
            xy(2) = y
          endif
        enddo
c
c       check the endpoints
c
        call chunksexps_fast(k, u, chunk, xcoefs, ycoefs)
        ifend = 0
c
        t0 = -1
        call legeexev(t0, x0, xcoefs, k-1)
        call legeexev(t0, y0, ycoefs, k-1)
        d0 = sqrt((targ(1) - x0)**2 + (targ(2) - y0)**2)
c
        t1 = 1
        call legeexev(t1, x1, xcoefs, k-1)
        call legeexev(t1, y1, ycoefs, k-1)
        d1 = sqrt((targ(1) - x1)**2 + (targ(2) - y1)**2)
c
        if (d0 .lt. dist) then
          tt = -1
          ifend = 1
          dist = d0
          xy(1) = x0
          xy(2) = y0
        endif
c
        if (d1 .lt. dist) then
          tt = 1
          ifend = 1
          dist = d1
          xy(1) = x1
          xy(2) = y1
        endif


c
c       if an endpoint is the closest, then return
c
        if (ifend .eq. 1) then
          call chunksexps_fast(k, u, der, xpcoefs, ypcoefs)
          call legeexev(tt, xp, xpcoefs, k-1)
          call legeexev(tt, yp, ypcoefs, k-1)
          dd = sqrt(xp**2+yp**2)
          xy_norm(1) = yp/dd
          xy_norm(2) = -xp/dd
          return
        endif
c
c       otherwise run newton
c
        call chunksexps_fast(k, u, der, xpcoefs, ypcoefs)
        call chunksexps_fast(k, u, der2, xppcoefs, yppcoefs)
c
        t0 = tt
        maxnewt = 200
        thresh = 1.0d-10
        thresh = 1.0d-9
        iextra = 3
        ifend = 0

        do i = 1,maxnewt
c
          call legeexev(t0, x, xcoefs, k-1)
          call legeexev(t0, y, ycoefs, k-1)
          call legeexev(t0, xp, xpcoefs, k-1)
          call legeexev(t0, yp, ypcoefs, k-1)
          call legeexev(t0, xpp, xppcoefs, k-1)
          call legeexev(t0, ypp, yppcoefs, k-1)
c
          dx = x - targ(1)
          dy = y - targ(2)
          d = dx**2 + dy**2
c    
          xp = xp*h
          yp = yp*h
          xpp = xpp*h*h
          ypp = ypp*h*h
          dp = 2*(dx*xp+dy*yp)
          dpp = 2*(xp**2 + dx*xpp + yp**2 + dy*ypp)
          t1 = t0 - dp/dpp
c
          if (t1 .gt. done) then
            t1 = (done+t0)/2
          endif
c
          if (t0 .lt. -done) then
            t1 = (-done+t0)/2
          endif
c
cccc          err = abs(t1 - t0)/abs(t1)
          err = abs(t1 - t0)
c
          if (err .lt. thresh) then
            ifend = ifend + 1
          endif
c
          t0 = t1
          xy(1) = x
          xy(2) = yx
          dist = d
          if (ifend .eq. iextra) goto 1100
c
        enddo
c
 1100 continue
c
        if (ifend .lt. iextra) then
          call prin2('newton bomb in chunk_dist1!*', x, 0)
          call prin2('t1 = *', t1, 1)
          call prin2('err = *', err, 1)
          call prin2('thresh = *', thresh, 1)
          stop
        endif
c
        if (err .gt. thresh/10) then
c          call prin2('WARNING: Large error in chunk_dist1!*', x, 0)
c          call prin2('t1 = *', t1, 1)
c          call prin2('err = *', err, 1)
c          call prin2('thresh = *', thresh, 1)
        endif
c
c       if here, then success - now return
c     
        call legeexev(t1, x, xcoefs, k-1)
        call legeexev(t1, y, ycoefs, k-1)
        dist = sqrt((targ(1) - x)**2 + (targ(2) - y)**2)
        xy(1) = x
        xy(2) = y
c
c       calculate the unit outward normal at t1
c

cc        call prin2('t1 in chunk_dist1=*',t1,1)
        call legeexev(t1, xp, xpcoefs, k-1)
        call legeexev(t1, yp, ypcoefs, k-1)
        dd = sqrt(xp**2+yp**2)
        xy_norm(1) = yp/dd
        xy_norm(2) = -xp/dd
c
        return
        end
c
c
c-------------------------------------------------------------c
c
        subroutine chunk_dist2(k,xnodes,chunk,xcoefs, ycoefs, xpcoefs,
     1      ypcoefs,xppcoefs,yppcoefs, h, targ, dist, xy, xy_norm)
        implicit real *8 (a-h,o-z)
        real *8 targ(2), xy(2), xy_norm(2)
        real *8 xcoefs(*), ycoefs(*) 
        real *8 xpcoefs(*), ypcoefs(*)
        real *8 xppcoefs(*), yppcoefs(*)
        real *8 chunk(2,*),xnodes(*)
c
c       unpackaged routine for computing the distance from chunk
c       to targ. This subroutine is the same as 
c       chunk_dist1. The only difference is that for this
c       subroutine we pass on the coefficients of the
c       legendre expansion for the coordinates and the
c       derivatives of the chunks as opposed to the point
c       values.
c
c       input:
c         k - number of points per chunk
c         xnodes - the k Legendre points
c         u - the matrix created by legewhts mapping points to coefs
c         chunk - points on the chunk at legendre nodes
c         der - derivatices w.r.t. "some" parameterization
c         der2 - 2nd derivatives w.r.t. the same parameterization
c         h - normalizing factor for the parameterization
c         targ - the distance from chunk ich to targ is desired
c
c       output:
c         dist - the distance from the chunk to targ
c         xy - the point lying on the chunk that is the closest, this
c             point is determined by newton
c         xy_norm - the outward normal, i.e. into the unbounded domain,
c             assuming that the chunk is parameterized counter-clockwise
c
c
        done = 1
c
c       find the closest legendre node
c
cccc        ifwhts = 0
cccc        call legewhts(k, xnodes, whts, ifwhts)

c
        dist = 1.0d15
c
        do i = 1,k
          x = chunk(1,i)
          y = chunk(2,i)
          d = sqrt((x-targ(1))**2 + (y-targ(2))**2)
          if (d .lt. dist) then
            tt = xnodes(i)
            dist = d
            xy(1) = x
            xy(2) = y
          endif
        enddo
c
c       check the endpoints
c
cccc        call chunksexps_fast(k, u, chunk, xcoefs, ycoefs)
        ifend = 0
c
        t0 = -1
        call legeexev(t0, x0, xcoefs, k-1)
        call legeexev(t0, y0, ycoefs, k-1)
        d0 = sqrt((targ(1) - x0)**2 + (targ(2) - y0)**2)
c
        t1 = 1
        call legeexev(t1, x1, xcoefs, k-1)
        call legeexev(t1, y1, ycoefs, k-1)
        d1 = sqrt((targ(1) - x1)**2 + (targ(2) - y1)**2)
c
        if (d0 .lt. dist) then
          tt = -1
          ifend = 1
          dist = d0
          xy(1) = x0
          xy(2) = y0
        endif
c
        if (d1 .lt. dist) then
          tt = 1
          ifend = 1
          dist = d1
          xy(1) = x1
          xy(2) = y1
        endif


c
c       if an endpoint is the closest, then return
c


cc        call prinf('ifend=*',ifend,1)
c
        if (ifend .eq. 1) then

cccc          call chunksexps_fast(k, u, der, xpcoefs, ypcoefs)
          call legeexev(tt, xp, xpcoefs, k-1)
          call legeexev(tt, yp, ypcoefs, k-1)
          dd = sqrt(xp**2+yp**2)
          xy_norm(1) = yp/dd
          xy_norm(2) = -xp/dd
          return
        endif
c
c       otherwise run newton
c
ccc        call chunksexps_fast(k, u, der, xpcoefs, ypcoefs)
ccc        call chunksexps_fast(k, u, der2, xppcoefs, yppcoefs)
c
        t0 = tt
        maxnewt = 200
        thresh = 1.0d-10
        thresh = 1.0d-9
        iextra = 3
        ifend = 0
cc        call prinf('k=*',k,1)

cc        call prin2('xcoefs=*',xcoefs,k)
cc        call prin2('ycoefs=*',ycoefs,k)
cc        call prin2('xpcoefs=*',xpcoefs,k)
cc        call prin2('ypcoefs=*',ypcoefs,k)
cc        call prin2('xppcoefs=*',xppcoefs,k)
cc        call prin2('yppcoefs=*',yppcoefs,k)

        do i = 1,maxnewt
c
          call legeexev(t0, x, xcoefs, k-1)
          call legeexev(t0, y, ycoefs, k-1)
          call legeexev(t0, xp, xpcoefs, k-1)
          call legeexev(t0, yp, ypcoefs, k-1)
          call legeexev(t0, xpp, xppcoefs, k-1)
          call legeexev(t0, ypp, yppcoefs, k-1)
c
          dx = x - targ(1)
          dy = y - targ(2)
          d = dx**2 + dy**2
c    
          xp = xp*h
          yp = yp*h
          xpp = xpp*h*h
          ypp = ypp*h*h
          dp = 2*(dx*xp+dy*yp)
          dpp = 2*(xp**2 + dx*xpp + yp**2 + dy*ypp)
          t1 = t0 - dp/dpp
c
          if (t1 .gt. done) then
            t1 = (done+t0)/2
          endif
c
          if (t0 .lt. -done) then
            t1 = (-done+t0)/2
          endif
c
cccc          err = abs(t1 - t0)/abs(t1)
          err = abs(t1 - t0)
c
          if (err .lt. thresh) then
            ifend = ifend + 1
          endif
c
          t0 = t1
          xy(1) = x
          xy(2) = yx
          dist = d
          if (ifend .eq. iextra) goto 1100
c
        enddo
c
 1100 continue
c
        if (ifend .lt. iextra) then
          call prin2('newton bomb in chunk_dist1!*', x, 0)
          call prin2('t1 = *', t1, 1)
          call prin2('dp = *',dp,1)
          call prin2('dpp = *',dpp,1)          
          call prin2('err = *', err, 1)
          call prin2('thresh = *', thresh, 1)
          stop
        endif
c
        if (err .gt. thresh/10) then
c          call prin2('WARNING: Large error in chunk_dist1!*', x, 0)
c          call prin2('t1 = *', t1, 1)
c          call prin2('err = *', err, 1)
c          call prin2('thresh = *', thresh, 1)
        endif
c
c       if here, then success - now return
c     
cc        call prin2('t1 in chunk_dist2=*',t1,1)
        call legeexev(t1, x, xcoefs, k-1)
        call legeexev(t1, y, ycoefs, k-1)
        dist = sqrt((targ(1) - x)**2 + (targ(2) - y)**2)
        xy(1) = x
        xy(2) = y
c
c       calculate the unit outward normal at t1
c
        call legeexev(t1, xp, xpcoefs, k-1)
        call legeexev(t1, yp, ypcoefs, k-1)
        dd = sqrt(xp**2+yp**2)
        xy_norm(1) = yp/dd
        xy_norm(2) = -xp/dd
c
        return
        end
c
c---------------------------------------------------------------c
        subroutine chunksexps(k, chunk, xcoefs, ycoefs)
        implicit real *8 (a-h,o-z)
        real *8 chunk(2,k), xcoefs(k), ycoefs(k)
        real *8 ts(1000), whts(1000), u(10000), v(10000)
        real *8 xs(1000), ys(1000)
c
c       this routine calculates the legendre expansions
c       of the points in chunk (chunk can be any 2,k vector
c       in fact)
c
        itype = 2
        call legeexps(itype, k, ts, u, v, whts)
c
        do i = 1,k
          xs(i) = chunk(1,i)
          ys(i) = chunk(2,i)
        enddo
c
        call chunkmatvec(k, k, u, xs, xcoefs)
        call chunkmatvec(k, k, u, ys, ycoefs)
c
        return
        end
c
c
c
c
c
        subroutine chunksexps_fast(k, u, chunk, xcoefs, ycoefs)
        implicit real *8 (a-h,o-z)
        real *8 chunk(2,k), xcoefs(k), ycoefs(k), u(k,k)
        real *8 ts(1000), whts(1000), v(10000)
        real *8 xs(1000), ys(1000)
c
c       this routine calculates the legendre expansions
c       of the points in chunk (chunk can be any 2,k vector
c       in fact)
c
cccc        itype = 2
cccc        call legeexps(itype, k, ts, u, v, whts)
c
        do i = 1,k
          xs(i) = chunk(1,i)
          ys(i) = chunk(2,i)
        enddo
c
        call chunkmatvec(k, k, u, xs, xcoefs)
        call chunkmatvec(k, k, u, ys, ycoefs)
c
        return
        end
c
c
c
c
c
        subroutine chunkdsdt(wgeo, dsdt)
        implicit double precision (a-h,o-z)
        double precision wgeo(1), dsdt(1)
c
c       calculates the value of dsdt at each point
c       on the chunk
c
        done = 1
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        ijk = 0
        ind = 0
        do i = 1,nch
          h = wgeo(ihs+i-1)
          do j = 1,k
            ind = ind+1
            dxdt = wgeo(iders+ijk)
            dydt = wgeo(iders+ijk+1)
            dsdt(ind) = sqrt(dxdt**2+dydt**2)*h
            ijk = ijk+2
          enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkrls(wgeo, rls, rltot)
        implicit double precision (a-h,o-z)
        double precision wgeo(1), rls(1), coefs(10000)
        double precision xs(10000), whts(10000)
        double precision u(100000), v(100000), coefs_out(10000)
c
        double precision, allocatable :: dsdt(:,:), rtemp(:,:)
        integer, allocatable :: adjs(:,:)
c
c       this routine returns in rls the stepwise arclengths at
c       each point - useful for plotting curvatures, densities, etc.
c
        done = 1
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        do i = 1,k*nch
          rls(i) = 0
        enddo
c
        itype = 2
        call legeexps(itype, k, xs, u, v, whts)

c
        ninre = sizeof(done)/sizeof(k)
        allocate(adjs(2,nch))
        call chunkcopy(nch*2/ninre, wgeo(iadjs), adjs)
c
        allocate(dsdt(k,nch))
        allocate(rtemp(k,nch))
        call chunkdsdt(wgeo, dsdt)

c
c       scan through the chunks, calcuate the legendre
c       series of dsdt, then integrate it on the interval [-1,x_i]
c
        do i = 1,nch
c
          call chunkmatvec(k, k, u, dsdt(1,i), coefs)
          call legeinte(coefs, k-1, coefs_out)
c
c         now eval coefs_out at each of the k legendre nodes
c
          do j = 1,k
            call legeexev(xs(j), rtemp(j,i), coefs_out, k)
          enddo
c
c         and evaluate the total arclength of the chunk
c
          call legeexev(done, rls(i), coefs_out, k)
c
        enddo

c
c       sum the arclengths of the chunks
c
        rltot = 0
        do i = 1,nch
c
          if (i .eq. 1) ichunk = 1
          if (i .ne. 1) ichunk = adjs(2,ichunk)
c
          call prinf('ichunk = *', ichunk, 1)
          call prin2('rtemp = *', rtemp(1,ichunk), k)
          do j = 1,k
            rtemp(j,ichunk) = rtemp(j,ichunk) + rltot
          enddo
c
          rltot = rltot + rls(ichunk)
c
        enddo
c
        call chunkcopy(k*nch, rtemp, rls)
c
        return
        end
c
c
c
c
c
        subroutine chunkcurvatures(wgeo, curvatures)
        implicit double precision (a-h,o-z)
        double precision wgeo(1), curvatures(1)
c
c       at each point on all chunks, calculate the signed curvature
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        ijk = 0
        ind = 0
        do i = 1,nch
          do j = 1,k
            dxdt = wgeo(iders+ijk)
            dydt = wgeo(iders+ijk+1)
            dxdt2 = wgeo(iders2+ijk)
            dydt2 = wgeo(iders2+ijk+1)
            dsdt = sqrt(dxdt**2+dydt**2)
            ijk = ijk+2
            ind = ind+1
            curvatures(ind) = (dxdt*dydt2-dydt*dxdt2)/dsdt**3
          enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunksplot(iw, itype, wgeo)
        implicit double precision (a-h,o-z)
        double precision wgeo(1)
        double precision, allocatable :: normals(:,:,:)
c
c       plots the chunks that are stored in wgeo
c
c       itype = 1   -   plot just the chunks
c       itype = 2   -   plot the chunks and the normals
c
        done = 1
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        if (itype .eq. 1) then
          idot = 1
          nnn = k*nch
          call zpyplot(iw, wgeo(ichunks), nnn, idot, 'The chunks*')
          return
        endif
c
        if (itype .eq. 2) then
c
           allocate(normals(2,k,nch))
           call chunknormals(wgeo, normals)
           ijk = 0
c
           sc = 1
           do i = 1,nch
             do j = 1,k
               normals(1,j,i) = wgeo(ichunks+ijk) + normals(1,j,i)/sc
               ijk = ijk+1
               normals(2,j,i) = wgeo(ichunks+ijk) + normals(2,j,i)/sc
               ijk = ijk+1
             enddo
           enddo
c
          idot = 1
          nnn = k*nch
          call zpyplot2(iw, wgeo(ichunks), nnn, idot, 
     1      normals, nnn, idot, 'The chunks and normals*')
c
        endif


c
        return
        end
c
c
c
c
c
        subroutine chunknormals(wgeo, normals)
        implicit double precision (a-h,o-z)
        double precision wgeo(1), normals(2,1)
c
c       this subroutine returns all the chunk normals at once,
c       assuming a counter-clockwise orientation
c
c       input:
c           wgeo - the geometry structure
c
c       output:
c           normals - the UNIT outward normals
c
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        ijk = 0
        do i = 1,nch
          do j = 1,k
            ijk = ijk+1
            dxdt = wgeo(iders+2*(ijk-1))
            dydt = wgeo(iders+2*(ijk-1)+1)
            dd = sqrt(dxdt**2 + dydt**2)
            normals(1,ijk) = dydt/dd
            normals(2,ijk) = -dxdt/dd
          enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkfunc_pack(eps, ifclosed, chsmall, ta, tb,
     1    funcurve, pars, k, nover, wgeo, lused)
        implicit double precision (a-h,o-z)
        double precision pars(1), wgeo(1)
        external funcurve
c
        double precision, allocatable :: chunks(:,:,:), ders(:,:,:)
        double precision, allocatable :: ders2(:,:,:), hs(:)
        integer, allocatable :: adjs(:,:)
c
c       this routine calls chunkfunc and then packs the
c       resulting arrays into a work array, wgeo
c
        maxpts = 10000000
        maxchunks = maxpts/k
c
        allocate(chunks(2,k,maxchunks))
        allocate(ders(2,k,maxchunks))
        allocate(ders2(2,k,maxchunks))
        allocate(hs(maxchunks))
        allocate(adjs(2,maxchunks))
c
        call chunkfunc(eps, ifclosed, chsmall, ta, tb, funcurve,
     1    pars, nover, k, nch, chunks, adjs, ders, ders2, hs)
c
        if (nch .gt. maxchunks) then
          call prinf('in chunkfunc_pack, nch = *', nch, 1)
          call prinf('in chunkfunc_pack, maxchunks = *', maxchunks, 1)
          stop
        endif
c
        call chunkpack(k, nch, chunks, adjs, ders, ders2,
     1    hs, wgeo, lused)
c
        return
        end
c
c
c
c
c
        subroutine chunkreverse(k,nch,chunks,adjs,
     1      ders,ders2,hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1)
c
c       this routine reverses the orientation of the curve described
c       by chunks, adjs, ders, ders2, and hs
c
c       reversing the orientation means:
c
c           chunks(,j,i) -> chunks(,k-j+1,i)
c           ders(,j,i)   -> -ders(,k-j+1,i)
c           ders2(,j,i)  -> ders2(,k-j+1,i)
c           adjs(1,i)   <-> adjs(2,i)
c
c       no chunk numbers are changed!! they are arbitrary to begin with
c
c
c       . . . first switch the order of the points on each chunk
c
        khalf=k/2
c
        do 4000 i=1,nch
c
        do 3400 j=1,khalf
c
        temp=chunks(1,j,i)
        chunks(1,j,i)=chunks(1,k-j+1,i)
        chunks(1,k-j+1,i)=temp
c
        temp=chunks(2,j,i)
        chunks(2,j,i)=chunks(2,k-j+1,i)
        chunks(2,k-j+1,i)=temp
c
        temp=ders(1,j,i)
        ders(1,j,i)=ders(1,k-j+1,i)
        ders(1,k-j+1,i)=temp
c
        temp=ders(2,j,i)
        ders(2,j,i)=ders(2,k-j+1,i)
        ders(2,k-j+1,i)=temp
c
        temp=ders2(1,j,i)
        ders2(1,j,i)=ders2(1,k-j+1,i)
        ders2(1,k-j+1,i)=temp
c
        temp=ders2(2,j,i)
        ders2(2,j,i)=ders2(2,k-j+1,i)
        ders2(2,k-j+1,i)=temp
c
 3400 continue

c
c       . . . flip the sign on the first derivative
c
        do 3600 j=1,k
        ders(1,j,i)=-ders(1,j,i)
        ders(2,j,i)=-ders(2,j,i)
 3600 continue
c
 4000 continue

c
c       now re-number the chunks from right to left, just transpose adjs
c
        do 4200 i=1,nch
        iii=adjs(1,i)
        adjs(1,i)=adjs(2,i)
        adjs(2,i)=iii
 4200 continue
c
        return
        end
c
c
c
c
c
        subroutine chunksort_pack(wgeo)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1)
c
c       sort the geometry *directly* in the packaged form
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1      iders2, ihs)
c
        call chunksort(k, nch, wgeo(ichunks), wgeo(iadjs),
     1      wgeo(iders), wgeo(iders2), wgeo(ihs))
c
        return
        end
c
c
c
c
c
        subroutine chunksort(k,nch,chunks,adjs,ders,ders2,hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1)
c
        real *8, allocatable :: c(:,:,:),d(:,:,:),d2(:,:,:),
     1      h(:)
c
c       this routine follows adjs all the way through and
c       reorders all the chunks so that the points are in a row,
c       positively oriented (assuming that each chunk by itself
c       is positively oriented)
c
c       this routine also assumes that chunks contains only one curve,
c       i.e. not two disjoint ones
c
        allocate(c(2,k,nch))
        allocate(d(2,k,nch))
        allocate(d2(2,k,nch))
        allocate(h(nch))
c
c       determine if the curve is closed or not...
c
        ifclosed = 1
        ileft = 0
        iright = 0
c
        do i = 1, nch
            if (adjs(1,i) .lt. 0) then
                ileft = i
                ifclosed = 0
            endif
            if (adjs(2,i) .lt. 0) then
                iright = i
                ifclosed = 0
            endif
        enddo
c
        icurrent=1
        if (ifclosed .eq. 0) icurrent = ileft
c
        do i = 1, nch
c
            if (i .ne. 1) icurrent=adjs(2,ilast)
            ilast=icurrent
c
            h(i)=hs(icurrent)
c
            do j = 1, k
                c(1,j,i)=chunks(1,j,icurrent)
                c(2,j,i)=chunks(2,j,icurrent)
                d(1,j,i)=ders(1,j,icurrent)
                d(2,j,i)=ders(2,j,icurrent)
                d2(1,j,i)=ders2(1,j,icurrent)
                d2(2,j,i)=ders2(2,j,icurrent)
            enddo
c
        enddo

c
c       copy them back over
c
        do 3600 i=1,nch
c
        hs(i)=h(i)
c
        do 3400 j=1,k
        chunks(1,j,i)=c(1,j,i)
        chunks(2,j,i)=c(2,j,i)
        ders(1,j,i)=d(1,j,i)
        ders(2,j,i)=d(2,j,i)
        ders2(1,j,i)=d2(1,j,i)
        ders2(2,j,i)=d2(2,j,i)
 3400 continue
 3600 continue

c
c       . . . now update adjs
c
        if (ifclosed .eq. 0) then
            adjs(1,1) = -1
            adjs(2,nch) = -1
        endif
c
        if (ifclosed .eq. 1) then
            adjs(1,1) = nch
            adjs(2,nch) = 1
        endif
c
        do i = 1, nch-1
            adjs(2,i) = i+1
            adjs(1,i+1) = i
        enddo
c
        return
        end
c
c
c
c
        subroutine chunksortab(k,nch,chunks,adjs,ders,ders2,hs,ab)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1),ab(2,1)
c
        real *8, allocatable :: c(:,:,:),d(:,:,:),d2(:,:,:),
     1      h(:),ab2(:,:)
c
c       this routine follows adjs all the way through and
c       reorders all the chunks so that the points are in a row,
c       positively oriented (assuming that each chunk by itself
c       is positively oriented)
c
c       this routine also assumes that chunks contains only one curve,
c       i.e. not two disjoint ones
c
        allocate(c(2,k,nch))
        allocate(d(2,k,nch))
        allocate(d2(2,k,nch))
        allocate(h(nch),ab2(2,nch))
c
c       determine if the curve is closed or not...
c
        ifclosed = 1
        ileft = 0
        iright = 0
c
        do i = 1, nch
            if (adjs(1,i) .lt. 0) then
                ileft = i
                ifclosed = 0
            endif
            if (adjs(2,i) .lt. 0) then
                iright = i
                ifclosed = 0
            endif
        enddo
c
        icurrent=1
        if (ifclosed .eq. 0) icurrent = ileft
c
        do i = 1, nch
c
            if (i .ne. 1) icurrent=adjs(2,ilast)
            ilast=icurrent
c
            h(i)=hs(icurrent)
            ab2(1,i) = ab(1,icurrent)
            ab2(2,i) = ab(2,icurrent)
c
            do j = 1, k
                c(1,j,i)=chunks(1,j,icurrent)
                c(2,j,i)=chunks(2,j,icurrent)
                d(1,j,i)=ders(1,j,icurrent)
                d(2,j,i)=ders(2,j,icurrent)
                d2(1,j,i)=ders2(1,j,icurrent)
                d2(2,j,i)=ders2(2,j,icurrent)
            enddo
c
        enddo

c
c       copy them back over
c
        do 3600 i=1,nch
c
        hs(i)=h(i)
        ab(1,i) = ab2(1,i)
        ab(2,i) = ab2(2,i)
c
        do 3400 j=1,k
        chunks(1,j,i)=c(1,j,i)
        chunks(2,j,i)=c(2,j,i)
        ders(1,j,i)=d(1,j,i)
        ders(2,j,i)=d(2,j,i)
        ders2(1,j,i)=d2(1,j,i)
        ders2(2,j,i)=d2(2,j,i)
 3400 continue
 3600 continue

c
c       . . . now update adjs
c
        if (ifclosed .eq. 0) then
            adjs(1,1) = -1
            adjs(2,nch) = -1
        endif
c
        if (ifclosed .eq. 1) then
            adjs(1,1) = nch
            adjs(2,nch) = 1
        endif
c
        do i = 1, nch-1
            adjs(2,i) = i+1
            adjs(1,i+1) = i
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkres(k,nch,chunks,ders,hs,sigma,errmax1,
     1      errmax2,tails)
        implicit real *8 (a-h,o-z)
        real *8 chunks(2,1,1),ders(2,1,1),hs(1),whts(1000),
     1      u(10000),v(10000),xs(1000),tails(1),xymid(10),
     2      cmass(10)
        real *8 sigma(k,1),coefs(1000)
c
c       this routine checks the resolution of a real-valued
c       function sigma along the curve chunks
c
c       input:
c
c         k - number of points per chunk
c         nch - number of chunks
c         chunks - chunk info, see chunkfunc for more info
c         ders - derivative info, see chunkfunc for more info
c         hs - scaling info, see chunkfunc for more info
c         sigma - complex valued function on chunks
c
c       output:
c
c         errmax1 - maximum tail size, scaled by the arclength of
c             the chunk
c         errmax2 - maximum relative tail size
c         tails - the size of the tails on each chunk, scaled by
c             the arclength of the chunk
c
c       NOTE: see below, but the size of the 'tail' is computed in
c       the l2 sense, and includes a varying number of coefficients
c       depending on the expansion size k
c
        done=1
c
        itype=2
        call legeexps(itype,k,xs,u,v,whts)
c
        if (k .eq. 1) nterms=1
        if (k .eq. 2) nterms=1
        if (k .eq. 3) nterms=1
        if (k .eq. 4) nterms=2
        if (k .eq. 5) nterms=2
        if (k .eq. 6) nterms=2
        if (k .eq. 7) nterms=2
        if (k .eq. 8) nterms=2
        if (k .eq. 9) nterms=2
        if (k .eq. 10) nterms=3
        if (k .eq. 11) nterms=3
        if (k .eq. 12) nterms=3
c
        if (k .gt. 12) nterms=2+k/12

c
c       loop over all chunks, and convert the values in
c       sigma to coefficients and compute the tails
c
        errmax1=-1
        errmax2=-1
c
        do 2800 i=1,nch
c
        call chunkmatvec(k,k,u,sigma(1,i),coefs)

c
c       compute the l2 norm of the coefficients
c
        rnorm=0
        do 1600 j=1,k
        rnorm=rnorm+abs(coefs(j))**2
 1600 continue
c
        rnorm=sqrt(rnorm)

c
c       compute the l2 norm of the tail
c
        rtail=0
        do 2000 j=1,nterms
        ioffset=j-1
        rtail=rtail+abs(coefs(k-ioffset))**2
 2000 continue
c
        rtail=sqrt(rtail)

c
c       weight the absolute tail value by the chunk size
c
        call chunksize(k,chunks(1,1,i),ders(1,1,i),hs(i),xymid,
     1      rad,cmass,crad,rltot)
c
        err1=rtail*rltot
        err2=rtail/rnorm
        tails(i)=err1
c
        if (err1 .gt. errmax1) errmax1=err1
        if (err2 .gt. errmax2) errmax2=err2
c
 2800 continue


c
        return
        end
c
c
c
c
c
        subroutine chunkzres(k,nch,chunks,ders,hs,sigma,errmax1,
     1      errmax2,tails)
        implicit real *8 (a-h,o-z)
        real *8 chunks(2,1,1),ders(2,1,1),hs(1),whts(1000),
     1      u(10000),v(10000),xs(1000),tails(1),xymid(10),
     2      cmass(10)
        complex *16 sigma(k,1),coefs(1000)
c
c       this routine checks the resolution of a complex-valued
c       function sigma along the curve chunks
c
c       input:
c
c         k - number of points per chunk
c         nch - number of chunks
c         chunks - chunk info, see chunkfunc for more info
c         ders - derivative info, see chunkfunc for more info
c         hs - scaling info, see chunkfunc for more info
c         sigma - complex valued function on chunks
c
c       output:
c
c         errmax1 - maximum tail size, scaled by the arclength of
c             the chunk
c         errmax2 - maximum relative tail size
c         tails - the size of the tails on each chunk, scaled by
c             the arclength of the chunk
c
c       NOTE: see below, but the size of the 'tail' is computed in
c       the l2 sense, and includes a varying number of coefficients
c       depending on the expansion size k
c
        done=1
c
        itype=2
        call legeexps(itype,k,xs,u,v,whts)
c
        if (k .eq. 1) nterms=1
        if (k .eq. 2) nterms=1
        if (k .eq. 3) nterms=1
        if (k .eq. 4) nterms=2
        if (k .eq. 5) nterms=2
        if (k .eq. 6) nterms=2
        if (k .eq. 7) nterms=2
        if (k .eq. 8) nterms=2
        if (k .eq. 9) nterms=2
        if (k .eq. 10) nterms=3
        if (k .eq. 11) nterms=3
        if (k .eq. 12) nterms=3
c
        if (k .gt. 12) nterms=2+k/12
cccc        call prinf('nterms = *', nterms, 1)
cccc        stop

c
c       loop over all chunks, and convert the values in
c       sigma to coefficients and compute the tails
c
        errmax1=-1
        errmax2=-1
c
        do i=1,nch
c
           call chunkmatvec2(k,k,u,sigma(1,i),coefs)
ccc           call prin2('coefs = *', coefs, 2*k)

c
c       compute the l2 norm of the coefficients
c     
           rnorm=0
           do j=1,k
              rnorm=rnorm+abs(coefs(j))**2
           enddo
c
           rnorm=sqrt(rnorm)
ccc           call prin2('rnorm = *', rnorm, 1)

c
c       compute the l2 norm of the tail
c
           rtail=0
           do j=1,nterms
              ioffset=j-1
              rtail=rtail+abs(coefs(k-ioffset))**2
ccc              call prin2('ce = *', coefs(k-ioffset), 2)
           enddo
c
           rtail=sqrt(rtail)
ccc           call prin2('rtail = *', rtail, 1)

c
c       weight the absolute tail value by the chunk size
c
           call chunksize(k,chunks(1,1,i),ders(1,1,i),hs(i),
     1          xymid,rad,cmass,crad,rltot)
c
           err1=rtail*rltot
           err2=rtail/rnorm*rltot

cccc           call prin2('err1 = *', err1, 1)
cccc           call prin2('err2 = *', err2, 1)
cccc           stop

           tails(i)=err1
c     
           if (err1 .gt. errmax1) errmax1=err1
           if (err2 .gt. errmax2) errmax2=err2
c     
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkwhts(k,nch,chunks,ders,hs,whts)
        implicit real *8 (a-h,o-z)
        real *8 chunks(2,k,1),ders(2,k,1),hs(1),whts(k,1),
     1      ts(1000),ws(1000)
c
c       this routine returns the smooth quadrature weights for
c       integrating functions on chunks - these weights 
c       contain contributions from ds/dt on the way out.
c
        done=1
        ifwhts=1
        call legewhts(k,ts,ws,ifwhts)
c
        do i = 1,nch
          do j = 1,k
            dsdt=sqrt(ders(1,j,i)**2+ders(2,j,i)**2)
            dsdt=dsdt*hs(i)
            whts(j,i)=ws(j)*dsdt
          enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkused(wgeo, lused)
        implicit real *8 (a-h,o-z)
        double precision wgeo(*)
c
c       return the length of memory used in wgeo
c
        lused = wgeo(8)
c
        return
        end
c
c
c
        
        subroutine chunkpack(k,nch,chunks,adjs,ders,ders2,
     1      hs,wgeo,lused)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        double precision chunks(2,k,1), ders(2,k,1), ders2(2,k,1)
        double precision hs(1), wgeo(1)
c
c       this routine packs all chunk info into one array, wpack
c
        d = 1.0d-0
        m = 5
        ninre = sizeof(d)/sizeof(m)
c
        wgeo(1)=k+.1d0
        wgeo(2)=nch+.1d0
c
        ichunks=11
        lchunks=2*k*nch
c
        iadjs=ichunks+lchunks+100
        ladjs=2*nch/ninre
c
        iders=iadjs+ladjs+100
        lders=2*k*nch
c
        iders2=iders+lders+100
        lders2=2*k*nch
c
        ihs=iders2+lders2+100
        lhs=nch
        lused=ihs+lhs+100
c
        wgeo(3)=ichunks+.1d0
        wgeo(4)=iadjs+.1d0
        wgeo(5)=iders+.1d0
        wgeo(6)=iders2+.1d0
        wgeo(7)=ihs+.1d0
        wgeo(8)=lused+.1d0
c
c       copy the data
c
        call chunkcopy(lchunks,chunks,wgeo(ichunks))
        call chunkcopy(ladjs,adjs,wgeo(iadjs))
        call chunkcopy(lders,ders,wgeo(iders))
        call chunkcopy(lders2,ders2,wgeo(iders2))
        call chunkcopy(lhs,hs,wgeo(ihs))
c
        return
        end
c
c
c
c
c
        subroutine chunkcopy(n,x,y)
        implicit real *8 (a-h,o-z)
        real *8 x(1),y(1)
c
c       copies x into y
c
        do 1400 i=1,n
          y(i)=x(i)
 1400 continue
c
        return
        end
c
c
        subroutine ichunkcopy(n,ix,iy)
        implicit real *8 (a-h,o-z)
        integer ix(*),iy(*)
c
c       copies ix into iy
c
        do 1400 i=1,n
          iy(i)=ix(i)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine chunksize(k, chunk, der, h, xymid, rad, cmass,
     1      crad, rltot)
        implicit real *8 (a-h,o-z)
        real *8 chunk(2,1),der(2,1),xs(1000),whts(1000),
     1      u(10000),v(10000),dsdt(1000),coefs(1000),
     2      ts(1000),coefs2(1000),xpts(1000),ypts(1000),
     3      xymid(1),ws(1000),cmass(1)
c
c       takes a discretized chunk and returns its midpoint (in
c       arclength) and total arclength
c
c       input:
c
c         k - number of nodes on the chunk
c         chunk - a (2,k) array containing x,y nodes
c         der - a (2,k) array containing dxdt, dydt where t is *some*
c             underlying parameterization
c         h - weighting factor to account for the arbitrary
c             parameterization (see chunkfunc)
c
c       output:
c
c         xymid - x,y coordinates of the midpoint in arclength
c         rad - the minimum radius of the ball centered at xymid that
c             could wholly contain the chunk
c         cmass - x,y coordinates of the center of mass of the
c             chunk
c         crad - the minimum radius of the ball centered at cmass that
c             could wholly contain the chunk
c         rltot - total arclength of the chunk
c
c       NOTE: this routine will break is k is larger than 100 !!!!
c
c
        done=1
c
c       first calculate the arclength - we are assuming that the chunk
c       has been computed before and resolved to "some" precision
c       using k nodes
c
        itype=2
        call legeexps(itype,k,xs,u,v,whts)
c
        rltot=0
        do 1000 i=1,k
          dsdt(i)=sqrt(der(1,i)**2+der(2,i)**2)*h
          rltot=rltot+dsdt(i)*whts(i)        
 1000 continue

c
c       use Newton to find t such that length(-1,t) = length(t,1) = rl/2
c
        call chunkmatvec(k,k,u,dsdt,coefs)
c
        t1=0
        rlhalf=rltot/2
        thresh=1.0d-10
        ifdone=0
c
        do 2000 ijk=1,1000

        do 1200 j=1,k
        ts(j)=-1+(t1+1)*(xs(j)+1)/2
        ws(j)=whts(j)/2*(t1+1)
 1200 continue

c
        rl1=0
        do 1400 i=1,k
        call legeexev(ts(i),val,coefs,k-1)
        rl1=rl1+val*ws(i)
 1400 continue
c
        call legeexev(t1,val,coefs,k-1)
        t2=t1-(rl1-rlhalf)/val

        err=rl1-rlhalf
cccc        call prin2('err in newton=*',err,1)
        if (abs(err) .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 3) goto 2100
c
        t1=t2
c
 2000 continue
 2100 continue
c
        if (ijk .ge. 1000) then
          call prinf('newton did not converge, chunksize!*',done,0)
          call prin2('t2=*',t2,1)
          call prin2('err=*',err,1)
          stop
        endif

c
c       chunk(t2) is now the midpoint - evaluate the expansions there
c
        do 3200 i=1,k
          xpts(i)=chunk(1,i)
          ypts(i)=chunk(2,i)
 3200 continue
c
        call chunkmatvec(k,k,u,xpts,coefs)
        call chunkmatvec(k,k,u,ypts,coefs2)
c
        call legeexev(t2,xmid,coefs,k-1)
        call legeexev(t2,ymid,coefs2,k-1)
c
        xymid(1)=xmid
        xymid(2)=ymid

c
c       compute the center of mass
c
        xc=0
        yc=0
        do 3600 i=1,k
          xc=xc+chunk(1,i)*dsdt(i)*whts(i)
          yc=yc+chunk(2,i)*dsdt(i)*whts(i)
 3600 continue
c
        xc=xc/rltot
        yc=yc/rltot
        cmass(1)=xc
        cmass(2)=yc

c
c       calculate various ball radii
c
        t0=-1
        t1=1
        call legeexev(t0,x0,coefs,k-1)
        call legeexev(t0,y0,coefs2,k-1)
        call legeexev(t1,x1,coefs,k-1)
        call legeexev(t1,y1,coefs2,k-1)
c
        r0=sqrt((x0-xmid)**2+(y0-ymid)**2)
        r1=sqrt((x1-xmid)**2+(y1-ymid)**2)
        rad=r0
        if (r1 .gt. r0) rad=r1
c
        cr0=sqrt((x0-xc)**2+(y0-yc)**2)
        cr1=sqrt((x1-xc)**2+(y1-yc)**2)
        crad=cr0
        if (cr1 .gt. cr0) crad=cr1
c
        do 4200 i=1,k
        r=sqrt((chunk(1,i)-xmid)**2+(chunk(2,i)-ymid)**2)
        cr=sqrt((chunk(1,i)-xc)**2+(chunk(2,i)-yc)**2)
        if (r .gt. rad) rad=r
        if (cr .gt. crad) crad=cr
 4200 continue

c
c       add a little buffer, probably not necessary...
c
        rad=rad*1.05d0
        crad=crad*1.05d0
c
        return
        end
c
c
c
c
        subroutine chunkcmf(k,chunk,der,h,xs,u,v,whts,xymid,rad)
        implicit real *8 (a-h,o-z)
        integer k
        real *8 chunk(2,k),der(2,k),xs(k),whts(k),
     1      u(k,k),v(k,k),dsdt(k),coefs(k),
     2      ts(k),coefs2(k),xpts(k),ypts(k),
     3      xymid(2),ws(k)
c
c       takes a discretized chunk and returns its midpoint (in
c       arclength) and total arclength
c
c       input:
c
c         k - number of nodes on the chunk
c         chunk - a (2,k) array containing x,y nodes
c         der - a (2,k) array containing dxdt, dydt where t is *some*
c             underlying parameterization
c         h - weighting factor to account for the arbitrary
c             parameterization (see chunkfunc)
c
c         xs - the order k gaussian nodes
c  
c         u - the k*k matrix converting the values of a polynomial
c             of order k-1 at k legendre nodes into the
c             coefficients of its legendre expansion
c
c         v - the k*k matrix converting the coefficinets
c             of a k-term legndre expansion into its values at 
c             the k legendre nodes
c
c         whts - the quadrature weights
c
c       output:
c
c         xymid - x,y coordinates of the midpoint in arclength
c         rad - the minimum radius of the ball centered at xymid that
c             could wholly contain the chunk
c
c
c
        done=1
c
c       first calculate the arclength - we are assuming that the chunk
c       has been computed before and resolved to "some" precision
c       using k nodes
c
        rltot=0
        do 1000 i=1,k
          dsdt(i)=sqrt(der(1,i)**2+der(2,i)**2)*h
          rltot=rltot+dsdt(i)*whts(i)        
 1000 continue

c
c       use Newton to find t such that length(-1,t) = length(t,1) = rl/2
c
        call chunkmatvec(k,k,u,dsdt,coefs)
c
        t1=0
        rlhalf=rltot/2
        thresh=1.0d-10
        ifdone=0
c
        do 2000 ijk=1,1000

        do 1200 j=1,k
        ts(j)=-1+(t1+1)*(xs(j)+1)/2
        ws(j)=whts(j)/2*(t1+1)
 1200 continue

c
        rl1=0
        do 1400 i=1,k
        call legeexev(ts(i),val,coefs,k-1)
        rl1=rl1+val*ws(i)
 1400 continue
c
        call legeexev(t1,val,coefs,k-1)
        t2=t1-(rl1-rlhalf)/val

        err=rl1-rlhalf
cccc        call prin2('err in newton=*',err,1)
        if (abs(err) .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 3) goto 2100
c
        t1=t2
c
 2000 continue
 2100 continue
c
        if (ijk .ge. 1000) then
          call prinf('newton did not converge, chunksize!*',done,0)
          call prin2('t2=*',t2,1)
          call prin2('err=*',err,1)
          stop
        endif

c
c       chunk(t2) is now the midpoint - evaluate the expansions there
c
        do 3200 i=1,k
          xpts(i)=chunk(1,i)
          ypts(i)=chunk(2,i)
 3200 continue
c
        call chunkmatvec(k,k,u,xpts,coefs)
        call chunkmatvec(k,k,u,ypts,coefs2)
c
        call legeexev(t2,xmid,coefs,k-1)
        call legeexev(t2,ymid,coefs2,k-1)
c
        xymid(1)=xmid
        xymid(2)=ymid

c       calculate various ball radii
c
        t0=-1
        t1=1
        call legeexev(t0,x0,coefs,k-1)
        call legeexev(t0,y0,coefs2,k-1)
        call legeexev(t1,x1,coefs,k-1)
        call legeexev(t1,y1,coefs2,k-1)
c
        r0=sqrt((x0-xmid)**2+(y0-ymid)**2)
        r1=sqrt((x1-xmid)**2+(y1-ymid)**2)
        rad=r0
        if (r1 .gt. r0) rad=r1
c
        do 4200 i=1,k
        r=sqrt((chunk(1,i)-xmid)**2+(chunk(2,i)-ymid)**2)
        if (r .gt. rad) rad=r
 4200 continue

c
c       add a little buffer, probably not necessary...
c
        rad=rad*1.05d0
c
        return
        end
c
        subroutine chunkpoly(eps,nverts,verts,ifclosed,ifdelete,
     1      nover,k,nch,chunks,adjs,ders,ders2,itypes,hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1),itypes(1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),
     1      hs(1),verts(2,1),
     2      xy0(10),xy1(10),ts(1000),whts(1000)
c
        double precision, allocatable :: ab(:,:), rls(:)
c
c       this routine returns chunks on a polygonal domain. note
c       that eps and nover behave differently in this routine than
c       they do in chunkfunc - as they should!
c
c       input:
c
c         eps - this is the desired size of the chunk near the corner
c         nverts - number of vertices (should be one more than the 
c             number of edges)
c         verts - number of vertices (if the curve is closed, then
c             verts(1) = verts(nverts) )
c         ifclosed - 0 if the curve is not closed, 1 if it is
c         ifdelete - set to 1 if you want to delete the chunk of
c             size eps next to the corner. setting ifdelete=1
c             will NOT delete an entire edge.
c         nover - this is the number of times to oversample the
c             polygon BEFORE the adaptive chunking begins, this is
c             in contrast to the post-adaptive chunking oversampling
c             performed in chunkfunc
c         k - number of nodes to put on each chunk
c
c       output:
c
c         nch - number of chunks created
c         chunks - x,y coordinates of each chunk
c         adjs - adjaceny information, note that ifdelete=1 then the
c             adjacent information for chunks near the corner will be -1
c         ders - values of dxdt, dydt for some parameterization t
c         ders2 - will be equal to 0 in this case
c         itypes - type of the chunk
c                      0 - middle segment, not touching a corner
c                      1 - the left end of the chunk is a corner
c                      2 - the right end of the chunk is a corner
c                      3 - both ends of the chunk are corners
c         hs - weighting parameters for each chunk to account for
c             arbitrary underlying parameterization (since things are
c             not necessarily in arclength)
c
c
        done=1
        pi=4*atan(done)
c
        maxpts = 10000000
        maxchunks = maxpts/k
c
        allocate(ab(6,maxchunks))
        allocate(rls(maxchunks))
c
c       start by oversampling each edge nover times
c
        nch=0
        nedges=nverts-1
        nover2=nover
        if (nover .le. 1) nover2=1
c
        do 3000 ijk=1,nedges
c
        xy0(1)=verts(1,ijk)
        xy0(2)=verts(2,ijk)
c
        xy1(1)=verts(1,ijk+1)
        xy1(2)=verts(2,ijk+1)
c
        dx=xy1(1)-xy0(1)
        dy=xy1(2)-xy0(2)
        rlen=sqrt(dx**2+dy**2)
        rlen5=rlen/nover2
c
        vecx=dx/rlen
        vecy=dy/rlen
c
        do 1600 i=1,nover2
c
        nch=nch+1
        ab(1,nch)=xy0(1)+(i-1)*rlen5*vecx
        ab(2,nch)=xy0(2)+(i-1)*rlen5*vecy
        ab(3,nch)=xy0(1)+i*rlen5*vecx
        ab(4,nch)=xy0(2)+i*rlen5*vecy
        rls(nch)=rlen5
c
        itypes(nch)=0
        if (i .eq. 1) itypes(nch)=1
        if (i .eq. nover2) itypes(nch)=2
        if (nover2 .eq. 1) itypes(nch)=3
c
 1600 continue
 3000 continue

c
c       compute adjacency information
c
        if (nch .eq. 1) then
            adjs(1,1)=-1
            adjs(2,1)=-1
            goto 4500
        endif
c
        do 3400 i=1,nch
c
        if ((i .eq. 1) .and. (ifclosed .ne. 1)) then
            adjs(1,i)=-1
            adjs(2,i)=i+1
            goto 3400
        endif
c
        if ((i .eq. 1) .and. (ifclosed .eq. 1)) then
            adjs(1,i)=nch
            adjs(2,i)=i+1
            goto 3400
        endif
c
        if ((i .eq. nch) .and. (ifclosed .ne. 1)) then
            adjs(1,i)=i-1
            adjs(2,i)=-1
            goto 3400
        endif
c
        if ((i .eq. nch) .and. (ifclosed .eq. 1)) then
            adjs(1,i)=i-1
            adjs(2,i)=1
            goto 3400
        endif
c
        adjs(1,i)=i-1
        adjs(2,i)=i+1
c
 3400 continue
 4500 continue

c
c       now run through the chunks that are not middle segments
c       and divide if they are larger than eps, update chunk info
c
cccc        ifdone=1
        maxiter=200
        do 7000 ijk=1,maxiter
c
        nchold=nch
        ifdone=1
c
        do 6400 i=1,nchold
c
        if (itypes(i) .eq. 0) goto 6400
c
c       check on whether to delete...
c
        ifd=1
        if (ifdelete .ne. 1) ifd=0
        if (itypes(i) .eq. 0) ifd=0
        if (itypes(i) .eq. 3) ifd=0
        if (rls(i)/2 .ge. eps) ifd=0
        if (ifd .eq. 0) goto 5100
c
c       if here, the chunk should be split and once it is split the
c       end should be deleted - no need to set ifdone=0 since this
c       chunk will be done
c
        if (itypes(i) .eq. 1) then
c
            itypes(i)=0
            xy0(1)=ab(1,i)
            xy0(2)=ab(2,i)
            xy1(1)=ab(3,i)
            xy1(2)=ab(4,i)
c
            rlsnew=rls(i)/2
            rls(i)=rlsnew
c
            ab(1,i)=(xy0(1)+xy1(1))/2
            ab(2,i)=(xy0(2)+xy1(2))/2
            i1=adjs(1,i)
            i2=adjs(2,i)
c
            if (i1 .gt. 0) adjs(2,i1)=-1
            adjs(1,i)=-1
c
        endif
c
        if (itypes(i) .eq. 2) then
c
            itypes(i)=0
            xy0(1)=ab(1,i)
            xy0(2)=ab(2,i)
            xy1(1)=ab(3,i)
            xy1(2)=ab(4,i)
c
            rlsnew=rls(i)/2
            rls(i)=rlsnew
c
            ab(3,i)=(xy0(1)+xy1(1))/2
            ab(4,i)=(xy0(2)+xy1(2))/2
            i1=adjs(1,i)
            i2=adjs(2,i)
c
            if (i2 .gt. 0) adjs(1,i2)=-1
            adjs(2,i)=-1
c
        endif
c
        goto 6400
c
 5100 continue
c
        if (rls(i) .lt. eps) goto 6400

c
c       if here, divide! the chunk is large enough so as to NOT be
c       deleted once split.
c
        ifdone=0
        if (itypes(i) .eq. 1) then
            itypes(i)=1
            itypes(nch+1)=0
        endif
c
        if (itypes(i) .eq. 2) then
            itypes(i)=0
            itypes(nch+1)=2
        endif
c
        if (itypes(i) .eq. 3) then
            itypes(i)=1
            itypes(nch+1)=2
        endif
c
        xy0(1)=ab(1,i)
        xy0(2)=ab(2,i)
        xy1(1)=ab(3,i)
        xy1(2)=ab(4,i)
c
        rlsnew=rls(i)/2
        rls(i)=rlsnew
        rls(nch+1)=rlsnew
c
        ab(3,i)=(xy0(1)+xy1(1))/2
        ab(4,i)=(xy0(2)+xy1(2))/2
c
        ab(1,nch+1)=(xy0(1)+xy1(1))/2
        ab(2,nch+1)=(xy0(2)+xy1(2))/2
        ab(3,nch+1)=xy1(1)
        ab(4,nch+1)=xy1(2)
c
        i1=adjs(1,i)
        i2=adjs(2,i)
c
        adjs(2,i)=nch+1
        adjs(1,i2)=nch+1
        adjs(1,nch+1)=i
        adjs(2,nch+1)=i2
c
        nch=nch+1
c

 6400 continue
c
        if (ifdone .eq. 1) goto 7100
c
 7000 continue
 7100 continue

        if (ijk .ge. maxiter) then
            call prinf('in chunkpoly, ijk=*',ijk,1)
            call prinf('did not finish! maxiter=*',maxiter,1)
            stop
        endif

c
c       now go through and fill in the legendre points on each chunk,
c       along with derivatives and scale factors hs
c
        ifwhts=1
        call legewhts(k,ts,whts,ifwhts)
c
        do 8400 i=1,nch
c
        slopex=ab(3,i)-ab(1,i)
        slopey=ab(4,i)-ab(2,i)

        do 7600 j=1,k       
c
        chunks(1,j,i)=ab(1,i)+slopex*(ts(j)+1)/2
        chunks(2,j,i)=ab(2,i)+slopey*(ts(j)+1)/2
c
        ders(1,j,i)=slopex/rls(i)
        ders(2,j,i)=slopey/rls(i)
c
        ders2(1,j,i)=0
        ders2(2,j,i)=0
c

 7600 continue
c
        hs(i)=rls(i)/2
c
 8400 continue
c
        return
        end
c
c
c
c
c
        subroutine chunksplit_pack(ich, wgeo, lused)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1)
c
        integer, allocatable :: adjs(:,:)
        double precision, allocatable :: chunks(:,:,:), ders(:,:,:)
        double precision, allocatable :: ders2(:,:,:), hs(:)
c
c       split chunk number ich which is stored in the work array
c       wgeo - wgeo AND lused are updated
c
        call chunkunpack1(wgeo, k, nch, ichunks, iadjs, iders,
     1    iders2, ihs)
c
        allocate(chunks(2,k,nch+10))
        allocate(adjs(2,nch+10))
        allocate(ders(2,k,nch+10))
        allocate(ders2(2,k,nch+10))
        allocate(hs(nch+10))
c
        call chunkunpack(wgeo, k, nch, chunks, adjs, ders, ders2, hs)
c
        call chunksplit1(ich, k, nch, chunks, adjs, ders,
     1    ders2, hs)
c
        call chunkpack(k, nch, chunks, adjs, ders, ders2,
     1    hs, wgeo, lused)
c
        return
        end
c
c
c
c
c
        subroutine chunksplit1(ich,k,nch,chunks,adjs,ders,
     1      ders2,hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1),
     1      coefs(1000),xs(10000),whts(10000),u(10000),
     2      v(10000),ts(1000),ws(1000),ts1(1000),ts2(1000),
     3      xymid(10),cmass(10),dsdt(1000),x(1000),y(1000),
     4      xp(1000),yp(1000),xpp(1000),ypp(1000),
     5      cx(1000),cy(1000),cxp(1000),cyp(1000),cxpp(1000),
     6      cypp(1000)
c
c       this routine takes the list of all chunks and splits one in
c       half with respect to arclength.
c
c       input:
c
c         ich - the chunk number to split
c         k - number of nodes per chunk
c         nch - total number of chunks
c         chunks - chunks array created by, for example, chunkfunc
c         adjs - adjacency information
c         ders - ders array created by, for example, chunkfunc
c         ders2 - ders2 array created by, for example, chunkfunc
c         hs - scaling factors created by, for example, chunkfunc
c
c       output:
c
c         nch - will be incremented by 1
c         chunks - updated with new chunk info
c         adjs - updated with new chunk info
c         ders - updated with new chunk info
c         ders2 - updated with new chunk info
c         hs - updated with new chunk info
c
c
        done=1
        itype=2
        call legeexps(itype,k,xs,u,v,whts)

c
c       first construct dsdt
c
        rltot=0
        do 1000 i=1,k
        dsdt(i)=sqrt(ders(1,i,ich)**2+ders(2,i,ich)**2)*hs(ich)
        rltot=rltot+dsdt(i)*whts(i)        
 1000 continue
c
        call chunkmatvec(k,k,u,dsdt,coefs)
cccc        call prin2('coefs for dsdt=*',coefs,k)
c
        t1=0
        rlhalf=rltot/2
        thresh=1.0d-10
        ifdone=0

c
c       use Newton to find t such that length(-1,t) = length(t,1) = rl/2
c
        do 2000 ijk=1,1000

        do 1200 j=1,k
        ts(j)=-1+(t1+1)*(xs(j)+1)/2
        ws(j)=whts(j)/2*(t1+1)
 1200 continue

c
        rl1=0
        do 1400 i=1,k
        call legeexev(ts(i),val,coefs,k-1)
        rl1=rl1+val*ws(i)
 1400 continue
c
        call legeexev(t1,val,coefs,k-1)
        t2=t1-(rl1-rlhalf)/val

        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 3) goto 2100
c
        t1=t2
c
 2000 continue
 2100 continue
c
        if (ijk .ge. 1000) then
            call prinf('newton did not converge! chunksplit1*',t1,0)
            call prinf('number of iterations=*',ijk,1)
            stop
        endif

c
c       now construct nodes on [-1,t2] and [t2,1]
c
        do 2600 i=1,k
        ts1(i)=-1+(t2+1)*(xs(i)+1)/2
        ts2(i)=t2+(1-t2)*(xs(i)+1)/2
 2600 continue

c
c       evaluate the new values of chunks, ders, ders2 and 
c       update nch, adjs, hs
c
        call chunk2to1(k,chunks(1,1,ich),x,y)
        call chunk2to1(k,ders(1,1,ich),xp,yp)
        call chunk2to1(k,ders2(1,1,ich),xpp,ypp)
c
        call chunkmatvec(k,k,u,x,cx)
        call chunkmatvec(k,k,u,y,cy)
        call chunkmatvec(k,k,u,xp,cxp)
        call chunkmatvec(k,k,u,yp,cyp)
        call chunkmatvec(k,k,u,xpp,cxpp)
        call chunkmatvec(k,k,u,ypp,cypp)
c
        i1=adjs(1,ich)
        i2=adjs(2,ich)
c
        do 3200 i=1,k
        call legeexev(ts1(i),chunks(1,i,ich),cx,k-1)
        call legeexev(ts1(i),chunks(2,i,ich),cy,k-1)
        call legeexev(ts1(i),ders(1,i,ich),cxp,k-1)
        call legeexev(ts1(i),ders(2,i,ich),cyp,k-1)
        call legeexev(ts1(i),ders2(1,i,ich),cxpp,k-1)
        call legeexev(ts1(i),ders2(2,i,ich),cypp,k-1)
 3200 continue        
c
        do 3600 i=1,k
        call legeexev(ts2(i),chunks(1,i,nch+1),cx,k-1)
        call legeexev(ts2(i),chunks(2,i,nch+1),cy,k-1)
        call legeexev(ts2(i),ders(1,i,nch+1),cxp,k-1)
        call legeexev(ts2(i),ders(2,i,nch+1),cyp,k-1)
        call legeexev(ts2(i),ders2(1,i,nch+1),cxpp,k-1)
        call legeexev(ts2(i),ders2(2,i,nch+1),cypp,k-1)
 3600 continue        
c
        hsold=hs(ich)
        hs(ich)=hsold*(t2+1)/2
        hs(nch+1)=hsold*(1-t2)/2
c
        adjs(1,ich)=i1
        adjs(2,ich)=nch+1

c
        adjs(1,nch+1)=ich
        adjs(2,nch+1)=i2
        nch=nch+1
c
        adjs(1,i2)=nch
c
        return
        end
    
        subroutine chunksplit1ab(ich,k,nch,chunks,adjs,ders,
     1      ders2,hs,ab)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1),ab(2,1),
     1      coefs(1000),xs(10000),whts(10000),u(10000),
     2      v(10000),ts(1000),ws(1000),ts1(1000),ts2(1000),
     3      xymid(10),cmass(10),dsdt(1000),x(1000),y(1000),
     4      xp(1000),yp(1000),xpp(1000),ypp(1000),
     5      cx(1000),cy(1000),cxp(1000),cyp(1000),cxpp(1000),
     6      cypp(1000)
c
c       this routine takes the list of all chunks and splits one in
c       half with respect to arclength.
c
c       input:
c
c         ich - the chunk number to split
c         k - number of nodes per chunk
c         nch - total number of chunks
c         chunks - chunks array created by, for example, chunkfunc
c         adjs - adjacency information
c         ders - ders array created by, for example, chunkfunc
c         ders2 - ders2 array created by, for example, chunkfunc
c         hs - scaling factors created by, for example, chunkfunc
c
c       output:
c
c         nch - will be incremented by 1
c         chunks - updated with new chunk info
c         adjs - updated with new chunk info
c         ders - updated with new chunk info
c         ders2 - updated with new chunk info
c         hs - updated with new chunk info
c
c
        done=1
        itype=2
        call legeexps(itype,k,xs,u,v,whts)

c
c       first construct dsdt
c
        rltot=0
        do 1000 i=1,k
        dsdt(i)=sqrt(ders(1,i,ich)**2+ders(2,i,ich)**2)*hs(ich)
        rltot=rltot+dsdt(i)*whts(i)        
 1000 continue
c
        call chunkmatvec(k,k,u,dsdt,coefs)
cccc        call prin2('coefs for dsdt=*',coefs,k)
c
        t1=0
        rlhalf=rltot/2
        thresh=1.0d-10
        ifdone=0

c
c       use Newton to find t such that length(-1,t) = length(t,1) = rl/2
c
        do 2000 ijk=1,1000

        do 1200 j=1,k
        ts(j)=-1+(t1+1)*(xs(j)+1)/2
        ws(j)=whts(j)/2*(t1+1)
 1200 continue

c
        rl1=0
        do 1400 i=1,k
        call legeexev(ts(i),val,coefs,k-1)
        rl1=rl1+val*ws(i)
 1400 continue
c
        call legeexev(t1,val,coefs,k-1)
        t2=t1-(rl1-rlhalf)/val

        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 3) goto 2100
c
        t1=t2
c
 2000 continue
 2100 continue
c
        if (ijk .ge. 1000) then
            call prinf('newton did not converge! chunksplit1*',t1,0)
            call prinf('number of iterations=*',ijk,1)
            stop
        endif

c
c       now construct nodes on [-1,t2] and [t2,1]
c
        do 2600 i=1,k
        ts1(i)=-1+(t2+1)*(xs(i)+1)/2
        ts2(i)=t2+(1-t2)*(xs(i)+1)/2
 2600 continue

c
c       evaluate the new values of chunks, ders, ders2 and 
c       update nch, adjs, hs
c
        call chunk2to1(k,chunks(1,1,ich),x,y)
        call chunk2to1(k,ders(1,1,ich),xp,yp)
        call chunk2to1(k,ders2(1,1,ich),xpp,ypp)
c
        call chunkmatvec(k,k,u,x,cx)
        call chunkmatvec(k,k,u,y,cy)
        call chunkmatvec(k,k,u,xp,cxp)
        call chunkmatvec(k,k,u,yp,cyp)
        call chunkmatvec(k,k,u,xpp,cxpp)
        call chunkmatvec(k,k,u,ypp,cypp)
c
        i1=adjs(1,ich)
        i2=adjs(2,ich)
c
        do 3200 i=1,k
        call legeexev(ts1(i),chunks(1,i,ich),cx,k-1)
        call legeexev(ts1(i),chunks(2,i,ich),cy,k-1)
        call legeexev(ts1(i),ders(1,i,ich),cxp,k-1)
        call legeexev(ts1(i),ders(2,i,ich),cyp,k-1)
        call legeexev(ts1(i),ders2(1,i,ich),cxpp,k-1)
        call legeexev(ts1(i),ders2(2,i,ich),cypp,k-1)
 3200 continue        
c
        do 3600 i=1,k
        call legeexev(ts2(i),chunks(1,i,nch+1),cx,k-1)
        call legeexev(ts2(i),chunks(2,i,nch+1),cy,k-1)
        call legeexev(ts2(i),ders(1,i,nch+1),cxp,k-1)
        call legeexev(ts2(i),ders(2,i,nch+1),cyp,k-1)
        call legeexev(ts2(i),ders2(1,i,nch+1),cxpp,k-1)
        call legeexev(ts2(i),ders2(2,i,nch+1),cypp,k-1)
 3600 continue        
c
        hsold=hs(ich)
        hs(ich)=hsold*(t2+1)/2
        hs(nch+1)=hsold*(1-t2)/2
c
        adjs(1,ich)=i1
        adjs(2,ich)=nch+1
        
        tt1 = ab(1,ich)
        tt2 = ab(2,ich)

        tt = tt1 + (t2 + 1)*(tt2-tt1)/2

        ab(1,ich) = tt1
        ab(2,ich) = tt

        ab(1,nch+1) = tt
        ab(2,nch+1) = tt2
c


c
        adjs(1,nch+1)=ich
        adjs(2,nch+1)=i2
        nch=nch+1
c
        adjs(1,i2)=nch
c
        return
        end
    
c        
c
        subroutine chunksplit2(ich,k,nch,chunks,adjs,ders,
     1      ders2,hs,ifsigma1,sigma1,ifsigma2,sigma2,ifsigma3,
     2      sigma3)
        implicit real *8 (a-h,o-z)
        integer ifsigma1,ifsigma2,ifsigma3
        integer adjs(2,1)
        real *8 chunks(2,k,1),ders(2,k,1),ders2(2,k,1),hs(1),
     1      coefs(1000),xs(10000),whts(10000),u(10000),
     2      v(10000),ts(1000),ws(1000),ts1(1000),ts2(1000),
     3      xymid(10),cmass(10),dsdt(1000),x(1000),y(1000),
     4      xp(1000),yp(1000),xpp(1000),ypp(1000),
     5      cx(1000),cy(1000),cxp(1000),cyp(1000),cxpp(1000),
     6      cypp(1000)

        complex *16 sigma1(k,1), sigma2(k,1), sigma3(k,1)
        real *8 sigma1r(k), sigma2r(k) ,sigma3r(k)
        real *8 sigma1c(k), sigma2c(k) ,sigma3c(k)
        real *8 s1r(k), s2r(k) ,s3r(k)
        real *8 s1c(k), s2c(k), s3c(k)
c
c       this routine takes the list of all chunks and splits one in
c       half with respect to arclength and resamples
c       the densities sigma1,sigma2 and sigma3 on the new
c       chunks
c
c       input:
c
c         ich - the chunk number to split
c         k - number of nodes per chunk
c         nch - total number of chunks
c         chunks - chunks array created by, for example, chunkfunc
c         adjs - adjacency information
c         ders - ders array created by, for example, chunkfunc
c         ders2 - ders2 array created by, for example, chunkfunc
c         hs - scaling factors created by, for example, chunkfunc
c         ifsigma1 - flag for resampling sigma1
c         sigma1 - density prescribed on the chunked domain
c         ifsigma2 - flag for resampling sigma2
c         sigma2 - density prescribed on the chunked domain
c         ifsigam3 - flag for resampling sigma3
c         sigma3 - density prescribed on the chunked domain
c         NOTE: It is assumed that all the densities
c         are of data type complex *16
c
c       output:
c
c         nch - will be incremented by 1
c         chunks - updated with new chunk info
c         adjs - updated with new chunk info
c         ders - updated with new chunk info
c         ders2 - updated with new chunk info
c         hs - updated with new chunk info
c         sigma1 - resampled density on the new chunks
c         sigma2 - resampled density on the new chunks
c         sigma3 - resampled density on the new chunks
c 
c
c
c
        done=1
        itype=2
        call legeexps(itype,k,xs,u,v,whts)

c
c       first construct dsdt
c
        rltot=0
        do 1000 i=1,k
        dsdt(i)=sqrt(ders(1,i,ich)**2+ders(2,i,ich)**2)*hs(ich)
        rltot=rltot+dsdt(i)*whts(i)        
 1000 continue
c
        call chunkmatvec(k,k,u,dsdt,coefs)
cccc        call prin2('coefs for dsdt=*',coefs,k)
c
        t1=0
        rlhalf=rltot/2
        thresh=1.0d-10
        ifdone=0

c
c       use Newton to find t such that length(-1,t) = length(t,1) = rl/2
c
        do 2000 ijk=1,1000

        do 1200 j=1,k
        ts(j)=-1+(t1+1)*(xs(j)+1)/2
        ws(j)=whts(j)/2*(t1+1)
 1200 continue

c
        rl1=0
        do 1400 i=1,k
        call legeexev(ts(i),val,coefs,k-1)
        rl1=rl1+val*ws(i)
 1400 continue
c
        call legeexev(t1,val,coefs,k-1)
        t2=t1-(rl1-rlhalf)/val

        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 3) goto 2100
c
        t1=t2
c
 2000 continue
 2100 continue
c
        if (ijk .ge. 1000) then
            call prinf('newton did not converge! chunksplit1*',t1,0)
            call prinf('number of iterations=*',ijk,1)
            stop
        endif

c
c       now construct nodes on [-1,t2] and [t2,1]
c
        do 2600 i=1,k
        ts1(i)=-1+(t2+1)*(xs(i)+1)/2
        ts2(i)=t2+(1-t2)*(xs(i)+1)/2
 2600 continue

c
c       evaluate the new values of chunks, ders, ders2 and 
c       update nch, adjs, hs, sigma1,sigma2,sigma3
c
        call chunk2to1(k,chunks(1,1,ich),x,y)
        call chunk2to1(k,ders(1,1,ich),xp,yp)
        call chunk2to1(k,ders2(1,1,ich),xpp,ypp)

        do i=1,k
           if(ifsigma1.eq.1) then
              sigma1r(i) = dreal(sigma1(i,ich))
              sigma1c(i) = dimag(sigma1(i,ich))
           endif

           if(ifsigma2.eq.1) then
              sigma2r(i) = dreal(sigma2(i,ich))
              sigma2c(i) = dimag(sigma2(i,ich))
           endif

           if(ifsigma3.eq.1) then
              sigma3r(i) = dreal(sigma3(i,ich))
              sigma3c(i) = dimag(sigma3(i,ich))
           endif

        enddo
c
        call chunkmatvec(k,k,u,x,cx)
        call chunkmatvec(k,k,u,y,cy)
        call chunkmatvec(k,k,u,xp,cxp)
        call chunkmatvec(k,k,u,yp,cyp)
        call chunkmatvec(k,k,u,xpp,cxpp)
        call chunkmatvec(k,k,u,ypp,cypp)

        if(ifsigma1.eq.1) then
           call chunkmatvec(k,k,u,sigma1r,s1r)
           call chunkmatvec(k,k,u,sigma1c,s1c)
        endif

        if(ifsigma2.eq.1) then
           call chunkmatvec(k,k,u,sigma2r,s2r)
           call chunkmatvec(k,k,u,sigma2c,s2c)
        endif
        
        if(ifsigma3.eq.1) then
           call chunkmatvec(k,k,u,sigma3r,s3r)
           call chunkmatvec(k,k,u,sigma3c,s3c)
        endif
        
        i1=adjs(1,ich)
        i2=adjs(2,ich)
c
        do 3200 i=1,k
        call legeexev(ts1(i),chunks(1,i,ich),cx,k-1)
        call legeexev(ts1(i),chunks(2,i,ich),cy,k-1)
        call legeexev(ts1(i),ders(1,i,ich),cxp,k-1)
        call legeexev(ts1(i),ders(2,i,ich),cyp,k-1)
        call legeexev(ts1(i),ders2(1,i,ich),cxpp,k-1)
        call legeexev(ts1(i),ders2(2,i,ich),cypp,k-1)

        if(ifsigma1.eq.1) then
           call legeexev(ts1(i),sigma1r(i),s1r,k-1)
           call legeexev(ts1(i),sigma1c(i),s1c,k-1)
           sigma1(i,ich) = dcmplx(sigma1r(i),sigma1c(i))
        endif

        if(ifsigma2.eq.1) then
           call legeexev(ts1(i),sigma2r(i),s2r,k-1)
           call legeexev(ts1(i),sigma2c(i),s2c,k-1)
           sigma2(i,ich) = dcmplx(sigma2r(i),sigma2c(i))
        endif

        if(ifsigma3.eq.1) then
           call legeexev(ts1(i),sigma3r(i),s3r,k-1)
           call legeexev(ts1(i),sigma3c(i),s3c,k-1)
           sigma3(i,ich) = dcmplx(sigma3r(i),sigma3c(i))
        endif

 3200 continue        
c
        do 3600 i=1,k
        call legeexev(ts2(i),chunks(1,i,nch+1),cx,k-1)
        call legeexev(ts2(i),chunks(2,i,nch+1),cy,k-1)
        call legeexev(ts2(i),ders(1,i,nch+1),cxp,k-1)
        call legeexev(ts2(i),ders(2,i,nch+1),cyp,k-1)
        call legeexev(ts2(i),ders2(1,i,nch+1),cxpp,k-1)
        call legeexev(ts2(i),ders2(2,i,nch+1),cypp,k-1)

        if(ifsigma1.eq.1) then
           call legeexev(ts2(i),sigma1r(i),s1r,k-1)
           call legeexev(ts2(i),sigma1c(i),s1c,k-1)
           sigma1(i,nch+1) = dcmplx(sigma1r(i),sigma1c(i))
        endif

        if(ifsigma2.eq.1) then
           call legeexev(ts2(i),sigma2r(i),s2r,k-1)
           call legeexev(ts2(i),sigma2c(i),s2c,k-1)
           sigma2(i,nch+1) = dcmplx(sigma2r(i),sigma2c(i))
        endif

        if(ifsigma3.eq.1) then
           call legeexev(ts2(i),sigma3r(i),s3r,k-1)
           call legeexev(ts2(i),sigma3c(i),s3c,k-1)
           sigma3(i,nch+1) = dcmplx(sigma3r(i),sigma3c(i))
        endif

 3600 continue        
c
        hsold=hs(ich)
        hs(ich)=hsold*(t2+1)/2
        hs(nch+1)=hsold*(1-t2)/2
c
        adjs(1,ich)=i1
        adjs(2,ich)=nch+1
c
        adjs(1,nch+1)=ich
        adjs(2,nch+1)=i2
        nch=nch+1
c
        adjs(1,i2)=nch
c
        return
        end
c
c
c
c
c
        subroutine chunkinterp2(k,chunk,kout,chunkout)
        implicit real *8 (a-h,o-z)
        real *8 chunk(2,k), chunkout(2,kout), xs(1000), ys(1000),
     1      xsout(1000),ysout(10000),coefs1(1000),coefs2(1000),
     2      ts(1000),whts(1000),ts2(1000),whts2(1000),
     3      u(100000),v(100000),u2(100000),v2(100000)
c
c       interpolate from k gaussian nodes to kout gaussian
c       nodes. kout should be less than 100!
c
c       input:
c           k - number of points on the original chunk
c           chunk - the chunk (or 2D function, however you want to
c               think about it)
c           kout - the number of points to upsample to
c
c       ouput:
c           chunkout - the up-sampled chunk (or 2D function) evaluated
c               at kout legendre nodes
c
c
        if (kout .gt. 100) then
            call prinf('kout if too big! kout = *', kout, 1)
            stop
        endif
c
        do 1800 i=1,k
        xs(i)=chunk(1,i)
        ys(i)=chunk(2,i)
 1800 continue
c
        itype=2
        call legeexps(itype,k,ts,u,v,whts)
        call legeexps(itype,kout,ts2,u2,v2,whts2)
c
        call chunkmatvec(k,k,u,xs,coefs1)
        call chunkmatvec(k,k,u,ys,coefs2)
c
        do 2200 i=k+1,kout
        coefs1(i)=0
        coefs2(i)=0
 2200 continue
c
        call chunkmatvec(kout,kout,v2,coefs1,xsout)
        call chunkmatvec(kout,kout,v2,coefs2,ysout)
c
        do 2600 i=1,kout
        chunkout(1,i)=xsout(i)
        chunkout(2,i)=ysout(i)
 2600 continue
c
        return
        end
c
c
c
c
        subroutine chunkinterp2f(k,chunk,uk,kout,chunkout,vkout)
        implicit real *8 (a-h,o-z)
        integer k,kout
        real *8 chunk(2,k), chunkout(2,kout), xs(k), ys(k),
     1      xsout(kout),ysout(kout),coefs1(kout),coefs2(kout),
     2      uk(k,k), vkout(kout,kout)
c
c       interpolate from k gaussian nodes to kout gaussian
c       nodes.
c
c       input:
c           k - number of points on the original chunk
c           chunk - the chunk (or 2D function, however you want to
c               think about it)
c           kout - the number of points to upsample to
c           uk - k*k matrix for converting the values of a polynomial
c                of order k-1 at k legendre nodes of its legendre
c                expansion
c           vkout - kout*kout matrix converting the coefficients of
c                   an koutterm-lengendre expansion into its values
c                   kout legendre nodes
c       ouput:
c           chunkout - the up-sampled chunk (or 2D function) evaluated
c               at kout legendre nodes
c
c
        do 1800 i=1,k
        xs(i)=chunk(1,i)
        ys(i)=chunk(2,i)
 1800 continue
c
        itype=2
c
        call chunkmatvec(k,k,uk,xs,coefs1)
        call chunkmatvec(k,k,uk,ys,coefs2)
c
        do 2200 i=k+1,kout
        coefs1(i)=0
        coefs2(i)=0
 2200 continue
c
        call chunkmatvec(kout,kout,vkout,coefs1,xsout)
        call chunkmatvec(kout,kout,vkout,coefs2,ysout)
c
        do 2600 i=1,kout
        chunkout(1,i)=xsout(i)
        chunkout(2,i)=ysout(i)
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine chunkfunc(eps, ifclosed, chsmall, ta, tb,
     1    funcurve, pars, nover, k, nch, chunks, adjs, ders,
     2    ders2, hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1),ifprocess(60000)
        double precision chunks(2,k,1),ders(2,k,1), ders2(2,k,1),
     1      xs(1000),ws(1000),xs2(1000),ws2(10000),u(10000),v(10000),
     2      u2(10000),v2(10000),ab(6,30000),hs(1),
     3      errs(1000),pars(1),
     4      errs0(1000)
        real *8 xcoefs(1000), ycoefs(1000), xpcoefs(1000)
        real *8 ypcoefs(1000), xpcoefs_out(1000)
        real *8 ypcoefs_out(1000)
        real *8 ch7(2,1000), der7(2,1000)
c
        double precision, allocatable :: fvals(:,:), coefs(:,:)
c
        external funcurve
c
c       using a user-defined subroutine funcurve, split up the
c       curve into chunks. the calling sequence of funcurve
c       should be
c
c           funcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
c
c       the routine will declare a chunk to be resolved when
c       chunks, ders, ders2 and dsdt = sqrt(ders(1,)**2 + ders(2,)**2)
c       have length 2k expansions with the last k coefficients
c       having relative-rmse less than eps.
c
c       NOTE: no effort whatsoever has been made to make this routine
c       efficient. it does what it does very deliberately, with some
c       redundancy.
c
c       NOTE 2: the routine returns chunks that are at most a factor
c       of two different in archlength than adjacent chunks
c
c       input:
c
c         eps - absolute precision to resolve the curve
c         ifclosed - switch telling the routine whether the curve is
c           closed or not, ifclosed=1 means a closed curve
c         chsmall - if the curve is NOT closed, then this is the maximum
c           size of the intervals on either end. this parameter
c           is ignored if ifclosed=0.
c         ta,tb - assume that curve is parameterized by t \in [ta,tb)
c         funcurve - see above
c         pars - parameter array to send to funcurve, see above
c         nover - post-process oversampling factor, nover <= 1 will result
c             in no change, nover=2 will split each chunk in half (with
c             respect to arclength!!)
c         k - number of legendre nodes per chunk
c
c       output:
c
c         nch - total number of chunks created
c         chunks - points on each chunk, dimensioned chunks(2,k,nch)
c         adjs - adjacency information, adjs(1,i) is the chunk before
c             chunk i and adjs(2,i) is the chunk after chunk i
c         ders - first derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         ders2 - second derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         hs - weight parameter to account for the arbitrary underlying
c             parameterization, used in arclength calculation and
c             subsequent integration
c
c
        done=1
        pi=4*atan(done)
c
        maxchunks=55000
        do i=1,maxchunks
          ifprocess(i)=0
        enddo

c
c       construct legendre nodes and weights, k and 2k of them, as well
c       as the interpolation/coefficients matrices
c
        itype=2
        call legeexps(itype,k,xs,u,v,ws)
        call legeexps(itype,2*k,xs2,u2,v2,ws2)

c
c       . . . start chunking
c
        ab(1,1)=ta
        ab(2,1)=tb
        nch=1
        ifdone=1
        adjs(1,1)=-1
        adjs(2,1)=-1
        nchnew=nch
c
        allocate(fvals(2*k,7))
        allocate(coefs(2*k,7))
c
        maxiter=10000
        do 5000 ijk = 1,maxiter

c
c       loop through all existing chunks, if resolved store, if not split
c
          ifdone=1
          do 4600 ich=1,nchnew
c
            if (ifprocess(ich) .eq. 1) goto 4600
            ifprocess(ich)=1
c
            a=ab(1,ich)
            b=ab(2,ich)
c
            do i=1,2*k
              t=a+(b-a)*(xs2(i)+1)/2
              call funcurve(t,pars,fvals(i,1),fvals(i,2),
     1          fvals(i,3),fvals(i,4),fvals(i,5),fvals(i,6))
              fvals(i,7)=sqrt(fvals(i,3)**2+fvals(i,4)**2)
            enddo
c
            do i=1,7
              call chunkmatvec(2*k,2*k,u2,fvals(1,i),coefs(1,i))
              errs0(i) = 0
              errs(i) = 0
            enddo
c
            do i=1,7
              do j=1,k
                errs0(i) = errs0(i) + coefs(j,i)**2
                errs(i) = errs(i) + coefs(k+j,i)**2
              enddo
            enddo
c
            rmsemax = -1.0d0
            do i=1,7
              errs(i) = sqrt(errs(i)/errs0(i)/k)
              if (errs(i) .gt. rmsemax) rmsemax=errs(i)
            enddo

c
c       . . . mark as processed and resolved if less than eps
c
            if (rmsemax .gt. eps) goto 2800
            goto 4600
            
c
 2800 continue

c
c       . . . if here, not resolved
c       divide - first update the adjacency list
c
            
            ifprocess(ich)=0
            ifdone=0
cccc            return
c
            if ((nch .eq. 1) .and. (ifclosed .gt. 0)) then
              adjs(1,nch)=2
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=1
            endif
c
            if ((nch .eq. 1) .and. (ifclosed .le. 0)) then
              adjs(1,nch)=-1
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=-1
            endif
c
            if (nch .gt. 1) then
              iold2=adjs(2,ich)
              adjs(2,ich)=nch+1
              if (iold2 .gt. 0) adjs(1,iold2)=nch+1
              adjs(1,nch+1)=ich
              adjs(2,nch+1)=iold2
            endif

c
c       now update the endpoints in ab
c
            ab(1,ich)=a
            ab(2,ich)=(a+b)/2
c
            nch=nch+1
            if (nch .gt. maxchunks) then
              call prinf('too many chunks in chunkfunc!*',done,0)
              call prinf('maxchunks=*',maxchunks,1)
              stop
            endif
c
            ab(1,nch)=(a+b)/2
            ab(2,nch)=b
c
 4500 continue
 4600 continue
c
          if ((ifdone .eq. 1) .and. (nchnew .eq. nch)) goto 5100
          nchnew=nch        
c
 5000 continue
 5100 continue

c
c       the curve should be resolved to precision eps now on
c       each interval ab(,i)
c       check the size of adjacent neighboring chunks - if off by a
c       factor of more than 2, split them as well. iterate until done.
c

        maxiter=1000
        do 9000 ijk=1,maxiter
c
          nchold=nch
          ifdone=1
          do 8600 i=1,nchold
c
            i1=adjs(1,i)
            i2=adjs(2,i)

c
c       calculate chunk lengths
c
            a=ab(1,i)
            b=ab(2,i)
            call chunklength(k,funcurve,pars,a,b,rlself)
c
            rl1=rlself
            rl2=rlself
c
            if (i1 .gt. 0) then
              a1=ab(1,i1)
              b1=ab(2,i1)
              call chunklength(k,funcurve,pars,a1,b1,rl1)
            endif
c
            if (i2 .gt. 0) then
              a2=ab(1,i2)
              b2=ab(2,i2)
              call chunklength(k,funcurve,pars,a2,b2,rl2)
            endif

c
c       only check if self is larger than either of adjacent blocks,
c       iterating a couple times will catch everything
c
            ifsplit=0
            sc = 2.05d0
            if (rlself .gt. sc*rl1) ifsplit=1
            if (rlself .gt. sc*rl2) ifsplit=1
            if (ifsplit .eq. 0) goto 8600
c
c       split chunk i now, and recalculate nodes, ders, etc
c
            ifdone=0
            a=ab(1,i)
            b=ab(2,i)
            ab2=(a+b)/2
c
            i1=adjs(1,i)
            i2=adjs(2,i)
c        
            adjs(1,i) = i1
            adjs(2,i) = nch+1
c
c       . . . first update nch+1
c
            adjs(1,nch+1) = i
            adjs(2,nch+1) = i2
c
c       . . . if there's an i2, update it
c
            if (i2 .gt. 0) then
              adjs(1,i2) = nch+1
            endif
c
            nch=nch+1
            if (nch .gt. maxchunks) then
                call prinf('too many chunks in chunkfunc!*',done,0)
                call prinf('maxchunks=*',maxchunks,1)
                stop
            endif
c
            ab(1,i)=a
            ab(2,i)=ab2
c
            ab(1,nch)=ab2
            ab(2,nch)=b
c
 8600 continue

c
          if (ifdone .ne. 0) goto 9100
c
 9000 continue
 9100 continue

c
c       go ahead and oversample by nover, updating
c       the adjacency information adjs along the way
c
        if (nover .le. 1) goto 6100

c
        do 6000 ijk=1,nover-1
c
        nchold=nch
        do 5600 i=1,nchold
c
        a=ab(1,i)
        b=ab(2,i)
c
c       find ab2 using newton such that 
c       len(a,ab2)=len(ab2,b)=half the chunk length
c
        call chunklength(k,funcurve,pars,a,b,rl)
c
        rlhalf=rl/2
        thresh=1.0d-8
        ifnewt=0
        ab0=(a+b)/2
c
        do 6600 iter=1,1000
c
        call chunklength(k,funcurve,pars,a,ab0,rl1)
        call funcurve(ab0,pars,x,y,dx,dy,ddx,ddy)
        dsdt=sqrt(dx**2+dy**2)
        ab1=ab0-(rl1-rlhalf)/dsdt
c
        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifnewt=ifnewt+1
c
        if (ifnewt .eq. 3) goto 6700
        ab0=ab1
c
 6600 continue
 6700 continue
c
        if (ifnewt .lt. 3) then
            call prin2('newton failed! interval not split.*',done,0)
            stop
        endif
c
        ab2=ab1
c
        i1=adjs(1,i)
        i2=adjs(2,i)
        adjs(2,i)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=i
        adjs(2,nch+1)=i2
c
        ab(1,i)=a
        ab(2,i)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
 5600 continue
c
 6000 continue
 6100 continue
c
        if (ifclosed .eq. 1) goto 7700
c
c       if the curve is open, check the dyadic refinement at the
c       ends, first find the end segments
c
        do 6400 i =1,nch
        if (adjs(1,i) .lt. 0) ileft=i
        if (adjs(2,i) .lt. 0) iright=i
 6400 continue
c
        a1=ab(1,ileft)
        b1=ab(2,ileft)
        a2=ab(1,iright)
        b2=ab(2,iright)
c
        call chunklength(k,funcurve,pars,a1,b1,rlleft)
        call chunklength(k,funcurve,pars,a2,b2,rlright)
c
c       . . . dyadically split the left segment
c
        do 7000 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('ifclosed = *', ifclosed, 1)
            call prinf('boom refining left end! ijk=*',ijk,1)
            call prinf('boom refining left end! rlleft=*',rlleft,1)
            stop
        endif
c
        if (rlleft .le. chsmall) goto 7100
        a=ab(1,ileft)
        b=ab(2,ileft)
        ab2=(a+b)/2
c
        i1=adjs(1,ileft)
        i2=adjs(2,ileft)
        adjs(2,ileft)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=ileft
        adjs(2,nch+1)=i2
c
        ab(1,ileft)=a
        ab(2,ileft)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
        call chunklength(k,funcurve,pars,a,ab2,rlleft)
c        
 7000 continue
 7100 continue

c
c       . . . dyadically split the right segment
c
        do 7400 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('boom refining right end! ijk=*',ijk,1)
            call prinf('boom right end! rlright=*',rlright,1)
            stop
        endif
c
        if (rlright .le. chsmall) goto 7500
        a=ab(1,iright)
        b=ab(2,iright)
        ab2=(a+b)/2
c
        i1=adjs(1,iright)
        i2=adjs(2,iright)
        adjs(2,iright)=nch+1
c
        adjs(1,nch+1)=iright
        adjs(2,nch+1)=i2
c
        ab(1,iright)=a
        ab(2,iright)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        iright=nch
        ab(1,nch)=ab2
        ab(2,nch)=b
c
        call chunklength(k,funcurve,pars,ab2,b,rlright)
c
 7400 continue
 7500 continue
c
 7700 continue

c
c       up to here, everything has been done in parameter space, [ta,tb]
c       . . . finally evaluate the k nodes on each chunk, along with 
c       derivatives and chunk lengths
c
        do i = 1, nch
c
          a=ab(1,i)
          b=ab(2,i)
          hs(i)=(b-a)/2
c
          do j = 1, k
            t=a+(b-a)*(xs(j)+1)/2
            call funcurve(t, pars, x, y, dx, dy, dx2, dy2)
            chunks(1,j,i) = x
            chunks(2,j,i) = y
            ders(1,j,i) = dx
            ders(2,j,i) = dy
            ders2(1,j,i) = dx2
            ders2(2,j,i) = dy2
          enddo
        enddo
c
        return
        end
c---------------------------------------------------------------------
        subroutine chunkfunc2(eps, ifclosed, chsmall, ta, tb,
     1    funcurve, pars, nover, k, nch, chunks, adjs, ders,
     2    ders2, hs, ab)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1),ifprocess(60000)
        double precision chunks(2,k,1),ders(2,k,1), ders2(2,k,1),
     1      xs(1000),ws(1000),xs2(1000),ws2(10000),u(10000),v(10000),
     2      u2(10000),v2(10000),ab(2,1),hs(1),
     3      errs(1000),pars(1),
     4      errs0(1000)
        real *8 xcoefs(1000), ycoefs(1000), xpcoefs(1000)
        real *8 ypcoefs(1000), xpcoefs_out(1000)
        real *8 ypcoefs_out(1000)
        real *8 ch7(2,1000), der7(2,1000)
c
        double precision, allocatable :: fvals(:,:), coefs(:,:)
c
        external funcurve
c
c       using a user-defined subroutine funcurve, split up the
c       curve into chunks. the calling sequence of funcurve
c       should be
c
c           funcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
c
c       the routine will declare a chunk to be resolved when
c       chunks, ders, ders2 and dsdt = sqrt(ders(1,)**2 + ders(2,)**2)
c       have length 2k expansions with the last k coefficients
c       having relative-rmse less than eps.
c
c       NOTE: no effort whatsoever has been made to make this routine
c       efficient. it does what it does very deliberately, with some
c       redundancy.
c
c       NOTE 2: the routine returns chunks that are at most a factor
c       of two different in archlength than adjacent chunks
c
c       input:
c
c         eps - absolute precision to resolve the curve
c         ifclosed - switch telling the routine whether the curve is
c           closed or not, ifclosed=1 means a closed curve
c         chsmall - if the curve is NOT closed, then this is the maximum
c           size of the intervals on either end. this parameter
c           is ignored if ifclosed=0.
c         ta,tb - assume that curve is parameterized by t \in [ta,tb)
c         funcurve - see above
c         pars - parameter array to send to funcurve, see above
c         nover - post-process oversampling factor, nover <= 1 will result
c             in no change, nover=2 will split each chunk in half (with
c             respect to arclength!!)
c         k - number of legendre nodes per chunk
c
c       output:
c
c         nch - total number of chunks created
c         chunks - points on each chunk, dimensioned chunks(2,k,nch)
c         adjs - adjacency information, adjs(1,i) is the chunk before
c             chunk i and adjs(2,i) is the chunk after chunk i
c         ders - first derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         ders2 - second derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         hs - weight parameter to account for the arbitrary underlying
c             parameterization, used in arclength calculation and
c             subsequent integration
c         ab - left and right end interval of chunked curve, 
c              dimensioned (2,nch)
c
c
        done=1
        pi=4*atan(done)
c
        maxchunks=55000
        do i=1,maxchunks
          ifprocess(i)=0
        enddo

c
c       construct legendre nodes and weights, k and 2k of them, as well
c       as the interpolation/coefficients matrices
c
        itype=2
        call legeexps(itype,k,xs,u,v,ws)
        call legeexps(itype,2*k,xs2,u2,v2,ws2)

c
c       . . . start chunking
c
        ab(1,1)=ta
        ab(2,1)=tb
        nch=1
        ifdone=1
        adjs(1,1)=-1
        adjs(2,1)=-1
        nchnew=nch
c
        allocate(fvals(2*k,7))
        allocate(coefs(2*k,7))
c
        maxiter=10000
        do 5000 ijk = 1,maxiter

c
c       loop through all existing chunks, if resolved store, if not split
c
          ifdone=1
          do 4600 ich=1,nchnew
c
            if (ifprocess(ich) .eq. 1) goto 4600
            ifprocess(ich)=1
c
            a=ab(1,ich)
            b=ab(2,ich)
c
            do i=1,2*k
              t=a+(b-a)*(xs2(i)+1)/2
              call funcurve(t,pars,fvals(i,1),fvals(i,2),
     1          fvals(i,3),fvals(i,4),fvals(i,5),fvals(i,6))
              fvals(i,7)=sqrt(fvals(i,3)**2+fvals(i,4)**2)
            enddo
c
            do i=1,7
              call chunkmatvec(2*k,2*k,u2,fvals(1,i),coefs(1,i))
              errs0(i) = 0
              errs(i) = 0
            enddo
c
            do i=1,7
              do j=1,k
                errs0(i) = errs0(i) + coefs(j,i)**2
                errs(i) = errs(i) + coefs(k+j,i)**2
              enddo
            enddo
c
            rmsemax = -1.0d0
            do i=1,7
              errs(i) = sqrt(errs(i)/errs0(i)/k)
              if (errs(i) .gt. rmsemax) rmsemax=errs(i)
            enddo

c
c       . . . mark as processed and resolved if less than eps
c
            if (rmsemax .gt. eps) goto 2800
            goto 4600
            
c
 2800 continue

c
c       . . . if here, not resolved
c       divide - first update the adjacency list
c
            
            ifprocess(ich)=0
            ifdone=0
cccc            return
c
            if ((nch .eq. 1) .and. (ifclosed .gt. 0)) then
              adjs(1,nch)=2
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=1
            endif
c
            if ((nch .eq. 1) .and. (ifclosed .le. 0)) then
              adjs(1,nch)=-1
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=-1
            endif
c
            if (nch .gt. 1) then
              iold2=adjs(2,ich)
              adjs(2,ich)=nch+1
              if (iold2 .gt. 0) adjs(1,iold2)=nch+1
              adjs(1,nch+1)=ich
              adjs(2,nch+1)=iold2
            endif

c
c       now update the endpoints in ab
c
            ab(1,ich)=a
            ab(2,ich)=(a+b)/2
c
            nch=nch+1
            if (nch .gt. maxchunks) then
              call prinf('too many chunks in chunkfunc!*',done,0)
              call prinf('maxchunks=*',maxchunks,1)
              stop
            endif
c
            ab(1,nch)=(a+b)/2
            ab(2,nch)=b
c
 4500 continue
 4600 continue
c
          if ((ifdone .eq. 1) .and. (nchnew .eq. nch)) goto 5100
          nchnew=nch        
c
 5000 continue
 5100 continue

c
c       the curve should be resolved to precision eps now on
c       each interval ab(,i)
c       check the size of adjacent neighboring chunks - if off by a
c       factor of more than 2, split them as well. iterate until done.
c

        maxiter=1000
        do 9000 ijk=1,maxiter
c
          nchold=nch
          ifdone=1
          do 8600 i=1,nchold
c
            i1=adjs(1,i)
            i2=adjs(2,i)

c
c       calculate chunk lengths
c
            a=ab(1,i)
            b=ab(2,i)
            call chunklength(k,funcurve,pars,a,b,rlself)
c
            rl1=rlself
            rl2=rlself
c
            if (i1 .gt. 0) then
              a1=ab(1,i1)
              b1=ab(2,i1)
              call chunklength(k,funcurve,pars,a1,b1,rl1)
            endif
c
            if (i2 .gt. 0) then
              a2=ab(1,i2)
              b2=ab(2,i2)
              call chunklength(k,funcurve,pars,a2,b2,rl2)
            endif

c
c       only check if self is larger than either of adjacent blocks,
c       iterating a couple times will catch everything
c
            ifsplit=0
            sc = 2.05d0
            if (rlself .gt. sc*rl1) ifsplit=1
            if (rlself .gt. sc*rl2) ifsplit=1
            if (ifsplit .eq. 0) goto 8600
c
c       split chunk i now, and recalculate nodes, ders, etc
c
            ifdone=0
            a=ab(1,i)
            b=ab(2,i)
            ab2=(a+b)/2
c
            i1=adjs(1,i)
            i2=adjs(2,i)
c        
            adjs(1,i) = i1
            adjs(2,i) = nch+1
c
c       . . . first update nch+1
c
            adjs(1,nch+1) = i
            adjs(2,nch+1) = i2
c
c       . . . if there's an i2, update it
c
            if (i2 .gt. 0) then
              adjs(1,i2) = nch+1
            endif
c
            nch=nch+1
            if (nch .gt. maxchunks) then
                call prinf('too many chunks in chunkfunc!*',done,0)
                call prinf('maxchunks=*',maxchunks,1)
                stop
            endif
c
            ab(1,i)=a
            ab(2,i)=ab2
c
            ab(1,nch)=ab2
            ab(2,nch)=b
c
 8600 continue

c
          if (ifdone .ne. 0) goto 9100
c
 9000 continue
 9100 continue

c
c       go ahead and oversample by nover, updating
c       the adjacency information adjs along the way
c
        if (nover .le. 1) goto 6100

c
        do 6000 ijk=1,nover-1
c
        nchold=nch
        do 5600 i=1,nchold
c
        a=ab(1,i)
        b=ab(2,i)
c
c       find ab2 using newton such that 
c       len(a,ab2)=len(ab2,b)=half the chunk length
c
        call chunklength(k,funcurve,pars,a,b,rl)
c
        rlhalf=rl/2
        thresh=1.0d-8
        ifnewt=0
        ab0=(a+b)/2
c
        do 6600 iter=1,1000
c
        call chunklength(k,funcurve,pars,a,ab0,rl1)
        call funcurve(ab0,pars,x,y,dx,dy,ddx,ddy)
        dsdt=sqrt(dx**2+dy**2)
        ab1=ab0-(rl1-rlhalf)/dsdt
c
        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifnewt=ifnewt+1
c
        if (ifnewt .eq. 3) goto 6700
        ab0=ab1
c
 6600 continue
 6700 continue
c
        if (ifnewt .lt. 3) then
            call prin2('newton failed! interval not split.*',done,0)
            stop
        endif
c
        ab2=ab1
c
        i1=adjs(1,i)
        i2=adjs(2,i)
        adjs(2,i)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=i
        adjs(2,nch+1)=i2
c
        ab(1,i)=a
        ab(2,i)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
 5600 continue
c
 6000 continue
 6100 continue
c
        if (ifclosed .eq. 1) goto 7700
c
c       if the curve is open, check the dyadic refinement at the
c       ends, first find the end segments
c
        do 6400 i =1,nch
        if (adjs(1,i) .lt. 0) ileft=i
        if (adjs(2,i) .lt. 0) iright=i
 6400 continue
c
        a1=ab(1,ileft)
        b1=ab(2,ileft)
        a2=ab(1,iright)
        b2=ab(2,iright)
c
        call chunklength(k,funcurve,pars,a1,b1,rlleft)
        call chunklength(k,funcurve,pars,a2,b2,rlright)
c
c       . . . dyadically split the left segment
c
        do 7000 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('ifclosed = *', ifclosed, 1)
            call prinf('boom refining left end! ijk=*',ijk,1)
            call prinf('boom refining left end! rlleft=*',rlleft,1)
            stop
        endif
c
        if (rlleft .le. chsmall) goto 7100
        a=ab(1,ileft)
        b=ab(2,ileft)
        ab2=(a+b)/2
c
        i1=adjs(1,ileft)
        i2=adjs(2,ileft)
        adjs(2,ileft)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=ileft
        adjs(2,nch+1)=i2
c
        ab(1,ileft)=a
        ab(2,ileft)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
        call chunklength(k,funcurve,pars,a,ab2,rlleft)
c        
 7000 continue
 7100 continue

c
c       . . . dyadically split the right segment
c
        do 7400 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('boom refining right end! ijk=*',ijk,1)
            call prinf('boom right end! rlright=*',rlright,1)
            stop
        endif
c
        if (rlright .le. chsmall) goto 7500
        a=ab(1,iright)
        b=ab(2,iright)
        ab2=(a+b)/2
c
        i1=adjs(1,iright)
        i2=adjs(2,iright)
        adjs(2,iright)=nch+1
c
        adjs(1,nch+1)=iright
        adjs(2,nch+1)=i2
c
        ab(1,iright)=a
        ab(2,iright)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        iright=nch
        ab(1,nch)=ab2
        ab(2,nch)=b
c
        call chunklength(k,funcurve,pars,ab2,b,rlright)
c
 7400 continue
 7500 continue
c
 7700 continue

c
c       up to here, everything has been done in parameter space, [ta,tb]
c       . . . finally evaluate the k nodes on each chunk, along with 
c       derivatives and chunk lengths
c
        do i = 1, nch
c
          a=ab(1,i)
          b=ab(2,i)
          hs(i)=(b-a)/2
c
          do j = 1, k
            t=a+(b-a)*(xs(j)+1)/2
            call funcurve(t, pars, x, y, dx, dy, dx2, dy2)
            chunks(1,j,i) = x
            chunks(2,j,i) = y
            ders(1,j,i) = dx
            ders(2,j,i) = dy
            ders2(1,j,i) = dx2
            ders2(2,j,i) = dy2
          enddo
        enddo
c
        return
        end
c-----------------------------------------------------------------------

        subroutine chunkfuncdens(eps, ifclosed, chsmall, ta, tb,
     1    funcurve,pars,nover,k0,nch0,nnh,aabb0,rdens1,rdens2, 
     2    k,nch,chunks,adjs, ders,ders2, hs)

        implicit real *8 (a-h,o-z)
        integer adjs(2,1),ifprocess(60000),nnh,k0,nch0
        double precision chunks(2,k,1),ders(2,k,1), ders2(2,k,1),
     1      xs(1000),ws(1000),xs2(1000),ws2(10000),u(10000),v(10000),
     2      u2(10000),v2(10000),ab(2,30000),hs(1),
     3      errs(1000),pars(1),
     4      errs0(1000)
       
        real *8 aabb0(2,1),rdens1(1),rdens2(1)
        real *8 xcoefs(1000), ycoefs(1000), xpcoefs(1000)
        real *8 ypcoefs(1000), xpcoefs_out(1000)
        real *8 ypcoefs_out(1000)
        real *8 ch7(2,1000), der7(2,1000)
c
        double precision, allocatable :: fvals(:,:), coefs(:,:)
c
        external funcurve
c
c       using a user-defined subroutine funcurve, split up the
c       curve into chunks. the calling sequence of funcurve
c       should be
c
c           funcurve(t,pars,x,y,dxdt,dydt,dxdt2,dydt2,k0,nch0,nnh,
c                    aabb0,rdens1,rdens2,f1,f2)
c
c       the routine will declare a chunk to be resolved when
c       chunks, ders, ders2, dsdt = sqrt(ders(1,)**2 + ders(2,)**2)
c       and the functions f1,f2
c       have length 2k expansions with the last k coefficients
c       having relative-rmse less than eps.
c
c
c       NOTE: no effort whatsoever has been made to make this routine
c       efficient. it does what it does very deliberately, with some
c       redundancy.
c
c       NOTE 2: the routine returns chunks that are at most a factor
c       of two different in archlength than adjacent chunks
c
c       input:
c
c         eps - absolute precision to resolve the curve
c         ifclosed - switch telling the routine whether the curve is
c           closed or not, ifclosed=1 means a closed curve
c         chsmall - if the curve is NOT closed, then this is the maximum
c           size of the intervals on either end. this parameter
c           is ignored if ifclosed=0.
c         ta,tb - assume that curve is parameterized by t \in [ta,tb)
c         funcurve - see above
c         pars - parameter array to send to funcurve, see above
c         nover - post-process oversampling factor, nover <= 1 will result
c             in no change, nover=2 will split each chunk in half (with
c             respect to arclength!!)
c         k - number of legendre nodes per chunk
c
c       output:
c
c         nch - total number of chunks created
c         chunks - points on each chunk, dimensioned chunks(2,k,nch)
c         adjs - adjacency information, adjs(1,i) is the chunk before
c             chunk i and adjs(2,i) is the chunk after chunk i
c         ders - first derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         ders2 - second derivative w.r.t. t on each chunk,
c             dimensioned chunks(2,k,nch)
c         hs - weight parameter to account for the arbitrary underlying
c             parameterization, used in arclength calculation and
c             subsequent integration
c
c
        done=1
        pi=4*atan(done)
c
        maxchunks=55000
        do i=1,maxchunks
          ifprocess(i)=0
        enddo


c
c       construct legendre nodes and weights, k and 2k of them, as well
c       as the interpolation/coefficients matrices
c
        itype=2
        call legeexps(itype,k,xs,u,v,ws)
        call legeexps(itype,2*k,xs2,u2,v2,ws2)

c
c       . . . start chunking
c
        ab(1,1)=ta
        ab(2,1)=tb
        nch=1
        ifdone=1
        adjs(1,1)=-1
        adjs(2,1)=-1
        nchnew=nch
c
        allocate(fvals(2*k,9))
        allocate(coefs(2*k,9))

c
        maxiter=10000
        do 5000 ijk = 1,maxiter

c
c       loop through all existing chunks, if resolved store, if not split
c
          ifdone=1
          do 4600 ich=1,nchnew

c
            if (ifprocess(ich) .eq. 1) goto 4600
            ifprocess(ich)=1

c
            a=ab(1,ich)
            b=ab(2,ich)


c
            do i=1,2*k
              t=a+(b-a)*(xs2(i)+1)/2


              call funcurve(t,pars,fvals(i,1),fvals(i,2),
     1          fvals(i,3),fvals(i,4),fvals(i,5),fvals(i,6),k0,nch0,
     2          nnh,aabb0,rdens1,rdens2,fvals(i,8),fvals(i,9))
              fvals(i,7)=sqrt(fvals(i,3)**2+fvals(i,4)**2)
            enddo
c
            do i=1,9
              call chunkmatvec(2*k,2*k,u2,fvals(1,i),coefs(1,i))
              errs0(i) = 0
              errs(i) = 0
            enddo
c
            do i=1,9
              do j=1,k
                errs0(i) = errs0(i) + coefs(j,i)**2
                errs(i) = errs(i) + coefs(k+j,i)**2
              enddo
            enddo
c
            rmsemax = -1.0d0
            do i=1,9
              errs(i) = sqrt(errs(i)/errs0(i)/k)
              if (errs(i) .gt. rmsemax) rmsemax=errs(i)
            enddo

c
c       . . . mark as processed and resolved if less than eps
c
            if (rmsemax .gt. eps) goto 2800
            goto 4600
            
c
 2800 continue

c
c       . . . if here, not resolved
c       divide - first update the adjacency list
c
            
            ifprocess(ich)=0
            ifdone=0
cccc            return
c
            if ((nch .eq. 1) .and. (ifclosed .gt. 0)) then
              adjs(1,nch)=2
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=1
            endif
c
            if ((nch .eq. 1) .and. (ifclosed .le. 0)) then
              adjs(1,nch)=-1
              adjs(2,nch)=2
              adjs(1,nch+1)=1
              adjs(2,nch+1)=-1
            endif
c
            if (nch .gt. 1) then
              iold2=adjs(2,ich)
              adjs(2,ich)=nch+1
              if (iold2 .gt. 0) adjs(1,iold2)=nch+1
              adjs(1,nch+1)=ich
              adjs(2,nch+1)=iold2
            endif

c
c       now update the endpoints in ab
c
            ab(1,ich)=a
            ab(2,ich)=(a+b)/2
c
            nch=nch+1
            if (nch .gt. maxchunks) then
              call prinf('too many chunks in chunkfunc!*',done,0)
              call prinf('maxchunks=*',maxchunks,1)
              stop
            endif
c
            ab(1,nch)=(a+b)/2
            ab(2,nch)=b
c
 4500 continue
 4600 continue
c
          if ((ifdone .eq. 1) .and. (nchnew .eq. nch)) goto 5100
          nchnew=nch        
c
 5000 continue
 5100 continue

c
c       the curve should be resolved to precision eps now on
c       each interval ab(,i)
c       check the size of adjacent neighboring chunks - if off by a
c       factor of more than 2, split them as well. iterate until done.
c

        maxiter=1000
        do 9000 ijk=1,maxiter
c
          nchold=nch
          ifdone=1
          do 8600 i=1,nchold
c
            i1=adjs(1,i)
            i2=adjs(2,i)

c
c       calculate chunk lengths
c
            a=ab(1,i)
            b=ab(2,i)
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a,b,rlself)
c
            rl1=rlself
            rl2=rlself
c
            if (i1 .gt. 0) then
              a1=ab(1,i1)
              b1=ab(2,i1)
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a1,b1,rl1)
            endif
c
            if (i2 .gt. 0) then
              a2=ab(1,i2)
              b2=ab(2,i2)
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a2,b2,rl2)
            endif

c
c       only check if self is larger than either of adjacent blocks,
c       iterating a couple times will catch everything
c
            ifsplit=0
            sc = 2.05d0
            if (rlself .gt. sc*rl1) ifsplit=1
            if (rlself .gt. sc*rl2) ifsplit=1
            if (ifsplit .eq. 0) goto 8600
c
c       split chunk i now, and recalculate nodes, ders, etc
c
            ifdone=0
            a=ab(1,i)
            b=ab(2,i)
            ab2=(a+b)/2
c
            i1=adjs(1,i)
            i2=adjs(2,i)
c        
            adjs(1,i) = i1
            adjs(2,i) = nch+1
c
c       . . . first update nch+1
c
            adjs(1,nch+1) = i
            adjs(2,nch+1) = i2
c
c       . . . if there's an i2, update it
c
            if (i2 .gt. 0) then
              adjs(1,i2) = nch+1
            endif
c
            nch=nch+1
            if (nch .gt. maxchunks) then
                call prinf('too many chunks in chunkfunc!*',done,0)
                call prinf('maxchunks=*',maxchunks,1)
                stop
            endif
c
            ab(1,i)=a
            ab(2,i)=ab2
c
            ab(1,nch)=ab2
            ab(2,nch)=b
c
 8600 continue

c
          if (ifdone .ne. 0) goto 9100
c
 9000 continue
 9100 continue

c
c       go ahead and oversample by nover, updating
c       the adjacency information adjs along the way
c
        if (nover .le. 1) goto 6100

c
        do 6000 ijk=1,nover-1
c
        nchold=nch
        do 5600 i=1,nchold
c
        a=ab(1,i)
        b=ab(2,i)
c
c       find ab2 using newton such that 
c       len(a,ab2)=len(ab2,b)=half the chunk length
c
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a,b,rl)
c
        rlhalf=rl/2
        thresh=1.0d-8
        ifnewt=0
        ab0=(a+b)/2
c
        do 6600 iter=1,1000
c
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a,ab0,rl1)
        call funcurve(ab0,pars,x,y,dx,dy,ddx,ddy,k0,nch0,nnh,aabb0,
     1                rdens1,rdens2,rrr1,rrr2)
        dsdt=sqrt(dx**2+dy**2)
        ab1=ab0-(rl1-rlhalf)/dsdt
c
        err=rl1-rlhalf
        if (abs(err) .lt. thresh) ifnewt=ifnewt+1
c
        if (ifnewt .eq. 3) goto 6700
        ab0=ab1
c
 6600 continue
 6700 continue
c
        if (ifnewt .lt. 3) then
            call prin2('newton failed! interval not split.*',done,0)
            stop
        endif
c
        ab2=ab1
c
        i1=adjs(1,i)
        i2=adjs(2,i)
        adjs(2,i)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=i
        adjs(2,nch+1)=i2
c
        ab(1,i)=a
        ab(2,i)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
 5600 continue
c
 6000 continue
 6100 continue
c
        if (ifclosed .eq. 1) goto 7700
c
c       if the curve is open, check the dyadic refinement at the
c       ends, first find the end segments
c
        do 6400 i =1,nch
        if (adjs(1,i) .lt. 0) ileft=i
        if (adjs(2,i) .lt. 0) iright=i
 6400 continue
c
        a1=ab(1,ileft)
        b1=ab(2,ileft)
        a2=ab(1,iright)
        b2=ab(2,iright)
c
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a1,b1,rlleft)
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a2,b2,rlright)
c
c       . . . dyadically split the left segment
c
        do 7000 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('ifclosed = *', ifclosed, 1)
            call prinf('boom refining left end! ijk=*',ijk,1)
            call prinf('boom refining left end! rlleft=*',rlleft,1)
            stop
        endif
c
        if (rlleft .le. chsmall) goto 7100
        a=ab(1,ileft)
        b=ab(2,ileft)
        ab2=(a+b)/2
c
        i1=adjs(1,ileft)
        i2=adjs(2,ileft)
        adjs(2,ileft)=nch+1
        if (i2 .gt. 0) adjs(1,i2)=nch+1
c
        adjs(1,nch+1)=ileft
        adjs(2,nch+1)=i2
c
        ab(1,ileft)=a
        ab(2,ileft)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        ab(1,nch)=ab2
        ab(2,nch)=b
c
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,a,ab2,rlleft)
c        
 7000 continue
 7100 continue

c
c       . . . dyadically split the right segment
c
        do 7400 ijk=1,1000
c
        if (ijk .gt. 100) then
            call prinf('boom refining right end! ijk=*',ijk,1)
            call prinf('boom right end! rlright=*',rlright,1)
            stop
        endif
c
        if (rlright .le. chsmall) goto 7500
        a=ab(1,iright)
        b=ab(2,iright)
        ab2=(a+b)/2
c
        i1=adjs(1,iright)
        i2=adjs(2,iright)
        adjs(2,iright)=nch+1
c
        adjs(1,nch+1)=iright
        adjs(2,nch+1)=i2
c
        ab(1,iright)=a
        ab(2,iright)=ab2
c
        nch=nch+1
        if (nch .gt. maxchunks) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            call prinf('maxchunks=*',maxchunks,1)
            stop
        endif
c
        iright=nch
        ab(1,nch)=ab2
        ab(2,nch)=b
c
            call chunklengthdens(k,funcurve,pars,k0,nch0,nnh,aabb0,
     1       rdens1,rdens2,ab2,b,rlright)
c
 7400 continue
 7500 continue
c
 7700 continue

c
c       up to here, everything has been done in parameter space, [ta,tb]
c       . . . finally evaluate the k nodes on each chunk, along with 
c       derivatives and chunk lengths

cccc        call prin2('ab=*',ab,2*nch)

        do i = 1, nch
c
          a=ab(1,i)
          b=ab(2,i)
          hs(i)=(b-a)/2
c
          do j = 1, k
            t=a+(b-a)*(xs(j)+1)/2
            call funcurve(t, pars, x, y, dx, dy, dx2, dy2,k0,nch0,nnh,
     1                    aabb0,rdens1,rdens2,rrr1,rrr2)

            chunks(1,j,i) = x
            chunks(2,j,i) = y
            ders(1,j,i) = dx
            ders(2,j,i) = dy
            ders2(1,j,i) = dx2
            ders2(2,j,i) = dy2
          enddo
        enddo

ccc        call prin2('hs=*',hs,nch)
c
        return
        end
c
c
c
c
c
        subroutine chunklength(k, funcurve, pars, a, b, rl)
        implicit real *8 (a-h,o-z)
        real *8 xs(1000),whts(1000),pars(1)
c
c       using a k-point gaussian quadrature, calculate the length
c       of the curve between a and b
c
        ifwhts=1
        call legewhts(k,xs,whts,ifwhts)
c
        rl=0
c
        do i=1,k
c
          t=a+(b-a)*(xs(i)+1)/2
          call funcurve(t,pars,x,y,xp,yp,xpp,ypp)
c
          dsdt=sqrt(xp**2+yp**2)*(b-a)/2
          rl=rl+dsdt*whts(i)
c
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunklengthdens(k,funcurve,pars,k0,nch0,nnh,
     1   ab0,rd1,rd2,a,b,rl)
        implicit real *8 (a-h,o-z)
        real *8 xs(1000),whts(1000),pars(1),rd1(1),rd2(1),ab0(2,1)
c
c       using a k-point gaussian quadrature, calculate the length
c       of the curve between a and b
c
        ifwhts=1
        call legewhts(k,xs,whts,ifwhts)
c
        rl=0
c
        do i=1,k
c
          t=a+(b-a)*(xs(i)+1)/2
          call funcurve(t,pars,x,y,xp,yp,xpp,ypp,k0,nch0,nnh,ab0,
     1                  rd1,rd2,rrr1,rrr2)
c
          dsdt=sqrt(xp**2+yp**2)*(b-a)/2
          rl=rl+dsdt*whts(i)
c
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunk_neighbors(ich, wgeo, ibefore, iafter)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(*)
c
c       returns the before and after neighbors of chunk ich
c
        iadjs=wgeo(4)
        call chunk_neighbors0(ich, wgeo(iadjs), ibefore, iafter)
c
        return
        end
c
c
        subroutine chunk_neighbors0(ich, adjs, ibefore, iafter)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
c
        ibefore = adjs(1,ich)
        iafter = adjs(2,ich)
        return
        end
c
c
c
c
c
        subroutine chunkunpack1(wgeo,k,nch,ichunks,iadjs,iders,
     1      iders2,ihs)
        implicit real *8 (a-h,o-z)
        real *8 wgeo(1)
c
c       returns the indices used to store chunks, ders, etc inside
c       wgeo
c
        k=wgeo(1)
        nch=wgeo(2)
        ichunks=wgeo(3)
        iadjs=wgeo(4)
        iders=wgeo(5)
        iders2=wgeo(6)
        ihs=wgeo(7)
c
        return
        end
c
c
c
c
c
        subroutine chunkunpack(wgeo,k,nch,chunks,adjs,ders,
     1      ders2,hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 wgeo(1),chunks(2,1),ders(2,1),ders2(2,1),hs(1)
c
c       unpacks data from the chunk storage array wgeo
c
        k=wgeo(1)
        nch=wgeo(2)
        ichunks=wgeo(3)
        iadjs=wgeo(4)
        iders=wgeo(5)
        iders2=wgeo(6)
        ihs=wgeo(7)
c
        lchunks=2*k*nch
        ladjs=2*nch
        lders=2*k*nch
        lders2=2*k*nch
        lhs=nch
c
c       copy the data
c
        call chunkcopy(lchunks,wgeo(ichunks),chunks)
        call chunkcopy(ladjs/2,wgeo(iadjs),adjs)
        call chunkcopy(lders,wgeo(iders),ders)
        call chunkcopy(lders2,wgeo(iders2),ders2)
        call chunkcopy(lhs,wgeo(ihs),hs)
c
        return
        end
c
c
c
c
c
        subroutine chunk2to1(n, xys, xs, ys)
        implicit real *8 (a-h,o-z)
        real *8 xys(2,1),xs(1),ys(1)
c
c       this routine pulls apart a 2 x n array into 2 separate
c       length n arrays
c
        do i=1,n
          xs(i)=xys(1,i)
          ys(i)=xys(2,i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkmatvec(m,n,a,x,y)
        implicit real *8 (a-h,o-z)
        real *8 a(m,n),x(1),y(1)
c
c       real matrix vector multiply
c
        do i=1,m
          dd=0
          do j=1,n
            dd=dd+a(i,j)*x(j)
          enddo
          y(i)=dd
        enddo
c
        return
        end
c
c
c
c
c
        subroutine chunkmatvec2(m,n,a,x,y)
        implicit real *8 (a-h,o-z)
        real *8 a(m,n)
        complex *16 x(1),y(1),dd
c
c       this routine multiplies a real matrix a by a complex
c       vector x, yielding a complex vector y
c
        do i=1,m
          dd=0
          do j=1,n
            dd=dd+a(i,j)*x(j)
          enddo
          y(i)=dd
        enddo
c
        return
        end
