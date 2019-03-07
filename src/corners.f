
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the code below contains several user-callable routines
c       for the construction of polygons and polygons with
c       smoothed corners. they are as follows:
c
c       corners_pack - 
c
c       chunkpolysmooth - given a set of vertices on a closed or
c           open curve, constructs a smoothed version of the geometry
c           using convolutional and adaptive discretization methods
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine corners_pack(ier, eps, widths, ibell, p1, p2, 
     1      i1, i2, nverts, verts, ifclosed, nover, k, wgeo, lused)
        implicit double precision (a-h,o-z)
        double precision widths(1), verts(2,1), wgeo(1)
c
        double precision, allocatable :: chunks(:,:,:), ders(:,:,:)
        double precision, allocatable :: hs(:), ders2(:,:,:)
        integer, allocatable :: adjs(:,:)
c
c       constructs a rounded polygon with all the packaged
c       featuers of chunkfunc_pack

c       input:
c         eps - 
c         widths - width to round at each corner
c         ibell - the type of rounding kernel to use
c           ibell = 1  <-> Gaussian
c           ibell = 2  <-> c_k (1-x**2)**k, where c_k is a normalizing constant
c         p1, p2 - real parameters for the kernel
c         i1, i2 - integer parameters for the kernel
c         nverts - number vertices, if ifclosed=1, then we expect
c           verts(*,1) = verts(*,nverts)
c         verts - (2,nverts) array of vertex locations, presumed to be
c           oriented counter-clockwise
c         ifclosed - 1 if closed, 0 if open
c         nover - oversampling parameter for each chunk
c         k - number of points per chunk
c
c       output:
c         ier - error code, 0 is good, not 0 is bad
c         wgeo - a work array storing the chunk information
c         lused - the length needed in wgeo in real *8 words
c
        maxchunks = 10000
        allocate(chunks(2,k,maxchunks))
        allocate(ders(2,k,maxchunks))
        allocate(ders2(2,k,maxchunks))
        allocate(hs(maxchunks))
        allocate(adjs(2,maxchunks))
c 
        call chunkpolysmooth(ier, eps, widths, ibell, p1, p2, 
     1    i1, i2, nverts, verts, ifclosed, nover, k, nch, chunks,
     2    adjs, ders, ders2, hs)
c
        if (nch .gt. maxchunks) then
          call prinf('nch = *', nch, 1)
          call prinf('error! bomb in corners_pack!*', maxchunks, 0)
        endif
c
c       now package up the data
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
        subroutine chunkpolysmooth(ier, eps, widths, ibell, p1, p2, 
     1      i1, i2, nverts, verts, ifclosed, nover, k, nch, chunks,
     2      adjs, ders, ders2, hs)
        implicit real *8 (a-h,o-z)
        integer adjs(2,1)
        real *8 widths(1), verts(2,1), chunks(2,k,1), ders(2,k,1),
     1      ders2(2,k,1), hs(1), xnodes(100), whts(100),
     2      pars(1000), xymid(10), cmass(10)
        real *8, allocatable :: evec(:,:)
        external fgauss, fbell
c
        allocate( evec(2,nverts) )
c
c       construct a polygon with rounded corners
c
c       input:
c
c         eps - the precision to which the curve should be resolved
c         widths - the size to cut out at each corner, should not be larger
c           than half of any adjacent line segment, should contain nverts
c           values, with widths(1)=widths(nverts) if a closed curve.
c           one should set widths(1)=widths(nverts)=0 if the curve is open.
c         ibell - the type of rounding to use
c             ibell=1    Gaussian
c             ibell=2    c_k (1-x**2)**k, where c_k is a normalizing constant
c             ibell=3    not implemented
c         p1,p2,i1,i2 - parameters for use in the bell
c             ibell=1    no parameters needed
c             ibell=2    i1=korder, the order of the bell
c         nverts - the number of vertices, if curve is closed then
c             we expect that verts(1) = verts(nverts)
c         verts - (2,nverts) array containing vertices in the plane, assumed
c             assumed to be oriented counterclockwise
c         ifclosed - set to 1 if the curve is closed, 0 if open
c         nover - number of times to oversample the chunks upon completion
c         k - number of gaussian nodes to put on each panel
c
c       output:
c
c         nch - number of chunks created
c         chunks - array of node locations on the chunks
c         adjs - adjancency information, adjs(1,i) is the chunk to the
c           chunk to the left of chunk i, adjs(2,i) to the right
c         ders - derivatives, scaled w.r.t. hs
c         ders2 - 2nd derivatives, scaled w.r.t. hs
c         hs - scaling parameters such that for t in [-1,1]
c
c       NOTE: in order to calculate dsdt relative to t in [-1,1], the
c         the scaling factors must be used:
c
c                   dsdt = sqrt(ders(1)**2+ders(2)**2)*hs
c
c
        done=1
        pi=4*atan(done)

c
c       check that widths is of proper size
c
        ier=0
        do 1400 i=1,nverts-1
c
        dx1=verts(1,i+1)-verts(1,i)
        dy1=verts(2,i+1)-verts(2,i)
        dx2=verts(1,i)-verts(1,i+1)
        dy2=verts(2,i)-verts(2,i+1)
c
        r1=sqrt(dx1**2+dy1**2)
        r2=sqrt(dx2**2+dy2**2)
c
        if (widths(i+1) .gt. r1/2) ier=2
        if (widths(i+1) .gt. r2/2) ier=2
c
 1400 continue
c
        if (ier .ne. 0) then
            call prinf('widths is too large, ier=*',ier,1)
            return
        endif

c
c       construct unit vectors which point along the edges,
c       will be useful
c
        do 1800 i=1,nverts-1
          dx1=verts(1,i+1)-verts(1,i)
          dy1=verts(2,i+1)-verts(2,i)
          r1=sqrt(dx1**2+dy1**2)
          evec(1,i)=dx1/r1
          evec(2,i)=dy1/r1
 1800 continue

c
c       construct all chunks and points, then scan them and
c       divide if necessary - start with the first segment, and
c       proceed - not the first corner
c
        ifwhts=1
        call legerts(ifwhts,k,xnodes,whts)
c
        if (ifclosed .eq. 0) then
            widths(1)=0
            widths(nverts)=0
        endif
c
        nch=0
        nedges=nverts-1
c

        do 5000 iseg=1,nedges
c
c       construct the single chunk on the middle part of the
c       edge
c
        x1=verts(1,iseg)+widths(iseg)*evec(1,iseg)
        y1=verts(2,iseg)+widths(iseg)*evec(2,iseg)
c
        x2=verts(1,iseg+1)-widths(iseg+1)*evec(1,iseg)
        y2=verts(2,iseg+1)-widths(iseg+1)*evec(2,iseg)
c
        nch=nch+1
c
        slopex=(x2-x1)
        slopey=(y2-y1)
        r=sqrt((x2-x1)**2+(y2-y1)**2)
        hs(nch)=r/2
c
        do 2200 j=1,k
c
        chunks(1,j,nch)=x1+(xnodes(j)+1)/2*slopex
        chunks(2,j,nch)=y1+(xnodes(j)+1)/2*slopey
c
        ders(1,j,nch)=slopex/r
        ders(2,j,nch)=slopey/r
c
        ders2(1,j,nch)=0
        ders2(2,j,nch)=0
 2200 continue
c
        if (iseg .eq. 1) then
          adjs(1,nch)=-1
          adjs(2,nch)=-1
          goto 3100
        endif


c
c       find the last segment created on a rounded corner
c       and update adjacencies
c
        do 2600 i=1,nch-1
        if (adjs(2,i) .lt. 0) then
            iright=i
            goto 2700
        endif
 2600 continue
 2700 continue
c
        adjs(2,iright)=nch
        adjs(1,nch)=iright
        adjs(2,nch)=-1
 3100 continue


        if ((ifclosed .ne. 1) .and. (iseg .eq. nedges)) goto 5100


c
c       . . . now construct a corner chunk
c
        x=verts(1,iseg+1)
        y=verts(2,iseg+1)
        w=widths(iseg+1)
c
        if (iseg .ne. nedges) then
            x3=verts(1,iseg+1)+w*evec(1,iseg+1)
            y3=verts(2,iseg+1)+w*evec(2,iseg+1)
        endif
c
        if (iseg .eq. nedges) then
            x3=verts(1,1)+w*evec(1,1)
            y3=verts(2,1)+w*evec(2,1)
        endif
c
        r23=sqrt((x3-x2)**2+(y3-y2)**2)
        b=sqrt(w*w-(r23/2)**2)
        a=-b/(r23/2)
        xmid=(x2+x3)/2
        ymid=(y2+y3)/2

c
c       calculate the bandwidth parameter h depending on which bell
c       specified by the user
c

        if (ibell .eq. 1) then
c
c           override the scale parameter so eps is 5.0 x 10^{-15}
c
            h = abs(b/a)/8
            goto 1235

cccc            call prin2('h = *', h, 1)
cccc            stop
c
            wd = r23/2
            thresh=eps/50
c
            h=wd/5
            do 1234 i=1,30
            h=sqrt(wd*wd/(-2*log(sqrt(2*pi)*thresh*h)))
            f1=1/sqrt(2*pi)/h*exp(-wd**2/2/h**2)
            if (abs(f1) .le. eps/10) goto 1235
 1234 continue
c
            call prinf('bell width did not converge!!! last h=*',h,1)
            call prin2('eps=*',eps,1)
            call prin2('f1=*',f1,1)
            stop
c
 1235 continue
c
        endif

c
        if (ibell .eq. 2) then
          h=abs(b/a)
          call prinf('ibell = *', ibell, 1)
          call prinf('not implemented and tested!!!!*', pi, 0)
          stop
        endif
c
c       create the chunk using either a gaussian bell or a finite bell
c
        ta=b/a
        tb=-b/a
        chsmall=1000
        nover5=1
        ifc=0
c
        
        if (ibell .eq. 1) then
c
            pars(1)=a
            pars(2)=b
cccc            pars(3)=h/2
            pars(3)=h
c
            call chunkfunc(eps,ifc,chsmall,ta,tb,fgauss,
     2          pars,nover5,k,nch5,chunks(1,1,nch+1),adjs(1,nch+1),
     3          ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
        endif
c

        if (ibell .eq. 2) then
c
            korder=i1
            pars(1)=a
            pars(2)=b
            pars(3)=h
            pars(4)=korder+0.1d0
c
            call chunkfunc(eps,ifc,chsmall,ta,tb,fbell,
     2          pars,nover5,k,nch5,chunks(1,1,nch+1),adjs(1,nch+1),
     3          ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
        endif

c
c       if left hand turn, reverse the orientation
c
        u1=x-x2
        u2=y-y2
        v1=x3-x
        v2=y3-y
c
        z=u1*v2-u2*v1
        if (z .ge. 0) ileft=1
        if (z .lt. 0) ileft=0
c
cccc        call prinf('test for left turn, ileft=*',ileft,1)
c
        if (ileft .eq. 1) then
            call chunkreverse(k,nch5,chunks(1,1,nch+1),adjs(1,nch+1),
     1          ders(1,1,nch+1),ders2(1,1,nch+1),hs(nch+1))
        endif

c
c       try new rotation scheme - calculate the vectors joining
c       the ends of the bump, and the ends of the polygon cut
c
        if (ileft .eq. 1) then
            yminusz1=x2-x3
            yminusz2=y2-y3
        endif
c
        if (ileft .eq. 0) then
            yminusz1=x3-x2
            yminusz2=y3-y2
        endif
c
        w1=abs(2*b/a)
        w2=0
c
        wlen=sqrt(w1**2+w2**2)
        yzlen=sqrt(yminusz1**2+yminusz2**2)
c
        phi=atan2(yminusz2,yminusz1)
cccc        if (ileft .eq. 0) phi=-phi

c
c       now rotate and translate the corner into place
c
        cphi=cos(phi)
        sphi=sin(phi)
c
        do i = nch+1,nch+nch5
          do j = 1,k
c
            chunks(2,j,i)=chunks(2,j,i)-b
c
            x7=chunks(1,j,i)
            y7=chunks(2,j,i)
            chunks(1,j,i)=cphi*x7-sphi*y7
            chunks(2,j,i)=sphi*x7+cphi*y7
c
            chunks(1,j,i)=chunks(1,j,i)+x
            chunks(2,j,i)=chunks(2,j,i)+y
c
            dx7=ders(1,j,i)
            dy7=ders(2,j,i)
            ders(1,j,i)=(cphi*dx7-sphi*dy7)
            ders(2,j,i)=(sphi*dx7+cphi*dy7)
c
            d2x7=ders2(1,j,i)
            d2y7=ders2(2,j,i)
            ders2(1,j,i)=cphi*d2x7-sphi*d2y7
            ders2(2,j,i)=sphi*d2x7+cphi*d2y7
c
          enddo
        enddo

c
        do i = nch+1,nch+nch5
          if (adjs(1,i) .gt. 0) adjs(1,i)=adjs(1,i)+nch
          if (adjs(2,i) .gt. 0) adjs(2,i)=adjs(2,i)+nch
        enddo
c
        do i = nch+1,nch+nch5
          if (adjs(1,i) .lt. 0) ileft=i
          if (adjs(2,i) .lt. 0) iright=i
        enddo
c
        adjs(1,ileft)=nch
        adjs(2,nch)=ileft
        nch=nch+nch5
c
 5000 continue
 5100 continue



c
c       if the curve is closed update the adjacency info
c
        if (ifclosed .ne. 1) goto 5700
c
        do i=1,nch
          if (adjs(1,i) .lt. 0) ileft=i
          if (adjs(2,i) .lt. 0) iright=i
        enddo
c
        adjs(1,ileft)=iright
        adjs(2,iright)=ileft
c
 5700 continue

c
c       now just check that no neighboring chunks are off by
c       more than a factor of two in arclength
c
        maxiter = 100000
        do ijk = 1,maxiter
c
          nchold = nch
          ifdone = 1
          dlen = -1
          id_chunk = -1
c
            do i = 1,nchold
c
c             find the longest chunk that requires splitting
c
              i1=adjs(1,i)
              i2=adjs(2,i)

c
c             calculate chunk lengths
c
              call chunksize(k,chunks(1,1,i),ders(1,1,i),hs(i),
     1          xymid,rad,cmass,crad,rlself)
c
              if (i1 .gt. 0) then
                call chunksize(k,chunks(1,1,i1),ders(1,1,i1),hs(i1),
     1            xymid,rad,cmass,crad,rl1)
              endif
c
              if (i2 .gt. 0) then
                call chunksize(k,chunks(1,1,i2),ders(1,1,i2),hs(i2),
     1            xymid,rad,cmass,crad,rl2)
              endif
c
c             only check if self is larger than either of adjacent blocks
c             iterating a couple times will catch everything
c
              ifsplit=0
              sc=2.05d0
cccc              sc=2
c
              if (i1 .gt. 0) then
                if (rlself .gt. sc*rl1) ifsplit=1
              endif
c
              if (i2 .gt. 0) then
                if (rlself .gt. sc*rl2) ifsplit=1
              endif
c
              if (ifsplit .ne. 0) then
                ifdone = 0
                if (rlself .gt. dlen) then
                  id_chunk = i
                  dlen = rlself
                endif
              endif
c
          enddo
c
          if (ifdone .eq. 1) goto 9100
c
c         if not done, split id_chunk
c
          call chunksplit1(id_chunk, k, nch, chunks, adjs, ders,
     1      ders2, hs)
c
        enddo
c
        call prinf('bomb! maxiter too small in smoothpoly!*', pi, 0)
        stop
c
 9100 continue

c        iw=77
c        itype=1
c        nnn = k*nch
c        call zpyplot(iw, chunks, nnn, itype, 'before oversampling*')
c        call prinf('before oversampling, nch = *', nch, 1)
c
cccc        stop

!!!!        call prin2('before oversampling, chunks = *', chunks, 2*nnn)

c
c       and finally oversample the curve by a factor of nover
c
        if (nover .le. 1) return
        niter=nover-1
c
        do i = 1,niter       
          nchold = nch
          do jj = 1,nchold
            call chunksplit1(jj, k, nch, chunks, adjs, ders,
     1        ders2, hs)
          enddo
        enddo
c
c        iw=78
c        itype=1
c        nnn = k*nch
c        call zpyplot(iw, chunks, nnn, itype, 'after oversampling*')
c        call prinf('after oversampling, nch = *', nch, 1)
c
cccc        call prin2('after oversampling, chunks = *', chunks, 2*nnn)
c        stop

        deallocate( evec )
c
        return
        end
c
c
c
c
c
        subroutine fgauss(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c       wrapper for the routine crn_fconvgauss
c
        a=pars(1)
        b=pars(2)
        h=pars(3)
c
        call crn_fconvgauss(t,a,b,h,val,der,der2)
        x=t
        y=val
c
        dxdt=1
        dydt=der
c
        dxdt2=0
        dydt2=der2
c
        return
        end
c
c
c
c
c
        subroutine crn_fconvgauss(x, a, b, h, val, der, der2)
        implicit real *8 (a-h,o-z)
c
c       this routine computes the convolution
c
c         ( a*abs(x)+b ) \star 1/(sqrt(2*pi)*h^2) exp(-x^2/(2*h^2))
c
c       this effectively smoothes off the corner from the abs(x) function
c
c
        done=1
        two=2
        pi=4*atan(done)
c
c       . . . formulas are computed via maple
c
        x2=x/sqrt(two)/h
        call qerrfun(x2,verf)
        val=a*x*verf+b+sqrt(two/pi)*a*h*exp(-x*x/two/h/h)
c
        fnorm=1/sqrt(2*pi)*exp(-x*x/2)
        der=a*verf
        der2=a*sqrt(two/pi)/h*exp(-x*x/two/h/h)
c
        return
        end
c
c
c
c
c
        subroutine fbell(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
        implicit real *8 (a-h,o-z)
        real *8 pars(1)
c
c       wrapper for the routine crn_fconvb1
c
        a=pars(1)
        b=pars(2)
        h=pars(3)
        k=pars(4)
c
        call crn_fconvb1(t,a,b,h,k,val,der,der2)
        x=t
        y=val
c
        dxdt=1
        dydt=der
c
        dxdt2=0
        dydt2=der2
c
        return
        end
c
c
c
c
c
        subroutine crn_fconvb1(x,a,b,h,k,val,der,der2)
        implicit real *8 (a-h,o-z)
        real *8 par1(100)
        external funcbe, funcbe1, funcbe2
c
c       computes the same convolution as in crn_fconvgauss
c       except using the bell computed in the subroutine bell1
c
c       the adaptive integration in this routine can be expensive...
c
        done=1

c
c       . . . use adaptive integration on the interval (x-h,x+h)
c
        rleft=x-h
        right=x+h
        zero=0
c
        m=16
        eps=1.0d-13
c
        par1(1)=x
        par1(2)=a
        par1(3)=b
        par1(4)=h
c
        if (rleft .ge. 0) then
            call adapgaus(ier,rleft,right,funcbe,par1,k,m,eps,
     1          rint1,maxrec,numint)
            call adapgaus(ier,rleft,right,funcbe1,par1,k,m,eps,
     1          rint2,maxrec,numint)
            call adapgaus(ier,rleft,right,funcbe2,par1,k,m,eps,
     1          rint3,maxrec,numint)
            val=rint1
            der=rint2
            der2=rint3
            return
        endif
c
        if (right .lt. 0) then
            call adapgaus(ier,rleft,right,funcbe,par1,k,m,eps,
     1          rint1,maxrec,numint)
            call adapgaus(ier,rleft,right,funcbe1,par1,k,m,eps,
     1          rint2,maxrec,numint)
            call adapgaus(ier,rleft,right,funcbe2,par1,k,m,eps,
     1          rint3,maxrec,numint)
            val=rint1
            der=rint2
            der2=rint3
            return
        endif

c
c       . . . if here, split into two pieces because of the abs
c       function's discontinuity
c
        call adapgaus(ier,rleft,zero,funcbe,par1,k,m,eps,
     1      rint1,maxrec,numint)
        call adapgaus(ier,rleft,zero,funcbe1,par1,k,m,eps,
     1      rint2,maxrec,numint)
        call adapgaus(ier,rleft,zero,funcbe2,par1,k,m,eps,
     1      rint3,maxrec,numint)
c
        call adapgaus(ier,zero,right,funcbe,par1,k,m,eps,
     1      rint4,maxrec,numint)
        call adapgaus(ier,zero,right,funcbe1,par1,k,m,eps,
     1      rint5,maxrec,numint)
        call adapgaus(ier,zero,right,funcbe2,par1,k,m,eps,
     1      rint6,maxrec,numint)
c
        val=rint1+rint4
        der=rint2+rint5
        der2=rint3+rint6
c
        return
        end
c
c
c
        function funcbe(t,par1,k)
        implicit real *8 (a-h,o-z)
        real *8 par1(1)
c
c       make sure that abs(x-t) .le. h !!!!! or all hell will break loose.
c
        x=par1(1)
        a=par1(2)
        b=par1(3)
        h=par1(4)
c
        xminust=x-t
        call bell1(xminust,h,k,val,der,der2)
        funcbe=(a*abs(t)+b)*val
c
        return
        end
c
c
c
        function funcbe1(t,par1,k)
        implicit real *8 (a-h,o-z)
        real *8 par1(1)
c
c       make sure that abs(x-t) .le. h !!!!! or all hell will break loose.
c
        x=par1(1)
        a=par1(2)
        b=par1(3)
        h=par1(4)
c
        xminust=x-t
        call bell1(xminust,h,k,val,der,der2)
        funcbe1=(a*abs(t)+b)*der
c
        return
        end
c
c
c
        function funcbe2(t,par1,k)
        implicit real *8 (a-h,o-z)
        real *8 par1(1)
c
c       make sure that abs(x-t) .le. h !!!!! or all hell will break loose.
c
        x=par1(1)
        a=par1(2)
        b=par1(3)
        h=par1(4)
c
        xminust=x-t
        call bell1(xminust,h,k,val,der,der2)
        funcbe2=(a*abs(t)+b)*der2
c
        return
        end
c
c
c
c
c
        subroutine bell1(x,h,k,val,der,der2)
        implicit real *8 (a-h,o-z)
c
c       this routine evaluates the bell
c
c           c_k * (1-(x/h)**2)**k * 1/h
c
c       along with its derivative (for no special reason).
c       c_k is chosen so that the integral is equal to 1, and
c       in fact,
c
c            c_k   =        \Gamma(k + 3/2) 
c                      -----------------------
c                       \Gamma(k + 1) * sqrt(pi)
c
c
        done=1
        pi=4*atan(done)
c
        if (abs(x) .gt. abs(h)) then
          val=0
          der=0
          der2=0
          return
        endif
c
c       evaluate 1/c_k
c
        t=k+1
        call gammanew_eval(t,gam1)
c
        t=k+1.5d0
        call gammanew_eval(t,gam2)
c
        ck=gam1*sqrt(pi)/gam2

        xh=x/h
        val=(1-xh**2)**k/h/ck
c
        der=-2*k*(1-xh**2)**(k-1)/h**3/ck*x
        der2=-2*k*(1-xh**2)**(k-1)/h**3/ck+
     1      4*k*(k-1)/(h**5)*(x**2)*(1-xh**2)**(k-2)
c
        return
        end
