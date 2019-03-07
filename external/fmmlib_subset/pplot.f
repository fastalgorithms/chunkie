c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code, and the beginning
c       of the python version of quaplot routines - and then some...
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine zpyplot(iw,zs,n,itype,title)
        implicit real *8 (a-h,o-z)
        real *8 zs(2,1)
        real *8, allocatable :: xs(:),ys(:)
        character *1 title(1)
c
        allocate(xs(n), ys(n))
c
        do 1600 i=1,n
        xs(i)=zs(1,i)
        ys(i)=zs(2,i)
 1600 continue
c
        call pyplot(iw,xs,ys,n,itype,title)
c
        return
        end
c
c
c
c
c
        subroutine zpyplot2(iw, zs, n, itype, zs2, n2, 
     1      itype2, title)
        implicit real *8 (a-h,o-z)
        real *8 zs(2,1),zs2(2,1)
        real *8, allocatable :: xs(:),ys(:),xs2(:),ys2(:)
        character *1 title(1)
c
        allocate(xs(n), ys(n), xs2(n2), ys2(n2))
c
        do 1600 i=1,n
        xs(i)=zs(1,i)
        ys(i)=zs(2,i)
 1600 continue
c
        do 1800 i=1,n2
        xs2(i)=zs2(1,i)
        ys2(i)=zs2(2,i)
 1800 continue
c
        call pyplot2(iw,xs,ys,n,itype,xs2,ys2,n2,itype2,title)
c
        return
        end
c
c
c
c
c
        subroutine zpyplot3(iw, zs, n, itype, zs2, n2, itype2,
     1      zs3, n3, itype3, title)
        implicit real *8 (a-h,o-z)
        real *8 zs(2,1), zs2(2,1), zs3(2,1)
        real *8, allocatable :: xs(:),ys(:),xs2(:),ys2(:)
        real *8, allocatable :: xs3(:), ys3(:)
        character *1 title(1)
c
        allocate(xs(n), ys(n), xs2(n2), ys2(n2))
        allocate(xs3(n), ys3(n))
c
        do i=1,n
          xs(i)=zs(1,i)
          ys(i)=zs(2,i)
        enddo
c
        do i=1,n2
          xs2(i)=zs2(1,i)
          ys2(i)=zs2(2,i)
        enddo
c
        do i=1,n3
          xs3(i) = zs3(1,i)
          ys3(i) = zs3(2,i)
        enddo
c
        call pyplot3(iw, xs, ys, n, itype, xs2, ys2, n2, itype2,
     1      xs3, ys3, n3, itype3, title)
c
        return
        end
c
c
c
c
c
        subroutine pyplot(iw,x,y,n,itype,title)
        implicit real *8 (a-h,o-z)
        real *8 x(1),y(1)
        character *1 a1,a10,file1(10),file11(9),title(1),
     1      temp(32)
        character *2 ls1
        character *9 file4
        character *10 file8
        character *32 title2
c
        equivalence (file1,file8), (file11,file4), (temp,title2)
c
c       plots the points determined by x,y using itype style
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)

        file8='plotiw.dat'
        file1(5)=a10
        file1(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1
c
c       print out the contents of the scripting file, plotiw.py
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c
        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c
        ls1='k.'
        if (itype .eq. 2) ls1='kx'
        if (itype .eq. 3) ls1='k-'
c
        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') ''
c
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a)') 'x = x.reshape(np.size(x)/2,2)'
        write(iun87,'(a)') 'x1=x[:,0]'
        write(iun87,'(a)') 'x2=x[:,1]'
c
        write(iun87,'(a,a,a)') 'pt.plot(x1,x2,"',ls1,'")'
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
cccc        write(iun87, '(a)') 'pt.axes().set_aspect("equal")'
        write(iun87,'(a)') 'pt.show()'

c
c       now print out the data file, plotiw.dat
c
        iun88=88
        open(unit=iun88,file=file8)
        do 1600 i=1,n
        write(iun88,*) x(i), y(i)
 1600 continue

c
        return
        end
c
c
c
c
c
        subroutine pyplot2(iw,x,y,n,itype,x2,y2,n2,itype2,title)
        implicit real *8 (a-h,o-z)
        real *8 x(1),y(1),x2(1),y2(1)
        character *1 a1,a10,file1(11),file2(11),file11(9),
     1      title(1),temp(32)
        character *2 ls1,ls2
        character *9 file4
        character *11 file8,file9
        character *32 title2
c
        equivalence (file1,file8), (file11,file4), (temp,title2),
     1      (file2,file9)
c
c       plots the points determined by x,y using itype style
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)
c
        file8='plotiw.dat1'
        file1(5)=a10
        file1(6)=a1
c
        file9='plotiw.dat2'
        file2(5)=a10
        file2(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1
c
c       print out the contents of the scripting file, plotiw.py
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c

        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c

        ls1='k.'
        if (itype .eq. 2) ls1='kx'
        if (itype .eq. 3) ls1='k-'
c
        ls2='g.'
        if (itype2 .eq. 2) ls2='gx'
        if (itype2 .eq. 3) ls2='g-'
c

        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') ''
c
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a)') 'x = x.reshape(np.size(x)/2,2)'
        write(iun87,'(a)') 'x1=x[:,0]'
        write(iun87,'(a)') 'x2=x[:,1]'
c
        write(iun87,'(a)') ''
        write(iun87,'(a,a,a)') 'y = np.loadtxt("',file9,'")'
        write(iun87,'(a)') 'y = y.reshape(np.size(y)/2,2)'
        write(iun87,'(a)') 'y1=y[:,0]'
        write(iun87,'(a)') 'y2=y[:,1]'
c
        write(iun87,'(a)') ''
        write(iun87,'(a,a,a)') 'pt.plot(x1,x2,"',ls1,'")'
        write(iun87,'(a,a,a)') 'pt.plot(y1,y2,"',ls2,'")'
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87, '(a)') 'pt.axes().set_aspect("equal")'
        write(iun87,'(a)') 'pt.show()'

c
c       now print out the data file, plotiw.dat
c
        iun88=888
        open(unit=iun88,file=file8)
        do 1600 i=1,n
        write(iun88,*) x(i), y(i)
 1600 continue
c
        iun89=899
        open(unit=iun89,file=file9)
        do 2600 i=1,n2
        write(iun89,*) x2(i), y2(i)
 2600 continue

c
        return
        end
c
c
c
c
c
        subroutine pyplot3(iw, x, y, n, itype, x2, y2, n2, itype2,
     1      x3, y3, n3, itype3, title)
        implicit real *8 (a-h,o-z)
        real *8 x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 a1,a10,file1(11),file2(11),file3(11),
     1      file11(9),title(1),temp(32)
        character *2 ls1,ls2,ls3
        character *9 file4
        character *11 file8,file9,file10
        character *32 title2
c
        equivalence (file1,file8), (file11,file4), (temp,title2),
     1      (file2,file9), (file3,file10)
c
c       plots the points determined by x,y using itype style
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)
c
        file8='plotiw.dat1'
        file1(5)=a10
        file1(6)=a1
c
        file9='plotiw.dat2'
        file2(5)=a10
        file2(6)=a1
c
        file10='plotiw.dat3'
        file3(5)=a10
        file3(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1
c
c       print out the contents of the scripting file, plotiw.py
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c

        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c

        ls1='k.'
        if (itype .eq. 2) ls1='kx'
        if (itype .eq. 3) ls1='k-'
c
        ls2='g.'
        if (itype2 .eq. 2) ls2='gx'
        if (itype2 .eq. 3) ls2='g-'
c
        ls3='r.'
        if (itype3 .eq. 2) ls2='rx'
        if (itype3 .eq. 3) ls2='r-'
c

        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') ''
c
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a)') 'x = x.reshape(np.size(x)/2,2)'
        write(iun87,'(a)') 'x1=x[:,0]'
        write(iun87,'(a)') 'x2=x[:,1]'
c
        write(iun87,'(a)') ''
        write(iun87,'(a,a,a)') 'y = np.loadtxt("',file9,'")'
        write(iun87,'(a)') 'y = y.reshape(np.size(y)/2,2)'
        write(iun87,'(a)') 'y1=y[:,0]'
        write(iun87,'(a)') 'y2=y[:,1]'
c
        write(iun87,'(a)') ''
        write(iun87,'(a,a,a)') 'z = np.loadtxt("',file10,'")'
        write(iun87,'(a)') 'z = z.reshape(np.size(z)/2,2)'
        write(iun87,'(a)') 'z1=z[:,0]'
        write(iun87,'(a)') 'z2=z[:,1]'
c
        write(iun87,'(a)') ''
        write(iun87,'(a,a,a)') 'pt.plot(x1,x2,"',ls1,'")'
        write(iun87,'(a,a,a)') 'pt.plot(y1,y2,"',ls2,'")'
        write(iun87,'(a,a,a)') 'pt.plot(z1,z2,"',ls3,'")'
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87, '(a)') 'pt.axes().set_aspect("equal")'
        write(iun87,'(a)') 'pt.show()'

c
c       now print out the data file, plotiw.dat
c
        iun88=888
        open(unit=iun88,file=file8)
        do 1600 i=1,n
        write(iun88,*) x(i), y(i)
 1600 continue
c
        iun89=899
        open(unit=iun89,file=file9)
        do 2600 i=1,n2
        write(iun89,*) x2(i), y2(i)
 2600 continue
c
        iun90=900
        open(unit=iun90,file=file10)
        do i=1,n3
          write(iun90,*) x3(i), y3(i)
        enddo

c
        return
        end
c
c
c
c
c
        subroutine pyimage(iw,m,n,vals,title)
        implicit real *8 (a-h,o-z)
        real *8 vals(m,n),target(2,m,n)
        character *1 a1,a10,file1(10),file11(9),title(1),
     1      temp(32), file1p(10)
        character *9 file4
        character *10 file8, file8p
        character *32 title2
c
        equivalence (file1,file8), (file11,file4), (temp,title2),
     1      (file1p,file8p)
c
c       output the data stored in vals to a file which python can
c       read and plot as a heat map using matplotlib
c
c       this routine requires a very special formatting of the inputs
c       target and vals - vals(1,1) is the value in the upper left
c       hand corner of the plot, vals(m,1) is the value in the lower
c       left hand corner, vals(1,n) is the value in the upper right
c       hand corner, and vals(m,n) is the value in the lower right 
c       hand corner.
c
c       input:
c
c         iw - file number, same as in quaplot
c         m - number of y values on the grid
c         n - number of x values on the grid
c         vals - 'z' value at the x,y grid points
c         title - the title of the plot, should be of the form 'title*',
c             and can only be a maximum of 32 characters long
c
c       note that on exit, two files will be created - plotiw.py 
c       and plotiw.dat
c
c
c       first construct the file names using iw
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)

        file8='plotiw.dat'
        file1(5)=a10
        file1(6)=a1
c
        file8p='plotiw.pdf'
        file1p(5)=a10
        file1p(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1
c
c       print out the contents of the scripting file, gniw
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c
        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue
c
        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') ''
        write(iun87,'(a)') 'pt.rc("font", size=16)'
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a,i5,a,i5,a)') 
     1      'x = x.reshape(', m, ',', n, ', order="F")'
        write(iun87,'(a)') 'pt.imshow(x, aspect="auto")'
        write(iun87,'(a)') 'cb=pt.colorbar(shrink=.9)'
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87,'(a,a,a)') 'pt.savefig("',file8p,'")'
        write(iun87,'(a)') 'pt.show()'

c
c       now print out the data file, plotiw.dat
c
        iun88=888
        open(unit=iun88,file=file8)
        do 1600 j=1,n
        do 1400 i=1,m
        write(iun88,*) vals(i,j)
 1400 continue
 1600 continue
c
        return
        end
c
c
c
c
c
        subroutine pyimage2(iw,m,n,vals,x,y,npts,x1,x2,
     1      y1,y2,title)
        implicit real *8 (a-h,o-z)
        real *8 vals(m,n),target(2,m,n),x(1),y(1)
        character *1 a1,a10,file1(11),file11(9),title(1),
     1      temp(32),file2(11),file1p(10)
        character *9 file4
        character *11 file8,file9
        character *10 file8p
        character *32 title2
c
        equivalence (file1,file8), (file11,file4), (temp,title2)
        equivalence (file2,file9), (file1p,file8p)
c
c       output the data stored in vals to a file which python can
c       read and plot as a heat map using matplotlib - this routine
c       also plots the scatterer, contains in x,y, as a filled
c       in polygon
c
c       this routine requires a very special formatting of the inputs
c       target and vals - vals(1,1) is the value in the upper left
c       hand corner of the plot, vals(m,1) is the value in the lower
c       left hand corner, vals(1,n) is the value in the upper right
c       hand corner, and vals(m,n) is the value in the lower right 
c       hand corner.
c
c       input:
c
c         iw - file number, same as in quaplot
c         m - number of y values on the grid
c         n - number of x values on the grid
c         vals - 'z' value at the x,y grid points
c         x,y - points on the scatterer
c         npts - length of x and y
c         title - the title of the plot, should be of the form 'title*',
c             and can only be a maximum of 32 characters long
c
c       note that on exit, two files will be created - plotiw.py 
c       and plotiw.dat
c
c
c       first construct the file names using iw
c
        i1=mod(iw,10)
        i10=(iw-i1)/10
        call int2char2(i1,a1)
        call int2char2(i10,a10)

        file8='plotiw.dat1'
        file1(5)=a10
        file1(6)=a1
c
        file8p='plotiw.pdf'
        file1p(5)=a10
        file1p(6)=a1
c
        file9='plotiw.dat2'
        file2(5)=a10
        file2(6)=a1
c
        file4='plotiw.py'
        file11(5)=a10
        file11(6)=a1
c
c       print out the contents of the scripting file, gniw
c
        iun87=877
        open(unit=iun87,file=file4)
c
        call quamesslen3(title,nchar)
c
        title2=''
        do 1200 i=1,nchar
        temp(i)=title(i)
 1200 continue

c
        write(iun87,'(a)') '#!/usr/bin/python'
        write(iun87,'(a)') 'import matplotlib.pyplot as pt'
        write(iun87,'(a)') 'import numpy as np'
        write(iun87,'(a)') ''
        write(iun87,'(a)') 'pt.rc("font", size=16)'
        write(iun87,'(a,a,a)') 'x = np.loadtxt("',file8,'")'
        write(iun87,'(a,i5,a,i5,a)') 
     1      'x = x.reshape(', m, ',', n, ', order="F")'
        write(iun87,'(a,f5.2,a,f5.2,a,f5.2,a,f5.2,a)')
     1      'pt.imshow(x, extent=[',x1,',',x2,',',y1,',',y2,'])'
        write(iun87,'(a)') 'cb=pt.colorbar(shrink=.9)'
        write(iun87,'(a,a,a)') 'pt.title("', title2, '")'
        write(iun87,'(a,a,a)') 'y = np.loadtxt("',file9,'")'
        write(iun87,'(a)') 'pt.fill(y[:,0],y[:,1],facecolor="w")'
        write(iun87,'(a,a,a)') 'pt.savefig("',file8p,'")'
        write(iun87,'(a)') 'pt.show()'

c
c       now print out the data file, plotiw.dat1 and plotiw.dat2
c
        iun88=888
        open(unit=iun88,file=file8)
        do 1600 j=1,n
        do 1400 i=1,m
        write(iun88,*) vals(i,j)
 1400 continue
 1600 continue
c
        call prinf('inside pyimage2, npts = *', npts, 1)
c
        iun89=889
        open(unit=iun89,file=file9)
        do 2600 i=1,npts
        write(iun89,*) x(i), y(i)
 2600 continue
c
        return
        end
c
c
c
c
c
        SUBROUTINE quamesslen3(MES,nchar)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c
        nchar=i1
c
        RETURN
        END
c
c
c
c
c
        subroutine int2char2(n,a)
        implicit real *8 (a-h,o-z)
        character *1 a
c
        if (n .eq. 0) a='0'
        if (n .eq. 1) a='1'
        if (n .eq. 2) a='2'
        if (n .eq. 3) a='3'
        if (n .eq. 4) a='4'
        if (n .eq. 5) a='5'
        if (n .eq. 6) a='6'
        if (n .eq. 7) a='7'
        if (n .eq. 8) a='8'
        if (n .eq. 9) a='9'
c
        return
        end
