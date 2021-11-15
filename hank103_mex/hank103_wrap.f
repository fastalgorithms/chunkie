        subroutine hank103_wrap(z,h0s,h1s,n)
        implicit real *8 (a-h,o-z)
        complex *16 z(n),h0s(n),h1s(n)

        ifexpon = 1

        do i=1,n
        call hank103(z(i),h0s(i),h1s(i),ifexpon)
        enddo

        return
        end
