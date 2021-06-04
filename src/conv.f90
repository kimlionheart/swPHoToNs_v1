!       interface to convolution by 2dcomp lib
!       Q. Wang 2017 Mar 9
subroutine get_local_size(global_size, vproc,&
                local_start, local_end, local_size)

use decomp_2d
use decomp_2d_fft
use MPI

implicit none

integer, dimension(2) :: vproc
integer, dimension(3) :: global_size, local_size, &
local_start,local_end

call decomp_2d_init(global_size(1),global_size(2),&
                global_size(3),vproc(1),vproc(2))

!        call decomp_2d_fft_init(PHYSICAL_IN_Z)
!        if(nrank.eq.0) print *,'nx=',global_size(1),', ny=',global_size(2), &
!        ',nz=',global_size(3), ', p_row=',vproc(1), ', p_col=',vproc(2)

        local_start = xstart -1 
        local_end   = xend - 1  
        local_size  = xsize

!        call decomp_2d_finalize
        return
end subroutine get_local_size


subroutine conv_pmonly(data, nside, param)
        use decomp_2d
        use decomp_2d_fft

        use MPI

        implicit none

        integer :: i, j, k, l,m,n, nhalf
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real*8, dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real*8,dimension(2) :: param
        !        real*8, dimension(2097152) :: data
        real*8, parameter :: M_PI = 3.141592653589793238462643383
        real(mytype) k2
        real(mytype) smooth, box, ismth2, gf, pref
        real(mytype), allocatable, dimension(:,:,:) :: in, in2
        complex(mytype), allocatable, dimension(:,:,:) :: out

        call decomp_2d_fft_init

        allocate(in(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))

!        write(*,*) xsize(1), xsize(2), xsize(3)
!        write(*,*) xstart(1), xstart(2), xstart(3)
!        write(*,*) xend(1), xend(2), xend(3)

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        in(i+xstart(1),j+xstart(2),k+xstart(3))=&
            data((i*xsize(2)+j)*xsize(3)+k+1)
        end do
        end do
        end do

        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

!        write(*,*) fft_start(1), fft_start(2), fft_start(3)
        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))


! compute r2c transform
        call decomp_2d_fft_3d(in,out)
        ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)
        ismth2 = 2*M_PI*smooth/(box)

        ismth2 = ismth2 * ismth2
        pref =  box * box / ( M_PI * nside(1) * nside(2) * nside(3) )

        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
                n = n - nside(3)
        endif
        
                do j=fft_start(2),fft_end(2),1
                m = j - 1
                if (m>nhalf) then
                        m = m - nside(2)
                endif
        
                        do, i=fft_start(1),fft_end(1),1
                        l = i - 1
                        if (l>nhalf) then
                                l = l - nside(1)
                        endif

                        k2 = REAL(n*n + m*m + l*l)

                        gf = pref/k2

                        if (l==0.and.m==0.and.n==0) then
                                gf = pref
                        endif

        
                        out(i,j,k) = out(i,j,k) * gf

                        end do
                end do
        end do

        ! compute c2r transform
        allocate(in2(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))

        !        call decomp_2d_fft_3d(out,in2)
        call decomp_2d_fft_3d(out,in2) 

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        data((i*xsize(2)+j)*xsize(3)+k+1) = &
                in2(i+xstart(1), j+xstart(2), k+xstart(3))
        end do
        end do
        end do

        deallocate(in,in2,out)
        call decomp_2d_fft_finalize
        !        write(*,*) smooth, box, M_PI, ismth2
        return
end subroutine conv_pmonly

subroutine convolution(data, nside, param)
        use decomp_2d
        use decomp_2d_fft

        use MPI

        implicit none

        integer :: i, j, k, l,m,n, nhalf
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real*8, dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real*8,dimension(2) :: param
        
        real*8, parameter :: M_PI = 3.1415926
!        real*8, parameter :: M_PI = 3.141592653589793238462643383
        real(mytype) k2
        real(mytype) smooth, box, ismth2, gf, pref
        real(mytype) fx, fy, fz, ff
        real(mytype), allocatable, dimension(:,:,:) :: in, in2
        complex(mytype), allocatable, dimension(:,:,:) :: out
        
        call decomp_2d_fft_init

        allocate(in(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))


        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
                in( i + xstart(1), j + xstart(2), k + xstart(3) )=&
                       data(((i)*xsize(2)+j)*xsize(3)+k+1)
!        in(i+xstart(1),j+xstart(2),k+xstart(3))=&
!                        data((i*xsize(2)+j)*xsize(3)+k+1)
        end do
        end do
        end do

        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))


        ! compute r2c transform
        call decomp_2d_fft_3d(in,out)

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)

!        write(*,*) smooth

        ismth2 = 2*M_PI*smooth/box
        ismth2 = ismth2 * ismth2
!        ismth2 = smooth/box
!        ismth2 = 2 * M_PI * ismth2 * ismth2
        pref =  box * box / ( M_PI * nside(1) * nside(2) * nside(3) )

!        write(*,*) smooth

        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
                n = n - nside(3)
        endif

            fz = M_PI*n/nside(3)
            fz = sin(fz)/fz

                if (n.EQ.0) then 
                        fz = 1
                endif        
            do j=fft_start(2),fft_end(2),1
            m = j - 1
            if (m>nhalf) then
                m = m - nside(2)
            endif

            fy = M_PI*m/nside(2)
            fy = sin(fy)/fy

                if (m.EQ.0) then 
                        fy = 1
                endif        
                do, i=fft_start(1),fft_end(1),1
                l = i - 1
                if (l>nhalf) then
                        l = l - nside(1)
                endif

                
                    fx = M_PI*l/nside(1)
                    fx = sin(fx)/fx
                if (l.EQ.0) then 
                        fx = 1
                endif        

                k2 = REAL(n*n + m*m + l*l)

                ff = 1.0/(fx*fy*fz)
                
!                ff = 1.0
                gf = pref*exp(-k2*ismth2)*ff*ff*ff*ff / k2

                if (l==0.and.m==0.and.n==0) then
                     gf = pref
                endif

                out(i,j,k) = out(i,j,k) * gf

               end do
            end do
        end do

        ! compute c2r transform
        allocate(in2(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))

        call decomp_2d_fft_3d(out,in2) 

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        data((i*xsize(2)+j)*xsize(3)+k+1) = &
                in2(i+xstart(1), j+xstart(2), k+xstart(3))
        end do
        end do
        end do

        deallocate(in,in2,out)
        call decomp_2d_fft_finalize
        return
end subroutine convolution



subroutine genicgrad2(data, nside, param )
        use decomp_2d
        use decomp_2d_fft

        use MPI

        implicit none

        integer :: i, j, k,dd, l,m,n,ll,mm,nn,ii,jj,kk, nhalf
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real*8, dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real*8,dimension(3) :: param
        real*8, parameter :: M_PI = 3.1415926
        real*8 k2, nor, kv
        real(mytype) smooth, box, ismth2, gf, pref, mag, delta, phase
        real(mytype) delta_re, delta_im
        real(mytype) fx, fy, fz, ff
        real(mytype), allocatable, dimension(:,:,:) :: in, in2
        complex(mytype), allocatable, dimension(:,:,:) :: out
        complex(mytype) tmps
        
        dd = param(3) 
        call decomp_2d_fft_init

        allocate(in(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))


        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
                in( i + xstart(1), j + xstart(2), k + xstart(3) )=&
                       data(((i)*xsize(2)+j)*xsize(3)+k+1)
        end do
        end do
        end do

        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))


        ! compute r2c transform
!        call decomp_2d_fft_3d(in,out)

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)     
        nor = 2*M_PI/(box)


        do k=fft_start(3),fft_end(3),1
           do j=fft_start(2),fft_end(2),1
                do, i=fft_start(1),fft_end(1),1
                        out(i,j,k) = 0
                end do
           end do
        end do


        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
                n = n-nside(3)
        endif

            do j=fft_start(2),fft_end(2),1
            m = j - 1
            if (m>nhalf) then
                m = m-nside(2)
            endif  


                do, i=fft_start(1),fft_end(1),1
                l = i - 1


                k2 = REAL(n*n + m*m + l*l)

                call ps_ic_k(k2*nor*nor, delta_re, delta_im)

                if (dd==1) then
                        kv=REAL(n*nor)
                endif

                if (dd==2) then
                        kv=REAL(m*nor)
                endif

                if (dd==3) then
                        kv=REAL(l*nor)
                endif

                gf = kv / ( k2  *nor *nor)
   
                if (l==0.and.m==0.and.n==0) then
                     gf = 0
                endif

                if (l>nhalf.and.m>nhalf.and.k>nhalf) then
                        gf = 0
                endif

                        out(i,j,k) = cmplx(-delta_re*gf, delta_im*gf)

               end do
            end do
        end do

        ! compute c2r transform
        allocate(in2(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))

        call decomp_2d_fft_3d(out,in2) 

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        data((i*xsize(2)+j)*xsize(3)+k+1) = &
                in2(i+xstart(1), j+xstart(2), k+xstart(3))
        end do
        end do
        end do

        deallocate(in,in2,out)
        call decomp_2d_fft_finalize
        return
end subroutine genicgrad2

subroutine densitykspace(data, nside, power, powcnt, param)
        use decomp_2d
        use decomp_2d_fft
        !        use glassman
        !        use decomp_2d_io

        use MPI

        implicit none

        integer :: i, j, k, l,m,n, nhalf, idx
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real*8, dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real*8,dimension(2) :: param
        real*8, dimension(64) :: power
        integer,dimension(64) :: powcnt

        !        real*8, dimension(2097152) :: data
        real*8, parameter :: M_PI = 3.141592653589793238462643383
        real(mytype) k2
        real(mytype) smooth, box, ismth2,pref
        real(mytype), allocatable, dimension(:,:,:) :: in
        complex(mytype), allocatable, dimension(:,:,:) :: out
        !        complex(mytype) :: ab = (1.,2.)
!        write(*,*) nside(1), nside(2), nside(3)

        !        ab = ab/2
        !        write(*,*) ab

!        write(*,*) xstart(1), xstart(2), xstart(3)
        call decomp_2d_fft_init

        allocate(in(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))

        !        write(*,*) xsize(1), xsize(2), xsize(3)
        !        write(*,*) xstart(1), xstart(2), xstart(3)
!        write(*,*) xend(1), xend(2), xend(3)

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        in(i+xstart(1),j+xstart(2),k+xstart(3))=&
                data((i*xsize(2)+j)*xsize(3)+k+1)
        end do
        end do
        end do

        !        do k=0,xsize(3)-1,1
        !        do j=0,xsize(2)-1,1
        !        do i=0,xsize(1)-1,1
        !        in(i+xstart(1),j+xstart(2),k+xstart(3))=&
!                data((i*xsize(2)+j)*xsize(3)+k+1)
        !        end do
        !        end do
        !        end do

call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

!        write(*,*) fft_start(1), fft_start(2), fft_start(3)
        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))


        ! compute r2c transform
        call decomp_2d_fft_3d(in,out)
        ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
        !        write(*,*) nside(1), nside(2), nside(3)
        !        write(*,*) fft_size(1), fft_size(2), fft_size(3)
!        write(*,*) fft_start(1), fft_start(2), fft_start(3)

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)
        ismth2 = 2*M_PI*smooth/(box)
        ismth2 = ismth2 * ismth2
        pref =  box * box / ( M_PI * nside(1) * nside(2) * nside(3) )
!        pref = -4*M_PI/(nside(1)*nside(2)*nside(3))

!        write(*,*) fft_start(1), fft_end(1)
        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
                n = n - nside(3)
        endif
        do j=fft_start(2),fft_end(2),1
        m = j - 1
        if (m>nhalf) then
                m = m - nside(2)
        endif
        do, i=fft_start(1),fft_end(1),1
        l = i - 1
        if (l>nhalf) then
                l = l - nside(1)
        endif

        k2 = REAL(n*n + m*m + l*l)

        idx = INT(sqrt(k2))

        if (idx < nside(1) ) then
        power(idx+1) = REAL(out(i,j,k))*REAL(out(i,j,k)) &
                + AIMAG(out(i,j,k))*AIMAG(out(i,j,k))
        powcnt(idx+1) = powcnt(idx+1)+1;
!                    write(*,*) idx, power(idx+1) , powcnt(idx+1), out(i,j,k)
        endif

        end do
        end do
        end do

        ! compute c2r transform
        !        allocate(in2(xstart(1):xend(1), &
                                !                xstart(2):xend(2),xstart(3):xend(3)))

        !        call decomp_2d_fft_3d(out,in2)
!        call decomp_2d_fft_3d(out,in2) 

        deallocate(in,out)
        call decomp_2d_fft_finalize
        !        write(*,*) smooth, box, M_PI, ismth2
        return
end subroutine densitykspace

