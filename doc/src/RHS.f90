module define_rhs
    !! This is a module to determine the right-hand side for different model problems
    use mpi
    use comm_variable
    implicit none

contains
    subroutine RHS_2DCloseOff(b,xx,yy)
        !! 2D close-off problem with Dirichlet boundary condition
        implicit none

        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: b
        real(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        integer :: i,j

        do j=1,ny
            do i=1,nx

                b(i,j)=5.d0*pi*pi*(sin(pi*xx(i,j))*sin(2.d0*pi*yy(i,j))) &
                    -k0*k0*(sin(pi*xx(i,j))*sin(2.d0*pi*yy(i,j))+1.d0)
                
                if (npx == 0 .and. i == 1) then
                    b(i,j)=(1.d0,0.d0)
                endif

                if (npx == npx0-1 .and. i == nx) then
                    b(i,j)=(1.d0,0.d0)
                endif

                if (npy == 0 .and. j == 1) then
                    b(i,j)=(1.d0,0.d0)
                endif

                if (npy == npy0-1 .and. j == ny) then
                    b(i,j)=(1.d0,0.d0)
                endif

            enddo
        enddo
    end subroutine RHS_2DCloseOff

    subroutine RHS_CenterSource2D(b,xx,yy)
        !! 2D constant wavenumber problem with a central point source
        implicit none

        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: b
        real(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        integer :: i,j

        do j=1,ny
            do i=1,nx
                !! A point source located at (0.5,0.5)
                if (0.5d0-hx/2.d0 <= xx(i,j) .and. xx(i,j) <= 0.5d0+hx/2.d0 &
              .and. 0.5d0-hy/2.d0 <= yy(i,j) .and. yy(i,j) <= 0.5d0+hy/2.d0) then
                    b(i,j) = 1.d0/hx/hy
                else
                    b(i,j) = 0.d0
                endif

            enddo
        enddo
    end subroutine RHS_CenterSource2D

    subroutine RHS_2DWedge(b,xx,yy)
        !! 2D Wedge problem
        implicit none

        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: b
        real(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        integer :: i,j

        do j=1,ny
            do i=1,nx
                !! A point source located at (300,0)
                if (300.d0-hx/2.d0 <= xx(i,j) .and. xx(i,j) <= 300.d0+hx/2.d0 &
                 .and. 0.d0-hy/2.d0 <= yy(i,j) .and. yy(i,j) <= 0.d0+hy/2.d0) then
                    b(i,j) = 1.d0/hx/hy
                else
                    b(i,j) = 0.d0
                endif

            enddo
        enddo
    end subroutine RHS_2DWedge

    subroutine RHS_marmousi(b,xx,yy)
        !! Marmousi problem
        implicit none

        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: b
        real(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        integer :: i,j

        do j=1,ny
            do i=1,nx
                !! A point source located at (6000,0)
                if (6000.d0-hx/2.d0 <= xx(i,j) .and. xx(i,j) <= 6000.d0+hx/2.d0 &
                 .and. 0.d0-hy/2.d0 <= yy(i,j) .and. yy(i,j) <= 0.d0+hy/2.d0) then
                    b(i,j) = 1.d0/hx/hy
                else
                    b(i,j) = 0.d0
                endif

            enddo
        enddo
    end subroutine RHS_marmousi

end module define_rhs
!=====================================