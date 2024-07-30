module analytical_sol
  !! This is a module to compute analytical solutions for some model problems
  use mpi
  use comm_variable
  implicit none

contains
  subroutine exact_2DCloseOff(u_ex,xx,yy)
      !! The analytical solution for 2D close-off problem
      implicit none

      complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u_ex
      real   (kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
      integer :: i,j

      do j=1,ny
          do i=1,nx
              u_ex(i,j) = dsin(pi*xx(i,j))*dsin(2.d0*pi*yy(i,j))+1.d0
          enddo
      enddo

  end subroutine exact_2DCloseOff

  subroutine exact_2DPointSource_1stSomer(u_ex,xx,yy)
    !! The analytical solution for constant wavenumber problem with central source point with Sommerfeld boundary condition
    implicit none

    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u_ex
    real   (kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
    integer :: i,j
    real (kind = realdp) :: r

    do j=1,ny
        do i=1,nx
            r = dsqrt((xx(i,j)-0.5d0)**2+(yy(i,j)-0.5d0)**2)
            if (r == 0.d0) then
              r=hx/2.d0
              u_ex(i,j) = cone/4.d0*(bessel_j0(k0*r)+cone*bessel_y0(k0*r))  !better that twice!
            else
              u_ex(i,j) = cone/4.d0*(bessel_j0(k0*r)+cone*bessel_y0(k0*r))
            endif
        enddo
    enddo

  end subroutine exact_2DPointSource_1stSomer

  subroutine exact_2DPointSource_Dirichlet(u_ex,xx,yy)
    !! The analytical solution for constant wavenumber problem with central source point with Dirichlet boundary condition
    ! Becareful it may take a lot of time to compute
    implicit none

    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u_ex
    real   (kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
    integer :: ii,jj

    do jj = 1, ny_global
      do ii = 1, nx_global
                u_ex = u_ex + 4*(sin(jj*pi*xx)*sin(ii*pi*yy)*sin(ii*pi*0.5d0)*sin(jj*pi*0.5d0)) &
                            /(jj**2*pi**2 + ii**2*pi**2 - k0**2)
      enddo
    enddo

  end subroutine exact_2DPointSource_Dirichlet

end module analytical_sol


!=====================================