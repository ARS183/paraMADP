module smoother
    !! A smoother module 
    use comm_variable
    use mpi_setup
    use define_BC
    use operators
    implicit none

contains
    subroutine Damp_Jacobi_smoother(x, rhs, ni, nj, hx_c, hy_c, kxy, kh2)
        !! This is a routine to apply damped Jacobbi smoother on coarse grid systems.

        implicit none
        integer, intent(in) :: ni,nj
            !! subdomain grid size on the current coarse level
        real(kind=realdp),intent(in) :: hx_c,hy_c
            !! coarse-grid space step size
        real   (kind=realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(in) :: kxy, kh2
            !! Wavenumber and (kh)^2
        complex(kind=realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP)  :: x
        complex(kind=realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP)  :: rhs

        !! local arguements
        complex(kind=realdp), allocatable, dimension(:,:) :: x_tmp 
        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
            !! Computational stencils
        integer :: k,sm_step,i,j

        real(kind=8) :: omega
            !! Relaxation factor, 0.8 by default

        sm_step=1
        allocate(x_tmp(1-LAP:ni+LAP,1-LAP:nj+LAP))

        omega=0.8d0
        x_tmp=(0.d0,0.d0)

        do k=1,sm_step
            call mg_check_xy2d(x,ni,nj) ! data communication, coarse level
            
            j = 1
            do i=1,ni
                call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                !==============================
                x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                    +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
            enddo
        
            j = nj
            do i=1,ni
                call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                !==============================
                x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                    +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
            enddo
        
            i = 1
            do j=2,nj-1
                call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                !==============================
                x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                    +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
            enddo
        
            i = ni
            do j=2,nj-1
                call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                !==============================
                x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                    +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
            enddo
            
            if (flag_BCs == 1) then
                do j=2,nj-1
                    do i=2,ni-1
                        call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                        call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                        !==============================
                        x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                            +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
                    enddo
                enddo
            elseif (flag_BCs == 2) then
                do j=2,nj-1
                    do i=2,ni-1
                        call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                        !==============================
                        x_tmp(i,j) = (rhs(i,j)-(ae*x(i+1,j)+aw*x(i-1,j) &
                                            +an*x(i,j+1)+as*x(i,j-1))/hx_c/hy_c)/(ap/hx_c/hy_c)
                    enddo
                enddo
            endif

            x=omega*x_tmp+(1-omega)*x
        enddo

        deallocate(x_tmp)

    end subroutine Damp_Jacobi_smoother

end module smoother
