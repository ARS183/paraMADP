module define_BC
    !! This is a module to deal with the computational stencils for the Dirichlet of Sommerfeld boundary conditions
    use comm_variable
    use wavenumber
    implicit none

contains
    subroutine if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy)
        !! Determine which BC, for specified coarse grid system
        implicit none

        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
        integer,intent(in) :: i,j,ni,nj
        real(kind = realdp),intent(in) :: hx_c,hy_c,kxy

        if (flag_BCs == 1) then
            call Dirichlet(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c)
        elseif (flag_BCs == 2) then
            call first_order_Neum(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy)
        else
            write(*,*) "There is no such a boundary conditions yet! Use first_order_Neum instead!!"
            call first_order_Neum(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy)
        endif
    end subroutine if_BCs

    subroutine Dirichlet(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c)
        !! Dirichlet, ap=h^2, which will devided by h^2 in the rountine of the operator
        !! When use this Dirichlet, you have to be carefull with the RHS respectively
        implicit none

        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
        integer,intent(in) :: i,j,ni,nj
        real(kind = realdp),intent(in) :: hx_c,hy_c

        if (npx == 0 .and. i==1) then
            ap = hx_c*hy_c + czero
            ae = 0.d0
            aw = 0.d0
            an = 0.d0
            as = 0.d0
        endif

        if (npx == npx0-1 .and. i==ni) then
            ap = hx_c*hy_c + czero
            ae = 0.d0
            aw = 0.d0
            an = 0.d0
            as = 0.d0
        endif

        if (npy == 0 .and. j==1) then
            ap = hx_c*hy_c + czero
            ae = 0.d0
            aw = 0.d0
            an = 0.d0
            as = 0.d0
        endif

        if (npy == npy0-1 .and. j==nj) then
            ap = hx_c*hy_c + czero
            ae = 0.d0
            aw = 0.d0
            an = 0.d0
            as = 0.d0
        endif

    end subroutine Dirichlet

    subroutine first_order_Neum(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy)
        !! Sommerfeld boundary conditions. Eliminate the ghost grid points by second-order discretization
        implicit none

        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
        integer,intent(in) :: i,j,ni,nj
        real(kind = realdp),intent(in) :: hx_c,hy_c,kxy

        !================================================
        if (npx == 0 .and. i==1) then
            aw = 0.d0
            ae = -2.d0
            ap = ap - 2.d0*kxy*hx_c*cone !Second-order discretization
        endif

        if (npx == npx0-1 .and. i==ni) then
            ae = 0.d0
            aw = -2.d0
            ap = ap - 2.d0*kxy*hx_c*cone!Second-order discretization
        endif

        if (npy == 0 .and. j==1) then
            as = 0.d0
            an = -2.d0
            ap = ap - 2.d0*kxy*hy_c*cone !Second-order discretization
        endif

        if (npy == npy0-1 .and. j==nj) then
            an = 0.d0
            as = -2.d0
            ap = ap - 2.d0*kxy*hy_c*cone !Second-order discretization
        endif

    end subroutine first_order_Neum

end module define_BC
!=====================================