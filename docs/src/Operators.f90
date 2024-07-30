module operators
    !! This is a module to implement some matrix-free operations
    use mpi
    use comm_variable
    use mpi_setup
    use wavenumber
    use define_BC
    implicit none

contains

    ! FUNCTION NORM =====================================================================
    real(kind = realdp) function norm(v)
        !! Compute the L2-norm of a complex variable, ONLY for the default grid system.
        ! Note that the dimension is not specified here, it means it calculate all the elements
        !  without exception. Thus, we must be careful about the input range of the vector v, as
        !  it was extended to (1-LAP:nx+LAP)
        implicit none
        real   (kind = realdp) :: norm_local
        complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP), intent(in) :: v
        integer :: j

        norm_local = 0.d0
        do j=1,ny
            norm_local = norm_local + real(dot_product(v(1:nx,j),v(1:nx,j)))
        enddo

        call MPI_Allreduce(norm_local,       &
                        norm,                &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

        norm = dsqrt(norm)
        return
    end function norm

    real(kind = realdp) function Realnorm(v)
        !! Compute the L2-norm of a real variable, ONLY for the default grid system.
        ! Note that the dimension is not specified here, it means it calculate all the elements
        !  without exception. Thus, we must be careful about the input range of the vector v, as
        !  it was extended to (1-LAP:nx+LAP)
        implicit none

        real(kind = realdp)                           :: norm_local
        real(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP), intent(in) :: v
        integer :: j

        norm_local = 0.d0
        do j=1,ny
            norm_local = norm_local + dot_product(v(1:nx,j),v(1:nx,j))
        enddo

        call MPI_Allreduce(norm_local,       &
                        Realnorm,            &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

        Realnorm = dsqrt(Realnorm)
        return
    end function Realnorm

    ! FUNCTION DOT_PROD =================================================================
    complex(kind = realdp) function dot_prod(v,w)
        !! Compute the do product of two complex variables, ONLY for the default grid system.
        ! Note that the dimension is not specified here, it means it calculate all the elements
        !  without exception. Thus, we must be careful about the input range of the vector v, as
        !  it was extended to (1-LAP:nx+LAP)
        implicit none

        complex(kind = realdp)                           :: dot_prod_local
        complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP), intent(in) :: v, w
        integer :: j

        dot_prod_local = czero
        dot_prod = czero

        do j= 1,ny
            dot_prod_local = dot_prod_local + dot_product(v(1:nx,j),w(1:nx,j))
        enddo

        call MPI_Allreduce(dot_prod_local,   &
                        dot_prod,            &
                        1,                   &
                        MPI_DOUBLE_COMPLEX,  &
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)
        return
    end function dot_prod

    ! FUNCTION NORM =====================================================================
    real(kind = realdp) function mg_norm(v,ni,nj)
        !! Compute the L2-norm of a complex variable, for a coarse-grid system.
        ! Note that the dimension is not specified here, it means it calculate all the elements
        !  without exception. Thus, we must be careful about the input range of the vector v, as
        !  it was extended to (1-LAP:ni+LAP)
        implicit none

        integer, intent(in) :: ni,nj
            !! ni, nj: grid size of the subdomain on the current coarse grid level
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(in) :: v
        real   (kind = realdp) :: norm_local
        integer :: j

        norm_local = 0.d0
        do j=1,nj
            norm_local = norm_local + real(dot_product(v(1:ni,j),v(1:ni,j)))
        enddo

        call MPI_Allreduce(norm_local,       &
                        mg_norm,             &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

        mg_norm = dsqrt(mg_norm)
        return
    end function mg_norm


    ! FUNCTION DOT_PROD =================================================================
    complex(kind = realdp) function mg_dot_prod(v,w,ni,nj)
        !! Compute the dot product of two complex variables, for a coarse-grid system.
        !! Note that the dimension is not specified here, it means it calculate all the elements
        !  without exception. Thus, we must be careful about the input range of the vector v, as
        !  it was extended to (1-LAP:ni+LAP)
        implicit none

        integer, intent(in) :: ni,nj
            !! ni, nj: grid size of the subdomain on the current coarse grid level
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(in) :: v, w
        complex(kind = realdp) :: dot_prod_local
        integer :: j

        dot_prod_local = czero
        mg_dot_prod = czero
        do j=1, nj
            dot_prod_local = dot_prod_local + dot_product(v(1:ni,j),w(1:ni,j))
        enddo

        call MPI_Allreduce(dot_prod_local,   &
                        mg_dot_prod,         &
                        1,                   &
                        MPI_DOUBLE_COMPLEX,  &
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)
        return

    end function mg_dot_prod
    !=============================================================================================================

    real(kind=realdp) function LOG2(x)
        !! Compute log2(x)
        implicit none
        real(kind=realdp), intent(in) :: x

        LOG2 = log(x) / log(2.d0)
    end function LOG2


    !=============================================================================================================
    !==                                             Linear operators                                            ==
    !=============================================================================================================

    subroutine Helmholtz2d_stencils(ap,an,as,aw,ae,hxhykxy2)
        !! Computational stencils for the Helmholtz operator, second-order central FD scheme
        implicit none

        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
        real(kind = realdp),intent(in) :: hxhykxy2
            !! The value of (kh)^2, which is computed and stored in advance.

        ap = 4.d0 - hxhykxy2 + czero
        ae = -1.d0
        aw = -1.d0
        an = -1.d0
        as = -1.d0
        
    end subroutine Helmholtz2d_stencils

    subroutine cslp2d_stencils(ap,an,as,aw,ae,hxhykxy2)
        !! Computational stencils for the CSLP operator
        implicit none

        complex(kind = realdp) :: ap
        real   (kind = realdp) :: ae,aw,an,as
        real(kind = realdp),intent(in) :: hxhykxy2
            !! The value of (kh)^2, which is computed and stored in advance

        ap = 4.d0 -(beta1 - beta2*cone)*hxhykxy2
        ae = -1.d0
        aw = -1.d0
        an = -1.d0
        as = -1.d0
        
    end subroutine cslp2d_stencils

    subroutine Helmholtz2d_BC(v_in,v_out)
        !! matrix-free Helmholtz operator, ONLY for the default gird system
        !! First deal with the boundary grid points and then the internals
        implicit none

        complex(kind = realdp),    dimension(1-LAP:nx+LAP,1-LAP:ny+LAP), intent(inout)   :: v_in
        complex(kind = realdp),    dimension(1-LAP:nx+LAP,1-LAP:ny+LAP), intent(inout)   :: v_out
        integer :: i,j
        real(kind = realdp) :: ae,aw,an,as
        complex(kind = realdp) :: ap

        v_out = czero
        call check_xy2d(v_in)
        
        j = 1
        do i=1,nx
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,nx,ny,hx,hy,wavenumber_k(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        j = ny
        do i=1,nx
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,nx,ny,hx,hy,wavenumber_k(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = 1
        do j=2,ny-1
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,nx,ny,hx,hy,wavenumber_k(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = nx
        do j=2,ny-1
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,nx,ny,hx,hy,wavenumber_k(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
        
        if (flag_BCs == 1) then
            do j=2,ny-1
                do i=2,nx-1
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
                    call if_BCs(ap,an,as,aw,ae,i,j,nx,ny,hx,hy,wavenumber_k(i,j)) !Actually no need anymore
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        elseif (flag_BCs == 2) then
            do j=2,ny-1
                do i=2,nx-1
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,kh_pow_2(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        endif

        v_out=v_out/hxhy

    end subroutine Helmholtz2d_BC

    subroutine Helmholtz2d_BC_mg(v_in,v_out,ni,nj,hx_c,hy_c,kxy,kh2)
        !! matrix-free Helmholtz operator, for the specified coarse gird systems
        implicit none

        integer,intent(in) :: ni,nj
        real(kind = realdp),intent(in) :: hx_c,hy_c
        real   (kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(in)      :: kxy, kh2
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(inout)   :: v_in
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(inout)   :: v_out
        
        real(kind = realdp) :: hxhy_c
        real(kind = realdp) :: ae,aw,an,as
        integer :: i,j
        complex(kind = realdp) :: ap

        if (hx_c /= hy_c) then
            write(*,*) "careful, hx not equals to hy, in helmholtzOP"
        endif
        hxhy_c = hx_c*hy_c

        v_out  = czero
        call mg_check_xy2d(v_in, ni, nj)

        j = 1
        do i=1,ni
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        j = nj
        do i=1,ni
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = 1
        do j=2,nj-1
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = ni
        do j=2,nj-1
            call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo

        if (flag_BCs == 1) then
            do j=2,nj-1
                do i=2,ni-1
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                    call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        elseif (flag_BCs == 2) then
            do j=2,nj-1
                do i=2,ni-1
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        endif

        v_out=v_out/hxhy_c

    end subroutine Helmholtz2d_BC_mg

    subroutine CSLP_OP_BC(v_in,v_out,ni,nj,hx_c,hy_c,kxy,kh2)
        !! matrix-free CSLP operator, for the specified gird systems
        implicit none

        integer        , intent(in) :: ni, nj
            !! ni,nj: grid size of subdomain on the current grid level
        real   (kind = realdp), intent(in) :: hx_c, hy_c
            !! hx_c,hy_c: space step of subdomain on the current grid level
        real   (kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(in) :: kxy,kh2
            !! kxy: wavenumber profile of subdomain on the current grid level,kh2=(kh)^2
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(inout) :: v_in
        complex(kind = realdp), dimension(1-LAP:ni+LAP,1-LAP:nj+LAP), intent(inout) :: v_out
        
        complex(kind = realdp) :: ap
        real   (kind = realdp) :: hxhy_c
        real   (kind = realdp) :: ae,aw,an,as
        integer :: i,j

        if (hx_c /= hy_c) then
            write(*,*) "careful, hx not equals to hy, in CSLP op"
        endif
        hxhy_c= hx_c*hy_c

        v_out = czero
        call mg_check_xy2d(v_in,ni,nj)

        j = 1
        do i=1,ni
            call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo

        j = nj
        do i=1,ni
            call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = 1
        do j=2,nj-1
            call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
    
        i = ni
        do j=2,nj-1
            call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
            call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
            !==============================
            v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                    +an*v_in(i,j+1)+as*v_in(i,j-1)
        enddo
        
        if (flag_BCs == 1) then
            do j=2,nj-1
                do i=2,ni-1
                    call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                    call if_BCs(ap,an,as,aw,ae,i,j,ni,nj,hx_c,hy_c,kxy(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        elseif (flag_BCs == 2) then
            do j=2,nj-1
                do i=2,ni-1
                    call cslp2d_stencils(ap,an,as,aw,ae,kh2(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                           +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            enddo
        endif

        v_out=v_out/hxhy_c

    end subroutine CSLP_OP_BC

    function Helm_Ax_nth(v_in,grid) 
        !! *ReD-Glk* Helmholtz operator for coarse levels
        implicit none

        type(Gridpara), intent(inout) :: grid
            !! Parameters of the current grid system
        complex(kind = realdp), intent(inout) :: v_in(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)
        complex(kind=realdp) :: Helm_Ax_nth(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)

        real(kind = realdp) :: ae,aw,an,as
        complex(kind = realdp) :: ap, coefA(-3:3,-3:3)
            !! computational stencils
        integer :: i,j,ia,ja
        integer(kind=4) :: level_i
        
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1

        Helm_Ax_nth   = czero
        v_in(1-LAP:0,:) = czero
        v_in(:,1-LAP:0) = czero
        v_in(grid%nx+1:grid%nx+LAP,:) = czero
        v_in(:,grid%ny+1:grid%ny+LAP) = czero

        call mg_check_xy2d(v_in, grid%nx, grid%ny)
        call GridBase_ExtrpltGhostBCs(v_in,grid)

        do j=1,grid%ny
            do i=1,grid%nx
                call Helm_Ax_nth_stencils(coefA,grid,i,j)
                do ja=-3,3
                    do ia=-3,3
                        Helm_Ax_nth(i,j)=Helm_Ax_nth(i,j)+coefA(ia,ja)*v_in(i+ia,j+ja)
                    enddo
                enddo
            enddo
        enddo

        if (flag_BCs == 1) then
            if (npy == 0) then
                Helm_Ax_nth(:,1)=v_in(:,1)*grid%hxhy
            endif
            if (npy == npy0-1) then
                Helm_Ax_nth(:,grid%ny)=v_in(:,grid%ny)*grid%hxhy
            endif
            if (npx == 0) then
                Helm_Ax_nth(1,:)=v_in(1,:)*grid%hxhy
            endif
            if (npx==npx0-1) then
                Helm_Ax_nth(grid%nx,:)=v_in(grid%nx,:)*grid%hxhy
            endif 

        elseif (flag_BCs ==2) then
        
            if (npy == 0) then
                j = 1
                do i=1,grid%nx
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    Helm_Ax_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npy == npy0-1) then
                j = grid%ny
                do i=1,grid%nx
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    Helm_Ax_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npx == 0) then
                i = 1
                do j=1,grid%ny
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    Helm_Ax_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npx==npx0-1) then
                i = grid%nx
                do j=1,grid%ny
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    Helm_Ax_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif
        endif

        Helm_Ax_nth=Helm_Ax_nth/grid%hxhy

    end function Helm_Ax_nth

    function CSLP_Mx_nth(v_in,grid) 
        !! *ReD-Glk* CSLP operator for coarse levels
        implicit none

        type(Gridpara), intent(inout) :: grid
        complex(kind = realdp), intent(inout) :: v_in(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)
        complex(kind=realdp) :: CSLP_Mx_nth(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)

        real(kind = realdp) :: ae,aw,an,as
        complex(kind = realdp) :: ap, coefA(-3:3,-3:3)
        integer :: i,j,ia,ja
        integer(kind=4) :: level_i
        
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1

        CSLP_Mx_nth   = czero
        v_in(1-LAP:0,:) = czero
        v_in(:,1-LAP:0) = czero
        v_in(grid%nx+1:grid%nx+LAP,:) = czero
        v_in(:,grid%ny+1:grid%ny+LAP) = czero

        call mg_check_xy2d(v_in, grid%nx, grid%ny)
        call GridBase_ExtrpltGhostBCs(v_in,grid)

        do j=1,grid%ny
            do i=1,grid%nx
                call CSLP_Mx_nth_stencils(coefA,grid,i,j)
                do ja=-3,3
                    do ia=-3,3
                        CSLP_Mx_nth(i,j)=CSLP_Mx_nth(i,j)+coefA(ia,ja)*v_in(i+ia,j+ja)
                    enddo
                enddo
            enddo
        enddo

        if (flag_BCs == 1) then

            if (npy == 0) then
                CSLP_Mx_nth(:,1)=v_in(:,1)*grid%hxhy
            endif
            if (npy == npy0-1) then
                CSLP_Mx_nth(:,grid%ny)=v_in(:,grid%ny)*grid%hxhy
            endif
            if (npx == 0) then
                CSLP_Mx_nth(1,:)=v_in(1,:)*grid%hxhy
            endif
            if (npx==npx0-1) then
                CSLP_Mx_nth(grid%nx,:)=v_in(grid%nx,:)*grid%hxhy
            endif 

        elseif (flag_BCs ==2) then
        
            if (npy == 0) then
                j = 1
                do i=1,grid%nx
                    call cslp2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    CSLP_Mx_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npy == npy0-1) then
                j = grid%ny
                do i=1,grid%nx
                    call cslp2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    CSLP_Mx_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npx == 0) then
                i = 1
                do j=1,grid%ny
                    call cslp2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    CSLP_Mx_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npx==npx0-1) then
                i = grid%nx
                do j=1,grid%ny
                    call cslp2d_stencils(ap,an,as,aw,ae,grid%wavenumber_kh_pow_2(i,j))
                    !==============================
                    CSLP_Mx_nth(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif
        endif

        CSLP_Mx_nth=CSLP_Mx_nth/grid%hxhy

    end function CSLP_Mx_nth

    subroutine GridBase_ExtrpltGhostBCs(v_in,grid)
        !! Extrapolation of a layer of ghost grid points
        implicit none

        type(Gridpara), intent(inout) :: grid
        complex(kind = realdp), intent(inout)  :: v_in(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)
        integer :: i,j
        
        if (flag_BCs == 1) then
            if (npx == 0) then
                do j=1-LAP,grid%ny+LAP
                    v_in(0,j)=2.d0*v_in(1,j)-v_in(2,j)
                    grid%wavenumber_k(0,j)=2.d0*grid%wavenumber_k(1,j)-grid%wavenumber_k(2,j)
                    grid%wavenumber_kh_pow_2(0,j)=grid%wavenumber_k(0,j)*grid%hxhy
                enddo
            endif
            if (npx==npx0-1) then
                do j=1-LAP,grid%ny+LAP
                    v_in(grid%nx+1,j)=2.d0*v_in(grid%nx,j)-v_in(grid%nx-1,j)
                    grid%wavenumber_k(grid%nx+1,j)=2.d0*grid%wavenumber_k(grid%nx,j)-grid%wavenumber_k(grid%nx-1,j)
                    grid%wavenumber_kh_pow_2(grid%nx+1,j)=grid%wavenumber_k(grid%nx+1,j)*grid%hxhy
                enddo
            endif
            if (npy == 0) then
                do i=1-LAP,grid%nx+LAP
                    v_in(i,0)=2.d0*v_in(i,1)-v_in(i,2)
                    grid%wavenumber_k(i,0)=2.d0*grid%wavenumber_k(i,1)-grid%wavenumber_k(i,2)
                    grid%wavenumber_kh_pow_2(i,0)=grid%wavenumber_k(i,0)*grid%hxhy
                enddo
            endif
            if (npy == npy0-1) then
                do i=1-LAP,grid%nx+LAP
                    v_in(i,grid%ny+1)=2.d0*v_in(i,grid%ny)-v_in(i,grid%ny-1)
                    grid%wavenumber_k(i,grid%ny+1)=2.d0*grid%wavenumber_k(i,grid%ny)-grid%wavenumber_k(i,grid%ny-1)
                    grid%wavenumber_kh_pow_2(i,grid%ny+1)=grid%wavenumber_k(i,grid%ny+1)*grid%hxhy
                enddo
            endif

            !=THE CORNER=
            if (npx == 0 .and. npy == 0) then
                v_in(0,0)=4.d0*v_in(1,1)-v_in(2,2)-v_in(0,2)-v_in(2,0)
                grid%wavenumber_k(0,0)=4.d0*grid%wavenumber_k(1,1)-grid%wavenumber_k(2,2)-grid%wavenumber_k(0,2)-grid%wavenumber_k(2,0)
                grid%wavenumber_kh_pow_2(0,0)=grid%wavenumber_k(0,0)*grid%hxhy
            endif

            if (npx == 0 .and. npy == npy0-1) then
                v_in(0,grid%ny+1)=4.d0*v_in(1,grid%ny)-v_in(2,grid%ny-1)-v_in(2,grid%ny+1)-v_in(0,grid%ny-1)
                grid%wavenumber_k(0,grid%ny+1)=4.d0*grid%wavenumber_k(1,grid%ny)-grid%wavenumber_k(2,grid%ny-1) &
                                           -grid%wavenumber_k(2,grid%ny+1)-grid%wavenumber_k(0,grid%ny-1)
                grid%wavenumber_kh_pow_2(0,grid%ny+1)=grid%wavenumber_k(0,grid%ny+1)*grid%hxhy
            endif

            if (npx==npx0-1 .and. npy == 0) then
                v_in(grid%nx+1,0)=4.d0*v_in(grid%nx,1)-v_in(grid%nx-1,2)-v_in(grid%nx-1,0)-v_in(grid%nx+1,2)
                grid%wavenumber_k(grid%nx+1,0)=4.d0*grid%wavenumber_k(grid%nx,1)-grid%wavenumber_k(grid%nx-1,2) &
                                           -grid%wavenumber_k(grid%nx-1,0)-grid%wavenumber_k(grid%nx+1,2)
                grid%wavenumber_kh_pow_2(grid%nx+1,0)=grid%wavenumber_k(grid%nx+1,0)*grid%hxhy
            endif

            if (npx==npx0-1 .and. npy == npy0-1) then
                v_in(grid%nx+1,grid%ny+1)=4.d0*v_in(grid%nx,grid%ny)-v_in(grid%nx-1,grid%ny-1) &
                                              -v_in(grid%nx-1,grid%ny+1)-v_in(grid%nx+1,grid%ny-1)
                grid%wavenumber_k(grid%nx+1,grid%ny+1)=4.d0*grid%wavenumber_k(grid%nx,grid%ny)-grid%wavenumber_k(grid%nx-1,grid%ny-1) &
                                               -grid%wavenumber_k(grid%nx-1,grid%ny+1)-grid%wavenumber_k(grid%nx+1,grid%ny-1)
                grid%wavenumber_kh_pow_2(grid%nx+1,grid%ny+1)=grid%wavenumber_k(grid%nx+1,grid%ny+1)*grid%hxhy
            endif

        elseif (flag_BCs == 2) then
            if (npx == 0) then
                do j=1-LAP,grid%ny+LAP
                    v_in(0,j)=2.d0*cone*grid%wavenumber_k(1,j)*grid%hx*v_in(1,j)+v_in(2,j)
                    grid%wavenumber_k(0,j)=0.d0
                    grid%wavenumber_kh_pow_2(0,j)=grid%wavenumber_k(0,j)*grid%hxhy
                enddo
            endif
            if (npx==npx0-1) then
                do j=1-LAP,grid%ny+LAP
                    v_in(grid%nx+1,j)=2.d0*cone*grid%wavenumber_k(grid%nx,j)*grid%hx* &
                                    v_in(grid%nx,j)+v_in(grid%nx-1,j)
                    grid%wavenumber_k(grid%nx+1,j)=0.d0
                    grid%wavenumber_kh_pow_2(grid%nx+1,j)=grid%wavenumber_k(grid%nx+1,j)*grid%hxhy
                enddo
            endif
            if (npy == 0) then
                do i=1-LAP,grid%nx+LAP
                    v_in(i,0)=2.d0*cone*grid%wavenumber_k(i,1)*grid%hy*v_in(i,1)+v_in(i,2)
                    grid%wavenumber_k(i,0)=0.d0
                    grid%wavenumber_kh_pow_2(i,0)=grid%wavenumber_k(i,0)*grid%hxhy
                    ! !!Another Ghost point seems no gains
                enddo
            endif
            if (npy == npy0-1) then
                do i=1-LAP,grid%nx+LAP
                    v_in(i,grid%ny+1)=2.d0*cone*grid%wavenumber_k(i,grid%ny)*grid%hy* &
                                    v_in(i,grid%ny)+v_in(i,grid%ny-1)
                    grid%wavenumber_k(i,grid%ny+1)=0.d0
                    grid%wavenumber_kh_pow_2(i,grid%ny+1)=grid%wavenumber_k(i,grid%ny+1)*grid%hxhy
                enddo
            endif

            !=THE CORNER=
            if (npx == 0 .and. npy == 0) then
            !   (-1, 3)--(0, 3)-|-(1, 3)--(2, 3)--(3, 3)
            !      |        |   |    |       |       |
            !   (-1, 2)--(0, 2)-|-(1, 2)--(2, 2)--(3, 2)
            !      |        |   |    |       |       |
            !   (-1, 1)--(0, 1)-|-(1, 1)--(2, 1)--(3, 1)
            !    __|________|___|____|_______|_______|__
            !   (-1, 0)--(0, 0)-|-(1, 0)--(2, 0)--(3, 0)
            !      |        |   |    |       |       |
            !   (-1,-1)--(0,-1)-|-(1,-1)--(2,-1)--(3,-1)
                v_in(0,0)=((2.d0*cone*grid%wavenumber_k(1,0)*grid%hx*v_in(1,0)+v_in(2,0)) &
                          +(2.d0*cone*grid%wavenumber_k(0,1)*grid%hy*v_in(0,1)+v_in(0,2)))/2.d0
                grid%wavenumber_k(0,0)=0.d0
                grid%wavenumber_kh_pow_2(0,0)=grid%wavenumber_k(0,0)*grid%hxhy
            endif

            if (npx == 0 .and. npy == npy0-1) then
            !   (-1,n+2)--(0,n+2)-|-(1,n+2)--(2,n+2)--(3,n+2)
            !      |        |     |    |       |       |
            !   (-1,n+1)--(0,n+1)-|-(1,n+1)--(2,n+1)--(3,n+1)
            !    __|________|_____|____|_______|_______|____
            !   (-1, n )--(0, n )-|-(1, n )--(2, n )--(3, n )
            !      |        |     |    |       |       |
            !   (-1,n-1)--(0,n-1)-|-(1,n-1)--(2,n-1)--(3,n-1)
            !      |        |     |    |       |       |
            !   (-1,n-2)--(0,n-2)-|-(1,n-2)--(2,n-2)--(3,n-2)
                v_in(0,grid%ny+1)=((2.d0*cone*grid%wavenumber_k(1,grid%ny+1)*grid%hx*v_in(1,grid%ny+1)+v_in(2,grid%ny+1)) &
                                  +(2.d0*cone*grid%wavenumber_k(0,grid%ny)*grid%hy*v_in(0,grid%ny)+v_in(0,grid%ny-1)))/2.d0
                grid%wavenumber_k(0,grid%ny+1)=0.d0
                grid%wavenumber_kh_pow_2(0,grid%ny+1)=grid%wavenumber_k(0,grid%ny+1)*grid%hxhy
            endif

            if (npx==npx0-1 .and. npy == 0) then
            !   (n-2, 3)--(n-1, 3)--( n , 3)--|--(n+1, 3)--(n+2, 3)
            !       |        |          |     |      |       |
            !   (n-2, 2)--(n-1, 2)--( n , 2)--|--(n+1, 2)--(n+2, 2)
            !       |        |          |     |      |       |
            !   (n-2, 1)--(n-1, 1)--( n , 1)--|--(n+1, 1)--(n+2, 1)
            !   ____|________|__________|_____|______|_______|______
            !       |        |          |     |      |       |
            !   (n-2, 0)--(n-1, 0)--( n , 0)--|--(n+1, 0)--(n+2, 0)
            !       |        |          |     |      |       |
            !   (n-2,-1)--(n-1,-1)--( n ,-1)--|--(n+1,-1)--(n+2,-1)

                v_in(grid%nx+1,0)=((2.d0*cone*grid%wavenumber_k(grid%nx,0)*grid%hx*v_in(grid%nx,0)+v_in(grid%nx-1,0)) &
                                  +(2.d0*cone*grid%wavenumber_k(grid%nx+1,1)*grid%hy*v_in(grid%nx+1,1)+v_in(grid%nx+1,2)))/2.d0
                grid%wavenumber_k(grid%nx+1,0)=0.d0
                grid%wavenumber_kh_pow_2(grid%nx+1,0)=grid%wavenumber_k(grid%nx+1,0)*grid%hxhy
            endif

            if (npx==npx0-1 .and. npy == npy0-1) then
            !   (n-2,n+2)--(n-1,n+2)--( n ,n+2)--|--(n+1,n+2)--(n+2,n+2)
            !       |          |          |      |      |          |
            !   (n-2,n+1)--(n-1,n+1)--( n ,n+1)--|--(n+1,n+1)--(n+2,n+1)
            !   ____|__________|__________|______|______|__________|______
            !       |          |          |      |      |          |
            !   (n-2, n )--(n-1, n )--( n , n )--|--(n+1, n )--(n+2, n )
            !       |          |          |      |      |          |
            !   (n-2,n-1)--(n-1,n-1)--( n ,n-1)--|--(n+1,n-1)--(n+2,n-1)
            !       |          |          |      |      |          |
            !   (n-2,n-2)--(n-1,n-2)--( n ,n-2)--|--(n+1,n-2)--(n+2,n-2)

                v_in(grid%nx+1,grid%ny+1)=((2.d0*cone*grid%wavenumber_k(grid%nx,grid%ny+1)*grid%hx*v_in(grid%nx,grid%ny+1)&
                                        +v_in(grid%nx-1,grid%ny+1)) &
                                        +(2.d0*cone*grid%wavenumber_k(grid%nx+1,grid%ny)*grid%hy*v_in(grid%nx+1,grid%ny)&
                                        +v_in(grid%nx+1,grid%ny-1)))/2.d0
                grid%wavenumber_k(grid%nx+1,grid%ny+1)=0.d0
                grid%wavenumber_kh_pow_2(grid%nx+1,grid%ny+1)=grid%wavenumber_k(grid%nx+1,grid%ny+1)*grid%hxhy
            endif
        else
            write(*,*) "There is no such a boundary conditions yet!!!"
            stop
        endif
        
    end subroutine GridBase_ExtrpltGhostBCs

    subroutine Helm_Ax_nth_stencils(Stcl_nth,grid,ic,jc)
        !! *ReD-Glk* computational stencils of the Helmholtz operator for differnet coarse-grid level
        implicit none

        type(Gridpara), intent(inout) :: grid
            !! parameters for the current coarse-grid system
        integer :: ic,jc
            !! Index for the wavenumber field

        integer :: ib,jb
        complex(kind = realdp) :: a(0:3,0:3),b(0:3,0:3)
            !! a: stencils for Laplace operator; b: stencils for wavenumber operator
        complex(kind = realdp) :: Stcl_nth(-3:3,-3:3)
            !! Resulting stencils for the Helmholtz operator
        integer(kind=4) :: level_i
            !! The current #level_i level
        
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1 
        Stcl_nth = czero
        if (level_i == 2) then
            a(:,0) =[complex(realdp) ::  3.82812500000d0,  0.21875000000d0, -0.38281250000d0, 0.d0 ]
            a(:,1) =[complex(realdp) ::  0.21875000000d0, -0.43750000000d0, -0.17187500000d0, 0.d0 ]
            a(:,2) =[complex(realdp) :: -0.38281250000d0, -0.17187500000d0, -0.01171875000d0, 0.d0 ]
            a(:,3) =[complex(realdp) ::  0.00000000000d0,  0.00000000000d0,  0.00000000000d0, 0.d0 ]

            b(:,0) =[complex(realdp) :: 1.196289062500d0, 0.478515625000d0, 0.017089843750d0, 0.d0 ]
            b(:,1) =[complex(realdp) :: 0.478515625000d0, 0.191406250000d0, 0.006835937500d0, 0.d0 ]
            b(:,2) =[complex(realdp) :: 0.017089843750d0, 0.006835937500d0, 0.000244140625d0, 0.d0 ]
            b(:,3) =[complex(realdp) :: 0.000000000000d0, 0.000000000000d0, 0.000000000000d0, 0.d0 ]
        elseif (level_i == 3) then
            a(:,0) =[complex(realdp) ::  11.2361450195313d0,   1.46577453613281d0,  -1.35073852539063d0,  -0.0456085205078125d0 ]
            a(:,1) =[complex(realdp) ::  1.46577453613281d0,  -1.12293624877930d0, -0.791053771972656d0,  -0.0220222473144531d0 ]
            a(:,2) =[complex(realdp) :: -1.35073852539063d0, -0.791053771972656d0, -0.125289916992188d0, -0.00203704833984375d0 ]
            a(:,3) =[complex(realdp) :: -0.0456085205078125d0, -0.0220222473144531d0, -0.00203704833984375d0, -1.14440917968750d-5 ]

            b(:,0) =[complex(realdp) :: 3.90293979644775d0, 1.84391236305237d0, 0.155307292938232d0, 0.000482320785522461d0 ]
            b(:,1) =[complex(realdp) :: 1.84391236305237d0, 0.871141493320465d0, 0.0733736753463745d0, 0.000227868556976318d0 ]
            b(:,2) =[complex(realdp) :: 0.155307292938232d0, 0.0733736753463745d0, 0.00618004798889160d0, 1.91926956176758d-5 ]
            b(:,3) =[complex(realdp) :: 0.000482320785522461d0, 0.000227868556976318d0, 1.91926956176758d-5, 5.96046447753906d-8 ]
        elseif (level_i == 4) then
            a(:,0) =[complex(realdp) ::  41.8592979907990d0, 6.16552120447159d0, -5.19245016574860d0, -0.230845034122467d0 ]
            a(:,1) =[complex(realdp) ::  6.16552120447159d0, -3.95059673488140d0, -3.20627358555794d0, -0.117296531796455d0 ]
            a(:,2) =[complex(realdp) :: -5.19245016574860d0, -3.20627358555794d0, -0.582719027996063d0, -0.0132198035717011d0 ]
            a(:,3) =[complex(realdp) :: -0.230845034122467d0, -0.117296531796455d0, -0.0132198035717011d0, -0.000154897570610046d0 ]

            b(:,0) =[complex(realdp) :: 14.9228199282661d0, 7.28212934662588d0, 0.703624110203236d0, 0.00486294622533023d0 ]
            b(:,1) =[complex(realdp) :: 7.28212934662588d0, 3.55357821617508d0, 0.343358815996908d0, 0.00237305037444457d0 ]
            b(:,2) =[complex(realdp) :: 0.703624110203236d0, 0.343358815996908d0, 0.0331764968577772d0, 0.000229292199946940d0 ]
            b(:,3) =[complex(realdp) :: 0.00486294622533023d0, 0.00237305037444457d0, 0.000229292199946940d0, 1.58470356836915d-6 ]
        elseif (level_i == 5) then
            a(:,0) =[complex(realdp) ::  164.563205060549d0, 24.9116344174836d0, -20.5528745003976d0, -0.972393697360531d0 ]
            a(:,1) =[complex(realdp) ::  24.9116344174836d0, -15.2913385086576d0, -12.8529095315607d0, -0.500045731023420d0 ]
            a(:,2) =[complex(realdp) :: -20.5528745003976d0, -12.8529095315607d0, -2.41139203752391d0, -0.0588705557165667d0 ]
            a(:,3) =[complex(realdp) :: -0.972393697360531d0, -0.500045731023420d0, -0.0588705557165667d0, -0.000785302079748362d0 ]

            b(:,0) =[complex(realdp) :: 59.0402948639394d0, 29.0317872028857d0, 2.89513443097121d0, 0.0230771133726648d0 ]
            b(:,1) =[complex(realdp) :: 29.0317872028857d0, 14.2757530282666d0, 1.42361969765568d0, 0.0113476710479858d0 ]
            b(:,2) =[complex(realdp) :: 2.89513443097121d0, 1.42361969765568d0, 0.141967505289585d0, 0.00113162282889334d0 ]
            b(:,3) =[complex(realdp) :: 0.0230771133726648d0, 0.0113476710479858d0, 0.00113162282889334d0, 9.02016432746677d-6 ]
        elseif (level_i == 6) then
            a(:,0) =[complex(realdp) ::  655.429688162392d0, 99.8838619755779d0, -81.9929805306983d0, -3.93873333857573d0 ]
            a(:,1) =[complex(realdp) ::  99.8838619755779d0, -60.6622323345425d0, -51.4359953570679d0, -2.03144743680355d0 ]
            a(:,2) =[complex(realdp) :: -81.9929805306983d0, -51.4359953570679d0, -9.72570502689905d0, -0.241711694433889d0 ]
            a(:,3) =[complex(realdp) :: -3.93873333857573d0, -2.03144743680355d0, -0.241711694433889d0, -0.00332380884970007d0 ]

            b(:,0) =[complex(realdp) :: 235.519105626998d0, 116.029743978018d0, 11.6606996686767d0, 0.0961113088190835d0 ]
            b(:,1) =[complex(realdp) :: 116.029743978018d0, 57.1626724369705d0, 5.74470590638153d0, 0.0473497490829131d0 ]
            b(:,2) =[complex(realdp) :: 11.6606996686767d0, 5.74470590638153d0, 0.577328605257279d0, 0.00475853160158357d0 ]
            b(:,3) =[complex(realdp) :: 0.0961113088190835d0, 0.0473497490829131d0, 0.00475853160158357d0, 3.92213772140715d-5 ]
        else !! Default 2nd-ordr stencils
            a(:,0) =[complex(realdp) ::   4.d0, -1.d0, 0.d0, 0.d0 ]
            a(:,1) =[complex(realdp) ::  -1.d0,  0.d0, 0.d0, 0.d0 ]
            a(:,2) =[complex(realdp) ::   0.d0,  0.d0, 0.d0, 0.d0 ]
            a(:,3) =[complex(realdp) ::   0.d0,  0.d0, 0.d0, 0.d0 ]

            b(:,0) =[complex(realdp) ::  1.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,1) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,2) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,3) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
        endif

        do jb=-3,3
            do ib=-3,3
                Stcl_nth(ib,jb) = a(ABS(ib),ABS(jb))-b(ABS(ib),ABS(jb))*grid%wavenumber_kh_pow_2(ic+ib,jc+jb)
            enddo
        enddo

    end subroutine Helm_Ax_nth_stencils

    subroutine CSLP_Mx_nth_stencils(Stcl_nth,grid,ic,jc)
        !! *ReD-Glk* computational stencils of the CSLP operator for differnet coarse-grid level
        implicit none

        type(Gridpara), intent(inout) :: grid
            !! parameters for the current coarse-grid system
        integer :: ic,jc
            !! Index for the wavenumber field

        integer :: ib,jb
        complex(kind = realdp) :: a(0:3,0:3),b(0:3,0:3)
            !! a: stencils for Laplace operator; b: stencils for wavenumber operator
        complex(kind = realdp) :: Stcl_nth(-3:3,-3:3)
            !! Resulting stencils for the CSLP operator
        integer(kind=4) :: level_i
            !! The current #level_i level
        
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1
        Stcl_nth = czero
        if (level_i == 2) then
            a(:,0) =[complex(realdp) ::  3.82812500000d0,  0.21875000000d0, -0.38281250000d0, 0.d0 ]
            a(:,1) =[complex(realdp) ::  0.21875000000d0, -0.43750000000d0, -0.17187500000d0, 0.d0 ]
            a(:,2) =[complex(realdp) :: -0.38281250000d0, -0.17187500000d0, -0.01171875000d0, 0.d0 ]
            a(:,3) =[complex(realdp) ::  0.00000000000d0,  0.00000000000d0,  0.00000000000d0, 0.d0 ]

            b(:,0) =[complex(realdp) :: 1.196289062500d0, 0.478515625000d0, 0.017089843750d0, 0.d0 ]
            b(:,1) =[complex(realdp) :: 0.478515625000d0, 0.191406250000d0, 0.006835937500d0, 0.d0 ]
            b(:,2) =[complex(realdp) :: 0.017089843750d0, 0.006835937500d0, 0.000244140625d0, 0.d0 ]
            b(:,3) =[complex(realdp) :: 0.000000000000d0, 0.000000000000d0, 0.000000000000d0, 0.d0 ]
        elseif (level_i == 3) then
            a(:,0) =[complex(realdp) ::  11.2361450195313d0,   1.46577453613281d0,  -1.35073852539063d0,  -0.0456085205078125d0 ]
            a(:,1) =[complex(realdp) ::  1.46577453613281d0,  -1.12293624877930d0, -0.791053771972656d0,  -0.0220222473144531d0 ]
            a(:,2) =[complex(realdp) :: -1.35073852539063d0, -0.791053771972656d0, -0.125289916992188d0, -0.00203704833984375d0 ]
            a(:,3) =[complex(realdp) :: -0.0456085205078125d0, -0.0220222473144531d0, -0.00203704833984375d0, -1.14440917968750d-5 ]

            b(:,0) =[complex(realdp) :: 3.90293979644775d0, 1.84391236305237d0, 0.155307292938232d0, 0.000482320785522461d0 ]
            b(:,1) =[complex(realdp) :: 1.84391236305237d0, 0.871141493320465d0, 0.0733736753463745d0, 0.000227868556976318d0 ]
            b(:,2) =[complex(realdp) :: 0.155307292938232d0, 0.0733736753463745d0, 0.00618004798889160d0, 1.91926956176758d-5 ]
            b(:,3) =[complex(realdp) :: 0.000482320785522461d0, 0.000227868556976318d0, 1.91926956176758d-5, 5.96046447753906d-8 ]
        elseif (level_i == 4) then
            a(:,0) =[complex(realdp) ::  41.8592979907990d0, 6.16552120447159d0, -5.19245016574860d0, -0.230845034122467d0 ]
            a(:,1) =[complex(realdp) ::  6.16552120447159d0, -3.95059673488140d0, -3.20627358555794d0, -0.117296531796455d0 ]
            a(:,2) =[complex(realdp) :: -5.19245016574860d0, -3.20627358555794d0, -0.582719027996063d0, -0.0132198035717011d0 ]
            a(:,3) =[complex(realdp) :: -0.230845034122467d0, -0.117296531796455d0, -0.0132198035717011d0, -0.000154897570610046d0 ]

            b(:,0) =[complex(realdp) :: 14.9228199282661d0, 7.28212934662588d0, 0.703624110203236d0, 0.00486294622533023d0 ]
            b(:,1) =[complex(realdp) :: 7.28212934662588d0, 3.55357821617508d0, 0.343358815996908d0, 0.00237305037444457d0 ]
            b(:,2) =[complex(realdp) :: 0.703624110203236d0, 0.343358815996908d0, 0.0331764968577772d0, 0.000229292199946940d0 ]
            b(:,3) =[complex(realdp) :: 0.00486294622533023d0, 0.00237305037444457d0, 0.000229292199946940d0, 1.58470356836915d-6 ]
        elseif (level_i == 5) then
            a(:,0) =[complex(realdp) ::  164.563205060549d0, 24.9116344174836d0, -20.5528745003976d0, -0.972393697360531d0 ]
            a(:,1) =[complex(realdp) ::  24.9116344174836d0, -15.2913385086576d0, -12.8529095315607d0, -0.500045731023420d0 ]
            a(:,2) =[complex(realdp) :: -20.5528745003976d0, -12.8529095315607d0, -2.41139203752391d0, -0.0588705557165667d0 ]
            a(:,3) =[complex(realdp) :: -0.972393697360531d0, -0.500045731023420d0, -0.0588705557165667d0, -0.000785302079748362d0 ]

            b(:,0) =[complex(realdp) :: 59.0402948639394d0, 29.0317872028857d0, 2.89513443097121d0, 0.0230771133726648d0 ]
            b(:,1) =[complex(realdp) :: 29.0317872028857d0, 14.2757530282666d0, 1.42361969765568d0, 0.0113476710479858d0 ]
            b(:,2) =[complex(realdp) :: 2.89513443097121d0, 1.42361969765568d0, 0.141967505289585d0, 0.00113162282889334d0 ]
            b(:,3) =[complex(realdp) :: 0.0230771133726648d0, 0.0113476710479858d0, 0.00113162282889334d0, 9.02016432746677d-6 ]
        elseif (level_i == 6) then
            a(:,0) =[complex(realdp) ::  655.429688162392d0, 99.8838619755779d0, -81.9929805306983d0, -3.93873333857573d0 ]
            a(:,1) =[complex(realdp) ::  99.8838619755779d0, -60.6622323345425d0, -51.4359953570679d0, -2.03144743680355d0 ]
            a(:,2) =[complex(realdp) :: -81.9929805306983d0, -51.4359953570679d0, -9.72570502689905d0, -0.241711694433889d0 ]
            a(:,3) =[complex(realdp) :: -3.93873333857573d0, -2.03144743680355d0, -0.241711694433889d0, -0.00332380884970007d0 ]

            b(:,0) =[complex(realdp) :: 235.519105626998d0, 116.029743978018d0, 11.6606996686767d0, 0.0961113088190835d0 ]
            b(:,1) =[complex(realdp) :: 116.029743978018d0, 57.1626724369705d0, 5.74470590638153d0, 0.0473497490829131d0 ]
            b(:,2) =[complex(realdp) :: 11.6606996686767d0, 5.74470590638153d0, 0.577328605257279d0, 0.00475853160158357d0 ]
            b(:,3) =[complex(realdp) :: 0.0961113088190835d0, 0.0473497490829131d0, 0.00475853160158357d0, 3.92213772140715d-5 ]
        else
            a(:,0) =[complex(realdp) ::   4.d0, -1.d0, 0.d0, 0.d0 ]
            a(:,1) =[complex(realdp) ::  -1.d0,  0.d0, 0.d0, 0.d0 ]
            a(:,2) =[complex(realdp) ::   0.d0,  0.d0, 0.d0, 0.d0 ]
            a(:,3) =[complex(realdp) ::   0.d0,  0.d0, 0.d0, 0.d0 ]

            b(:,0) =[complex(realdp) ::  1.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,1) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,2) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
            b(:,3) =[complex(realdp) ::  0.d0,  0.d0, 0.d0, 0.d0 ]
        endif

        do jb=-3,3
            do ib=-3,3
                Stcl_nth(ib,jb) = a(ABS(ib),ABS(jb))-(beta1 - beta2*cone)*b(ABS(ib),ABS(jb))*grid%wavenumber_kh_pow_2(ic+ib,jc+jb)
            enddo
        enddo

    end subroutine CSLP_Mx_nth_stencils

end module operators