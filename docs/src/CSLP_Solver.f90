module CSLP_Solver
    !! This is a module to solve the inverse of CSLP approximately
    !! By multigrid methods or Krylov iterations
    use mpi
    use comm_variable
    use mpi_setup
    use wavenumber, only:wavenumber_k,kh_pow_2
    use operators
    use smoother
    implicit none

    type, public    :: GridSystem
        !! A data type that collect some info of a (coarse) grid system 
        !! It is different the type Gridpara, which does not consist any variable like u, rhs, res, and etc.
        integer(kind=4)     :: nxc_global, nyc_global
            !! Global grid size, in x and y directions respectively.
        integer(kind=4)     :: nxc, nyc
            !! Local grid size of a subdomain, in x and y directions respectively.
        integer, dimension(0:npMax - 1) :: ic_offset  
            !! ic_offset: an array that contains the index offset for MPI ranks in x direction   
        integer, dimension(0:npMax - 1) :: jc_offset  
            !! jc_offset: an array that contains the index offset for MPI ranks in y direction  
        integer, dimension(0:npMax - 1) :: ic_nn  
            !! ic_nn: an array that contains the number of grid points for MPI ranks in x direction   
        integer, dimension(0:npMax - 1) :: jc_nn  
            !! jc_nn: an array that contains the number of grid points for MPI ranks in y direction
        real(kind = realdp)        :: hxc, hyc, hxhyc
            !! space step, hxhyc==hx*hy==hx^2

        complex(kind = realdp), allocatable, dimension(:,:)     :: u_c    
            !! Solution variable in subdomain
        complex(kind = realdp), allocatable, dimension(:,:)     :: rhs_c    
            !! RHS variable in subdomain
        complex(kind = realdp), allocatable, dimension(:,:)     :: res_c    
            !! Residual variable in subdomain
        real(kind = realdp),    allocatable, dimension(:,:)     :: kxy_c,kh2_c    
            !! wavenumber in subdomain
    end type GridSystem

    contains

    function MGCSLP_invMx(b_in)  
        !! Multigrid-Based CSLP, ONLY starts from the default (finest) grid sysyem
        implicit none
        
        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: b_in
        complex(kind = realdp) :: MGCSLP_invMx(1-LAP:nx+LAP,1-LAP:ny+LAP)
            !! MGCSLP_invMx = M_h^(-1)(b_h)

        type(GridSystem) :: mg_global
    
        MGCSLP_invMx = czero

        call finestgrid_define(mg_global)
        mg_global%rhs_c=b_in

        if (MG_flag == 0 .or. MG_flag == 1) then
            !! One V-cycle
            if (mod(mg_global%nxc_global-1,2) == 0 .and. mod(mg_global%nyc_global-1,2) == 0) then
                call V_cycle(mg_global)
            endif
        elseif (MG_flag == 12) then
            !! Two V-cycle
            if (mod(mg_global%nxc_global-1,2) == 0 .and. mod(mg_global%nyc_global-1,2) == 0) then
                call V_cycle(mg_global)
                call V_cycle(mg_global)
            endif
        elseif (MG_flag == 3) then
            !! One F-cycle
            if (mod(mg_global%nxc_global-1,2) == 0 .and. mod(mg_global%nyc_global-1,2) == 0) then
                call F_cycle(mg_global)
            endif
        else
            write(*,*) "No such a multigrid method yet!!"
            stop
        endif

        MGCSLP_invMx=mg_global%u_c
        call grid_destroy(mg_global)

        return
    end function MGCSLP_invMx

    function MGCSLP_invMHx(b_in, grid)  
        !! Multigrid-based CSLP, starts from the a specified grid system
        implicit none
        type(Gridpara) :: grid
            !! A data type that collect some parameters of a system
        complex(kind = realdp),dimension(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)  :: b_in
        complex(kind = realdp),dimension(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP) :: MGCSLP_invMHx
            !! MGCSLP_invMHx = M_H^(-1)(b_H)
        
        type(GridSystem) :: basegridsys
    
        MGCSLP_invMHx = czero

        call gridsys_define(grid, basegridsys)
        basegridsys%rhs_c=b_in

        if (MG_flag == 0 .or. MG_flag == 1) then
            !! One V-cycle
            if (mod(basegridsys%nxc_global-1,2) == 0 .and. mod(basegridsys%nyc_global-1,2) == 0) then
                call V_cycle(basegridsys)
            endif
        elseif (MG_flag == 12) then
            !! Two V-cycle
            if (mod(basegridsys%nxc_global-1,2) == 0 .and. mod(basegridsys%nyc_global-1,2) == 0) then
                call V_cycle(basegridsys)
                call V_cycle(basegridsys)
            endif
        elseif (MG_flag == 3) then
            !! One F-cycle
            if (mod(basegridsys%nxc_global-1,2) == 0 .and. mod(basegridsys%nyc_global-1,2) == 0) then
                call F_cycle(basegridsys)
            endif
        else
            write(*,*) "No such a multigrid method yet!!"
            stop
        endif

        MGCSLP_invMHx=basegridsys%u_c
        call grid_destroy(basegridsys)

        return
    end function MGCSLP_invMHx

    function KrylovCSLP_invMHx(b_in, grid)  
        !! Krylov-based CSLP, starts from the a specified grid sysyem
        implicit none
        type(Gridpara) :: grid
            !! A data type that collect some parameters of a system
        complex(kind = realdp),dimension(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)  :: b_in
        complex(kind = realdp),dimension(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)  :: KrylovCSLP_invMHx
        
        type(GridSystem) :: basegridsys
        integer(kind=4) :: C_it
            !! The maximum number of iterations for the Krylov iterations

        KrylovCSLP_invMHx = czero

        C_it = ceiling(6.d0*sqrt(sqrt(dble(grid%nx_global*grid%ny_global))))


        call gridsys_define(grid, basegridsys)
        basegridsys%rhs_c=b_in

        !! Choose Bi-CGSATB or GMRES by uncomment or comment 
        ! call ReD_Glk_CSLP_bicgstab(basegridsys, grid, C_it) 
        call ReD_Glk_CSLP_gmres(basegridsys, grid, C_it)

        KrylovCSLP_invMHx=basegridsys%u_c
        call grid_destroy(basegridsys)

        return
    end function KrylovCSLP_invMHx

    subroutine finestgrid_define(mg_finest)
        !! This is a routine that define the finest grid system from the default (finest) grid parameters
        implicit none
        type(GridSystem), intent(inout)  :: mg_finest

        mg_finest%nxc_global = nx_global
        mg_finest%nyc_global = ny_global
        mg_finest%hxc        = hx
        mg_finest%hyc        = hy
        mg_finest%hxhyc      = hxhy
        mg_finest%nxc        = nx
        mg_finest%nyc        = ny
        mg_finest%ic_offset  = i_offset
        mg_finest%jc_offset  = j_offset
        mg_finest%ic_nn      = i_nn
        mg_finest%jc_nn      = j_nn

        allocate(mg_finest%u_c(1-LAP:mg_finest%nxc+LAP,1-LAP:mg_finest%nyc+LAP))
        allocate(mg_finest%rhs_c(1-LAP:mg_finest%nxc+LAP,1-LAP:mg_finest%nyc+LAP))
        allocate(mg_finest%res_c(1-LAP:mg_finest%nxc+LAP,1-LAP:mg_finest%nyc+LAP))
        allocate(mg_finest%kxy_c(1-LAP:mg_finest%nxc+LAP,1-LAP:mg_finest%nyc+LAP))
        allocate(mg_finest%kh2_c(1-LAP:mg_finest%nxc+LAP,1-LAP:mg_finest%nyc+LAP))

        mg_finest%u_c   = (0.d0,0.d0)
        mg_finest%rhs_c = (0.d0,0.d0)
        mg_finest%res_c = (0.d0,0.d0)
        mg_finest%kxy_c = wavenumber_k
        mg_finest%kh2_c = kh_pow_2

    end subroutine finestgrid_define

    subroutine gridsys_define(grid, basegridsys)
        !! This is a routine that define a GridSystem type from the grid parameters of a data type Gridpara
        !! One can find the differnece of GridSystem and Gridpara 
        implicit none
        type(Gridpara)   :: grid
        type(GridSystem) :: basegridsys

        basegridsys%nxc_global = grid%nx_global
        basegridsys%nyc_global = grid%ny_global
        basegridsys%hxc        = grid%hx
        basegridsys%hyc        = grid%hy
        basegridsys%hxhyc      = grid%hxhy
        basegridsys%nxc        = grid%nx
        basegridsys%nyc        = grid%ny
        basegridsys%ic_offset  = grid%i_offset
        basegridsys%jc_offset  = grid%j_offset
        basegridsys%ic_nn      = grid%i_nn
        basegridsys%jc_nn      = grid%j_nn

        allocate(basegridsys%u_c(1-LAP:basegridsys%nxc+LAP,1-LAP:basegridsys%nyc+LAP))
        allocate(basegridsys%rhs_c(1-LAP:basegridsys%nxc+LAP,1-LAP:basegridsys%nyc+LAP))
        allocate(basegridsys%res_c(1-LAP:basegridsys%nxc+LAP,1-LAP:basegridsys%nyc+LAP))
        allocate(basegridsys%kxy_c(1-LAP:basegridsys%nxc+LAP,1-LAP:basegridsys%nyc+LAP))
        allocate(basegridsys%kh2_c(1-LAP:basegridsys%nxc+LAP,1-LAP:basegridsys%nyc+LAP))

        basegridsys%u_c   = (0.d0,0.d0)
        basegridsys%rhs_c = (0.d0,0.d0)
        basegridsys%res_c = (0.d0,0.d0)
        basegridsys%kxy_c = grid%wavenumber_k
        basegridsys%kh2_c = grid%wavenumber_kh_pow_2

    end subroutine gridsys_define

    subroutine coarsegrid_create(mg_coarse,mg_fine)
        !! This is a routine that define a coarse grid system mg_coarse from a fine grid system mg_fine
        implicit none

        type(GridSystem), intent(inout)  :: mg_coarse,mg_fine

        integer(kind=4) :: k 

        mg_coarse%nxc_global =(mg_fine%nxc_global-1)/2+1
        mg_coarse%nyc_global =(mg_fine%nyc_global-1)/2+1

        mg_coarse%hxc=slx/(mg_coarse%nxc_global-1)
        mg_coarse%hyc=sly/(mg_coarse%nyc_global-1)
        mg_coarse%hxhyc=mg_coarse%hxc*mg_coarse%hyc

        if (mod(mg_fine%ic_offset(npx),2) == 0) then
            if (mod(mg_fine%nxc,2) == 0) then
                mg_coarse%nxc=mg_fine%nxc / 2
            else
                mg_coarse%nxc=(mg_fine%nxc - 1) / 2
            endif
        else
            if (mod(mg_fine%nxc,2) == 0) then
                mg_coarse%nxc=mg_fine%nxc / 2
            else
                mg_coarse%nxc=(mg_fine%nxc + 1) / 2
            endif
        endif

        if (mod(mg_fine%jc_offset(npy),2) == 0) then
            if (mod(mg_fine%nyc,2) == 0) then
                mg_coarse%nyc=mg_fine%nyc / 2
            else
                mg_coarse%nyc=(mg_fine%nyc - 1) / 2
            endif
        else
            if (mod(mg_fine%nyc,2) == 0) then
                mg_coarse%nyc=mg_fine%nyc / 2
            else
                mg_coarse%nyc=(mg_fine%nyc + 1) / 2
            endif
        endif


        do k=0,npx0-1
            if (mod(mg_fine%ic_offset(k),2) == 0) then

                mg_coarse%ic_offset(k) = (mg_fine%ic_offset(k)+1+1) / 2

                if (mod(mg_fine%ic_nn(k),2) == 0) then
                    mg_coarse%ic_nn(k)=mg_fine%ic_nn(k) / 2
                else
                    mg_coarse%ic_nn(k)=(mg_fine%ic_nn(k) - 1) / 2
                endif

            else

                mg_coarse%ic_offset(k) = (mg_fine%ic_offset(k)+1) / 2

                if (mod(mg_fine%ic_nn(k),2) == 0) then
                    mg_coarse%ic_nn(k)=mg_fine%ic_nn(k) / 2
                else
                    mg_coarse%ic_nn(k)=(mg_fine%ic_nn(k) + 1) / 2
                endif

            endif
        enddo

        do k=0,npy0-1
            if (mod(mg_fine%jc_offset(k),2) == 0) then

                mg_coarse%jc_offset(k) = (mg_fine%jc_offset(k)+1+1) / 2

                if (mod(mg_fine%jc_nn(k),2) == 0) then
                    mg_coarse%jc_nn(k)=mg_fine%jc_nn(k) / 2
                else
                    mg_coarse%jc_nn(k)=(mg_fine%jc_nn(k) - 1) / 2
                endif
            else

                mg_coarse%jc_offset(k) = (mg_fine%jc_offset(k)+1) / 2

                if (mod(mg_fine%jc_nn(k),2) == 0) then
                    mg_coarse%jc_nn(k)=mg_fine%jc_nn(k) / 2
                else
                    mg_coarse%jc_nn(k)=(mg_fine%jc_nn(k) + 1) / 2
                endif
            endif
        enddo

        allocate(mg_coarse%u_c(1-LAP:mg_coarse%nxc+LAP,1-LAP:mg_coarse%nyc+LAP))
        allocate(mg_coarse%rhs_c(1-LAP:mg_coarse%nxc+LAP,1-LAP:mg_coarse%nyc+LAP))
        allocate(mg_coarse%res_c(1-LAP:mg_coarse%nxc+LAP,1-LAP:mg_coarse%nyc+LAP))
        allocate(mg_coarse%kxy_c(1-LAP:mg_coarse%nxc+LAP,1-LAP:mg_coarse%nyc+LAP))
        allocate(mg_coarse%kh2_c(1-LAP:mg_coarse%nxc+LAP,1-LAP:mg_coarse%nyc+LAP))

        mg_coarse%u_c   = (0.d0,0.d0)
        mg_coarse%rhs_c = (0.d0,0.d0)
        mg_coarse%res_c = (0.d0,0.d0)
        mg_coarse%kxy_c = 0.d0
        mg_coarse%kh2_c = 0.d0

    end subroutine coarsegrid_create

    subroutine grid_destroy(mg_dm)
        !! This is a routine that deallocate the arrays of a grid system
        implicit none

        type(GridSystem), intent(inout)  :: mg_dm

        deallocate(mg_dm%u_c, mg_dm%rhs_c, mg_dm%res_c, mg_dm%kxy_c, mg_dm%kh2_c)

    end subroutine grid_destroy

    subroutine restriction(mg_coarse,mg_fine)
        !! This is routine that perform full-weight resctriction from fine to coarse grid system, mainly
        !! residual (res) on fine grid --> right-hand side (rhs) on coarse grid,
        !! wavenumber on fine grid --> wavenumber on the coarse grid
        !! Based on the relationship of the index between the fine and coarse grid
        implicit none

        type(GridSystem), intent(inout)  :: mg_coarse,mg_fine

        integer :: ii,jj,iif,jjf,icH_global,jcH_global,ifh_global,jfh_global

        call mg_check_xy2d(mg_fine%res_c,mg_fine%nxc,mg_fine%nyc)
        call mg_checkreal_xy2d(mg_fine%kxy_c,mg_fine%nxc,mg_fine%nyc)

        do jj=1,mg_coarse%nyc
            do ii=1,mg_coarse%nxc
                icH_global=mg_coarse%ic_offset(npx)-1+ii
                jcH_global=mg_coarse%jc_offset(npy)-1+jj
                ifh_global=2*icH_global-1
                jfh_global=2*jcH_global-1
                iif=ifh_global-(mg_fine%ic_offset(npx)-1)
                jjf=jfh_global-(mg_fine%jc_offset(npy)-1)


                mg_coarse%rhs_c(ii,jj)=(4.d0*mg_fine%res_c(iif,jjf) &
                                       +2.d0*mg_fine%res_c(iif-1,jjf)+2.d0*mg_fine%res_c(iif+1,jjf) &
                                       +2.d0*mg_fine%res_c(iif,jjf-1)+2.d0*mg_fine%res_c(iif,jjf+1) &
                                       +1.d0*mg_fine%res_c(iif-1,jjf-1)+1.d0*mg_fine%res_c(iif-1,jjf+1) &
                                       +1.d0*mg_fine%res_c(iif+1,jjf-1)+1.d0*mg_fine%res_c(iif+1,jjf+1))/16.d0

                mg_coarse%kxy_c(ii,jj)=(4.d0*mg_fine%kxy_c(iif,jjf) &
                                       +2.d0*mg_fine%kxy_c(iif-1,jjf)+2.d0*mg_fine%kxy_c(iif+1,jjf) &
                                       +2.d0*mg_fine%kxy_c(iif,jjf-1)+2.d0*mg_fine%kxy_c(iif,jjf+1) &
                                       +1.d0*mg_fine%kxy_c(iif-1,jjf-1)+1.d0*mg_fine%kxy_c(iif-1,jjf+1) &
                                       +1.d0*mg_fine%kxy_c(iif+1,jjf-1)+1.d0*mg_fine%kxy_c(iif+1,jjf+1))/16.d0

                if (npx == 0 .and. ii==1) then
                    mg_coarse%rhs_c(ii,jj)=mg_fine%res_c(iif,jjf)
                    mg_coarse%kxy_c(ii,jj)=mg_fine%kxy_c(iif,jjf)
                endif

                if (npx == npx0-1 .and. ii==mg_coarse%nxc) then
                    mg_coarse%rhs_c(ii,jj)=mg_fine%res_c(iif,jjf)
                    mg_coarse%kxy_c(ii,jj)=mg_fine%kxy_c(iif,jjf)
                endif

                if (npy == 0 .and. jj==1) then
                    mg_coarse%rhs_c(ii,jj)=mg_fine%res_c(iif,jjf)
                    mg_coarse%kxy_c(ii,jj)=mg_fine%kxy_c(iif,jjf)
                endif

                if (npy == npy0-1 .and. jj==mg_coarse%nyc) then
                    mg_coarse%rhs_c(ii,jj)=mg_fine%res_c(iif,jjf)
                    mg_coarse%kxy_c(ii,jj)=mg_fine%kxy_c(iif,jjf)
                endif

            enddo
        enddo

        if (npx == 0 .and. npy == 0 ) then
            mg_coarse%rhs_c(1,1)=mg_fine%res_c(1,1)
            mg_coarse%kxy_c(1,1)=mg_fine%kxy_c(1,1)
        endif

        if (npx == npx0-1 .and. npy == 0) then
            mg_coarse%rhs_c(mg_coarse%nxc,1)=mg_fine%res_c(mg_fine%nxc,1)
            mg_coarse%kxy_c(mg_coarse%nxc,1)=mg_fine%kxy_c(mg_fine%nxc,1)
        endif

        if (npx == 0 .and. npy == npy0-1 ) then
            mg_coarse%rhs_c(1,mg_coarse%nyc)=mg_fine%res_c(1,mg_fine%nyc)
            mg_coarse%kxy_c(1,mg_coarse%nyc)=mg_fine%kxy_c(1,mg_fine%nyc)
        endif

        if (npx == npx0-1 .and. npy == npy0-1) then
            mg_coarse%rhs_c(mg_coarse%nxc,mg_coarse%nyc)=mg_fine%res_c(mg_fine%nxc,mg_fine%nyc)
            mg_coarse%kxy_c(mg_coarse%nxc,mg_coarse%nyc)=mg_fine%kxy_c(mg_fine%nxc,mg_fine%nyc)
        endif

        mg_coarse%kh2_c=mg_coarse%kxy_c*mg_coarse%kxy_c*mg_coarse%hxhyc

    end subroutine restriction

    subroutine prolongation_en_correct(mg_coarse,mg_fine)
        !! This is routine that perform bilinear interpolation and correction from coarse to fine grid system, mainly
        !! solution (u) on coarse grid --> correction (e_h) for solution (u) of fine grid
        !! Based on the relationship of the index between the fine and coarse grid
        implicit none

        type(GridSystem), intent(inout)  :: mg_coarse,mg_fine
        complex(kind = realdp),allocatable, dimension(:,:) :: e_h
            !! Correction
        integer :: ii,jj,iif,jjf,icH_global,jcH_global,ifh_global,jfh_global

        allocate(e_h(1-LAP:mg_fine%nxc+LAP,1-LAP:mg_fine%nyc+LAP))
        e_h=(0.d0,0.d0)

        call mg_check_xy2d(mg_coarse%u_c,mg_coarse%nxc,mg_coarse%nyc)

        do jj=1,mg_coarse%nyc
            do ii=1,mg_coarse%nxc
                icH_global=mg_coarse%ic_offset(npx)-1+ii
                jcH_global=mg_coarse%jc_offset(npy)-1+jj
                ifh_global=2*icH_global-1
                jfh_global=2*jcH_global-1
                iif=ifh_global-(mg_fine%ic_offset(npx)-1)
                jjf=jfh_global-(mg_fine%jc_offset(npy)-1)

                e_h(iif,jjf) = mg_coarse%u_c(ii,jj)

                e_h(iif,jjf+1) = (mg_coarse%u_c(ii,jj) + mg_coarse%u_c(ii,jj+1))/2.d0
                e_h(iif,jjf-1) = (mg_coarse%u_c(ii,jj) + mg_coarse%u_c(ii,jj-1))/2.d0
                e_h(iif+1,jjf) = (mg_coarse%u_c(ii,jj) + mg_coarse%u_c(ii+1,jj))/2.d0
                e_h(iif-1,jjf) = (mg_coarse%u_c(ii,jj) + mg_coarse%u_c(ii-1,jj))/2.d0

                e_h(iif+1,jjf+1) = (mg_coarse%u_c(ii+1,jj) + mg_coarse%u_c(ii,jj+1) &
                                 + mg_coarse%u_c(ii+1,jj+1) + mg_coarse%u_c(ii,jj))/4.d0
                e_h(iif+1,jjf-1) = (mg_coarse%u_c(ii+1,jj) + mg_coarse%u_c(ii,jj-1) &
                                 + mg_coarse%u_c(ii+1,jj-1) + mg_coarse%u_c(ii,jj))/4.d0
                e_h(iif-1,jjf+1) = (mg_coarse%u_c(ii-1,jj) + mg_coarse%u_c(ii,jj+1) &
                                 + mg_coarse%u_c(ii-1,jj+1) + mg_coarse%u_c(ii,jj))/4.d0
                e_h(iif-1,jjf-1) = (mg_coarse%u_c(ii-1,jj) + mg_coarse%u_c(ii,jj-1) &
                                 + mg_coarse%u_c(ii-1,jj-1) + mg_coarse%u_c(ii,jj))/4.d0

            enddo
        enddo

        mg_fine%u_c(1:mg_fine%nxc,1:mg_fine%nyc)=mg_fine%u_c(1:mg_fine%nxc,1:mg_fine%nyc)+e_h(1:mg_fine%nxc,1:mg_fine%nyc)

        deallocate(e_h)
    end subroutine prolongation_en_correct

    recursive subroutine V_cycle(mg_fine)
        !! A classic multigrid V-cycle
        implicit none

        type(GridSystem), intent(inout):: mg_fine
        type(GridSystem) :: mg_coarse
        complex(kind = realdp), allocatable, dimension(:,:)  :: op_Mx

        allocate(op_Mx(1-LAP:mg_fine%nxc+LAP,1-LAP:mg_fine%nyc+LAP))
        op_Mx=(0.d0, 0.d0)

        call Damp_Jacobi_smoother(mg_fine%u_c,mg_fine%rhs_c,mg_fine%nxc,mg_fine%nyc, &
                                  mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        call CSLP_OP_BC(mg_fine%u_c,op_Mx,mg_fine%nxc,mg_fine%nyc, &
                        mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        mg_fine%res_c = mg_fine%rhs_c - op_Mx
        op_Mx=(0.d0, 0.d0)
        deallocate(op_Mx)

        call coarsegrid_create(mg_coarse,mg_fine)
        call restriction(mg_coarse,mg_fine)

        if (MG_flag == 0) then
            !! A classic Two-Cycle
            if (cslp_mg_miter > 1) then !! GMRES or Bi-CGSTAB can be chose to solve the coarse-grid problem
                call mg_fullgmres(mg_coarse)
            else 
                call mg_bicgstab(mg_coarse)
            endif
        else
            if (mod(mg_coarse%nxc_global-1,2) == 0 .and. mg_coarse%nxc_global > nx_min  .and. &
                mod(mg_coarse%nyc_global-1,2) == 0 .and. mg_coarse%nyc_global > nx_min) then
                call V_cycle(mg_coarse)
            else
                if (cslp_mg_miter > 1) then !! GMRES or Bi-CGSTAB can be chose to solve the coarsest-grid problem
                    call mg_fullgmres(mg_coarse)
                else 
                    call mg_bicgstab(mg_coarse)
                endif
            endif
        endif

        call prolongation_en_correct(mg_coarse,mg_fine)

        call Damp_Jacobi_smoother(mg_fine%u_c,mg_fine%rhs_c,mg_fine%nxc,mg_fine%nyc, &
                                  mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        call grid_destroy(mg_coarse)

    end subroutine V_cycle

    recursive subroutine F_cycle(mg_fine)
        !! This is a classic F-cycle multigrid 
        implicit none

        type(GridSystem), intent(inout):: mg_fine
        type(GridSystem) :: mg_coarse
        complex(kind = realdp), allocatable, dimension(:,:)  :: op_Mx

        call Damp_Jacobi_smoother(mg_fine%u_c,mg_fine%rhs_c,mg_fine%nxc,mg_fine%nyc, &
                                  mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)
        
        allocate(op_Mx(1-LAP:mg_fine%nxc+LAP,1-LAP:mg_fine%nyc+LAP))
        op_Mx=(0.d0, 0.d0)
        call CSLP_OP_BC(mg_fine%u_c,op_Mx,mg_fine%nxc,mg_fine%nyc, &
                        mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        mg_fine%res_c = mg_fine%rhs_c - op_Mx
        op_Mx=(0.d0, 0.d0)
        deallocate(op_Mx)

        call coarsegrid_create(mg_coarse,mg_fine)
        call restriction(mg_coarse,mg_fine)

        if (mod(mg_coarse%nxc_global-1,2) == 0 .and. mg_coarse%nxc_global > nx_min  .and. &
            mod(mg_coarse%nyc_global-1,2) == 0 .and. mg_coarse%nyc_global > nx_min) then
            call F_cycle(mg_coarse)
        else
            if (cslp_mg_miter > 1) then
                call mg_fullgmres(mg_coarse)
            else 
                call mg_bicgstab(mg_coarse)
            endif
        endif

        call prolongation_en_correct(mg_coarse,mg_fine)

        call Damp_Jacobi_smoother(mg_fine%u_c,mg_fine%rhs_c,mg_fine%nxc,mg_fine%nyc, &
                                  mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        call grid_destroy(mg_coarse)

        allocate(op_Mx(1-LAP:mg_fine%nxc+LAP,1-LAP:mg_fine%nyc+LAP))
        op_Mx=(0.d0, 0.d0)
        call CSLP_OP_BC(mg_fine%u_c,op_Mx,mg_fine%nxc,mg_fine%nyc, &
                        mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        mg_fine%res_c = mg_fine%rhs_c - op_Mx
        op_Mx=(0.d0, 0.d0)
        deallocate(op_Mx)

        call coarsegrid_create(mg_coarse,mg_fine)
        call restriction(mg_coarse,mg_fine)

        if (mod(mg_coarse%nxc_global-1,2) == 0 .and. mg_coarse%nxc_global > nx_min  .and. &
            mod(mg_coarse%nyc_global-1,2) == 0 .and. mg_coarse%nyc_global > nx_min) then
                call V_cycle(mg_coarse)
            else
                if (cslp_mg_miter > 1) then
                    call mg_fullgmres(mg_coarse)
                else 
                    call mg_bicgstab(mg_coarse)
                endif
        endif

        call prolongation_en_correct(mg_coarse,mg_fine)

        call Damp_Jacobi_smoother(mg_fine%u_c,mg_fine%rhs_c,mg_fine%nxc,mg_fine%nyc, &
                                  mg_fine%hxc,mg_fine%hyc,mg_fine%kxy_c,mg_fine%kh2_c)

        call grid_destroy(mg_coarse)


    end subroutine F_cycle


    subroutine mg_fullgmres(mg_solve)
        !! This is a rountine to solve the (specified) coarsest-grid CSLP system approximately by using GMRES
        !! The CSLP is define by ReD-O2 scheme
        implicit none

        !Subroutine arguments -------------------------------------------------------
        type(GridSystem), intent(inout) :: mg_solve
            !! The current coarsest grid system
    
        !Local arguments ------------------------------------------------------------
        real(kind = realdp):: Rerror
            !! Relative residual
        integer :: iter
            !! The number of iterations
        integer                        :: k,ki,j
        real(kind = realdp)                   :: b_norm, res_norm
        complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
        complex(kind = realdp), allocatable, dimension(:,:)   :: Au0, H  
        complex(kind = realdp), allocatable, dimension(:,:,:) :: V
    
      ! Subroutine content ---------------------------------------------------------
        allocate(sn(cslp_mg_miter+1), cs(cslp_mg_miter+1), beta(cslp_mg_miter+1))
        allocate(Au0(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(H(cslp_mg_miter+1,cslp_mg_miter))
        allocate(V(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP,cslp_mg_miter+1))
    
        b_norm = mg_norm(mg_solve%rhs_c,mg_solve%nxc,mg_solve%nyc)
    
        Au0  = (0.d0,0.d0)
        sn   = (0.d0,0.d0)
        cs   = (0.d0,0.d0)
        V    = (0.d0,0.d0)
        H    = (0.d0,0.d0)
        beta = (0.d0,0.d0)
    
        call CSLP_OP_BC(mg_solve%u_c,Au0,mg_solve%nxc,mg_solve%nyc, &
                        mg_solve%hxc,mg_solve%hyc,mg_solve%kxy_c,mg_solve%kh2_c) ! Compute Ax

        mg_solve%res_c = mg_solve%rhs_c - Au0  !r=b-Ax
        res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)  ! ||r||

        if (res_norm == 0.) then
            if (my_id .eq. 0 ) then
                write(*,"(A)") "   CSLP RHS = 0 , Error = 0 "
            endif
        else
            Rerror   = res_norm / b_norm  ! scaled residual
            beta(1) = res_norm       ! beta(1)=||r0||
            V(:,:,1)  = mg_solve%res_c / res_norm  !  This is V(:,1) i.e. v1 in the algorithm
        
            k=0
            do j = 1,cslp_mg_miter
                k=k+1   !!Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
                call mg_arnoldi(V,H,k,mg_solve%nxc,mg_solve%nyc,mg_solve%hxc,mg_solve%hyc,&
                                mg_solve%kxy_c,mg_solve%kh2_c)
        
                call mg_apply_givens_rotation(H, cs, sn, k)
        
                beta(k+1) = -sn(k)*beta(k)
                beta(k)   =  conjg(cs(k))*beta(k)
        
                Rerror = CDABS(beta(k+1))/b_norm
        
                if (Rerror<1.d-8) then  !! The default tolerance is 1E-08
                    exit
                end if
            enddo
        
            if (my_id .eq. 0 ) then
             write(*,"(A,I9,A,E16.9)") "  Final MGsolve GMRES Iter.    ", k, "     Rel. res =", Rerror
            endif
        
            call mg_back_substitute(H,beta,k)
        
        
            do ki = 1,k
                mg_solve%u_c = mg_solve%u_c + beta(ki)*V(:,:,ki)
            enddo
        endif
    
        deallocate(sn, cs, beta)
        deallocate(Au0)!, res,u0
        deallocate(H)
        deallocate(V)
    
        iter = k
    
    end subroutine mg_fullgmres
      !====================================================================================
    
    subroutine mg_arnoldi(V,H,k,ni,nj,hx_c,hy_c,kxy,kh2)
        !! Arnoldi precess
        implicit none
    
      ! Subroutine arguments -------------------------------------------------------
        integer, intent(in) :: ni,nj,k
        real(kind = realdp),intent(in) :: hx_c,hy_c
    
        complex(kind = realdp),    dimension(1-LAP:ni+LAP,1-LAP:nj+LAP,cslp_mg_miter+1), intent(inout) :: V
        complex(kind = realdp),    dimension(cslp_mg_miter+1,cslp_mg_miter),                  intent(inout) :: H
        real   (kind = realdp),    dimension(1-LAP:ni+LAP,1-LAP:nj+LAP),            intent(in)    :: kxy, kh2
    
      ! Local arguments ------------------------------------------------------------
        integer                    :: i
    
      ! Subroutine content ---------------------------------------------------------
        call CSLP_OP_BC(V(:,:,k),V(:,:,k+1),ni,nj,hx_c,hy_c,kxy,kh2) 
            !!  w=A*v_i
    
        do i=1,k
    
            H(i, k)  = mg_dot_prod(V(:,:,i), V(:,:,k+1),ni,nj)     !! Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
            V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
    
        end do
    
        H(k+1, k) = mg_norm(V(:,:,k+1),ni,nj)
        V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)
    
    end subroutine mg_arnoldi
    
      ! SUBROUTINE APPLY_GIVENS_ROTATION ====================================================
    subroutine mg_apply_givens_rotation(H, cs, sn, k)
        !! This is routine to perform apply_givens_rotation
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        integer, intent(in) :: k
        complex(kind = realdp),    dimension(cslp_mg_miter+1,cslp_mg_miter), intent(inout)   :: H
        complex(kind = realdp),    dimension(cslp_mg_miter+1),          intent(inout)   :: cs, sn

        ! Local arguments ------------------------------------------------------------
        complex(kind = realdp)  :: temp
        integer  :: i

        ! Subroutine content ---------------------------------------------------------

        do i=1,k-1
            temp     =  conjg(cs(i))*H(i,k) + conjg(sn(i))*H(i+1,k)
            H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
            H(i,k)   =  temp
        end do

        if (H(k,k)==0.) then
            cs(k) = 0.0d0
            sn(k) = 1.0d0
        else
            temp  = CDSQRT((CDABS(H(k,k)))**2 + (H(k+1,k))**2)
            cs(k) = H(k,k) / temp
            sn(k) = H(k+1,k) / temp
        end if

        H(k,k)   = conjg(cs(k))*H(k,k) + conjg(sn(k))*H(k+1,k)
        H(k+1,k) = (0.0d0,0.0d0)
    
    end subroutine mg_apply_givens_rotation
    
    
    ! SUBROUTINE BACK_SUBSTITUTION =========================================================
    subroutine mg_back_substitute(H,beta,k)
        !! This is routine that performs back substitute
        implicit none
        integer :: i
        integer, intent(in) :: k
        complex(kind = realdp),    dimension(cslp_mg_miter+1,cslp_mg_miter), intent(in)    :: H
        complex(kind = realdp),    dimension(cslp_mg_miter+1),          intent(inout) :: beta

        beta(k) = beta(k)/H(k,k)

        do i=k-1,1,-1
            beta(i) = (beta(i) - sum(H(i,i+1:k)*beta(i+1:k)))/H(i,i)
        end do
    end subroutine mg_back_substitute

    !Algorithm Bi-CSGTAB=========================================================================
    subroutine mg_bicgstab(mg_solve, maximum_iterations)
        !! This is a rountine to solve the (specified) coarsest-grid CSLP system approximately by using Bi-CGSTAB
        !! The CSLP is define by ReD-O2 scheme
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(GridSystem), intent(inout) :: mg_solve

        ! optional arguments
        integer, optional, intent(in)   :: maximum_iterations 
            !! User-specified maximum number of iterations
    
        integer :: iter
            !! iter: the number of iterations
        integer :: maxit
            !! maxit: the maximum number of iterations, default is 3000, the default tolerance is 1E-08
        complex(kind = realdp), allocatable, dimension(:,:) :: res_hat, Au0, vi, p, s, t
        complex(kind=realdp)                                :: rho_cg, rho0_cg
        complex(kind=realdp)                                :: alpha_cg, omega_cg, beta_cg
        real(kind=realdp)                                   :: res_norm, b_norm       

        maxit = 3000
        if ( present(maximum_iterations) ) maxit = maximum_iterations

        allocate(res_hat(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(vi(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(p(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(s(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(t(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(Au0(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))

        
        rho_cg   = (1.d0,0.d0)
        alpha_cg = (1.d0,0.d0)
        omega_cg = (1.d0,0.d0)  

        vi    = (0.d0,0.d0)
        p     = (0.d0,0.d0)
        s     = (0.d0,0.d0)
        t     = (0.d0,0.d0)
        Au0   = (0.d0,0.d0)

        call CSLP_OP_BC(mg_solve%u_c,Au0,mg_solve%nxc,mg_solve%nyc, &
                        mg_solve%hxc,mg_solve%hyc,mg_solve%kxy_c,mg_solve%kh2_c)
        mg_solve%res_c = mg_solve%rhs_c - Au0  !r=b-Ax
        
        res_hat = mg_solve%res_c

        res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
        b_norm   = mg_norm(mg_solve%rhs_c,mg_solve%nxc,mg_solve%nyc)
        
        if (res_norm == 0.) then
            if (my_id .eq. 0 ) then
                write(*,"(A)") "   CSLP RHS = 0 , Error = 0 "
            endif
        else
            iter = 0 
            do while(res_norm /= 0.d0 .and. res_norm > 1.d-8*b_norm .and. iter < maxit) 
                rho0_cg = rho_cg
                rho_cg  = mg_dot_prod(res_hat,mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
                if(rho_cg == czero) then
                    write(*,*) "OMG, rho=0!"
                    stop
                endif

                beta_cg = (rho_cg/rho0_cg)*(alpha_cg/omega_cg)
                
                p       = mg_solve%res_c + beta_cg * (p - omega_cg*vi)

                call CSLP_OP_BC(p,vi,mg_solve%nxc,mg_solve%nyc, &
                                mg_solve%hxc,mg_solve%hyc,mg_solve%kxy_c,mg_solve%kh2_c)

                alpha_cg= rho_cg / mg_dot_prod(res_hat,vi,mg_solve%nxc,mg_solve%nyc)

                if(alpha_cg == czero) then
                    write(*,*) "OMG, alpha_cg=0!"
                    stop
                endif

                s = mg_solve%res_c - alpha_cg*vi
                
                if (mg_norm(s,mg_solve%nxc,mg_solve%nyc) < 1.d-8*b_norm) then
                    mg_solve%u_c = mg_solve%u_c + alpha_cg*p

                    Au0 = czero
                    call CSLP_OP_BC(mg_solve%u_c,Au0,mg_solve%nxc,mg_solve%nyc, &
                                    mg_solve%hxc,mg_solve%hyc,mg_solve%kxy_c,mg_solve%kh2_c)
                    mg_solve%res_c = mg_solve%rhs_c - Au0

                    res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)

                    Exit
                endif
                
                call CSLP_OP_BC(s,t,mg_solve%nxc,mg_solve%nyc, &
                                mg_solve%hxc,mg_solve%hyc,mg_solve%kxy_c,mg_solve%kh2_c)
            
                omega_cg  = mg_dot_prod(t,s,mg_solve%nxc,mg_solve%nyc) &
                        /mg_dot_prod(t,t,mg_solve%nxc,mg_solve%nyc)

                if(omega_cg == czero) then
                    write(*,*) "OMG, omega_cg=0!"
                    stop
                endif

                mg_solve%u_c   = mg_solve%u_c + alpha_cg*p + omega_cg*s
                
                mg_solve%res_c = s - omega_cg*t
                
                res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
                
                iter = iter + 1

            end do !-------> END OF While LOOP
            
            if (my_id .eq. 0 ) then
                write(*,"(A,I5,A,E16.9)") "            Final MGsolve BiCGSTAB Iter.    ", iter, "     Rel. res =", res_norm/b_norm
            end if
        endif
        
        deallocate(vi ,p, s, t, res_hat, Au0) 

    end subroutine mg_bicgstab

    subroutine ReD_Glk_CSLP_bicgstab(mg_solve, grid, maximum_iterations)
        !! This is a rountine to solve a coarse grid CSLP system approximately by using Bi-CGSTAB. 
        !! The CSLP operator is defined by ReD-Glk scheme. 
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(GridSystem), intent(inout) :: mg_solve
            !! The current coarse grid system
        type(Gridpara) :: grid
            !! Parameters of current coarse grid

        ! optional arguments
        integer, optional, intent(in)   :: maximum_iterations 

        integer :: iter
            !! iter: the number of iterations
        integer :: maxit
            !! maxit: the maximum number of iterations, default is 1000
        complex(kind = realdp), allocatable, dimension(:,:) :: res_hat, Au0, vi, p, s, t
        complex(kind=realdp)                                :: rho_cg, rho0_cg
        complex(kind=realdp)                                :: alpha_cg, omega_cg, beta_cg
        real(kind=realdp)                                   :: res_norm, b_norm, s_norm   

        !------------------------END PARAMETER AND VARIABLE-----------------------------!  

        maxit = 1000
        if ( present(maximum_iterations) ) maxit = maximum_iterations

        allocate(res_hat(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(vi(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(p(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(s(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(t(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(Au0(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))

        rho_cg   = (1.d0,0.d0)
        alpha_cg = (1.d0,0.d0)
        omega_cg = (1.d0,0.d0)  
        
        vi    = (0.d0,0.d0)
        p     = (0.d0,0.d0)
        s     = (0.d0,0.d0)
        t     = (0.d0,0.d0)
        Au0   = (0.d0,0.d0)

        Au0=CSLP_Mx_nth(mg_solve%u_c,grid) ! Ax
            !! CSLP operator by ReD-Glk 
        mg_solve%res_c = mg_solve%rhs_c - Au0  !r=b-Ax
        
        res_hat = mg_solve%res_c

        res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
        b_norm   = mg_norm(mg_solve%rhs_c,mg_solve%nxc,mg_solve%nyc)
        
        if (res_norm == 0.) then
            if (my_id .eq. 0 ) then
                write(*,"(A)") "   CSLP RHS = 0 , Error = 0 "
            endif
        else
            iter = 0 
            do while(res_norm /= 0.d0 .and. res_norm > cslp_mg_tol*b_norm .and. iter < maxit) 
                rho0_cg = rho_cg
                rho_cg  = mg_dot_prod(res_hat,mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
                if(rho_cg == czero) then
                    write(*,*) "OMG, rho=0!"
                    stop
                endif

                beta_cg = (rho_cg/rho0_cg)*(alpha_cg/omega_cg)
                
                p       = mg_solve%res_c + beta_cg * (p - omega_cg*vi)

                vi = CSLP_Mx_nth(p,grid)

                alpha_cg= rho_cg / mg_dot_prod(res_hat,vi,mg_solve%nxc,mg_solve%nyc)

                if(alpha_cg == czero) then
                    write(*,*) "OMG, alpha_cg=0!"
                    stop
                endif

                s = mg_solve%res_c - alpha_cg*vi
                s_norm = mg_norm(s,mg_solve%nxc,mg_solve%nyc)
                if ( s_norm < cslp_mg_tol*1.d-1*b_norm) then
                    mg_solve%u_c = mg_solve%u_c + alpha_cg*p

                    Au0 = czero
                    Au0 = CSLP_Mx_nth(mg_solve%u_c,grid)
                    mg_solve%res_c = mg_solve%rhs_c - Au0

                    res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)

                    Exit
                endif

                t = CSLP_Mx_nth(s,grid)
            
                omega_cg  = mg_dot_prod(t,s,mg_solve%nxc,mg_solve%nyc) &
                        /mg_dot_prod(t,t,mg_solve%nxc,mg_solve%nyc)

                if(omega_cg == czero) then
                    write(*,*) "OMG, omega_cg=0!"
                    stop
                endif

                mg_solve%u_c   = mg_solve%u_c + alpha_cg*p + omega_cg*s
                
                mg_solve%res_c = s - omega_cg*t
                
                res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)
                
                iter = iter + 1

            end do !-------> END OF While LOOP
            
            if (my_id .eq. 0 ) then
                write(*,"(A,I5,A,E16.9)") "            Final CSLPsolve BiCGSTAB Iter.    ", iter, "     Rel. res =", res_norm/b_norm
            end if
        endif
        
        deallocate(vi ,p, s, t, res_hat, Au0)

    end subroutine ReD_Glk_CSLP_bicgstab     


    subroutine ReD_Glk_CSLP_gmres(mg_solve, grid, maximum_iterations)
        !! This is a rountine to solve a coarse grid CSLP system approximately by using GMRES
        !! The CSLP operator is defined by ReD-Glk scheme
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(GridSystem), intent(inout) :: mg_solve
        type(Gridpara) :: grid

        ! optional arguments
        integer, optional, intent(in)   :: maximum_iterations 
            !! User-specified maximum number of iterations,default is 300

        ! Local arguments ------------------------------------------------------------
        real(kind = realdp):: Rerror
        integer :: iter,maxit,k,ki,i,j
        real(kind = realdp) :: b_norm, res_norm
        complex(kind = realdp)  :: temp
        complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
        complex(kind = realdp), allocatable, dimension(:,:)   :: Au0, H
        complex(kind = realdp), allocatable, dimension(:,:,:) :: V
        integer(kind=4) :: level_i
        
        ! optional present
        maxit = 200
        if ( present(maximum_iterations) ) maxit = maximum_iterations
    
        ! Subroutine content ---------------------------------------------------------
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1

        allocate(sn(maxit+1), cs(maxit+1), beta(maxit+1))
        allocate(Au0(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP))
        allocate(H(maxit+1,maxit))
        allocate(V(1-LAP:mg_solve%nxc+LAP,1-LAP:mg_solve%nyc+LAP,maxit+1))
    
        b_norm = mg_norm(mg_solve%rhs_c,mg_solve%nxc,mg_solve%nyc)
    
        Au0  = czero
        sn   = czero
        cs   = czero
        V    = czero
        H    = czero
        beta = czero
    
        Au0=CSLP_Mx_nth(mg_solve%u_c,grid) ! Compute Ax
            !! CSLP operator by ReD-Glk 

        mg_solve%res_c = mg_solve%rhs_c - Au0  !r=b-Ax
        res_norm = mg_norm(mg_solve%res_c,mg_solve%nxc,mg_solve%nyc)  ! ||r||

        Au0 = mg_solve%u_c ! storge the initial u instead

        if (res_norm == 0.) then
            if (my_id .eq. 0 ) then
                write(*,"(A)") "   CSLP RHS = 0 , Error = 0 "
            endif
        else
            Rerror   = res_norm / b_norm  ! scaled error
            beta(1) = res_norm       ! beta(1)=||r0||
            V(:,:,1)  = mg_solve%res_c / res_norm  !  This is V(:,1) i.e. v1 in the algorithm
        
            k=0
            do j = 1,maxit
                k=k+1   !!Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
                !---------------------START arnoldi process----------------------------------
                V(:,:,k+1) = CSLP_Mx_nth(V(:,:,k), grid)

                do i=1,k
                    H(i, k)  = mg_dot_prod(V(:,:,i), V(:,:,k+1),mg_solve%nxc,mg_solve%nyc)     !! Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
                    V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
                end do
    
                H(k+1, k) = mg_norm(V(:,:,k+1),mg_solve%nxc,mg_solve%nyc)
                V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)
                !---------------------END arnoldi process----------------------------------
        

                !---------------------START givens_rotation----------------------------------
                do i=1,k-1
                    temp     =  conjg(cs(i))*H(i,k) + conjg(sn(i))*H(i+1,k)
                    H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
                    H(i,k)   =  temp
                end do
                
                if (H(k,k)==0) then
                    cs(k) = 0.0d0
                    sn(k) = 1.0d0
                else
                    temp  = CDSQRT((CDABS(H(k,k)))**2 + (H(k+1,k))**2)
                    cs(k) = H(k,k) / temp
                    sn(k) = H(k+1,k) / temp
                end if
                
                H(k,k)   = conjg(cs(k))*H(k,k) + conjg(sn(k))*H(k+1,k)
                H(k+1,k) = (0.0d0,0.0d0)
                !---------------------END givens_rotation----------------------------------
        
                
                beta(k+1) = -sn(k)*beta(k)
                beta(k)   =  conjg(cs(k))*beta(k)
        
                Rerror = CDABS(beta(k+1))/b_norm
        
                if (Rerror<cslp_mg_tol) then ! The tolerance for the coarse-grid CSLP solver is set in the Input file
                    exit
                end if
            enddo
            
            if (my_id .eq. 0 ) then
                write(*,"(A,I9,A,E16.9)") "  Final CSLPsolve GMRES Iter.    ", k, "     Rel. res =", Rerror
            endif
        
            !!---call back_substitute(H,beta,k)---------------------------------
            beta(k) = beta(k)/H(k,k)
            do i=k-1,1,-1
                beta(i) = (beta(i) - sum(H(i,i+1:k)*beta(i+1:k)))/H(i,i)
            end do
            !!------------------------------------------------------------------
            mg_solve%u_c = czero ! make use of the memmory, use as a temp variable here
            do ki = 1,k
                mg_solve%u_c = mg_solve%u_c+ beta(ki)*V(:,:,ki)
            enddo
            mg_solve%u_c = mg_solve%u_c+Au0
        endif
    
        deallocate(sn, cs, beta)
        deallocate(Au0)
        deallocate(H)
        deallocate(V)
    
        iter = k

    end subroutine ReD_Glk_CSLP_gmres 

end module CSLP_Solver