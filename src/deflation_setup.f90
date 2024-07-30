module deflaion_setup
    !! This is an integrated module for deflation preconditioning 
    use mpi
    use comm_variable
    use mpi_setup
    use wavenumber, only:wavenumber_k, kh_pow_2
    use operators
    use CSLP_Solver
    implicit none

    type, public  :: TwoGrids
        !! Type definition for a two-grid system. "f" in variable names stands for fine grid, "c" in variable names stands for coarse grid
        integer(kind=4)     :: nxf_global, nyf_global
        integer(kind=4)     :: nxf, nyf
        integer(kind=4)     :: nxc_global, nyc_global
        integer(kind=4)     :: nxc, nyc

        real(kind = realdp)        :: hxf, hyf, hxhyf
        real(kind = realdp)        :: hxc, hyc, hxhyc

        integer, dimension(0:npMax - 1) :: if_offset, jf_offset, if_nn, jf_nn
        integer, dimension(0:npMax - 1) :: ic_offset, jc_offset, ic_nn, jc_nn

        real(kind = realdp), allocatable, dimension(:,:) :: kxy_f, kxy_c, kh2_f, kh2_c
        
    end type TwoGrids

    interface TwoGrids
        procedure :: FromFine2Coarse
    end interface TwoGrids

contains

    function DEF_Px(x)
        !! Identifier for different deflation methods, ONLY for an input array from the default (finest) grid system
        implicit none
        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: x
        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: DEF_Px

        type(Gridpara) :: finestgrid
        type(TwoGrids) :: f2c

        call default_gridpara(finestgrid)
        f2c = TwoGrids(finestgrid)

        select case (M_flag)
  
            case (2)
            DEF_Px = P_DEFx(x,f2c) !! A basic deflation method
    
            case (3) 
            DEF_Px = P_ADEF1x(x,f2c) !! Adapted deflation method, including higher-order deflation vectors
    
            case (4)
            DEF_Px = P_TLKMx(x,f2c) !! Two-level Krylov Method

            case (5)
            DEF_Px = MultiLevelADP_Px(x,finestgrid) !! Multilevel deflation method
    
            case default
            DEF_Px = P_ADEF1x(x,f2c)
         
        end select

        deallocate(f2c%kxy_f,f2c%kh2_f,f2c%kxy_c,f2c%kh2_c)
        deallocate(finestgrid%wavenumber_k,finestgrid%wavenumber_kh_pow_2)

    end function DEF_Px

    subroutine default_gridpara(finestgrid)
        !! This is a rountine to define the default finest Gridpara from the common grid parameters of the whole project 
        implicit none

        type (Gridpara), intent(inout) :: finestgrid

        finestgrid%nx_global = nx_global
        finestgrid%ny_global = ny_global
        finestgrid%hx        = hx
        finestgrid%hy        = hy
        finestgrid%hxhy      = hx*hy
        finestgrid%nx        = nx
        finestgrid%ny        = ny
        finestgrid%i_offset  = i_offset
        finestgrid%j_offset  = j_offset
        finestgrid%i_nn      = i_nn
        finestgrid%j_nn      = j_nn

        allocate(finestgrid%wavenumber_k(1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(finestgrid%wavenumber_kh_pow_2(1-LAP:nx+LAP,1-LAP:ny+LAP))
        finestgrid%wavenumber_k = wavenumber_k
        finestgrid%wavenumber_kh_pow_2 = kh_pow_2

    end subroutine default_gridpara

    type(TwoGrids) function FromFine2Coarse(OneGrid)
        !! A procedure to create a two-grid system from a given fine Gridpara type
        implicit none

        type(Gridpara), intent(inout) :: OneGrid

        integer(kind=4) :: k

        FromFine2Coarse%nxf_global = OneGrid%nx_global
        FromFine2Coarse%nyf_global = OneGrid%ny_global
        FromFine2Coarse%hxf        = OneGrid%hx
        FromFine2Coarse%hyf        = OneGrid%hy
        FromFine2Coarse%hxhyf      = OneGrid%hxhy
        FromFine2Coarse%nxf        = OneGrid%nx
        FromFine2Coarse%nyf        = OneGrid%ny
        FromFine2Coarse%if_offset  = OneGrid%i_offset
        FromFine2Coarse%jf_offset  = OneGrid%j_offset
        FromFine2Coarse%if_nn      = OneGrid%i_nn
        FromFine2Coarse%jf_nn      = OneGrid%j_nn

        allocate(FromFine2Coarse%kxy_f(1-LAP:FromFine2Coarse%nxf+LAP,1-LAP:FromFine2Coarse%nyf+LAP))
        FromFine2Coarse%kxy_f  = OneGrid%wavenumber_k

        allocate(FromFine2Coarse%kh2_f(1-LAP:FromFine2Coarse%nxf+LAP,1-LAP:FromFine2Coarse%nyf+LAP))
        FromFine2Coarse%kh2_f = OneGrid%wavenumber_kh_pow_2


        FromFine2Coarse%nxc_global =(OneGrid%nx_global-1)/2+1
        FromFine2Coarse%nyc_global =(OneGrid%ny_global-1)/2+1

        FromFine2Coarse%hxc=slx/(FromFine2Coarse%nxc_global-1)
        FromFine2Coarse%hyc=sly/(FromFine2Coarse%nyc_global-1)
        FromFine2Coarse%hxhyc = FromFine2Coarse%hxc*FromFine2Coarse%hyc

        if (mod(OneGrid%i_offset(npx),2) == 0) then
            if (mod(OneGrid%nx,2) == 0) then
                FromFine2Coarse%nxc=OneGrid%nx / 2
            else
                FromFine2Coarse%nxc=(OneGrid%nx - 1) / 2
            endif
        else
            if (mod(OneGrid%nx,2) == 0) then
                FromFine2Coarse%nxc=OneGrid%nx / 2
            else
                FromFine2Coarse%nxc=(OneGrid%nx + 1) / 2
            endif
        endif

        if (mod(OneGrid%j_offset(npy),2) == 0) then
            if (mod(OneGrid%ny,2) == 0) then
                FromFine2Coarse%nyc=OneGrid%ny / 2
            else
                FromFine2Coarse%nyc=(OneGrid%ny - 1) / 2
            endif
        else
            if (mod(OneGrid%ny,2) == 0) then
                FromFine2Coarse%nyc=OneGrid%ny / 2
            else
                FromFine2Coarse%nyc=(OneGrid%ny + 1) / 2
            endif
        endif


        do k=0,npx0-1
            if (mod(OneGrid%i_offset(k),2) == 0) then

                FromFine2Coarse%ic_offset(k) = (OneGrid%i_offset(k)+1+1) / 2

                if (mod(OneGrid%i_nn(k),2) == 0) then
                    FromFine2Coarse%ic_nn(k)=OneGrid%i_nn(k) / 2
                else
                    FromFine2Coarse%ic_nn(k)=(OneGrid%i_nn(k) - 1) / 2
                endif

            else

                FromFine2Coarse%ic_offset(k) = (OneGrid%i_offset(k)+1) / 2

                if (mod(OneGrid%i_nn(k),2) == 0) then
                    FromFine2Coarse%ic_nn(k)=OneGrid%i_nn(k) / 2
                else
                    FromFine2Coarse%ic_nn(k)=(OneGrid%i_nn(k) + 1) / 2
                endif

            endif
        enddo

        do k=0,npy0-1
            if (mod(OneGrid%j_offset(k),2) == 0) then

                FromFine2Coarse%jc_offset(k) = (OneGrid%j_offset(k)+1+1) / 2

                if (mod(OneGrid%j_nn(k),2) == 0) then
                    FromFine2Coarse%jc_nn(k)=OneGrid%j_nn(k) / 2
                else
                    FromFine2Coarse%jc_nn(k)=(OneGrid%j_nn(k) - 1) / 2
                endif
            else

                FromFine2Coarse%jc_offset(k) = (OneGrid%j_offset(k)+1) / 2

                if (mod(OneGrid%j_nn(k),2) == 0) then
                    FromFine2Coarse%jc_nn(k)=OneGrid%j_nn(k) / 2
                else
                    FromFine2Coarse%jc_nn(k)=(OneGrid%j_nn(k) + 1) / 2
                endif
            endif
        enddo

        allocate(FromFine2Coarse%kxy_c(1-LAP:FromFine2Coarse%nxc+LAP,1-LAP:FromFine2Coarse%nyc+LAP))
        call wavenumber_FWrestriction(FromFine2Coarse) !! The wavenumber of the coarse level is restricted from the fine level

        allocate(FromFine2Coarse%kh2_c(1-LAP:FromFine2Coarse%nxc+LAP,1-LAP:FromFine2Coarse%nyc+LAP))
        FromFine2Coarse%kh2_c = 0.d0 !This is for safe
        FromFine2Coarse%kh2_c = FromFine2Coarse%kxy_c * FromFine2Coarse%kxy_c * FromFine2Coarse%hxhyc

        return
    end function FromFine2Coarse

    subroutine wavenumber_FWrestriction(f2c)
        !! This routine obtains the wavenumber of the coarse level by full-weight restriction from the fine level
        implicit none

        type(TwoGrids), intent(inout)  :: f2c

        integer :: ii,jj,iif,jjf,icH_global,jcH_global,ifh_global,jfh_global

        call mg_checkreal_xy2d(f2c%kxy_f,f2c%nxf,f2c%nyf)
        call mg_checkreal_xy2d(f2c%kh2_f,f2c%nxf,f2c%nyf)
        ! Before adding the next line, MP2BC1_nx65k40_pre3_ReD-Glk_DirAveGhost_np4 is not consistent even diverge.
        f2c%kxy_c = 0.d0 !This is the most significant line

        do jj=1,f2c%nyc
            do ii=1,f2c%nxc
                icH_global=f2c%ic_offset(npx)-1+ii
                jcH_global=f2c%jc_offset(npy)-1+jj
                ifh_global=2*icH_global-1
                jfh_global=2*jcH_global-1
                iif=ifh_global-(f2c%if_offset(npx)-1)
                jjf=jfh_global-(f2c%jf_offset(npy)-1)

                f2c%kxy_c(ii,jj) = ( 4.d0*f2c%kxy_f(iif,jjf) &
                                    +2.d0*f2c%kxy_f(iif-1,jjf)+2.d0*f2c%kxy_f(iif+1,jjf) &
                                    +2.d0*f2c%kxy_f(iif,jjf-1)+2.d0*f2c%kxy_f(iif,jjf+1) &
                                    +1.d0*f2c%kxy_f(iif-1,jjf-1)+1.d0*f2c%kxy_f(iif-1,jjf+1) &
                                    +1.d0*f2c%kxy_f(iif+1,jjf-1)+1.d0*f2c%kxy_f(iif+1,jjf+1))/16.d0

                if (npx == 0 .and. ii==1) then
                    f2c%kxy_c(ii,jj)=f2c%kxy_f(iif,jjf)
                endif

                if (npx == npx0-1 .and. ii==f2c%nxc) then
                    f2c%kxy_c(ii,jj)=f2c%kxy_f(iif,jjf)
                endif

                if (npy == 0 .and. jj==1) then
                    f2c%kxy_c(ii,jj)=f2c%kxy_f(iif,jjf)
                endif

                if (npy == npy0-1 .and. jj==f2c%nyc) then
                    f2c%kxy_c(ii,jj)=f2c%kxy_f(iif,jjf)
                endif

            enddo
        enddo

        if (npx == 0 .and. npy == 0 ) then
            f2c%kxy_c(1,1)=f2c%kxy_f(1,1)
        endif

        if (npx == npx0-1 .and. npy == 0) then
            f2c%kxy_c(f2c%nxc,1)=f2c%kxy_f(f2c%nxf,1)
        endif

        if (npx == 0 .and. npy == npy0-1 ) then
            f2c%kxy_c(1,f2c%nyc)=f2c%kxy_f(1,f2c%nyf)
        endif

        if (npx == npx0-1 .and. npy == npy0-1) then
            f2c%kxy_c(f2c%nxc,f2c%nyc)=f2c%kxy_f(f2c%nxf,f2c%nyf)
        endif

        call mg_checkreal_xy2d(f2c%kxy_c,f2c%nxc,f2c%nyc)

    end subroutine wavenumber_FWrestriction

    type(Gridpara) function CoarseGridpara(f2c)
        !! A procedure to create a Gridpara type of the coarse grid system from a given a two-grid system
        implicit none

        type(TwoGrids), intent(inout) :: f2c

        CoarseGridpara%nx_global = f2c%nxc_global
        CoarseGridpara%ny_global = f2c%nyc_global
        CoarseGridpara%hx        = f2c%hxc
        CoarseGridpara%hy        = f2c%hyc
        CoarseGridpara%hxhy      = f2c%hxhyc
        CoarseGridpara%nx        = f2c%nxc
        CoarseGridpara%ny        = f2c%nyc
        CoarseGridpara%i_offset  = f2c%ic_offset
        CoarseGridpara%j_offset  = f2c%jc_offset
        CoarseGridpara%i_nn      = f2c%ic_nn
        CoarseGridpara%j_nn      = f2c%jc_nn

        allocate(CoarseGridpara%wavenumber_k(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        CoarseGridpara%wavenumber_k = f2c%kxy_c

        allocate(CoarseGridpara%wavenumber_kh_pow_2(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        CoarseGridpara%wavenumber_kh_pow_2 = f2c%kh2_c

    end function CoarseGridpara

    function ZTx(x,f2c) 
        !! Restriction of a variable, by full-weight restriction or higher-order restriction depending on the laylers of overlapping grid points used.
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        complex(kind=realdp) :: ZTx(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        integer :: ii,jj,iif,jjf,icH_global,jcH_global,ifh_global,jfh_global

        x(1-LAP:0,:) = czero
        x(:,1-LAP:0) = czero

        x(f2c%nxf+1:f2c%nxf+LAP,:) = czero
        x(:,f2c%nyf+1:f2c%nyf+LAP) = czero

        call mg_check_xy2d(x,f2c%nxf,f2c%nyf)

        ZTx=czero

        if(LAP == 1) then
            do jj=1,f2c%nyc
                do ii=1,f2c%nxc
                    icH_global=f2c%ic_offset(npx)-1+ii
                    jcH_global=f2c%jc_offset(npy)-1+jj
                    ifh_global=2*icH_global-1
                    jfh_global=2*jcH_global-1
                    iif=ifh_global-(f2c%if_offset(npx)-1)
                    jjf=jfh_global-(f2c%jf_offset(npy)-1)


                    ZTx(ii,jj)=(4.d0*x(iif,jjf) &
                                +2.d0*x(iif-1,jjf)+2.d0*x(iif+1,jjf) &
                                +2.d0*x(iif,jjf-1)+2.d0*x(iif,jjf+1) &
                                +1.d0*x(iif-1,jjf-1)+1.d0*x(iif-1,jjf+1) &
                                +1.d0*x(iif+1,jjf-1)+1.d0*x(iif+1,jjf+1))/4.d0

                    if (npx == 0 .and. ii==1) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npx == npx0-1 .and. ii==f2c%nxc) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npy == 0 .and. jj==1) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npy == npy0-1 .and. jj==f2c%nyc) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                enddo
            enddo
        else
            do jj=1,f2c%nyc
                do ii=1,f2c%nxc
                    icH_global=f2c%ic_offset(npx)-1+ii
                    jcH_global=f2c%jc_offset(npy)-1+jj
                    ifh_global=2*icH_global-1
                    jfh_global=2*jcH_global-1
                    iif=ifh_global-(f2c%if_offset(npx)-1)
                    jjf=jfh_global-(f2c%jf_offset(npy)-1)


                    ZTx(ii,jj)=( 1.d0*x(iif-2,jjf+2) + 4.d0*x(iif-1,jjf+2) + 6.d0*x(iif,jjf+2) &
                               + 4.d0*x(iif+1,jjf+2) + 1.d0*x(iif+2,jjf+2) &
                               + 4.d0*x(iif-2,jjf+1) + 16.d0*x(iif-1,jjf+1) + 24.d0*x(iif,jjf+1) &
                               + 16.d0*x(iif+1,jjf+1) + 4.d0*x(iif+2,jjf+1) &
                               + 6.d0*x(iif-2,jjf)   + 24.d0*x(iif-1,jjf)   + 36.d0*x(iif,jjf)   &
                               + 24.d0*x(iif+1,jjf)   + 6.d0*x(iif+2,jjf) &
                               + 4.d0*x(iif-2,jjf-1) + 16.d0*x(iif-1,jjf-1) + 24.d0*x(iif,jjf-1) &
                               + 16.d0*x(iif+1,jjf-1) + 4.d0*x(iif+2,jjf-1) &
                               + 1.d0*x(iif-2,jjf-2) + 4.d0*x(iif-1,jjf-2) + 6.d0*x(iif,jjf-2) &
                               + 4.d0*x(iif+1,jjf-2) + 1.d0*x(iif+2,jjf-2))/64.d0

                    if (npx == 0 .and. ii==1) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npx == npx0-1 .and. ii==f2c%nxc) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npy == 0 .and. jj==1) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                    if (npy == npy0-1 .and. jj==f2c%nyc) then
                        ZTx(ii,jj)=x(iif,jjf)
                    endif

                enddo
            enddo
        endif

        if (npx == 0 .and. npy == 0 ) then
            ZTx(1,1)=x(1,1)
        endif

        if (npx == npx0-1 .and. npy == 0) then
            ZTx(f2c%nxc,1)=x(f2c%nxf,1)
        endif

        if (npx == 0 .and. npy == npy0-1 ) then
            ZTx(1,f2c%nyc)=x(1,f2c%nyf)
        endif

        if (npx == npx0-1 .and. npy == npy0-1) then
            ZTx(f2c%nxc,f2c%nyc)=x(f2c%nxf,f2c%nyf)
        endif

    end function ZTx

    function Zx(x,f2c) 
        !! Interpolation of a variable from coarse to fine, by bilinear or higher-order interpolation depending on the laylers of overlapping grid points used.
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        complex(kind=realdp) :: Zx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        integer :: ii,jj,iif,jjf,icH_global,jcH_global,ifh_global,jfh_global

        Zx=czero

        x(1-LAP:0,:) = czero
        x(:,1-LAP:0) = czero

        x(f2c%nxc+1:f2c%nxc+LAP,:) = czero
        x(:,f2c%nyc+1:f2c%nyc+LAP) = czero

        call mg_check_xy2d(x,f2c%nxc,f2c%nyc)
        if (LAP == 1) then
            do jj=1,f2c%nyc
                do ii=1,f2c%nxc
                    icH_global=f2c%ic_offset(npx)-1+ii
                    jcH_global=f2c%jc_offset(npy)-1+jj
                    ifh_global=2*icH_global-1
                    jfh_global=2*jcH_global-1
                    iif=ifh_global-(f2c%if_offset(npx)-1)
                    jjf=jfh_global-(f2c%jf_offset(npy)-1)

                    Zx(iif,jjf) = x(ii,jj)

                    Zx(iif,jjf+1) = (x(ii,jj) + x(ii,jj+1))/2.d0
                    Zx(iif,jjf-1) = (x(ii,jj) + x(ii,jj-1))/2.d0
                    Zx(iif+1,jjf) = (x(ii,jj) + x(ii+1,jj))/2.d0
                    Zx(iif-1,jjf) = (x(ii,jj) + x(ii-1,jj))/2.d0

                    Zx(iif+1,jjf+1) = (x(ii+1,jj) + x(ii,jj+1) &
                                    + x(ii+1,jj+1) + x(ii,jj))/4.d0
                    Zx(iif+1,jjf-1) = (x(ii+1,jj) + x(ii,jj-1) &
                                    + x(ii+1,jj-1) + x(ii,jj))/4.d0
                    Zx(iif-1,jjf+1) = (x(ii-1,jj) + x(ii,jj+1) &
                                    + x(ii-1,jj+1) + x(ii,jj))/4.d0
                    Zx(iif-1,jjf-1) = (x(ii-1,jj) + x(ii,jj-1) &
                                    + x(ii-1,jj-1) + x(ii,jj))/4.d0

                enddo
            enddo
        else
            do jj=1,f2c%nyc
                do ii=1,f2c%nxc
                    icH_global=f2c%ic_offset(npx)-1+ii
                    jcH_global=f2c%jc_offset(npy)-1+jj
                    ifh_global=2*icH_global-1
                    jfh_global=2*jcH_global-1
                    iif=ifh_global-(f2c%if_offset(npx)-1)
                    jjf=jfh_global-(f2c%jf_offset(npy)-1)

                    Zx(iif,jjf) = (1.d0*x(ii-1,jj-1)+6.d0*x(ii,jj-1)+1.d0*x(ii+1,jj-1) &
                                  +6.d0*x(ii-1,jj)+36.d0*x(ii,jj)+6.d0*x(ii+1,jj)      &
                                  +1.d0*x(ii-1,jj+1)+6.d0*x(ii,jj+1)+1.d0*x(ii+1,jj+1))*0.015625d0

                    Zx(iif,jjf+1) = (1.d0*x(ii-1,jj)   + 6.d0*x(ii,jj)   + 1.d0*x(ii+1,jj) &
                                   + 1.d0*x(ii-1,jj+1) + 6.d0*x(ii,jj+1) + 1.d0*x(ii+1,jj+1))*0.0625d0

                    Zx(iif,jjf-1) = (1.d0*x(ii-1,jj)   + 6.d0*x(ii,jj)   + 1.d0*x(ii+1,jj) &
                                   + 1.d0*x(ii-1,jj-1) + 6.d0*x(ii,jj-1) + 1.d0*x(ii+1,jj-1))*0.0625d0

                    Zx(iif+1,jjf) = (1.d0*x(ii,jj-1)   + 6.d0*x(ii,jj)   + 1.d0*x(ii,jj+1) &
                                   + 1.d0*x(ii+1,jj-1) + 6.d0*x(ii+1,jj) + 1.d0*x(ii+1,jj+1))*0.0625d0

                    Zx(iif-1,jjf) = (1.d0*x(ii,jj-1)   + 6.d0*x(ii,jj)   + 1.d0*x(ii,jj+1) &
                                   + 1.d0*x(ii-1,jj-1) + 6.d0*x(ii-1,jj) + 1.d0*x(ii-1,jj+1))*0.0625d0

                    Zx(iif+1,jjf+1) = (x(ii+1,jj) + x(ii,jj+1) &
                                    + x(ii+1,jj+1) + x(ii,jj))*0.25d0
                    Zx(iif+1,jjf-1) = (x(ii+1,jj) + x(ii,jj-1) &
                                    + x(ii+1,jj-1) + x(ii,jj))*0.25d0
                    Zx(iif-1,jjf+1) = (x(ii-1,jj) + x(ii,jj+1) &
                                    + x(ii-1,jj+1) + x(ii,jj))*0.25d0
                    Zx(iif-1,jjf-1) = (x(ii-1,jj) + x(ii,jj-1) &
                                    + x(ii-1,jj-1) + x(ii,jj))*0.25d0


                    if (npx == 0 .and. ii==1) then
                        Zx(iif,jjf) = (x(ii,jj-1) + 6.d0*x(ii,jj) + x(ii,jj+1))*0.125d0
                        Zx(iif,jjf+1) = (x(ii,jj) + x(ii,jj+1))*0.5d0
                        Zx(iif,jjf-1) = (x(ii,jj) + x(ii,jj-1))*0.5d0
                    endif

                    if (npx == npx0-1 .and. ii==f2c%nxc) then
                        Zx(iif,jjf) = (x(ii,jj-1) + 6.d0*x(ii,jj) + x(ii,jj+1))*0.125d0
                        Zx(iif,jjf+1) = (x(ii,jj) + x(ii,jj+1))*0.5d0
                        Zx(iif,jjf-1) = (x(ii,jj) + x(ii,jj-1))*0.5d0
                    endif

                    if (npy == 0 .and. jj==1) then
                        Zx(iif,jjf) = (x(ii-1,jj) + 6.d0*x(ii,jj) + x(ii+1,jj))*0.125d0
                        Zx(iif+1,jjf) = (x(ii+1,jj) + x(ii,jj))*0.5d0
                        Zx(iif-1,jjf) = (x(ii-1,jj) + x(ii,jj))*0.5d0
                    endif

                    if (npy == npy0-1 .and. jj==f2c%nyc) then
                        Zx(iif,jjf) = (x(ii-1,jj) + 6.d0*x(ii,jj) + x(ii+1,jj))*0.125d0
                        Zx(iif+1,jjf) = (x(ii+1,jj) + x(ii,jj))*0.5d0
                        Zx(iif-1,jjf) = (x(ii-1,jj) + x(ii,jj))*0.5d0
                    endif

                enddo
            enddo

            if (npx == 0 .and. npy == 0 ) then
                Zx(1,1)=x(1,1)
            endif
    
            if (npx == npx0-1 .and. npy == 0) then
                Zx(f2c%nxf,1)=x(f2c%nxc,1)
            endif
    
            if (npx == 0 .and. npy == npy0-1 ) then
                Zx(1,f2c%nyf)=x(1,f2c%nyc)
            endif
    
            if (npx == npx0-1 .and. npy == npy0-1) then
                Zx(f2c%nxf,f2c%nyf)=x(f2c%nxc,f2c%nyc)
            endif
        endif

    end function Zx

    function ZTZx(x,f2c) 
        !! A routine that first perform restriction and then interpolation
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind=realdp) :: ZTZx(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        complex(kind = realdp), allocatable, dimension(:,:) :: temp
        allocate(temp(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))

        temp = czero
        ZTZx=czero

        temp = Zx(x,f2c)
        ZTZx = ZTx(temp,f2c)

        deallocate(temp)

    end function ZTZx

    function Helm_Ahx(x,f2c)
        !! A function that performs the Helmholtz operator on the fine grid of the two-grid system, by ReD-O2 method
        implicit none
        !input
        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        !output
        complex(kind=realdp) :: Helm_Ahx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        
        
        Helm_Ahx = czero
        call Helmholtz2d_BC_mg(x,Helm_Ahx,f2c%nxf,f2c%nyf,f2c%hxf,f2c%hyf,f2c%kxy_f,f2c%kh2_f)

    end function

    function Helm_A2hx(x,f2c)
        !! A function that performs the Helmholtz operator on the coarse grid of the two-grid system, by different methods
        implicit none
        !input
        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        !output
        complex(kind=realdp) :: Helm_A2hx(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        Helm_A2hx = czero

        A2h_flag = 2

        select case (A2h_flag)
  
            case (1) !! ReD-O2
            call Helmholtz2d_BC_mg(x,Helm_A2hx,f2c%nxc,f2c%nyc,f2c%hxc,f2c%hyc,f2c%kxy_c,f2c%kh2_c)

            case (2) !! ReD-Glk
            call Helmholtz2d_ReD_Glk(x,Helm_A2hx,f2c)

            case (3) !! ReD-cmpO4
            call Helmholtz2d_O4cmpct(x,Helm_A2hx,f2c)

            case default
            call Helmholtz2d_BC_mg(x,Helm_A2hx,f2c%nxc,f2c%nyc,f2c%hxc,f2c%hyc,f2c%kxy_c,f2c%kh2_c)
     
        end select

    end function Helm_A2hx

    function Ex(x,f2c) 
        !! Coarse-grid operation for two-level deflation
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind=realdp) :: Ex(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        type(Gridpara) :: CurrentGrid
        complex(kind=realdp), allocatable, dimension(:,:)  :: temp
        integer(kind=4)  :: re_A2h
            !! Flag for using Re-descretization approach (1) or not (0)

        CurrentGrid = CoarseGridpara(f2c)
        Ex   = czero
        re_A2h = 1
        if (re_A2h == 0) then !! Straight-forward Galerkin coarsening approach
            allocate(temp(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
            temp = czero
            temp = Zx(x,f2c)          !temp = Z x
            temp = Helm_Ahx(temp,f2c) !temp = A_h Z x
            Ex   = ZTx(temp,f2c)      !Ex = (Z^T A_h Z) x
            deallocate(temp)
        else
            if(M_flag == 4) then 
                Ex = Helm_A2hx(x,f2c)
                Ex = MGCSLP_invMHx(Ex,CurrentGrid)
                Ex = ZTZx(Ex,f2c)
            else
                Ex = Helm_A2hx(x,f2c)
            endif
        endif

        deallocate(CurrentGrid%wavenumber_k,CurrentGrid%wavenumber_kh_pow_2)

    end function Ex

    function invEy(y,f2c)
        !! Invert the coarse-grid operator for two-level delfation method, by using GMRES or Bi-CGSTAB
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: y(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind=realdp) :: invEy(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        invEy = czero
        if (def_mg_miter > 1) then !! determined by the specified maximum number of iterations on the coarse grid
            call DEF_fullgmres(y,invEy,f2c)
        else 
            call DEF_bicgstab(y,invEy,f2c)
        endif

    end function

    recursive function MultiLevel_invEy(y,f2c)
        !! A recursive function to invert the coarse-level operators for multilevel delfation methods, by using preconditioned FGMRES
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: y(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind=realdp) :: MultiLevel_invEy(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        integer(kind=4) :: level_i
        
        level_i=int(LOG2(dble((nx_global-1)/(f2c%nxc_global-1))))+1
            
        MultiLevel_invEy = czero
        if (level_i == 2 .and. Sv_L2 < 9.9d0) then
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[START]"
            endif
            call DEF_prefgmres(y,MultiLevel_invEy,f2c,level=level_i, rtol=dble(Sv_L2)*1.d-1)
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[FINISH]"
            endif
        elseif (level_i == 3 .and. Sv_L3 < 9.9d0) then
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[START]"
            endif
            call DEF_prefgmres(y,MultiLevel_invEy,f2c,level=level_i, rtol=dble(Sv_L3)*1.d-1)
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[FINISH]"
            endif
        elseif (level_i == 4 .and. Sv_L4 < 9.9d0) then
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[START]"
            endif
            call DEF_prefgmres(y,MultiLevel_invEy,f2c,level=level_i, rtol=dble(Sv_L4)*1.d-1)
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[FINISH]"
            endif
        else
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[START]"
            endif
            call DEF_prefgmres(y,MultiLevel_invEy,f2c,maximum_iterations=1,level=level_i)
            if (my_id == 0) then
                write(*,*) "Smoothing the coarser level", level_i, "........................[FINISH]"
            endif
        endif
    end function

    function Qx(x,f2c)
        !! Perform y=Qx in deflation definition, where Q=ZE^(-1)Z^T
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind = realdp) :: Qx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        complex(kind=realdp), allocatable, dimension(:,:)  :: temp

        allocate(temp(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        temp = czero
        Qx   = czero

        temp = ZTx(x,f2c) !temp = Z^T x
        temp = invEy(temp,f2c)  ! temp = E^{-1} Z^T x
        Qx   = Zx(temp,f2c) ! Qx = (Z E^{-1} Z^T) x

        deallocate(temp)

    end function Qx

    function Px(x,f2c)
        !! Perform y=Px in deflation definition, where P=I-AQ
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind = realdp) :: Px(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind=realdp), allocatable, dimension(:,:)  :: temp

        allocate(temp(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        temp = czero
        Px   = czero

        temp=Qx(x,f2c)
        Px = x - Helm_Ahx(temp,f2c)

        deallocate(temp)

    end function Px

    function P_DEFx(x,f2c)
        !! Deflation preconditioning, P = (I-AQ)+Q
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind = realdp) :: P_DEFx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        complex(kind=realdp), allocatable, dimension(:,:)  :: temp

        allocate(temp(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        temp = czero
        P_DEFx = czero

        temp=Qx(x,f2c)
        P_DEFx = x - Helm_Ahx(temp,f2c)+temp

        deallocate(temp)

    end function P_DEFx

    function P_ADEF1x(x,f2c)
        !! Adapted Deflation Preconditioning, P = M^(-1)(I-AQ)+Q, including higher-order deflation (ADP) if LAP > 1
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind = realdp) :: P_ADEF1x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        complex(kind=realdp), allocatable, dimension(:,:)  :: temp

        allocate(temp(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        temp = czero
        P_ADEF1x = czero

        temp=Qx(x,f2c)
        P_ADEF1x = x - Helm_Ahx(temp,f2c)
        P_ADEF1x = MGCSLP_invMx(P_ADEF1x)
        P_ADEF1x = P_ADEF1x + temp

        deallocate(temp)

    end function P_ADEF1x

    recursive function MultiLevelADP_Px(x,grid)
        !! Multilevel deflation preconditioning (MADP) 
        implicit none
        type(Gridpara), intent(inout) :: grid
        complex(kind = realdp) :: x(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)
        complex(kind = realdp) :: MultiLevelADP_Px(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP)
        
        type(TwoGrids) :: present_f2c
        complex(kind=realdp), allocatable, dimension(:,:)  :: tmp_Qx
        complex(kind=realdp), allocatable, dimension(:,:)  :: tmp_OnCoarse
        integer(kind=4) :: level_i
        
        level_i=int(LOG2(dble((nx_global-1)/(grid%nx_global-1))))+1

        if(level_i > def_nlevel - 1) then
            write(*,*) "The coarsening goes too far!!"
            stop
        else
            !coarsen
            if (mod(grid%nx_global-1,2) /= 0 .or. mod(grid%ny_global-1,2) /= 0) then
                write(*,*) "Coarsen exceed!"
                stop
            endif
            present_f2c = TwoGrids(grid)

            allocate(tmp_Qx(1-LAP:grid%nx+LAP,1-LAP:grid%ny+LAP))
            allocate(tmp_OnCoarse(1-LAP:present_f2c%nxc+LAP,1-LAP:present_f2c%nyc+LAP))
            tmp_Qx = czero
            tmp_OnCoarse = czero
            MultiLevelADP_Px = czero
            
            tmp_OnCoarse = ZTx(x,present_f2c) !Z^T x
            tmp_OnCoarse = MultiLevel_invEy(tmp_OnCoarse, present_f2c)  ! solve y = E^{-1} (Z^T x)
            tmp_Qx = Zx(tmp_OnCoarse,present_f2c) ! Z E^{-1} Z^T x

            MultiLevelADP_Px = x - Helm_Ax_nth(tmp_Qx,grid)
            
            if (level_i < 3) then 
                !! Multigrid-based CSLP on the finest and second level
                MultiLevelADP_Px = MGCSLP_invMHx(MultiLevelADP_Px, grid)
            else
                !! Krylov-based CSLP on the finest and second level
                MultiLevelADP_Px = KrylovCSLP_invMHx(MultiLevelADP_Px, grid)
            endif
            
            MultiLevelADP_Px = MultiLevelADP_Px + tmp_Qx

            deallocate(tmp_Qx, tmp_OnCoarse)
            deallocate(present_f2c%kxy_f,present_f2c%kh2_f,present_f2c%kxy_c,present_f2c%kh2_c)
        endif

    end function MultiLevelADP_Px

    function P_TLKMx(x,f2c)
        !! Two-Level Krylov Method, P = [(I-M^(-1)AQ`)+Q`]M^(-1), where Q' is defined based on M^(-1)A
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout) :: x(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)
        complex(kind = realdp) :: P_TLKMx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP)

        complex(kind=realdp), allocatable, dimension(:,:)  :: invMx,Q_invMx, invMA_Q_invMx

        allocate(invMx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        allocate(Q_invMx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        allocate(invMA_Q_invMx(1-LAP:f2c%nxf+LAP,1-LAP:f2c%nyf+LAP))
        invMx = czero
        Q_invMx = czero
        invMA_Q_invMx = czero

        P_TLKMx = czero
        invMx   = MGCSLP_invMx(x) !! M^(-1)x
        Q_invMx = Qx(invMx,f2c) !! QM^(-1)x

        invMA_Q_invMx = Helm_Ahx(Q_invMx,f2c) !! AQM^(-1)x
        invMA_Q_invMx = MGCSLP_invMx(invMA_Q_invMx) !! M^(-1)AQM^(-1)x

    
        P_TLKMx = invMx - invMA_Q_invMx + Q_invMx !! [(I-M^(-1)AQ)+Q]M^(-1)

        deallocate(invMx,Q_invMx,invMA_Q_invMx)

    end function P_TLKMx

    !================================================================================================
    subroutine DEF_fullgmres(y,x,f2c)
        !! A (CSLP preconditioned) GMRES solver for the coarse-grid problem in two-level deflation method
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: y(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind = realdp), intent(inout)  :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        
    
        ! Local arguments ------------------------------------------------------------
        type(Gridpara) :: CurrentGrid
        real(kind = realdp):: Rerror
        integer :: iter
        integer :: k,ki,j
        real(kind = realdp)                   :: b_norm, res_norm, Mres_norm, Mb_norm, MRerror
        complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
        complex(kind = realdp), allocatable, dimension(:,:)   :: res, Mres, Mb, H  !, u0
        complex(kind = realdp), allocatable, dimension(:,:,:) :: V
    
        ! Subroutine content ---------------------------------------------------------
        allocate(sn(def_mg_miter+1), cs(def_mg_miter+1), beta(def_mg_miter+1))
        allocate(res(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP), Mres(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(Mb(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(H(def_mg_miter+1,def_mg_miter))
        allocate(V(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP,def_mg_miter+1))
    
        CurrentGrid = CoarseGridpara(f2c)
        b_norm = mg_norm(y,f2c%nxc,f2c%nyc)
        
        if(M2h_flag == 1) then
            Mb = MGCSLP_invMHx(y,CurrentGrid)
        else
            Mb = y
        endif
        Mb_norm = mg_norm(Mb,f2c%nxc,f2c%nyc)
    
        res  = czero
        Mres = czero
        sn   = czero
        cs   = czero
        V    = czero
        H    = czero
        beta = czero

        res = y - Ex(x,f2c)  !r=b-Ax
        if(M2h_flag == 1) then
            Mres = MGCSLP_invMHx(res,CurrentGrid)
        else
            Mres = res
        endif

        Mres_norm = mg_norm(Mres,f2c%nxc,f2c%nyc)  ! ||r||
        MRerror   = Mres_norm / Mb_norm  ! scaled error?
        beta(1) = Mres_norm       ! beta(1)=||r0||

        V(:,:,1)  = Mres / Mres_norm  !  This is V(:,1) i.e. v1 in the algorithm
    
        k=0
        do j = 1,def_mg_miter
            k=k+1   !!Be careful!!, after the whole iteration without achieving eps, then the value of j will be "def_mg_miter+1".So we need a k.
            call def_arnoldi(V,H,k,f2c,CurrentGrid)
    
            call def_apply_givens_rotation(H, cs, sn, k)
    
            beta(k+1) = -sn(k)*beta(k)
            beta(k)   =  conjg(cs(k))*beta(k)
    
            MRerror = CDABS(beta(k+1))/Mb_norm
    
            if (MRerror < def_rtol) then
                exit
            end if
        enddo
    
        call def_back_substitute(H,beta,k)
    
        do ki = 1,k
            x = x + beta(ki)*V(:,:,ki)
        enddo

        res = y - Ex(x,f2c)
        res_norm = mg_norm(res,f2c%nxc,f2c%nyc)
        Rerror   = res_norm / b_norm

        if (my_id .eq. 0 ) then
         write(*,"(A,I9,A,E16.9)") "   Final DEF coarse solve GMRES Iter.    ", k, "    Error    ", Rerror
        endif
    
        deallocate(sn, cs, beta)
        deallocate(res,Mres)
        deallocate(Mb)
        deallocate(H)
        deallocate(V)
        deallocate(CurrentGrid%wavenumber_k,CurrentGrid%wavenumber_kh_pow_2)
    
        iter = k
    
    end subroutine DEF_fullgmres
    !====================================================================================
    subroutine def_arnoldi(V,H,k,f2c,CurrentGrid)
        !! Arnoldi prosess of DEF_fullgmres
        implicit none
    
        ! Subroutine arguments -------------------------------------------------------
        type(TwoGrids), intent(inout) :: f2c
        type(Gridpara) :: CurrentGrid
        complex(kind = realdp), dimension(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP,def_mg_miter+1), intent(inout) :: V
        complex(kind = realdp), dimension(def_mg_miter+1,def_mg_miter), intent(inout) :: H
        integer, intent(in) :: k

    
        ! Local arguments ------------------------------------------------------------
        integer :: i
    
        ! Subroutine content ---------------------------------------------------------
        V(:,:,k+1) = Ex(V(:,:,k),f2c)

        if (M2h_flag == 1) then
            V(:,:,k+1) = MGCSLP_invMHx(V(:,:,k+1),CurrentGrid)
        endif

        do i=1,k
    
            H(i, k)  = mg_dot_prod(V(:,:,i), V(:,:,k+1),f2c%nxc,f2c%nyc)     !! Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
            V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
    
        end do
    
        H(k+1, k) = mg_norm(V(:,:,k+1),f2c%nxc,f2c%nyc)
        V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)
    
    end subroutine def_arnoldi
    
    ! SUBROUTINE APPLY_GIVENS_ROTATION ====================================================
    subroutine def_apply_givens_rotation(H, cs, sn, k)
        !! Apply givens rotation of DEF_fullgmres
        implicit none
        
        ! Subroutine arguments -------------------------------------------------------
        integer, intent(in) :: k
        complex(kind = realdp), dimension(def_mg_miter+1,def_mg_miter), intent(inout) :: H
        complex(kind = realdp), dimension(def_mg_miter+1), intent(inout) :: cs, sn
        
        ! Local arguments ------------------------------------------------------------
        complex(kind = realdp)  :: temp
        integer  :: i
        
        ! Subroutine content ---------------------------------------------------------
        do i=1,k-1
            temp     =  conjg(cs(i))*H(i,k) + conjg(sn(i))*H(i+1,k)
            H(i+1,k) = -sn(i)*H(i,k) + cs(i)*H(i+1,k)
            H(i,k)   =  temp
        end do
        
        if (H(k,k)==czero) then
            cs(k) = 0.0d0
            sn(k) = 1.0d0
        else
            temp  = CDSQRT((CDABS(H(k,k)))**2 + (H(k+1,k))**2)
            cs(k) = H(k,k) / temp
            sn(k) = H(k+1,k) / temp
        end if
        
        H(k,k)   = conjg(cs(k))*H(k,k) + conjg(sn(k))*H(k+1,k)
        H(k+1,k) = (0.0d0,0.0d0)
    
    end subroutine def_apply_givens_rotation
    
    ! SUBROUTINE BACK_SUBSTITUTION =========================================================
    subroutine def_back_substitute(H,beta,k)
        !! Perform back substitute in of DEF_fullgmres
        implicit none
        integer :: i
        integer, intent(in) :: k
        complex(kind = realdp),    dimension(def_mg_miter+1,def_mg_miter), intent(in)    :: H
        complex(kind = realdp),    dimension(def_mg_miter+1),          intent(inout) :: beta
        
        beta(k) = beta(k)/H(k,k)
        
        do i=k-1,1,-1
            beta(i) = (beta(i) - sum(H(i,i+1:k)*beta(i+1:k)))/H(i,i)
        end do
    
    end subroutine def_back_substitute

    !Algorithm Bi-CSGTAB=========================================================================
    subroutine DEF_bicgstab(y,x,f2c)
        !! A (CSLP preconditioned) Bi-CGSTAB solver for the coarse-grid problem in two-level deflation method
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: y(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind = realdp), intent(inout)  :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
    
        integer :: iter
        type(Gridpara) :: CurrentGrid

        complex(kind = realdp), allocatable, dimension(:,:) :: res_hat, res, vi, p, invMp, s, invMs, t
        complex(kind=realdp)                                :: rho_cg, rho0_cg
        complex(kind=realdp)                                :: alpha_cg, omega_cg, beta_cg
        real(kind=realdp)                                   :: res_norm, b_norm       

        !------------------------END PARAMETER AND VARIABLE-----------------------------!  

        allocate(res(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(res_hat(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(vi(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(p(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(invMp(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(s(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(invMs(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(t(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))

        
        rho_cg   = (1.d0,0.d0)
        alpha_cg = (1.d0,0.d0)
        omega_cg = (1.d0,0.d0)  

        res   = czero
        vi    = czero
        p     = czero
        invMp = czero
        s     = czero
        invMs = czero
        t     = czero

        CurrentGrid = CoarseGridpara(f2c)
        res = y - Ex(x,f2c)
        res_hat = res

        res_norm = mg_norm(res,f2c%nxc,f2c%nyc)
        b_norm   = mg_norm(y,f2c%nxc,f2c%nyc)
        
        iter = 0 
        do while(res_norm /= 0.d0 .and. res_norm > def_rtol*b_norm .and. iter < 3000) 
            rho0_cg = rho_cg
            rho_cg  = mg_dot_prod(res_hat,res,f2c%nxc,f2c%nyc)
            if(rho_cg == czero) then
                write(*,*) "OMG, rho=0!"
                stop
            endif

            beta_cg = (rho_cg/rho0_cg)*(alpha_cg/omega_cg)
            
            p       = res + beta_cg * (p - omega_cg*vi)

            if (M2h_flag == 1) then
                invMp = MGCSLP_invMHx(p,CurrentGrid)
            else
                invMp = p
            endif

            vi = Ex(invMp,f2c)
            alpha_cg= rho_cg / mg_dot_prod(res_hat,vi,f2c%nxc,f2c%nyc)

            if(alpha_cg == czero) then
                write(*,*) "OMG, alpha_cg=0!"
                stop
            endif

            s       = res - alpha_cg*vi
            
            if (mg_norm(s,f2c%nxc,f2c%nyc) < def_rtol*b_norm) then
                x = x + alpha_cg*invMp

                res = y - Ex(x,f2c) 

                res_norm = mg_norm(res,f2c%nxc,f2c%nyc)

                Exit
            endif
            
            if (M2h_flag == 1) then
                invMs = MGCSLP_invMHx(s,CurrentGrid)
            else 
                invMs = s
            endif
            t = Ex(invMs,f2c)
        
            omega_cg  = mg_dot_prod(t,s,f2c%nxc,f2c%nyc) &
                    /mg_dot_prod(t,t,f2c%nxc,f2c%nyc)

            if(omega_cg == czero) then
                write(*,*) "OMG, omega_cg=0!"
                stop
            endif

            x   = x + alpha_cg*invMp + omega_cg*invMs
            
            res = s - omega_cg*t
            
            res_norm = mg_norm(res,f2c%nxc,f2c%nyc)

            iter = iter + 1

        end do !-------> END OF While LOOP

        if (my_id .eq. 0 ) then
            write(*,"(A,I9,A,E16.9)") "   Final DEF coarse solve BiCGSTAB Iter.    ", iter, "     Rel. res.=", res_norm/b_norm
        end if
        
        deallocate(res, vi ,p, invMp, s, invMs, t, res_hat)
        deallocate(CurrentGrid%wavenumber_k, CurrentGrid%wavenumber_kh_pow_2)

    end subroutine DEF_bicgstab     

    ! SUBROUTINE Flexible GMRES with right Deflated precondition ==================================================================
    recursive subroutine DEF_prefgmres(y,x,f2c,maximum_iterations,level,rtol)
        !! A recursive flexible GMRES with right deflation preconditioning
        implicit none

        ! Subroutine arguments -------------------------------------------------------
        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: y(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind = realdp), intent(inout)  :: x(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        
        ! optional arguments
        integer, optional, intent(in) :: maximum_iterations 
        integer, optional, intent(in) :: level
        real(kind = realdp), optional, intent(in) :: rtol

        ! Local arguments ------------------------------------------------------------
        type(Gridpara) :: CurrentGrid
        integer :: k,ki,i,j,iter,maxit,li
        real(kind = realdp) :: Rerror, res_tol 
        real(kind = realdp) :: b_norm, res_norm
        complex(kind = realdp)  :: temp
        complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
        complex(kind = realdp), allocatable, dimension(:,:)   :: res, H
        complex(kind = realdp), allocatable, dimension(:,:,:) :: V, Z

        ! Subroutine content ---------------------------------------------------------
        maxit = 55
        if ( present(maximum_iterations) ) maxit = maximum_iterations

        li = 2
        if ( present(level) ) li = level

        res_tol = 1.d-1
        if ( present(rtol) ) res_tol = rtol

        allocate(sn(maxit+1), cs(maxit+1), beta(maxit+1))
        allocate(res(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP))
        allocate(H(maxit+1,maxit))
        allocate(V(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP,maxit+1))
        allocate(Z(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP,maxit+1))

        CurrentGrid = CoarseGridpara(f2c)
        b_norm = mg_norm(y,f2c%nxc,f2c%nyc)

        res  = czero
        sn   = czero
        cs   = czero
        V    = czero
        H    = czero
        beta = czero

        res = y - Helm_Ax_nth(x,CurrentGrid)  !r=y-Ax

        res_norm = mg_norm(res,f2c%nxc,f2c%nyc)  ! ||r||
        Rerror   = res_norm / b_norm  

        beta(1) = res_norm       ! beta(1)=||r0||
        V(:,:,1)  = res / res_norm  !  This is V(:,1) i.e. v1 in the algorithm

        k=0
        do j=1,maxit
            k=k+1   !!Be careful!!, after the whole iteration without achieving eps, then the value of j will be "maxit+1".So we need a k.
            
            !call Prearnoldi(V, H, k)----------------------------------------
            if (li  < def_nlevel) then
                Z(:,:,k) = MultiLevelADP_Px(V(:,:,k),CurrentGrid)
            elseif (li  == def_nlevel) then
                Z(:,:,k) = KrylovCSLP_invMHx(V(:,:,k),CurrentGrid)
            else 
                write(*,*) "Warning!! Level exceeds!!!"
                stop
            endif
            V(:,:,k+1) = Helm_Ax_nth(Z(:,:,k), CurrentGrid)

            do i=1,k
                H(i, k)  = mg_dot_prod(V(:,:,i), V(:,:,k+1),f2c%nxc,f2c%nyc)     !! Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
                V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
            end do

            H(k+1, k) = mg_norm(V(:,:,k+1),f2c%nxc,f2c%nyc)
            V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)
            !-----------------------------------------------------------------

            !call apply_givens_rotation(H, cs, sn, k)-----------------------
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
            !-----------------------------------------------------------------


            beta(k+1) = -sn(k)*beta(k)
            beta(k)   =  conjg(cs(k))*beta(k)

            Rerror = CDABS(beta(k+1))/b_norm

            if (Rerror < res_tol) then
                exit
            end if

        end do

        !----call back_substitute(H,beta,k)---------------------------------
        beta(k) = beta(k)/H(k,k)
        do i=k-1,1,-1
            beta(i) = (beta(i) - sum(H(i,i+1:k)*beta(i+1:k)))/H(i,i)
        end do
        !------------------------------------------------------------------

        do ki = 1,k
            x = x + beta(ki)*Z(:,:,ki)
        enddo

        !calculate the final Res
        res = y - Helm_Ax_nth(x,CurrentGrid)  !r=y-Ax
        res_norm = mg_norm(res,f2c%nxc,f2c%nyc)
        Rerror   = res_norm / b_norm

        if (my_id .eq. 0 ) then
            write(*,"(A,I2,A,I5,A,E14.9)") "       Level", li, ": Final DEF coarse solve FGMRES Iter. ", k, "    Rres ", Rerror
        endif

        deallocate(sn, cs, beta)
        deallocate(res)
        deallocate(H)
        deallocate(V)
        deallocate(Z)
        deallocate(CurrentGrid%wavenumber_k,CurrentGrid%wavenumber_kh_pow_2)

        iter = k
    end subroutine DEF_prefgmres

    subroutine Helmholtz2d_O4cmpct(v_in,v_out,f2c)
        !! A coarse-grid Helmholtz operator for two-level deflation, descretized by a classic compact fourth-order FD scheme
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: v_in(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind = realdp), intent(inout)  :: v_out(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        
        real(kind = realdp) :: cmpctAs,cmpctAc
        integer :: i,j
        complex(kind = realdp) :: cmpctA0

        if (f2c%hxc /= f2c%hyc) then
            write(*,*) "careful, hx not equals to hy, in helmholtzOP"
        endif
        v_out = czero

        v_in(1-LAP:0,:) = czero
        v_in(:,1-LAP:0) = czero
        v_in(f2c%nxc+1:f2c%nxc+LAP,:) = czero
        v_in(:,f2c%nyc+1:f2c%nyc+LAP) = czero

        if (flag_BCs == 1) then
            if (npy == 0) then
                do i=1,f2c%nxc
                    v_in(i,0)=v_in(i,1)
                enddo
            endif
            if (npy == npy0-1) then
                do i=1,f2c%nxc
                    v_in(i,f2c%nyc+1)=v_in(i,f2c%nyc)
                enddo
            endif
            if (npx == 0) then
                do j=1,f2c%nyc
                    v_in(0,j)=v_in(1,j)
                enddo
            endif
            if (npx==npx0-1) then
                do j=1,f2c%nyc
                    v_in(f2c%nxc+1,j)=v_in(f2c%nxc,j)
                enddo
            endif

            if (npx == 0 .and. npy == 0) then
                v_in(0,0)=(v_in(0,1)+v_in(1,0))/2.d0
            endif

            if (npx == 0 .and. npy == npy0-1) then
                v_in(0,f2c%nyc+1)=(v_in(0,f2c%nyc)+v_in(1,f2c%nyc+1))/2.d0
            endif

            if (npx==npx0-1 .and. npy == 0) then
                v_in(f2c%nxc+1,0)=(v_in(f2c%nxc,0)+v_in(f2c%nxc+1,1))/2.d0
            endif

            if (npx==npx0-1 .and. npy == npy0-1) then
                v_in(f2c%nxc+1,f2c%nyc+1)=(v_in(f2c%nxc,f2c%nyc+1)+v_in(f2c%nxc+1,f2c%nyc))/2.d0
            endif

        elseif (flag_BCs == 2) then
            if (npy == 0) then
                do i=1,f2c%nxc
                    v_in(i,0)=2.d0*cone*f2c%kxy_c(i,1)*f2c%hyc*(1.d0-f2c%kh2_c(i,1)/6.d0)*v_in(i,1)+v_in(i,2)
                enddo
            endif
            if (npy == npy0-1) then
                do i=1,f2c%nxc
                    v_in(i,f2c%nyc+1)=2.d0*cone*f2c%kxy_c(i,f2c%nyc)*f2c%hyc* &
                                    (1.d0-f2c%kh2_c(i,f2c%nyc)/6.d0)*v_in(i,f2c%nyc)+v_in(i,f2c%nyc-1)
                enddo
            endif
            if (npx == 0) then
                do j=1,f2c%nyc
                    v_in(0,j)=2.d0*cone*f2c%kxy_c(1,j)*f2c%hxc*(1.d0-f2c%kh2_c(1,j)/6.d0)*v_in(1,j)+v_in(2,j)
                enddo
            endif
            if (npx==npx0-1) then
                do j=1,f2c%nyc
                    v_in(f2c%nxc+1,j)=2.d0*cone*f2c%kxy_c(f2c%nxc,j)*f2c%hxc* &
                                    (1.d0-f2c%kh2_c(f2c%nxc,j)/6.d0)*v_in(f2c%nxc,j)+v_in(f2c%nxc-1,j)
                enddo
            endif

            if (npx == 0 .and. npy == 0) then
                v_in(0,0)=(v_in(0,1)+v_in(1,0))/2.d0
            endif

            if (npx == 0 .and. npy == npy0-1) then
                v_in(0,f2c%nyc+1)=(v_in(0,f2c%nyc)+v_in(1,f2c%nyc+1))/2.d0
            endif

            if (npx==npx0-1 .and. npy == 0) then
                v_in(f2c%nxc+1,0)=(v_in(f2c%nxc,0)+v_in(f2c%nxc+1,1))/2.d0
            endif

            if (npx==npx0-1 .and. npy == npy0-1) then
                v_in(f2c%nxc+1,f2c%nyc+1)=(v_in(f2c%nxc,f2c%nyc+1)+v_in(f2c%nxc+1,f2c%nyc))/2.d0
            endif


        else
            write(*,*) "There is no such a boundary conditions yet!!!"
            stop
        endif
        !------------------------------------------------------------------------------------------------<

        call mg_check_xy2d(v_in, f2c%nxc, f2c%nyc)

        do j=1,f2c%nyc
            do i=1,f2c%nxc
                cmpctA0 = 10.d0/3.d0-f2c%kh2_c(i,j)*(2.d0/3.d0+1.d0/36.d0)
                cmpctAs = -2.d0/3.d0-f2c%kh2_c(i,j)*(1.d0/12.d0-1.d0/72.d0)
                cmpctAc = -1.d0/6.d0-f2c%kh2_c(i,j)*(1.d0/144.d0)
                v_out(i,j) = cmpctA0*v_in(i,j)+cmpctAs*(v_in(i,j+1)+v_in(i+1,j)+v_in(i,j-1)+v_in(i-1,j)) &
                            +cmpctAc*(v_in(i+1,j+1)+v_in(i+1,j-1)+v_in(i-1,j-1)+v_in(i-1,j+1))
            enddo
        enddo

        v_out=v_out/f2c%hxhyc

        if (flag_BCs == 1) then
            if (npy == 0) then
                v_out(:,1)=v_in(:,1)
            endif
            if (npy == npy0-1) then
                v_out(:,f2c%nyc)=v_in(:,f2c%nyc)
            endif
            if (npx == 0) then
                v_out(1,:)=v_in(1,:)
            endif
            if (npx==npx0-1) then
                v_out(f2c%nxc,:)=v_in(f2c%nxc,:)
            endif 
        endif

    end subroutine Helmholtz2d_O4cmpct

    subroutine Helmholtz2d_ReD_Glk(v_in,v_out,f2c)
        !! A coarse-grid Helmholtz operator for two-level deflation, descretized by ReD-GLK FD scheme
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: v_in(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        complex(kind = realdp), intent(inout)  :: v_out(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)
        
        real(kind = realdp) :: ae,aw,an,as
        integer :: i,j,ia,ja
        complex(kind = realdp) :: ap,a(-2:2,-2:2),b(-2:2,-2:2)

        if (f2c%hxc /= f2c%hyc) then
            write(*,*) "careful, hx not equals to hy, in helmholtzOP"
        endif
        v_out = czero

        v_in(1-LAP:0,:) = czero
        v_in(:,1-LAP:0) = czero
        v_in(f2c%nxc+1:f2c%nxc+LAP,:) = czero
        v_in(:,f2c%nyc+1:f2c%nyc+LAP) = czero

        call mg_check_xy2d(v_in, f2c%nxc, f2c%nyc)
        call ExtrpltGhostBCs(v_in,f2c) !! fill in a layer of ghost grid points by boundary conditions

        do j=1,f2c%nyc
            do i=1,f2c%nxc
                call ReD_Glk_stencils(a,b,f2c,i,j)
                !==============================
                do ja=-2,2
                    do ia=-2,2
                        v_out(i,j)=v_out(i,j)+a(ia,ja)*v_in(i+ia,j+ja)
                    enddo
                enddo
            enddo
        enddo

        if (flag_BCs == 1) then

            if (npy == 0) then
                v_out(:,1)=v_in(:,1)*f2c%hxhyc
            endif
            if (npy == npy0-1) then
                v_out(:,f2c%nyc)=v_in(:,f2c%nyc)*f2c%hxhyc
            endif
            if (npx == 0) then
                v_out(1,:)=v_in(1,:)*f2c%hxhyc
            endif
            if (npx==npx0-1) then
                v_out(f2c%nxc,:)=v_in(f2c%nxc,:)*f2c%hxhyc
            endif 

        elseif (flag_BCs ==2) then
            !! The boundary keeps 2nd-order finite difference scheme using five points stencils 
            if (npy == 0) then
                j = 1
                do i=1,f2c%nxc
                    call Helmholtz2d_stencils(ap,an,as,aw,ae,f2c%kh2_c(i,j))
                    !==============================
                    v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                            +an*v_in(i,j+1)+as*v_in(i,j-1)
                enddo
            endif

            if (npy == npy0-1) then
                j = f2c%nyc
                    do i=1,f2c%nxc
                        call Helmholtz2d_stencils(ap,an,as,aw,ae,f2c%kh2_c(i,j))
                        !==============================
                        v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                                +an*v_in(i,j+1)+as*v_in(i,j-1)
                    enddo
            endif

            if (npx == 0) then
                i = 1
                    do j=1,f2c%nyc
                        call Helmholtz2d_stencils(ap,an,as,aw,ae,f2c%kh2_c(i,j))
                        !==============================
                        v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                                +an*v_in(i,j+1)+as*v_in(i,j-1)
                    enddo
            endif

            if (npx==npx0-1) then
                i = f2c%nxc
                    do j=1,f2c%nyc
                        call Helmholtz2d_stencils(ap,an,as,aw,ae,f2c%kh2_c(i,j))
                        !==============================
                        v_out(i,j)=ap*v_in(i,j)+ae*v_in(i+1,j)+aw*v_in(i-1,j) &
                                                +an*v_in(i,j+1)+as*v_in(i,j-1)
                    enddo
            endif
        endif

        v_out=v_out/f2c%hxhyc

    end subroutine Helmholtz2d_ReD_Glk

    subroutine ReD_Glk_stencils(a,b,f2c,ic,jc)
        !! *ReD-Glk* computational stencils of the Helmholtz operator for coarse-grid level in two-deflation method
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp) :: a(-2:2,-2:2)
            !! a: Laplace operator
        complex(kind = realdp) :: b(-2:2,-2:2)
            !! b: wavenumber operator
        integer :: ic,jc,ib,jb

        a(:,-2)=[complex(realdp) :: -0.0117187500d0, -0.17187500d0, -0.382812500d0, -0.17187500d0, -0.0117187500d0 ]
        a(:,-1)=[complex(realdp) :: -0.1718750000d0, -0.43750000d0,  0.218750000d0, -0.43750000d0, -0.1718750000d0 ]
        a(:,0) =[complex(realdp) :: -0.3828125000d0,  0.21875000d0,  3.828125000d0,  0.21875000d0, -0.3828125000d0 ]
        a(:,1) =[complex(realdp) :: -0.1718750000d0, -0.43750000d0,  0.218750000d0, -0.43750000d0, -0.1718750000d0 ]
        a(:,2) =[complex(realdp) :: -0.0117187500d0, -0.17187500d0, -0.382812500d0, -0.17187500d0, -0.0117187500d0 ]

        b(:,-2)=[complex(realdp) :: 0.000244140625d0, 0.006835937500d0, 0.017089843750d0, 0.006835937500d0, 0.000244140625d0 ]
        b(:,-1)=[complex(realdp) :: 0.006835937500d0, 0.191406250000d0, 0.478515625000d0, 0.191406250000d0, 0.006835937500d0 ]
        b(:,0) =[complex(realdp) :: 0.017089843750d0, 0.478515625000d0, 1.196289062500d0, 0.478515625000d0, 0.017089843750d0 ]
        b(:,1) =[complex(realdp) :: 0.006835937500d0, 0.191406250000d0, 0.478515625000d0, 0.191406250000d0, 0.006835937500d0 ]
        b(:,2) =[complex(realdp) :: 0.000244140625d0, 0.006835937500d0, 0.017089843750d0, 0.006835937500d0, 0.000244140625d0 ]

        do jb=-2,2
            do ib=-2,2
                b(ib,jb) = b(ib,jb)*f2c%kh2_c(ic+ib,jc+jb)
            enddo
        enddo
        a=a-b
        
    end subroutine ReD_Glk_stencils

    subroutine ExtrpltGhostBCs(v_in,f2c)
        !! Extrapolate a layer of ghost grid points based on the boundary conditions 
        implicit none

        type(TwoGrids), intent(inout) :: f2c
        complex(kind = realdp), intent(inout)  :: v_in(1-LAP:f2c%nxc+LAP,1-LAP:f2c%nyc+LAP)

        integer :: i,j
        
        if (flag_BCs == 1) then
            if (npx == 0) then
                do j=1-LAP,f2c%nyc+LAP
                    v_in(0,j)=2.d0*v_in(1,j)-v_in(2,j)
                    f2c%kxy_c(0,j)=0.d0 !! The zero padding of wavenumber is ok in practical
                    f2c%kh2_c(0,j)=f2c%kxy_c(0,j)*f2c%hxhyc
                enddo
            endif
            if (npx==npx0-1) then
                do j=1-LAP,f2c%nyc+LAP
                    v_in(f2c%nxc+1,j)=2.d0*v_in(f2c%nxc,j)-v_in(f2c%nxc-1,j)
                    f2c%kxy_c(f2c%nxc+1,j)=0.d0
                    f2c%kh2_c(f2c%nxc+1,j)=f2c%kxy_c(f2c%nxc+1,j)*f2c%hxhyc
                enddo
            endif
            if (npy == 0) then
                do i=1-LAP,f2c%nxc+LAP
                    v_in(i,0)=2.d0*v_in(i,1)-v_in(i,2)
                    f2c%kxy_c(i,0)=0.d0
                    f2c%kh2_c(i,0)=f2c%kxy_c(i,0)*f2c%hxhyc
                enddo
            endif
            if (npy == npy0-1) then
                do i=1-LAP,f2c%nxc+LAP
                    v_in(i,f2c%nyc+1)=2.d0*v_in(i,f2c%nyc)-v_in(i,f2c%nyc-1)
                    f2c%kxy_c(i,f2c%nyc+1)=0.d0
                    f2c%kh2_c(i,f2c%nyc+1)=f2c%kxy_c(i,f2c%nyc+1)*f2c%hxhyc
                enddo
            endif

            !=THE EDGES=
            if (npx == 0 .and. npy == 0) then
                v_in(0,0)=4.d0*v_in(1,1)-v_in(2,2)-v_in(0,2)-v_in(2,0)
                f2c%kxy_c(0,0)=0.d0
                f2c%kh2_c(0,0)=f2c%kxy_c(0,0)*f2c%hxhyc
            endif

            if (npx == 0 .and. npy == npy0-1) then
                v_in(0,f2c%nyc+1)=4.d0*v_in(1,f2c%nyc)-v_in(2,f2c%nyc-1)-v_in(2,f2c%nyc+1)-v_in(0,f2c%nyc-1)
                f2c%kxy_c(0,f2c%nyc+1)=0.d0
                f2c%kh2_c(0,f2c%nyc+1)=f2c%kxy_c(0,f2c%nyc+1)*f2c%hxhyc
            endif

            if (npx==npx0-1 .and. npy == 0) then
                v_in(f2c%nxc+1,0)=4.d0*v_in(f2c%nxc,1)-v_in(f2c%nxc-1,2)-v_in(f2c%nxc-1,0)-v_in(f2c%nxc+1,2)
                f2c%kxy_c(f2c%nxc+1,0)=0.d0
                f2c%kh2_c(f2c%nxc+1,0)=f2c%kxy_c(f2c%nxc+1,0)*f2c%hxhyc
            endif

            if (npx==npx0-1 .and. npy == npy0-1) then
                v_in(f2c%nxc+1,f2c%nyc+1)=4.d0*v_in(f2c%nxc,f2c%nyc)-v_in(f2c%nxc-1,f2c%nyc-1) &
                                              -v_in(f2c%nxc-1,f2c%nyc+1)-v_in(f2c%nxc+1,f2c%nyc-1)
                f2c%kxy_c(f2c%nxc+1,f2c%nyc+1)=0.d0
                f2c%kh2_c(f2c%nxc+1,f2c%nyc+1)=f2c%kxy_c(f2c%nxc+1,f2c%nyc+1)*f2c%hxhyc
            endif

        elseif (flag_BCs == 2) then
            if (npx == 0) then
                do j=1-LAP,f2c%nyc+LAP
                    v_in(0,j)=2.d0*cone*f2c%kxy_c(1,j)*f2c%hxc*v_in(1,j)+v_in(2,j)
                    f2c%kxy_c(0,j)=0.d0
                    f2c%kh2_c(0,j)=f2c%kxy_c(0,j)*f2c%hxhyc
                enddo
            endif
            if (npx==npx0-1) then
                do j=1-LAP,f2c%nyc+LAP
                    v_in(f2c%nxc+1,j)=2.d0*cone*f2c%kxy_c(f2c%nxc,j)*f2c%hxc* &
                                    v_in(f2c%nxc,j)+v_in(f2c%nxc-1,j)
                    f2c%kxy_c(f2c%nxc+1,j)=0.d0
                    f2c%kh2_c(f2c%nxc+1,j)=f2c%kxy_c(f2c%nxc+1,j)*f2c%hxhyc
                enddo
            endif
            if (npy == 0) then
                do i=1-LAP,f2c%nxc+LAP
                    v_in(i,0)=2.d0*cone*f2c%kxy_c(i,1)*f2c%hyc*v_in(i,1)+v_in(i,2)
                    f2c%kxy_c(i,0)=0.d0
                    f2c%kh2_c(i,0)=f2c%kxy_c(i,0)*f2c%hxhyc
                    ! !!Another Ghost point seems no gains
                enddo
            endif
            if (npy == npy0-1) then
                do i=1-LAP,f2c%nxc+LAP
                    v_in(i,f2c%nyc+1)=2.d0*cone*f2c%kxy_c(i,f2c%nyc)*f2c%hyc* &
                                    v_in(i,f2c%nyc)+v_in(i,f2c%nyc-1)
                    f2c%kxy_c(i,f2c%nyc+1)=0.d0
                    f2c%kh2_c(i,f2c%nyc+1)=f2c%kxy_c(i,f2c%nyc+1)*f2c%hxhyc
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

                v_in(0,0)=((2.d0*cone*f2c%kxy_c(1,0)*f2c%hxc*v_in(1,0)+v_in(2,0)) &
                          +(2.d0*cone*f2c%kxy_c(0,1)*f2c%hyc*v_in(0,1)+v_in(0,2)))/2.d0
                f2c%kxy_c(0,0)=0.d0
                f2c%kh2_c(0,0)=f2c%kxy_c(0,0)*f2c%hxhyc
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

                v_in(0,f2c%nyc+1)=((2.d0*cone*f2c%kxy_c(1,f2c%nyc+1)*f2c%hxc*v_in(1,f2c%nyc+1)+v_in(2,f2c%nyc+1)) &
                                  +(2.d0*cone*f2c%kxy_c(0,f2c%nyc)*f2c%hyc*v_in(0,f2c%nyc)+v_in(0,f2c%nyc-1)))/2.d0
                f2c%kxy_c(0,f2c%nyc+1)=0.d0
                f2c%kh2_c(0,f2c%nyc+1)=f2c%kxy_c(0,f2c%nyc+1)*f2c%hxhyc

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

                v_in(f2c%nxc+1,0)=((2.d0*cone*f2c%kxy_c(f2c%nxc,0)*f2c%hxc*v_in(f2c%nxc,0)+v_in(f2c%nxc-1,0)) &
                                  +(2.d0*cone*f2c%kxy_c(f2c%nxc+1,1)*f2c%hyc*v_in(f2c%nxc+1,1)+v_in(f2c%nxc+1,2)))/2.d0
                f2c%kxy_c(f2c%nxc+1,0)=0.d0
                f2c%kh2_c(f2c%nxc+1,0)=f2c%kxy_c(f2c%nxc+1,0)*f2c%hxhyc

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

                v_in(f2c%nxc+1,f2c%nyc+1)=((2.d0*cone*f2c%kxy_c(f2c%nxc,f2c%nyc+1)*f2c%hxc*v_in(f2c%nxc,f2c%nyc+1)&
                                        +v_in(f2c%nxc-1,f2c%nyc+1)) &
                                        +(2.d0*cone*f2c%kxy_c(f2c%nxc+1,f2c%nyc)*f2c%hyc*v_in(f2c%nxc+1,f2c%nyc)&
                                        +v_in(f2c%nxc+1,f2c%nyc-1)))/2.d0
                f2c%kxy_c(f2c%nxc+1,f2c%nyc+1)=0.d0
                f2c%kh2_c(f2c%nxc+1,f2c%nyc+1)=f2c%kxy_c(f2c%nxc+1,f2c%nyc+1)*f2c%hxhyc
                
            endif
        else
            write(*,*) "There is no such a boundary conditions yet!!!"
            stop
        endif
        
    end subroutine ExtrpltGhostBCs

end module deflaion_setup
