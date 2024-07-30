module read_setup
    !! This is a module to read the input parameters for the solver settings
    use mpi
    use comm_variable

    implicit none
    integer, parameter :: NparaMax = 50          
        !! the maximun number of integer input parameter
    integer, parameter :: RparaMax = 50          
        !! the maximun number of real input parameter

contains
    subroutine read_parameter()
        !! read the input parameters for the solver settings
        implicit none
        integer,             allocatable, dimension(:) :: Nparameters
            !! Integer input parameter
        real(kind = realdp), allocatable, dimension(:) :: Rparameters
            !! Real input parameter

        character(len=64) :: arg
            !! Input filename

        allocate(Nparameters(1:NparaMax))
        allocate(Rparameters(1:RparaMax))

        Nparameters=0
        Rparameters=0.d0

        !! Read by Rank 0 and then broadcast to the other ranks
        if (my_id .eq. 0) then
            call get_command_argument(1, arg)
                !! Get the input filename

            if (LEN_TRIM(arg) == 0) then
                !! The Default input filename is Input/Helmholtz.in
                arg='Input/Helmholtz.in'
            endif

            open(100,file=trim(arg))

            read(100,*)
            read(100,*)
            read(100,*)nx_global,ny_global,LAP
            read(100,*)
            read(100,*)npx0,npy0
            read(100,*)
            read(100,*)slx,sly
            read(100,*)
            read(100,*)Iperiodic_X, Iperiodic_Y
            read(100,*)
            read(100,*)i_case, flag_BCs
            read(100,*)
            read(100,*)m_iter,Algorithm,Irestart,eps
            read(100,*)
            read(100,*)k_case, k0, freq
            read(100,*)
            read(100,*)M_flag, M2h_flag, beta1, beta2
            read(100,*)
            read(100,*)MG_flag, def_mg_miter, cslp_mg_miter, nx_min
            read(100,*)
            read(100,*)def_rtol, cslp_mg_tol
            read(100,*)
            read(100,*)def_nlevel, Sv_L2, Sv_L3, Sv_L4

            close(100)
        endif

        Nparameters(1) = nx_global
        Nparameters(2) = ny_global
        Nparameters(3) = npx0
        Nparameters(4) = npy0
        Nparameters(5) = Iperiodic_X
        Nparameters(6) = Iperiodic_Y
        Nparameters(7) = m_iter
        Nparameters(8) = Irestart
        Nparameters(9) = M_flag
        Nparameters(10)= MG_flag
        Nparameters(11)= def_mg_miter
        Nparameters(12)= nx_min
        Nparameters(13)= i_case
        Nparameters(14)= flag_BCs
        Nparameters(15)= Algorithm
        Nparameters(16)= k_case
        Nparameters(17)= M2h_flag
        Nparameters(18)= cslp_mg_miter
        Nparameters(19)= LAP
        Nparameters(20)= def_nlevel
        Nparameters(21)= Sv_L2
        Nparameters(22)= Sv_L3
        Nparameters(23)= Sv_L4

        call MPI_BCAST(Nparameters,NparaMax,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        Rparameters(1) = slx
        Rparameters(2) = sly
        Rparameters(3) = eps
        Rparameters(4) = k0
        Rparameters(5) = beta1
        Rparameters(6) = beta2
        Rparameters(7) = freq
        Rparameters(8) = def_rtol
        Rparameters(9) = cslp_mg_tol

        call MPI_BCAST(Rparameters,RparaMax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        nx_global       = Nparameters(1)
        ny_global       = Nparameters(2)
        npx0            = Nparameters(3)
        npy0            = Nparameters(4)
        Iperiodic_X     = Nparameters(5)
        Iperiodic_Y     = Nparameters(6)
        m_iter          = Nparameters(7)
        Irestart        = Nparameters(8)
        M_flag          = Nparameters(9)
        MG_flag         = Nparameters(10)
        def_mg_miter    = Nparameters(11)
        nx_min          = Nparameters(12)
        i_case          = Nparameters(13)
        flag_BCs        = Nparameters(14)
        Algorithm       = Nparameters(15)
        k_case          = Nparameters(16)
        M2h_flag        = Nparameters(17)
        cslp_mg_miter   = Nparameters(18)
        LAP             = Nparameters(19)
        def_nlevel      = Nparameters(20)
        Sv_L2           = Nparameters(21)
        Sv_L3           = Nparameters(22)
        Sv_L4           = Nparameters(23)

        slx             = Rparameters(1)
        sly             = Rparameters(2)
        eps             = Rparameters(3)
        k0              = Rparameters(4)
        beta1           = Rparameters(5)
        beta2           = Rparameters(6)
        freq            = Rparameters(7)
        def_rtol        = Rparameters(8)
        cslp_mg_tol     = Rparameters(9)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        deallocate(Nparameters, Rparameters)

    end subroutine read_parameter

end module read_setup