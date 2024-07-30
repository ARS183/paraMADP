program helmholtz_2d
  !! Parallel proconditioned Krylov slover for 2D Helmholtz equation
  use mpi
  use comm_variable
  use read_setup
  use mpi_setup
  use define_grid
  use define_rhs
  use analytical_sol
  use wavenumber
  use solvers
  use write_data
  use user_module
  use idrs_module
  

  implicit none
  !Arguements
  type(matrix)                  :: A
    !! An empty matrix type for calling idrs solver
  type(preconditioner)          :: M1
    !! a preconditioner type for calling idrs solver. Since the IDR(s) solver part is ported, its input and output parameters are somewhat different from other solvers.
  complex(kind = realdp), allocatable,dimension(:,:) :: u
    !! arrays for the solution
  complex(kind = realdp), allocatable,dimension(:,:) :: u_ex
    !! arrays for exact solution
  complex(kind = realdp), allocatable,dimension(:,:) :: u_err
    !! arrays for the error of solution
  complex(kind = realdp), allocatable,dimension(:,:) :: b
    !! arrays for the right-hand side
  real(kind = realdp),    allocatable,dimension(:,:) :: xx,yy
    !! coordinates
  real(kind=realdp), allocatable :: resvec(:) 
    !! Residual vector for idrs solver
  real(kind = realdp) :: time_start  
    !! time_start: variable for runtime measurement start  
  real(kind = realdp) :: time_end  
    !! time_end: variable for runtime measurement end  
  real(kind = realdp) :: u_err_max  
    !! u_err_max: local max|u_err|=max|u-u_ex|  
  real(kind = realdp) :: Aerror  
    !! Aerror: global max|u_err|=max|u-u_ex|  
  real(kind = realdp) :: Rerror  
    !! Rerror: relative residual ||b-Au||/||b||  
  real(kind = realdp) :: max_k_local  
    !! max_k_local: local maximum wavenumber  
  real(kind = realdp) :: max_k_global  
    !! max_k_global: global maximum wavenumber
  integer :: iter,k,i,j
    !! number of iterations and loop indexs
  integer :: flag
    !! an output flag for idrs solvers

  ! Initilize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np_size,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

  ! read parameters from Input/Helmtoltz.in
  call read_parameter()

  ! Assess if uniform grid size
  if (i_case == 3 .or. i_case == 3) then
    if (dble(ny_global-1)/dble(nx_global-1) /= sly/slx) then
      if (my_id == 0) then
        write(*,*) 'I cannot handle non-equal gird so far!!!'
      endif
      stop
    endif
  endif

  ! compute grid size
  hy   = sly/dble(ny_global-1)
  hx   = slx/dble(nx_global-1)
  hxhy = hx*hy

  ! 2D domain patition and MPI setting
  call part2d()
  ! variables alllocation and initialization
  allocate(u(1-LAP:nx+LAP,1-LAP:ny+LAP),b(1-LAP:nx+LAP,1-LAP:ny+LAP), u_ex(1-LAP:nx+LAP,1-LAP:ny+LAP))
  allocate(xx(1-LAP:nx+LAP,1-LAP:ny+LAP), yy(1-LAP:nx+LAP,1-LAP:ny+LAP))
  allocate(u_err(1:nx,1:ny))
  allocate(resvec(m_iter+1))
  u     = czero
  b     = czero
  u_ex  = czero

  xx        = 0.d0
  yy        = 0.d0
  Rerror    = 0.d0
  u_err     = 0.d0
  u_err_max = 0.d0
  Aerror    = 0.d0
  iter  = 0

  M1%ni = nx
  M1%nj = ny
  M1%hx_c = hx
  M1%hy_c = hy

  ! define 2D coordinates for different model problems 
  if (i_case == 1) then
    call define_uniform_grid(xx,yy)
  elseif (i_case == 2) then
    call define_uniform_grid(xx,yy)
  elseif (i_case == 3) then
    call define_wedge_grid(xx,yy)
  elseif (i_case == 4) then
    call define_marmousi_grid(xx,yy)
  else
    write(*,*) "No such a case so far!"
    stop
  endif
  
  ! Analytical Solutions for two model problems
  if (i_case == 1) then
    call exact_2DCloseOff(u_ex,xx,yy)
  elseif (i_case == 2) then
    if (flag_BCs ==1) then
      ! Be careful that it will take too much time for large grid size
      call exact_2DPointSource_Dirichlet(u_ex,xx,yy)
    elseif (flag_BCs == 2) then
      call exact_2DPointSource_1stSomer(u_ex,xx,yy)
    else
      write(*,*) "No such a BC!"
      stop
    endif
  endif

  ! Define right-hand side 
  if (i_case == 1) then
    call RHS_2DCloseOff(b,xx,yy)
  elseif (i_case == 2) then
    call RHS_CenterSource2D(b,xx,yy)
  elseif (i_case == 3) then
      call RHS_2DWedge(b,xx,yy)
  elseif (i_case == 4) then
      call RHS_marmousi(b,xx,yy)
  else
    write(*,*) "No such a case so far!"
    stop
  endif

  ! Determine constant-wavenumber or non-constant wavenumber problems
  if (k_case == 0) then
    call Const_K()
  else
    if (i_case == 3) then
      call wavenumber_k_Wedge(xx,yy)
    elseif (i_case == 4) then
      call read_wavenumber_k_marmousi()
    endif
  endif

  ! Determine maximum wavenumber 
  max_k_local = maxval(wavenumber_k(1:nx,1:ny))
  indf_level = 2
  call MPI_ALLREDUCE(max_k_local, max_k_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  ! Calculate the indefinite coarse level approxiamtely 
  indf_level = int(LOG2(1.5d0/(max_k_global*hx)))+1
  if (my_id == 0) then
    write(*,*) "kh = ", max_k_global*hx, "indf_level = ", indf_level
  end if

  ! Assess the input number of levels for multilevel defiation 
  if (M_flag == 5) then
    !Samller shift based on dimensionless wavenumber. Uncomment if only Krylov iterations are used for CSLP approxiamtion 
    !beta2 = -1.d0/sqrt(max_k_global*max_k_global*slx*sly)
    if (mod(nx_global-1,2**(def_nlevel-1)) /= 0 .or. mod(ny_global-1,2**(def_nlevel-1)) /= 0) then
      write(*,*) "Cannot go to such a deep level!"
    endif
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Creat an output file ended with iters.plt to record the residual and iterations
  if (my_id == 0 ) then
    call system('mkdir -p Output')
    write(logname,"('Output/MP',I1,'BC',I1,'nx',I5.5,'k',I5.5,     &
                   & 'thm',I2.2,'pre',I1,'nldef',I1, &
                   & 'cslptol',I2.2,'deftol',I2.2,'np',I3.3,'iters.plt')") &
                   & i_case,flag_BCs,nx_global,int(k0),            &
                   & Algorithm,M_flag,def_nlevel,     &
                   & int(-LOG10(cslp_mg_tol)),int(-LOG10(def_rtol)),npx0*npy0
    open(1234, file=trim(logname), status='REPLACE')
    write(1234,*) 'variables=iter, res'
    close(1234)
  end if

  ! Start to call Krylov subspace solvers
  time_start=MPI_WTIME()
  ! GMRES 
  if (Algorithm == 1) then
    if (M_flag == 0) then ! no precondition
      if (Irestart == 1) then
        call restartgmres(b,u,Rerror,iter)
      else
        call fullgmres(b,u,Rerror,iter)
      endif
    else !left precondition GMRES
      if (Irestart == 1) then
        call Pre_restartgmres(b,u,Rerror,iter)
      else
        call Pre_fullgmres(b,u,Rerror,iter)
      endif
    endif 
    
  ! right preconditioned GMRES
  elseif (Algorithm == 2) then
    call full_pgmres(b,u,Rerror,iter)
  ! preconditioned fgmres
  elseif (Algorithm == 3) then
    call pfgmres(b,u,Rerror,iter)
  ! preconditioned GCR
  elseif (Algorithm == 4) then
    call full_pgcr(b,u,Rerror,iter)
  ! IDR(s) 
  elseif (Algorithm == 21) then
    u = idrs( A, b, M1, 1, eps, m_iter, 'biortho ', flag, Rerror, iter, resvec=resvec )
  elseif (Algorithm == 24) then
    u = idrs( A, b, M1, 4, eps, m_iter, 'biortho ', flag, Rerror, iter, resvec=resvec )
  elseif (Algorithm == 28) then
    u = idrs( A, b, M1, 8, eps, m_iter, 'biortho ', flag, Rerror, iter, resvec=resvec )
  elseif (Algorithm == 31) then
    u = idrs( A, b, M1, 1, eps, m_iter, 'bicgstab', flag, Rerror, iter, resvec=resvec )
  endif

  time_end=MPI_WTIME()
  ! Exit from solver
  ! write the number of iterations, relative residual and runtime
  if (my_id == 0 ) then
    write(*,"('Iterations  ', I6, '[Exit!!!]')") iter
    write(*,"('Relative residual=',E16.9)") Rerror
    write(*,"('Time=',E16.9)") time_end-time_start
  end if

  ! If idrs, write the residual verus iterations.
  if (Algorithm == 21 .or. Algorithm == 24 .or. Algorithm == 28 .or. Algorithm == 31) then
    if (my_id == 0 ) then
      open(1234, file=trim(logname),position='APPEND',status='OLD')
      do k = 1, iter
        write(1234, "(I9,E17.9)") k, resvec(k)
      enddo
      close(1234)
    end if
  endif

  ! Calculate the absolute error max|u-u_ex|
  ! (if there is no analytical solutions, u_ex=0)
  u_err=u_ex(1:nx,1:ny)-u(1:nx,1:ny)
  u_err_max=maxval(abs(u_err))

  call MPI_Reduce(u_err_max,           &
                  Aerror,              &
                  1,                   &
                  MPI_DOUBLE_PRECISION,&
                  MPI_MAX,             &
                  0,                   &
                  mpi_comm_world,      &
                  ierr)

  ! Output a .dat file, including the number of iterations, runtime, absolute error and relative residual
  if(my_id == 0) then
    write(output_name,"('Output/MP',I1,'BC',I1,'nx',I5.5,'k',I5.5,     &
                       & 'thm',I2.2,'pre',I1,'nldef',I1, &
                       & 'cslptol',I2.2,'deftol',I2.2,'np',I3.3,'.dat')")&
                       & i_case,flag_BCs,nx_global,int(k0),            &
                       & Algorithm,M_flag, def_nlevel,      &
                       & int(-LOG10(cslp_mg_tol)),int(-LOG10(def_rtol)),npx0*npy0
    open(12,file=trim(output_name),status='unknown')
    write(12,"('Total Iterations=',I6)") iter
    write(12,"('Time=',E16.9)") time_end-time_start
    write(12,"('Absolute Error=',E16.9)") Aerror
    write(12,"('Relative Error=',E16.9)") Rerror
    close(12)
  endif
  
  ! Write the solution
  call write_data_whole(xx,yy,u)

  !The field of wavenumber can be obtained as follow
  ! call write_real_data_whole(xx,yy,wavenumber_k)

  ! Deallocate wavenumber 
  call wavenumber_k_destroy()
  ! Deallocate the variables
  deallocate(u,b,u_ex,xx,yy,u_err)  
  ! Finalize MPI
  call mpi_finalize(ierr)

end program helmholtz_2d
