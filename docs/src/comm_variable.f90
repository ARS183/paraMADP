module comm_variable
   !! This is a module that define some common variables for the whole project
   use mpi
   implicit none

   !=====================math parameter===================
   integer, parameter :: realdp = 8             
      !! precision for the whole program
   real(kind=realdp), parameter :: pi = 3.141592653589793238d0
   real(kind=realdp), parameter :: e = 2.7182818285d0


   !========= program parameter=================
   integer, parameter :: npMax = 400 
      !! the maximun number of partitions in one dimension

   complex(kind = realdp), parameter :: cone  = (0.d0,1.d0)
   complex(kind = realdp), parameter :: czero = (0.d0,0.d0)
      !! Constants for complex numbers

   character (len = 100) :: filename
   character (len = 100) :: output_name
   character (len = 100) :: logname
      !! Strings for file names

   !=========solver parameter================
   real(kind = realdp) :: eps
      !! eps: tolerance for the outer iterations
   real(kind = realdp) :: def_rtol
      !! def_rtol: tolerance for the coarse-grid solver in two-level deflation method
   real(kind = realdp) :: cslp_mg_tol
      !! cslp_mg_tol: tolerance for Krylov based CSLP solver on coarse-grid levels
   integer :: m_iter
      !! m_iter: maximum number of outer iterations
   integer :: Algorithm
      !! Algorithm: solver identifier
   integer :: Irestart
      !! Irestart: Restart GMRES (1) or full GMRES (0) 
   integer :: M_flag
      !! M_flag: Flag for preconditioner type
   integer :: M2h_flag
      !! M2h_flag: In two-level deflation method, flag for using CSLP for the coarse level (1) or not (0)
   integer :: A2h_flag
      !! A2h_flag: In two-level deflation method, flag for the ReD method used for the coarse level operator
   integer :: def_mg_miter
      !! def_mg_miter: maximum number of iterations for deflation on the coarse levels
   integer :: cslp_mg_miter
      !! cslp_mg_miter: maximum number of iterations for solving CSLP by multigrid method on the coarsest level
   integer :: MG_flag
      !! MG_flag: flag for the type of multigrid method
   integer :: nx_min
      !! nx_min: minimum number of grid points in x-direction for the coarsest grid, it is used to specified the grid size of the coarsest grid when coarsening
   integer :: def_nlevel   
      !! def_nlevel: Number of levels in multilevel deflation method
   integer :: indf_level
      !! indf_level: The level that the linear system becomes indefinite
   integer :: Sv_L2
      !! Sv_L2: In multilevel method, a single-digit number that sets the tolerance for the second coarse grid iterations, ex. 1 means 1E-1
   integer :: Sv_L3   
      !! Sv_L3: In multilevel method, a single-digit number that sets the tolerance for the third coarse grid iterations, ex. 3 means 3E-1
   integer :: Sv_L4   
      !! Sv_L4: In multilevel method, a single-digit number that sets the tolerance for the fourth coarse grid iterations, ex. 10 means one iteration
   integer ::  npx0
      !! npx0: Number of partitions in x-direction
   integer ::  npy0
      !! npy0: Number of partitions in y-direction
   integer ::  npx
      !! npx: The coordinate/order of the current partition in x-direction
   integer ::  npy
      !! npy: The coordinate/order of the current partition in y-direction
   integer ::  nx_global
      !! nx_global: Global number of grid points in x-direction
   integer ::  ny_global
      !! ny_global: Global number of grid points in y-direction
   integer ::  nx
      !! nx: Local number of grid points of the current subdomain in x-direction
   integer ::  ny
      !! ny: Local number of grid points of the current subdomain in y-direction
   integer :: ID_XM1  
      !! ID_XM1: Identifier for neighboring partition in the -x direction  
   integer :: ID_XP1  
      !! ID_XP1: Identifier for neighboring partition in the +x direction  
   integer :: ID_YM1  
      !! ID_YM1: Identifier for neighboring partition in the -y direction  
   integer :: ID_YP1  
      !! ID_YP1: Identifier for neighboring partition in the +y direction  
   integer :: Iperiodic_X  
      !! Iperiodic_X: Periodicity flag in x-direction  
   integer :: Iperiodic_Y  
      !! Iperiodic_Y: Periodicity flag in y-direction
   real(kind = realdp) :: slx  
      !! slx: Physical length of the computational domain in x direction  
   real(kind = realdp) :: sly  
      !! sly: Physical length of the computational domain in y direction  
   real(kind = realdp) :: hx  
      !! hx: Grid spacing in x-direction  
   real(kind = realdp) :: hy  
      !! hy: Grid spacing in y-direction  
   real(kind = realdp) :: hxhy  
      !! hxhy: Product of grid spacings hx * hy, for this project hx must equal to hy, so hxhy=h^2
   integer, dimension(0:npMax - 1) :: i_offset  
      !! i_offset: an array that contains the index offset for MPI ranks in x direction   
   integer, dimension(0:npMax - 1) :: j_offset  
      !! j_offset: an array that contains the index offset for MPI ranks in y direction  
   integer, dimension(0:npMax - 1) :: i_nn  
      !! i_nn: an array that contains the number of grid points for MPI ranks in x direction  
   integer, dimension(0:npMax - 1) :: j_nn  
      !! j_nn: an array that contains the number of grid points for MPI ranks in y direction
   integer :: LAP
      !! The number of overlapping layers. The higher order deflation will be turned on once LAP > 1
   integer :: i_case
      !! Case identifier
   real(kind = realdp) :: freq  
      !! freq: the frequency for Helmholtz problem  
   real(kind = realdp) :: k0  
      !! k0: the wavenumber for constant-wavenumber model problem  
   real(kind = realdp) :: beta1  
      !! beta1: the real term for CSLP, usually 1.d0  
   real(kind = realdp) :: beta2  
      !! beta2: the complex shift for CSLP
   integer :: k_case
      !! k_case: Case identifier for constant (0) or non-constant (1) wavenumber
   integer :: flag_BCs
      !! flag_BCs: Boundary conditions flag
   integer :: ierr  
      !! ierr: MPI error code  
   integer :: np_size  
      !! np_size: the total number of MPI processes  
   integer :: my_id  
      !! my_id: MPI rank ID
   integer::STATUS(MPI_STATUS_SIZE), request1, request2
      !! MPI status array

   type, public  :: Gridpara
      !! Type definition for grid parameters
      integer(kind=4)     :: nx_global, ny_global
      integer(kind=4)     :: nx, ny
      real(kind = realdp) :: hx, hy, hxhy
      integer, dimension(0:npMax - 1) :: i_offset, j_offset, i_nn, j_nn
      real(kind = realdp), allocatable, dimension(:,:) :: wavenumber_k  
         !! wavenumber_k: the wavenumber profile of the domain  
      real(kind = realdp), allocatable, dimension(:,:) :: wavenumber_kh_pow_2  
         !! wavenumber_kh_pow_2 = (kh)^2
   end type Gridpara

end module comm_variable