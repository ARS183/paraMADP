! IDRS Induced Dimension Reduction method
!
!   x = IDRS( A, b, M, s,                             (mandatory)
!           [ tolerance, maximum_iterations, variant, (optional input scalars)
!             flag, relres, iterations,               (optional output scalars)
!             x0, U0, omega,                          (optional input vectors)
!             resvec, H] )                            (optional output vectors)
!
!   solves the system of linear equations A*x=b for x with the method IDRS
!   The dimension of the matrix A is NxN 
!   The dimension of the right-hand-side b is NxNRHS 
!   (NRHS: number of rhs vectors)
!
!   IDRS uses a user defined TYPE MATRIX, which is a structure that contains
!   all the parameters that define the system matrix A
!   The operator * must be overloaded to perform the operation A*v
!   To this end overloading must be defined for (real and complex) user 
!   routines that perform the matrix-vector multiplication. These
!   routines should be included in the module user_module that is used in the
!   module idrs_module.f90.
!
!   IDRS uses a user defined TYPE PRECONDITIONER, which is a structure that 
!   contains all the parameters that define the preconditioning matrix M
!   Overloading of the operator / is used to perform the preconditioning 
!   operation v/M1. To this end overloading must be defined for (real and 
!   complex) user routines that perform the preconditiong operations. These 
!   routines should be included in the module user_module that is used in the
!   module idrs_module.f90.
!
!   The precision of the complex and real variables are defined by the 
!   parameters rp (real precision)and cp (complex precision). These parameters 
!   should be set in the module user_module.
!
!   Input and output variables of IDRS:
!
!   Function output:
!        x: REAL (kind=rp) or COMPLEX (kind=cp) vector of dimension NxNRHS
!         
!   Function parameters:
!     Input, required:
!        A: TYPE MATRIX (user defined), A defines the system matrix 
!        b: REAL (kind=rp) or COMPLEX (kind=cp) matrix of dimension NxNRHS,
!           b is the matrix of right-hand-side vectors
!        M: TYPE PRECONDITIONER, M defines the preconditioner
!        s: INTEGER (must be > 0), s defines the dimension of the shadow space 
!           in IDRS, and by this the depth of the recursions 
!
!     Input, optional parameters:
!        tolerance: REAL (kind=rp). IDRS is terminated if |b-Ax|/|b| < tolerance
!           Default: tolerance = 1e-6
!        maximum_iterations: INTEGER. IDRS stops of the number of iterations 
!           exceeds maximum_iterations. 
!           Default: maximum_iterations = min(2*N,1000)
!        variant: CHARACTER*8. Selects the specific IDR-variant. 
!           Possible variants are:
!           variant = 'biortho ': IDR(s) with random shadow vectors,
!                                 'maintaining convergence' for computing omega
!           variant = 'minsync ': Piecewise constant shadow vectors, 
!                                 'minimal residual' for computing omega
!           variant = 'bicgstab': s must be set to 1. 
!                                 Initial residual is the shadow vector,
!                                 'minimal residual' for omega. 
!                                 This gives an algorithm that is
!                                 mathematically the same as BiCGSTAB
!           Default: variant = 'biortho '
!
!     Output, optional parameters:
!        flag: INTEGER. Indicates convergence condition
!           flag = 0: required tolerance satisfied
!           flag = 1: no convergence to the required tolerance within maximum
!                     number of iterations
!           flag = 2: final residual norm above tolerance
!           flag = 3: one of the iteration parameters became zero, 
!                     causing break down
!        relres: REAL (kind=rp), |b-Ax|/|b|. The norms are trace norms.
!        iterations: INTEGER. The number of iterations that has been performed 
!                             to compute the solution
!
!     Input, optional arrays:
!        x0: REAL (kind=rp) or COMPLEX (kind=cp) matrix of dimension NxNRHS,
!           initial guess for the solution. Default: x0 = 0.
!        U0: REAL (kind=rp) or COMPLEX (kind=cp) array of dimension NxNRHSxs,
!           Initial search space. Default: no initial search space.
!        omega: REAL (kind=rp) or COMPLEX (kind=cp) array of dimension n_omega,
!           where n_omega is the number of user defined omega-parameters.
!           The parameters are used cyclicly. User defined omega's make it 
!           possible to define and examine new IDR(s) variants.
!
!     Output, optional arrays:
!        resvec: REAL (kind=rp), array of dimension maximum_iterations+1
!           resvec contains for every iteration the norm of the residual
!        H: REAL (kind=rp) or COMPLEX (kind=cp) matrix of dimension 
!           (nritz+1) x nritz, in which nritz is the number of ritz values
!           that can be computed. Note that in s+1 iterations s additional 
!           ritz values can be computed.  H is an extended upper Hessenberg 
!           matrix of which the upper bandwidth is s. It can be used for 
!           analysis purposes, the eigenvalues of the leading submatrices of 
!           H are ritzvalues
!
!   Acknowledgement: Duncan van der Heul, Reinaldo Astudillo, Jan de Gier and 
!                    Marielba Rojas are gratefully acknowledged for their 
!                    advise on many aspects of this software. 
!
!   The software is distributed without any warranty.
!
!   Martin van Gijzen
!   Copyright (c) July 2015
!

module idrs_module
   !! IDRS Induced Dimension Reduction method
   use mpi
   use comm_variable
   use read_setup
   use mpi_setup
   use user_module

   implicit none

   private
   public :: IDRS, TRACE_DOT, FROB_NORM, P_DOT
   !public :: CIDRS, RIDRS, TRACE_DOT, FROB_NORM, P_DOT

   interface TRACE_DOT
      module procedure CTRACE_DOT, RTRACE_DOT, RCTRACE_DOT
   end interface

   interface P_DOT
      module procedure CP_DOT, RP_DOT
   end interface

   interface FROB_NORM
      module procedure CFROB_NORM, RFROB_NORM
   end interface

   interface IDRS
      module procedure CIDRS, RIDRS
   end interface

contains

   ! Trace inner product of complex matrices
   function CTRACE_DOT(v,w)
      complex(kind=realdp), intent(in)      :: v(1-LAP:nx+LAP,1-LAP:ny+LAP), w(1-LAP:nx+LAP,1-LAP:ny+LAP)
      complex(kind=realdp)                  :: CTRACE_DOT_local
      complex(kind=realdp)                  :: CTRACE_DOT
      
      CTRACE_DOT_local = sum( conjg(v(1:nx,1:ny))*w(1:nx,1:ny) )
      call MPI_Allreduce(CTRACE_DOT_local,   &
                        CTRACE_DOT,          &
                        1,                   &
                        MPI_DOUBLE_COMPLEX,  &
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

   end function CTRACE_DOT

   ! Trace inner product of real matrices
   function RTRACE_DOT(v,w)
      real(kind=realdp), intent(in)      :: v(1-LAP:nx+LAP,1-LAP:ny+LAP), w(1-LAP:nx+LAP,1-LAP:ny+LAP)
      real(kind=realdp)                  :: RTRACE_DOT_local
      real(kind=realdp)                  :: RTRACE_DOT

      RTRACE_DOT_local = sum( v(1:nx,1:ny)*w(1:nx,1:ny) )
      
      call MPI_Allreduce(RTRACE_DOT_local,   &
                        RTRACE_DOT,          &
                        1,                   &
                        MPI_DOUBLE_PRECISION,&
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

   end function RTRACE_DOT

   ! Trace inner product of real and complex matrices
   function RCTRACE_DOT(v,w)
      real(kind=realdp), intent(in)      :: v(1-LAP:nx+LAP,1-LAP:ny+LAP)
      complex(kind=realdp), intent(in)   :: w(1-LAP:nx+LAP,1-LAP:ny+LAP)
      complex(kind=realdp)               :: RCTRACE_DOT
      complex(kind=realdp)               :: RCTRACE_DOT_local

      RCTRACE_DOT_local = sum( v(1:nx,1:ny)*w(1:nx,1:ny) )
      call MPI_Allreduce(RCTRACE_DOT_local,  &
                        RCTRACE_DOT,         &
                        1,                   &
                        MPI_DOUBLE_COMPLEX,  &
                        MPI_SUM,             &
                        MPI_COMM_WORLD,      &
                        ierr)

   end function RCTRACE_DOT

   ! P inner product of complex matrices
   function CP_DOT(P,R0,w,s)
      integer                                   :: s
      real(kind=realdp),    allocatable, intent(in) :: P(:,:,:)!(1-LAP:nx+LAP,1-LAP:ny+LAP,s) 
      complex(kind=realdp), allocatable, intent(in) :: R0(:,:)
      complex(kind=realdp), intent(in)              :: w(1-LAP:nx+LAP,1-LAP:ny+LAP)
      complex(kind=realdp)                          :: v(s), CP_DOT(s)
      complex(kind=realdp)                          :: vi_local(s)
      integer                                   :: i, j, N, low(s), up(s), step, nrhs

      if ( allocated(P) ) then
         ! Biortho: P has orthogonal random numbers
         !write(*,*) "Hi, Do you call me? CP_DOT"
         do i = 1,s
            vi_local(i) = sum( P(1:nx,1:ny,i)*w(1:nx,1:ny) )
         end do
         call MPI_Allreduce(vi_local,           &
                           v,                   &
                           s,                   &
                           MPI_DOUBLE_COMPLEX,  &
                           MPI_SUM,             &
                           MPI_COMM_WORLD,      &
                           ierr)

      else if ( allocated(R0) ) then
         ! BiCGSTAB: shadow vector equal to initial residual
         vi_local(1) = sum( conjg(R0(1:nx,1:ny))*w(1:nx,1:ny) ) !no!!! sum( conjg(R0)*w ) no!!!
         !made a big mistake here, only found after done numerous tests, only not converge for MP-4 nx=2209,np>13?!!
         call MPI_Allreduce(vi_local(1),        &
                           v(1),                &
                           1,                   &
                           MPI_DOUBLE_COMPLEX,  &
                           MPI_SUM,             &
                           MPI_COMM_WORLD,      &
                           ierr)

      else
         ! Minsync: P is piecewise constant 
         N    = size(w,1)
         nrhs = size(w,2)
         step = N/s
         low(1) = 1
         do i = 1,s-1
            low(i+1) = i*step+1
            up(i) = i*step
         end do
         up(s) = N

         do i = 1,s
            v(i)  = 0.
            do j = 1, nrhs
               v(i) = v(i) + sum( w(low(i):up(i),j) )
            end do
         end do
      end if
      
      CP_DOT = v

   end function CP_DOT


   ! P inner product of real matrices
   function RP_DOT(P,R0,w,s)
      integer                              :: s
      real(kind=realdp),allocatable,intent(in) :: P(:,:,:)!(1-LAP:nx+LAP,1-LAP:ny+LAP,s)
      real(kind=realdp),allocatable,intent(in) :: R0(:,:) 
      real(kind=realdp), intent(in)            :: w(1-LAP:nx+LAP,1-LAP:ny+LAP)
      real(kind=realdp)                        :: v(s), RP_DOT(s)
      real(kind=realdp)                        :: vi_local(s)
      integer                              :: i, j, N, nrhs, low(s), up(s), step

      if ( allocated(P) ) then
         ! Biortho: P has orthogonal random numbers
         !write(*,*) "Hi, Do you call me? RP_DOT"
         do i = 1,s
            vi_local(i) = sum( P(1:nx,1:ny,i)*w(1:nx,1:ny) )
         end do
         call MPI_Allreduce(vi_local,           &
                           v,                   &
                           s,                   &
                           MPI_DOUBLE_PRECISION,&
                           MPI_SUM,             &
                           MPI_COMM_WORLD,      &
                           ierr)
      else if (allocated(R0) ) then
         ! Bi-CGSTAB: shadow vector equal to initial residual
         vi_local(1) = sum( R0(1:nx,1:ny)*w(1:nx,1:ny) ) 
         !!made a big mistake here, only found after done numerous tests,
         call MPI_Allreduce(vi_local(1),        &
                           v(1),                &
                           1,                   &
                           MPI_DOUBLE_PRECISION,&
                           MPI_SUM,             &
                           MPI_COMM_WORLD,      &
                           ierr)
      else
         ! Minsync: P has piecewise constant columns
         N = size(w,1)
         nrhs = size(w,2)
         step = N/s
         low(1) = 1
         up(s) = N
         do i = 1,s-1
            low(i+1) = i*step+1
            up(i) = i*step
         end do

         do i = 1,s
            v(i) = 0.
            do j = 1, nrhs
               v(i) = sum( w(low(i):up(i),j) )
            end do
         end do
      end if
         
      RP_DOT = v

   end function RP_DOT 

   
   function CFROB_NORM(v)
   !! Frobenius norm of complex matrix
      complex(kind=realdp), intent(in)      :: v(1-LAP:nx+LAP,1-LAP:ny+LAP)
      real(kind=realdp)                     :: CFROB_NORM
      real(kind=realdp)                     :: CFROB_NORM_local

      CFROB_NORM_local = sum( conjg(v(1:nx,1:ny))*v(1:nx,1:ny) ) 
      call MPI_Allreduce(CFROB_NORM_local,       &
                         CFROB_NORM,             &
                         1,                      &
                         MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,                &
                         MPI_COMM_WORLD,         &
                         ierr)
      CFROB_NORM = sqrt(CFROB_NORM)

   end function CFROB_NORM


   
   function RFROB_NORM(v)
   !! Frobenius norm of real matrix
      real(kind=realdp), intent(in)      :: v(1-LAP:nx+LAP,1-LAP:ny+LAP)
      real(kind=realdp)                  :: RFROB_NORM_local
      real(kind=realdp)                  :: RFROB_NORM

      RFROB_NORM_local =  sum( v(1:nx,1:ny)*v(1:nx,1:ny) ) 
      call MPI_Allreduce(RFROB_NORM_local,       &
                         RFROB_NORM,             &
                         1,                      &
                         MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,                &
                         MPI_COMM_WORLD,         &
                         ierr)
      RFROB_NORM = sqrt(RFROB_NORM)

   end function RFROB_NORM


   function CIDRS( A, b, M1, s, &                            ! required
                     tolerance, maximum_iterations, variant, & ! optional input
                     flag, relres, iterations, &               ! optional output
                     x0, U0, omega, resvec, H )                ! optional arrays

      IMPLICIT NONE

      ! Required input parameters:
      type(matrix), intent(in)               :: A               ! system matrix
      complex(kind=realdp),allocatable,intent(in):: b(:,:)!(1-LAP:nx+LAP,1-LAP:ny+LAP)          ! system rhs
      !complex(kind=realdp), allocatable,dimension(:,:) :: b
      type(preconditioner), intent(in)       :: M1              ! preconditioner
      integer, intent(in)                    :: s               ! s parameter

      ! Solution:
      complex(kind=realdp)                       :: CIDRS(size(b,1),size(b,2))

      ! Optional input parameters:
      real(kind=realdp), optional, intent(in)    :: tolerance 
      integer, optional, intent(in)          :: maximum_iterations 
      character(len=8), optional, intent(in) :: variant     

      ! Optional output parameters:
      integer, optional, intent(out)         :: flag        
      real(kind=realdp), optional, intent(out)   :: relres     
      integer, optional, intent(out)         :: iterations

      ! Optional input arrays:
      complex(kind=realdp), optional, intent(in) :: x0(:,:)  
      complex(kind=realdp), optional, intent(in) :: U0(:,:,:) 
      complex(kind=realdp), optional, intent(in) :: omega(:) 

      ! Optional output arrays:
      real(kind=realdp), optional, intent(out)   :: resvec(:)   !Jchen, change it to relresvec
      complex(kind=realdp), optional, intent(out):: H(:,:)   
      
      ! Local arrays:
      real(kind=realdp),    allocatable          :: P(:,:,:) 
      complex(kind=realdp), allocatable          :: R0(:,:) 

      complex(kind=realdp), allocatable, dimension(:)     :: f,mu,alpha,beta,gamma
      complex(kind=realdp), allocatable, dimension(:,:)   :: x,r,v,t,M
      complex(kind=realdp), allocatable, dimension(:,:,:) :: G,U
      complex(kind=realdp)                       :: om, tr    
      real(kind=realdp)                          :: nr, nt, rho, kappa

      ! complex(kind=realdp)                       :: x(size(b,1),size(b,2))
      ! complex(kind=realdp)                       :: G(size(b,1),size(b,2),s)
      ! complex(kind=realdp)                       :: U(size(b,1),size(b,2),s)
      ! complex(kind=realdp)                       :: r(size(b,1),size(b,2)) 
      ! complex(kind=realdp)                       :: v(size(b,1),size(b,2))   
      ! complex(kind=realdp)                       :: t(size(b,1),size(b,2))  
      ! complex(kind=realdp)                       :: M(s,s), f(s), mu(s)
      ! complex(kind=realdp)                       :: alpha(s), beta(s), gamma(s)
      
      !include 'idrs_body.f90'
      !------------------------------idrs_body.f90-------------------------------------
      ! Declarations:
      integer               :: n                 ! dimension of the system
      integer               :: nrhs              ! Number of RHS-vectors
      integer               :: maxit             ! maximum number of iterations
      integer               :: method            ! which IDR(s) variant?
      real(kind=realdp)         :: tol               ! actual tolerance
      integer               :: info              ! convergence indicator
      logical               :: out_flag          ! store flag
      logical               :: out_relres        ! store relres
      logical               :: out_iterations    ! store number of iterations
      logical               :: inispace          ! initial search space
      logical               :: user_omega        ! user defined omega present
      integer               :: n_omega           ! number of user defined omega's
      logical               :: out_resvec        ! store residual norms
      logical               :: out_H             ! store iteration parameters in H
      integer               :: nritz             ! Number of wanted ritz values

      integer               :: iter              ! number of iterations
      integer               :: ii                ! inner iterations index
      integer               :: jj                ! G-space index
      real(kind=realdp)         :: normb, normr, tolb! for tolerance check
      integer               :: i,j,k,l           ! loop counters

      ! Problem size:
      n    = size(b,1)
      ! Number of right-hand-side vectors:
      nrhs = size(b,2)

      ! allocate local arrays:
      allocate(f(s), mu(s), alpha(s), beta(s), gamma(s))
      allocate(x(1-LAP:nx+LAP,1-LAP:ny+LAP) , r(1-LAP:nx+LAP,1-LAP:ny+LAP) , v(1-LAP:nx+LAP,1-LAP:ny+LAP) )
      allocate(t(1-LAP:nx+LAP,1-LAP:ny+LAP) , M(s,s))
      allocate(G(1-LAP:nx+LAP,1-LAP:ny+LAP,s), U(1-LAP:nx+LAP,1-LAP:ny+LAP,s))

      
      ! Check optional input parameters:
      if ( present(tolerance) ) then
         if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
         tol = tolerance 
      else
         tol = 1e-6
      endif

      maxit=min(2*n,1000)
      if ( present(maximum_iterations) ) maxit = maximum_iterations 
      
      method = 1 ! biortho   
      if ( present(variant) ) then
         if ( variant == 'minsync' ) then
            method = 2
         elseif ( variant == 'bicgstab' ) then 
            method = 3
         endif
      endif

      ! Initialize the output variables 
      out_flag       = present(flag)
      if ( out_flag )       flag = -1 
      out_relres     = present(relres)
      if ( out_relres)      relres = 1. 
      out_iterations = present(iterations)
      if ( out_iterations ) iterations = 0 

      ! Check optional input arrays:
      x = 0.
      if ( present(x0) ) x = x0
         
      U = 0.
      inispace =  present(U0)
      if ( inispace ) U = U0

      user_omega = present(omega)
      if ( user_omega ) then
         n_omega = size(omega)
      end if

      ! Check output arrays
      out_resvec     = present(resvec)
      if ( out_resvec ) then
         if ( maxit+1 > size(resvec) ) &
         stop "Length of vector with residual norms too small, should be maxit+1"
      end if

      out_H = present(H)
      if ( out_H ) then
         nritz = size(H,1)-1
         if ( size(H,2) /= nritz ) &
            stop "Second dimension of H incompatible, with first"
         H = 0.
      end if

      ! compute initial residual, set absolute tolerance
      
      normb = FROB_NORM(b)
      tolb = tol * normb
      r = b - A*x
      normr = FROB_NORM(r)
      if ( out_resvec ) resvec(1)= normr/normb
      
      ! check if the initial solution is not already a solution within the prescribed
      ! tolerance
      if (normr <= tolb) then      
         if ( out_iterations ) iterations = 0               
         if ( out_flag )       flag  = 0
         if ( out_relres )     relres = normr/normb
         return
      end if

      ! Define P and kappa (depending on the method)
      if ( method == 1 ) then
         !allocate(P(n,nrhs,s))
         allocate( P(1-LAP:nx+LAP,1-LAP:ny+LAP,s) )
         call RANDOM_SEED
         call RANDOM_NUMBER(P)
         do j = 1,s
            do k = 1,j-1
               alpha(k) = TRACE_DOT( P(:,:,k),P(:,:,j) )
               P(:,:,j) = P(:,:,j) - alpha(k)*P(:,:,k)
            end do
            !P(:,:,j) = P(:,:,j)/FROB_NORM(P(:,:,j))
            P(:,:,j) = P(:,:,j)/RFROB_NORM(P(:,:,j))
         end do
         kappa = 0.7
      elseif ( method == 2 ) then
      ! P is piecewise constant, minimum residual for omega
         kappa = 0.
      elseif ( method == 3 ) then
         if ( s /= 1 ) stop "s=1 is required for variant bicgstab"
         !allocate(R0(n,nrhs))
         allocate( R0(1-LAP:nx+LAP,1-LAP:ny+LAP) )
         R0 = r
         kappa = 0.
      endif

      ! Initialize local variables:
      M = 0.
      om = 1.
      iter = 0
      info = -1
      jj = 0
      ii = 0

      
      ! This concludes the initialisation phase

      ! Main iteration loop, build G-spaces:
      
      do while (  info < 0 )  ! start of iteration loop
      
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Generate s vectors in G_j
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

         ! New right-hand side for small system:
         f = P_DOT( P, R0, r, s )

         do k=1, s

            ! Update inner iteration counter
            ii = ii + 1

            ! Compute new v
            v = r 
            if ( jj > 0 ) then

               ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
               do i = k,s
                  gamma(i) = f(i)
                  do j = k, i-1
                     gamma(i) = gamma(i) - M(i,j)*gamma(j)
                  end do
                  gamma(i) = gamma(i)/M(i,i)
                  v = v - gamma(i)*G(:,:,i)
               end do

               ! Compute new U(:,:,k)
               t = om*(v/M1)
               do i = k,s
                  t = t + gamma(i)*U(:,:,i)
               end do
               U(:,:,k) = t

               ! Compute Hessenberg matrix?
               if ( out_H .and. ii <= nritz ) &
                  H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)

            else if ( .not. inispace ) then

               ! Updates for the first s iterations (in G_0):
               U(:,:,k) = v/M1

            end if

            ! Compute new G(:,:,k), G(:,:,k) is in space G_j
            G(:,:,k) = A*U(:,:,k)
         
            ! Bi-Orthogonalise the new basis vectors: 
            mu = P_DOT( P, R0, G(:,:,k), s )
            do i = 1,k-1
               alpha(i) = mu(i)
               do j = 1, i-1
                  alpha(i) = alpha(i) - M(i,j)*alpha(j)
               end do
               alpha(i) = alpha(i)/M(i,i)
               G(:,:,k) = G(:,:,k) - G(:,:,i)*alpha(i)
               U(:,:,k) = U(:,:,k) - U(:,:,i)*alpha(i)
               mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
            end do
            M(k:s,k) = mu(k:s)

            ! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz .and. k  > 1 ) &
               H(ii-k+1:ii-1,ii) =  alpha(1:k-1)/beta(1:k-1)

            ! Break down?
            if ( abs(M(k,k)) <= tiny(tol) ) then
               info = 3
               exit
            end if

            ! Make r orthogonal to p_i, i = 1..k, update solution and residual 
            beta(k) = f(k)/M(k,k)
            r = r - beta(k)*G(:,:,k)
            x = x + beta(k)*U(:,:,k)

            ! New f = P'*r (first k  components are zero)
            if ( k < s ) then
               f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
            end if

            ! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz ) then     
               H(ii,ii) = 1./beta(k)
               l = max(1,ii-s)
               H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
               H(l:ii+1,ii)   = H(l:ii+1,ii)/om
            end if

            ! Check for convergence
            normr = FROB_NORM(r)
            iter = iter + 1

            if(my_id == 0) then
               write(*,"(I9,E17.9)")  iter, normr
            endif  

            if ( out_resvec ) resvec(iter + 1) = normr/normb
            if ( normr < tolb ) then
               info = 0
               exit
            elseif ( iter == maxit ) then
               info = 1
               exit
            end if 

         end do ! Now we have computed s+1 vectors in G_j
         if ( info >= 0 )  then
            exit
         end if

         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Compute first residual in G_j+1
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

         ! Update G-space counter
         jj = jj + 1

         ! Compute first residual in G_j+1
         ! Note: r is already perpendicular to P so v = r
         ! Preconditioning:
         v = r/M1
         t = A*v

         ! Computation of a new omega
         if ( user_omega ) then
            i = mod(jj,n_omega)
            if ( i == 0 ) i = n_omega
            om = omega(i)
         elseif ( kappa == 0. ) then

            ! Minimal residual (same as in Bi-CGSTAB):
            om = TRACE_DOT(t,r)/TRACE_DOT(t,t)
         else

            ! 'Maintaining the convergence':
            nr = FROB_NORM(r)
            nt = FROB_NORM(t)
            tr = TRACE_DOT(t,r)
            rho = abs(tr/(nt*nr))
            om=tr/(nt*nt)
            if ( rho < kappa ) then
               om = om*kappa/rho
            end if
         end if
         if ( abs(om) <= epsilon(tol) ) then 
            info = 3
            exit
         end if 

         ! Update solution and residual
         r = r - om*t 
         x = x + om*v 

         ! Check for convergence
         normr =FROB_NORM(r)
         iter = iter + 1

         if(my_id == 0) then
            write(*,"(I9,E17.9)")  iter, normr
         endif  

         if ( out_resvec ) resvec(iter + 1) = normr/normb
         if ( normr < tolb ) then
            info = 0
         elseif ( iter == maxit ) then
            info = 1
         end if 
        
      end do ! end of while loop

      ! Set output parameters
      r = b - A*x
      normr = FROB_NORM(r)

      if ( info == 0 .and. normr > tolb ) info = 2
      if ( out_iterations ) iterations = iter
      if ( out_relres )     relres=normr/normb
      if ( out_flag )       flag = info
      
      !------------------------------idrs_body.f90-------------------------------------
      cidrs = x

      ! deallocate local arrays:
      deallocate(f, mu, alpha, beta, gamma)
      deallocate(x, r, v)
      deallocate(t, M)
      deallocate(G, U)
      
   end function CIDRS

   function RIDRS( A, b, M1, s, &                            ! required
                  tolerance, maximum_iterations, variant, & ! optional input
                  flag, relres, iterations, &               ! optional output
                  x0, U0, omega, H, resvec )                ! optional arrays

      IMPLICIT NONE

      ! Required input and output parameters:
      type(matrix), intent(in)               :: A               ! system matrix
      real(kind=realdp),allocatable,intent(in):: b(:,:)
      !real(kind=realdp), intent(in)              :: b(1-LAP:nx+LAP,1-LAP:ny+LAP)         ! system rhs
      !real(kind=realdp), allocatable,dimension(:,:) :: b
      type(preconditioner), intent(in)       :: M1              ! preconditioner
      integer, intent(in)                    :: s               ! s parameter

      ! Solution
      real(kind=realdp)                          :: RIDRS(size(b,1),size(b,2))

      ! Optional input parameters:
      real(kind=realdp), optional, intent(in)    :: tolerance
      integer, optional, intent(in)          :: maximum_iterations
      character(len=8), optional, intent(in) :: variant

      ! Optional output parameters:
      integer, optional, intent(out)         :: flag
      real(kind=realdp), optional, intent(out)   :: relres
      integer, optional, intent(out)         :: iterations 

      ! Optional input arrays:
      real(kind=realdp), optional, intent(in)    :: x0(:,:)
      real(kind=realdp), optional, intent(in)    :: U0(:,:,:)
      real(kind=realdp), optional, intent(in)    :: omega(:)

      ! Optional output arrays
      real(kind=realdp), optional, intent(out)   :: resvec(:)
      real(kind=realdp), optional, intent(out)   :: H(:,:)  
      
      ! Local arrays:
      real(kind=realdp), allocatable             :: P(:,:,:) 
      real(kind=realdp), allocatable             :: R0(:,:) 
      real(kind=realdp)                          :: om, kappa           
      real(kind=realdp)                          :: nr, nt, tr, rho 

      real(kind=realdp), allocatable, dimension(:)     :: f,mu,alpha,beta,gamma
      real(kind=realdp), allocatable, dimension(:,:)   :: x,r,v,t,M
      real(kind=realdp), allocatable, dimension(:,:,:) :: G,U

      ! real(kind=realdp)                          :: x(size(b,1),size(b,2))
      ! real(kind=realdp)                          :: G(size(b,1),size(b,2),s)
      ! real(kind=realdp)                          :: U(size(b,1),size(b,2),s)
      ! real(kind=realdp)                          :: r(size(b,1),size(b,2))            
      ! real(kind=realdp)                          :: v(size(b,1),size(b,2))         
      ! real(kind=realdp)                          :: t(size(b,1),size(b,2))        
      ! real(kind=realdp)                          :: M(s,s), f(s), mu(s)
      ! real(kind=realdp)                          :: alpha(s), beta(s), gamma(s) 
      
      !include 'idrs_body.f90'
      !----------------------------------idrs_body.f90---------------------------------------!
      ! Declarations:
      integer               :: n                 ! dimension of the system
      integer               :: nrhs              ! Number of RHS-vectors
      integer               :: maxit             ! maximum number of iterations
      integer               :: method            ! which IDR(s) variant?
      real(kind=realdp)         :: tol               ! actual tolerance
      integer               :: info              ! convergence indicator
      logical               :: out_flag          ! store flag
      logical               :: out_relres        ! store relres
      logical               :: out_iterations    ! store number of iterations
      logical               :: inispace          ! initial search space
      logical               :: user_omega        ! user defined omega present
      integer               :: n_omega           ! number of user defined omega's
      logical               :: out_resvec        ! store residual norms
      logical               :: out_H             ! store iteration parameters in H
      integer               :: nritz             ! Number of wanted ritz values

      integer               :: iter              ! number of iterations
      integer               :: ii                ! inner iterations index
      integer               :: jj                ! G-space index
      real(kind=realdp)         :: normb, normr, tolb! for tolerance check
      integer               :: i,j,k,l           ! loop counters

      ! Problem size:
      n    = size(b,1)
      ! Number of right-hand-side vectors:
      nrhs = size(b,2)

      ! allocate local arrays:
      allocate(f(s), mu(s), alpha(s), beta(s), gamma(s))
      allocate(x(1-LAP:nx+LAP,1-LAP:ny+LAP) , r(1-LAP:nx+LAP,1-LAP:ny+LAP) , v(1-LAP:nx+LAP,1-LAP:ny+LAP) )
      allocate(t(1-LAP:nx+LAP,1-LAP:ny+LAP) , M(s,s))
      allocate(G(1-LAP:nx+LAP,1-LAP:ny+LAP,s), U(1-LAP:nx+LAP,1-LAP:ny+LAP,s))
      
      ! Check optional input parameters:
      if ( present(tolerance) ) then
         if ( tolerance < 0 ) stop "Illegal value for parameter tolerance"
         tol = tolerance 
      else
         tol = 1e-6
      endif

      maxit=min(2*n,1000)
      if ( present(maximum_iterations) ) maxit = maximum_iterations 
      
      method = 1 ! biortho   
      if ( present(variant) ) then
         if ( variant == 'minsync' ) then
            method = 2
         elseif ( variant == 'bicgstab' ) then 
            method = 3
         endif
      endif

      ! Initialize the output variables 
      out_flag       = present(flag)
      if ( out_flag )       flag = -1 
      out_relres     = present(relres)
      if ( out_relres)      relres = 1. 
      out_iterations = present(iterations)
      if ( out_iterations ) iterations = 0 

      ! Check optional input arrays:
      x = 0.
      if ( present(x0) ) x = x0
         
      U = 0.
      inispace =  present(U0)
      if ( inispace ) U = U0

      user_omega = present(omega)
      if ( user_omega ) then
         n_omega = size(omega)
      end if

      ! Check output arrays
      out_resvec     = present(resvec)
      if ( out_resvec ) then
         if ( maxit+1 > size(resvec) ) &
         stop "Length of vector with residual norms too small, should be maxit+1"
      end if

      out_H = present(H)
      if ( out_H ) then
         nritz = size(H,1)-1
         if ( size(H,2) /= nritz ) &
            stop "Second dimension of H incompatible, with first"
         H = 0.
      end if

      ! compute initial residual, set absolute tolerance
      
      normb = FROB_NORM(b)
      tolb = tol * normb
      r = b - A*x
      normr = FROB_NORM(r)
      if ( out_resvec ) resvec(1)= normr/normb
      
      ! check if the initial solution is not already a solution within the prescribed
      ! tolerance
      if (normr <= tolb) then      
         if ( out_iterations ) iterations = 0               
         if ( out_flag )       flag  = 0
         if ( out_relres )     relres = normr/normb
         return
      end if

      ! Define P and kappa (depending on the method)
      if ( method == 1 ) then
         allocate( P(n,nrhs,s) )
         call RANDOM_SEED
         call RANDOM_NUMBER(P)
         do j = 1,s
            do k = 1,j-1
               alpha(k) = TRACE_DOT( P(:,:,k),P(:,:,j) )
               P(:,:,j) = P(:,:,j) - alpha(k)*P(:,:,k)
            end do
            ! P(:,:,j) = P(:,:,j)/FROB_NORM(P(:,:,j))
            P(:,:,j) = P(:,:,j)/RFROB_NORM(P(:,:,j))
         end do
         kappa = 0.7
      elseif ( method == 2 ) then
      ! P is piecewise constant, minimum residual for omega
         kappa = 0.
      elseif ( method == 3 ) then
         if ( s /= 1 ) stop "s=1 is required for variant bicgstab"
         allocate( R0(n,nrhs) )
         R0 = r
         kappa = 0.
      endif

      ! Initialize local variables:
      M = 0.
      om = 1.
      iter = 0
      info = -1
      jj = 0
      ii = 0

      
      ! This concludes the initialisation phase

      ! Main iteration loop, build G-spaces:
      
      do while (  info < 0 )  ! start of iteration loop
      
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Generate s vectors in G_j
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

         ! New right-hand side for small system:
         f = P_DOT( P, R0, r, s )

         do k=1, s

            ! Update inner iteration counter
            ii = ii + 1

            ! Compute new v
            v = r 
            if ( jj > 0 ) then

               ! Solve small system (Note: M is lower triangular) and make v orthogonal to P:
               do i = k,s
                  gamma(i) = f(i)
                  do j = k, i-1
                     gamma(i) = gamma(i) - M(i,j)*gamma(j)
                  end do
                  gamma(i) = gamma(i)/M(i,i)
                  v = v - gamma(i)*G(:,:,i)
               end do

               ! Compute new U(:,:,k)
               t = om*(v/M1)
               do i = k,s
                  t = t + gamma(i)*U(:,:,i)
               end do
               U(:,:,k) = t

               ! Compute Hessenberg matrix?
               if ( out_H .and. ii <= nritz ) &
                  H(ii-s:ii-k,ii)   = -gamma(k:s)/beta(k:s)

            else if ( .not. inispace ) then

               ! Updates for the first s iterations (in G_0):
               U(:,:,k) = v/M1

            end if

            ! Compute new G(:,:,k), G(:,:,k) is in space G_j
            G(:,:,k) = A*U(:,:,k)
         
            ! Bi-Orthogonalise the new basis vectors: 
            mu = P_DOT( P, R0, G(:,:,k), s )
            do i = 1,k-1
               alpha(i) = mu(i)
               do j = 1, i-1
                  alpha(i) = alpha(i) - M(i,j)*alpha(j)
               end do
               alpha(i) = alpha(i)/M(i,i)
               G(:,:,k) = G(:,:,k) - G(:,:,i)*alpha(i)
               U(:,:,k) = U(:,:,k) - U(:,:,i)*alpha(i)
               mu(k:s)  = mu(k:s)  - M(k:s,i)*alpha(i)
            end do
            M(k:s,k) = mu(k:s)

            ! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz .and. k  > 1 ) &
               H(ii-k+1:ii-1,ii) =  alpha(1:k-1)/beta(1:k-1)

            ! Break down?
            if ( abs(M(k,k)) <= tiny(tol) ) then
               info = 3
               exit
            end if

            ! Make r orthogonal to p_i, i = 1..k, update solution and residual 
            beta(k) = f(k)/M(k,k)
            r = r - beta(k)*G(:,:,k)
            x = x + beta(k)*U(:,:,k)

            ! New f = P'*r (first k  components are zero)
            if ( k < s ) then
               f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k)
            end if

            ! Compute Hessenberg matrix?
            if ( out_H .and. ii <= nritz ) then     
               H(ii,ii) = 1./beta(k)
               l = max(1,ii-s)
               H(l+1:ii+1,ii) = (H(l+1:ii+1,ii) - H(l:ii,ii))
               H(l:ii+1,ii)   = H(l:ii+1,ii)/om
            end if

            ! Check for convergence
            normr = FROB_NORM(r)
            iter = iter + 1

            if(my_id == 0) then
               write(*,"(I9,E17.9)")  iter, normr
            endif  

            if ( out_resvec ) resvec(iter + 1) = normr/normb
            if ( normr < tolb ) then
               info = 0
               exit
            elseif ( iter == maxit ) then
               info = 1
               exit
            end if 

         end do ! Now we have computed s+1 vectors in G_j
         if ( info >= 0 )  then
            exit
         end if

         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Compute first residual in G_j+1
         !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

         ! Update G-space counter
         jj = jj + 1

         ! Compute first residual in G_j+1
         ! Note: r is already perpendicular to P so v = r
         ! Preconditioning:
         v = r/M1
         t = A*v

         ! Computation of a new omega
         if ( user_omega ) then
            i = mod(jj,n_omega)
            if ( i == 0 ) i = n_omega
            om = omega(i)
         elseif ( kappa == 0. ) then

            ! Minimal residual (same as in Bi-CGSTAB):
            om = TRACE_DOT(t,r)/TRACE_DOT(t,t)
         else

            ! 'Maintaining the convergence':
            nr = FROB_NORM(r)
            nt = FROB_NORM(t)
            tr = TRACE_DOT(t,r)
            rho = abs(tr/(nt*nr))
            om=tr/(nt*nt)
            if ( rho < kappa ) then
               om = om*kappa/rho
            end if
         end if
         if ( abs(om) <= epsilon(tol) ) then 
            info = 3
            exit
         end if 

         ! Update solution and residual
         r = r - om*t 
         x = x + om*v 

         ! Check for convergence
         normr =FROB_NORM(r)
         iter = iter + 1

         if(my_id == 0) then
            write(*,"(I9,E17.9)")  iter, normr
         endif  

         if ( out_resvec ) resvec(iter + 1) = normr/normb
         if ( normr < tolb ) then
            info = 0
         elseif ( iter == maxit ) then
            info = 1
         end if 
         
      end do ! end of while loop

      ! Set output parameters
      r = b - A*x
      normr = FROB_NORM(r)


      if ( info == 0 .and. normr > tolb ) info = 2
      if ( out_iterations ) iterations = iter
      if ( out_relres )     relres=normr/normb
      if ( out_flag )       flag = info
      
      !----------------------------------idrs_body.f90---------------------------------------!

      RIDRS = x
      
      ! deallocate local arrays:
      deallocate(f, mu, alpha, beta, gamma)
      deallocate(x, r, v)
      deallocate(t, M)
      deallocate(G, U)
      
      end function RIDRS

   end module idrs_module
