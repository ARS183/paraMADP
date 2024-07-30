module user_module
   !! A module for IDR(s) to make data type and operators compatible
   use mpi
   use comm_variable
   use operators
   use CSLP_Solver
   use deflaion_setup
   !use preconditioned
   implicit none


   type matrix
      !! This is an empty matrix type
   end type matrix

   ! Overload * to define the matrix-vector multiplication using the matrix type

   INTERFACE OPERATOR(*)
      module procedure rmatvec, cmatvec
   END INTERFACE

   ! Define the preconditioner type 

   type preconditioner
      integer                          :: ni, nj
      real   (kind = realdp)               :: hx_c, hy_c
   end type preconditioner

   ! Overload / to define the preconditioner operation using the preconditioner type

   INTERFACE OPERATOR(/)
      module procedure rprecon, cprecon
   END INTERFACE


   contains

   function rmatvec( A, v ) 
      !! This is an empty matrix-vector multiplication for real-type input v
      type(matrix), intent(in)        :: A
      real(kind=realdp), intent(in)       :: v(:,:)
      real(kind=realdp)                   :: w(size(v,1),size(v,2))
      real(kind=realdp)                   :: rmatvec(size(v,1),size(v,2))
      
      w = v
      rmatvec = w

   end function rmatvec

   function cmatvec( A, v ) 
      !! This is complex matrix-vector multiplication
      type(matrix), intent(in)        :: A
      complex(kind=realdp), intent(in)    :: v(:,:)
      complex(kind=realdp)                :: w(size(v,1),size(v,2))
      complex(kind=realdp)                :: v_tmp(size(v,1),size(v,2))
      complex(kind=realdp)                :: cmatvec(size(v,1),size(v,2))
      
      v_tmp = v
      call Helmholtz2d_BC(v_tmp, w)
      cmatvec = w

   end function cmatvec

   function rprecon( v, M1 )
      !! This is an empty preconditioner for real-type input v
      implicit none

      type(preconditioner), intent(in)              :: M1
      real(kind=realdp), dimension(:,:), intent(in)     :: v
      real(kind=realdp), dimension(size(v,1),size(v,2)) :: rprecon
      
      rprecon = v

   end function rprecon

   function cprecon( v, M1 )
      !! This is an preconditioner for complex-type input v, i.e. M1^(-1)v
      !! We can manully call different preconditioners here. Here only use multigrid-based CLSP 
      type(preconditioner), intent(in)      :: M1
      complex(kind=realdp), intent(in)          :: v(:,:)
      complex(kind=realdp)                      :: w(size(v,1),size(v,2))
      complex(kind=realdp)                      :: cprecon(size(v,1),size(v,2))
      
      w = v
      cprecon = MGCSLP_invMx(w)

   end function cprecon

end module user_module
