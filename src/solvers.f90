module solvers
  !! Serveral GMRES-type Krylov solvers. Only applies to the (default) finest grid level
  use mpi
  use comm_variable
  use read_setup
  use operators
  use CSLP_Solver
  use deflaion_setup
  implicit none

contains

  !====================================================================================
  subroutine fullgmres(b,u,Rerror,iter)
    !! Full GMRES without precondition
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                        :: k,ki,j
    real(kind = realdp)                   :: b_norm, res_norm
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res, u0, Au0, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP), u0(1-LAP:nx+LAP,1-LAP:ny+LAP), Au0(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)  ! what if b is real?

    res  = (0.d0,0.d0)
    Au0  = (0.d0,0.d0)
    sn   = (0.d0,0.d0)
    cs   = (0.d0,0.d0)
    V    = (0.d0,0.d0)
    H    = (0.d0,0.d0)
    beta = (0.d0,0.d0)

    u0=u

    call Helmholtz2d_BC(u0,Au0) ! matrix-vec multiplication
    res = b - Au0  !r=b-Ax
    res_norm = norm(res)  ! ||r||

    Rerror   = res_norm / b_norm  ! scaled error
    beta(1) = res_norm       ! beta(1)=||r0||
    V(:,:,1)  = res / res_norm  !  This is V(:,1) i.e. v1 in the algorithm
    
    k=0
    do j=1,m_iter
      k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
      call arnoldi(V, H, k)
      call apply_givens_rotation(H, cs, sn, k)

      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  conjg(cs(k))*beta(k)

      Rerror = CDABS(beta(k+1))/b_norm

      if (my_id == 0 ) then
        write(*,"(I9,E17.9)")  k, Rerror
        open(1234, file=trim(logname),position='APPEND',status='OLD')
        write(1234, "(I9,E17.9)") k, Rerror
        close(1234)
      end if

      if (Rerror<eps) then
        exit
      end if
    end do

    call back_substitute(H,beta,k)

    do ki = 1,k
      u = u + beta(ki)*V(:,:,ki)
    enddo

    iter = k

    deallocate(sn, cs, beta)
    deallocate(res, u0, Au0)
    deallocate(H)
    deallocate(V)

  end subroutine fullgmres
  !====================================================================================

  !====================================================================================
  subroutine restartgmres(b,u,Rerror,iter_total)
    !! Restart GMRES
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter_total

    !-Local arguments ------------------------------------------------------------
    integer                        :: k, ki, j, Out_iter
    real(kind = realdp)                   :: b_norm, res0_norm, res_norm, pRerror
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res0, res, u0, Au0, Au, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res0(1-LAP:nx+LAP,1-LAP:ny+LAP), res(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(u0(1-LAP:nx+LAP,1-LAP:ny+LAP), Au0(1-LAP:nx+LAP,1-LAP:ny+LAP), Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)

    Out_iter=0
    Rerror=1.d0
    do while (Rerror > eps)
      Out_iter = Out_iter+1

      res  = (0.d0,0.d0)
      Au0  = (0.d0,0.d0)
      sn   = (0.d0,0.d0)
      cs   = (0.d0,0.d0)
      V    = (0.d0,0.d0)
      H    = (0.d0,0.d0)
      beta = (0.d0,0.d0)
      u0=u

      call Helmholtz2d_BC(u0,Au0) !-compute Ax
      res0 = b - Au0  !r=b-Ax
      res0_norm = norm(res0)  ! ||r||
      beta(1) = res0_norm       ! beta(1)=||r0||
      V(:,:,1)  = res0 / res0_norm  !  This is V(:,1) i.e. v1 in the algorithm

      k=0
      do j=1,m_iter
        k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
        call arnoldi(V, H, k)

        call apply_givens_rotation(H, cs, sn, k)

        beta(k+1) = -sn(k)*beta(k)
        beta(k)   =  conjg(cs(k))*beta(k)

        pRerror = CDABS(beta(k+1))/b_norm
      end do

      call back_substitute(H,beta,k)
      do ki = 1,k
        u = u + beta(ki)*V(:,:,ki)
      enddo

      ! compute Ax
      call Helmholtz2d_BC(u,Au)
      res = b - Au  !r=b-Ax
      res_norm = norm(res)  ! ||r||
      Rerror = res_norm / b_norm

      if (my_id == 0 ) then
        write(*,"(A,I9,A,E14.6,A,E14.6)") "   Outer Iteration   ", Out_iter, "   pError   ", pRerror, "   Error   ", Rerror
      end if

    enddo

    deallocate(sn, cs, beta)
    deallocate(res0, res, u0, Au0, Au)
    deallocate(H)
    deallocate(V)

    iter_total = Out_iter*m_iter
  end subroutine restartgmres
  !====================================================================================

  !====================================================================================
  subroutine Pre_fullgmres(b,u,Rerror,iter)
    !! Full GMRES with left preconditioned
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                        :: k,ki,j
    real(kind = realdp)                   :: b_norm, Mres_norm, Mb_norm, MRerror, res_norm
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res, Mres, Mb, u0, Au0, Au, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP), Mres(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Mb(1-LAP:nx+LAP,1-LAP:ny+LAP), u0(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Au0(1-LAP:nx+LAP,1-LAP:ny+LAP), Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b) ! ||b||_2
    Mb = Precond_x(b)  ! apply preconditioner, P^(-1)b

    Mb_norm = norm(Mb)

    res  = (0.d0,0.d0)
    Mres = (0.d0,0.d0)
    Au0  = (0.d0,0.d0)
    sn   = (0.d0,0.d0)
    cs   = (0.d0,0.d0)
    V    = (0.d0,0.d0)
    H    = (0.d0,0.d0)
    beta = (0.d0,0.d0)

    u0=u

    call Helmholtz2d_BC(u0,Au0) ! compute Ax
    res = b - Au0  !r=b-Ax
    Mres = Precond_x(res) ! M^{-1}r  !Mres = res

    Mres_norm = norm(Mres)  ! ||r||
    MRerror   = Mres_norm / Mb_norm  ! scaled error

    beta(1) = Mres_norm       ! beta(1)=||r0||
    V(:,:,1)  = Mres / Mres_norm  !  This is V(:,1) i.e. v1 in the algorithm

    k=0
    do j=1,m_iter
      k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
      call Prearnoldi(V, H, k)

      call apply_givens_rotation(H, cs, sn, k)

      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  conjg(cs(k))*beta(k)

      MRerror = CDABS(beta(k+1))/Mb_norm

      if (my_id == 0 ) then
        write(*,"(I9,E17.9)")  k, MRerror
        open(1234, file=trim(logname),position='APPEND',status='OLD')
        write(1234, "(I9,E17.9)") k, MRerror
        close(1234)
      end if

      if (MRerror < eps) then
        exit
      end if

    end do

    call back_substitute(H,beta,k)

    do ki = 1,k
      u = u + beta(ki)*V(:,:,ki)
    enddo

    call Helmholtz2d_BC(u,Au)

    res = b - Au
    res_norm = norm(res)
    Rerror   = res_norm / b_norm

    deallocate(sn, cs, beta)
    deallocate(res, Mres, Mb, u0, Au0, Au)
    deallocate(H)
    deallocate(V)

    iter = k
  end subroutine Pre_fullgmres
  !====================================================================================

  !====================================================================================
  subroutine full_pgmres(b,u,Rerror,iter)
    !! Full GMRES with right precondition
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                        :: k,ki,j,i
    real(kind = realdp)                   :: b_norm, res_norm
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res, z, u0, Au, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(u0(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(z(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)

    res  = (0.d0,0.d0)
    sn   = (0.d0,0.d0)
    cs   = (0.d0,0.d0)
    V    = (0.d0,0.d0)
    H    = (0.d0,0.d0)
    beta = (0.d0,0.d0)

    u0=u

    call Helmholtz2d_BC(u,Au) ! compute Ax
    res = b - Au  !r=b-Ax

    res_norm = norm(res)  ! ||r||
    Rerror   = res_norm / b_norm  

    beta(1) = res_norm       ! beta(1)=||r0||
    V(:,:,1)  = res / res_norm  !  This is V(:,1) i.e. v1 in the algorithm

    k=0
    do j=1,m_iter
      k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
      
      z = Precond_x(V(:,:,k))

      call Helmholtz2d_BC(z,V(:,:,k+1))

      do i=1,k
          H(i, k)  = dot_prod(V(:,:,i), V(:,:,k+1))     !-Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
          V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
      end do

      H(k+1, k) = norm(V(:,:,k+1))
      V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)


      call apply_givens_rotation(H, cs, sn, k)

      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  conjg(cs(k))*beta(k)

      Rerror = CDABS(beta(k+1))/b_norm

      if (my_id == 0 ) then
        write(*,"(I9,E17.9)")  k, Rerror
        open(1234, file=trim(logname),position='APPEND',status='OLD')
        write(1234, "(I9,E17.9)") k, Rerror
        close(1234)
      end if

      if (Rerror < eps) then
        exit
      end if

    end do

    call back_substitute(H,beta,k)

    u=(0.d0,0.d0)
    do ki = 1,k
      u = u + beta(ki)*V(:,:,ki)
    enddo
    u=Precond_x(u)
    u=u+u0

    call Helmholtz2d_BC(u,Au)

    res = b - Au
    res_norm = norm(res)
    Rerror   = res_norm / b_norm

    deallocate(sn, cs, beta)
    deallocate(res, z, u0, Au)
    deallocate(H)
    deallocate(V)

    iter = k
  end subroutine full_pgmres
  !====================================================================================

  !====================================================================================
  subroutine Pre_restartgmres(b,u,Rerror,iter)
    !! Preconditioned restart GMRES 
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real(kind = realdp), intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                        :: k,ki,j, Out_iter
    real(kind = realdp)                   :: b_norm, Mres_norm, Mb_norm, MRerror, res_norm
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res, Mres, Mb, u0, Au0, Au, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP), Mres(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Mb(1-LAP:nx+LAP,1-LAP:ny+LAP), u0(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Au0(1-LAP:nx+LAP,1-LAP:ny+LAP), Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)  ! what if b is real?
    Mb=Precond_x(b)
    Mb_norm = norm(Mb)

    Out_iter=0
    Rerror=1.d0
    MRerror=1.d0
    do while (MRerror .gt. eps)
      Out_iter = Out_iter+1

      res  = (0.d0,0.d0)
      Mres  = (0.d0,0.d0)
      Au0  = (0.d0,0.d0)
      sn   = (0.d0,0.d0)
      cs   = (0.d0,0.d0)
      V    = (0.d0,0.d0)
      H    = (0.d0,0.d0)
      beta = (0.d0,0.d0)

      u0=u

      call Helmholtz2d_BC(u0,Au0)
      Mres = b - Au0  !r=b-Ax
      Mres = Precond_x(Mres)
      Mres_norm = norm(Mres)  ! ||r||
      MRerror   = Mres_norm / Mb_norm  ! scaled error?
      beta(1) = Mres_norm       ! beta(1)=||r0||
      V(:,:,1)  = Mres / Mres_norm  !  This is V(:,1) i.e. v1 in the algorithm

      k=0
      do j=1,m_iter 
          !-Here m_iter is the iteration to restart
          k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
          call Prearnoldi(V, H, k)

          call apply_givens_rotation(H, cs, sn, k)

          beta(k+1) = -sn(k)*beta(k)
          beta(k)   =  conjg(cs(k))*beta(k)

          MRerror = CDABS(beta(k+1))/Mb_norm
      enddo

      call back_substitute(H,beta,k)

      do ki = 1,k
      u = u + beta(ki)*V(:,:,ki)
      enddo

      !call Helmholtz2d(u,Au)
      call Helmholtz2d_BC(u,Au)
      res = b - Au
      res_norm = norm(res)
      Rerror   = res_norm / b_norm

      if (my_id == 0 ) then
        write(*,"(A,I9,A,E14.6,A,E14.6)") "   Outer Iteration   ", Out_iter, "   MError   ", MRerror, "   Error   ", Rerror
      end if
      
    enddo

    deallocate(sn, cs, beta)
    deallocate(res, Mres, Mb, u0, Au0, Au)
    deallocate(H)
    deallocate(V)
    iter = Out_iter*m_iter
  end subroutine Pre_restartgmres
  !====================================================================================

  !================================================================
  subroutine arnoldi(V, H, k)
    !! ARNOLDI precess
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    integer, intent(in)    :: k
    complex(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1), intent(inout) :: V
    complex(kind = realdp), dimension(m_iter+1,m_iter),                    intent(inout) :: H

    integer :: i

    !-Subroutine content ---------------------------------------------------------

    !- w=A*v_i
    call Helmholtz2d_BC(V(:,:,k),V(:,:,k+1))

    do i=1,k
        H(i, k)  = dot_prod(V(:,:,i), V(:,:,k+1))     !-Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
        V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
    end do

    H(k+1, k) = norm(V(:,:,k+1))
    V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)

  end subroutine arnoldi
  !================================================================
  
  subroutine Prearnoldi(V, H, k)
    !! Preconditioned ARNOLDI precess
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    integer,                                                   intent(in)    :: k
    complex(kind = realdp),    dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1),   intent(inout) :: V
    complex(kind = realdp),    dimension(m_iter+1,m_iter),                      intent(inout) :: H

    integer                    :: i

    !-Subroutine content ---------------------------------------------------------
    !- w=A*v_i
    call Helmholtz2d_BC(V(:,:,k),V(:,:,k+1))
    !-Precondition
    V(:,:,k+1) = Precond_x(V(:,:,k+1))

    do i=1,k
        H(i, k)  = dot_prod(V(:,:,i), V(:,:,k+1))     !-Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
        V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
    end do

    H(k+1, k) = norm(V(:,:,k+1))
    V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)

  end subroutine Prearnoldi

  !================================================================
  subroutine apply_givens_rotation(H, cs, sn, k)
    !! APPLY GIVENS ROTATION 
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),    dimension(m_iter+1,m_iter), intent(inout)   :: H
    integer,                         intent(in)      :: k
    complex(kind = realdp),    dimension(m_iter+1),   intent(inout)   :: cs, sn

    !-Local arguments ------------------------------------------------------------
    complex(kind = realdp)    :: temp
    integer         :: i

    !-Subroutine content ---------------------------------------------------------
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

  end subroutine apply_givens_rotation


  !===============================================================
  subroutine back_substitute(H,beta,k)
    !! Apply BACK SUBSTITUTION 
    implicit none

    integer :: i
    integer, intent(in) :: k
    complex(kind = realdp), dimension(m_iter+1,m_iter), intent(in)    :: H
    complex(kind = realdp), dimension(m_iter+1),        intent(inout) :: beta

    beta(k) = beta(k)/H(k,k)

    do i=k-1,1,-1
      beta(i) = (beta(i) - sum(H(i,i+1:k)*beta(i+1:k)))/H(i,i)
    end do

  end subroutine back_substitute


  !====================================================================================
  subroutine full_pgcr(b,u,Rerror,iter)
    !! Right preconditioned GCR 
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                 :: k,i,j
    complex(kind = realdp)  :: alpha, beta
    real(kind = realdp)     :: b_norm, res_norm, vv_norm

    complex(kind = realdp), allocatable, dimension(:,:)   :: res, Au
    complex(kind = realdp), allocatable, dimension(:,:,:) :: vv,ss

    !-Subroutine content ---------------------------------------------------------
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(vv(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))
    allocate(ss(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)

    res  = (0.d0,0.d0)
    vv   = (0.d0,0.d0)
    ss   = (0.d0,0.d0)

    call Helmholtz2d_BC(u,Au) ! compute Ax
    res = b - Au  !r=b-Ax
    res_norm = norm(res)  ! ||r||
    Rerror   = res_norm / b_norm  ! scaled error?
    
    k=0
    do i=1,m_iter
      k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
      
      ss(:,:,i) = Precond_x(res)
      call Helmholtz2d_BC(ss(:,:,i),vv(:,:,i))

      do j=1,i-1
        alpha = dot_prod(vv(:,:,j),vv(:,:,i))
        ss(:,:,i) = ss(:,:,i) - alpha*ss(:,:,j)
        vv(:,:,i) = vv(:,:,i) - alpha*vv(:,:,j)
      end do

      vv_norm = norm(vv(:,:,i))
      ss(:,:,i)=ss(:,:,i)/vv_norm
      vv(:,:,i)=vv(:,:,i)/vv_norm
      beta=dot_prod(vv(:,:,i),res)
      u   = u + beta*ss(:,:,i)
      res = res - beta*vv(:,:,i)

      res_norm = norm(res)
      Rerror   = res_norm / b_norm
      if (my_id == 0 ) then
        write(*,"(I9,E17.9)")  k, Rerror
        open(1234, file=trim(logname),position='APPEND',status='OLD')
        write(1234, "(I9,E17.9)") k, Rerror
        close(1234)
      end if

      if (Rerror < eps) then
        exit
      end if

    end do

    call Helmholtz2d_BC(u,Au)
    res = b - Au
    res_norm = norm(res)
    Rerror   = res_norm / b_norm

    deallocate(res, Au)
    deallocate(vv)
    deallocate(ss)

    iter = k
  end subroutine full_pgcr
  !====================================================================================


  !====================================================================================
  subroutine pfgmres(b,u,Rerror,iter)
    !! Flexible GMRES with right precondition 
    implicit none

    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),   dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(in)    :: b
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP),intent(inout) :: u
    real   (kind = realdp),intent(out) :: Rerror
    integer, intent(inout) :: iter

    !-Local arguments ------------------------------------------------------------
    integer                        :: k,ki,j,i
    real(kind = realdp)                   :: b_norm, res_norm
    complex(kind = realdp), allocatable, dimension(:)     :: sn, cs, beta
    complex(kind = realdp), allocatable, dimension(:,:)   :: res, u0, Au, H
    complex(kind = realdp), allocatable, dimension(:,:,:) :: V, Z

    !-Subroutine content ---------------------------------------------------------
    allocate(sn(m_iter+1), cs(m_iter+1), beta(m_iter+1))
    allocate(res(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(u0(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(Au(1-LAP:nx+LAP,1-LAP:ny+LAP))
    allocate(H(m_iter+1,m_iter))
    allocate(V(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))
    allocate(Z(1-LAP:nx+LAP,1-LAP:ny+LAP,m_iter+1))

    b_norm = norm(b)

    res  = (0.d0,0.d0)
    sn   = (0.d0,0.d0)
    cs   = (0.d0,0.d0)
    V    = (0.d0,0.d0)
    H    = (0.d0,0.d0)
    beta = (0.d0,0.d0)

    u0=u

    call Helmholtz2d_BC(u,Au) ! compute Ax
    res = b - Au  !r=b-Ax

    res_norm = norm(res)  ! ||r||
    Rerror   = res_norm / b_norm  

    beta(1) = res_norm       ! beta(1)=||r0||
    V(:,:,1)  = res / res_norm  !  This is V(:,1) i.e. v1 in the algorithm

    k=0
    do j=1,m_iter
      k=k+1   !-Be careful!!, after the whole iteration without achieving eps, then the value of j will be "m_iter+1".So we need a k.
      
      Z(:,:,k) = Precond_x(V(:,:,k))
      call Helmholtz2d_BC(Z(:,:,k),V(:,:,k+1))

      do i=1,k
          H(i, k)  = dot_prod(V(:,:,i), V(:,:,k+1))     !-Attention: h_(i,j)=(w,v_i)=v^H*w, for complex value, so the code should be dot_product(v_i,w)
          V(:,:,k+1) = V(:,:,k+1) - H(i,k) * V(:,:,i)
      end do

      H(k+1, k) = norm(V(:,:,k+1))
      V(:,:,k+1)  = V(:,:,k+1) / H(k+1, k)


      call apply_givens_rotation(H, cs, sn, k)

      beta(k+1) = -sn(k)*beta(k)
      beta(k)   =  conjg(cs(k))*beta(k)

      Rerror = CDABS(beta(k+1))/b_norm

      if (my_id == 0 ) then
        write(*,"(I9,E17.9)")  k, Rerror
        open(1234, file=trim(logname),position='APPEND',status='OLD')
        write(1234, "(I9,E17.9)") k, Rerror
        close(1234)
      end if

      if (Rerror < eps) then
        exit
      end if

    end do

    call back_substitute(H,beta,k)

    do ki = 1,k
      u = u + beta(ki)*Z(:,:,ki)
    enddo

    call Helmholtz2d_BC(u,Au)
    res = b - Au
    res_norm = norm(res)
    Rerror   = res_norm / b_norm

    deallocate(sn, cs, beta)
    deallocate(res, u0, Au)
    deallocate(H)
    deallocate(V)
    deallocate(Z)

    iter = k
  end subroutine pfgmres
  !====================================================================================


  function Precond_x(x)
    !! This is a routine to select which preconditioner is applied, Precond_x = P^(-1)x
    implicit none
    !-Subroutine arguments -------------------------------------------------------
    complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: x
    complex(kind = realdp) :: Precond_x(1-LAP:nx+LAP,1-LAP:ny+LAP)

    select case (M_flag)
   
      case (1) 
        !-Multigrid based CSLP 
        Precond_x = MGCSLP_invMx(x) 

      case (2:5)
        !-Deflation preconditioning
        Precond_x = DEF_Px(x)

      case default
        !-No precondition
        Precond_x = x
       
    end select

  end function Precond_x

end module solvers