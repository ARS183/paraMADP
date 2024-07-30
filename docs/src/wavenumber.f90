module wavenumber
    !! This module is used to determine the wavenumber of the computational domain, as well as (kh)^2
    use mpi
    use comm_variable
    use mpi_setup
    implicit none

    real(kind = realdp), allocatable, dimension(:,:)     :: wavenumber_k
    real(kind = realdp), allocatable, dimension(:,:)     :: kh_pow_2
    ! If the wavenumber_k is not initialized to zero, some unexpected values will occur in parallel computing. 
    ! It will leads to crash. Ex. MP3 nx2305f40 np192

    contains

    subroutine Const_K()
        !! constant wavenumber, determine by input variable k0
        implicit none

        allocate(wavenumber_k(1-LAP:nx+LAP,1-LAP:ny+LAP))
        allocate(kh_pow_2(1-LAP:nx+LAP,1-LAP:ny+LAP))

        wavenumber_k = 0.d0
        kh_pow_2 = 0.d0

        wavenumber_k = k0
        kh_pow_2 = k0*k0*hxhy

        !! data communication with neighbouring subdomains
        call mg_checkreal_xy2d(wavenumber_k, nx, ny)
        call mg_checkreal_xy2d(kh_pow_2, nx, ny)

    end subroutine Const_K

    subroutine wavenumber_k_Wedge(xx,yy)
        !! Wavenumber for Wedge model problem
        implicit none

        real(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        real(kind = realdp) :: c1, c2, c3
            !! Three acoustis velocity for three layers
        integer :: i,j

        allocate(wavenumber_k(1-LAP:nx+LAP,1-LAP:ny+LAP))
        wavenumber_k = 0.d0

        c1 = 2000.d0 ! m/s
        c2 = 1500.d0 ! m/s
        c3 = 3000.d0 ! m/s

        do j=1,ny
            do i=1,nx
                if (yy(i,j) > -1.d0/6.d0*xx(i,j)-400.d0) then
                    !! The top layer
                    wavenumber_k(i,j) = 2.d0*pi*freq/c1
                elseif (yy(i,j) < 1.d0/3.d0*xx(i,j)-800.d0) then
                    !! The bottom layer
                    wavenumber_k(i,j) = 2.d0*pi*freq/c3
                else
                    !! The middle layer
                    wavenumber_k(i,j) = 2.d0*pi*freq/c2
                endif
            enddo
        enddo

        allocate(kh_pow_2(1-LAP:nx+LAP,1-LAP:ny+LAP))
        kh_pow_2 = 0.d0
        kh_pow_2 = wavenumber_k*wavenumber_k*hxhy

        call mg_checkreal_xy2d(wavenumber_k, nx, ny)
        call mg_checkreal_xy2d(kh_pow_2, nx, ny)

    end subroutine wavenumber_k_Wedge

    subroutine read_wavenumber_k_marmousi()
        !! Read the velocity profile of Marmousi problem with specified grid size and calculate the wavenumber
        implicit none

        real(kind=realdp),allocatable,dimension(:,:) :: vp_global
            !! Global velocity profile 
        real(kind=realdp),allocatable,dimension(:,:) :: temp_vp
            !! Local velocity profile for subdomains
        integer i,j,i_global,j_global,ii,jj,k3
        real(kind=realdp) :: temp

        allocate(wavenumber_k(1-LAP:nx+LAP,1-LAP:ny+LAP))
        wavenumber_k = 0.d0

        allocate(vp_global(1:nx_global,1:ny_global))

        if (my_id .eq. 0) then
            ! Rank 0: read the velocity data from files starting with "vp_" in Input directory, with specified grid size
            write(filename,"('Input/vp_nx',I5.5,'ny',I4.4,'.plt')") nx_global, ny_global
            open(77,file=trim(filename),status='old')

            read(77,*)
            read(77,*)

                do jj=1,ny_global
                    do ii=1,nx_global
                            read(77,*) temp, temp, vp_global(ii,jj)
                    enddo
                enddo

            close(77)
            !========must remember that the data you read is velocity, ===================!
            !========don't forget to transfer them to wavenumber within every processor===!
            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i

                    if (k3 .eq. 0) then
                        ! Rank 0
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                wavenumber_k(ii,jj)  = 2.d0*pi*freq/vp_global(i_global,j_global) 
                            enddo
                        enddo
                    else
                        ! Send to the other ranks
                        allocate(temp_vp(i_nn(i),j_nn(j)))
                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                temp_vp(ii,jj) = vp_global(i_global,j_global)
                            enddo
                        enddo
                        call MPI_SEND(temp_vp,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,9200,MPI_COMM_WORLD,ierr)
                        deallocate(temp_vp)
                    endif
                enddo
            enddo
        else
            ! The other ranks:  receive data and compute wavenumber
            allocate(temp_vp(i_nn(npx),j_nn(npy)))
            call MPI_RECV(temp_vp,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,9200,MPI_COMM_WORLD,status,ierr)
            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    wavenumber_k(ii,jj)  = 2.d0*pi*freq/temp_vp(ii,jj)
                enddo
            enddo
            deallocate(temp_vp)
        endif

    deallocate(vp_global)

    allocate(kh_pow_2(1-LAP:nx+LAP,1-LAP:ny+LAP))
    kh_pow_2 = 0.d0
    kh_pow_2 = wavenumber_k*wavenumber_k*hxhy

    call mg_checkreal_xy2d(wavenumber_k, nx, ny)
    call mg_checkreal_xy2d(kh_pow_2, nx, ny)

    end subroutine read_wavenumber_k_marmousi

    subroutine wavenumber_k_destroy()
        !! Deallocate the wavenumber field in the end
        implicit none

        deallocate(wavenumber_k,kh_pow_2)

    end subroutine wavenumber_k_destroy

end module