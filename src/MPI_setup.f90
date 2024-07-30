module mpi_setup
    !! Module for setting up MPI partition and data exchange
    use mpi
    use comm_variable
    implicit none

    contains
    integer function my_mod1(na,nb)
        !! A function to calculate the modulo of two integers na and nb
        implicit none
        integer :: na,nb

        if(na<0) then
            my_mod1=na+nb
        else if (na>nb-1) then
            my_mod1=na-nb
        else
            my_mod1=na
        endif
        return
    end function my_mod1

    !************************************************************************************************
    subroutine part2d()
        !! Partitions the 2D domain for parallel processing and calculates offsets and sizes for each partition
        implicit none
        integer :: k,ka,npx1_x,npx1_y,npx2_x,npx2_y,npy1_x,npy1_y,npy2_x,npy2_y

        if (np_size /= npx0*npy0) then
            if (my_id == 0) then
                write(*,*)'Number of Total Procs is ',np_size
                write(*,"('Number of Total procs is not equal to npx0*npy0=',I5,' !!!')")npx0*npy0
            endif
            stop
        endif

        !! the ID of partion along 2 directions
        npx=mod(my_id,npx0)
        npy=mod(my_id,npx0*npy0)/npx0

        nx=nx_global/npx0
        ny=ny_global/npy0
        if(npx < mod(nx_global,npx0)) then
            nx=nx+1
        endif
        if(npy < mod(ny_global,npy0)) then
            ny=ny+1
        endif

        do k=0,npx0-1
            ka=min(k,mod(nx_global,npx0))
            i_offset(k)=int(nx_global/npx0)*k+ka+1
            i_nn(k)=nx_global/npx0
            if(k < mod(nx_global,npx0)) then
                i_nn(k)=i_nn(k)+1
            endif
        enddo

        do k=0,npy0-1
            ka=min(k,mod(ny_global,npy0))
            j_offset(k)=int(ny_global/npy0)*k+ka+1
            j_nn(k)=ny_global/npy0
            if(k < mod(ny_global,npy0)) then
                j_nn(k)=j_nn(k)+1
            endif
        enddo

        npx1_x=my_mod1(npx-1,npx0)
        npx2_x=my_mod1(npx+1,npx0)
        npx1_y=npy
        npx2_y=npy
        ID_XM1=npx1_y*npx0+npx1_x    !! -1 proc in x-direction
        ID_XP1=npx2_y*npx0+npx2_x    !! +1 proc in x-direction
        if(Iperiodic_X == 0 .and. npx == 0) then
            ID_XM1=MPI_PROC_NULL     !! if not periodic, 0 node donot send mesg to npx0-1 node
        endif
        if(Iperiodic_X == 0 .and. npx == npx0-1) then
            ID_XP1=MPI_PROC_NULL    !! if not periodic, npx0-1 node donot send mesg to 0 node
        endif

        npy1_x=npx
        npy2_x=npx
        npy1_y=my_mod1(npy-1,npy0)
        npy2_y=my_mod1(npy+1,npy0)
        ID_YM1=npy1_y*npx0+npy1_x
        ID_YP1=npy2_y*npx0+npy2_x
        if(Iperiodic_Y== 0 .and. npy == 0) then
            ID_YM1=MPI_PROC_NULL
        endif
        if(Iperiodic_Y== 0 .and. npy == npy0-1) then
            ID_YP1=MPI_PROC_NULL
        endif

        call MPI_barrier(MPI_COMM_WORLD,ierr)

    end subroutine part2d

    !************************************************************************************************
    subroutine check_x2d(f)
        !! Exchanges data in the x direction, ONLY for the default (finest) grid system
        implicit none

        integer i,j,k1
        complex(kind = realdp):: f(1-LAP:nx+LAP,1-LAP:ny+LAP)
        complex(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        complex(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2

        allocate(tmp_send1(LAP*ny), tmp_send2(LAP*ny))
        allocate(tmp_recv1(LAP*ny), tmp_recv2(LAP*ny))

        k1=0
        do i=1,LAP
            do j=1,ny
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(nx-LAP+i,j)
            enddo
        enddo

        if (npx == 0) then
            k1=0
            do i=1,LAP
                do j=1,ny
                    k1=k1+1
                    tmp_send1(k1)=f(i+1,j)
                enddo
            enddo
        endif

        if (npx == npx0-1) then
            k1=0
            do i=1,LAP
                do j=1,ny
                    k1=k1+1
                    tmp_send2(k1)=f(nx-LAP+i-1,j)
                enddo
            enddo
        endif


        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_COMPLEX, ID_XM1, 9021, &
                            tmp_recv2, k1,  MPI_DOUBLE_COMPLEX, ID_XP1, 9021, MPI_COMM_WORLD,Status,ierr)
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_COMPLEX, ID_XP1, 8021, &
                                tmp_recv1, k1,  MPI_DOUBLE_COMPLEX, ID_XM1, 8021, MPI_COMM_WORLD,Status,ierr)

        if(ID_XM1 /= MPI_PROC_NULL) then
                k1=0
                do i=1,LAP
                    do j=1,ny
                        k1=k1+1
                        f(i-LAP,j)=tmp_recv1(k1)
                    enddo
                enddo
        endif

        if(ID_XP1 /= MPI_PROC_NULL) then
            k1=0
            do i=1,LAP
                do j=1,ny
                k1=k1+1
                f(nx+i,j)=tmp_recv2(k1)
                enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine check_x2d

    !************************************************************************************************
    subroutine check_y2d(f)
        !! Exchanges data in the y direction, ONLY for the default (finest) grid system
        implicit none

        integer i,j,k1 
        complex(kind = realdp):: f(1-LAP:nx+LAP,1-LAP:ny+LAP)
        complex(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        complex(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2

        allocate(tmp_send1(LAP*(nx+2*LAP)), tmp_send2(LAP*(nx+2*LAP)))
        allocate(tmp_recv1(LAP*(nx+2*LAP)), tmp_recv2(LAP*(nx+2*LAP)))

        k1=0
        do j=1,LAP
            do i=1-LAP,nx+LAP
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(i,ny+j-LAP)
            enddo
        enddo

        if (npy == 0) then
            k1=0
            do j=1,LAP
                do i=1-LAP,nx+LAP
                    k1=k1+1
                    tmp_send1(k1)=f(i,j+1)
                enddo
            enddo
        endif

        if (npy == npy0-1) then
            k1=0
            do j=1,LAP
                do i=1-LAP,nx+LAP
                    k1=k1+1
                    tmp_send2(k1)=f(i,ny+j-LAP-1)
                enddo
            enddo
        endif

        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_COMPLEX, ID_YM1, 9022, &
                            tmp_recv2, k1,  MPI_DOUBLE_COMPLEX, ID_YP1, 9022, MPI_COMM_WORLD,Status,ierr)
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_COMPLEX, ID_YP1, 022, &
                            tmp_recv1, k1,  MPI_DOUBLE_COMPLEX, ID_YM1, 022, MPI_COMM_WORLD,Status,ierr)

        if(ID_YM1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,nx+LAP
                        k1=k1+1
                        f(i,j-LAP)=tmp_recv1(k1)
                    enddo
            enddo
        endif
        if(ID_YP1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,nx+LAP    !!no2
                        k1=k1+1
                        f(i,ny+j)=tmp_recv2(k1)
                    enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine check_y2d
    !************************************************************************************************

    subroutine check_xy2d(f)
        !! Exchanges data first in the x direction and then y direction, ONLY for the default (finest) grid system
        implicit none
        complex(kind = realdp) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP)

        call check_x2d(f)
        call check_y2d(f)

    end subroutine check_xy2d
    !************************************************************************************************

    subroutine mg_check_x2d(f,ni,nj)
        !! Exchanges data in the x direction, for a specified-size grid system
        implicit none

        integer,intent(in) :: ni,nj
        complex(kind = realdp):: f(1-LAP:ni+LAP,1-LAP:nj+LAP)
        complex(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        complex(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2
        integer i,j,k1 

        allocate(tmp_send1(LAP*nj),tmp_send2(LAP*nj))
        allocate(tmp_recv1(LAP*nj),tmp_recv2(LAP*nj))

        k1=0
        do i=1,LAP
            do j=1,nj
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(ni-LAP+i,j)
            enddo
        enddo

        if (npx == 0) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    tmp_send1(k1)=f(i+1,j)
                enddo
            enddo
        endif

        if (npx == npx0-1) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    tmp_send2(k1)=f(ni-LAP+i-1,j)
                enddo
            enddo
        endif

        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_COMPLEX, ID_XM1, 9421, &
                            tmp_recv2, k1,  MPI_DOUBLE_COMPLEX, ID_XP1, 9421, MPI_COMM_WORLD,Status,ierr)   !!no2
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_COMPLEX, ID_XP1, 421, &
                                tmp_recv1, k1,  MPI_DOUBLE_COMPLEX, ID_XM1, 421, MPI_COMM_WORLD,Status,ierr)   !!no2

        if(ID_XM1 /= MPI_PROC_NULL) then
                k1=0
                do i=1,LAP
                    do j=1,nj
                        k1=k1+1
                        f(i-LAP,j)=tmp_recv1(k1)
                    enddo
                enddo
        endif


        if(ID_XP1 /= MPI_PROC_NULL) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    f(ni+i,j)=tmp_recv2(k1)
                enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine mg_check_x2d

    !************************************************************************************************

    subroutine mg_check_y2d(f,ni,nj)
        !! Exchanges data in the y direction, for a specified-size grid system
        implicit none

        integer,intent(in) :: ni,nj
        complex(kind = realdp):: f(1-LAP:ni+LAP,1-LAP:nj+LAP)
        complex(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        complex(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2
        integer i,j,k1 

        allocate(tmp_send1(LAP*(ni+2*LAP)),tmp_send2(LAP*(ni+2*LAP)))
        allocate(tmp_recv1(LAP*(ni+2*LAP)),tmp_recv2(LAP*(ni+2*LAP)))

        k1=0
        do j=1,LAP
            do i=1-LAP,ni+LAP
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(i,nj+j-LAP)
            enddo
        enddo

        if (npy == 0) then
            k1=0
            do j=1,LAP
                do i=1-LAP,ni+LAP
                    k1=k1+1
                    tmp_send1(k1)=f(i,j+1)
                enddo
            enddo
        endif

        if (npy == npy0-1) then
            k1=0
            do j=1,LAP
                do i=1-LAP,ni+LAP
                    k1=k1+1
                    tmp_send2(k1)=f(i,nj+j-LAP-1)
                enddo
            enddo
        endif

        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_COMPLEX, ID_YM1, 9422, &
                            tmp_recv2, k1,  MPI_DOUBLE_COMPLEX, ID_YP1, 9422, MPI_COMM_WORLD,Status,ierr)   !!no2
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_COMPLEX, ID_YP1, 8422, &
                            tmp_recv1, k1,  MPI_DOUBLE_COMPLEX, ID_YM1, 8422, MPI_COMM_WORLD,Status,ierr)   !!no2

        if(ID_YM1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,ni+LAP
                        k1=k1+1
                        f(i,j-LAP)=tmp_recv1(k1)
                    enddo
            enddo
        endif

        if(ID_YP1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,ni+LAP
                        k1=k1+1
                        f(i,nj+j)=tmp_recv2(k1)
                    enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine mg_check_y2d

    !************************************************************************************************

    subroutine mg_check_xy2d(f,ni,nj)
        !! Exchanges data first in the x direction and then y direction, for a specified-size grid system
        implicit none
        integer,intent(in) :: ni,nj
        complex(kind = realdp) :: f(1-LAP:ni+LAP,1-LAP:nj+LAP)

        call mg_check_x2d(f,ni,nj)
        call mg_check_y2d(f,ni,nj)

    end subroutine mg_check_xy2d

    !************************************************************************************************

    subroutine mg_checkreal_x2d(f,ni,nj)
        !! Exchanges REAL data in the x direction, for a specified-size grid system
        implicit none

        integer,intent(in) :: ni,nj
        real(kind = realdp):: f(1-LAP:ni+LAP,1-LAP:nj+LAP)
        real(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        real(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2
        integer i,j,k1 

        allocate(tmp_send1(LAP*nj),tmp_send2(LAP*nj))
        allocate(tmp_recv1(LAP*nj),tmp_recv2(LAP*nj))

        k1=0
        do i=1,LAP
            do j=1,nj
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(ni-LAP+i,j)
            enddo
        enddo

        if (npx == 0) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    tmp_send1(k1)=f(i+1,j)
                enddo
            enddo
        endif

        if (npx == npx0-1) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    tmp_send2(k1)=f(ni-LAP+i-1,j)
                enddo
            enddo
        endif

        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_PRECISION, ID_XM1, 9421, &
                            tmp_recv2, k1,  MPI_DOUBLE_PRECISION, ID_XP1, 9421, MPI_COMM_WORLD,Status,ierr)   !!no2
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_PRECISION, ID_XP1, 421, &
                                tmp_recv1, k1,  MPI_DOUBLE_PRECISION, ID_XM1, 421, MPI_COMM_WORLD,Status,ierr)   !!no2

        if(ID_XM1 /= MPI_PROC_NULL) then
                k1=0
                do i=1,LAP
                    do j=1,nj
                        k1=k1+1
                        f(i-LAP,j)=tmp_recv1(k1)
                    enddo
                enddo
        endif


        if(ID_XP1 /= MPI_PROC_NULL) then
            k1=0
            do i=1,LAP
                do j=1,nj
                    k1=k1+1
                    f(ni+i,j)=tmp_recv2(k1)
                enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine mg_checkreal_x2d

    !************************************************************************************************

    subroutine mg_checkreal_y2d(f,ni,nj)
        !! Exchanges REAL data in the y direction, for a specified-size grid system
        implicit none

        integer,intent(in) :: ni,nj
        real(kind = realdp):: f(1-LAP:ni+LAP,1-LAP:nj+LAP)
        real(kind = realdp), allocatable, dimension(:) :: tmp_send1, tmp_send2
        real(kind = realdp), allocatable, dimension(:) :: tmp_recv1, tmp_recv2
        integer i,j,k1 

        allocate(tmp_send1(LAP*(ni+2*LAP)),tmp_send2(LAP*(ni+2*LAP)))
        allocate(tmp_recv1(LAP*(ni+2*LAP)),tmp_recv2(LAP*(ni+2*LAP)))

        k1=0
        do j=1,LAP
            do i=1-LAP,ni+LAP
                k1=k1+1
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(i,nj+j-LAP)
            enddo
        enddo

        if (npy == 0) then
            k1=0
            do j=1,LAP
                do i=1-LAP,ni+LAP
                    k1=k1+1
                    tmp_send1(k1)=f(i,j+1)
                enddo
            enddo
        endif

        if (npy == npy0-1) then
            k1=0
            do j=1,LAP
                do i=1-LAP,ni+LAP
                    k1=k1+1
                    tmp_send2(k1)=f(i,nj+j-LAP-1)
                enddo
            enddo
        endif

        call MPI_Sendrecv(tmp_send1, k1,  MPI_DOUBLE_PRECISION, ID_YM1, 9422, &
                            tmp_recv2, k1,  MPI_DOUBLE_PRECISION, ID_YP1, 9422, MPI_COMM_WORLD,Status,ierr)   !!no2
        call MPI_Sendrecv(tmp_send2, k1,  MPI_DOUBLE_PRECISION, ID_YP1, 8422, &
                            tmp_recv1, k1,  MPI_DOUBLE_PRECISION, ID_YM1, 8422, MPI_COMM_WORLD,Status,ierr)   !!no2

        if(ID_YM1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,ni+LAP
                        k1=k1+1
                        f(i,j-LAP)=tmp_recv1(k1)
                    enddo
            enddo
        endif

        if(ID_YP1 /= MPI_PROC_NULL) then
            k1=0
            do j=1,LAP
                    do i=1-LAP,ni+LAP
                        k1=k1+1
                        f(i,nj+j)=tmp_recv2(k1)
                    enddo
            enddo
        endif

        deallocate(tmp_send1, tmp_send2)
        deallocate(tmp_recv1, tmp_recv2)

    end subroutine mg_checkreal_y2d

    !************************************************************************************************

    subroutine mg_checkreal_xy2d(f,ni,nj)
        !! Exchanges REAL data first in the x direction and then y direction, for a specified-size grid system
        implicit none
        integer,intent(in) :: ni,nj
        real(kind = realdp) :: f(1-LAP:ni+LAP,1-LAP:nj+LAP)

        call mg_checkreal_x2d(f,ni,nj)
        call mg_checkreal_y2d(f,ni,nj)

    end subroutine mg_checkreal_xy2d

end module mpi_setup