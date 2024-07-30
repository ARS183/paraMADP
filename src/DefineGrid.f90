module define_grid
    !! This is a module to determine the grid coordinate. 
    !! First determined by Rank 0 and then distributed to the other ranks. 
    !! For now, only rectangular domain with uniform space step size, i.e. hx==hy. 
    use mpi
    use comm_variable

    implicit none

contains
    subroutine define_uniform_grid(phys1,phys2)
        !! A square computational domain with a unit side length
        implicit none
        real(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: phys1,phys2
            !! Local coordinates in x and y directions

        real(kind = realdp), allocatable, dimension(:,:) :: phys1_global,phys2_global
            !! Global coordinates in x and y directions
        real(kind = realdp), allocatable, dimension(:,:) :: temp_phys1,temp_phys2
            !! Temp local coordinates in x and y directions for data sending
        integer :: i,j,ka,i_global,j_global,ii,jj

        hy=sly/dble(ny_global-1)
        hx=slx/dble(nx_global-1)

        if (my_id .eq. 0) then
            allocate(phys1_global(1:nx_global,1:ny_global), phys2_global(1:nx_global,1:ny_global))

            do j=1,ny_global
                do i=1,nx_global
                    phys1_global(i,j)=0.d0+dble(i-1)*hx
                    phys2_global(i,j)=0.d0+dble(j-1)*hy
                enddo
            enddo

            do j=0,npy0-1
                do i=0,npx0-1
                    ka=j*npx0+i

                    if (ka .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                phys1(ii,jj)=phys1_global(i_global,j_global)
                                phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                    else
                        allocate(temp_phys1(i_nn(i),j_nn(j)),temp_phys2(i_nn(i),j_nn(j)))

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                temp_phys1(ii,jj)=phys1_global(i_global,j_global)
                                temp_phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                        call MPI_SEND(temp_phys1,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,1,MPI_COMM_WORLD,ierr)
                        call MPI_SEND(temp_phys2,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,2,MPI_COMM_WORLD,ierr)

                        deallocate(temp_phys1,temp_phys2)

                    endif
                enddo
            enddo

            deallocate(phys1_global, phys2_global)

        else !(my_id .eq. 0)

            allocate(temp_phys1(i_nn(npx),j_nn(npy)),temp_phys2(i_nn(npx),j_nn(npy)))

            call MPI_RECV(temp_phys1,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,1,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(temp_phys2,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,2,MPI_COMM_WORLD,status,ierr)

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    phys1(ii,jj)=temp_phys1(ii,jj)
                    phys2(ii,jj)=temp_phys2(ii,jj)
                enddo
            enddo

            deallocate(temp_phys1,temp_phys2)

        endif 

    end subroutine define_uniform_grid


    subroutine define_wedge_grid(phys1,phys2)
        !! Wedge model problem
        implicit none
        real(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: phys1,phys2

        real(kind = realdp), allocatable, dimension(:,:) :: phys1_global,phys2_global
        real(kind = realdp), allocatable, dimension(:,:) :: temp_phys1,temp_phys2
        integer :: i,j,ka,i_global,j_global,ii,jj

        hy=sly/dble(ny_global-1)
        hx=slx/dble(nx_global-1)

        if (my_id .eq. 0) then
            allocate(phys1_global(1:nx_global,1:ny_global), phys2_global(1:nx_global,1:ny_global))

            do j=1,ny_global
                do i=1,nx_global
                    phys1_global(i,j)=0.d0+dble(i-1)*hx
                    phys2_global(i,j)=-sly+dble(j-1)*hy
                enddo
            enddo

            do j=0,npy0-1
                do i=0,npx0-1
                    ka=j*npx0+i

                    if (ka .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                phys1(ii,jj)=phys1_global(i_global,j_global)
                                phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                    else
                        allocate(temp_phys1(i_nn(i),j_nn(j)),temp_phys2(i_nn(i),j_nn(j)))

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                temp_phys1(ii,jj)=phys1_global(i_global,j_global)
                                temp_phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                        call MPI_SEND(temp_phys1,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,1,MPI_COMM_WORLD,ierr)
                        call MPI_SEND(temp_phys2,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,2,MPI_COMM_WORLD,ierr)

                        deallocate(temp_phys1,temp_phys2)

                    endif
                enddo
            enddo

            deallocate(phys1_global, phys2_global)

        else !(my_id .eq. 0)

            allocate(temp_phys1(i_nn(npx),j_nn(npy)),temp_phys2(i_nn(npx),j_nn(npy)))

            call MPI_RECV(temp_phys1,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,1,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(temp_phys2,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,2,MPI_COMM_WORLD,status,ierr)

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    phys1(ii,jj)=temp_phys1(ii,jj)
                    phys2(ii,jj)=temp_phys2(ii,jj)
                enddo
            enddo

            deallocate(temp_phys1,temp_phys2)

        endif !(my_id .eq. 0)

    end subroutine define_wedge_grid


    subroutine define_marmousi_grid(phys1,phys2)
        !! Marmousi problem
        implicit none
        real(kind = realdp), dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: phys1,phys2

        real(kind = realdp), allocatable, dimension(:,:) :: phys1_global,phys2_global
        real(kind = realdp), allocatable, dimension(:,:) :: temp_phys1,temp_phys2
        integer :: i,j,ka,i_global,j_global,ii,jj

        hy=sly/dble(ny_global-1)
        hx=slx/dble(nx_global-1)

        if (my_id .eq. 0) then
            allocate(phys1_global(1:nx_global,1:ny_global), phys2_global(1:nx_global,1:ny_global))

            do j=1,ny_global
                do i=1,nx_global
                    phys1_global(i,j)=0.d0+dble(i-1)*hx
                    phys2_global(i,j)=-sly+dble(j-1)*hy
                enddo
            enddo

            do j=0,npy0-1
                do i=0,npx0-1
                    ka=j*npx0+i

                    if (ka .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                phys1(ii,jj)=phys1_global(i_global,j_global)
                                phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                    else
                        allocate(temp_phys1(i_nn(i),j_nn(j)),temp_phys2(i_nn(i),j_nn(j)))

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                temp_phys1(ii,jj)=phys1_global(i_global,j_global)
                                temp_phys2(ii,jj)=phys2_global(i_global,j_global)
                            enddo
                        enddo

                        call MPI_SEND(temp_phys1,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,1,MPI_COMM_WORLD,ierr)
                        call MPI_SEND(temp_phys2,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,  &
                            ka,2,MPI_COMM_WORLD,ierr)

                        deallocate(temp_phys1,temp_phys2)

                    endif
                enddo
            enddo

            deallocate(phys1_global, phys2_global)

        else !(my_id .eq. 0)

            allocate(temp_phys1(i_nn(npx),j_nn(npy)),temp_phys2(i_nn(npx),j_nn(npy)))

            call MPI_RECV(temp_phys1,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,1,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(temp_phys2,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,  &
                0,2,MPI_COMM_WORLD,status,ierr)

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    phys1(ii,jj)=temp_phys1(ii,jj)
                    phys2(ii,jj)=temp_phys2(ii,jj)
                enddo
            enddo

            deallocate(temp_phys1,temp_phys2)

        endif !(my_id .eq. 0)

    end subroutine define_marmousi_grid

end module define_grid

