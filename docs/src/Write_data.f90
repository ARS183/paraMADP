module write_data
    !! A module for solutions output. 
    !! The idea is to collect the data from the other ranks to rank 0, and then write into a tecplot file (.plt). 
    ! There are three "repeated" variants so far:
    ! one for the numerical solution, one for the exact/analytical solution, one for data with real type (wavenumber)
    ! To improve: an version that can handle different tpyes of input data and can manually specify the output filename.
    use mpi
    use comm_variable
    implicit none

contains
    subroutine write_data_whole(xx,yy,u)
        !! Write the solutions
        implicit none
        integer :: i,j,k3,i_global,j_global,ii,jj
        real(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) ::u

        real(kind = realdp),   allocatable,dimension(:,:) :: xx_global,yy_global
        complex(kind = realdp),allocatable,dimension(:,:) :: u_global

        real(kind = realdp),   allocatable,dimension(:,:) :: temp_xx,temp_yy
        complex(kind = realdp),allocatable,dimension(:,:) :: temp_u

        allocate(xx_global(1:nx_global,1:ny_global),yy_global(1:nx_global,1:ny_global))
        allocate(u_global(1:nx_global,1:ny_global))

        if (my_id .eq. 0 ) then

            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                xx_global(i_global,j_global)=xx(ii,jj)
                                yy_global(i_global,j_global)=yy(ii,jj)
                            enddo
                        enddo

                    else
                        allocate(temp_xx(i_nn(i),j_nn(j)),temp_yy(i_nn(i),j_nn(j)))

                        call MPI_RECV(temp_xx,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,1,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(temp_yy,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,2,MPI_COMM_WORLD,status,ierr)

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                xx_global(i_global,j_global)=temp_xx(ii,jj)
                                yy_global(i_global,j_global)=temp_yy(ii,jj)
                            enddo
                        enddo

                        deallocate(temp_xx,temp_yy)
                    endif
                enddo
            enddo



            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                u_global(i_global,j_global)    = u(ii,jj)
                            enddo
                        enddo
                    else
                        allocate(temp_u(i_nn(i),j_nn(j)))
                        call MPI_RECV(temp_u,i_nn(i)*j_nn(j),MPI_DOUBLE_COMPLEX,k3,4,MPI_COMM_WORLD,status,ierr)
                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                u_global(i_global,j_global)    =    temp_u(ii,jj)
                            enddo
                        enddo
                        deallocate(temp_u)
                    endif
                enddo
            enddo

            !===output plt for whole-area
            write(filename,"('Output/MP',I1,'BC',I1,'nx',I5.5,'k',I5.5,     &
                           & 'thm',I2.2,'pre',I1,'nldef',I1,                &
                           & 'np',I3.3,'.plt')")                            &
                           & i_case,flag_BCs,nx_global,int(k0),             &
                           & Algorithm,M_flag,def_nlevel,                   &
                           & npx0*npy0
            open(70,file=trim(filename))
            write(70,*)'variables=x,y,u'
            write(70,*)'zone i=',nx_global,' j=',ny_global
            !write(70,*)'solutiontime=',iter*dt,'strandid=',strandid

            do jj=1,ny_global
                do ii=1,nx_global
                    write(70,*) xx_global(ii,jj),yy_global(ii,jj),real(u_global(ii,jj))
                enddo
            enddo

            close(70)

        else  !if (my_id .eq. 0 )

            allocate(temp_xx(i_nn(npx),j_nn(npy)),temp_yy(i_nn(npx),j_nn(npy)), &
                             temp_u(i_nn(npx),j_nn(npy)))

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_xx(ii,jj)=xx(ii,jj)
                    temp_yy(ii,jj)=yy(ii,jj)
                enddo
            enddo

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_u(ii,jj)=u(ii,jj)
                enddo
            enddo

            call MPI_SEND(temp_xx,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_yy,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_u, i_nn(npx)*j_nn(npy),MPI_DOUBLE_COMPLEX,0,4,MPI_COMM_WORLD,ierr)

            deallocate(temp_xx,temp_yy,temp_u)

        endif  !if (my_id .eq. 0 )
        deallocate(xx_global, yy_global)
        deallocate(u_global)

    end subroutine write_data_whole

    subroutine write_exact_data_whole(xx,yy,u)
        !! Write the analytical solutions
        implicit none
        integer :: i,j,k3,i_global,j_global,ii,jj
        real(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        complex(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) ::u

        real(kind = realdp),   allocatable,dimension(:,:) :: xx_global,yy_global
        complex(kind = realdp),allocatable,dimension(:,:) :: u_global

        real(kind = realdp),   allocatable,dimension(:,:) :: temp_xx,temp_yy
        complex(kind = realdp),allocatable,dimension(:,:) :: temp_u

        allocate(xx_global(1:nx_global,1:ny_global),yy_global(1:nx_global,1:ny_global))
        allocate(u_global(1:nx_global,1:ny_global))

        if (my_id .eq. 0 ) then

            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                xx_global(i_global,j_global)=xx(ii,jj)
                                yy_global(i_global,j_global)=yy(ii,jj)
                            enddo
                        enddo

                    else
                        allocate(temp_xx(i_nn(i),j_nn(j)),temp_yy(i_nn(i),j_nn(j)))

                        call MPI_RECV(temp_xx,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,10,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(temp_yy,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,20,MPI_COMM_WORLD,status,ierr)

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                xx_global(i_global,j_global)=temp_xx(ii,jj)
                                yy_global(i_global,j_global)=temp_yy(ii,jj)
                            enddo
                        enddo

                        deallocate(temp_xx,temp_yy)
                    endif
                enddo
            enddo



            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                u_global(i_global,j_global)    = u(ii,jj)
                            enddo
                        enddo
                    else
                        allocate(temp_u(i_nn(i),j_nn(j)))
                        call MPI_RECV(temp_u,i_nn(i)*j_nn(j),MPI_DOUBLE_COMPLEX,k3,40,MPI_COMM_WORLD,status,ierr)
                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                u_global(i_global,j_global)    =    temp_u(ii,jj)
                            enddo
                        enddo
                        deallocate(temp_u)
                    endif
                enddo
            enddo

            !===output plt for whole-area
            write(filename,"('Output/MP',I1,'_BC',I1, '_exact_nx',I5.5,'k',I5.5,'thm',I2.2,'np',I3.3,'.plt')") &
                            int(i_case),int(flag_BCs),int(nx_global),int(k0),Algorithm,npx0*npy0
            open(71,file=trim(filename))
            write(71,*)'variables=x,y,u'
            write(71,*)'zone i=',nx_global,' j=',ny_global
            !write(71,*)'solutiontime=',iter*dt,'strandid=',strandid

            do jj=1,ny_global
                do ii=1,nx_global
                    write(71,*) xx_global(ii,jj),yy_global(ii,jj),real(u_global(ii,jj))
                enddo
            enddo

            close(71)

        else  !if (my_id .eq. 0 )

            allocate(temp_xx(i_nn(npx),j_nn(npy)),temp_yy(i_nn(npx),j_nn(npy)), &
                             temp_u(i_nn(npx),j_nn(npy)))

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_xx(ii,jj)=xx(ii,jj)
                    temp_yy(ii,jj)=yy(ii,jj)
                enddo
            enddo

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_u(ii,jj)=u(ii,jj)
                enddo
            enddo

            call MPI_SEND(temp_xx,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,10,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_yy,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,20,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_u, i_nn(npx)*j_nn(npy),MPI_DOUBLE_COMPLEX,0,40,MPI_COMM_WORLD,ierr)

            deallocate(temp_xx,temp_yy,temp_u)

        endif  !if (my_id .eq. 0 )
        deallocate(xx_global, yy_global)
        deallocate(u_global)

    end subroutine write_exact_data_whole

    subroutine write_real_data_whole(xx,yy,u)
        !! Write the real wavenumber/variables
        implicit none
        integer :: i,j,k3,i_global,j_global,ii,jj
        real(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
        real(kind = realdp),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) ::u

        real(kind = realdp),   allocatable,dimension(:,:) :: xx_global,yy_global
        real(kind = realdp),allocatable,dimension(:,:) :: u_global

        real(kind = realdp),   allocatable,dimension(:,:) :: temp_xx,temp_yy
        real(kind = realdp),allocatable,dimension(:,:) :: temp_u

        allocate(xx_global(1:nx_global,1:ny_global),yy_global(1:nx_global,1:ny_global))
        allocate(u_global(1:nx_global,1:ny_global))

        if (my_id .eq. 0 ) then

            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                xx_global(i_global,j_global)=xx(ii,jj)
                                yy_global(i_global,j_global)=yy(ii,jj)
                            enddo
                        enddo

                    else
                        allocate(temp_xx(i_nn(i),j_nn(j)),temp_yy(i_nn(i),j_nn(j)))

                        call MPI_RECV(temp_xx,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,10,MPI_COMM_WORLD,status,ierr)
                        call MPI_RECV(temp_yy,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,20,MPI_COMM_WORLD,status,ierr)

                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                xx_global(i_global,j_global)=temp_xx(ii,jj)
                                yy_global(i_global,j_global)=temp_yy(ii,jj)
                            enddo
                        enddo

                        deallocate(temp_xx,temp_yy)
                    endif
                enddo
            enddo



            do j=0,npy0-1
                do i=0,npx0-1
                    k3=j*npx0+i
                    if (k3 .eq. 0) then
                        do jj=1,ny
                            do ii=1,nx
                                i_global=i_offset(i)-1+ii
                                j_global=j_offset(j)-1+jj
                                u_global(i_global,j_global)    = u(ii,jj)
                            enddo
                        enddo
                    else
                        allocate(temp_u(i_nn(i),j_nn(j)))
                        call MPI_RECV(temp_u,i_nn(i)*j_nn(j),MPI_DOUBLE_PRECISION,k3,40,MPI_COMM_WORLD,status,ierr)
                        do jj=1,j_nn(j)
                            do ii=1,i_nn(i)
                                i_global=i_offset(i)+ii-1
                                j_global=j_offset(j)+jj-1
                                u_global(i_global,j_global)    =    temp_u(ii,jj)
                            enddo
                        enddo
                        deallocate(temp_u)
                    endif
                enddo
            enddo

            !===output plt for whole-area
            write(filename,"('Output/MP',I1,'_BC',I1, '_wn_nx',I5.5,'k',I5.5,'thm',I2.2,'np',I3.3,'.plt')") &
                            int(i_case),int(flag_BCs),int(nx_global),int(k0),Algorithm,npx0*npy0
            open(71,file=trim(filename))
            write(71,*)'variables=x,y,k'
            write(71,*)'zone i=',nx_global,' j=',ny_global

            do jj=1,ny_global
                do ii=1,nx_global
                    write(71,*) xx_global(ii,jj),yy_global(ii,jj),real(u_global(ii,jj))
                enddo
            enddo

            close(71)

        else  !if (my_id .eq. 0 )

            allocate(temp_xx(i_nn(npx),j_nn(npy)),temp_yy(i_nn(npx),j_nn(npy)), &
                             temp_u(i_nn(npx),j_nn(npy)))

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_xx(ii,jj)=xx(ii,jj)
                    temp_yy(ii,jj)=yy(ii,jj)
                enddo
            enddo

            do jj=1,j_nn(npy)
                do ii=1,i_nn(npx)
                    temp_u(ii,jj)=u(ii,jj)
                enddo
            enddo

            call MPI_SEND(temp_xx,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,10,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_yy,i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,20,MPI_COMM_WORLD,ierr)
            call MPI_SEND(temp_u, i_nn(npx)*j_nn(npy),MPI_DOUBLE_PRECISION,0,40,MPI_COMM_WORLD,ierr)

            deallocate(temp_xx,temp_yy,temp_u)

        endif  !if (my_id .eq. 0 )
        deallocate(xx_global, yy_global)
        deallocate(u_global)

    end subroutine write_real_data_whole

end module write_data
