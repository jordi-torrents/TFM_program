module interaction
  use def_variables
contains
  subroutine build_cell_lists()
    integer :: particle, cell
    header = 0
    do particle=1,N
      cell = int(pos(particle,1)) + int_L*int(pos(particle,2)) + 1
      cell_list(particle) = header(cell)
      header(cell) = particle
    end do
  end subroutine

  subroutine all_neighbours_direction(nbr_direction)
    real(8) :: nbr_vel(N,2), nbr_direction(N)
    integer :: cell, i, j, k
    call build_cell_lists()
    nbr_vel=unitary_vel

    do cell=1,N_cells
      i = header(cell)
      do while (i>0)
        j = header(cell)
        do while (j>0)
          if (i < j) then
            nbr_vel(i,:) = nbr_vel(i,:) + unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + unitary_vel(i,:)
          end if
          j = cell_list(j)
        end do
        i = cell_list(i)
      end do

      do k=1,4
        i = header(cell)
        do while (i>0)
          j = header(list_nearest_nbr_cells(cell,k))
          do while (j>0)
            nbr_vel(i,:) = nbr_vel(i,:) + unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + unitary_vel(i,:)
            j = cell_list(j)
          end do
          i = cell_list(i)
        end do
      end do
    end do

    nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
  end subroutine


  subroutine all_neighbours_direction2()
    real(8) :: nbr_vel(N), nbr_dir(N,2), theta_i
    integer :: cell, i, j, k, counter(N)
    call build_cell_lists()
    nbr_dir=unitary_vel
    nbr_vel=vel
    counter=1

    do cell=1,N_cells
      i = header(cell)
      do while (i>0)
        j = header(cell)
        do while (j>0)
          if (i < j) then
            nbr_vel(i) = nbr_vel(i) + vel(j)
            nbr_vel(j) = nbr_vel(j) + vel(i)
            nbr_dir(i,:) = nbr_dir(i,:) + unitary_vel(j,:)
            nbr_dir(j,:) = nbr_dir(j,:) + unitary_vel(i,:)
            counter(i)=counter(i)+1
            counter(j)=counter(j)+1
          end if
          j = cell_list(j)
        end do
        i = cell_list(i)
      end do

      do k=1,4
        i = header(cell)
        do while (i>0)
          j = header(list_nearest_nbr_cells(cell,k))
          do while (j>0)
            nbr_vel(i) = nbr_vel(i) + vel(j)
            nbr_vel(j) = nbr_vel(j) + vel(i)
            nbr_dir(i,:) = nbr_dir(i,:) + unitary_vel(j,:)
            nbr_dir(j,:) = nbr_dir(j,:) + unitary_vel(i,:)
            counter(i)=counter(i)+1
            counter(j)=counter(j)+1
            j = cell_list(j)
          end do
          i = cell_list(i)
        end do
      end do
    end do
    do i=1,N
      theta_i=atan2(nbr_dir(i,2),nbr_dir(i,1))+nu*pi*2.d0*(r1279()-0.5)
      unitary_vel(i,:)=(/cos(theta_i),sin(theta_i)/)
      vel(i)=max(vel(i),nbr_vel(i)/dble(counter(i)))
    end do
  end subroutine

end module
