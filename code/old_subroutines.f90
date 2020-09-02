subroutine vicsek_steps2(isteps)
  integer :: cell, particle, i
  do j=1,isteps

    call neighbours_direction()
    do i=1,N
      randoms(i) = r1279()
    end do
    theta(:) = nbr_direction(:)+noise_cnst*(randoms(:)-0.5)
    vel(:,1) = dcos(theta(:))
    vel(:,2) = dsin(theta(:))

    do i=1,number_inner_cells
      cell=inner_cells(i)
      particle=header(cell)
      do while (particle>0)
        pos(particle,:) = pos(particle,:) + dt*vel(particle,:)
        particle = cell_list(particle)
      end do
    end do

    do i=1,number_outer_cells
      cell=outer_cells(i)
      particle=header(cell)
      do while (particle>0)
        pos(particle,:)=modulo(pos(particle,:) + dt*vel(particle,:),L)
        particle = cell_list(particle)
      end do
    end do

  end do
end subroutine

subroutine neighbours_direction2()
  real(8) :: dxy(2), nbr_vel(N,2)
  integer :: cell, i, j, k, l
  call build_cell_lists()
  nbr_vel=vel

  do l=1,number_inner_cells
    cell = inner_cells(l)

    i = header(cell)
    do while (i>0)
      j = header(cell)
      do while (j>0)
        if (i < j) then
          if (sum((pos(i,:)-pos(j,:))**2)<threshold) then
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
          end if
        end if
        j = cell_list(j)
      end do
      i = cell_list(i)
    end do

    do k=1,4
      i = header(cell)
      do while (i>0)
        j = header(nbr_cell_list(cell,k))
        do while (j>0)
          if (sum((pos(i,:)-pos(j,:))**2)<threshold) then
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
          end if
          j = cell_list(j)
        end do
        i = cell_list(i)
      end do
    end do
  end do

  do l=1,number_outer_cells
    cell = outer_cells(l)

    i = header(cell)
    do while (i>0)
      j = header(cell)
      do while (j>0)
        if (i < j) then
          if (sum((pos(i,:)-pos(j,:))**2)<threshold) then
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
          end if
        end if
        j = cell_list(j)
      end do
      i = cell_list(i)
    end do

    do k=1,4
      i = header(cell)
      do while (i>0)
        j = header(nbr_cell_list(cell,k))
        do while (j>0)
          dxy = (/PBC_dist(pos(i,1)-pos(j,1)),PBC_dist(pos(i,2)-pos(j,2))/)
          if (sum(dxy**2)<threshold) then
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
          end if
          j = cell_list(j)
        end do
        i = cell_list(i)
      end do
    end do
  end do

  nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
end subroutine

subroutine neighbours_direction3()
  real(8) :: dxy(2), nbr_vel(N,2), shift(2)
  integer :: cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, cell, nbr_cell, i, j
  call build_cell_lists()
  nbr_vel=vel

  do cell_X=0,int_L-1
    do cell_Y=0,int_L-1
      cell = cell_X + int_L*cell_Y + 1
      do nbr_cell_X=cell_X-1,cell_X+1
        do nbr_cell_Y=cell_Y-1,cell_Y+1

          if (nbr_cell_X < 0) then
            shift(1) = -L
            nbr_cell = nbr_cell_X + int_L
          else if (nbr_cell_X >= int_L) then
            shift(1) =  L
            nbr_cell = nbr_cell_X - int_L
          else
            shift(1) = 0.d0
            nbr_cell = nbr_cell_X
          end if

          if (nbr_cell_Y < 0) then
            shift(2) = -L
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y + int_L) + 1
          else if (nbr_cell_Y >= int_L) then
            shift(2) =  L
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y - int_L) + 1
          else
            shift(2) = 0.d0
            nbr_cell = nbr_cell + int_L*nbr_cell_Y + 1
          end if

          i = header(cell)
          do while (i.ne.0)
            j = header(nbr_cell)
            do while (j.ne.0)
              if (i < j) then
                dxy = pos(i,:)-pos(j,:)-shift(:)
                if (sum(dxy**2)<threshold) then
                  nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
                  nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
                end if
              end if
              j = cell_list(j)
            end do
            i = cell_list(i)
          end do
        end do
      end do
    end do
  end do
  nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
end subroutine

subroutine neighbours_direction4()
  real(8) :: dxy(2), nbr_vel(N,2)
  ! call build_cell_lists()
  nbr_vel=vel
  do i=1,n-1
    do j=i+1,n
      dxy = (/PBC_dist(pos(i,1)-pos(j,1)),&
             &PBC_dist(pos(i,2)-pos(j,2))/)
      if (sum(dxy**2)<threshold) then
        nbr_vel(i,:) = nbr_vel(i,:) + vel(j,:)
        nbr_vel(j,:) = nbr_vel(j,:) + vel(i,:)
      end if
    end do
  end do
  nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
end subroutine

function PBC_dist(x)
  real(8) :: PBC_dist, x
  PBC_dist = x - int(x/L2)*L
end function PBC_dist


  subroutine all_neighbours_direction(nbr_direction)
    real(8) :: nbr_vel(N,2), nbr_pos(N,2), nbr_direction(N), dist(2)
    integer :: cell, i, j, k, counter(N), integer0
    call build_cell_lists()
    nbr_vel=unitary_vel
    nbr_pos=0.d0
    counter=0

    do integer0=1,number_inner_cells
      cell=inner_cells(integer0)
      i = header(cell)
      do while (i>0)
        j = header(cell)
        do while (j>0)
          if (i < j) then
            nbr_pos(i,:) = nbr_pos(i,:) + pos(j,:)+unitary_vel(j,:)
            nbr_pos(j,:) = nbr_pos(j,:) + pos(i,:)+unitary_vel(i,:)
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
            nbr_pos(i,:) = nbr_pos(i,:) + pos(j,:)+unitary_vel(j,:)
            nbr_pos(j,:) = nbr_pos(j,:) + pos(i,:)+unitary_vel(i,:)
            counter(i)=counter(i)+1
            counter(j)=counter(j)+1
            j = cell_list(j)
          end do
          i = cell_list(i)
        end do
      end do
    end do

    do integer0=1,number_outer_cells
      cell=outer_cells(integer0)
      i = header(cell)
      do while (i>0)
        j = header(cell)
        do while (j>0)
          if (i < j) then
            nbr_pos(i,:) = nbr_pos(i,:) + pos(j,:)+unitary_vel(j,:)
            nbr_pos(j,:) = nbr_pos(j,:) + pos(i,:)+unitary_vel(i,:)
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
            dist=pos(i,:)-pos(j,:)
            ! if (dist(1)>3.d0) dist(1)=dist(1)-L
            ! if (dist(2)>3.d0) dist(2)=dist(2)-L
            ! if (dist(1)<-3.d0) dist(1)=dist(1)+L
            ! if (dist(2)<-3.d0) dist(2)=dist(2)+L
            ! dist(1)=dist(1)-dble(int(2.d0*dist(1)/L))*L
            dist=dist-dble(int(2.d0*dist/L))*L
            ! if (any(abs(dist)>5.d0)) print*, pos(i,:),pos(j,:),(pos(i,:)-pos(j,:)),dist
            nbr_pos(i,:) = nbr_pos(i,:) + pos(i,:)-dist+unitary_vel(j,:)
            nbr_pos(j,:) = nbr_pos(j,:) + pos(j,:)+dist+unitary_vel(i,:)
            counter(i)=counter(i)+1
            counter(j)=counter(j)+1
            j = cell_list(j)
          end do
          i = cell_list(i)
        end do
      end do
    end do
    do i=1,N
      if (counter(i)>0) then
        nbr_pos(i,:) = nbr_pos(i,:)/dble(counter(i))-pos(i,:)
      else
        nbr_pos(i,:) = unitary_vel(i,:)
      end if
    end do
    nbr_direction=atan2(nbr_pos(:,2),nbr_pos(:,1))
  end subroutine

  subroutine weighted_all_neighbours_direction(nbr_direction, nbr_parti_vel)
    real(8) :: nbr_vel(N,2), nbr_direction(N), nbr_parti_vel(N)
    integer :: cell, i, j, k, counter(N)
    call build_cell_lists()
    do i=1,N
      nbr_vel(i,:)=vel(i)*unitary_vel(i,:)
      nbr_parti_vel(i)=vel(i)
    end do
    counter=1
    do cell=1,N_cells
      i = header(cell)
      do while (i>0)
        j = header(cell)
        do while (j>0)
          if (i < j) then
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j)*unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i)*unitary_vel(i,:)
            nbr_parti_vel(i)=nbr_parti_vel(i)+vel(j)
            nbr_parti_vel(j)=nbr_parti_vel(j)+vel(i)
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
            nbr_vel(i,:) = nbr_vel(i,:) + vel(j)*unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + vel(i)*unitary_vel(i,:)
            nbr_parti_vel(i)=nbr_parti_vel(i)+vel(j)
            nbr_parti_vel(j)=nbr_parti_vel(j)+vel(i)
            counter(i)=counter(i)+1
            counter(j)=counter(j)+1
            j = cell_list(j)
          end do
          i = cell_list(i)
        end do
      end do
    end do
    nbr_parti_vel=nbr_parti_vel/dble(counter)
    nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
  end subroutine

  subroutine selective_neighbours_direction(to_be_updated, nbr_direction)
    logical :: to_be_updated(N)
    real(8) :: nbr_vel(N,2), nbr_direction(N)
    integer :: cell, i, j, k
    call build_cell_lists()

    do i=1,N
      if (to_be_updated(i)) then
        nbr_vel(i,:)=0.d0
        cell = where_is(i)
        do k=1,9
          j = header(list_nearest_nbr_cells(cell,k))
          do while (j>0)
            nbr_vel(i,:) = nbr_vel(i,:) + unitary_vel(j,:)
            j = cell_list(j)
          end do
        end do
        nbr_direction(i)=atan2(nbr_vel(i,2),nbr_vel(i,1))
      end if
    end do
  end subroutine

  subroutine weighted_selective_neighbours_direction(to_be_updated, nbr_direction)
    logical :: to_be_updated(N)
    real(8) :: nbr_vel(N,2), nbr_direction(N)
    integer :: cell, i, j, k
    call build_cell_lists()

    do i=1,N
      if (to_be_updated(i)) then
        nbr_vel(i,:)=0.d0
        cell = where_is(i)
        do k=1,9
          j = header(list_nearest_nbr_cells(cell,k))
          do while (j>0)
            nbr_vel(i,:) = nbr_vel(i,:) + weight(j)*unitary_vel(j,:)
            j = cell_list(j)
          end do
        end do
        nbr_direction(i)=atan2(nbr_vel(i,2),nbr_vel(i,1))
      end if
    end do
  end subroutine


  subroutine one_neighbours_direction(to_be_updated, nbr_direction, nbr_velocity)
    real(8) :: nbr_vel(2), nbr_vel2, nbr_direction, nbr_velocity
    integer :: cell, cell_j, j, k, to_be_updated, counter
    counter=0
    nbr_vel=0.d0
    nbr_vel2=0.d0
    cell = modulo(floor(pos(to_be_updated,1)),int_L) + int_L*modulo(floor(pos(to_be_updated,2)),int_L) + 1

    do k=1,25
      j = header(list_2nd_nbr_cells(cell,k))
      do while (j>0)
        cell_j=modulo(floor(pos(j,1)),int_L) + int_L*modulo(floor(pos(j,2)),int_L) + 1
        if (ANY( list_nearest_nbr_cells(cell,:)==cell_j )) then
          nbr_vel(:) = nbr_vel(:) + unitary_vel(j,:)
          nbr_vel2 = nbr_vel2 + vel(j)
          counter=counter+1
        end if
        j = cell_list(j)
      end do
    end do
    nbr_velocity=nbr_vel2/dble(counter)
    nbr_direction=atan2(nbr_vel(2),nbr_vel(1))
  end subroutine

  subroutine build_cell_lists()
    integer particle, cell
    header = 0
    do particle=1,N
      cell = modulo(floor(pos(particle,1)),int_L) + int_L*modulo(floor(pos(particle,2)),int_L) + 1
      where_is(particle) = cell
      cell_list(particle) = header(cell)
      header(cell) = particle
    end do
  end subroutine


    subroutine asinc_vicsek(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt

        to_be_updated=timekeeper<0
        call selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,1) = cos(theta_i)
            unitary_vel(i,2) = sin(theta_i)
            timekeeper(i) = timekeeper(i)+1.d0
          end if
        end do

        pos(:,1) = modulo(pos(:,1) + dt*mean_flight*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*mean_flight*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_V_levy(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt

        to_be_updated=timekeeper<0
        call selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,1) = cos(theta_i)
            unitary_vel(i,2) = sin(theta_i)
            vel(i) = mean_flight*(A*r1279()+B)**EXP
            timekeeper(i) = timekeeper(i)+1.d0
          end if
        end do
        pos(:,1) = modulo(pos(:,1) + dt*vel*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*vel*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine sinc_levy_fixed_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      noise_cnst = nu*pi*2.d0

      do j=1,isteps
        call weighted_all_neighbours_direction(nbr_direction)
        do i=1,N
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
        end do

        pos(:,1) = modulo(pos(:,1) + vel*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + vel*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine sinc_levy_changing_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      noise_cnst = nu*pi*2.d0

      do j=1,isteps
        call weighted_all_neighbours_direction(nbr_direction)
        do i=1,N
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
            weight(i) = (A*r1279()+B)**EXP
            vel(i)=mean_flight*weight(i)
        end do

        pos(:,1) = modulo(pos(:,1) + vel*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + vel*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_V_levy_fixed_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt
        to_be_updated=timekeeper<0
        call weighted_selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
            timekeeper(i) = timekeeper(i)+1.d0
          end if
        end do

        pos(:,1) = modulo(pos(:,1) + dt*vel*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*vel*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_T_levy_fixed_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt
        to_be_updated=timekeeper<0
        call weighted_selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
            timekeeper(i) = timekeeper(i)+weight(i)
          end if
        end do

        pos(:,1) = modulo(pos(:,1) + dt*mean_flight*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*mean_flight*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_V_levy_changing_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt
        to_be_updated=timekeeper<0
        call weighted_selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
            weight(i) = (A*r1279()+B)**EXP
            vel(i) = mean_flight*weight(i)
            timekeeper(i) = timekeeper(i)+1.d0
          end if
        end do

        pos(:,1) = modulo(pos(:,1) + dt*vel*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*vel*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_T_levy_changing_weight(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0

      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt
        to_be_updated=timekeeper<0
        call weighted_selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
            weight(i) = (A*r1279()+B)**EXP
            timekeeper(i) = timekeeper(i)+weight(i)
          end if
        end do

        pos(:,1) = modulo(pos(:,1) + dt*mean_flight*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + dt*mean_flight*unitary_vel(:,2),L)
      end do
    end subroutine

    subroutine asinc_T_levy(isteps)
      real(8) :: noise_cnst, theta_i, nbr_direction(N)
      logical :: to_be_updated(N)
      noise_cnst = nu*pi*2.d0
      pos = modulo(pos,L)
      do j=1,nint(isteps/dt)
        timekeeper=timekeeper-dt

        to_be_updated=timekeeper<0
        call selective_neighbours_direction(to_be_updated, nbr_direction)
        do i=1,N
          if (to_be_updated(i)) then
            theta_i = nbr_direction(i)+noise_cnst*(r1279()-0.5)
            unitary_vel(i,1) = cos(theta_i)
            unitary_vel(i,2) = sin(theta_i)
            timekeeper(i) = timekeeper(i)+(A*r1279()+B)**EXP
          end if
        end do
        pos(:,1) = pos(:,1) + dt*mean_flight*unitary_vel(:,1)
        pos(:,2) = pos(:,2) + dt*mean_flight*unitary_vel(:,2)
      end do
      t=t+dble(isteps)
    end subroutine

    subroutine Guillespie2(isteps)
      real :: r
      real(8) :: CDF(N), nbr_direction, old_lambda, delta_t, theta_j,noise_cnst, nbr_velocity, acc
      integer :: every
      acc=1.0d0
      every=N
      noise_cnst = nu*pi*2.0d0
      CDF(1)=lambda(1)
      do j=2,N
        CDF(j)=CDF(j-1)+lambda(j)
      end do
      do i=0,isteps-1
        if (mod(i,every)==0) then
          pos = modulo(pos,L)
          call build_cell_lists()
        end if
        r=r1279()*real(CDF(N))
        do j=1,N
          if (CDF(j)>r) exit
        end do

        call one_neighbours_direction(j, nbr_direction, nbr_velocity)
        ! vel(j) = mean_flight
        ! vel(j) = mean_flight*(r1279()+0.5)
        ! vel(j) = mean_flight*(A*r1279()+B)**EXP
        ! vel(j) = nbr_velocity
        vel(j)=vel(j)+0.3d0
        theta_j = nbr_direction+noise_cnst*(r1279()-0.5)
        unitary_vel(j,:) = (/cos(theta_j),sin(theta_j)/)
        tau(j)=max(vel(j)-nbr_velocity,0.1d0)

        old_lambda=lambda(j)
        lambda(j)=1.d0/tau(j)
        CDF(j:N)=CDF(j:N)-old_lambda+lambda(j)
        delta_t=-log(r1279()+1.e-10)/CDF(N)
        t=t+delta_t
        timekeeper=timekeeper+delta_t
        A=timekeeper(j)
        timekeeper(j)=0.d0


        pos(:,1) = pos(:,1) + delta_t*0.5d0*(vel+max(vel-delta_t*acc,0.d0))*unitary_vel(:,1)
        pos(:,2) = pos(:,2) + delta_t*0.5d0*(vel+max(vel-delta_t*acc,0.d0))*unitary_vel(:,2)
        vel=max(vel-delta_t*acc,0.d0)

      end do
    end subroutine

    subroutine Guillespie(isteps)
      real :: r
      real(8) :: CDF(N), nbr_direction, old_lambda, delta_t, theta_j,noise_cnst, nbr_velocity
      integer :: every
      every=10*N
      noise_cnst = nu*pi*2.d0
      CDF(1)=lambda(1)
      do j=2,N
        CDF(j)=CDF(j-1)+lambda(j)
      end do
      do i=0,isteps*N-1
        if (mod(i,every)==0) then
          pos = modulo(pos,L)
          call build_cell_lists()
        end if
        r=r1279()*real(CDF(N))
        do j=1,N
          if (CDF(j)>r) exit
        end do

        delta_t=-log(r1279()+1.e-10)/CDF(N)



        call one_neighbours_direction(j, nbr_direction, nbr_velocity)
        theta_j = nbr_direction+noise_cnst*(r1279()-0.5)
        unitary_vel(j,:) = (/cos(theta_j),sin(theta_j)/)



        pos(:,1) = pos(:,1) + delta_t*exp(-internal_t(:))*unitary_vel(:,1)
        pos(:,2) = pos(:,2) + delta_t*exp(-internal_t(:))*unitary_vel(:,2)
        internal_t=internal_t+delta_t
        t=t+delta_t

        if (internal_t(j)>1.d0) then
          internal_t(j)=0.d0
          ! call fire(j)
        end if

        tau(j)=0.75d0*mean_flight**2/(exp(-internal_t(j))**2)
        old_lambda=lambda(j)
        lambda(j)=1.d0/tau(j)
        CDF(j:N)=CDF(j:N)-old_lambda+lambda(j)

      end do
    end subroutine

    subroutine sinc_vicsek(isteps) ! dt=1.0
      real(8) :: noise_cnst, theta(N), nbr_direction(N), nbr_vel(N)
      noise_cnst = nu*pi*2.d0
      do j=1,isteps
        call weighted_all_neighbours_direction(nbr_direction, nbr_vel)
        do i=1,N
          theta(i) = nbr_direction(i)+noise_cnst*(r1279()-0.5)
          energy(i)=energy(i)+mean_flight*mean_flight-vel(i)*vel(i)

          vel(i)=

          ! if (nbr_vel(i)>sqrt(energy(i))) then
          !   vel(i)=sqrt(energy(i))
          ! else if (energy(i)>1.d0) then
          !   vel(i)=10.d0*mean_flight
          ! else
          !   vel(i)=nbr_vel(i)
          ! end if

          ! ratio=sqrt(energy(i))/nbr_vel(i)
          ! x_min=((2.d0-alpha)/(1.d0-alpha))*((ratio**(1.d0-alpha)-1.d0)/(ratio**(2.d0-alpha)-1.d0))
          ! A=(ratio*x_min)**(1.d0-alpha)-x_min**(1.d0-alpha)
          ! B=x_min**(1.d0-alpha)
          ! vel(i)=nbr_vel(i)*(A*r1279()+B)**EXP

          ! if (nbr_vel(i)>sqrt(energy(i))) then
          !   vel(i)=nbr_vel(i)*(1.d0+log(r1279()+1.e-10))
          !   else
          !   vel(i)=nbr_vel(i)-sqrt(energy(i))*log(r1279()+1.e-10)
          ! end if

          ! energy(i)=energy(i)+mean_flight*mean_flight-vel(i)*vel(i)
          ! if (energy(i)<0.d0) print*, t,'fuk',energy(i),vel(i)
          ! if (nbr_vel(i)>sqrt(energy(i))) then
          !   vel(i)=nbr_vel(i)-abs(sqrt(-2.d0*energy(i)*log(r1279()))*cos(2.d0*pi*r1279()))
          !   else
          !   vel(i)=nbr_vel(i)+abs(sqrt(-2.d0*energy(i)*log(r1279()))*cos(2.d0*pi*r1279()))
          ! end if
          ! energy(i)=energy(i)-vel(i)*vel(i)

        end do
        unitary_vel(:,1) = cos(theta)
        unitary_vel(:,2) = sin(theta)
        pos(:,1) = modulo(pos(:,1) + vel(:)*unitary_vel(:,1),L)
        pos(:,2) = modulo(pos(:,2) + vel(:)*unitary_vel(:,2),L)
        t=t+1.d0
      end do

    end subroutine
