module init_module
  use def_variables
  use integrators
  use read_input
  use various
contains
  subroutine initialize()
    integer :: i, j, io, nlines
    real(8) :: dr, r_nbr(Nnbr)

    time_counter=0
    if (Scan_noise) then
      nlines = 0
      open(153, file = 'nus.dat', status='OLD', action='READ')
      do
        read(153,*,iostat=io)
        if (io/=0) exit
        nlines = nlines + 1
      end do
      close(153)
      allocate(nu_array(nlines))
      open(153, file = 'nus.dat', status='OLD', action='READ')
      do i=1,size(nu_array)
        read(153,*) nu_array(i)
      end do
      close(153)
    else
      allocate(nu_array(1))
      nu_array(1)=nu
    end if

    open(unit=stat_output_unit, file=trim(stat_output_name))
    write(stat_output_unit,*) 'Date : '//date(7:8)//'-'//date(5:6)//'-'//date(1:4)&
    //' Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6), ' | Mode: ',mode_str
    close(stat_output_unit)

    do j=1,size(nu_array)
      nu=nu_array(j)
      call open_output_files()

      if (pols_active) then
        write(pol_output_unit,'(A34)') '# Date : '//date(7:8)//'-'//date(5:6)//'-'//date(1:4)&
        //' Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)
        write(pol_output_unit,'(A5,i7,A6,F5.3,A10,I3,A10,I6,A7,A20)') '# N =',N,' nu = ',nu,&
        ' N_steps =',N_steps,' N_reset =',N_reset,' MODE =',mode_str
      end if

      ! if (msd_active) then
      !   write(msd_output_unit,'(A34)') '# Date : '//date(7:8)//'-'//date(5:6)//'-'//date(1:4)&
      !   //' Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)
      ! end if

      if (nbrs_active) then
        write(nbr_output_unit,'(A34)') '# Date : '//date(7:8)//'-'//date(5:6)//'-'//date(1:4)&
        //' Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)
        write(nbr_output_unit,'(A5,i7,A6,F5.3,A10,I3,A10,I6,A7,A20)') '# N =',N,' nu = ',nu,&
        ' N_steps =',N_steps,' N_reset =',N_reset,' MODE =',mode_str
        dr=(L2)/dble(Nnbr)
        do i=1,Nnbr
          r_nbr(i)=i*dr
        end do
        write(nbr_output_unit,'(A1,f8.2,100f9.3)') '#',r_nbr
      end if
      call close_output_files()
    end do
  end subroutine

  subroutine set_geometry()
    integer :: cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell
    integer :: nbr_X(9), nbr_Y(9), k
    nbr_X=(/-1,0,1,-1,0,1,-1, 0, 1/)
    nbr_Y=(/ 1,1,1, 0,0,0,-1,-1,-1/)

    do cell_X=0,int_L-1
      do cell_Y=0,int_L-1
        cell = cell_X + int_L*cell_Y + 1
        do k=1,9
          nbr_cell_X=cell_X+nbr_X(k)
          nbr_cell_Y=cell_Y+nbr_Y(k)

          if (nbr_cell_X < 0) then
            nbr_cell = nbr_cell_X + int_L
          else if (nbr_cell_X >= int_L) then
            nbr_cell = nbr_cell_X - int_L
          else
            nbr_cell = nbr_cell_X
          end if

          if (nbr_cell_Y < 0) then
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y + int_L) + 1
          else if (nbr_cell_Y >= int_L) then
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y - int_L) + 1
          else
            nbr_cell = nbr_cell + int_L*nbr_cell_Y + 1
          end if

          list_nearest_nbr_cells(cell,k) = nbr_cell
        end do
      end do
    end do
  end subroutine

  subroutine reset_system()
    real :: r1279
    real(8) :: theta_i
    configuration_counter=0
    int_time=0
    timer=0.d0
    vel=0.d0
    do i=1,N
      pos(i,:) = (/r1279(), r1279()/)*L
      theta_i = r1279()*2.0*pi
      unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
    end do
  end subroutine

  subroutine prepare_systems()
    integer :: i
    logical :: file_exists
    character(2) :: nu_str
    call system('mkdir -p '//trim(folder)//'/last_state')
    ! print*, 'Preparing system'

    do i=1,size(nu_array)
      nu=nu_array(i)
      write(nu_str,'(i2.2)') nint(nu*100.d0)
      inquire(file=trim(trim(folder)//'/last_state/nu'//nu_str//'.csv'), exist=file_exists)

      if ((.not.file_exists).or.(ignore_last_config)) then
        call reset_system()
        nu=0.0
        call integrate_simple_vicsek(N_reset)
        nu=nu_array(i)
        if (mode_int==0) call integrate_simple_vicsek(N_reset)
        if (mode_int==1) call integrate_levy_behaviour(N_reset)
        if (mode_int==2) call integrate_burstandcoast(N_reset)
        int_t=0


        call write_last_configuration()
        open(unit=stat_output_unit, file=trim(stat_output_name), status='OLD', action='WRITE', position='APPEND')
        write(stat_output_unit,'(A24,f6.3)') 'System prepared with nu=',nu
        close(stat_output_unit)
      end if
    end do
  end subroutine
end module
