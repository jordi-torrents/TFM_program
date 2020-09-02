program levy_program
  use def_variables ;  use init_module
  use various       ;  use read_input
  use integrators   ;  use measurements
  ! integer :: i
  call cpu_time(start)
  call date_and_time(date,hour)

  call open_input()
  call read_parameters()
  call setr1279(seed)
  call set_geometry()

  call initialize()
  call prepare_systems()
  do  iteration=1,N_iterations
    call cpu_time(start_iteration)

    if (mode_int==0) then

      do i2=1,size(nu_array)
        nu=nu_array(i2)

        call read_last_configuration()
        call open_output_files()

        do i1=1,N_measures
          call integrate_simple_vicsek(N_steps)
          ! call Guillespie(N_steps)

          if (pols_active) write(pol_output_unit,'(f17.15)') sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N)
          if (nbrs_active) call gdr_nbr_step()
          ! if (msd_active) write(msd_output_unit,'(f9.1,E15.6)') t, sum(abs_pos**2)/dble(N)
          if (write_config_active) call write_configuration()
        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)
    else if (mode_int==1) then
      do i2=1,size(nu_array)
        nu=nu_array(i2)

        call read_last_configuration()
        call open_output_files()

        do i1=1,N_measures
          call integrate_levy_behaviour(N_steps)

          if (pols_active) write(pol_output_unit,'(f17.15)') sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N)
          if (nbrs_active) call gdr_nbr_step()
          if (write_config_active) call write_configuration()
        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)
    else if (mode_int==1) then
      do i2=1,size(nu_array)
        nu=nu_array(i2)

        call read_last_configuration()
        call open_output_files()

        do i1=1,N_measures
          call integrate_burstandcoast(N_steps)

          if (pols_active) write(pol_output_unit,'(f17.15)') sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N)
          if (nbrs_active) call gdr_nbr_step()
          if (write_config_active) call write_configuration()
        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)
    end if
  end do
  call open_output_files()
  ! do i=1,600
  !   write(nbr_output_unit,*) -0.01+real(i)*0.02, speed_hist(i)
  ! end do


      ! call reset_system()
      ! call initialize()
      ! call Guillespie(N_reset)
      ! open(unit=555,file=trim(folder)//'.dat')
      ! t=0.d0
      ! do i1=1,N_measures
      !   call Guillespie(1)
      !   mean=sum(vel)/dble(N)
      !   write(555,*) sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N),&
      !   mean*sqrt(pi)/3.d0!,sqrt(sum((vel-mean)**2)/dble(N-1))*sqrt(pi)/3.d0, &
      !   ! vel(1)*sqrt(pi)/3.d0,vel(2)*sqrt(pi)/3.d0,vel(3)*sqrt(pi)/3.d0
      !   ! call write_configuration()
      ! end do


  call cpu_time(finish)
  write(*,*) 'Total CPU time:',finish-start,'s'
  open(unit=stat_output_unit, file=trim(stat_output_name),status='OLD', action='WRITE', position='APPEND')
  write(stat_output_unit,*) 'Total CPU time:',finish-start,'s'
  close(stat_output_unit, status='DELETE')
end program levy_program
