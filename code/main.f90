program levy_program
  use def_variables ;  use init_module
  use various       ;  use read_input
  use integrators
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

    if (mode_int==0) then ! SIMPLE VICSEK MODEL MODE

      do i2=1,size(eta_array)
        eta=eta_array(i2)

        call read_last_configuration()
        call open_output_files()
        print*, 'alehop'

        do i1=1,N_measures
          call integrate_simple_vicsek(N_steps)
          call measurements_and_prints()
        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)

    else if (mode_int==1) then  ! LÉVY FLIGHTS MODE

      do i2=1,size(eta_array)
        eta=eta_array(i2)

        call read_last_configuration()
        call open_output_files()

        do i1=1,N_measures
          call integrate_levy_behaviour(N_steps)

          call measurements_and_prints()
        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)

    else if (mode_int==1) then   ! BURST-AND-COAST MODE

      do i2=1,size(eta_array)
        eta=eta_array(i2)

        call read_last_configuration()
        call open_output_files()

        do i1=1,N_measures
          call integrate_burstandcoast(N_steps)

          call measurements_and_prints()

        end do
        call close_output_files()
        call write_last_configuration()
      end do
      call print_stat(iteration,start_iteration)
    end if
  end do
  call open_output_files()

  call cpu_time(finish)
  open(unit=stat_output_unit, file=trim(stat_output_name),status='OLD', action='WRITE', position='APPEND')
  write(stat_output_unit,*) 'Total CPU time:',finish-start,'s'
  print*, 'Total CPU time:',finish-start,'s'
  close(stat_output_unit, status='DELETE')
end program levy_program
