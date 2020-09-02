module various
  use def_variables
contains
  subroutine print_stat(i, start_iteration)
    real :: start_iteration, finish_iteration

    call date_and_time(date,hour)
    call cpu_time(finish_iteration)
    open(unit=stat_output_unit, file=trim(stat_output_name), status='OLD', action='WRITE', position='APPEND')
    write(stat_output_unit,'(A10,I4,A4,I4,A22,f7.2,A17)') 'Iteration ',i,' of ',N_iterations,&
            ' | Iteration time (s):', finish_iteration-start_iteration, ' | Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)
    close(stat_output_unit)
  end subroutine

  subroutine open_output_files()
    character(2) :: nu_str
    character(50):: pol_output_name, nbr_output_name!, msd_output_name
    write(nu_str,'(i2.2)') nint(nu*100.d0)
    pol_output_name = trim(folder)//'/pol/nu'//nu_str//'.csv'
    nbr_output_name = trim(folder)//'/nbr/nu'//nu_str//'.csv'
    ! msd_output_name = trim(folder)//'/msd/nu'//nu_str//'.csv'
    if (pols_active) open(unit=pol_output_unit, file=trim(pol_output_name), position='APPEND')
    if (nbrs_active) open(unit=nbr_output_unit, file=trim(nbr_output_name), position='APPEND')
    ! if (msd_active)  open(unit=msd_output_unit, file=trim(msd_output_name), position='APPEND')
  end subroutine

  subroutine close_output_files()
    if (pols_active) close(pol_output_unit)
    if (nbrs_active) close(nbr_output_unit)
    if (pols_active) close(msd_output_unit)
  end subroutine

  subroutine read_last_configuration()
    character(2) :: nu_str
    write(nu_str,'(i2.2)') nint(nu*100.d0)
    open(unit=542, file=trim(trim(folder)//'/last_state/nu'//nu_str//'.csv'))
    do i=1,N
      read(542,*) pos(i,:), unitary_vel(i,:), vel(i), abs_pos(i,:)
    end do
      read(542,*) t
    close(542)
  end subroutine

  subroutine write_last_configuration()
    character(2) :: nu_str
    write(nu_str,'(i2.2)') nint(nu*100.d0)
    open(unit=542, file=trim(trim(folder)//'/last_state/nu'//nu_str//'.csv'))
    do i=1,N
      write(542,*) pos(i,:), unitary_vel(i,:), vel(i), abs_pos(i,:)
    end do
      write(542,*) t
    close(542)
  end subroutine
end module
