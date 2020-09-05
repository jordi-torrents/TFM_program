module various
  use def_variables
contains

  subroutine measurements_and_prints()
    integer :: i, nbr(Nnbr), indx_nbr
    real(8) :: r, dr, dxy(2), mean
    character(4) :: time_str ; character(2) :: eta_str

    if (print_polarization) write(pol_output_unit,'(f17.15)') sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N)

    if (print_speed) then
        mean=sum(vel)/dble(N)
        write(speed_output_unit,*) mean*sqrt(pi)/3.d0, sqrt(sum((vel-mean)**2)/dble(N-1))*sqrt(pi)/3.d0
    end if

    if (print_GNF) then
      nbr=0
      dr=(L2)/dble(Nnbr)
      do i=1,N
          dxy = pos(i,:)-(/L2,L2/)
          r=sqrt(sum(dxy**2))
          indx_nbr = int(r/dr) + 1
          if (indx_nbr<=Nnbr) nbr(indx_nbr) = nbr(indx_nbr) + 1
      end do
      do i=Nnbr,1,-1
        nbr(i)=sum(nbr(1:i))
      end do
      write(nbr_output_unit,'(1000I6)') nbr(:)
    end if

    if (print_configuration) then
      write(time_str,'(i4.4)') configuration_counter
      write(eta_str,'(i2.2)') nint(eta*100.d0)
      open(unit=823,file=trim(trim(folder)//'/configuration/eta'//eta_str//'/'//time_str//'.csv'))
      do i=1,N
        write(823,'(2f8.3)') pos(i,:)
      end do
      close(823)
      configuration_counter = configuration_counter+1
    end if
  end subroutine

  subroutine print_stat(i, start_iteration)
    real :: start_iteration, finish_iteration
    call date_and_time(date,hour)
    call cpu_time(finish_iteration)
    open(unit=stat_output_unit, file=trim(stat_output_name), status='OLD', action='WRITE', position='APPEND')
    write(stat_output_unit,'(A10,I4,A4,I4,A22,f7.2,A17)') 'Iteration ',i,' of ',N_iterations,&
            ' | Iteration time (s):', finish_iteration-start_iteration, ' | Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)

    print'(A10,I4,A4,I4,A22,f7.2,A17)', 'Iteration ',i,' of ',N_iterations,&
            ' | Iteration time (s):', finish_iteration-start_iteration, ' | Hour: '//hour(1:2)//':'//hour(3:4)//':'//hour(5:6)
    close(stat_output_unit)
  end subroutine

  subroutine open_output_files()
    character(2) :: eta_str
    character(50):: pol_output_name, nbr_output_name, speed_output_name
    write(eta_str,'(i2.2)') nint(eta*100.d0)
    pol_output_name = trim(folder)//'/polarization/eta'//eta_str//'.csv'
    nbr_output_name = trim(folder)//'/GNF/eta'//eta_str//'.csv'
    speed_output_name = trim(folder)//'/speed/eta'//eta_str//'.csv'
    if (print_polarization) open(unit=pol_output_unit, file=trim(pol_output_name), position='APPEND')
    if (print_GNF) open(unit=nbr_output_unit, file=trim(nbr_output_name), position='APPEND')
    if (print_speed)  open(unit=speed_output_unit, file=trim(speed_output_name), position='APPEND')
  end subroutine

  subroutine close_output_files()
    if (print_polarization) close(pol_output_unit)
    if (print_GNF) close(nbr_output_unit)
    if (print_speed) close(speed_output_unit)
  end subroutine

  subroutine read_last_configuration()
    character(2) :: eta_str
    write(eta_str,'(i2.2)') nint(eta*100.d0)
    open(unit=542, file=trim(trim(folder)//'/last_state/eta'//eta_str//'.csv'))
    do i=1,N
      read(542,*) pos(i,:), unitary_vel(i,:), vel(i)
    end do
    read(542,*) int_time, timer
    close(542)
  end subroutine

  subroutine write_last_configuration()
    character(2) :: eta_str
    write(eta_str,'(i2.2)') nint(eta*100.d0)
    open(unit=542, file=trim(trim(folder)//'/last_state/eta'//eta_str//'.csv'))
    do i=1,N
      write(542,*) pos(i,:), unitary_vel(i,:), vel(i)
    end do
    write(542,*) int_time, timer
    close(542)
  end subroutine
end module
