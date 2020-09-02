module read_input
  use def_variables
contains
  subroutine open_input()
    character(24) :: fName
    integer :: fStat
    call get_command_argument(1,fName, status=fStat)
    if (fStat /= 0) then
      print*,'Failed at reading input file. Exitting program...'
      call exit()
    end if
    open(unit=un_input,file=trim(fName), status='old')
  end subroutine open_input

  subroutine read_parameters()
    character(100) :: original_folder
    read(un_input,*) folder
    read(un_input,*) stat_output_name
    read(un_input,*) mode_int
    read(un_input,*) Scan_noise
    read(un_input,*) nu
    read(un_input,*) alpha
    read(un_input,*) mean_flight
    read(un_input,*) burst_amplitude
    read(un_input,*) N_measures
    read(un_input,*) N_steps
    read(un_input,*) N_reset
    read(un_input,*) N_iterations
    read(un_input,*) Nnbr
    read(un_input,*) pols_active
    read(un_input,*) nbrs_active
    read(un_input,*) write_config_active
    ! read(un_input,*) msd_active
    read(un_input,*) ignore_last_config
    read(un_input,*) int_L
    read(un_input,*) rho
    read(un_input,*) seed
    close(un_input)

    if ((trim(stat_output_name).eq.'None').or.&
        (trim(stat_output_name).eq.'none').or.&
        (trim(stat_output_name).eq.'NONE'))&
        stat_output_name=trim(folder)//'_'//hour(1:6)//'.stat'

    rho=rho*pi/9.d0
    mean_flight=mean_flight*3.d0/sqrt(pi)
    burst_amplitude=burst_amplitude*3.d0/sqrt(pi)

    L=dble(int_L)
    N=nint(L*L*rho)
    N_cells=int_L*int_L
    L2=0.5d0*L
    print*, 'N=',N


    original_folder=trim(folder)
    call read_levy_parameters()

    if (mode_int==0) mode_str='VicsekM' ! Classic Vicsek Model
    if (mode_int==1) mode_str='_LevyM_' ! Vicsek Model with Lévy behaviour
    if (mode_int==2) mode_str='bur-coa' ! burst-and-coast model

    folder=trim(folder)//'/'//mode_str
    print*, 'MODE: ',mode_str


    if (pols_active) call system('mkdir -p '//trim(folder)//'/pol')
    if (nbrs_active) call system('mkdir -p '//trim(folder)//'/nbr')
    ! if (msd_active)  call system('mkdir -p '//trim(folder)//'/msd')
    if (write_config_active) call system('mkdir -p '//trim(folder)//'/config')
    allocate(pos(N,2))
    allocate(abs_pos(N,2))
    allocate(unitary_vel(N,2))
    allocate(vel(N))
    allocate(header(N_cells))
    allocate(cell_list(N))
    allocate(list_nearest_nbr_cells(N_cells,9))
    allocate(speed_hist(600))
  end subroutine read_parameters

  subroutine read_levy_parameters()
    character(2) :: alpha_str
    write(alpha_str,'(i2)') nint(10*alpha)
    open(123, file='levy_parameters/A'//alpha_str//'.dat', status='old')
    read(123,*)
    read(123,*)
    read(123,*) alpha, A, B
    EXPONENT=1.d0/(1.d0-alpha)
    close(123)
  end subroutine read_levy_parameters

end module read_input
