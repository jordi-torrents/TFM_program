module def_variables
  implicit none
  real(4) :: start, finish
  real(8) :: pi=dacos(-1.d0), rho, L, L2, eta, dt, alpha, A, B, EXPONENT, mean_flight, timer, burst_amplitude
  real(8), allocatable :: pos(:,:), unitary_vel(:,:), eta_array(:), weight(:), vel(:)
  integer :: seed, un_input=56, N, N_steps, N_measures, N_iterations, N_reset, Nnbr, N_cells, int_L, int_time, mode_int
  integer :: pol_output_unit=0, speed_output_unit=1, nbr_output_unit=2, stat_output_unit=3, msd_output_unit=4, configuration_counter
  integer, allocatable :: header(:), cell_list(:), list_nearest_nbr_cells(:,:), speed_hist(:)
  character(100) :: folder, stat_output_name; character(7) :: mode_str; character(8) :: date ; character(10) :: hour
  logical :: print_polarization, print_GNF, scan_noise, print_configuration, ignore_last_config, print_speed
contains
end module
