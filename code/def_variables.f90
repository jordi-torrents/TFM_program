module def_variables
  implicit none
  real(4) :: start, finish
  real(8) :: pi=dacos(-1.d0), rho, L, L2, nu, dt, alpha, A, B, EXPONENT, mean_flight, timer, burst_amplitude
  real(8), allocatable :: pos(:,:),unitary_vel(:,:), nu_array(:),weight(:),vel(:), abs_pos(:,:)
  integer :: seed, un_input=56, N, N_steps, N_measures, N_iterations, N_reset, Nnbr, N_cells, int_L, int_time, mode_int
  integer :: pol_output_unit=0, nbr_output_unit=1, stat_output_unit=3, msd_output_unit=4, configuration_counter
  integer, allocatable :: header(:), cell_list(:), list_nearest_nbr_cells(:,:), speed_hist(:)
  character(100) :: folder, stat_output_name; character(7) :: mode_str; character(8) :: date ; character(10) :: hour
  logical :: pols_active, nbrs_active, scan_noise, write_config_active, ignore_last_config
contains
end module
