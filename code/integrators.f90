module integrators
  use def_variables
  use interaction
contains

  subroutine integrate_simple_vicsek(isteps)
    real(8) :: noise_cnst, theta(N), nbr_direction(N)
    noise_cnst = eta*pi*2.d0
    do j=1,isteps
      call all_neighbours_direction(nbr_direction)
      do i=1,N
        theta(i) = nbr_direction(i)+noise_cnst*(r1279()-0.5)
      end do
      unitary_vel(:,1) = cos(theta)
      unitary_vel(:,2) = sin(theta)

      pos(:,1) = modulo(pos(:,1) + mean_flight*unitary_vel(:,1),L)
      pos(:,2) = modulo(pos(:,2) + mean_flight*unitary_vel(:,2),L)
    end do
    int_time=int_time+isteps
  end subroutine

  subroutine integrate_levy_behaviour(isteps)
    real(8) :: noise_cnst, theta(N), nbr_direction(N)
    noise_cnst = eta*pi*2.d0
    do j=1,isteps
      call all_neighbours_direction(nbr_direction)
      do i=1,N
        theta(i) = nbr_direction(i)+noise_cnst*(r1279()-0.5)
        vel(i) = mean_flight*(A*r1279()+B)**EXPONENT
      end do
      unitary_vel(:,1) = cos(theta)
      unitary_vel(:,2) = sin(theta)

      pos(:,1) = modulo(pos(:,1) + vel*unitary_vel(:,1),L)
      pos(:,2) = modulo(pos(:,2) + vel*unitary_vel(:,2),L)
    end do
    int_time=int_time+isteps
  end subroutine

  subroutine integrate_burstandcoast(isteps)
    integer :: particle_i
    do j=1,isteps
      do while(timer<0.d0)
        particle_i = int(r1279()*dble(N))+1
        vel(particle_i)=vel(particle_i)+burst_amplitude
        timer=timer-log(r1279()+1.e-10)*500.d0
      end do
      call all_neighbours_direction_and_speed()
      vel=vel*0.998d0
      pos(:,1) = modulo(pos(:,1) + vel(:)*unitary_vel(:,1),L)
      pos(:,2) = modulo(pos(:,2) + vel(:)*unitary_vel(:,2),L)
      timer=timer-1.d0
    end do
    int_time=int_time+isteps
  end subroutine

end module
