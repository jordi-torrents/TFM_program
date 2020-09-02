module measurements
  use def_variables
contains

  subroutine gdr_nbr_step()
    integer :: i, nbr(Nnbr), indx_nbr
    real(8) :: r, dr, dxy(2)

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
  end subroutine

  subroutine speed_hist_step()
    integer :: i, indx_vel
    do i=1,N
      indx_vel=int(vel(i)/0.02d0)+1
      if (indx_vel<=600) speed_hist(indx_vel) = speed_hist(indx_vel) + 1
    end do
  end subroutine

  subroutine write_configuration()
    integer :: i
    character(4) :: time_str
    write(time_str,'(i4.4)') configuration_counter
    open(unit=823,file=trim(trim(folder)//'/config/'//time_str//'.csv'))
    do i=1,N
      write(823,'(2f8.3)') pos(i,:)
    end do
    close(823)
    configuration_counter = configuration_counter+1
  end subroutine
end module
