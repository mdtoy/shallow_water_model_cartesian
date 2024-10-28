module output

!-----------------------------------------------------------------
!   This module outputs the model data in ascii and/or binary 
!   form.
!-----------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters
use prognostics
use step


implicit none


contains



subroutine initial_output (tau)

implicit none

real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

integer :: nt,i,j
integer :: ntm

ntm = int(tau_duration*3600._dbl_kind/(dt*out_freq))+1

write (31) ntm, im, jm
write (31) (tau + nt*tau_duration/(ntm-1),    nt = 0,ntm-1)
write (31) (x_h(i),                  i = 1,im)
write (31) (x_u(i),                  i = 1,im)
write (31) (y_h(j),                  j = 1,jm)
write (31) (y_v(j),                  j = 1,jm)


do j = 1,jm
   write (31) (hs(i,j),        i = 1,im)
end do


end subroutine initial_output



subroutine output_data(stp_cnt,tau)

implicit none


!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
integer (kind = int_kind) :: stp_cnt         ! Number of time step
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

!----------------------------------------------------------------------------
! LOCAL
!----------------------------------------------------------------------------
integer :: i,j




do j = 1, jm
   write (31) (u(i,j,n4),           i = 1,im)
   write (31) (v(i,j,n4),           i = 1,im)
   write (31) (zeta(i,j),           i = 1,im)
   write (31) (pv(i,j),             i = 1,im)
   write (31) (h_star(i,j,n4),      i = 1,im)
   write (31) (h(i,j),              i = 1,im)
end do





call calc_gmeans(tau)

print *, "Output data has been written."




end subroutine output_data





subroutine calc_gmeans (tau)
! Calculates and outputs global mass-weighted means of various variables

implicit none

!----------------------------------------------------------------------------
! INTENT IN 
!----------------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau    ! Time in hours

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: i,j,k

real (kind = dbl_kind) ::                                            &
      h_star_bar,       & ! area-weighted mean mass (m)
      h_star_temp,      &
      ke_bar,           & ! mass-weighted mean kinetic energy (J/kg)
      geop_bar,         & ! mass-weighted mean geopotential energy (J/kg)
      zeta_bar,         & ! mass-weighted mean relative vorticity (s^-1)
      pv_bar,           & ! mass-weighted mean potential vort. (Pa^-1 s^-1)
      pot_enstr_bar,    & ! mass-weighted mean potential enstrophy (Pa^-2 s^-2)
      total_energy_bar    ! mass-weighted mean total energy (J/kg)


! Initialize
h_star_bar = c0
ke_bar = c0
geop_bar = c0
zeta_bar = c0
pv_bar = c0
pot_enstr_bar = c0
total_energy_bar = c0


! Calculate area-weighted mass
do j = 1,jm
   do i = 1,im
      h_star_bar = h_star_bar + h_star(i,j,n4)
   end do
end do
h_star_bar = h_star_bar / (im*jm)



! Calculate global mass-weighted kinetic energy, geopotential energy,
! potential vorticity and enstrophy
do j = 1,jm
   do i = 1,im
      ke_bar = ke_bar + h_star(i,j,n4)*ke_horiz(i,j)
      geop_bar = geop_bar + h_star(i,j,n4)*grav*( p5*h_star(i,j,n4)+ &
                    hs(i,j) )
      h_star_temp = p25*( h_star(i,j,n4) + h_star(im1(i),j,n4) +     &
                     h_star(im1(i),jm1(j),n4) + h_star(i,jm1(j),n4) )
      zeta_bar = zeta_bar + h_star_temp*zeta(i,j)
      pv_bar = pv_bar + h_star_temp*pv(i,j)
      pot_enstr_bar = pot_enstr_bar + h_star_temp*p5*pv(i,j)**2
   end do
end do
ke_bar = ke_bar / ( h_star_bar*im*jm )
geop_bar = geop_bar / ( h_star_bar*im*jm )
zeta_bar = zeta_bar / ( h_star_bar*im*jm )
pv_bar = pv_bar / ( h_star_bar*im*jm )
pot_enstr_bar = pot_enstr_bar / ( h_star_bar*im*jm )


! Calculate global mass-weighted mean total energy (i.e. geopotential
! and kinetic energy)
total_energy_bar = geop_bar + ke_bar




! Write global means to file
write (45, "(F24.10,F24.6,7(ES24.14E3))" )                             &
             tau, tau*3600._dbl_kind, h_star_bar, geop_bar, ke_bar,    &
             total_energy_bar, zeta_bar, pv_bar, pot_enstr_bar



print "(A9,F24.13,A26,F24.13,A19,F24.13)", "ke_bar =", ke_bar,         &
         "total_energy_bar =", total_energy_bar, "h_star_bar =", h_star_bar

end subroutine calc_gmeans





end module output
