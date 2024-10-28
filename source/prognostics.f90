module prognostics

!-----------------------------------------------------------------
!
!   Contains arrays related to the prognostic variables, and
!       those diagnostic variables derived from the prognostic
!       variables.  Variables are initialized to zero here.
!       Surface geopotential also declared here.
!
!-----------------------------------------------------------------


use kinds
use model_parameters
use physical_parameters


implicit none
save


!
! Declare prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,ntprog) ::          &
                     u, v,   &   ! horizontal velocity components (m/s)
                     h_star      ! thickness of fluid (m) (h-hs)

!
! Declare tendencies of prognostic variables
!
real (kind = dbl_kind), dimension(im,jm,nttend) ::          &
                     u_f,     &   ! d/dt (u)   (m/s^2)
                     v_f,     &   ! d/dt (v)   (m/s^2)
                     h_star_f     ! d/dt (h_star)   (m/s)

!
! Declare diagnostic variables
!
real (kind = dbl_kind), dimension(im,jm) ::                 &
                     ke_horiz,  &   ! contribution to kinetic energy
                                    ! from horizontal velocity (J/kg)
                     zeta,      &   ! relative vorticity (s^-1)
                     pv,        &   ! vertical component of potential 
                                    ! vorticity (m^2 eta/kg/s)
                     h              ! height of free surface (m)



!
! Declare topographic height
!
real (kind = dbl_kind), dimension (im,jm) ::                &
                     hs             ! topographic height (m)






contains


!======================================================================
! BEGINNING OF INIT_PROGNOSTICS
!======================================================================

subroutine init_prognostics

implicit none

! initialize prognostic arrays
u = c0
v = c0
h_star = c0

u_f = c0
v_f = c0
h_star_f = c0


end subroutine init_prognostics

!======================================================================
! END OF INIT_PROGNOSTICS
!======================================================================


end module prognostics
