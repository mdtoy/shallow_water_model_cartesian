module step

!-----------------------------------------------------------------
!   This module updates the diagnostic variables and then time-
!   steps the prognostic variables.
!   Euler forward is used for the first time step, followed by
!   Adams-Bashforth 2nd order for the second time step, and
!   subsequently by Adams-Bashforth 3rd order time stepping.
!   The "traditional" HPGF horizontal discretization is used.
!-----------------------------------------------------------------

use kinds
use physical_parameters
use model_parameters
use prognostics
use momentum_tendency

implicit none
save


integer :: n3,  &                  ! index for time step n values
                                   ! i.e. previous time step
                                   
           n4                      ! index for time step n+1 values
                                   ! i.e. current time step

integer :: n3_f, n2_f, n1_f        ! index for {n, n-1, n-2} tendencies




contains


!======================================================================
! BEGINNING OF UPDATE_DIAGNOSTICS
!======================================================================

subroutine update_diagnostics ( tau, w1, w2, w3 )

!----------------------------------------------------------------------
! PURPOSE:
!   Updates diagnostic variables and calls subroutines to calculate
!   tendencies of prognostic variables in preparation for
!   time-stepping the prognostic variables.
!----------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
real (kind = dbl_kind), intent(in) :: tau         ! time in hours
real (kind = dbl_kind), intent(in) :: w1, w2, w3  ! Time stepping
                                                  ! weights

!-------------------------------------------------------------------
! LOCAL
!-------------------------------------------------------------------

!
! Declare mass flux variables to be used in the "3rd-order" Takacs
! advection schemes for continuity and theta advection
!
real (kind = dbl_kind), dimension(im,jm) ::                          &
           F_x,          &  ! Mass flux in x-direction
                            ! (colocated with u-points)
           F_y              ! Mass flux in y-direction
                            ! (colocated with v-points)







!
! Get tendency of prognostic mass variable h_star
!
F_x(:,:) = u(:,:,n4)*p5*(h_star(:,:,n4)+h_star(im1(:),:,n4))
F_y(:,:) = v(:,:,n4)*p5*(h_star(:,:,n4)+h_star(:,jm1(:),n4))

h_star_f(:,:,n3_f) =  - invdx*(F_x(ip1(:),:)-F_x(:,:)) -             &
                        invdy*(F_y(:,jp1(:))-F_y(:,:))



! Calculate height of free surface
h(:,:) = hs(:,:) + h_star(:,:,n4)


!
! Get tendencies of prognostic variables u and v
!
call get_uf_vf ( u(:,:,n4),v(:,:,n4),F_x,F_y,h,h_star(:,:,n4),       &
                 pv,zeta,ke_horiz,u_f(:,:,n3_f),v_f(:,:,n3_f) )







end subroutine update_diagnostics

!======================================================================
! END OF UPDATE_DIAGNOSTICS
!======================================================================




!======================================================================
! BEGINNING OF STEP_DYNAMICS
!======================================================================

subroutine step_dynamics ( step_count, w1, w2, w3 )

!----------------------------------------------------------------------
! PURPOSE:
!   Performs the dynamics time stepping using the Adams-Bashworth 3rd.
!   order scheme.
!----------------------------------------------------------------------

implicit none

!-------------------------------------------------------------------
! INTENT IN 
!-------------------------------------------------------------------
integer (kind = int_kind), intent(in) :: step_count

real (kind = dbl_kind), intent(in) :: w1, w2, w3  ! Time stepping
                                                  ! weights




! Advance prognostic time step indices
n4 = mod(step_count+1,2) + 1
n3 = mod(step_count,2) + 1


! Step prognostic variables
u(:,:,n4) = u(:,:,n3) + dt * ( w3*u_f(:,:,n3_f) + w2*u_f(:,:,n2_f) + &
                               w1*u_f(:,:,n1_f) )

v(:,:,n4) = v(:,:,n3) + dt * ( w3*v_f(:,:,n3_f) + w2*v_f(:,:,n2_f) + &
                               w1*v_f(:,:,n1_f) )

h_star(:,:,n4) = h_star(:,:,n3) + dt * ( w3*h_star_f(:,:,n3_f) +     &
                    w2*h_star_f(:,:,n2_f) + w1*h_star_f(:,:,n1_f) )



! Advance tendency time step indices
n3_f = mod(step_count+2,3) + 1
n2_f = mod(step_count+1,3) + 1
n1_f = mod(step_count,3) + 1


end subroutine step_dynamics

!======================================================================
! END OF STEP_DYNAMICS
!======================================================================




end module step
