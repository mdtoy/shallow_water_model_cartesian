module model_parameters

!-----------------------------------------------------------------
!    This module specifies the grid size, time step length, and
!    run duration.
!-----------------------------------------------------------------


use kinds
use physical_parameters


implicit none
save

!
! Set grid dimensions
!
integer, parameter :: im = 1,    &    ! number of x-direction grid points
                      jm = 80          ! number of y-direction grid points

real (kind=dbl_kind), parameter :: &
                      dx = 100.E+03_dbl_kind,     &   ! x-direction grid spacing (m)
                      dy = 100.E+03_dbl_kind          ! y-direction grid spacing (m)


!
! Set time step length
!
real (kind=dbl_kind), parameter :: dt = 30._dbl_kind      ! time step length (s)


!
! Set termination time
!
real (kind=dbl_kind), parameter ::  &
                       tau_duration = 2400._dbl_kind        ! run duration (hours)


!
! Set frequency of output
! (i.e. number of timesteps per output calls)
!
integer, parameter :: out_freq = 120


real (kind=dbl_kind), parameter :: &
                      invdx  = 1.0_dbl_kind/dx,     &   ! inverse dx
                      invdy  = 1.0_dbl_kind/dy          ! inverse dy


integer, parameter ::   &
           ntprog = 2,  &   ! no. of time slots needed for prognostic variables
           nttend = 3       ! no. of time slots needed for tendency variables



! Declare horizontal grid point "+n/-n" indices:
!  ip1 represents "i + 1", im1 represents "i - 1"
!  ip2 represents "i + 2", im2 represents "i - 2"
!  Periodic boundary conditions are accounted for in that grid point
!  index zero has the value im, and index im+1 gets the value 1, etc.
!  Ditto the above for "j".
!  Indices are set in subroutine initialize_model.
integer, dimension (1:im) :: ip1, im1, ip2, im2
integer, dimension (1:jm) :: jp1, jm1, jp2, jm2



!
! Specify weighting factors for time-stepping
!

! Euler forward factors
real (kind=dbl_kind), parameter ::    &
          w3_ef = 1.0_dbl_kind,       &
          w2_ef = 0.0_dbl_kind,       &
          w1_ef = 0.0_dbl_kind

! Adams-Bashworth 2nd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab2 =  1.5_dbl_kind,     &
          w2_ab2 = -0.5_dbl_kind,     &
          w1_ab2 =  0.0_dbl_kind

! Adams-Bashworth 3rd. order factors
real (kind=dbl_kind), parameter ::    &
          w3_ab3 =  1.91666666666666666666666666666666_dbl_kind,     &
          w2_ab3 = -1.33333333333333333333333333333333_dbl_kind,     &
          w1_ab3 =  0.41666666666666666666666666666666_dbl_kind


real (kind = dbl_kind), dimension(jm) ::                          &
                   f_cor       ! Coriolis parameter (s^-1)


! Domain coordinate (x-y) values
real (kind=dbl_kind), dimension (1:im) :: x_h, x_u
real (kind=dbl_kind), dimension (1:jm) :: y_h, y_v






end module model_parameters
