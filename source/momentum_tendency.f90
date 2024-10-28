module momentum_tendency

!-----------------------------------------------------------------------
! PURPOSE: Calculates the tendencies of u and v.
!-----------------------------------------------------------------------

use kinds
use model_parameters
use physical_parameters


implicit none
save



contains



!======================================================================
! BEGINNING OF GET_UF_VF
!======================================================================

subroutine get_uf_vf ( u, v, F_x, F_y, h, h_star, pv, zeta, ke_horiz,   &
                       u_f, v_f )

!---------------------------------------------------------------------------
! PURPOSE:
!   Computes the time tendency of the x and y components of velocity,
!   i.e. u and v respectively.
!---------------------------------------------------------------------------

implicit none

!---------------------------------------------------------------------------
! INTENT IN
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(in) ::            &
          u, v           ! horiz. velocity components (m/s)

real (kind=dbl_kind), dimension(im,jm), intent(in) ::            &
          F_x, F_y,   &  ! Mass fluxes in x and y directions
          h,    &        ! height of free surface (m)
          h_star         ! thickness of fluid (m) (h-hs)

!---------------------------------------------------------------------------
! INTENT OUT
!---------------------------------------------------------------------------
real (kind=dbl_kind), dimension(im,jm), intent(out) ::           &
          pv,      &     ! vert. component of pot. vorticity (m-1 s-1)
          zeta,    &     ! relative vorticity (s^-1)
          ke_horiz       ! contribution to kinetic energy
                         ! from horizontal velocity (J/kg)

real (kind=dbl_kind), dimension(im,jm), intent(out) ::           &
          u_f, v_f       ! tendency of u and v (m/s^2)

!---------------------------------------------------------------------------
! LOCAL
!---------------------------------------------------------------------------
integer :: j

real (kind=dbl_kind), parameter ::                                   &
          inv24 = c1/24.00000_dbl_kind

real (kind=dbl_kind), dimension(im,jm) ::                            &
          abs_vort,        &                  ! absolute vorticity (s^-1)
          alfa, beta, gamm, delt, epsln, fi   ! linear combinations of pv





!---------------------------------------------------------------------------
! Calculate the (PV k) cross (m times velocity) term.       (term 1 of 5)
! Potential enstrophy and energy conserving scheme of Arakawa and 
! Lamb (1981) is used.
!---------------------------------------------------------------------------


! calculate potential vorticity
zeta(:,:) = invdy * ( u(:,jm1(:)) - u(:,:) )  +                      &
            invdx * ( v(:,:) - v(im1(:),:) )
do j = 1, jm
   abs_vort(:,j) = f_cor(j) + zeta(:,j)
   pv(:,j) = abs_vort(:,j) / ( p25*(h_star(:,j) + h_star(im1(:),j) + &
                      h_star(im1(:),jm1(j)) + h_star(:,jm1(j))) )
end do



! calculate linear combinations of pv ( eqn. (3.34) of AL (1981) )
alfa(:,:)  = inv24 * ( c2*pv(ip1(:),jp1(:)) + pv(:,jp1(:)) +  &
                       c2*pv(:,:) + pv(ip1(:),:) )
beta(:,:)  = inv24 * ( pv(:,jp1(:)) + c2*pv(im1(:),jp1(:)) +  &
                       pv(im1(:),:) + c2*pv(:,:) )
gamm(:,:)  = inv24 * ( c2*pv(:,jp1(:)) + pv(im1(:),jp1(:)) +  &
                       c2*pv(im1(:),:) + pv(:,:) )
delt(:,:)  = inv24 * ( pv(ip1(:),jp1(:)) + c2*pv(:,jp1(:)) +  &
                       pv(:,:) + c2*pv(ip1(:),:) )
epsln(:,:) = inv24 * ( pv(ip1(:),jp1(:)) + pv(:,jp1(:)) -     &
                       pv(:,:) - pv(ip1(:),:) )
fi(:,:)    = inv24 * (-pv(ip1(:),jp1(:)) + pv(:,jp1(:)) +     &
                       pv(:,:) - pv(ip1(:),:) ) 


!
! calculate u and v tendencies -- term 1 of 5 ( see eqns. (3.5) and (3.6)
!                                 of AL (1981) )

u_f(:,:) = alfa(:,:)*F_y(:,jp1(:)) + beta(:,:)*F_y(im1(:),jp1(:)) +  &
             gamm(:,:)*F_y(im1(:),:) + delt(:,:)*F_y(:,:) -          &
             epsln(:,:)*F_x(ip1(:),:) + epsln(im1(:),:)*F_x(im1(:),:)
v_f(:,:) = -gamm(ip1(:),:)*F_x(ip1(:),:) - delt(:,:)*F_x(:,:) -      &
              alfa(:,jm1(:))*F_x(:,jm1(:)) -                         &
              beta(ip1(:),jm1(:))*F_x(ip1(:),jm1(:)) -               &
              fi(:,:)*F_y(:,jp1(:)) + fi(:,jm1(:))*F_y(:,jm1(:))



!---------------------------------------------------------------------------
! Add contribution of the horiz. gradient of kinetic energy.   (term 2 of 5)
! And also divergence damping!!!
!---------------------------------------------------------------------------


! Calculate contribution to kinetic energy
! from the horizontal velocity  ( see eqn (3.41) of AL (1981) )
! Note:  expect SICK to result from use of this K.E.
ke_horiz(:,:) = p5 * (p5*u(:,:)**2 + p5*u(ip1(:),:)**2) +            &
                p5 * (p5*v(:,:)**2 + p5*v(:,jp1(:))**2)

! "SICK-PROOF" K.E.
! ke_horiz(:,:) = inv8 * ( p25*(u(:,jp1(:))+u(:,jm1(:)))**2 +               &
!                  u(:,:)**2 + p25*(u(ip1(:),jp1(:))+u(ip1(:),jm1(:)))**2 + &
!                   u(ip1(:),:)**2 )   +                                    &
!                   p25 * ( v(:,:)**2 + v(:,jp1(:))**2 )

! Add contribution of horiz. gradient of K.E.
u_f(:,:) = u_f(:,:) - invdx*(ke_horiz(:,:)-ke_horiz(im1(:),:))
v_f(:,:) = v_f(:,:) - invdy*(ke_horiz(:,:)-ke_horiz(:,jm1(:)))



!---------------------------------------------------------------------------
! Add contributions due to vert. advection of horiz. momentum and
! subgrid-scale turbulent momentum flux.            (terms 3 & 4 of 5)
! XXXXXXXXXXXXXXX NOT APPLICABLE TO SHALLOW WATER EQUATIONS XXXXXXXXXXXX
!---------------------------------------------------------------------------


!---------------------------------------------------------------------------
! Add contribution of the horizontal pressure gradient force.  (term 5 of 5)
!---------------------------------------------------------------------------

u_f(:,:) = u_f(:,:) - invdx * grav * ( h(:,:) - h(im1(:),:) )
v_f(:,:) = v_f(:,:) - invdy * grav * ( h(:,:) - h(:,jm1(:)) )




end subroutine get_uf_vf

!======================================================================
! END OF GET_UF_VF
!======================================================================



end module momentum_tendency
