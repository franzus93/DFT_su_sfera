!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!					!!!!!!!!!
!!!!!!! 	MODULO_DFT		!!!!!!!!!
!!!!!!!					!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! questo moduto contiene le subrutines e le     !
! funzioni necessarie a completare il programma !
! in python in modo che sia abbastanza veloce   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!La subrutine "convolve" prende le trasformate  !
!di potenziale e distribuzione e moltiplica i   !
!coefficienti tra di loro.			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine convolve( pot_, distro_, convo_,BETA, L_MAX)
  
  implicit none
  
  integer, intent(in) :: L_MAX
  real(8), intent(in) :: BETA
  real(8), intent(in),dimension(0:L_MAX) :: pot_
  real(8), intent(in),dimension(0:1, 0:L_MAX, 0:L_MAX) :: distro_ 
  real(8), intent(inout),dimension(0:1, 0:L_MAX, 0:L_MAX) :: convo_ !f2py intent(in,out) :: convo_
  
  real(8) :: four_pi
  real(8) :: pot_l
  integer :: u, l, m
  
  four_pi = 16.*atan(1.0)*BETA*BETA
  
  do u = 0, 1
       do l = 0, L_MAX
            pot_l = pot_(l)*sqrt(four_pi/dble(2*l+1))           
            do m = 0, l                                         
                 convo_(u,l,m) = distro_(u,l,m)*pot_l
            end do
       end do
  end do  

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine convolve_product calcola la somma !
!del prodotto fra i coefficienti di pot e due   !
!campi eventualmente diversi			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine convolve_product(pot_, distroA_, distroB_, convo, BETA, L_MAX)

  implicit none
  
  integer, intent(in) :: L_MAX
  real(8), intent(in) :: BETA
  real(8), intent(in), dimension(0:L_MAX) :: pot_
  real(8), intent(in), dimension(0:1, 0:L_MAX, 0:L_MAX) :: distroA_, distroB_
  real(8), intent(inout), dimension(0:0) :: convo !f2py intent(in,out) :: convo
  
  real(8) :: four_pi
  real(8) :: pot_l
  integer :: u, l, m
  
  convo(0) = 0.
  four_pi = 16.*atan(1.0)*BETA*BETA

  do u = 0, 1
       do l = 0, L_MAX
            pot_l = pot_(l)*sqrt(four_pi/dble(2*l+1))           
            do m = 0, l                                         
                 convo(0) = convo(0) + distroA_(u,l,m)*pot_l*distroB_(u,l,m)
            end do
       end do
  end do  

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine integrate calcola l'integrale del !
!campo che le viene dato in pasto, sulla sfera  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine integrate(field, integral, N)
  
  implicit none
  
  integer, intent(in) :: N
  real(8), intent(in), dimension(0:N-1, 0:2*N-1) :: field
  real(8), intent(inout), dimension(0:0) :: integral !f2py intent(in,out) :: integral
  
  real(8) :: pi = 4.*atan(1.0)
  real(8) :: piN, sint
  integer :: t,p
  
  integral(0) = 0.
  
  piN = pi/dble(N)
  
  do t = 0, N-1
       sint = piN*piN*sin(piN*t)
       do p = 0, 2*N-1
            integral(0) = integral(0) + sint*field(t,p)
       end do
  end do

end subroutine
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine compute_grad_id calcola la parte  !
!ideale del gradiente del gran potenziale       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_grad_id( distro, grad_id, TRIAL_RHO, BETA, MU_EX, N)
  
  implicit none
  
  integer, intent(in) :: N
  real(8), intent(in) :: TRIAL_RHO, BETA, MU_EX
  real(8), intent(in), dimension(0:N-1, 0:2*N-1) :: distro
  real(8), intent(inout), dimension(0:N-1, 0:2*N-1) :: grad_id !f2py intent(in,out) :: grad_id
  
  integer :: t,p
  real(8) :: lBMU
  
  lBMU = log(TRIAL_RHO) + BETA*MU_EX
  
  do t = 0, N-1
       do p = 0, 2*N-1
            grad_id(t,p) = log(distro(t,p)) - lBMU
       end do
  end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine compute_wz calcola w(0) a partire !
!dai coefficienti di w_(l) del potenziale       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_wz( pot_, wz, L_MAX)
  
  implicit none
  
  integer, intent(in) :: L_MAX
  real(8), intent(inout), dimension(0:0) :: wz !f2py intent(in,out) :: wz
  real(8), intent(in), dimension(0:L_MAX) :: pot_
  
  integer :: l
  real(8) :: four_pi = 16.*atan(1.0)
  
  wz(0) = 0.
  
  do l=0, L_MAX
      wz(0) = wz(0) + sqrt( dble(2*l+1)/four_pi )*pot_(l)
  end do
  
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine compute_precon calcola i precon-  !
!-ditioners.			                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_precon( distro, grad, wz, precon, BETA, N)

  implicit none
  
  integer, intent(in) :: N
  real(8), intent(in) :: wz, BETA
  real(8), intent(in), dimension(0:N-1, 0:2*N-1) :: distro, grad
  real(8), intent(inout), dimension(0:N-1, 0:2*N-1) :: precon  !f2py intent(in,out) :: precon
  
  real(8) :: pi = 4.*atan(1.0)
  real(8) :: dt, wBpiN, wBpiNsint
  integer :: t, p
  
  dt = pi/dble(N)
  
  wBpiN = wz*pi*pi*BETA/dble(N*N)
  
  do t = 0, N-1
       wBpiNsint = wBpiN*sin(pi*dble(t)/dble(N))
       do p = 0, 2*N-1
            precon(t,p) = distro(t,p)*grad(t,p)/( 1. + wBpiNsint*distro(t,p) )
       end do
  end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine calcola zeta, che serve per       !
!calcolare gli psi                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
subroutine compute_zeta( precon, oldcon, zeta, N)

  implicit none
  
  integer, intent(in) :: N
  real(8), intent(in), dimension(0:N-1, 0:2*N-1) :: precon, oldcon
  real(8), intent(inout), dimension(0:0) :: zeta  !f2py intent(in,out) :: zeta
  
  real(8) :: sum_oldcon_sq = 0., sum_oldnew = 0.
  integer :: t,p
  
  
  do t = 0, N-1
       do p = 0, 2*N-1
          sum_oldcon_sq = sum_oldcon_sq + oldcon(t,p)**2
          sum_oldnew    = sum_oldnew + precon(t,p)*( precon(t,p) - oldcon(t,p) )
       end do
  end do
  
  zeta(0) = sum_oldnew/sum_oldcon_sq
  
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!la subrutine update_distro aggiorna la distro  !
!controllando che non diventi mai negativa	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_distro( distro, psi, step, EPSILONN, N)

  implicit none
  
  integer, intent(in) :: N
  real(8), intent(in) :: step, EPSILONN
  real(8), intent(in), dimension(0:N-1, 0:2*N-1) :: psi
  real(8), intent(inout), dimension(0:N-1, 0:2*N-1) :: distro !f2py intent(in,out) :: distro
  
  integer :: t,p
  
  do t = 0, N-1
       do p = 0, 2*N-1
            distro(t,p) = distro(t,p) - step*psi(t,p)
            if ( distro(t,p) .le. 0. ) then
                 distro(t,p) = EPSILONN
            end if
       end do
  end do
  
end subroutine



