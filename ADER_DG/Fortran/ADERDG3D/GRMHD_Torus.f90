! ***************************************************************************
#ifdef TORUS
module settings
  IMPLICIT NONE
  integer,parameter :: nx=1,ny=1,nz=300,ngc=3
  integer,parameter :: nx1=1-ngc,nx2=nx+ngc
  integer,parameter :: ny1=1-ngc,ny2=ny+ngc
  integer,parameter :: nz1=1-ngc,nz2=nz+ngc
  integer,parameter :: nv=8,nu=5
  integer,parameter :: krho=1,kvx=2,kvy=3,kvz=4,kp=5,kbx=6,kby=7,kbz=8

end module settings

module initial
  use settings
  IMPLICIT NONE
  real,dimension(nv,nx1:nx2,ny1:ny2,nz1:nz2) ::  vinit,uinit,vold
  real,dimension(3,nx1:nx2,ny1:ny2,nz1:nz2)  ::  other
  real,dimension(nx1:nx2,ny1:ny2,nz1:nz2)    ::  aphi, aphi0, e_field
  real,dimension(nx1:nx2,ny1:ny2,nz1:nz2)    ::  aphix, aphiz
end module initial

module packet_1
  IMPLICIT NONE
  real              :: gam,d,e,s2,b2,sb2,w
end module packet_1

module packet_2
  IMPLICIT NONE
  real        :: gam,sb2,b2,s2,d,e,eps
end module packet_2

module packet_3
  IMPLICIT NONE
  real    ::    gam,d,e,s2,b2,sb2
end module packet_3

module packet_4
  IMPLICIT NONE
  real    ::    d,s,s2,b2,sb2
end module packet_4

module isen
  IMPLICIT NONE
  real    :: u5_isen
    end module isen

! *************************************************************************** 
    
    MODULE DISC_PARAMETER
       USE TypesDef, ONLY : aom, rho_floor,p_floor
        IMPLICIT NONE
        REAL :: r_hrz , gamma,gamma1 
      real :: gamma_th, gamma1_th
      real :: gamma_cd, gamma1_cd
      REAL :: lapse, gp, gm, aphi
      REAL :: g_contr(3,3), g_cov(3,3) ,shift(3), Kex(3,3),phi_tmp 
  real               :: rho,vr,vtheta,vphi,p,br,btheta,bphi
  INTEGER :: iNewton
  INTEGER :: MAXNewton = 15
      real               :: x1, x2, x3,pi,g00,g03,g33,ut_in,tmp
     REAL :: Lcap,eta,w_c,Mcap
     REAL :: tol
      REAL :: rho_c, rho_atmo, p_atmo,r_in,r,ut 
  !real               :: zbrent_c_bis,zbrent_rho_bis
  !real               :: zbrent_rin_bis, zbrent_rout_bis
     
  REAL               :: x(3) 
  !integer :: inew, iz_pos, indeces ( 3)
  !INTEGER :: exact_quartic
  !integer :: eos, rec_flag_x, rec_flag_y, rec_flag_z, rec_flag, atmo

  !real :: Omega,R_star,LC,domeg,Prd,Omega_0
  !real :: tdyn,omeg_f,e_phi
  !real :: omg,dzomg,JJ,vanish,dtp,r_hrz
  !real :: norm0_bx,norm0_bz,norm0_vx,norm0_vy,norm0_vz,norm0_rho
  !real :: coordinates(3), op(9)
  !real :: v_kick, time, M_bh, thick, rho_c, T_c, r_brem, rhou
  !real :: rho_atmo, p_atmo
  !real :: r_border, rho_jump, inject_vel, bphi_par
  !real, parameter ::   p_floor   = 1.0e-30
  !real, parameter ::   rho_floor = 1.0e-30
  !
  !real :: conv_dist, cdelta_x = pi / 1.e+6
     !REAL :: delta_z  
     !REAL ::  xmin, xmax !
     !REAL ::  ymin, ymax !ibcxmin ibcxmax=2
     !REAL ::  zmin
     REAL, PARAMETER :: zmax = 60. ! RMAX
     REAL :: tc,c1,c2  
  real               :: delta_x,delta_z 
  real               :: alpha_c,r_c,ct,st
  real               :: test,fact,w,orbital_vel,pot,t,h,rr,pot_c,p_c,rr_c
  real               :: pot_in,omega_c,rho_atm
  real               :: gamma2,nmn
  real               :: glf,f,df,dt,ur,rcri,vc2,urc,pc
  real               :: pm_c,  pd,rhod
  real               :: cc1
  real               :: vr_down, vtheta_down, vphi_down
  real               :: sigma,zz,zoz,gtt_x,gtt_z,gtp_x,gtp_z,gpp_x,gpp_z
  real               :: bx, by, bz
  real               :: Z1, Z2, rms, rmb, espo, rcusp 
  real               :: pot_cusp, omega_kep, elk 
  real               :: r_out
  real               :: ut_cusp, al_mb, al_ms, utdown, utup
  
  integer            :: ix,iy,iz,n 
  
  logical,parameter :: pe0=.FALSE.,prl=.false.    ! pe0 is for the output
  
  
  INTEGER, PARAMETER :: coordnts = 1  ! 0 for BL coordinates, 1 for KS coord  
  INTEGER, PARAMETER :: atmo = 2      ! 1: Michel atmosphere; 2: Static (Landau); 3: Hawley (quasi-geodetic) atmosphere,  4: Uniform atmosphere
  INTEGER, PARAMETER :: b_pol_type = 0  !0: Wals solution, 1: Isodensity magnetic potential: for iso. stsaggered potential is needed

  
      REAL            :: par_al  = 3.8
      REAL, PARAMETER :: Delta_W = -1e-3   !-0.8
      !REAL, PARAMETER :: rhoc  = 1.146875622133909D-004
      REAL, PARAMETER :: rhoc  = 1D-08      ! solo per la Michel atmosphere 
      REAL, PARAMETER :: t_atm = 1.0   ! only for Landau atmosphere
      REAL, PARAMETER :: jump_l = 1e-6  ! only for Landau atmosphere
      !REAL, PARAMETER :: strc = 1.0
      REAL, PARAMETER :: kpol = 4.96*1e-2    ! politropic constant
      !REAL :: tmax = 1.0
      REAL            :: sigma_c = 0.1   ! beta_c  for the magnetic field    <====  
      REAL, PARAMETER :: b0 = 0.0   
      REAL, PARAMETER :: prtb = 0.0      ! it looks to be a radial constant velocity
      REAL, PARAMETER :: rho_atmo4 = 1e-10       ! 1e-6*rhoc   !rho_floor
      REAL, PARAMETER :: p_atmo4 =0.5*rho_atmo4  !p_floor
  
  CHARACTER(len=16), PARAMETER :: limit ='yes'
    
 END MODULE DISC_PARAMETER
    
subroutine init_disk_isentropic_par(METRIC_loc)
  USE DISC_PARAMETER
  USE TypesDef, ONLY : EQN,xL,xR, b0_torus_par, al_torus_par  
  implicit none

! rho_c is the central density of the disc. You should set
!  rho_c   = 1.1402329E-04 ! to have the model in kerr_A1a.par
!  rho_c   = 1.420174E-05  ! this is the model in ./disk_swallowing
 
  REAL :: zbrent_c_bis
  REAL :: zbrent_rout_bis
  REAL :: zbrent_rin_bis
  REAL :: al_loc  
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 
  
  par_al  = al_torus_par 
  al_loc  = par_al 
  sigma_c = b0_torus_par
  pi=Dacos(-1.)

! Griddelta_x = pi / 1.e+6
     delta_z = 0.5
     
     !xmin= xL(2)   ; xmax=xR(2) ; ! ibcxmin=2; ibcxmax=2
     !ymin= xL(3)          ; ymax= xR(3)     ;! ibcymin=1; ibcymax=1
     !zmin= xL(1) ; ! zmax= xR(1)     ; !ibczmin=5; ibczmax=6

  !
  IF ( coordnts == 0) THEN  ! Boyer Lindquist coordinates

        r_hrz  = 1.0d0 + sqrt(1.0d0**2 - aom**2)

  ELSE                      ! KS   coordinates

        !r_hrz  = 1.0
        r_hrz  = 1.0d0 + sqrt(1.0d0**2 - aom**2)

  ENDIF


! Initial conditions

  gamma  = EQN%gamma !4.0 / 3.0
  gamma1 = gamma / ( gamma - 1.0)
  gamma2 = 1.0 / ( gamma - 1.0)

  gamma_th = 4./3.; gamma1_th = gamma_th / ( gamma_th - 1.0)
!  gamma_th = 2.0; gamma1_th = gamma_th / ( gamma_th - 1.0)
  gamma_cd = 4./3.; gamma1_cd = gamma_cd / ( gamma_cd - 1.0)

  espo  = 1.0d0/3.0d0
  Z1    = 1.0d0 + (1.0d0 - aom*aom)**espo*((1.0d0 + aom)**espo + &
       & (1.0d0 - aom)**espo)
  Z2    = sqrt(3.0d0*aom*aom + Z1**2)
  
  rmb   = 2.0d0*(1.0d0 - aom/2.0d0 + sqrt(1.0d0 - aom))
  rms   = (3.0d0 + Z2 - sqrt((3.0d0 - Z1)*(3.0d0 + Z1 + 2.0d0*Z2)))

  !if ( pe0) then
  !   write(*,*)'r_hrz,rmb,rms',r_hrz,rmb,rms
  !endif

! Initial conditions for Michel's inflow: fix kpol and rhoc
  n    = 3   !from Olindo was ng     = 1.0/(EQN%gamma - 1.0) =3 only if gamma = 4/3
  vc2  = ( ( n + 1) * kpol * rhoc** ( gamma - 1.0))/  &
       & ( n * ( 1.0 + ( n + 1) * kpol * rhoc** ( gamma - 1.0)))
  tc   = n * vc2 / ( ( 1.0 + n) * ( 1.0 - n * vc2))
  pc   = rhoc * tc
  urc  = sqrt ( vc2 / ( 1.0 + 3.0 * vc2))
  rcri = 1.0 / ( 2.0 * urc**2)
  c1   = urc * ( tc / kpol)**n * rcri**2
  c2   = ( 1.0 + ( 1.0 + n) * tc)**2* ( 1.0 - 2.0 / rcri + urc**2)

  !if ( pe0) then
  !    write(*,*)'zmin,zmax,delta_z',zmin,zmax,delta_z  
  !    write(*,*)'kpol',kpol
  !    write(*,*)'rcri',rcri
  !    write(*,*)'c1,c2',c1,c2
  !endif


  ! Parameters of the disc

  ! Evaluate l_mb
  x1 = rmb
  x2 = pi / 2.0
  x3 = 0.0d0
  !call system_metric( x2, x3, x1,(/1,2,3/))   ! COMPUTE THE METRIC 
  
  
	! Prepare quantities:
    x = (/x1,x2,x3 /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp  )
    !
  
  
  
  !g33     =  g_cov ( 2) ! this is g_{phi,phi}
  !g03     = g_cov ( 2) * beta ( 2) + g_cov ( 4) * beta ( 3)
  !g00     = -alpha**2 + g_cov ( 1) * beta ( 1)**2 + &
  !     &                g_cov ( 2) * beta ( 2)**2 + &
  !     &                g_cov ( 3) * beta ( 3)**2 + &
  !     &          2.0 * g_cov ( 4) * beta ( 2) * beta ( 3)
  
  !
  g33     = g_cov(3,3) ! g_cov ( 2) ! this is g_{phi,phi}
  g03     = g_cov(3,3) * shift (3) + g_cov(1,3) * shift (1)
  g00     = -lapse**2 + g_cov(2,2) * shift(2)**2 + &
       &                g_cov(3,3) * shift (3)**2 + &
       &                g_cov(1,1) * shift (1)**2 + &
       &          2.0 * g_cov(1,3) * shift (3)*shift (1)
  
  
  omega_kep  = 1.0 / ( rmb**1.5 + aom)
  al_mb = - ( g03 + omega_kep * g33) / ( g00 + omega_kep * g03)


  ! Evaluate l_ms
  x1 = rms
  x2 = pi / 2.0
  x3 = 0.0d0
  !call system_metric( x2, x3, x1,(/1,2,3/))
	! Prepare quantities:
    x = (/x1,x2,x3 /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
    !
  
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)
  !
  omega_kep  = 1.0 / ( rms**1.5 + aom)
  al_ms = - ( g03 + omega_kep * g33) / ( g00 + omega_kep * g03)

  !if ( pe0) then
  !   write(*,*)'el',al
  !   write(*,*)'el_mb',al_mb
  !   write(*,*)'el_ms',al_ms
  !endif

!!$  if ( pe0) then
!!$     if ( al < al_ms .or. al > al_mb) then
!!$        
!!$        write(*,*)'I suggest a value of al between al_ms, al_mb'
!!$        write(*,*)'continue anyway?'
!!$        pause
!!$
!!$     endif
!!$  endif

  !! Plot the Keplerian angular momentum and potential 
  !    open(unit=56,file='elk_in.dat',status='unknown')
  !    open(unit=74,file='w_eq_in.dat',status='unknown') 
!  !
!  do iz = 1, nz
!
!     r = z ( iz)
!
!     !call system_metric(pi/2.,0.,r,(/1,2,3/))
!     
!    x = (/r ,pi/2.,0./)
!    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi )
!     
!     
!     g33     = g_cov (3,3) ! this is g_{phi,phi}
!     g03     = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
!     g00     = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
!          &                g_cov (3,3) * shift (3)**2 + &
!          &                g_cov (1,1) * shift (1)**2 + &
!          &          2.0 * g_cov (1,3) * shift (3) * shift (1)
!     !
!     omega_kep  = 1.0 / ( r**1.5 + aom)
!     elk        = - ( g03 + omega_kep * g33) / ( g00 + omega_kep * g03)
!
!     omega = - ( g03 + g00 * al) / ( g33 + g03 * al)
!
!     tmp  = ( -1.0 / ( g00 + 2.0 * omega * g03 + omega**2 * g33))
!
!     utup     = sqrt ( abs ( tmp))
!     utdown   = - 1.0 / ( 1.0 - omega * al) / utup
!     pot      = log ( abs( utdown))
!!     write(*,*)'W1',pot
!
!!!$     rr       = g03**2 - g00*g33
!!!$     pot      = 0.5 * log ( abs ( rr / ( g00 * al**2 + 2.0 * g03 * al + g33)))
!!!$!     pot      = 0.5 * log ( rr / ( g00 * al**2 + 2.0 * g03 * al + g33))
!!!$     write(*,*)'W2',pot
!!!$     write(*,*)'----------'
!     if ( pe0) then
!        write ( 56, *) r, elk, iz
!        write ( 74, *) r, pot, iz
!     endif
!  enddo
!  close ( 56)


  ! Find the position of the first keplerian point
  IF ( coordnts == 0) x1    = r_hrz + 1.0d-9
  IF ( coordnts == 1) x1    = 2.01 ! empirical
  x2    = rms
  tol   = 1.d-7
  rcusp = zbrent_c_bis ( x1, x2, al_loc, tol,METRIC_loc)
  if ( pe0) write(*,*)'rcusp(as first keplerian point)=',rcusp


  ! Compute the potential at rcusp
  !call system_metric(pi/2.,0.,rcusp,(/1,2,3/))
  
    x = (/rcusp ,pi/2.,0./)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
     
      
  g33     = g_cov (3,3) ! this is g_{phi,phi}
  g03     = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)!recall g03=beta_phi 
  g00     = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &                g_cov (3,3) * shift (3)**2 + &
       &                g_cov (1,1) * shift (1)**2 + &
       &          2.0 * g_cov (1,3) * shift (3) * shift (1)
  rr       = g03**2 - g00 * g33
  pot_cusp = 0.5 * log ( abs ( rr / ( g00 * al_loc**2 + 2.0 * g03 * al_loc + g33)))
  ut_cusp  = exp ( pot_cusp)
  if ( pe0) write(*,*)'pot_cusp',pot_cusp   

  ! Find the position of the second keplerian point
  x1  = rms
  x2  = 60.0d0
  tol = 1.d-7
  r_c = zbrent_c_bis ( x1, x2, al_loc, tol,METRIC_loc)
  omega_c = 1.0 / ( r_c**1.5 + aom)
  if ( pe0) write(*,*)'r_c',r_c


  ! Compute potential at centre
  !call system_metric(pi/2.,0.,r_c,(/1,2,3/))
  
  
    x = (/r_c ,pi/2.,0./)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
     
     
  alpha_c = lapse
  !
  g33     = g_cov (3,3) ! this is g_{phi,phi}
  g03     = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)!recall g03=beta_phi 
  g00     = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &                g_cov (3,3) * shift (3)**2 + &
       &                g_cov (1,1) * shift (1)**2 + &
       &          2.0 * g_cov (1,3) * shift (3) * shift (1)
  !
  rr_c    = g03**2 - g00 * g33
  pot_c   = 0.5 * log ( rr_c / ( g00*al_loc**2 + 2.0*g03*al_loc + g33))
  if ( pe0) write(*,*)'pot_c',pot_c   

  ut_in  = ut_cusp * exp ( Delta_W)
  pot_in = log ( ut_in)


  IF ( al_loc < al_ms .or. al_loc > al_mb) then

     ! Since there is no cusp, it is useless to compute
     ! W_cusp and Delta_W
     r_in  = 40.0

     !call system_metric(pi/2.,0.,r_in,(/1,2,3/))
     ! 
     x = (/r_in ,pi/2.,0./)
     CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
     !
      alpha_c = lapse
      !
      g33     = g_cov (3,3) ! this is g_{phi,phi}
      g03     = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)!recall g03=beta_phi 
      g00     = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &                g_cov (3,3) * shift (3)**2 + &
       &                g_cov (1,1) * shift (1)**2 + &
       &          2.0 * g_cov (1,3) * shift (3) * shift (1)
     rr_c    = g03**2 - g00 * g33
     pot_in  = 0.5 * log ( rr_c / ( g00*al_loc**2 + 2.0*g03*al_loc + g33))

     x1    = r_in + 1.0
     x2    = zmax
     tol   = 1.d-7
     r_out = zbrent_rout_bis( x1, x2, tol, r_in, al_loc,METRIC_loc)       

  ELSE 
      ! Compute the external radius
      if ( Delta_W > 0.0d0) then
         ! Overflow
         if ( pe0) then
            write(*,*)'roche lobe overflow!'
            write(*,*)'pot_in (from energy gap)',pot_in
         endif
         r_in  = rcusp
         x1    = r_in
         x2    = 50.
         tol   = 1.d-7
         r_out = zbrent_rin_bis ( x1, x2, tol, ut_in, al_loc,METRIC_loc)        

      elseif (Delta_W.eq.0.0d0) then
         ! Filling Roche lobe exaclty!
         if ( pe0) then
            write(*,*)'filling roche lobe exactly!'
            write(*,*)'pot_in (from energy gap)',pot_in
         endif
         r_in  = rcusp
         x1    = r_c
         x2    = 50.
         tol   = 1.d-5
         r_out = zbrent_rout_bis( x1, x2, tol, r_in, al_loc,METRIC_loc)       

      else
         ! Torus inside Roche lobe
         if ( pe0) then
            write(*,*)'torus inside roche lobe!'
            write(*,*)'pot_in (from energy gap)',pot_in
         endif
         x1    = rcusp
         x2    = r_c
         tol   = 1.0d-5
         r_in  = zbrent_rin_bis ( x1, x2, tol, ut_in, al_loc,METRIC_loc)        
         x1    = r_c
         x2    = 50.
         tol   = 1.d-3
         r_out = zbrent_rout_bis ( x1, x2, tol, r_in, al_loc,METRIC_loc)        

      endif

  ENDIF


  ! sigma = pmag/pgas = 1/beta: sigma = 0 for B = 0
  eta     = gamma

  rho_c = ((gamma - 1.0)/(kpol*gamma) * ( exp ( pot_in - pot_c) - 1.0d0 ))**  &
       &   ( 1.0d0 / ( gamma - 1.0))

  p_c = kpol * rho_c**gamma

  pm_c  = p_c*sigma_c

  if (pe0) write(*,*)'p_c',p_c
  if (pe0) write(*,*)'pm_c',pm_c

  w_c   = rho_c + gamma*p_c/(gamma - 1.0)

  rr_c  = g03**2 - g00*g33
  Mcap  = pm_c / (( w_c**gamma) * rr_c**( gamma - 1.d0))

  if (pe0) write(*,*) 'Lcap,Mcap',rr_c,Mcap
  ! The definition of  
  !    rr_c = alpha**2*g33 
  ! is equivalent to 
  !    rr_c = g03**2 - g00*g33
  ! only in BL coordinates. To be general, we use here the most general:

  if (pe0) write(*,*) 'rr_c',rr_c
  pot_c   = 0.5 * log ( rr_c / ( g00 * al_loc**2 + 2.0 * g03 * al_loc + g33))
  ! pot_c   = 0.5 * log ( rr_c / ( g00*al**2+2.*g03*al+g33))
  ! pot_c   = log(alpha/sqrt((1.+beta(2)*al)**2-alpha**4*al**2/rr_c))
  
  if (pe0) write(*,*)'pot_c',pot_c

  !nmn  = tmax
  !tdyn = 2.0 * pi / omega_c  
  !tmax = nmn * tdyn
  !
  !if (pe0) then
  !    write(*,*)'pot_in',pot_in
  !    write(*,*)'pot_in - pot_cusp',pot_in - pot_cusp
  !    write(*,*)'r_in',r_in
  !    write(*,*)'r_out',r_out
  !    write(*,*)'r_in [Km]', r_in*n_solar*1.4766
  !    write(*,*)'r_out [Km]',r_out*n_solar*1.4766
  !
  !    write(*,*)'rho_c[geo],rho_c[cgs]',real(rho_c),real(rho_c*conv_rho)
  !    print*,'l_0 = ',al,'P_c = ',2.*pi/omega_c
  !    print*,'D W = ',pot_in-pot_c,'sig = ',sigma_c
  !    print*,'tdyn[geo],tdyn[sec]',tdyn,tdyn*conv_time
  !    print*,'tmax = ',tmax,'tmax/P = ',nmn
  !endif

!
  continue
  !
    END SUBROUTINE init_disk_isentropic_par


    
SUBROUTINE init_disk_isentropic(xp, V, METRIC_loc) ! , tmax)
  USE DISC_PARAMETER
  USE typesDef, ONLY: aom, Mbh
  implicit none
  REAL :: V(8)
  REAL :: xp(3) 
  real :: zbrent_c_bis,zbrent_rho_bis
  real               :: zbrent_rin_bis, zbrent_rout_bis, al_loc  
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 
  
  !
  al_loc = par_al 
  
     r  = xp(1) !z ( iz)
     ct = cos (xp(2)) ! xp0 ( ixp))
     st = sin (xp(2)) !xp0 ( ixp))
     
     !call system_metric( x0 ( ix), y ( iy), z ( iz), (/1,2,3/)) 
    CALL METRIC_loc( xp, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
      
    
  ! ==============================================================
  ! Loop 1
  ! This is just to compute A_phi in case of b_pol_type=1,
  ! since this is based on rho computation.

  ! Aphi is staggered in both r and theta, centered in phi

    IF ( b_pol_type == 1) THEN   ! Isodensity magnetic potential 
 
        PRINT *, 'Isodensity magnetic potential NOT IMPLEMENTED!'
        STOP

    r  = xp(1) !z0 ( iz)
    ct = cos (xp(2))  
    st = sin (xp(2))  
 
    CALL METRIC_loc( xp, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
      

    if ( atmo == 1) then
        ! -------------------------------
        ! Michel's atmosphere

        t=tc
        do
           ur = c1/(r**2*(t/kpol)**n)
           f  = (1. + (1 + n)*t)**2*(1.0 - 2.0/r + ur**2) - c2
           df=2*(1+n)*(1.+(1+n)*t)*(1.-2./r+ur**2)-2*n*ur**2/t*(1.+(1+n)*t)**2
           dt = -f/df
           if (abs(dt)<1.e-18) exit
           t = t + dt
        end do
        rho = rhoc*(t/tc)**n
        p   = rho*t  ! or p = kpol*rho**gamma (it is the same)
        ut  = sqrt(g_cov(1,1) + (g_cov(1,1)*ur)**2)
        glf = lapse*ut
        vr  = -ur/glf
        vtheta = 0.0
        vphi   = 0.0
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    elseif ( atmo == 2) then
        ! -------------------------------
        ! Static (Landau) atmosphere

        rho_atm= jump_l * rho_c
        h      = (1.0 + gamma1*t_atm) * alpha_c/lapse
        t      = (h - 1.0) / gamma1
        rho    = rho_atm * ( t / t_atm)**gamma2
        p      = rho*t
        vr     = 0.0
        vtheta = 0.0
        vphi   = 0.0
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

        rho_atmo = rho
        p_atmo   = p

    elseif ( atmo == 3) then
        ! -------------------------------
        ! Hawley (quasi-geodetic) atmosphere
        cc1    = 1.0d-9
        vtheta = 0.0
        vphi   = 0.0
        vr     = -(1.d0/lapse)*sqrt(2.d0/r) / g_cov(1,1)
        rho    = cc1/(r*r*sqrt(2.d0/r))
        p      = kpol*rho**gamma
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    elseif ( atmo == 4) then 
        ! -------------------------------
        ! Uniform atmosphere
        rho_atmo = rho_atmo4  ! this is used in system_cons_to_prim
        p_atmo   = p_atmo4    ! this is used in system_cons_to_prim

        vtheta = 0.0
        vphi   = 0.0
        vr     = 0.0
        rho    = rho_atmo4
    !    p      = kpol*rho**gamma
        p      = p_atmo4
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    endif
! End of "if atmo = "


    ! -------------------------------
    ! Now build the torus
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)
    rr    = g03**2 - g00*g33
    pot   = 0.5 * log ( abs( rr / ( g00 * al_loc**2 + 2.0 * g03 * al_loc + g33)))

    orbital_vel = - ( g03 + g00 * al_loc) / ( g33 + g03 * al_loc)
    test        = g00+2.*orbital_vel*g03+orbital_vel**2*g33
    
    Lcap  = rr


    if (r > r_in .and. pot < pot_in .and. test < 0.) then ! Disk
       
       fact = sigma_c*(rr/rr_c)**(gamma-1.)

       x1   = 1.0d-25
       x2   = 1.0d+25
       tol  = 1.d-16
       rhod = zbrent_rho_bis ( x1, x2, tol, eta, pot, pot_in, Lcap, Mcap)
       pd   = kpol * rhod**gamma

       IF ( limit=='yes') THEN

!          if ( rhod > rho) then
          if ( pd > p) then
             rho    = rhod
             p      = pd
             vr     = shift (1) / lapse
             vtheta = 0.0

             vphi = (orbital_vel + shift(3))/lapse
             bphi = sqrt ( 2.0 * p * fact / g33)
          end if

       ELSE

          rho    = rhod
          p      = pd
          vr     = shift (1) / lapse
          vtheta = 0.0

          vphi = (orbital_vel + shift(3))/lapse
          bphi = sqrt ( 2.0 * p * fact / g33)

       ENDIF
       
    end if


    ! Compute vector potential for poloidal magnetic field:


    IF ( b_pol_type == 0) THEN  ! Wald m.f. 
       aphi = b0 / 2.0d0 * xp(1)**2* ( sin (xp(2)))**2 
    ELSE                        ! Isodensity magnetic potential 
       aphi = b0 * max( rho / rho_c - 0.2d0, 0.0d0) 
    ENDIF



 !end do
 !end do
 !end do
     
    ENDIF


! Loop 2

!======================================================================
! Here we set the staggered m.f.
!======================================================================
!----------------------------------------------
! Isodensity magnetic potential

  IF ( b_pol_type == 1) THEN  ! 

 
        PRINT *, 'Isodensity magnetic potential NOT IMPLEMENTED!'
    STOP
    
    !V(6) = - ( aphi ( ix, iy, iz + 1) - aphi ( ix, iy, iz)) * &
    !        &             ddz ( iz)
    !
    !
    !
    !V(8) = ( aphi ( ix + 1, iy, iz) - aphi ( ix, iy, iz)) * &
    !        &             ddx ( ix)

 

  ENDIF   




!=============================================================
! Loop 3
! This is the true loop for building initial conditions 

    if ( atmo == 1) then
        ! -------------------------------
        ! Michel's atmosphere

        t=tc
        DO iNewton = 1, MAXNEWTON 
           ur = c1/(r**2*(t/kpol)**n)
           f  = (1. + (1 + n)*t)**2*(1.0 - 2.0/r + ur**2) - c2
           df=2*(1+n)*(1.+(1+n)*t)*(1.-2./r+ur**2)-2*n*ur**2/t*(1.+(1+n)*t)**2
           dt = -f/df
           if (abs(dt)<1.e-10) exit
           t = t + dt
        end do
        rho = rhoc*(t/tc)**n
        p   = rho*t  ! or p = kpol*rho**gamma (it is the same)
        ut  = sqrt(g_cov(1,1) + (g_cov(1,1)*ur)**2)
        glf = lapse*ut
        vr  = -ur/glf
        vtheta = 0.0
        vphi   = 0.0
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    elseif ( atmo == 2) then
        ! -------------------------------
        ! Static (Landau) atmosphere

        rho_atm= jump_l * rho_c
        h      = (1.0 + gamma1*t_atm) * alpha_c/lapse
        t      = (h - 1.0) / gamma1
        rho    = rho_atm * ( t / t_atm)**gamma2
        p      = rho*t
        vr     = 0.0
        vtheta = 0.0
        vphi   = 0.0
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    elseif ( atmo == 3) then
        ! -------------------------------
        ! Hawley (quasi-geodetic) atmosphere
        cc1    = 1.0d-10
        vtheta = 0.0
        vphi   = 0.0
        vr     = -(1.d0/lapse)*sqrt(2.d0/r) / g_cov(1,1)
        rho    = cc1/(r*r*sqrt(2.d0/r))
        p      = kpol*rho**gamma
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    elseif ( atmo == 4) then
        ! -------------------------------
        ! Uniform atmosphere
        rho_atmo = rho_atmo4  ! this is used in system_cons_to_prim
        p_atmo   = p_atmo4    ! this is used in system_cons_to_prim

        vtheta = 0.0
        vphi   = 0.0
        vr     = 0.0
        rho    = rho_atmo4
    !    p      = kpol*rho**gamma
        p      = p_atmo4
        br     = 0.0
        btheta = 0.0
        bphi   = 0.0

    endif
! End of "if atmo = "


    ! -------------------------------
    ! Now build the torus
      g33 = g_cov (3,3) ! this is g_{phi,phi}
      g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
      g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
           &            g_cov (3,3) * shift (3)**2 + &
           &            g_cov (1,1) * shift (1)**2 + &
           &      2.0 * g_cov (1,3) * shift (3) * shift (1)
    rr    = g03**2 - g00*g33
    pot   = 0.5 * log ( abs( rr / ( g00 * al_loc**2 + 2.0 * g03 * al_loc + g33)))

    orbital_vel = - ( g03 + g00 * al_loc) / ( g33 + g03 * al_loc)
    test        = g00+2.*orbital_vel*g03+orbital_vel**2*g33
    
    Lcap  = rr


    if (r > r_in .and. pot < pot_in .and. test < 0.) then ! Disk
       
       fact = sigma_c*(rr/rr_c)**(gamma-1.)

       x1   = 1.0d-25
       x2   = 1.0d+25
       tol  = 1.d-16
        !IF(FLAGtmp) THEN
        !    open(unit=123,file='TMP.dat'  ,status='unknown') 
        !    write(123,'( E16.6, E16.6, E16.6,E16.6,E16.6,E16.6,E16.6,E16.6)') x1, x2, tol, eta, pot, pot_in, Lcap, Mcap
        !    CLOSE(123)
        !    STOP
        !ENDIF
       rhod = zbrent_rho_bis( x1, x2, tol, eta, pot, pot_in, Lcap, Mcap)
       pd   = kpol * rhod**gamma

       IF ( limit=='yes') THEN

!          if ( rhod > rho) then
          if ( pd > p) then
             rho    = rhod
             p      = pd
             vr     = shift (1) / lapse + prtb
             vtheta = 0.0

             vphi = (orbital_vel + shift(3))/lapse
             bphi = sqrt ( 2.0 * p * fact / g33)
          end if

       ELSE

          rho    = rhod
          p      = pd
          vr     = shift (1) / lapse + prtb
          vtheta = 0.0

          vphi = (orbital_vel + shift(3))/lapse
          bphi = sqrt ( 2.0 * p * fact / g33)

       ENDIF
       
    end if

!!!    ! Introduce a poloidal magnetic field:
!!!
!!!
!!!    IF ( b_pol_type == 0) THEN  ! Wald m.f.
!!!
!!!
!!!        ! compute metric term derivatives
!!!        sigma = r**2 + aom**2*ct**2
!!!        zz    = 2.0*r/sigma
!!!        zoz   = (r**2 - aom**2*ct**2)/sigma**2
!!!    
!!!        gtt_z = -2.0*zoz
!!!        gtp_z = -aom*st**2*gtt_z
!!!        gpp_z = 2.0*st**2*(r - aom**2*st**2*zoz)
!!!    
!!!        gtt_x = 2.0*zz/sigma*aom**2*st*ct
!!!        gtp_x = -2.0*aom*zz*st*ct*(r**2 + aom**2)/sigma
!!!        gpp_x = 2.0*st*ct*(sigma + aom**2*st**2*(1.0 + 2.0*zz) + &
!!!             &             aom**4*st**4/sigma*zz)
!!!    
!!!        ! compute the magnetic field
!!!        if (aom == 0.0) then
!!!       
!!!           bx  = -b0 / gp * r   *st**2
!!!           bz  =  b0 / gp * r**2*ct*st
!!!       
!!!        else
!!!       
!!!           bx =  -b0 * gm * ( 0.5 * gpp_z + aom * gtp_z)
!!!           bz =   b0 * gm * ( 0.5 * gpp_x + aom * gtp_x)
!!!       
!!!        endif
!!!    
!!!
!!!    ELSE  ! Isodensity magnetic potential
!!!        PRINT *, 'Isodensity magnetic potential NOT IMPLEMENTED!'
!!!       STOP
!!!
!!!       !bx = .5*(u(kbx,ix-1,iy,iz)+u(kbx,ix,iy,iz))
!!!       !bz = .5*(u(kbz,ix,iy,iz-1)+u(kbz,ix,iy,iz))
!!!       !
!!!       !bx = gm * bx
!!!       !bz = gm * bz
!!!
!!!    ENDIF
!!!! End of 'IF ( b_pol_type == 0) THEN'
!!!
!!!
!!!    br     = bz
!!!    btheta = bx
!!!    
    !!*************************************  END TORUS AND MAG. FIELD
    !IF(FLAGtmp) THEN
    !    open(unit=123,file='TMP.dat'  ,status='unknown')
    !    write(123,'( E16.6, E16.6, E16.6,E16.6,E16.6)') xp, rho , rhod
    !    write(123,'( E16.6, E16.6, E16.6,E16.6,E16.6,E16.6)') g_cov(1:3,1),g_cov(2:3,2),g_cov(3,3)
    !    CLOSE(123)
    !    STOP
    !ENDIF
    
    ! Now set velocity to covariant:
    vr_down     = g_cov (1,1) * vr     + g_cov (1,3) * vphi
    vtheta_down = g_cov (2,2) * vtheta 
    vphi_down   = g_cov (3,3) * vphi   + g_cov (1,3) * vr

    V(1:8) = (/ rho, vr_down, vtheta_down, vphi_down, p, br , btheta, bphi/)    

!!$    ! Keep contravariant velocity
!!$    vv=(/ rho,vtheta,vphi,vr,p,btheta,bphi,br /)

    !call system_prim_to_cons ( vv, uu)
    
    
    !u (1:nu,ix,iy,iz) = uu
    !v0(1:nv,ix,iy,iz) = vv
    
    !IF ( ny == 1) THEN
    !   ! The original one, for runs in (theta, r)
    !
    !   u(kby,ix,0,iz) = bphi * gp
    !   u(kby,ix,1,iz) = bphi * gp
    !   
    !ELSE
    !   
    !   u(kby,ix,0,iz) = bphi * gp
    !   u(kby,ix,iy,iz) = bphi * gp
    !   
    !ENDIF


!end do
!end do
!end do









300 format(1x,(i9,1x),1(E23.9,1x))
301 format(E23.9)
350 format(1x,(i9,1x),2(E23.9,1x))

end subroutine init_disk_isentropic

! ***************************************************************************

FUNCTION zbrent_rho_bis(x1,x2,tol,eta,W_c,W_in,Lcap,Mcap)
 ! USE nrutil, ONLY : nrerror
  !use system
  IMPLICIT NONE
  REAL, INTENT(IN) :: x1,x2,tol,W_c,W_in,eta
  REAL :: zbrent_rho_bis,fo,Lcap,Mcap
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  INTERFACE
     subroutine func_rho_bis(rho,eta,W_c,W_in,Lcap,Mcap,f_out)
       !use system
       Implicit none
       real          :: rho,f_out,eta
       real          :: tmp1,tmp3,Lcap
       real          :: W_c,W_in,Mcap,gamma_mo
       Intent(IN) :: rho,eta,W_in,W_c,Lcap,Mcap
       Intent(OUT):: f_out
     end subroutine func_rho_bis
  END INTERFACE
!  write(*,*)'params',beta_c,eta,W_c,W_in,Lcap,Mcap
  a=x1
  b=x2
  call func_rho_bis(a,eta,W_c,W_in,Lcap,Mcap,fo)
  fa = fo
  call func_rho_bis(b,eta,W_c,W_in,Lcap,Mcap,fo)
  fb = fo
!!$  write(*,*)'a,fa',real(a),real(fa)
!!$  write(*,*)'b,fb',real(b),real(fb)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0))  then
     write(*,*)'root must be bracketed for zbrent_rho_bis'
     stop
  endif
  c=b
  fc=fb
  do iter=1,ITMAX
     if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_rho_bis=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     call func_rho_bis(b,eta,W_c,W_in,Lcap,Mcap,fo)
     fb = fo
  end do
!  call nrerror('zbrent_rho_bis: exceeded maximum iterations')
  zbrent_rho_bis=b
END FUNCTION zbrent_rho_bis



subroutine func_rho_bis(rho,eta,W_c,W_in,Lcap,Mcap,f_out)
  !use system
  USE DISC_PARAMETER, ONLY: kpol
  USE TypesDef, ONLY : EQN
  Implicit none
  real          :: rho,f_out,eta
  real          :: tmp1,tmp3,Lcap
  real          :: W_c,W_in,Mcap,gamma_mo,gamma
  Intent(IN) :: rho,eta,W_in,W_c,Lcap,Mcap
  Intent(OUT):: f_out

!  write(*,*)'from func_rho_bis:params',sigma_c,eta,W_c,W_in
  gamma = EQN%gamma
  gamma_mo = gamma - 1.0

  tmp1   = 1.0d0 + kpol*gamma/gamma_mo*rho**gamma_mo

  ! tmp3 is (h*rho)
  tmp3  = rho + gamma/gamma_mo*kpol*rho**gamma
  
  f_out = W_c - W_in + log(tmp1) +                   &
       &  eta/(eta - 1.0d0)*Mcap*Lcap**(eta - 1.0d0)*tmp3**(eta - 1.0d0)

    end subroutine func_rho_bis





    
    
    
    
    
FUNCTION zbrent_c_bis ( x1, x2, el0, tol,METRIC_loc)
  USE recipies_mod, ONLY : nrerror
  IMPLICIT NONE
  REAL, INTENT(IN) :: x1, x2, el0, tol
  REAL :: zbrent_c_bis,fo
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  INTERFACE
     subroutine func_c_bis ( xg, el0, f_out,METRIC_loc)
       !Use system
       Implicit none
       real          :: xg, el0, f_out
       real          :: Omega_k, el_k
       real          :: x1, x2, x3, pi
       real          :: g00, g03, g33
       
       Intent(IN) :: xg, el0
       Intent(OUT):: f_out
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

     end subroutine func_c_bis
  END INTERFACE
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

  a=x1
  b=x2
  call func_c_bis(a, el0, fo,METRIC_loc)
  fa = fo
  call func_c_bis(b, el0, fo,METRIC_loc)
  fb = fo
!!$  write(*,*)'a,fa',a,fa
!!$  write(*,*)'b,fb',b,fb
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
     write(*,*)'root must be bracketed for zbrent_c_bis'
     stop
  endif
  c=b
  fc=fb
  do iter=1,ITMAX
     if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_c_bis=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
           q=(q-1.0)*(r-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     call func_c_bis(b, el0, fo,METRIC_loc)
     fb = fo
  end do
  call nrerror('zbrent_c_bis: exceeded maximum iterations')
  zbrent_c_bis=b
END FUNCTION zbrent_c_bis



subroutine func_c_bis ( xg, el0, f_out,METRIC_loc)

  !Use system
  !Use DISC_PARAMETER
  USE TypesDef, ONLY : aom
  Implicit none
  REAL          :: lapse, gp, gm, shift(3),Kex(3,3),g_cov(3,3),g_contr(3,3),phi_tmp
  real          :: xg, el0, f_out
  real          :: Omega_k, el_k
  real          :: x1, x2, x3, pi
  real          :: g00, g03, g33
  REAL          :: x(3)
  Intent(IN) :: xg, el0
  Intent(OUT):: f_out
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

  pi = acos ( -1.0)

  x1 = xg
  x2 = pi / 2.0
  x3 = 0.0d0

  !call system_metric ( x2, x3, x1, (/1,2,3/))

    x = (/x1,x2,x3 /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
    !
  
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)

  Omega_k  = 1.0 / ( xg**1.5 + aom)

  el_k  = - ( g03 + g33 * Omega_k) / ( g00 + g03 * Omega_k)

  f_out = el_k - el0


end subroutine func_c_bis


FUNCTION zbrent_rin_bis( x1, x2, tol, ut_in,el0,METRIC_loc)
  USE recipies_mod, ONLY : nrerror
  IMPLICIT NONE
  REAL, INTENT(IN) :: x1,x2,tol,ut_in,el0
  REAL :: zbrent_rin_bis,fo
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL :: a,b,c,d,e,fa,fb,fc,p,q,r_1,s,tol1,xm
  INTERFACE
     subroutine func_rin_bis (xg, el0, ut_in, f_out,METRIC_loc)
       !Use system
       Implicit none
       real          :: xg,f_out
       real          :: aom,el0,Delta_W
       real          :: g_tt,g_ft,g_ff
       real          :: ut_in,utx,tmp
       real          :: g00, g03, g33
       Intent(IN) :: xg, el0, ut_in
       Intent(OUT):: f_out
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

     end subroutine func_rin_bis
     
  END INTERFACE
  
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

  a=x1
  b=x2
  call func_rin_bis( a, el0, ut_in, fo,METRIC_loc)
  fa = fo
  call func_rin_bis( b, el0, ut_in, fo,METRIC_loc)
  fb = fo
!!$  write(*,*)'a,fa',a,fa
!!$  write(*,*)'b,fb',b,fb
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
       call nrerror('root must be bracketed for zbrent_rin_bis')
  c=b
  fc=fb
  do iter=1,ITMAX
     if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_rin_bis=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r_1=fb/fc
           p=s*(2.0*xm*q*(q-r_1)-(b-a)*(r_1-1.0))
           q=(q-1.0)*(r_1-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     call func_rin_bis ( b, el0, ut_in, fo,METRIC_loc)
     fb = fo
  end do
  call nrerror('zbrent_rin_bis: exceeded maximum iterations')
  zbrent_rin_bis=b
END FUNCTION zbrent_rin_bis

subroutine func_rin_bis (xg, el0, ut_in, f_out,METRIC_loc)
  !Use system
  !Use DISC_PARAMETER
  USE TypesDef, ONLY : aom
  Implicit none
  REAL          :: lapse, gp, gm, shift(3),Kex(3,3),g_cov(3,3),g_contr(3,3),phi_tmp
  real              :: xg
  real          :: f_out
  real          :: el0
  real          :: g_tt,g_ft,g_ff,pi
  real          :: ut_in,utx,tmp
  real          :: g00, g03, g33
  Intent(IN) :: xg, el0, ut_in
  Intent(OUT):: f_out
  REAL          :: x(3)
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

!!$  write(*,*)'from r_in'
!!$  write(*,*)xg,cs,cs2,sn,sn2,M_bh,aom,el0,ut_in,f_out
  pi = acos ( -1.0)

  !call system_metric ( pi / 2.,0.,xg,(/1,2,3/))
  
    x = (/xg,pi / 2.,0. /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
    !
  
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)

  utx   = sqrt((g03**2 - g33*g00) / (g33 + 2.0d0*el0*g03 + el0**2*g00))


  f_out = utx - ut_in


end subroutine func_rin_bis


FUNCTION zbrent_rout_bis ( x1, x2, tol, r_in2, el0,METRIC_loc)
  USE recipies_mod, ONLY : nrerror
  !Use DISC_PARAMETER
  IMPLICIT NONE
  REAL, INTENT(IN) :: x1,x2,tol,r_in2,el0
  REAL :: zbrent_rout_bis,fo
  INTEGER, PARAMETER :: ITMAX=100
  REAL, PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL    :: a,b,c,d,e,fa,fb,fc,p,q,r_1,s,tol1,xm
  INTERFACE
     subroutine func_rout_bis( xg, r_in2, el0, f_out,METRIC_loc)
       !Use system
       Implicit none
       real          :: xg, f_out
       real          :: r_in2, el0
       real          :: ut_in,ut,tmp
       real          :: g00, g03, g33
       Intent(IN) :: xg, r_in2, el0
       Intent(OUT):: f_out
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 
     end subroutine func_rout_bis
  END INTERFACE
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 
  a=x1
  b=x2
  call func_rout_bis ( a, r_in2, el0, fo,METRIC_loc)
  fa = fo
  call func_rout_bis ( b, r_in2, el0, fo,METRIC_loc)
  fb = fo
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
       call nrerror('root must be bracketed for zbrent_rout_bis')
  c=b
  fc=fb
  do iter=1,ITMAX
     if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a
        fc=fa
        d=b-a
        e=d
     end if
     if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     end if
     tol1=2.0*EPS*abs(b)+0.5*tol
     xm=0.5*(c-b)
     if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_rout_bis=b
        RETURN
     end if
     if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa
        if (a == c) then
           p=2.0*xm*s
           q=1.0-s
        else
           q=fa/fc
           r_1=fb/fc
           p=s*(2.0*xm*q*(q-r_1)-(b-a)*(r_1-1.0))
           q=(q-1.0)*(r_1-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q
        p=abs(p)
        if (2.0*p  <  min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        end if
     else
        d=xm
        e=d
     end if
     a=b
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
     call func_rout_bis ( b, r_in2, el0, fo,METRIC_loc)
     fb = fo
  end do
  call nrerror('zbrent_rout_bis: exceeded maximum iterations')
  zbrent_rout_bis=b
END FUNCTION zbrent_rout_bis


subroutine func_rout_bis( xg, r_in2, el0, f_out,METRIC_loc)
  !Use DISC_PARAMETER
  USE TypesDef, ONLY : aom
  Implicit none
  REAL          :: lapse, gp, gm, shift(3),Kex(3,3),g_cov(3,3),g_contr(3,3),phi_tmp
  real          :: f_out, el0
  real          :: r_in2, xg
  real          :: ut_in,ut_1,tmp
  real          :: g00, g03, g33, pi
  Intent(IN) :: xg, r_in2, el0
  Intent(OUT):: f_out
  REAL          :: x(3)
   INTERFACE
      SUBROUTINE METRIC_loc(xc,lapse,gp,gm,shift,Kex,g_cov,g_contr,phi)
#ifdef GRMHD
          USE typesDef, ONLY: aom, Mbh, P_eps
#endif  
          REAL, DIMENSION(3), intent(IN) :: xc
          REAL                           :: lapse,gp,gm, phi 
          REAL,dimension(3)              :: shift
          REAL,dimension(3,3)            :: Kex, g_cov, g_contr 
      END SUBROUTINE METRIC_loc 
   END INTERFACE 

  pi = acos ( -1.0)
  
  !call system_metric ( pi / 2.,0.,r_in,(/1,2,3/))
      x = (/r_in2,pi / 2.,0. /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
    !
  
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)

  ut_in   = sqrt((g03**2 - g33*g00)/               &
       &    (g33 + 2.0d0*el0*g03 + el0**2*g00))


 ! call system_metric ( pi / 2.,0.,xg,(/1,2,3/))
      x = (/xg,pi / 2.,0. /)
    CALL METRIC_loc( x, lapse, gp, gm, shift, Kex, g_cov, g_contr, phi_tmp )
    !
  
  g33 = g_cov (3,3) ! this is g_{phi,phi}
  g03 = g_cov (3,3) * shift (3) + g_cov (1,3) * shift (1)
  g00 = -lapse**2 + g_cov (2,2) * shift (2)**2 + &
       &            g_cov (3,3) * shift (3)**2 + &
       &            g_cov (1,1) * shift (1)**2 + &
       &      2.0 * g_cov (1,3) * shift (3) * shift (1)

  ut_1      = sqrt((g03**2 - g33*g00)/               &
       &    (g33 + 2.0d0*el0*g03 + el0**2*g00))



  f_out = ut_1 - ut_in


    end subroutine func_rout_bis
    
    
    
    
    
    
!
!
!!!$SUBROUTINE mnewt_cons_to_prim(dim,ntrial,x,tolx,tolf)
!!!$  USE nrtype
!!!$!  USE nr, ONLY : lubksb,ludcmp
!!!$  IMPLICIT NONE
!!!$  INTEGER(I4B), INTENT(IN) :: ntrial,dim
!!!$  REAL, INTENT(IN) :: tolx,tolf
!!!$  REAL, DIMENSION(dim), INTENT(INOUT) :: x
!!!$  INTERFACE
!!!$     SUBROUTINE usrfun(x,fvec,fjac)
!!!$       Use gammie_params
!!!$       USE nrtype
!!!$       IMPLICIT NONE
!!!$       REAL, DIMENSION(dim), INTENT(IN) :: x
!!!$       REAL, DIMENSION(dim), INTENT(OUT) :: fvec
!!!$       REAL, DIMENSION(dim,dim), INTENT(OUT) :: fjac
!!!$     END SUBROUTINE usrfun
!!!$  END INTERFACE
!!!$  INTEGER(I4B) :: i
!!!$  INTEGER(I4B), DIMENSION(dim) :: indx
!!!$  REAL :: d
!!!$  REAL, DIMENSION(dim) :: fvec,p
!!!$  REAL, DIMENSION(dim,dim) :: fjac
!!!$  do  i=1,ntrial
!!!$     call usrfun(x,fvec,fjac)
!!!$!     write(*,*)'x,fevc,fjac',x,fvec,fjac
!!!$     if (sum(abs(fvec)) <= tolf) RETURN
!!!$     p=-fvec
!!!$     call ludcmp(fjac,indx,d)
!!!$     call lubksb(fjac,indx,p)
!!!$     x=x+p
!!!$     if (sum(abs(p)) <= tolx) RETURN
!!!$  end do
!!!$END SUBROUTINE mnewt_cons_to_prim

    
#endif    