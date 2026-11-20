module mod_collisionProducts
  use iso_fortran_env, only: int32, real64
  use mod_particles,   only: ParticleSet
  use mod_RNG,         only: ran2

  implicit none
  private
  public :: apply_collision_products

  real(real64), parameter :: QE_ABS = 1.602176634e-19_real64
  real(real64), parameter :: PI_R8  = 3.1415926535897932384626433832795_real64

contains

  subroutine apply_collision_products( &
      part, iproc, ptype, ip, target_type, target_vx, target_vy, target_vz, &
      c_ind, col_info, sig_Eex, mass, Nm, Pcoll, iseed)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: iproc, ptype, ip
    integer(int32), intent(in) :: target_type
    real(real64), intent(in) :: target_vx, target_vy, target_vz
    integer(int32), intent(in) :: c_ind
    integer(int32), intent(in) :: col_info(:,:)
    real(real64), intent(in) :: sig_Eex(:,:), mass(:), Nm(:)
    real(real64), intent(inout) :: Pcoll(:,:)
    integer(int32), intent(inout) :: iseed

    integer(int32) :: n_re, n_by, rt

    n_re = col_info(c_ind,1)
    n_by = col_info(c_ind,2)
    rt   = col_info(c_ind,2+n_re+n_by+1)

    select case (rt)
    case (1_int32, 2_int32, 3_int32, 5_int32)
      call apply_non_charge_exchange_products( &
          part, iproc, ptype, ip, target_type, target_vx, target_vy, target_vz, &
          c_ind, col_info, sig_Eex, mass, Nm, Pcoll, iseed)

    case (4_int32)
      ! Charge exchange: the ion inherits the neutral's sampled velocity.
      ! The resulting fast neutral is released into the background (not tracked).
      call apply_charge_exchange_products( &
          part, iproc, ptype, ip, &
          target_vx, target_vy, target_vz, &
          mass, Nm, Pcoll)

    case default
      error stop "apply_collision_products: unknown reaction type"
    end select

  end subroutine apply_collision_products


  subroutine apply_non_charge_exchange_products( &
      part, iproc, ptype, ip, target_type, target_vx, target_vy, target_vz, &
      c_ind, col_info, sig_Eex, mass, Nm, Pcoll, iseed)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: iproc, ptype, ip
    integer(int32), intent(in) :: target_type
    real(real64), intent(in) :: target_vx, target_vy, target_vz
    integer(int32), intent(in) :: c_ind
    integer(int32), intent(in) :: col_info(:,:)
    real(real64), intent(in) :: sig_Eex(:,:), mass(:), Nm(:)
    real(real64), intent(inout) :: Pcoll(:,:)
    integer(int32), intent(inout) :: iseed

    integer(int32) :: n_re, n_by, rt
    integer(int32) :: i_by, btype, ib
    integer(int32) :: flag_reactant1_kept
    real(real64) :: vx1, vy1, vz1
    real(real64) :: vx2, vy2, vz2
    real(real64) :: vr, mu, sum_mass, vx_cm, vy_cm, vz_cm
    real(real64) :: Eth, dE_heavy, Erel_after
    real(real64) :: th, th_add, phi, cos_th, sin_th
    real(real64) :: cos_th_s, phi_s
    real(real64) :: ex, ey, ez, ex1, ey1, ez1
    real(real64) :: vp, v2old, v2new
    real(real64) :: rnd1, rnd2

    n_re = col_info(c_ind,1)
    n_by = col_info(c_ind,2)
    rt   = col_info(c_ind,2+n_re+n_by+1)

    vx1 = part(ptype,iproc)%pv(4,ip)
    vy1 = part(ptype,iproc)%pv(5,ip)
    vz1 = part(ptype,iproc)%pv(6,ip)

    vx2 = target_vx
    vy2 = target_vy
    vz2 = target_vz

    mu = abs(mass(ptype))*abs(mass(target_type)) / &
         (abs(mass(ptype)) + abs(mass(target_type)))

    vr = sqrt((vx1-vx2)**2 + (vy1-vy2)**2 + (vz1-vz2)**2)

    vx_cm = (abs(mass(ptype))*vx1 + abs(mass(target_type))*vx2) / &
            (abs(mass(ptype)) + abs(mass(target_type)))
    vy_cm = (abs(mass(ptype))*vy1 + abs(mass(target_type))*vy2) / &
            (abs(mass(ptype)) + abs(mass(target_type)))
    vz_cm = (abs(mass(ptype))*vz1 + abs(mass(target_type))*vz2) / &
            (abs(mass(ptype)) + abs(mass(target_type)))

    Eth      = sig_Eex(c_ind,1)
    dE_heavy = sig_Eex(c_ind,2)

    Erel_after = 0.5_real64*mu*vr*vr - Eth*QE_ABS


    if (Erel_after < 0.0_real64) return

    sum_mass = 0.0_real64
    do i_by = 1_int32, n_by
      btype = col_info(c_ind,2+n_re+i_by)
      if (mass(btype) > 0.0_real64) then
        sum_mass = sum_mass + 1.0_real64 / abs(mass(btype))
      end if
    end do

    if (sum_mass <= 0.0_real64) return

    rnd1 = ran2(iseed)
    rnd2 = ran2(iseed)

    th_add = 2.0_real64*PI_R8/real(n_by,real64)
    cos_th = 1.0_real64 - 2.0_real64*rnd1
    th     = acos(cos_th)
    phi    = 2.0_real64*PI_R8*rnd2

    rnd1 = ran2(iseed)
    rnd2 = ran2(iseed)

    cos_th_s = 1.0_real64 - 2.0_real64*rnd1
    phi_s    = 2.0_real64*PI_R8*rnd2

    flag_reactant1_kept = 0_int32

    do i_by = 1_int32, n_by



      btype = col_info(c_ind,2+n_re+i_by)
			

      if (mass(btype) <= 0.0_real64) cycle

      cos_th = cos(th)
      sin_th = sin(th)
      th = th + th_add

      ex = cos_th
      ey = sin_th*sin(phi)
      ez = sin_th*cos(phi)

      call scatter_vector(ex1, ey1, ez1, ex, ey, ez, cos_th_s, phi_s)

      if (i_by == 1_int32 .and. btype == ptype) then
        ib = ip
        flag_reactant1_kept = 1_int32

      else if ((rt == 1_int32 .or. rt == 3_int32) .and. btype == ptype) then
        ib = ip
        flag_reactant1_kept = 1_int32

      else if ((rt == 1_int32 .or. rt == 3_int32) .and. btype == target_type) then
        ! Neutral target is not tracked as a particle in this module.
        ! Skip updating it.
        cycle

      else
        call append_particle_like_position(part(ptype,iproc), part(btype,iproc), ip, ib)
      end if

      v2old = part(btype,iproc)%pv(4,ib)**2 + &
              part(btype,iproc)%pv(5,ib)**2 + &
              part(btype,iproc)%pv(6,ib)**2

      vp = sqrt(2.0_real64*(Erel_after/sum_mass) / abs(mass(btype))**2)

      part(btype,iproc)%pv(4,ib) = vx_cm + vp*ex1
      part(btype,iproc)%pv(5,ib) = vy_cm + vp*ey1
      part(btype,iproc)%pv(6,ib) = vz_cm + vp*ez1

      if (dE_heavy > 0.0_real64 .and. abs(mass(btype)) > 10.0_real64*abs(mass(ptype))) then
        call add_isotropic_energy(part(btype,iproc)%pv(4,ib), &
                                  part(btype,iproc)%pv(5,ib), &
                                  part(btype,iproc)%pv(6,ib), &
                                  dE_heavy, mass(btype), iseed)
      end if

      v2new = part(btype,iproc)%pv(4,ib)**2 + &
              part(btype,iproc)%pv(5,ib)**2 + &
              part(btype,iproc)%pv(6,ib)**2

      Pcoll(btype,iproc) = Pcoll(btype,iproc) + &
          0.5_real64*Nm(btype)*mass(btype)*(v2new - v2old)

    end do

    if (flag_reactant1_kept == 0_int32) then
      call mark_particle_dead(part(ptype,iproc), ip)
    end if

    ! Reactant 2 is the neutral target in this module.
    ! It is represented by ni0/ngas, not by a tracked ParticleSet,
    ! so we do not kill or decrement it here.

    ! Pcoll(ptype,iproc) = Pcoll(ptype,iproc) - Nm(ptype)*Eth*QE_ABS

  end subroutine apply_non_charge_exchange_products


  subroutine append_particle_like_position(src, dst, ip_src, ip_new)
    type(ParticleSet), intent(in) :: src
    type(ParticleSet), intent(inout) :: dst
    integer(int32), intent(in) :: ip_src
    integer(int32), intent(out) :: ip_new

    ip_new = dst%n + 1_int32

    if (ip_new > size(dst%pv, 2)) then
      error stop "append_particle_like_position: particle capacity exceeded"
    end if

    dst%n = ip_new

    dst%pv(1,ip_new) = src%pv(1,ip_src)
    dst%pv(2,ip_new) = src%pv(2,ip_src)
    dst%pv(3,ip_new) = src%pv(3,ip_src)

    dst%pv(4,ip_new) = 0.0_real64
    dst%pv(5,ip_new) = 0.0_real64
    dst%pv(6,ip_new) = 0.0_real64

  end subroutine append_particle_like_position


  subroutine mark_particle_dead(p, ip)
    type(ParticleSet), intent(inout) :: p
    integer(int32), intent(in) :: ip

    ! Prefer your real flag_dead if ParticleSet has one.
    ! For now, remove by swap-with-last.
    if (ip < p%n) then
      p%pv(1,ip)  = p%pv(1,p%n)
      p%pv(2,ip)  = p%pv(2,p%n)
      p%pv(3,ip)  = p%pv(3,p%n)
      p%pv(4,ip) = p%pv(4,p%n)
      p%pv(5,ip) = p%pv(5,p%n)
      p%pv(6,ip) = p%pv(6,p%n)
    end if

    ! Vacated slot: clear its flags (slots past n are kept zeroed).
    if (allocated(p%flag_dead)) p%flag_dead(p%n) = 0
    if (allocated(p%flag_cex))  p%flag_cex(p%n)  = 0_int32

    p%n = p%n - 1_int32

  end subroutine mark_particle_dead


  subroutine add_isotropic_energy(vx, vy, vz, dE_eV, m, iseed)
    real(real64), intent(inout) :: vx, vy, vz
    real(real64), intent(in) :: dE_eV, m
    integer(int32), intent(inout) :: iseed

    real(real64) :: dv, ex, ey, ez

    if (dE_eV <= 0.0_real64) return
    if (m <= 0.0_real64) return

    dv = sqrt(2.0_real64*dE_eV*QE_ABS/abs(m))

    call isotropic_unit_vector(ex, ey, ez, iseed)

    vx = vx + dv*ex
    vy = vy + dv*ey
    vz = vz + dv*ez

  end subroutine add_isotropic_energy


  subroutine isotropic_unit_vector(ex, ey, ez, iseed)
    real(real64), intent(out) :: ex, ey, ez
    integer(int32), intent(inout) :: iseed

    real(real64) :: r1, r2, costh, sinth, phi

    r1 = ran2(iseed)
    r2 = ran2(iseed)

    costh = 2.0_real64*r1 - 1.0_real64
    sinth = sqrt(max(0.0_real64, 1.0_real64 - costh*costh))
    phi = 2.0_real64*PI_R8*r2

    ex = sinth*cos(phi)
    ey = sinth*sin(phi)
    ez = costh

  end subroutine isotropic_unit_vector


  subroutine scatter_vector(vx1, vy1, vz1, vx, vy, vz, costheta, phi)
    real(real64), intent(out) :: vx1, vy1, vz1
    real(real64), intent(in) :: vx, vy, vz, costheta, phi

    real(real64) :: sintheta, cosphi, sinphi, v, vv

    sintheta = sqrt(max(1.0_real64 - costheta**2, 0.0_real64))
    sinphi = sin(phi)
    cosphi = cos(phi)
    v = sqrt(vx**2 + vy**2 + vz**2)

    if (v == 0.0_real64) then
      vx1 = 0.0_real64
      vy1 = 0.0_real64
      vz1 = 0.0_real64
      return
    end if

    if (abs(vy) > abs(vz)) then
      vv = sqrt(vx**2 + vy**2)
      vx1 = vx*costheta + (vy*v*sinphi + vx*vz*cosphi)/vv*sintheta
      vy1 = vy*costheta + (-vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
      vz1 = vz*costheta - vv*cosphi*sintheta
    else
      vv = sqrt(vx**2 + vz**2)
      vx1 = vx*costheta + (vz*v*sinphi - vy*vx*cosphi)/vv*sintheta
      vy1 = vy*costheta + vv*cosphi*sintheta
      vz1 = vz*costheta - (vx*v*sinphi + vy*vz*cosphi)/vv*sintheta
    end if

  end subroutine scatter_vector


  subroutine apply_charge_exchange_products( &
      part, iproc, ptype, ip, &
      target_vx, target_vy, target_vz, &
      mass, Nm, Pcoll)

    type(ParticleSet), intent(inout) :: part(:,:)
    integer(int32), intent(in) :: iproc, ptype, ip
    real(real64), intent(in) :: target_vx, target_vy, target_vz
    real(real64), intent(in) :: mass(:), Nm(:)
    real(real64), intent(inout) :: Pcoll(:,:)

    real(real64) :: v2old, v2new

    v2old = part(ptype,iproc)%pv(4,ip)**2 + &
            part(ptype,iproc)%pv(5,ip)**2 + &
            part(ptype,iproc)%pv(6,ip)**2

    part(ptype,iproc)%pv(4,ip) = target_vx
    part(ptype,iproc)%pv(5,ip) = target_vy
    part(ptype,iproc)%pv(6,ip) = target_vz

    v2new = target_vx**2 + target_vy**2 + target_vz**2

    Pcoll(ptype,iproc) = Pcoll(ptype,iproc) + &
        0.5_real64 * Nm(ptype) * abs(mass(ptype)) * (v2new - v2old)

  end subroutine apply_charge_exchange_products

end module mod_collisionProducts