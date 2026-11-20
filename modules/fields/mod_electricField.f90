module mod_electricField
  use iso_fortran_env, only: real64, int32
  implicit none
  private
  public :: calc_Efield_modular
  public :: calc_Efield_energy_conserving

contains

  subroutine calc_Efield_modular(n, h, phi, E, bcnd)
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    real(real64),   intent(in) :: phi(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    integer :: ix, iy, iz

    ! Parallel zero (60 MB at ITER size; was a single-threaded array
    ! assignment in front of the already-parallel loops below).
    !$omp parallel do private(iz) schedule(static)
    do iz = 0, n(3)+2
      E(:,:,:,iz) = 0.0_real64
    end do
    !$omp end parallel do

    !
    ! Interior points
    !
    !$OMP PARALLEL
    !$OMP DO
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1

          if (bcnd(ix,iy,iz) <= 0) then

            ! Second-order centered differences
            E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
            E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
            E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))

          end if

        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !
    ! Boundary conditions
    !
    !$OMP PARALLEL
    !$OMP DO
    do iz = 1, n(3)+1
      do iy = 1, n(2)+1
        do ix = 1, n(1)+1

          if (bcnd(ix,iy,iz) == -1) goto 10

          !
          ! Neumann BCs
          !
          if (bcnd(ix,iy,iz) == -2) then
            ! YZ plane, LHS only

            E(1,1,iy,iz) = 0.0_real64

            if (iy > 1 .and. iy < n(2)+1 .and. &
                iz > 1 .and. iz < n(3)+1) then
              E(2,1,iy,iz) = -(phi(1,iy+1,iz) - phi(1,iy-1,iz)) / (2.0_real64*h(2))
              E(3,1,iy,iz) = -(phi(1,iy,iz+1) - phi(1,iy,iz-1)) / (2.0_real64*h(3))
              goto 10
            end if
          end if

          !
          ! Periodic BCs
          !
          if (bcnd(ix,iy,iz) == 0 .or. bcnd(ix,iy,iz) == -2) then

            ! iy = 0,1,ny+1,ny+2 planes
            if (iy == 1) then
              E(2,ix,1,iz)       = -(phi(ix,2,iz) - phi(ix,n(2),iz)) / (2.0_real64*h(2))
              E(2,ix,0,iz)       = E(2,ix,n(2),iz)
              E(2,ix,n(2)+1,iz)  = E(2,ix,1,iz)
              E(2,ix,n(2)+2,iz)  = E(2,ix,2,iz)

              E(1,ix,0,iz)       = E(1,ix,n(2),iz)
              E(1,ix,n(2)+2,iz)  = E(1,ix,2,iz)

              E(3,ix,0,iz)       = E(3,ix,n(2),iz)
              E(3,ix,n(2)+2,iz)  = E(3,ix,2,iz)
            end if

            ! iz = 0,1,nz+1,nz+2 planes
            if (iz == 1) then
              E(3,ix,iy,1)       = -(phi(ix,iy,2) - phi(ix,iy,n(3))) / (2.0_real64*h(3))
              E(3,ix,iy,0)       = E(3,ix,iy,n(3))
              E(3,ix,iy,n(3)+1)  = E(3,ix,iy,1)
              E(3,ix,iy,n(3)+2)  = E(3,ix,iy,2)

              E(1,ix,iy,0)       = E(1,ix,iy,n(3))
              E(1,ix,iy,n(3)+2)  = E(1,ix,iy,2)

              E(2,ix,iy,0)       = E(2,ix,iy,n(3))
              E(2,ix,iy,n(3)+2)  = E(2,ix,iy,2)
            end if

            goto 10
          end if

          !
          ! Walls
          !
          if (bcnd(ix,iy,iz) >= 1) then

            ! West wall
            if (bcnd(ix+1,iy,iz) <= 0) then
              E(1,ix,iy,iz) = 2.0_real64*E(1,ix+1,iy,iz) - E(1,ix+2,iy,iz)
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! East wall
            if (bcnd(ix-1,iy,iz) <= 0) then
              E(1,ix,iy,iz) = 2.0_real64*E(1,ix-1,iy,iz) - E(1,ix-2,iy,iz)
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! South wall
            if (bcnd(ix,iy+1,iz) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = 2.0_real64*E(2,ix,iy+1,iz) - E(2,ix,iy+2,iz)
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! North wall
            if (bcnd(ix,iy-1,iz) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = 2.0_real64*E(2,ix,iy-1,iz) - E(2,ix,iy-2,iz)
              E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz-1)) / (2.0_real64*h(3))
            end if

            ! Bottom wall
            if (bcnd(ix,iy,iz+1) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = 2.0_real64*E(3,ix,iy,iz+1) - E(3,ix,iy,iz+2)
            end if

            ! Top wall
            if (bcnd(ix,iy,iz-1) <= 0) then
              E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix-1,iy,iz)) / (2.0_real64*h(1))
              E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy-1,iz)) / (2.0_real64*h(2))
              E(3,ix,iy,iz) = 2.0_real64*E(3,ix,iy,iz-1) - E(3,ix,iy,iz-2)
            end if

          end if

10        continue

        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

  end subroutine calc_Efield_modular


  subroutine calc_Efield_energy_conserving(n, h, phi, E, bcnd, flag_pbc, flag_pbcz)
    ! Face-centered E for the explicit energy-conserving (EC-PIC) scheme
    ! (Powis & Kaganovich, Phys. Plasmas 31, 023901 (2024), Eq. 5):
    !   E(1,ix,iy,iz) = E_x at the +x face of node ix, i.e. physically at
    !                   (x_ix + h(1)/2, y_iy, z_iz) - same array shape/index
    !                   convention as calc_Efield_modular, reinterpreted.
    ! and cyclically for E(2,...)/E(3,...) at +y/+z faces. Paired with
    ! gather_E_energy_conserving (mod_particleMover.f90).
    !
    ! Unlike calc_Efield_modular, no one-sided wall extrapolation is needed:
    ! a face value only ever reads the two nodal potentials either side of
    ! it, and phi is already well-defined at every node touched by a wall
    ! (interior-solved, fixed at its Dirichlet wall value, or - for a
    ! dielectric node - set to its accumulated-charge potential by
    ! apply_dielectric_bc_to_phi, mod_state.f90, before this is called) - x
    ! is never periodic in this codebase (always wall/Neumann), so the
    ! x-direction (E(1,...)'s own values, and every component's x-ghost
    ! planes) never needs special-casing regardless of flag_pbc/flag_pbcz
    ! below. bcnd is accepted (matching calc_Efield_modular's signature so
    ! callers can select between the two uniformly) but not used - a
    ! dielectric node's phi already carries the right value by the time it
    ! gets here, so flag_die needs no branch. flag_nmn is still rejected by
    ! the push_scheme guard in mod_simulation.f90's init(), so there is no
    ! Neumann case to special-case here.
    !
    ! flag_pbc/flag_pbcz (y/z periodic - see mod_generateBoundary.f90;
    ! flag_pbcz=1 always implies flag_pbc=1 too, never z-only) trigger the
    ! ghost fix-up pass below. Node n(2)+1 (resp. n(3)+1) is a real,
    ! independently-solved node the Poisson solver keeps equal to node 1
    ! (the periodic seam) - only iy=0/n(2)+2 (resp. iz=0/n(3)+2) are true
    ! unpopulated ghosts. E_x/E_z's own formulas never difference across y,
    ! so they're already correct at n(2)+1 for free; only their y-ghost
    ! planes need a periodic copy. E_y's own formula *does* difference
    ! across y, so its face values at iy=0 and iy=n(2)+1 (which would
    ! otherwise reach into the true ghosts) need overriding too. Mirrors
    ! calc_Efield_modular's E(2,ix,0,iz)=E(2,ix,n(2),iz)/
    ! E(2,ix,n(2)+1,iz)=E(2,ix,1,iz) ghost-copy pattern, adapted from node
    ! to face indexing. z mirrors this exactly when flag_pbcz==1.
    integer(int32), intent(in) :: n(3)
    real(real64),   intent(in) :: h(3)
    real(real64),   intent(in) :: phi(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32), intent(in) :: flag_pbc, flag_pbcz

    integer :: ix, iy, iz

    E = 0.0_real64

    !$OMP PARALLEL
    !$OMP DO
    do iz = 0, n(3)+2
      do iy = 0, n(2)+2
        do ix = 0, n(1)+1
          E(1,ix,iy,iz) = -(phi(ix+1,iy,iz) - phi(ix,iy,iz)) / h(1)
        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !$OMP PARALLEL
    !$OMP DO
    do iz = 0, n(3)+2
      do iy = 0, n(2)+1
        do ix = 0, n(1)+2
          E(2,ix,iy,iz) = -(phi(ix,iy+1,iz) - phi(ix,iy,iz)) / h(2)
        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !$OMP PARALLEL
    !$OMP DO
    do iz = 0, n(3)+1
      do iy = 0, n(2)+2
        do ix = 0, n(1)+2
          E(3,ix,iy,iz) = -(phi(ix,iy,iz+1) - phi(ix,iy,iz)) / h(3)
        end do
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    !
    ! Periodic ghost fix-up (see header comment) - must run after the three
    ! loops above, since it overwrites specific ghost-plane entries they
    ! computed from unpopulated ghost phi.
    !
    if (flag_pbc == 1_int32) then
      !$OMP PARALLEL
      !$OMP DO
      do iz = 0, n(3)+2
        do ix = 0, n(1)+2
          E(2,ix,0,iz)      = E(2,ix,n(2),iz)
          E(2,ix,n(2)+1,iz) = E(2,ix,1,iz)

          E(1,ix,0,iz)      = E(1,ix,n(2),iz)
          E(1,ix,n(2)+2,iz) = E(1,ix,2,iz)

          E(3,ix,0,iz)      = E(3,ix,n(2),iz)
          E(3,ix,n(2)+2,iz) = E(3,ix,2,iz)
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL
    end if

    if (flag_pbcz == 1_int32) then
      !$OMP PARALLEL
      !$OMP DO
      do iy = 0, n(2)+2
        do ix = 0, n(1)+2
          E(3,ix,iy,0)      = E(3,ix,iy,n(3))
          E(3,ix,iy,n(3)+1) = E(3,ix,iy,1)

          E(1,ix,iy,0)      = E(1,ix,iy,n(3))
          E(1,ix,iy,n(3)+2) = E(1,ix,iy,2)

          E(2,ix,iy,0)      = E(2,ix,iy,n(3))
          E(2,ix,iy,n(3)+2) = E(2,ix,iy,2)
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL
    end if

  end subroutine calc_Efield_energy_conserving

end module mod_electricField
