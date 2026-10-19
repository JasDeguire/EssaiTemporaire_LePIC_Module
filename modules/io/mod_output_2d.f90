module mod_output_2d
  use iso_fortran_env, only: real64, int32

  implicit none
  private

  public :: write_plane_xy_scalar, write_plane_xz_scalar, write_plane_yz_scalar
  public :: write_plane_xy_scalar_2d, write_plane_xz_scalar_2d, write_plane_yz_scalar_2d
  public :: write_density_planes, write_scalar_planes, write_vector_component_planes
  public :: write_plane_xy_flux_2d, write_plane_xz_flux_2d, write_plane_yz_flux_2d

contains

  integer(int32) pure function npoints_1_to_np1(ncell, every) result(npts)
    integer(int32), intent(in) :: ncell, every
    npts = ((ncell + 1_int32 - 1_int32) / every) + 1_int32
  end function npoints_1_to_np1

  subroutine write_plane_xy_scalar(filename, f, n, iz_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), iz_plane, every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: ix, iy
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(2), every)
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(ix,iy,iz_plane), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_scalar


  subroutine write_plane_xz_scalar(filename, f, n, iy_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), iy_plane, every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: ix, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(ix,iy_plane,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_scalar


  subroutine write_plane_yz_scalar(filename, f, n, ix_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), ix_plane, every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: iy, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(2), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(ix_plane,iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_scalar


  subroutine write_density_planes(np, n, ptype, ix_plane, iy_plane, iz_plane, every, prefix)
    ! Writes one species from np(:,:,:,ptype) directly.
    ! Avoids allocating/copying a full 3D temporary array.
    character(len=*), intent(in) :: prefix
    integer(int32),   intent(in) :: n(3), ptype
    integer(int32),   intent(in) :: ix_plane, iy_plane, iz_plane, every
    real(real64),     intent(in) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,*)

    character(len=256) :: fxy, fxz, fyz

    write(fxy,'(a,"_xy.mco")') trim(prefix)
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    write(fyz,'(a,"_yz.mco")') trim(prefix)

    call write_plane_xy_species(fxy, np, n, ptype, iz_plane, every)
    call write_plane_xz_species(fxz, np, n, ptype, iy_plane, every)
    call write_plane_yz_species(fyz, np, n, ptype, ix_plane, every)
  end subroutine write_density_planes


  subroutine write_scalar_planes(f, n, ix_plane, iy_plane, iz_plane, every, prefix)
    character(len=*), intent(in) :: prefix
    integer(int32),   intent(in) :: n(3), ix_plane, iy_plane, iz_plane, every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2,0:n(3)+2)

    character(len=256) :: fxy, fxz, fyz

    write(fxy,'(a,"_xy.mco")') trim(prefix)
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    write(fyz,'(a,"_yz.mco")') trim(prefix)

    call write_plane_xy_scalar(fxy, f, n, iz_plane, every)
    call write_plane_xz_scalar(fxz, f, n, iy_plane, every)
    call write_plane_yz_scalar(fyz, f, n, ix_plane, every)
  end subroutine write_scalar_planes


  subroutine write_vector_component_planes(E, n, comp, ix_plane, iy_plane, iz_plane, every, prefix)
    ! Writes E(comp,:,:,:) directly.
    ! Avoids Ex3/Ey3/Ez3 allocation and full 3D copies in mod_simulation.
    character(len=*), intent(in) :: prefix
    integer(int32),   intent(in) :: n(3), comp
    integer(int32),   intent(in) :: ix_plane, iy_plane, iz_plane, every
    real(real64),     intent(in) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)

    character(len=256) :: fxy, fxz, fyz

    write(fxy,'(a,"_xy.mco")') trim(prefix)
    write(fxz,'(a,"_xz.mco")') trim(prefix)
    write(fyz,'(a,"_yz.mco")') trim(prefix)

    call write_plane_xy_component(fxy, E, n, comp, iz_plane, every)
    call write_plane_xz_component(fxz, E, n, comp, iy_plane, every)
    call write_plane_yz_component(fyz, E, n, comp, ix_plane, every)
  end subroutine write_vector_component_planes


  subroutine write_plane_xy_species(filename, np, n, ptype, iz_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), ptype, iz_plane, every
    real(real64),     intent(in) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,*)
    integer(int32) :: ix, iy
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(2), every)
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( np(ix,iy,iz_plane,ptype), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_species


  subroutine write_plane_xz_species(filename, np, n, ptype, iy_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), ptype, iy_plane, every
    real(real64),     intent(in) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,*)
    integer(int32) :: ix, iz
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( np(ix,iy_plane,iz,ptype), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_species


  subroutine write_plane_yz_species(filename, np, n, ptype, ix_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), ptype, ix_plane, every
    real(real64),     intent(in) :: np(0:n(1)+2,0:n(2)+2,0:n(3)+2,*)
    integer(int32) :: iy, iz
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(2), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( np(ix_plane,iy,iz,ptype), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_species


  subroutine write_plane_xy_component(filename, E, n, comp, iz_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), comp, iz_plane, every
    real(real64),     intent(in) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: ix, iy
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(2), every)
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( E(comp,ix,iy,iz_plane), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_component


  subroutine write_plane_xz_component(filename, E, n, comp, iy_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), comp, iy_plane, every
    real(real64),     intent(in) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: ix, iz
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( E(comp,ix,iy_plane,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_component


  subroutine write_plane_yz_component(filename, E, n, comp, ix_plane, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), comp, ix_plane, every
    real(real64),     intent(in) :: E(3,0:n(1)+2,0:n(2)+2,0:n(3)+2)
    integer(int32) :: iy, iz
    integer :: u
    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(2), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( E(comp,ix_plane,iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_component


  subroutine write_plane_xy_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(2)+2)
    integer(int32) :: ix, iy
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(2), every)
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(ix,iy), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_scalar_2d


  subroutine write_plane_xz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: f(0:n(1)+2,0:n(3)+2)
    integer(int32) :: ix, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(ix,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_scalar_2d


  subroutine write_plane_yz_scalar_2d(filename, f, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: f(0:n(2)+2,0:n(3)+2)
    integer(int32) :: iy, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(2), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( f(iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_scalar_2d

  ! Legacy j<i>_*.mco: a flux-magnitude block followed by a "vector" block
  ! of the in-plane flow angle (legacy write_data.f90: atan2 of the two
  ! in-plane mean-velocity components) - two blocks, one file.
  subroutine write_plane_xy_flux_2d(filename, mag, ang, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: mag(0:n(1)+2,0:n(2)+2), ang(0:n(1)+2,0:n(2)+2)
    integer(int32) :: ix, iy
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(2), every)
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( mag(ix,iy), ix=1,n(1)+1,every )
    end do
    write(u,'(a)') 'vector'
    do iy = n(2)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( ang(ix,iy), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xy_flux_2d


  subroutine write_plane_xz_flux_2d(filename, mag, ang, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: mag(0:n(1)+2,0:n(3)+2), ang(0:n(1)+2,0:n(3)+2)
    integer(int32) :: ix, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(1), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( mag(ix,iz), ix=1,n(1)+1,every )
    end do
    write(u,'(a)') 'vector'
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( ang(ix,iz), ix=1,n(1)+1,every )
    end do
    close(u)
  end subroutine write_plane_xz_flux_2d


  subroutine write_plane_yz_flux_2d(filename, mag, ang, n, every)
    character(len=*), intent(in) :: filename
    integer(int32),   intent(in) :: n(3), every
    real(real64),     intent(in) :: mag(0:n(2)+2,0:n(3)+2), ang(0:n(2)+2,0:n(3)+2)
    integer(int32) :: iy, iz
    integer :: u

    open(newunit=u, file=filename, status='replace', action='write')
    write(u,*) npoints_1_to_np1(n(2), every), npoints_1_to_np1(n(3), every)
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( mag(iy,iz), iy=1,n(2)+1,every )
    end do
    write(u,'(a)') 'vector'
    do iz = n(3)+1, 1, -every
      write(u,'(*(es18.10,1x))') ( ang(iy,iz), iy=1,n(2)+1,every )
    end do
    close(u)
  end subroutine write_plane_yz_flux_2d

end module mod_output_2d
