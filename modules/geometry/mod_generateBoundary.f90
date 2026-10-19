module mod_generateBoundary
  use iso_fortran_env, only: real64
  use mod_utils
  use mod_constants, only: pi
  implicit none
  private
  public :: BoundaryInputs, BoundaryOutputs, generate_boundary, GeometrySegments

  ! ----------------------------
  ! Pack inputs that were globals
  ! ----------------------------
  type :: BoundaryInputs
    integer :: flag_restart = 0
    integer :: flag_convP   = 0
    integer :: flag_grd     = 0

    ! grid holes
    integer :: ind_g = 0
    integer :: nhy   = 0
    integer :: nhz   = 0
    real(real64) :: Lhy = 0.0_real64
    real(real64) :: Lhz = 0.0_real64

    ! secondary emission / injection knobs
    real(real64) :: gam_sec = 0.0_real64
    integer      :: opt_inj = 0
    integer      :: igrid_sec = 0
    integer      :: write_mco = 1
  end type BoundaryInputs

  ! -----------------------------
  ! Pack outputs that were globals
  ! -----------------------------
  type :: BoundaryOutputs
    ! geometry (meters at the end)
    real(real64) :: xmax = 0.0_real64, ymax = 0.0_real64, zmax = 0.0_real64

    ! grid “window” geometry (meters at the end)
    real(real64) :: Lgy  = 0.0_real64, Lgz  = 0.0_real64
    real(real64) :: Sg   = 0.0_real64
    real(real64) :: xg1  = 0.0_real64
    integer      :: ixg  = 0

    ! flags
    integer :: flag_circxh = 0
    integer :: flag_pbc    = 0
    integer :: flag_pbcz   = 0
    integer :: flag_nmn    = 0
    integer :: flag_die    = 0

    ! dielectric bookkeeping
    integer :: ig_die(5) = 0

    ! secondary emission info
    real(real64) :: zg_sec(2) = 0.0_real64
    integer      :: n_cath = 0
    integer      :: dir_sec = 1

    ! output decimation
    integer :: every = 1
  end type BoundaryOutputs



  type :: GeometrySegments
    real(real64) :: xmax = 0.0_real64, ymax = 0.0_real64, zmax = 0.0_real64
    integer      :: nseg = 0
    real(real64), allocatable :: xl(:), yl(:), zl(:), xr(:), yr(:), zr(:)
    integer,      allocatable :: ind(:)
  end type GeometrySegments

contains

  subroutine generate_boundary(u,n,h,bcnd,V,ngrid,dtype,xl_pow,xr_pow,mpi_rank, inp, outp, geo)
    implicit none

    ! ---- core arrays (legacy layout) ----
    integer,      intent(in)    :: n(3), ngrid, mpi_rank
    integer,      intent(inout) :: bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: u   (0:n(1)+2,0:n(2)+2,0:n(3)+2)
    real(real64), intent(inout) :: h(3)
    real(real64), intent(inout) :: V(ngrid)
    integer,      intent(inout) :: dtype(0:ngrid)
    type(GeometrySegments), intent(in) :: geo

    ! heating deposition window along x
    real(real64), intent(in)    :: xl_pow, xr_pow

    type(BoundaryInputs),  intent(in)    :: inp
    type(BoundaryOutputs), intent(inout) :: outp

    ! ---- locals (same as legacy) ----
    integer :: ig,ix,iy,iz,ind,ixl,ixr,iyl,iyr,izl,izr,iyh,izh
    integer :: igl,igr,cnt_hy,cnt_hz,cnt_yz_planes,ind_val
    integer :: flag_circz, flag_circx, flag_circx_tmp
    real(real64) :: xl,xr,yl,yr,zl,zr,x,y,z,yd,zd,R,Sh
    ! -------------------
    ! Initialization block
    ! -------------------
    iyh = 0
    izh = 0
    outp%Lgy  = 0.0_real64
    outp%Lgz  = 0.0_real64
    bcnd = 0
    u    = 0.0_real64

    outp%flag_pbc    = 0
    outp%flag_pbcz   = 0
    outp%flag_nmn    = 0
    outp%flag_die    = 0
    outp%flag_circxh = 0

    cnt_yz_planes = 0
    outp%ig_die    = 0
    outp%zg_sec    = 0.0_real64
    outp%n_cath    = 0
    outp%dir_sec   = 1
    R = 0.0_real64

    outp%every = 1
    dtype = 0
    if (maxval(n) >  512) outp%every = 2
    if (maxval(n) > 1024) outp%every = 4
    if (maxval(n) > 2048) outp%every = 8

    ! -------------------
    ! Read potential values
    ! -------------------
    if (inp%flag_restart == 1 .and. inp%flag_convP == 1) then
      open(41,file='./Backups/DATA.BAK/Vgrd.bak',form='UNFORMATTED')
      read(41) V(1:ngrid)
      close(41)
    end if
90  continue
    close(10)

    ! -------------------
    ! Read geometry.inp
    ! -------------------
    outp%xmax = geo%xmax
    outp%ymax = geo%ymax
    outp%zmax = geo%zmax

    if (outp%xmax < 0.0_real64) outp%flag_nmn = 1
    outp%xmax = abs(outp%xmax)

    flag_circx_tmp = 0
    if (outp%ymax < 0.0_real64) flag_circx_tmp = 1
    outp%ymax = abs(outp%ymax)

    flag_circz = 0
    if (outp%zmax < 0.0_real64) flag_circz = 1
    outp%zmax = abs(outp%zmax)

    h(1) = outp%xmax / n(1)
    h(2) = outp%ymax / n(2)
    h(3) = outp%zmax / n(3)

    ! -------------------
    ! Draw walls and segments
    ! -------------------
    do ig=1,geo%nseg
      xl = geo%xl(ig); yl = geo%yl(ig); zl = geo%zl(ig)
      xr = geo%xr(ig); yr = geo%yr(ig); zr = geo%zr(ig)
      ind = geo%ind(ig)

      dtype(abs(ind)) = 1

      if (flag_circx_tmp == 1) then
        flag_circx = 1
      else
        flag_circx = 0
      end if

      if (yl < 0.0_real64) then
        if (flag_circx_tmp == 1) flag_circx = 0
        if (flag_circx_tmp == 0) flag_circx = 1
      end if
      yl = abs(yl)

      if (flag_circx == 1) then
        R = min((yr-yl)/2.0_real64, (zr-zl)/2.0_real64)
        if (xl <= xl_pow .and. xr >= xr_pow) outp%flag_circxh = 1
      end if

      if (flag_circz == 1) R = min((xr-xl)/2.0_real64, (yr-yl)/2.0_real64)

      if (ind == 0) outp%flag_pbc = 1

      if (xr < 0.0_real64) then
        xr = abs(xr)
        if (cnt_yz_planes == 0) then
          dtype(abs(ind)) = 3
        else
          dtype(abs(ind)) = 4
        end if
        cnt_yz_planes = cnt_yz_planes + 1
        outp%flag_die = 1
        outp%ig_die(dtype(abs(ind))) = int(xl/h(1)) + 1
      end if

      if (yr < 0.0_real64) then
        yr = abs(yr)
        dtype(abs(ind)) = 2
        outp%flag_die = 1
        outp%ig_die(2) = int(yl/h(2)) + 1
      end if

      if (zr < 0.0_real64) then
        zr = abs(zr)
        dtype(abs(ind)) = -dtype(abs(ind))
        outp%flag_pbc = 1
        outp%flag_pbcz = 1
      end if

      if (abs(ind) > ngrid) then
        if (mpi_rank == 0) then
          print*, 'Insufficient number of wall labels found in file boundary.inp'
          print*, 'please correct ...'
        end if
        call stop_calculation
      end if

      ixl = int(xl/h(1)) + 1; if (ixl <= 1) ixl = 0
      ixr = int(xr/h(1)) + 1; if (ixr >= n(1)) ixr = n(1)+2

      iyl = int(yl/h(2)) + 1; if (iyl <= 1) iyl = 0
      iyr = int(yr/h(2)) + 1; if (iyr >= n(2)) iyr = n(2)+2

      izl = int(zl/h(3)) + 1; if (izl <= 1) izl = 0
      izr = int(zr/h(3)) + 1; if (izr >= n(3)) izr = n(3)+2

      ! secondary emission along Oz
      outp%dir_sec = 1
      if ( (inp%gam_sec > 0.0_real64 .or. abs(inp%opt_inj) == 4) .and. inp%igrid_sec == abs(ind) ) then
        if (zr /= zl) then
          print*, 'Warning: secondary particle emission model only along (Oz), please correct...'
          call stop_calculation
        end if

        if (zl == 0.0_real64) then
          outp%zg_sec(1) = zr
          outp%n_cath = outp%n_cath + 1
        end if

        if (zr == outp%zmax) then
          outp%zg_sec(2) = zl
          outp%dir_sec = -1
          outp%n_cath = outp%n_cath + 1
        end if
      end if

      ! grid index where holes are drawn
      if (abs(ind) == inp%ind_g) then
        igl   = ixl
        igr   = ixr
        outp%Lgy = min((iyr-iyl)*h(2), outp%ymax)
        outp%Lgz = min((izr-izl)*h(3), outp%zmax)
        outp%ixg = nint(real(igl+igr,real64)/2.0_real64)
        outp%xg1 = (ixl-1)*h(1)
      end if

      if (ind == 0) then
        ixl = max(ixl,2)
        ixr = min(ixr,n(1))
      end if

      if (ind >= 0) then
        !$OMP PARALLEL
        !$OMP DO
        do iz=0,n(3)+2
          do iy=0,n(2)+2
            do ix=ixl,ixr
              bcnd(ix,iy,iz)=ind
              if(ind>0) u(ix,iy,iz)=V(ind)
            end do
          end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
      end if

      if (ind >= 0) then
        ixl=max(ixl,2); ixr=min(ixr,n(1))
        iyl=max(iyl,2); iyr=min(iyr,n(2))
        izl=max(izl,2); izr=min(izr,n(3))
      end if

      !$OMP PARALLEL
      !$OMP DO
      do iz=izl,izr
        z=(iz-1)*h(3)
        do iy=iyl,iyr
          y=(iy-1)*h(2)
          do ix=ixl,ixr
            x=(ix-1)*h(1)
            if (ind >= 0) then
              if(flag_circz==1 .and. ind>0 .and. ((x-outp%xmax/2.0_real64)**2 + (y-outp%ymax/2.0_real64)**2) > R**2) goto 40
              if(flag_circx==1 .and. ind>0 .and. ((y-outp%ymax/2.0_real64)**2 + (z-outp%zmax/2.0_real64)**2) > R**2) goto 40
              bcnd(ix,iy,iz) = -1
              u(ix,iy,iz)    = 0.0_real64
40            continue
            else
              if(flag_circz==1) then
                if(((x-outp%xmax/2.0_real64)**2 + (y-outp%ymax/2.0_real64)**2) > R**2) goto 45
              end if
              if(flag_circx==1) then
                if(((y-outp%ymax/2.0_real64)**2 + (z-outp%zmax/2.0_real64)**2) > R**2) goto 45
              end if
              bcnd(ix,iy,iz) = abs(ind)
              if(abs(ind)>0) u(ix,iy,iz)=V(abs(ind))
45            continue
            end if
          end do
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL

    end do

50  if (outp%flag_nmn == 1) then
      bcnd(0:1,:,:) = -2
      u(0:1,:,:)    = 0.0_real64
    end if

    if (outp%flag_pbc == 1) then
      !$OMP PARALLEL
      !$OMP DO
      do iy=0,n(2)+2
        do ix=0,n(1)+2
          ig = bcnd(ix,iy,1)
          if (dtype(ig) < 0) then
            if (bcnd(ix,iy,2) == -1) then
              bcnd(ix,iy,0:1) = 0
              bcnd(ix,iy,n(3)+1:n(3)+2) = 0
            end if
          end if
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP END PARALLEL
    end if

    dtype = abs(dtype)

    ! ---- holes ----
    if (inp%flag_grd == 0) goto 100

    iyh = nint(inp%Lhy/2.0_real64/h(2)) + 1
    izh = nint(inp%Lhz/2.0_real64/h(3)) + 1

    if (inp%Lhy > outp%ymax) stop "Lhy > ymax"
    if (inp%Lhz > outp%zmax) stop "Lhz > zmax"

    if (inp%Lhy == inp%Lhz) then
      R  = inp%Lhy/2.0_real64
      Sh = pi*R**2
    else
      Sh = inp%Lhy*inp%Lhz
    end if

    do cnt_hz=1,inp%nhz
      zd  = -(inp%nhz-1)*3.0_real64*inp%Lhz/4.0_real64 + (cnt_hz-1)*3.0_real64*inp%Lhz/2.0_real64 + outp%zmax/2.0_real64
      izl = nint(zd/h(3)) + 1 - izh
      izr = izl + 2*izh
      if (izl < 1) izl = 1
      if (izr > n(3)+1) izr = n(3)+1

      do cnt_hy=1,inp%nhy
        yd  = -(inp%nhy-1)*3.0_real64*inp%Lhy/4.0_real64 + (cnt_hy-1)*3.0_real64*inp%Lhy/2.0_real64 + outp%ymax/2.0_real64
        iyl = nint(yd/h(2)) + 1 - iyh
        iyr = iyl + 2*iyh
        if (iyl < 1) iyl = 1
        if (iyr > n(2)+1) iyr = n(2)+1

        do iz=izl,izr
          do iy=iyl,iyr
            y=(iy-1)*h(2)
            z=(iz-1)*h(3)
            ind_val = -1
            if (iy==1 .or. iy==n(2)+1 .or. iz==1 .or. iz==n(3)+1) ind_val = 0
            do ix=igl,igr
              if (inp%Lhy == inp%Lhz) then
                if (((y-yd)**2 + (z-zd)**2) <= R**2) then
                  bcnd(ix,iy,iz) = ind_val
                  u(ix,iy,iz)    = 0.0_real64
                end if
              else
                bcnd(ix,iy,iz) = ind_val
                u(ix,iy,iz)    = 0.0_real64
              end if
            end do
          end do
        end do

      end do
    end do

100 continue

    ! Convert units to meters
    h        = h * 1.0e-2_real64
    outp%xmax = outp%xmax * 1.0e-2_real64
    outp%ymax = outp%ymax * 1.0e-2_real64
    outp%zmax = outp%zmax * 1.0e-2_real64

    outp%xg1 = outp%xg1 * 1.0e-2_real64
    outp%Lgy = outp%Lgy * 1.0e-2_real64
    outp%Lgz = outp%Lgz * 1.0e-2_real64

    if (inp%flag_grd == 1) then
        outp%Sg = outp%Lgy*outp%Lgz - (inp%nhy-1)*(inp%nhz-1)*Sh*1.0e-4_real64
    else
        outp%Sg = outp%Lgy*outp%Lgz
    end if


    outp%zg_sec = outp%zg_sec * 1.0e-2_real64

    close(10)

    ! write 2D files
    if (mpi_rank == 0 .and. inp%write_mco == 1) then
      izl = n(3)/2 + 1
      open(12,file='./Output/Output_2D/bcnd_xy.mco')
      write(12,*) n(1)/outp%every, n(2)/outp%every
      do iy=n(2)+1,1,-1*outp%every
        write(12,101) (bcnd(ix,iy,izl), ix=1,n(1)+1,outp%every)
101     format(800(i3,1x))
      end do
      close(12)

      open(12,file='./Output/Output_2D/bcnd_xz.mco')
      write(12,*) n(1)/outp%every, n(3)/outp%every
      do iz=n(3)+1,1,-1*outp%every
        write(12,101) (bcnd(ix,n(2)/2+1,iz), ix=1,n(1)+1,outp%every)
      end do
      close(12)

      ix = n(1)/2 + 1
      if (inp%flag_grd == 1) ix = outp%ixg
      open(12,file='./Output/Output_2D/bcnd_yz.mco')
      write(12,*) n(2)/outp%every, n(3)/outp%every
      do iz=n(3)+1,1,-1*outp%every
        write(12,101) (bcnd(ix,iy,iz), iy=1,n(2)+1,outp%every)
      end do
      close(12)


      izl = n(3)/2 + 1
      open(12,file='./Output/Output_2D/phi_xy.mco')
      write(12,*) n(1)/outp%every, n(2)/outp%every
102   format(800(f6.2,1x))
      do iy=n(2)+1,1,-1*outp%every
        write(12,102) (u(ix,iy,izl), ix=1,n(1)+1,outp%every)
        
      end do
      close(12)

      open(12,file='./Output/Output_2D/phi_xz.mco')
      write(12,*) n(1)/outp%every, n(3)/outp%every
      do iz=n(3)+1,1,-1*outp%every
        write(12,102) (u(ix,n(2)/2+1,iz), ix=1,n(1)+1,outp%every)
      end do
      close(12)

      ix = n(1)/2 + 1
      if (inp%flag_grd == 1) ix = outp%ixg
      open(12,file='./Output/Output_2D/phi_yz.mco')
      write(12,*) n(2)/outp%every, n(3)/outp%every
      do iz=n(3)+1,1,-1*outp%every
        write(12,102) (u(ix,iy,iz), iy=1,n(2)+1,outp%every)
      end do
      close(12)

    end if

  end subroutine generate_boundary

end module mod_generateBoundary
