# LePIC+ GUIs: input editor and Output/ viewers

## Without tkinter: the browser version

    python3 GUI/web_gui.py            # then open http://localhost:8765

Needs only numpy and matplotlib. It serves the same three tools (Inputs, Output 2D, Output 3D) as a web page on
127.0.0.1 of the machine where LePIC+ lives:

* **VS Code Remote-SSH**: VS Code forwards the port by itself; click the link it shows (or the Ports panel).
* **Plain ssh**: `ssh -L 8765:localhost:8765 <host>`, then open http://localhost:8765 on your own computer.

Options: `--port N`, `--root <LePIC+ folder>` (only folders inside it can be opened; default: the folder above GUI/).
It edits one case at a time (single user) with the same parser, checks and save logic as the Tk editor
(`pic_inputs.py`); help texts are shared through `gui_help.py`. Stop it with Ctrl-C.

The Tk versions below are unchanged and still work where tkinter is available.

Adapted from the LePIC3D_DD 3.4.1 editor (`GUI_and_tutorial/`) to the input format read by LePIC+
(`modules/config/mod_readConditions.f90`, `modules/geometry/mod_boundary.f90`, `mod_generateBoundary.f90`).

    python3 GUI/pic_input_gui.py                       # opens ./input_dir
    python3 GUI/pic_input_gui.py input_dir/Cyl_coil_RFant

Requires Python 3 with tkinter (RHEL/Alma/Rocky: `sudo dnf install python3-tkinter`; Ubuntu: `sudo apt install python3-tk`).
numpy is optional (geometry preview and 3D view).

`python3 GUI/pic_inputs.py <case folder> ...` runs the headless parser/validator (round trip + checks), no tkinter needed.

## What differs from the LePIC3D_DD editor

* Line 1 holds the reactions file **and the number of time steps**. LePIC+ runs `nsteps` steps; `t_max` is read but unused.
  The Summary tab shows `nsteps x dt`.
* RF modes `No / Electrode / Wrapped / Planar / Coil` (legacy letters D, A, P, C are understood).
  `Electrode` is flagged: LePIC+ stops on it. For the antennas the editor checks what the code checks
  (f0 > 0, initial amplitude > 0, Pabs > 0, antenna radius > 0) and shows the RF update period 1/(f0 dt) and c/wpe.
* Optional values: skin depth `delta` (3rd value of the Pabs line; negative = self-consistent), `E0_RF`
  (the coil current I in amperes for `Coil`) and `R_win` (Planar/Coil). Clearing one removes it from the line.
* **Coil mode**: the turns of `coil.inp` (N, then `a(cm) x(cm) w`) are edited in the Emission & RF tab and saved
  with the case. Each turn must satisfy x < max(0, x_min of the heating region), as in `mod_state.f90`.
* The last two lines are `Coulomb collisions (no/yes)` and `push scheme (momentum/energy)`. A LePIC3D_DD file
  ("Save 3D" line + MC/EC flag) is converted on load and flagged for review.
* A negative number of wall labels is valid (no 3D time averages) and is a check-box.
* Multigrid levels: LePIC+ always uses the file value; each of nx, ny, nz must be divisible by 2^(levels-1).
* MPI rank count is read from `runModular.sh` (`#SBATCH --nodes x --ntasks-per-node`) to check that nz is divisible.
* Restart paths are `Backups/DATA.BAK/` and `Backups/SAV_DATA/`. `runModular.sh` rewrites the restart line itself.
* **Copy to input_dir/**: LePIC+ only reads the files directly in `input_dir/`. This button copies the open case
  (and its own reactions file, if it has one) there; the replaced files go to `input_dir/gui_backups/`.

Save writes only the lines you changed; previous files go to `<case>/gui_backups/`. "Review changes" shows the diff first.

# Output viewers

Adapted from the LePIC3D_DD 3.4.1 viewers (`plot_2D_gui.py`, `plot_3D_param_gui.py`) to `Output/`.
Both need numpy, matplotlib and tkinter. Run them from the run directory (the one holding `Output/` and `input_dir/`):

    python3 GUI/plot_2D_gui.py            # opens ./Output
    python3 GUI/plot_3D_gui.py            # opens ./Output/phi_n_3D.dat, cell size from input_dir/geometry.inp

## plot_2D_gui.py: plane maps and time series

* Maps from `Output/Output_2D/it<step>_<name>_<plane>.mco`. A **Snapshot** box selects the save; "latest"
  follows a running job while auto-update is on, "initial" shows the start-up maps (n<k>, phi1, E).
* Per species: n, T, mean velocities u<k>x/y/z, |flux| j<k> with streamlines (from the angle block of the file),
  source / sink rates. Fields: phi, Ex/Ey/Ez, dr = dh/lambda_De, bcnd; |B|, Bx/By/Bz.
* Time series: Iw<k> (A) and Pw<k> (W) per wall label, Ptot.dat, Mtot.dat, phi_Te_ne_cntr.dat (probe)
  and skin_depth.dat (RF antenna: <ne>, c/wpe, delta, E0 or coil current, P_RF).
* Records repeated after a restart (runModular.sh chains redo the steps after the last backup) are removed on
  reading: the last record of each time is kept.
* Also still reads LePIC2D `DATA/` and 3-D code `DATA/DATA_2D/` folders.

## plot_3D_gui.py: the 3D maps

* Reads `Output/phi_n_3D.dat`: time-averaged phi and the density of species min(2, charged species), i.e. the
  first ion, written when the 3D averages are on (number of wall labels > 0) and overwritten at every save.
* The file has no cell size: it is taken from `geometry.inp` next to the conditions.inp given in the Files box.
  Point that box at the conditions.inp **of the run that wrote the file**.
* E is computed from phi (central differences, masked where n = 0). Slices XY / XZ / YZ and a 3D view.
* Also still reads `DATA/3D_param*.dat` of the 3-D code.
