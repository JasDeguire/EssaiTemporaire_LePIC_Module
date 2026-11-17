# -*- coding: utf-8 -*-
"""
gui_help.py  --  help texts of the LePIC+ input editor (shared by the Tk and the web front ends)
"""
import re


# =========================================================================
#  Help texts   key -> (title, code variable, text)
#  Checked against the LePIC+ modules (mod_readConditions, mod_simParams,
#  mod_generateBoundary, mod_restart, mod_state, mod_coilField ...).
# =========================================================================
HELP = {
    'rname': ('Reactions file', 'rname',
              'Data file (in input_dir/) that defines the species (electrons, ions, neutrals, optionally [H-]) '
              'and every collision / reaction cross section used by the Monte-Carlo collision model. '
              'At most 20 characters, no spaces. Line 1 also holds the number of time steps.'),
    'nsteps': ('Number of time steps', 'nsteps',
               'Length of the run in time steps (second value of line 1, must be > 0). LePIC+ stops after '
               'nsteps steps, not at t_max. A restarted run executes nsteps MORE steps from the backup.'),
    'Te': ('Electron temperature', 'Ti(1)  [eV]',
           'Reference electron temperature. It sets the electron thermal velocity vte = sqrt(2 e Te / me), '
           'hence the time step dt = kt * min(dx,dy,dz) / vte, the Debye length, and the velocities of the '
           'electrons that are loaded or injected. With power-controlled heating the temperature used by '
           'the heating step is recomputed during the run; when particle injection is on, the heating step '
           'keeps this value.'),
    'n0': ('Plasma density', 'n0  [m^-3]',
           'Target plasma density. It fixes the weight of one macro-particle: Nm = n0 * dx*dy*dz / np_cell, '
           'the Debye length and the plasma frequency (both computed with the scaled permittivity).'),
    'ngas': ('Neutral gas density', 'ngas  [m^-3]',
             'Background neutral density used by the collision model and for the weight of neutral '
             'macro-particles.'),
    'np_cell': ('Macro-particles per cell', 'np_cell',
                'Number of macro-particles per cell for each species when the plasma is loaded. Larger = '
                'less noise but slower. Tick "Start without initial plasma" (stored as a negative number) '
                'to begin with an empty box; the particle weight still uses |np_cell|.'),
    'k_eps0': ('Permittivity scaling', 'k_eps0',
               'The code uses eps0 -> k_eps0 * eps0 in Poisson\'s equation. Values > 1 lower the plasma '
               'frequency and enlarge the Debye length (a standard trick to relax the time-step and '
               'grid constraints); results must then be interpreted with that scaling. Cannot be 0.\n'
               'The RF skin depth c/wpe is always computed with the TRUE eps0.\n'
               'Negative value: in LePIC+ the only effect is that a restart reloads the electrode potentials '
               'from Backups/DATA.BAK/Vgrd.bak (the iterative power/current/voltage loop of the legacy code '
               'is not ported).'),
    'tmax': ('Simulated time', 'tmax  [s]',
             'Read but NOT used by LePIC+: the run length is the number of time steps on line 1. '
             'The Summary tab shows nsteps x dt.'),
    'kt': ('Time-step coefficient', 'kt',
           'dt = kt * min(dx,dy,dz) / vte, i.e. the fraction of a cell an average electron crosses per step. '
           'It must also keep wpe*dt small (explicit leap-frog is unstable for wpe*dt > 2, and usually kept '
           'below ~0.2-0.5 for accuracy). kt also selects how often heating is applied '
           '(every 4 steps if kt >= 0.1, 20 if 0.05-0.1, 40 below 0.05).'),
    'nx': ('Cells along x', 'n(1)',
           'Number of grid CELLS along x (not an exponent - the comment in older files is misleading). '
           'dx = xmax / nx. Nodes are numbered 1 ... nx+1. Larger grids need more memory: arrays scale as '
           'nx*ny*nz. If you run with several MPI processes, nz must be divisible by their number.'),
    'ny': ('Cells along y', 'n(2)', 'Number of grid cells along y; dy = ymax / ny.'),
    'nz': ('Cells along z', 'n(3)',
           'Number of grid cells along z; dz = zmax / nz. The Poisson domain is split along z between MPI '
           'processes, so nz must be divisible by the number of MPI processes (the N of mpirun -np N, or the SLURM task count). '
           'The editor reads that number from runModular.sh (next to input_dir/) and checks nz.'),
    'ng': ('Multigrid levels', 'ng',
           'Number of levels of the V-cycle multigrid Poisson solver. LePIC+ always uses this value. '
           'nx, ny and nz must each be divisible by 2^(levels-1) (the coarse grids are n/2^(level-1)); '
           'LePIC+ does not check this at start-up.'),
    'eps': ('Poisson-solver precision', 'eps',
            'Residual threshold at which the multigrid solver stops. Smaller = more accurate but slower.'),
    'omega': ('SOR relaxation parameter', 'omega',
              'Over-relaxation factor of the red-black SOR smoother. Must be in [1, 2); typically 1.4-1.8. '
              'Too close to 2 can make the solver diverge ("PDE solver is diverging").'),
    'xl_pow': ('Heating / injection region', 'xl_pow ... zr_pow  [cm]',
               'Box (in cm) used by two features.\n'
               '- Electron heating: applied to electrons whose x lies between x_min and x_max.\n'
               '- Particle injection (options 1-4): new particles are created uniformly inside the whole box '
               '(x, y and z limits); with the cosine option the profile follows cos() over the box.\n'
               'In round (cylindrical) geometries y_max also defines the heated radius '
               'R = y_max - ymax/2 around the box axis.\n'
               'Limits larger than the box are clipped to it.'),
    'yr_pow': ('Heating box - y_max', 'yr_pow  [cm]',
               'Upper y limit of the heating/injection box. In a cylindrical geometry it defines the heated '
               'radius R = y_max - ymax/2. Tick "annular" (stored as a negative number) to heat only the '
               'electrons OUTSIDE that radius instead of inside it (e.g. RF-driver heating at the wall).'),
    'x_load': ('Plasma loading region', 'x_load  [cm]',
               'The initial plasma is loaded only where x <= this value (clipped to xmax). It also marks the '
               'position of the plasma electrode for the negative-ion / charge-exchange bookkeeping.'),
    'Pabs': ('Absorbed power', 'Pabs  [W]',
             'Power-controlled electron heating: at each heating call the temperature of the re-sampled '
             'Maxwellian is computed from the power balance Te = 2/3 (dEk + Pabs/nu_h) / (Nm e Nh). '
             'A negative value switches this mode off. Needs a heating frequency >= 0. '
             'With an RF antenna (Wrapped, Planar, Coil) this is the target absorbed power: the drive amplitude '
             '(E0, or the coil current I) is rescaled every RF period by sqrt(Pabs / P_RF).'),
    'delta_heat': ('Heating skin depth', 'delta_heat  [cm]',
                   'Optional third value of the Pabs line. 0 or absent = uniform Maxwellian heating in the heating '
                   'region. > 0 = heating weighted by exp((r-R)/delta) from the wall radius R. Tick '
                   '"self-consistent" (stored negative) to recompute delta = c/wpe every 1/f0 (needs f0 > 0 on the '
                   'RF line). Not used with an RF antenna.'),
    'nu_h': ('Heating frequency', 'nu_h  [s^-1]',
             'Rate at which electrons inside the heating region are re-drawn from a Maxwellian. '
             'A negative value turns electron heating off completely. The probability per heating call, '
             'nu_h * (heating interval) * dt, must not exceed 1 - otherwise the code stops '
             '("Warning, nu*dt>1").'),
    'I_inj': ('Injected current', 'I_inj  [A]',
              'Current of electrons injected inside the heating/injection box. For options 1 and 2 each '
              'injection event creates a quasi-neutral set of particles (electrons plus ions, in the '
              'proportions of the reactions file); options 3 and 4 create electrons only. '
              'A negative value switches injection off. Not allowed together with secondary emission > 0. '
              'For option 2 only the sign matters (the amount is whatever replaces the ions lost at the walls).'),
    'opt_inj': ('Injection mode', 'opt_inj',
                '1 = fixed current I_inj, particles uniform in the box.\n'
                '2 = as many injection events as ions lost at the walls since the last call (steady particle '
                'number; I_inj only has to be positive).\n'
                '3 = electron beam along +z (shifted Maxwellian, drift energy = emission temperature TH).\n'
                '4 = electron emission from a cathode plate (half-Maxwellian, temperature TH); the cathode is the '
                'wall with the "emitting electrode label", which must be a flat plate normal to z.\n'
                'Tick "cosine profile" (negative option) for a cos(x)cos(y) injection profile instead of '
                'uniform.'),
    'nB': ('Magnetic-field entries', 'nB',
           'Number of field contributions (0-5). They are summed. Tick "positive ions are not magnetized" '
           '(stored as a negative nB) to apply the field to electrons and negative ions only.'),
    'nBx': ('Magnetic-field map grid', 'n_B(1:3)',
            'Number of cells of the grid on which the field is tabulated. 1 1 1 means a spatially uniform '
            'field. If any value is <= 1 the code resets all three to 1. If the grid differs from the main '
            'grid the field is interpolated.'),
    'B.file': ('Entry type', 'B_file',
               'Not used (s) = the entry is ignored.  Analytical profile (n) = one of the built-in models.  '
               'Field map files (y) = read mag_field_files/<name>Bx.dat, <name>By.dat and <name>Bz.dat '
               '(lines "ix iy iz B") and multiply by the scale factor.'),
    'B.name': ('Model / file prefix', 'B_name',
               'For analytical entries the FIRST LETTER selects the model: g = Gaussian filter field, '
               'P = PG-current filter field (ELISE fit), e = cusp (extraction) magnets, t = Hall-thruster '
               'magnetic circuit. For map files it is the file-name prefix.'),
    'B.scale': ('Scale factor', 'B_scale', 'Multiplies the values read from the map files (map entries only).'),
    'B.B0': ('Field amplitude', 'B0  [G]',
             'Gauss. Gaussian: peak value. PG current: multiplies the built-in 1D profile. Cusp: amplitude. '
             'Hall thruster: peak value at x0; a NEGATIVE value gives the pseudo-1D configuration '
             '(Bx = 0, By uniform along y).'),
    'B.dir': ('Field direction', 'B_info',
              'Bx, By or Bz (Gaussian and PG models only). The cusp model produces Bx and Bz, the thruster model '
              'Bx and By.'),
    'B.dL': ('Width / length scale', 'dL  [cm]',
             'Gaussian: standard deviation sigma of B = B0 exp(-(x-x0)^2 / 2 sigma^2) - use a huge value for a '
             'uniform field. Cusp: distance d between the magnets. Hall thruster: radial distance "a" of the '
             'magnets.'),
    'B.x0': ('Position x0', 'x0  [cm]',
             'Gaussian: centre. Cusp: position of the magnets (the field decays as exp(-pi (x0-x)/d)). '
             'Hall thruster: position of the field maximum.'),
    'B.y0': ('Position y0', 'y0  [cm]', 'Read but not used by the built-in models.'),
    'B.z0': ('Position z0', 'z0  [cm]', 'Cusp model: z position in the middle of the magnets. Unused otherwise.'),
    'ngrid': ('Number of wall labels', 'ngrid',
              'Number of electrode / wall labels. The code reads this many potentials from boundary.inp '
              '(missing ones are 0 V), and every label used in geometry.inp (|label|) must be <= this number. '
              'Labels are numbered 1 ... ngrid; label 0 means "periodic". Tick "no 3D time averages" (stored '
              'negative) to skip the time-averaged 3D potential/density maps.'),
    'nsav': ('Saving frequency', 'nsav  [time steps]',
             'Every nsav steps the code writes its time-averaged output. Tick "keep accumulating averages" '
             '(negative value) so that the averages are not reset after each save.'),
    'nbak': ('Backup frequency', 'nbak  [time steps]',
             'Every nbak steps a restart backup is written (Backups/DATA.BAK/ and DATA.BAK2/ alternately). Tick "do not time-average the saved '
             'data" (negative value) to save instantaneous snapshots instead of averages.'),
    'cnt_plt': ('Number of snapshots', 'cnt_plt',
                'Saves a sequence of cnt_plt profile snapshots, evenly spaced in the time window [t1, t2] '
                '(spacing = |t2-t1| / cnt_plt, rounded to a multiple of the saving period). Negative = no '
                'sequence. Must not be 0 (division by zero).'),
    't1': ('Snapshot window start', 't1  [s]', 'Start of the time window of the snapshot sequence.'),
    't2': ('Snapshot window end', 't2  [s]', 'End of the time window of the snapshot sequence.'),
    'pdf': ('Velocity distribution output', 'ans, ptype, yl ... zr, n_rg',
            'Writes the velocity distribution of one species for particles inside the listed regions: '
            'x ranges per region, one common y and z range. The same regions are also used to tabulate '
            'collision counts. At most 5 regions.'),
    'ptype_pdf': ('Species number', 'ptype_pdf', 'Species index as listed in the reactions file (1 = electrons).'),
    'restart': ('Restart mode', 'ans',
                'New run = start from scratch.  Continue = read Backups/DATA.BAK/ written by a run with the SAME '
                'number of MPI ranks.  Import = read Backups/SAV_DATA/ files written with a DIFFERENT number of '
                'ranks (then fill in the layout of those files below). Note: runModular.sh rewrites this line '
                '(threads, ranks, and "yes 1." for chained segments).'),
    'omp_max': ('OpenMP threads of the backup', 'omp_rank_max',
                'Only used when importing external backups: number of OpenMP threads of the run that wrote the '
                'files (<= 100). It does NOT set how many threads this run uses - that is set in the run script '
                '(OMP_NUM_THREADS).'),
    'mpi_max': ('MPI ranks of the backup', 'mpi_rank_max',
                'Only used when importing external backups: number of MPI ranks of the run that wrote the files; '
                'must be a multiple of the number of ranks used now. It does NOT set the number of processes '
                '(mpirun -np in the run script).'),
    'np_dup': ('Particle duplication / removal', 'np_dup',
               'At restart only. 1 = unchanged. A whole number N >= 2 duplicates every macro-particle N times '
               '(lower noise). A fraction between 0 and 1 randomly keeps that fraction of the particles. '
               'Ignored for a new run. Use whole numbers when above 1.'),
    'holes': ('Grid apertures', 'ans, Lhy, Lhz, nhy, nhz, ind_g',
              'Drills an array of nhy x nhz apertures of size Lhy x Lhz (cm) in the wall with the given label '
              '(a disc when Lhy = Lhz, otherwise a rectangular slit). Apertures are spaced 1.5 x their size '
              'apart and centred in the box. The YZ output plane is taken through this grid. The label '
              'also identifies the plasma-electrode / anode surface used by the negative-ion emission and '
              'the thruster mode.'),
    'jne': ('Negative-ion flux', 'jne  [mA/cm^2]',
            'Current density of negative ions emitted from the surface of the grid selected above (apertures '
            'excluded), at temperature TH. 0 = off. Requires a species called [H-] in the reactions file.'),
    'THm': ('Emission temperature / energy', 'THm  [eV]',
            'Used by three mechanisms: temperature of the emitted negative ions; temperature of secondary '
            'electrons; drift energy (e-beam) or temperature (emissive cathode) of injection options 3 and 4. '
            'Must be > 0 when one of them is used.'),
    'num_grd': ('Grid for negative-ion statistics', 'num_grd',
                'Label of the wall on which the impact characteristics of negative ions are recorded.'),
    'Ca': ('Dielectric capacitance', 'Ca/S  [F/m^2]',
           'Surface capacitance of dielectric walls. The charge accumulated on a dielectric cell is converted '
           'into its potential through this capacitance. Only relevant if the geometry has dielectric surfaces. '
           'A negative value (restart only) tells the code not to reload the dielectric potentials from the '
           'backup; |Ca| is used.'),
    'gam_sec': ('Secondary-electron emission', 'gam_sec',
                'Mean number of electrons released when a positive ion hits the emitting electrode (they leave '
                'with temperature TH, along z). Negative = off. Needs the electrode to be a flat plate normal to '
                'z, and cannot be combined with injected current > 0. When the iterative mode (k_eps0 < 0) is '
                'on and no power/current is set, the electrode voltage is adjusted automatically.'),
    'igrid_sec': ('Emitting electrode label', 'igrid_sec',
                  'Label of the electrode that emits secondary electrons or acts as emissive cathode.'),
    'rf': ('RF mode', 'rf_mode',
           'No = no RF.\n'
           'Electrode = RF voltage on an electrode (not implemented in LePIC+: the code stops).\n'
           'Wrapped = inductive coil wrapped around the driver wall: E_theta peaks at the radius '
           'R = y_max - ymax/2 of the heating region and decays inward over the skin depth.\n'
           'Planar = flat coil under a dielectric window at x = x_min of the heating region: E_theta decays '
           'axially away from the window; radial profile r/R_win inside R_win, R_win/r outside.\n'
           'Coil = same window, but E_theta AND the RF magnetic field come from the real turns listed in '
           'coil.inp (exact loop field, shielded by exp(-(x-x0)/delta)); every species uses the Boris push.\n'
           'With any antenna the skin depth is c/wpe from the mean density, current injection and Maxwellian '
           'heating are switched off, and the amplitude is rescaled every RF period to reach Pabs.'),
    'phi0': ('Voltage amplitude 1', 'phi0_RF  [V]', 'First RF voltage amplitude (Electrode mode, not ported).'),
    'f0': ('Frequency', 'f0_RF  [Hz]',
           'RF frequency. With an antenna the amplitude and the skin depth are updated every 1/f0 '
           '(= 1/(f0 dt) steps), so f0 must be > 0.'),
    'phi1': ('Voltage amplitude 2', 'phi1_RF  [V]', 'Second RF voltage amplitude (Electrode mode, not ported).'),
    'f1': ('Frequency 2', 'f1_RF  [Hz]', 'Second RF frequency (Electrode mode, not ported).'),
    'E0': ('Initial drive amplitude', 'E0_RF  [V/m, or A in Coil mode]',
           'Starting value of the antenna drive, then rescaled every RF period toward Pabs. Wrapped / Planar: '
           'peak field E0 (V/m). Coil: coil current amplitude I (A, peak - divide by sqrt 2 for rms). '
           'Must be > 0.'),
    'R_win': ('Window radius', 'R_win  [cm]',
              'Planar / Coil only. Radius of the disk (r <= R_win, inside the heating slab) over which the mean '
              'electron density that sets delta = c/wpe is taken. Planar: also the radius where the radial '
              'profile peaks. 0 or absent = use y_max - ymax/2 of the heating region.'),
    'coil': ('Coil turns', 'input_dir/coil.inp',
             'Coil mode only. One concentric circular loop per line, coaxial with x: radius a (cm), axial '
             'position x (cm) and relative current w (1 for a series coil). Every turn must be behind the window: '
             'x < max(0, x_min of the heating region), or the code stops (the field is singular on a wire). '
             'The vacuum field per ampere is tabulated at start-up and written to Output/coil_field.dat.'),
    'coulomb': ('Coulomb collisions', 'flag_coulomb',
                'Binary Coulomb collisions between charged particles (yes/no).'),
    'push': ('Push scheme', 'push_scheme',
             'momentum = node-centred E, same trilinear shape for deposit and gather (default). '
             'energy = explicit energy-conserving EC-PIC (Powis & Kaganovich 2024): face-centred E, nearest-grid-'
             'point gather normal to each component; tolerates dx >> lambda_D without grid heating. Not compatible '
             'with the Neumann boundary at x = 0. Must be written in lower case.'),
    'g.xmax': ('Box size x', 'xmax  [cm]',
               'Box length along x in cm. Tick "zero-gradient boundary at x = 0" (stored as a negative xmax) for '
               'a Neumann condition on the left side instead of a wall.'),
    'g.ymax': ('Box size y', 'ymax  [cm]',
               'Box length along y in cm. Tick "round openings along x" (negative ymax) to draw the wall '
               'openings as discs (cylinders along x, centred on the box axis) instead of rectangles; individual '
               'segments can flip this with the "toggle round" option.'),
    'g.zmax': ('Box size z', 'zmax  [cm]',
               'Box length along z in cm. Tick "round objects along z" (negative zmax) to draw discs in the '
               'xy plane (cylinders along z).'),
    'g.seg': ('Wall segment', 'x<,y<,z<, x>,y>,z>, label',
              'Each line draws one region; later segments overwrite earlier ones.\n'
              'WALL WITH OPENING (label > 0): the whole slab x< ... x> becomes a wall at the potential of the '
              'label, EXCEPT the window y<..y> x z<..z> which stays plasma. Typical use: the box walls and '
              'channel/driver walls.\n'
              'SOLID BLOCK (label < 0 in the file): the region x<..x> y<..y> z<..z> is filled with metal/dielectric '
              'of label |label| (plates, electrodes, grids).\n'
              'PERIODIC (label 0): the transverse faces wrap around.\n'
              'Options: dielectric in an x = const plane (x> stored negative, first one is type 3, second 4), '
              'dielectric in a y = const plane (y> negative, type 2), periodic in z (z> negative), toggle round '
              '(y< negative). A flag needs a non-zero value to carry its sign.'),
    'b.V': ('Electrode potential', 'V(label)  [V]',
            'Potential of the wall / electrode with this label. For dielectric surfaces it is the initial '
            'potential. For the emitting electrode in iterative mode it is a starting value (on a restart '
            'with k_eps0 < 0 the code reads it from DATA.BAK/Vgrd.bak instead).'),
}


def help_key(key):
    if key is None:
        return None
    if re.match(r'^B\d+\.', key):
        return 'B.' + key.split('.', 1)[1]
    if re.match(r'^g\.seg\d+', key) or key.startswith('g.f.'):
        return 'g.seg'
    if re.match(r'^b\.V\d+', key):
        return 'b.V'
    if key in ('nBy', 'nBz'):
        return 'nBx'
    if key in ('ny', 'nz', 'nx'):
        return key
    if key in ('yl_pow', 'zl_pow', 'zr_pow', 'xr_pow'):
        return 'xl_pow'
    if key in ('yl_rg', 'yr_rg', 'zl_rg', 'zr_rg', 'n_rg'):
        return 'pdf'
    if key in ('Lhy', 'Lhz', 'nhy', 'nhz', 'ind_g'):
        return 'holes'
    if re.match(r'^coil(\d+\.\d+|\.n)$', key):
        return 'coil'
    return key
