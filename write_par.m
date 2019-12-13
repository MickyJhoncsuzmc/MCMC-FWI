function write_par(simulation_dir, nx, nz, nt, dt, xmin, xmax, zmin, zmax, nrec, rec_x, rec_z, sx, sz, f0, rho, v_p, v_s)

%% Source file
fid = fopen([simulation_dir,'/Standard_run/DATA/SOURCE'], 'w');
fprintf(fid, '#-----------------------------------------------------------------------------\n'); 
fprintf(fid, 'source_surf                     = .false.         # source inside the medium, or source automatically moved exactly at the surface by the solver\n'); 
fprintf(fid, 'xs                              = %f             # source location x in meters\n',sx(1)); 
fprintf(fid, 'zs                              = %f             # source location z in meters (zs is ignored if source_surf is set to true, it is replaced with the topography height)\n',sz(1)); 
fprintf(fid, 'source_type                     = 1              \n'); 
fprintf(fid, 'time_function_type              = 1              \n'); 
fprintf(fid, 'name_of_source_file             = /home/haipeng/Software/specfem2d/DATA/source1.dat                            \n'); 
fprintf(fid, 'burst_band_width                = 0.             # Only for option 9 : band width of the burst                 \n'); 
fprintf(fid, 'f0                              = %f             # dominant source frequency (Hz) if not Dirac or Heaviside    \n', f0); 
fprintf(fid, 'tshift                          = 0.0            # time shift when multi sources (if one source, must be zero) \n'); 
fprintf(fid, 'anglesource                     = 0.             # angle of the source (for a force only); for a plane wave, this is the incidence angle; for moment tensor sources this is unused\n'); 
fprintf(fid, 'Mxx                             = 1.             # Mxx component (for a moment tensor source only)   \n'); 
fprintf(fid, 'Mzz                             = 1.             # Mzz component (for a moment tensor source only)   \n'); 
fprintf(fid, 'Mxz                             = 0.             # Mxz component (for a moment tensor source only)   \n'); 
fprintf(fid, 'factor                          = 1.d10          # amplification factor                              \n'); 
fclose(fid);

%% Model file
fid = fopen([simulation_dir,'/Standard_run/DATA/topo_model.dat'], 'w');                                   
fprintf(fid, '# number of interfaces                                                                 \n'); 
fprintf(fid, ' 2                                                                                     \n'); 
fprintf(fid, '#                                                                                      \n'); 
fprintf(fid, '# for each interface below, we give the number of points and then x,z for each point   \n'); 
fprintf(fid, '#                                                                                      \n'); 
fprintf(fid, '# interface number 1 (bottom of the mesh)                                              \n'); 
fprintf(fid, ' 2                                                                                     \n'); 
fprintf(fid, ' %f %f                                                                                 \n', xmin, zmin); 
fprintf(fid, ' %f %f                                                                                 \n', xmax, zmin); 
fprintf(fid, '# interface number 2 (topography, top of the mesh)                                     \n'); 
fprintf(fid, ' 2                                                                                     \n'); 
fprintf(fid, ' %f %f                                                                                 \n', xmin, zmax);  
fprintf(fid, ' %f %f                                                                                 \n', xmax, zmax); 
fprintf(fid, '#                                                                                      \n'); 
fprintf(fid, '# for each layer, we give the number of spectral elements in the vertical direction    \n'); 
fprintf(fid, '#                                                                                      \n'); 
fprintf(fid, '# layer number 1 (bottom layer)                                                        \n'); 
fprintf(fid, '%d                                                                                     \n', nz); 
fclose(fid);



%% Par_file
fid = fopen([simulation_dir,'/Standard_run/DATA/Par_file'], 'w');
fprintf(fid, '#-----------------------------------------------------------------------------\n'); 
fprintf(fid, '#-----------------------------------------------------------------------------\n'); 
fprintf(fid, '#                                                                             \n'); 
fprintf(fid, '# simulation input parameters                                                 \n'); 
fprintf(fid, '#                                                                             \n'); 
fprintf(fid, '#-----------------------------------------------------------------------------\n'); 
fprintf(fid, '                                                                              \n'); 
fprintf(fid, '# title of job                                                                \n'); 
fprintf(fid, 'title                           = FWI Simulation                              \n'); 

fprintf(fid, '# forward or adjoint simulation                                               \n');
fprintf(fid, '# 1 = forward, 2 = adjoint, 3 = both simultaneously                           \n');
fprintf(fid, '# note: 2 is purposely UNUSED (for compatibility with the numbering of our 3D codes)\n');
fprintf(fid, 'SIMULATION_TYPE                 = 1                                           \n');
fprintf(fid, '# 0 = regular wave propagation simulation, 1/2/3 = noise simulation           \n');
fprintf(fid, 'NOISE_TOMOGRAPHY                = 0                                           \n');
fprintf(fid, '# save the last frame, needed for adjoint simulation                          \n');
fprintf(fid, 'SAVE_FORWARD                    = .false.                                     \n');
fprintf(fid, '# parameters concerning partitioning\n');
fprintf(fid, 'NPROC                           = 1              # number of processes        \n');
fprintf(fid, 'partitioning_method             = 3              # SCOTCH = 3, ascending order (very bad idea) = 1\n');

fprintf(fid, '# number of control nodes per element (4 or 9)\n');
fprintf(fid, 'ngnod                           = 9                                           \n');

fprintf(fid, '# time step parameters                                                        \n');
fprintf(fid, '# total number of time steps                                                  \n');
fprintf(fid, 'NSTEP                           = %d                                          \n', nt);
fprintf(fid, '# duration of a time step (see section "How to choose the time step" of the manual for how to do this)\n');
fprintf(fid, 'DT                              = %f                                          \n', dt);

fprintf(fid, '# time stepping                                                               \n');
fprintf(fid, '# 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta), 3 = classical RK4 4th-order 4-stage Runge-Kutta\n');
fprintf(fid, 'time_stepping_scheme            = 1                                           \n');

fprintf(fid, '# axisymmetric (2.5D) or Cartesian planar (2D) simulation                     \n');
fprintf(fid, 'AXISYM                          = .false.                                     \n');

fprintf(fid, '# set the type of calculation (P-SV or SH/membrane waves)                     \n');
fprintf(fid, 'P_SV                            = .true.                                      \n');

fprintf(fid, '# set to true to use GPUs\n');
fprintf(fid, 'GPU_MODE                        = .true.                                     \n');

fprintf(fid, '# creates/reads a binary database that allows to skip all time consuming setup steps in initialization\n');
fprintf(fid, '# 0 = does not read/create database    \n');
fprintf(fid, '# 1 = creates database                 \n');
fprintf(fid, '# 2 = reads database                   \n');
fprintf(fid, 'setup_with_binary_database      = 0    \n');
 
fprintf(fid, '# available models                                        \n');
fprintf(fid, '#   default       - define model using nbmodels below     \n');
fprintf(fid, '#   ascii         - read model from ascii database file   \n');
fprintf(fid, '#   binary        - read model from binary databse file   \n');
fprintf(fid, '#   binary_voigt  - read Voigt model from binary database file \n');
fprintf(fid, '#   external      - define model using define_external_model subroutine \n');
fprintf(fid, '#   gll           - read GLL model from binary database file \n');
fprintf(fid, '#   legacy        - read model from model_velocity.dat_input \n');
fprintf(fid, 'MODEL                           = default                    \n');

fprintf(fid, '# Output the model with the requested type, does not save if turn to default or .false.\n');
fprintf(fid, '# (available output formats: ascii,binary,gll,legacy)                                  \n');
fprintf(fid, 'SAVE_MODEL                      = default                                              \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# attenuation                                                                 \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# attenuation parameters\n');
fprintf(fid, 'ATTENUATION_VISCOELASTIC        = .false.        \n');
fprintf(fid, 'ATTENUATION_VISCOACOUSTIC       = .false.        \n');

fprintf(fid, '# for viscoelastic or viscoacoustic attenuation\n');
fprintf(fid, 'N_SLS                           = 3              \n');
fprintf(fid, 'ATTENUATION_f0_REFERENCE        = 5.196          \n');
fprintf(fid, 'READ_VELOCITIES_AT_f0           = .false.        \n');
fprintf(fid, 'USE_SOLVOPT                     = .false.        \n');

fprintf(fid, '# for poroelastic attenuation                    \n');
fprintf(fid, 'ATTENUATION_PORO_FLUID_PART     = .false.        \n');
fprintf(fid, 'Q0_poroelastic                  = 1              \n');
fprintf(fid, 'freq0_poroelastic               = 10             \n');

fprintf(fid, '# to undo attenuation and/or PMLs for sensitivity kernel calculations or forward runs with SAVE_FORWARD                             \n');
fprintf(fid, '# use the flag below. It performs undoing of attenuation and/or of PMLs in an exact way for sensitivity kernel calculations         \n');
fprintf(fid, '# but requires disk space for temporary storage, and uses a significant amount of memory used as buffers for temporary storage.     \n');
fprintf(fid, '# When that option is on the second parameter indicates how often the code dumps restart files to disk (if in doubt, use something between 100 and 1000).\n');
fprintf(fid, 'UNDO_ATTENUATION_AND_OR_PML     = .false.   \n');
fprintf(fid, 'NT_DUMP_ATTENUATION             = 500       \n');

fprintf(fid, '# Instead of reconstructing the forward wavefield, this option reads it from the disk using asynchronous I/O.   \n');
fprintf(fid, '# Outperforms conventional mode using a value of NSTEP_BETWEEN_COMPUTE_KERNELS high enough.                     \n');
fprintf(fid, 'NO_BACKWARD_RECONSTRUCTION      = .false.                                                                       \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# sources                                                                     \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# source parameters\n');
fprintf(fid, 'NSOURCES                        = 1              \n');
fprintf(fid, 'force_normal_to_surface         = .false.        \n');

fprintf(fid, '# use an existing initial wave field as source or start from zero (medium initially at rest)\n');
fprintf(fid, 'initialfield                    = .false.        \n');
fprintf(fid, 'add_Bielak_conditions_bottom    = .false.        \n');
fprintf(fid, 'add_Bielak_conditions_right     = .false.        \n');
fprintf(fid, 'add_Bielak_conditions_top       = .false.        \n');
fprintf(fid, 'add_Bielak_conditions_left      = .false.        \n');

fprintf(fid, '# acoustic forcing\n');
fprintf(fid, 'ACOUSTIC_FORCING                = .false.        \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# receivers                                                                   \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# receiver set parameters for recording stations (i.e. recording points)\n');
fprintf(fid, 'seismotype                      = 4              # record 1=displ 2=veloc 3=accel 4=pressure 5=curl of displ 6=the fluid potential  \n');

fprintf(fid, '# subsampling of the seismograms to create smaller files (but less accurately sampled in time)                                      \n');
fprintf(fid, 'subsamp_seismos                 = 1                                                                                                 \n');

fprintf(fid, '# so far, this option can only be used if all the receivers are in acoustic elements                                                \n');
fprintf(fid, 'USE_TRICK_FOR_BETTER_PRESSURE   = .false.                                                                                           \n');

fprintf(fid, '# every how many time steps we save the seismograms                                                                                 \n');
fprintf(fid, '# (costly, do not use a very small value; if you use a very large value that is larger than the total number                        \n');
fprintf(fid, '#  of time steps of the run, the seismograms will automatically be saved once at the end of the run anyway)                         \n');
fprintf(fid, 'NSTEP_BETWEEN_OUTPUT_SEISMOS    = 1000                                                                                              \n');

fprintf(fid, '# use this t0 as earliest starting time rather than the automatically calculated one                                                \n');
fprintf(fid, 'USER_T0                         = 0.0d0                                                                                             \n');

fprintf(fid, '# seismogram formats                             \n');
fprintf(fid, 'save_ASCII_seismograms          = .true.         \n');
fprintf(fid, 'save_binary_seismograms_single  = .true.         \n');
fprintf(fid, 'save_binary_seismograms_double  = .false.        \n');
fprintf(fid, 'SU_FORMAT                       = .false.        \n');
fprintf(fid, '# use an existing STATION file found in ./DATA or create a new one from the receiver positions below in this Par_file   \n');
fprintf(fid, 'use_existing_STATIONS           = .false.        \n');

fprintf(fid, '# number of receiver sets (i.e. number of receiver lines to create below)  \n');
fprintf(fid, 'nreceiversets                   = 1              \n');

fprintf(fid, '# orientation\n');
fprintf(fid, 'anglerec                        = 0.d0           # angle to rotate components at receivers                               \n');
fprintf(fid, 'rec_normal_to_surface           = .false.        # base anglerec normal to surface (external mesh and curve file needed) \n');

fprintf(fid, '# first receiver set (repeat these 6 lines and adjust nreceiversets accordingly)       \n');
fprintf(fid, 'nrec                            = %d             # number of receivers                 \n',nrec);
fprintf(fid, 'xdeb                            = %d             # first receiver x in meters          \n',rec_x(1));
fprintf(fid, 'zdeb                            = %d             # first receiver z in meters          \n',rec_z(1));
fprintf(fid, 'xfin                            = %d             # last receiver x in meters (ignored if only one receiver)   \n', rec_x(end));
fprintf(fid, 'zfin                            = %d             # last receiver z in meters (ignored if only one receiver)   \n', rec_z(end));
fprintf(fid, 'record_at_surface_same_vertical = .false.         # receivers inside the medium or at the surface (z values are ignored if this is set to true, they are replaced with the topography height)\n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# adjoint kernel outputs                                                      \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# save sensitivity kernels in ASCII format (much bigger files, but compatible with current GMT scripts) or in binary format  \n');
fprintf(fid, 'save_ASCII_kernels              = .false.                                                                                    \n');

fprintf(fid, '# since the accuracy of kernel integration may not need to respect the CFL, this option permits to save computing time, and memory with UNDO_ATTENUATION_AND_OR_PML mode \n');
fprintf(fid, 'NSTEP_BETWEEN_COMPUTE_KERNELS   = 1                                           \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# boundary conditions                                                         \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# Perfectly Matched Layer (PML) boundaries       \n');
fprintf(fid, '# absorbing boundary active or not               \n');
fprintf(fid, 'PML_BOUNDARY_CONDITIONS         = .false.         \n');
fprintf(fid, 'NELEM_PML_THICKNESS             = 3              \n');
fprintf(fid, 'ROTATE_PML_ACTIVATE             = .false.        \n');
fprintf(fid, 'ROTATE_PML_ANGLE                = 30.            \n');
fprintf(fid, '# change the four parameters below only if you know what you are doing; they change the damping profiles inside the PMLs       \n');
fprintf(fid, 'K_MIN_PML                       = 1.0d0 # from Gedney page 8.11                                                                \n');
fprintf(fid, 'K_MAX_PML                       = 1.0d0       \n');
fprintf(fid, 'damping_change_factor_acoustic  = 0.5d0       \n');
fprintf(fid, 'damping_change_factor_elastic   = 1.0d0       \n');
fprintf(fid, 'PML_PARAMETER_ADJUSTMENT        = .false.     \n');

fprintf(fid, '# Stacey ABC       \n');
fprintf(fid, 'STACEY_ABSORBING_CONDITIONS     = .true.       \n');

fprintf(fid, '# periodic boundaries                           \n');
fprintf(fid, 'ADD_PERIODIC_CONDITIONS         = .false.       \n');
fprintf(fid, 'PERIODIC_HORIZ_DIST             = 4000.d0       \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# velocity and density models                                                 \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, 'nbmodels                        = %d                                          \n',nx * nz);
fprintf(fid, '# available material types (see user manual for more information)             \n');
fprintf(fid, '#   acoustic:    model_number 1 rho Vp 0  0 0 QKappa 9999 0 0 0 0 0 0 (for QKappa use 9999 to ignore it)             \n');
fprintf(fid, '#   elastic:     model_number 1 rho Vp Vs 0 0 QKappa Qmu 0 0 0 0 0 0 (for QKappa and Qmu use 9999 to ignore them)    \n');
kk = 0;
rho = flipud(rho);
v_p = flipud(v_p);
v_s = flipud(v_s);

for ii = 1 : nx
    for jj = 1 : nz
        kk = kk + 1;
        fprintf(fid, '%d 1 %f %f %f 0 0 9999 9999 0 0 0 0 0 0    \n', kk, rho(jj, ii), v_p(jj,ii), v_s(jj,ii));
    end
end


fprintf(fid, '# external tomography file                                \n');
fprintf(fid, 'TOMOGRAPHY_FILE                 = ./DATA/tomo_file.xyz    \n');

fprintf(fid, '# use an external mesh created by an external meshing tool or use the internal mesher    \n');
fprintf(fid, 'read_external_mesh              = .false.                                                \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# PARAMETERS FOR EXTERNAL MESHING                                             \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# data concerning mesh, when generated using third-party app (more info in README)                                    \n');
fprintf(fid, '# (see also absorbing_conditions above)                                                                               \n');
fprintf(fid, 'mesh_file                       = ./DATA/mesh_file          # file containing the mesh                                \n');
fprintf(fid, 'nodes_coords_file               = ./DATA/nodes_coords_file  # file containing the nodes coordinates                   \n');
fprintf(fid, 'materials_file                  = ./DATA/materials_file     # file containing the material number for each element    \n');
fprintf(fid, 'free_surface_file               = ./DATA/free_surface_file  # file containing the free surface                        \n');
fprintf(fid, 'axial_elements_file             = ./DATA/axial_elements_file   # file containing the axial elements if AXISYM is true \n');
fprintf(fid, 'absorbing_surface_file          = ./DATA/absorbing_surface_file   # file containing the absorbing surface             \n');
fprintf(fid, 'acoustic_forcing_surface_file   = ./DATA/MSH/Surf_acforcing_Bottom_enforcing_mesh   # file containing the acoustic forcing surface     \n');
fprintf(fid, 'absorbing_cpml_file             = ./DATA/absorbing_cpml_file   # file containing the CPML element numbers                              \n');
fprintf(fid, 'tangential_detection_curve_file = ./DATA/courbe_eros_nodes  # file containing the curve delimiting the velocity model                  \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# PARAMETERS FOR INTERNAL MESHING                                             \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# file containing interfaces for internal mesh                                               \n');
fprintf(fid, 'interfacesfile                  = ./topo_model.dat                                           \n');
fprintf(fid, '# geometry of the model (origin lower-left corner = 0,0) and mesh description                \n');
fprintf(fid, 'xmin                            = %f             # abscissa of left side of the model        \n', xmin);
fprintf(fid, 'xmax                            = %f             # abscissa of right side of the model       \n', xmax);
fprintf(fid, 'nx                              = %d             # number of elements along X                \n', nx);

fprintf(fid, '# absorbing boundary parameters (see absorbing_conditions above)      \n');
fprintf(fid, 'absorbbottom                    = .true.      \n');
fprintf(fid, 'absorbright                     = .true.      \n');
fprintf(fid, 'absorbtop                       = .false.     \n');
fprintf(fid, 'absorbleft                      = .true.      \n');

fprintf(fid, '# define the different regions of the model in the (nx,nz) spectral-element mesh      \n');
fprintf(fid, 'nbregions                       = %d                                                  \n', nx * nz);
fprintf(fid, '# format of each line: nxmin nxmax nzmin nzmax material_number                        \n');
kk = 0;
for ii = 1 : nx
    for jj = 1 : nz
        kk = kk + 1;
        fprintf(fid, '%d  %d  %d  %d  %d \n', ii, ii, jj, jj, kk);
    end
end

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# display parameters                                                          \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');

fprintf(fid, '# every how many time steps we display information about the simulation (costly, do not use a very small value)                        \n');
fprintf(fid, 'NSTEP_BETWEEN_OUTPUT_INFO       = 1000                                                                                                 \n');

fprintf(fid, '# meshing output                                                                                                                       \n');
fprintf(fid, 'output_grid_Gnuplot             = .false.        # generate a GNUPLOT file containing the grid, and a script to plot it                \n');
fprintf(fid, 'output_grid_ASCII               = .false.        # dump the grid in an ASCII text file consisting of a set of X,Y,Z points or not      \n');

fprintf(fid, '# to plot total energy curves, for instance to monitor how CPML absorbing layers behave;  \n');
fprintf(fid, '# should be turned OFF in most cases because a bit expensive                              \n');
fprintf(fid, 'OUTPUT_ENERGY                   = .false.                                                 \n');

fprintf(fid, '# every how many time steps we compute energy (which is a bit expensive to compute)       \n');
fprintf(fid, 'NTSTEP_BETWEEN_OUTPUT_ENERGY    = 100                                                     \n');

fprintf(fid, '# Compute the field int_0^t v^2 dt for a set of GLL points and write it to file. Use      \n');
fprintf(fid, '# the script utils/visualisation/plotIntegratedEnergyFile.py to watch. It is refreshed at the same time than the seismograms      \n');
fprintf(fid, 'COMPUTE_INTEGRATED_ENERGY_FIELD = .false.                                                 \n');

fprintf(fid, '#-----------------------------------------------------------------------------\n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '# movies/images/snaphots                                                      \n');
fprintf(fid, '#                                                                             \n');
fprintf(fid, '#-----------------------------------------------------------------------------\n');


fprintf(fid, '# every how many time steps we draw JPEG or PostScript pictures of the simulation                                \n');
fprintf(fid, '# and/or we dump results of the simulation as ASCII or binary files (costly, do not use a very small value)      \n');
fprintf(fid, 'NSTEP_BETWEEN_OUTPUT_IMAGES     = 200                                                                            \n');

fprintf(fid, '# minimum amplitude kept in for the JPEG and PostScript snapshots; amplitudes below that are muted               \n');
fprintf(fid, 'cutsnaps                        = 1.                   \n');

fprintf(fid, '#### for JPEG color images ####                        \n');
fprintf(fid, 'output_color_image              = .true.               \n'); 
fprintf(fid, 'imagetype_JPEG                  = 2                    \n'); 
fprintf(fid, 'factor_subsample_image          = 1.0d0                \n'); 
fprintf(fid, 'USE_CONSTANT_MAX_AMPLITUDE      = .false.              \n'); 
fprintf(fid, 'CONSTANT_MAX_AMPLITUDE_TO_USE   = 1.17d4               \n');
fprintf(fid, 'POWER_DISPLAY_COLOR             = 0.30d0               \n'); 
fprintf(fid, 'DRAW_SOURCES_AND_RECEIVERS      = .true.               \n'); 
fprintf(fid, 'DRAW_WATER_IN_BLUE              = .true.               \n'); 
fprintf(fid, 'USE_SNAPSHOT_NUMBER_IN_FILENAME = .false.              \n'); 

fprintf(fid, '#### for PostScript snapshots ####                     \n');
fprintf(fid, 'output_postscript_snapshot      = .false.              \n'); 
fprintf(fid, 'imagetype_postscript            = 1                    \n');
fprintf(fid, 'meshvect                        = .true.               \n');
fprintf(fid, 'modelvect                       = .false.              \n'); 
fprintf(fid, 'boundvect                       = .true.               \n'); 
fprintf(fid, 'interpol                        = .true.               \n'); 
fprintf(fid, 'pointsdisp                      = 6                    \n');
fprintf(fid, 'subsamp_postscript              = 1                    \n');
fprintf(fid, 'sizemax_arrows                  = 1.d0                 \n');
fprintf(fid, 'US_LETTER                       = .false.              \n'); 

fprintf(fid, '#### for wavefield dumps ####                          \n'); 
fprintf(fid, 'output_wavefield_dumps          = .false.         # output wave field to a text file (creates very big files)                 \n');
fprintf(fid, 'imagetype_wavefield_dumps       = 1              # display 1=displ vector 2=veloc vector 3=accel vector 4=pressure           \n');
fprintf(fid, 'use_binary_for_wavefield_dumps  = .false.        # use ASCII or single-precision binary format for the wave field dumps      \n');

fprintf(fid, '#-----------------------------------------------------------      \n');
fprintf(fid, 'NUMBER_OF_SIMULTANEOUS_RUNS     = 1                               \n');
fprintf(fid, 'BROADCAST_SAME_MESH_AND_MODEL   = .true.                          \n');
fclose(fid);




end 
