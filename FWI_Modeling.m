%  Full Waveform Inversion acquisition
%  By Haipeng Li               
%  at USTC, 18 Nov. 2019        
%  haipengl@mail.ustc.edu.cn    

clear all;close all;clc;
addpath(genpath(pwd),'-begin');


simulation_dir = '/data1/FWI_Project/data/observe_data_simulation';

[nx, nz, dt, nt, f0, xmin, xmax, zmin, zmax, rec_x, rec_z, rec_n, src_x, src_z, src_n] = acquisition();
[ rho, vp, vs ] = velocity_model( nx, nz);

tic;
write_par(simulation_dir, nx, nz, nt, dt, xmin, xmax, zmin, zmax, rec_n, rec_x, rec_z, src_x, src_z, f0, rho, vp, vs);
cd data/observe_data_simulation/
run_simulation(src_n, src_x(1), src_x(end), src_z(1), src_z(end), f0, simulation_dir );

cd ../../../
toc




