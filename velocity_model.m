function [ rho, vp, vs ] = velocity_model( nx, nz)

%vel = load([model_dir,'/vel.mat']);

rho = ones(nz, nx) * 2000;
vp = ones(nz, nx) * 2500;
rho(30:40,70:90) = 2200;
vp(30:40,70:90) = 3000;

% vs = vp/sqrt(3);
vs = vp * 0;


end 