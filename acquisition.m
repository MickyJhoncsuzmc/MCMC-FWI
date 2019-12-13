function [nx, nz, dt, nt, f0, xmin, xmax, zmin, zmax, rec_x, rec_z, rec_n, src_x, src_z, src_n]=acquisition()

% 25 m per element
nx = 160;
nz = 80;
nt = 3000;
dt = 0.0005;
f0 = 10;  

xmin = 0; xmax = 4000;
zmin = 0; zmax = 2000;

rec_x = [500 : 50 : 3500];
rec_z = (zmax - 100) * ones(size(rec_x));
rec_n = size(rec_x,2);

src_x = [500 : 5 : 3500];
src_z = (zmax - 100) * ones(size(src_x)); 
src_n = size(src_x,2);

end 