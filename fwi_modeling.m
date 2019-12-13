clear all;close all;clc;

%% load the velocity model
% model_dir='marmousi2_151_401_10m';
model_dir='simple_51_101_10m';
load(['./model/',model_dir,'/vel.mat']);
x=0:dx:(nx-1)*dx;x=x/1000;
z=0:dx:(nz-1)*dx;z=z/1000;
figure;subplot(211);imagesc(x,z,vel);
title('true velocity');ylabel('Z (km)');

load(['./model/',model_dir,'/acquisition.mat']);

% parallel_init;
% seis=zeros(nt,ng,ns);

tic;
parfor is=1:ns
    display(['Modeling, is=',num2str(is),' ns=',num2str(ns)]);
    seis=a2d_mod_abc24(vel,nbc,dx,nt,dt,s,sx(is),sz(is),gx,gz,isFS);
    write_bin(['./data/',model_dir,'/csg_raw_',num2str(is),'.bin'],seis);
end
toc;

load(['./model/',model_dir,'/v0.mat']);
subplot(212);imagesc(x,z,vel0);
title('initial velocity');xlabel('X (km)');ylabel('Z (km)');

tic;
parfor is=1:ns
    display(['Modeling, is=',num2str(is),' ns=',num2str(ns)]);
    seis=a2d_mod_abc24(vel0,nbc,dx,nt,dt,s,sx(is),sz(is),gx,gz,isFS);
    write_bin(['./data/',model_dir,'/csg_dir_',num2str(is),'.bin'],seis);
end
toc;

tic;
parfor is=1:ns
    csg=read_bin(['./data/',model_dir,'/csg_raw_',num2str(is),'.bin'],nt,ng);
    dir=read_bin(['./data/',model_dir,'/csg_dir_',num2str(is),'.bin'],nt,ng);
    refl=csg-dir;
    write_bin(['./data/',model_dir,'/csg_refl_',num2str(is),'.bin'],refl);
end
toc;
% parallel_stop;
