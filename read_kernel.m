M = 160;
N = 80;
fid = fopen('proc000000_alpha_kernel.bin', 'rb');
% [data, count] = fread(fid, [M,N],'float');
[data, count] = fread(fid, inf,'int');

fclose(fid);
imagesc(data);