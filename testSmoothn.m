pivfile = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\output\Well__E_007_Tritc\piv\deformations_pass1.bin';
[x,y,t,u,v] = readPIV_bin(pivfile);
uc = u;
vc = v;
ug = gpuArray(uc);
vg = gpuArray(vc);
% GPU is about 4x faster with individual vector components but the function breaks if input vectors {u,v}
tic
zc = smoothn(uc,'robust');
toc
tic
zg = smoothn(ug,'robust');
toc
tic
zg = smoothn(vg,'robust');
toc

% Does it make a big difference between running u,v separate vs together as
% vector?
% tic
% usc = smoothn(uc,'robust');
% toc
% tic
% vsc = smoothn(vc,'robust');
% toc
% tic
% zsc = smoothn({uc,vc},'robust');
% toc

scl = 0.2;
figure(1); clf
t = tiledlayout(1,2);
ax1 = nexttile; 
quiver(uc(:,:,25)*scl,vc(:,:,25)*scl,0,'r'); 
axis image
ax2 = nexttile; 
quiver(zsc{1}(:,:,25)*scl,zsc{2}(:,:,25)*scl,0,'b');
axis image
linkaxes([ax1,ax2])




