function [TFM_results] = TFM_step(rootdir,alias,cfg_data)

% Load PIV
Ngrids = length(cfg_data.Deformation);
pivdir = fullfile(rootdir,'output',alias,'piv');
PIV_file = fullfile(pivdir,['smooth_deformations_pass',num2str(Ngrids),'.bin']);
E = cfg_data.TFM.GelStiffness*1e3;
h = cfg_data.TFM.GelThickness;
fcalx = cfg_data.TFM.CalibrationFactor;
[xvec,yvec,tvec,U,V,Xdrift,Ydrift] = readPIV_bin(PIV_file);
Um = nanmean(U,4);
Vm = nanmean(V,4);
Ntimepoints = length(tvec);
nx = length(xvec);
ny = length(yvec);
nx0 = nx + mod(nx,2);
ny0 = ny + mod(ny,2);

[Xgrid,Ygrid] = meshgrid(xvec,yvec);
Nyq = 2;
TX = nan(ny0,nx0,Ntimepoints);
TY = TX;
parfor itimepoint = 1 : Ntimepoints
    u = Um(:,:,itimepoint);
    v = Vm(:,:,itimepoint);
    [tx,ty,nxx,nxy,nyy,ut0,vt0,X0,Y0,Lx,Ly]...
        = mytrac_msm_2D_rik_from_dp(E,h,h,nx0,ny0,Xgrid,Ygrid,u,v,fcalx,Nyq);
    TX(:,:,itimepoint) = tx;
    TY(:,:,itimepoint) = ty;

end

TFM_xvec = linspace(xvec(1),xvec(end),nx0);
TFM_yvec = linspace(yvec(1),yvec(end),ny0);
TFM_tvec = tvec;
TFM_repvec = 1;
[NY,NX,NF,NR] = size(TX);

% Output
TFM_dir = fullfile(rootdir,'output',alias,'tfm');
if ~exist(TFM_dir,'dir'); mkdir(TFM_dir); end

fid2 = fopen(fullfile(TFM_dir,'traction_stresses.bin'),'w');
fwrite(fid2,'v3','uchar');
fwrite(fid2,single([NX,NY,NF,NR]),'single');
fwrite(fid2,single(TFM_xvec),'single');
fwrite(fid2,single(TFM_yvec),'single');
fwrite(fid2,single(TFM_tvec),'single');
fwrite(fid2,single(TFM_repvec),'single');
fwrite(fid2,single(Xdrift),'single');
fwrite(fid2,single(Ydrift),'single');
fwrite(fid2,single(TX),'single');
fwrite(fid2,single(TY),'single');
fclose(fid2);
