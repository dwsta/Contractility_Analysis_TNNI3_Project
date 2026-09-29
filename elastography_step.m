function elastography_step(rootdir,alias,cfg_data)
% Flags for elastography
runNH = 1;
runLinE = 1;

% Load PIV
Ngrids = length(cfg_data.Deformation);
pivdir = fullfile(rootdir,'output',alias,'piv');
filpiv = ['smooth_deformations_pass',num2str(Ngrids),'.bin'];
[xvec,yvec,tvec,U,V,Xdrift,Ydrift] = readPIV_bin(fullfile(pivdir,filpiv));
U = U*cfg_data.TFM.CalibrationFactor;
V = V*cfg_data.TFM.CalibrationFactor;

[X,Y] = meshgrid(xvec,yvec);

% Load peak frames (used to compute elastography
analysisdir = fullfile(rootdir,'output',alias,'whole-ROI_analysis');
frames = load(fullfile(analysisdir,'peak_frames.txt'));
inds = find(ismember(tvec',frames));
% Run elastography
outdir = fullfile(rootdir,'output',alias,'tfm');
% [E,nu,G,K,LinE,NH] = func_runElastographyVideos(X,Y,U,V,inds,runLinE,runNH);
[E,nu,G,K,LinE,NH] = func_Pades_runElastographyVideos_rik(X,Y,U,V,inds,runLinE,runNH);

t = table;
t.alias=repmat(alias,[numel(inds),1]);
t.peakID = [1:height(t)]';
t.E = E';
t.nu = nu';
t.G = G';
t.K = K';
elastofilename = fullfile(outdir,[alias,'_elastography_metrics.csv']);
writetable(t,elastofilename)

elastoSolvers = struct;
elastoSolvers.LinE = LinE;
elastoSolvers.NH   = NH;

str = jsonencode(elastoSolvers);
elastoSolversfname = fullfile(outdir,[alias,'_elastography_convergence.json']);
fid = fopen(elastoSolversfname,'w');
fprintf(fid,'%s',str);
fclose(fid);

mean_E = nanmean(E);
mean_nu = nanmean(nu);

% Write out monolayer stresses for all timepoints

for i = 1 : size(U,3)
    deltaS = (xvec(2) - xvec(1))*cfg_data.TFM.CalibrationFactor;

    % Derivtive of Strain matrix
%     [dUdX1,dUdX2] = Func_1stDer_PadesScheme(squeeze(U(:,:,i)),deltaS,deltaS);
%     [dVdX1,dVdX2] = Func_1stDer_PadesScheme(squeeze(V(:,:,i)),deltaS,deltaS);
    dUdX1 = dfdx_pade_2D_rik(squeeze(U(:,:,i)),deltaS);
    dUdX2 = dfdy_pade_2D_rik(squeeze(U(:,:,i)),deltaS);
    dVdX1 = dfdx_pade_2D_rik(squeeze(V(:,:,i)),deltaS);
    dVdX2 = dfdy_pade_2D_rik(squeeze(V(:,:,i)),deltaS);
    

    Exx = dUdX1; Eyy = dVdX2; Exy = (dUdX2+dVdX1)./2;
    

    SXX(:,:,i) = (mean_E./(1-mean_nu.^2)) .* (Exx + mean_nu .* Eyy);
    SYY(:,:,i) = (mean_E./(1-mean_nu.^2)) .* (Eyy + mean_nu .* Exx);
    SXY(:,:,i) = (mean_E./(1+mean_nu)) .* (Exy);
end


MSM_repvec = 1;
[NY,NX,NF,NR] = size(SXX);

fid2 = fopen(fullfile(outdir,'monolayer_stresses.bin'),'w');
fwrite(fid2,'v3','uchar');
fwrite(fid2,single([NX,NY,NF,NR]),'single');
fwrite(fid2,single(xvec),'single');
fwrite(fid2,single(yvec),'single');
fwrite(fid2,single(tvec),'single');
fwrite(fid2,single(MSM_repvec),'single');
fwrite(fid2,single(Xdrift),'single');
fwrite(fid2,single(Ydrift),'single');
fwrite(fid2,single(SXX),'single');
fwrite(fid2,single(SYY),'single');
fwrite(fid2,single(SXY),'single');
fclose(fid2);

