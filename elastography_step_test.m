%% Run the analysis
rootdir = 'G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\contractility_run_20230318_174730';
alias = 'B003_r0004_c0005_beads_interleaved';
cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));
elastography_step(rootdir,alias,cfg_data);

outdir = fullfile(rootdir,'output',alias,'tfm');
fname = fullfile(outdir,'monolayer_stresses.bin');
[x,y,t,SXX,SYY,SXY] = readMSM_bin(fname);

%% Plot
pframes = load(fullfile(rootdir,'output',alias,'whole-ROI_analysis','peak_frames.txt'));
indx = pframes(1);
[xt,yt,tt,Tx,Ty,Xdrif,Ydrift] = readPIV_bin(fullfile(rootdir,'output',alias,'tfm','traction_stresses.bin'));
[X,Y] = meshgrid(xt,yt);

figure(1)
imagesc(SXX(:,:,indx)+SYY(:,:,indx),'XData',x,'YData',y); hold on
quiver(X,Y,Tx(:,:,indx),Ty(:,:,indx))

