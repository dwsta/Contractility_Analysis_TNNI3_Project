rootdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\contractility_run_20221122_211235';
alias = 'Fitc_2';
cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));

Ngrids = length(cfg_data.Deformation);
Ngrids = 1;
pivdir = fullfile(rootdir,'output',alias,'piv');
filpiv = ['deformations_pass',num2str(Ngrids),'.bin'];
ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
[xvec,yvec,tvec,U,V] = readPIV_bin(fullfile(pivdir,filpiv));
U = U(:,:,1:50);
[U1] = smoothn_rik(U,0.001,'robust');
U3 = smoothn_rik(U,0.0691);
