rootdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter';
jobname = 'jobfile.csv';

cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));
t = readJobFile2(jobname, rootdir);
iexp = 1;
alias = t(iexp,:).Alias{1};

% Taken from: elastography_step(rootdir,alias,cfg_data)
% Load PIV
Ngrids = length(cfg_data.Deformation);
pivdir = fullfile(rootdir,'output',alias,'piv');
filpiv = ['smooth_deformations_pass',num2str(Ngrids),'.bin'];
[xvec,yvec,tvec,U,V] = readPIV_bin(fullfile(pivdir,filpiv));
U = U*cfg_data.TFM.CalibrationFactor;
V = V*cfg_data.TFM.CalibrationFactor;

[X,Y] = meshgrid(xvec,yvec);


% Load peak frames (used to compute elastography)
analysisdir = fullfile(rootdir,'output',alias,'whole-ROI_analysis');
peakframesfile =fullfile(analysisdir,'peak_frames.txt');
if ~exist(peakframesfile,'file')
    output_filename= fullfile(rootdir,'foo.csv');
    output_file_connection = fopen(output_filename,'w');
    whole_ROI_analysis(rootdir,alias,cfg_data,output_file_connection)
end
frames = load(peakframesfile);

inds = find(ismember(tvec',frames));
% Run elastography
outdir = fullfile(rootdir,'output',alias,'tfm');
% [E,nu,G,K,LinE,NH] = func_runElastographyVideos(X,Y,U,V,inds,runLinE,runNH);
runLinE = 1;
runNH = 1;
[E,nu,G,K,LinE,NH] = func_Pades_runElastographyVideos_rik(X,Y,U,V,inds,runLinE,runNH);

