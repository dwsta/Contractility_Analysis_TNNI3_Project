rootdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\contractility_run_20221122_211235';
alias = 'Fitc_2';
cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));
cfg_data.TFM.DoTFM=0;
% Load cell segmentation
maskdir = fullfile(rootdir,'output',alias,'mask');
filmask= dir([maskdir,filesep,'*_cp_masks.png']);
mask = imread(fullfile(filmask.folder,filmask.name));
% props = regionprops(mask,'Area','Centroid','Orientation','Circularity','MajoraxisLength','MinoraxisLength','Perimeter');
fnameout = 'test_cell_by_cell_TFM.pdf';
[trace_matrix,frames] = compute_traces_cell_by_cell(cfg_data, rootdir, alias,mask);
tvec = frames /cfg_data.OtherParameters.FrameRate;
save_traces(trace_matrix,[],tvec,fnameout,cfg_data)
