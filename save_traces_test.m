pivdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\output\Well__E_007_Tritc\piv';
filpiv = 'smooth_deformations_pass4.bin';
ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
[xvec,yvec,tvec,U,V] = readPIV_bin(fullfile(pivdir,filpiv));
% Remove the measurements at reference frames, as it can introduce errors
% later in the analysis
vecind = load(fullfile(pivdir,'vecind.txt'));
[tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V

tvec(idx)=[];
U(:,:,idx,:)=[];
V(:,:,idx,:)=[];

raw_time = tvec;
[X,Y] = meshgrid(xvec,yvec);
div = divergence_rik(X,Y,U,V);
raw_signal = squeeze(mean(sqrt(div.^2),[1 2],'omitnan'));


% [xvec,yvec,frames,U,V] = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\contractility_run_20221105_184302\output\Well__E_002_Tritc\piv'
cfg_data = loadJsonConfig('G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\contractility_run_20221105_184302\config.json');

fnameout = 'print_test.pdf';

[peakIDs,frame_interp,signal_interp] = split_peaks(raw_time,raw_signal);

% p = plot_multicolor_rik(tvec,trace_matrix,peakIDs)
save_traces(signal_interp,peakIDs,frame_interp,fnameout,cfg_data)

[metrics] = measure_peaks(peakIDs,frame_interp,signal_interp)

