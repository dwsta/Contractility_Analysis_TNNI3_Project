%% Load data and open output file connection
pivdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\output\Well__E_007_Tritc\piv';
filpiv = 'smooth_deformations_pass4.bin';
outcsv = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\results.csv';
ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
[xvec,yvec,tvec,U,V] = readPIV_bin(fullfile(pivdir,filpiv));
% Remove the measurements at reference frames, as it can introduce errors
% later in the analysis
vecind = load(fullfile(pivdir,'vecind.txt'));
[tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
tvec(idx)=[];
U(:,:,idx,:)=[];
V(:,:,idx,:)=[];
output_file_connection = fopen(outcsv,'w');
headers = {'Alias','T peak','T rise','T fall','TCT',...
    'D peak','D valley','D high','D low','AUC','Power', 'CR', 'RR',...
    'Trace','Time','Phase','Average Trace','dPA_dt','Time for Average'};
fprintf(output_file_connection,'%s\n',strjoin(headers,','));
%% Compute divergence
raw_time = tvec;
[X,Y] = meshgrid(xvec,yvec);
div = divergence_rik(X,Y,U,V);
raw_signal = squeeze(mean(sqrt(div.^2),[1 2],'omitnan'));
alias = 'An alias';
%%
computeParameters_phase_avg(pivdir,alias,raw_time,raw_signal,output_file_connection)
fclose(output_file_connection)