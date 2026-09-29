function cell_by_cell_analysis(rootdir,alias,cfg_data,output_file_connection)

analysisdir = fullfile(rootdir,'output',alias,'cell-by-cell_analysis');
mkdir(analysisdir)

% Load cell segmentation
maskdir = fullfile(rootdir,'output',alias,'mask');

%DWS - in the case of just providing a mask, just get
%the file copied as "mask_*.tif"
%if cellpose, get the cellpose output of "*cp_masks.png"
if cfg_data.SingleCell.UseCellpose
    filmask= dir([maskdir,filesep,'*_cp_masks.png']);
else
    filmask= dir([maskdir,filesep,'mask_*.tif']);
end
mask = imread(fullfile(filmask.folder,filmask.name));
props = regionprops(mask,'Area','Centroid','Orientation','Circularity','MajoraxisLength','MinoraxisLength','Perimeter');
Tprops = struct2table(props);
[trace_matrix,frames] = compute_traces_cell_by_cell(cfg_data, rootdir, alias, mask);

if cfg_data.TFM.DoTFM
    fnameout = fullfile(analysisdir,[alias,'_msm_trace_raw.pdf']);
    fnameout2 = fullfile(analysisdir,[alias,'_msm_trace_peaks.pdf']);
else
    fnameout = fullfile(analysisdir,[alias,'_div_trace_raw.pdf']);
    fnameout2 = fullfile(analysisdir,[alias,'_div_trace_peaks.pdf']);
end

tvec= (frames-min(frames)) / cfg_data.OtherParameters.FrameRate;
save_traces(trace_matrix,[],tvec,fnameout,cfg_data)
Ntraces = size(trace_matrix,1);

for itrace = 1 : Ntraces
[peakIDs,time,signal] = split_peaks(tvec,trace_matrix(itrace,:));
trace_matrix_interp(itrace,:) = signal;
peakIDs_matrix_interp(itrace,:) = peakIDs;
metrics = measure_peaks(peakIDs,time,signal);
t = struct2table(metrics);

means = varfun(@nanmean, t(:,2:end), 'InputVariables', @isnumeric);
sds = varfun(@nanstd, t(:,2:end), 'InputVariables', @isnumeric);

t.alias = repmat(alias,[height(t),1]);
t.cellID = repmat(itrace,[height(t),1]);
t = [t(:,end-1:end) t(:,1:end-2)];
if itrace == 1 
    bigT = t;
else
    bigT = [bigT;t];
end
signal_str = ['{', sprintf('%f;',signal) ,'}'];
time_str   = ['{', sprintf('%f;',time)   ,'}'];
peakID_str = ['{', sprintf('%d;',peakIDs)   ,'}'];
props_str = arrayfun(@num2str,Tprops{itrace,:},'UniformOutput',false);
means_str = arrayfun(@num2str,means{1,:},'UniformOutput',false);
sds_str  = arrayfun(@num2str,sds{1,:},'UniformOutput',false);
stringout = strjoin([alias, num2str(itrace) , props_str ,signal_str , time_str , peakID_str, means_str,sds_str],',');
fprintf(output_file_connection,'\n%s',stringout);

end
save_traces(trace_matrix_interp,peakIDs_matrix_interp,time,fnameout2,cfg_data)

peaksfilename = fullfile(analysisdir,[alias,'_peak_metrics.csv']);
writetable(bigT,peaksfilename)
