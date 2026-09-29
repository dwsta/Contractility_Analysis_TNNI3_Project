function register_peak_contraction_frames(rootdir,alias,time,signal,frames)
[peakIDs,time_interp,signal_interp] = split_peaks(time,signal);

analysisdir = fullfile(rootdir,'output',alias,'whole-ROI_analysis');
mkdir(analysisdir)
% Make a file with the frames at peak contracion (to be used in
% elastography)
fname_peaks = fullfile(analysisdir,'peak_frames.txt');
fid = fopen(fname_peaks,'w');
peak_frames = [];
nPeaks = max(peakIDs);
for iPeak = 1 : nPeaks
    inPeak = peakIDs == iPeak;
    [~,loc] = max(inPeak.*signal_interp);
    [~,ind] = min(abs(time-time_interp(loc)));
    peak_frames(iPeak) = frames(ind);
end
fprintf(fid,'%d ',peak_frames);
fclose(fid);
