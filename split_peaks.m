% Copyright (c) 2022, Ricardo Serrano
% All rights reserved.
function [peakID,time,signal] = split_peaks(raw_time,raw_signal,varargin)
% Aim: Find start of beats as "kinks" in the signal. We will do that by
% finding big differences in slope of locally fitted lines (with least-squares) 
% using a stencil that takes the point and points forward (controlled by stencil).
%
% Note: At the beginning of the signal, there may be a portion of peak, but if no 
% kink is observed, these piece of peak will receive the value (peakID==0) and it
% won't be analyzed, as it is incomplete.
stencil = 2;

% If peaks are negative, flip the sign
if nanmean(raw_signal(:))<0 
    raw_signal = -raw_signal;
    sign_flag = -1;
else
    sign_flag = 1;
end

% The contractility algorithm may have skipped frames of low contractility,
% generating a non-uniform time signal. Interpolate the signal to match
% highest frequency sampling in the signal (i.e. minimum sampling period)
dx = min(diff(raw_time));
time = min(raw_time):dx:max(raw_time);

nanpt = isnan(raw_signal);
if all(nanpt) 
    peakID = zeros(size(time));
    signal = nan(size(time));
    return
else
    raw_signal(nanpt) = [];
    raw_time(nanpt)=[];
end

signal = interp1(raw_time,raw_signal,time,'spline');

x = time;
y = signal;

%DWS edit - make this callable

old_method = 0;

disp(nargin);

if nargin>2
    if varargin{1}=="findkink"
        old_method = 1;
    end
end

if old_method == 1
    xbar = stencil*dx/2; % For uniform time resolution and stencil size, this is constant.The interpolation above ensures this uniformity.
    ybar = movmean(signal,[0 stencil]);
    ybar = ybar(stencil:end-stencil-1); % discard starting and ending points, where either trailing or forward fits won't have enough data
    Npts = length(ybar);
    x_minus_xbarlead = (dx*[0:stencil])-xbar; % True for  equally spaced data
    for irow = 1:Npts
        num(irow) = sum(x_minus_xbarlead.*(y(irow+stencil:irow+2*stencil)-ybar(irow)),2);
    end
    den = sum(x_minus_xbarlead.^2,2);
    slope = num./den';
    ordinate = ybar - slope.*xbar;
    
    [kink_peakval,kink_locs] = findpeaks(slope,'MinPeakHeight',0.25*max(slope));
    peakID = cumsum(ismember(1:length(y),kink_locs-stencil));
else
    disp("DWS_method");
    %this is based on the following assumptions:
    % 1) peaks are upright
    % 2) all peaks start at a local minima and end at a local minima
    % 3) peaks cannot start at index 1 or max index
    % 4) peaks have a prominence of at least 50% of the max-min of the
    % signal (this will break with photobleaching)
    min_peak_time = 0.3; %200 bpm = 60/200 = 0.3 - maximum possible rate allowed

    min_peak_frames = ceil(min_peak_time/dx); %minimum distance between peaks


    % alternate calculation (experimental)
    % perform a cross correlation
    % find the first peak (shortest possible repetitive signal)
    % choose half of that or the maximum chosen above
    cross_corr = xcorr(signal-mean(signal),'biased');
    cross_corr = cross_corr(length(signal):length(cross_corr)); %only choose second half
    %find lag locations
    [dummy,xcor_loc] = findpeaks(cross_corr);
    % first peak will be the minimum lag for a repetitive signal
    % choose half of that or the maximum above (200bpm default)
    % whichever is larger
    min_peak_frames = max(min_peak_frames,ceil(xcor_loc(1)/2));
    
    max_signal = max(signal);
    min_signal = min(signal);
    max_amplitude = max_signal-min_signal;
    
    prominence_threshold = max_amplitude * 0.5;
    
    %find peaks that are at least min_peak_frames apart and
    %prominence_threshold prominent
    [pks,maxlocs,w,p] = findpeaks(signal,"MinPeakDistance",min_peak_frames);
   
    %filter out low peaks
    %minimum height will be 1/5 of the max prominence (arbitrary cutoff)
    %originally thought of median but ran into situation
    %where a lot of small peaks were counted
    max_prominence = max(p);
    maxlocs = maxlocs(find(p>max_prominence/5));
    % 
    % figure;
    % findpeaks(signal,"MinPeakDistance",min_peak_frames,'Annotate','extents');
    % p
    % median(p)

    %create vector of minima
    minlocs = NaN([length(maxlocs)-1 1]);
    
    %each maximum should be surrounded by two minima which delimit the peak
    for i = 1:length(maxlocs)-1
        [m,mloc] = min(signal(maxlocs(i):maxlocs(i+1)));
        minlocs(i) = maxlocs(i)-1+mloc; %-1 because index starts at 1, so if the minimum were at the left end, want minlocs(i) = maxlocs(i)
    end
    
    %left side of plot
    [m,leftmin] = min(signal(1:maxlocs(1)));
    
    if leftmin == 1
        maxlocs = maxlocs(2:length(maxlocs)); % get rid of incomplete peak on left
    else
        minlocs = [leftmin;minlocs]; % add minimum on left
    end
    
    %right side of plot
    [m,rightmin] = min(signal(maxlocs(length(maxlocs)):length(signal)));
    rightmin = rightmin + maxlocs(length(maxlocs))-1;
    
    if rightmin == length(signal)
        maxlocs = maxlocs(1:length(maxlocs)-1); %get rid of incomplete peak on right
    else
        minlocs = [minlocs;rightmin]; %add minimum on right
    end
    
    %peak_to_peak = [diff(maxlocs);NaN]; %last one will always have no length
    
    peakID = cumsum(ismember(1:length(signal),minlocs)); %note first frames will be 0 and not a peak
    peakID(peakID==max(peakID)) = 0; %mark end as not a peak
end
signal = sign_flag*signal;
% plot(x,y,'r'); hold on
% plot(x(kink_locs),y(kink_locs),'ko');
% plot_multicolor_rik(x,y,peakID)
