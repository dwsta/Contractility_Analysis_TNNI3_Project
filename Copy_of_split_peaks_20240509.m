% Copyright (c) 2022, Ricardo Serrano
% All rights reserved.
function [peakID,time,signal] = split_peaks(raw_time,raw_signal)
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
signal = sign_flag*signal;
% plot(x,y,'r'); hold on
% plot(x(kink_locs),y(kink_locs),'ko');
% plot_multicolor_rik(x,y,peakID)
