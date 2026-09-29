function [metrics] = measure_peaks(peakIDs,time_interp,signal_interp)

nPeaks = max(peakIDs);
metrics = struct;
if nPeaks == 0
    metrics.peakID = NaN;
    metrics.amplitude = NaN;
    metrics.duration = NaN;
    metrics.rise_time = NaN;
    metrics.fall_time = NaN;
    metrics.pw90 = NaN;
    metrics.pw50 = NaN;
    metrics.pw30 = NaN;
    metrics.peak_value = NaN;
    metrics.baseline = NaN;
    metrics.valley = NaN;
end

%DWS edit: save peak-to-peak
peak_x = NaN(nPeaks,1);

for iPeak = 1:nPeaks
    try
        y = signal_interp(peakIDs==iPeak);
        % If peaks are negative, flip the sign
        if nanmean(y(:))<0
            y = -y;
            sign_flag = -1;
        else
            sign_flag = 1;
        end

        x = time_interp(peakIDs==iPeak);
        %DWS_edit save start time before setting to 0
        startx = x(1);

        x = x-x(1);
        [value,loc] = max(y);

        % Fit a parabola given 3 points (maximum, and 2 neighbors)
        x1 = x(loc-1);
        x2 = x(loc);
        x3 = x(loc+1);
        y1 = y(loc-1);
        y2 = y(loc);
        y3 = y(loc+1);

        den = (x1-x2)*(x1-x3)*(x2-x3);
        a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / den;
        b = (x1^2*(y2-y3) + x3^2*(y1-y2) + x2^2*(y3-y1)) / den;
        c = (x2^2*(x3*y1 - x1*y3) + x2*(x1^2*y3 - x3^2*y1) + x1*x3*(x3-x1)*y2)/den;

        % The peak value and time at which peak occurs
        xpeak = -b/(2*a);
        ypeak = -b^2/(4*a) + c;
        valley = quantile(y(x<=xpeak),0.05); % A bit more robust than min(y)
        baseline = (ypeak-valley)*0.1+valley; % 10% amplitude of the peak
        amplitude = ypeak - baseline;

        %DWSedit save interpolated peak time
        peak_x(iPeak) = xpeak + startx;

        % Let's find the start and end of the peak by finding the closest points
        % where the signal crosses baseline value.
        % The baseline crossing before peak.time will be the start of the peak
        % the one after peak.time will be the end of the peak.

        upstroke_x = x(x<=xpeak);
        upstroke_y = y(x<=xpeak);
        upstroke_perc = (1-(upstroke_y-baseline)/amplitude)*100;
        downstroke_x = x(x>=xpeak);
        downstroke_y = y(x>=xpeak);
        downstroke_perc = (1-(downstroke_y-baseline)/amplitude)*100;

        up90_time = interp1(upstroke_perc, upstroke_x, 90);
        down90_time = interp1(downstroke_perc, downstroke_x, 90);
        rise_time = xpeak-up90_time;
        fall_time = down90_time-xpeak;
        pw90 = down90_time-up90_time;
        up50_time = interp1(upstroke_perc, upstroke_x, 50);
        down50_time = interp1(downstroke_perc, downstroke_x, 50);
        pw50 = down50_time - up50_time;
        up30_time = interp1(upstroke_perc, upstroke_x, 30);
        down30_time = interp1(downstroke_perc, downstroke_x, 30);
        pw30 = down30_time - up30_time;

        % Output
        metrics(iPeak).peakID = iPeak;
        metrics(iPeak).amplitude = sign_flag*amplitude;
        %metrics(iPeak).duration = max(x)-min(x);
        %DWSEdit
        metrics(iPeak).duration = NaN; %set to NaN for now
        metrics(iPeak).rise_time = rise_time;
        metrics(iPeak).fall_time = fall_time;
        metrics(iPeak).pw90 = pw90;
        metrics(iPeak).pw50 = pw50;
        metrics(iPeak).pw30 = pw30;
        metrics(iPeak).peak_value = ypeak;
        metrics(iPeak).baseline = baseline;
        metrics(iPeak).valley = valley;
    catch ME
        metrics(iPeak).peakID = iPeak;
        metrics(iPeak).amplitude = NaN;
        metrics(iPeak).duration = NaN;
        metrics(iPeak).rise_time = NaN;
        metrics(iPeak).fall_time = NaN;
        metrics(iPeak).pw90 = NaN;
        metrics(iPeak).pw50 = NaN;
        metrics(iPeak).pw30 = NaN;
        metrics(iPeak).peak_value = NaN;
        metrics(iPeak).baseline = NaN;
        metrics(iPeak).valley = NaN;
    end

end

%DWSEdit
%calculate peak to peak interval
try
    peak_to_peak = diff(peak_x);

    for iPeak = 1:(nPeaks-1)
        metrics(iPeak).duration = peak_to_peak(iPeak);
    end

    metrics(nPeaks).duration = NaN;
catch ME
    disp("problem with calculating duration");
end

end