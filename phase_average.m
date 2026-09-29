% Copyright (c) 2022, Ricardo Serranoddiv
% All rights reserved.param_out

function phase_average(cfg_data,pivdir,alias,frames,raw_signal,output_file_connection)
try
% headers = {'Alias','T peak','T rise','T fall','TCT',...
%     'D peak','D valley','D high','D low','AUC','Power', 'CR', 'RR',...
%     'Trace','Time','Frame,'Phase','Average Trace','dPA_dt','Time for Average'};
raw_time = frames / cfg_data.OtherParameters.FrameRate;

% The contractility algorithm may have skipped frames of low contractility,
% generating a non-uniform time signal. Interpolate the signal to match
% highest frequency sampling in the signal (i.e. minimum sampling period)
min_sampling_time = min(diff(raw_time));
time = min(raw_time):min_sampling_time:max(raw_time);
signal = interp1(raw_time,raw_signal,time,'spline');

smooth_signal = smooth(signal,0.005,'rloess');

discard = round(4*0.005*length(smooth_signal));    % To remove end effects
smooth_signal = smooth_signal(discard:end-discard);
[acor,lag] = xcorr(smooth_signal-mean(smooth_signal),'biased');

[acor_pks,lag_locs] = findpeaks(acor,'minpeakheight',max(acor)/5);

% lag_locs = lag_locs+discard;                % Shift the location of lag to match the
acor_pks = acor_pks(lag(lag_locs)>=0);
lag_locs = lag_locs(lag(lag_locs)>=0);

% [acor_pks,I] = sort(acor_pks,'descend');
% lag_locs = lag_locs(I);
if length(lag_locs) > 1
    [t,div,phase_time,s_avg,t_peak] = phase_average(signal,time,raw_signal,raw_time,lag,lag_locs,smooth_signal,pivdir);
else
    div = signal;
    t = time;
    phase_time = raw_time;
    s_avg = nan;
    t_peak = nan;
end


%
%   Find crossing over d_low and determine if it is the beginning of a
%   contraction (end of relaxed state) or the end of a relaxation
%   (beginning of relaxed state)
%

a = div(1:end-1)-d_low;
b = div(2:end)-d_low;
cross_low_event_indx = find(a.*b<0);
cross_low_event_indx(cross_low_event_indx==1)=[];
cross_low_event_indx(cross_low_event_indx==length(div))=[];
contraction_start_time = [];
relaxation_end_time = [];
for ievent = 1:length(cross_low_event_indx)
    indx_cross = cross_low_event_indx(ievent);
    stencil = [indx_cross-3:indx_cross+3];
    stencil = stencil(stencil>1&stencil<length(div)); % Ensure the stencil does not go out of borders
    div_interpolation_segment = div(stencil);
    time_interpolation_segment = t(stencil);
    time_event = interp1(div_interpolation_segment,time_interpolation_segment,d_low);
    if mean(diff(div_interpolation_segment))>0 % contraction
        contraction_start_time = [contraction_start_time,time_event];
        
    else
        relaxation_end_time = [relaxation_end_time,time_event];
    end
    
end

%
%   Find crossing over d_high and determine if it is the end of a
%   contraction or the beginning of a relaxation
%
a = div(1:end-1)-d_high;
b = div(2:end)-d_high;

cross_high_event_indx = find(a.*b<0) ;


contraction_end_time = [];
relaxation_start_time = [];
for ievent = 1:length(cross_high_event_indx)
    indx_cross = cross_high_event_indx(ievent);
    stencil = [indx_cross-3:indx_cross+3];
    stencil = stencil(stencil>0&stencil<length(div)); % Ensure the stencil does not take points out of bounds
    div_interpolation_segment = div(stencil);
    time_interpolation_segment = t(stencil);
    
    time_event = interp1(div_interpolation_segment,time_interpolation_segment,d_high);
    if mean(diff(div_interpolation_segment))>0 % contraction
        contraction_end_time = [contraction_end_time,time_event];
    else
        relaxation_start_time = [relaxation_start_time,time_event];
    end
    
end

[~,indx_max_div]=max(div);
time_max_div = t(indx_max_div);

% contraction happens before maximum
contraction_start_time(contraction_start_time>time_max_div)=[];
contraction_end_time(contraction_end_time>time_max_div)=[];

% relaxation after maximum
relaxation_start_time(relaxation_start_time<time_max_div)=[];
relaxation_end_time(relaxation_end_time<time_max_div)=[];

if range_rik(contraction_start_time)<0.2; contraction_start_time= mean(contraction_start_time); end
if range_rik(contraction_end_time)<0.2;   contraction_end_time  = mean(contraction_end_time); end
if range_rik(relaxation_start_time)<0.2; relaxation_start_time= mean(relaxation_start_time); end
if range_rik(relaxation_end_time)<0.2; relaxation_end_time= mean(relaxation_end_time); end



t_rise = contraction_end_time-contraction_start_time;
t_fall = relaxation_end_time-relaxation_start_time;


if length(contraction_start_time) == length(relaxation_end_time)
    total_contraction_time = relaxation_end_time-contraction_start_time;
else
    total_contraction_time = nan;
end

if length(t_rise)~= 1
    t_rise = nan;
end
if length(t_fall)~= 1
    t_fall = nan;
end

param_out{3} = num2str(t_rise);
param_out{4} = num2str(t_fall);
param_out{5} = num2str(total_contraction_time);

% Area Under Curve and Power (3/14/2017)
auc = trapz(t,div)-d_valley*(t(end)-t(1));
pow = auc/t_peak;

param_out{10} = num2str(auc);
param_out{11} = num2str(pow);


% figure(2);clf
set(0,'CurrentFigure',figHandles(2));clf

plot(t,s_avg','.');
hold on
plot(t,div,'LineWidth',3)

line([xlim],[d_peak, d_peak],'Color','k','LineStyle','-.')
line([xlim],[d_low, d_low],'Color','k','LineStyle',':')
line([xlim],[d_high, d_high],'Color','k','LineStyle',':')
line([xlim],[d_valley, d_valley],'Color','k','LineStyle','-.')

axis([t(1) t(end) , ylim*1.1])

yr = range_rik(ylim);

text(t(end)*0.03,d_peak+0.03*yr,sprintf('D peak = %1.4e',d_peak));
text(t(end)*0.03,d_high+0.03*yr,sprintf('D high = %1.4e',d_high));
text(t(end)*0.03,d_low+0.03*yr,sprintf('D low = %1.4e',d_low));
text(t(end)*0.03,d_valley+0.03*yr,sprintf('D valley = %1.4e',d_valley));
text(t(end)*0.03,0.7*yr,sprintf('T peak = %4.4f',t_peak));
text(t(end)*0.03,0.6*yr,sprintf('T rise = %4.4f',t_rise));
text(t(end)*0.03,0.5*yr,sprintf('T fall = %4.4f',t_fall));
text(t(end)*0.03,0.4*yr,sprintf('TCT = %4.4f',total_contraction_time));
text(t(end)*0.03,0.3*yr,sprintf('AUC = %4.4f',auc));
text(t(end)*0.03,0.2*yr,sprintf('Power = %4.4f',pow));


set(figHandles(2), 'PaperPositionMode', 'manual');
set(figHandles(2), 'PaperUnits', 'points');
set(figHandles(2), 'PaperPosition', [1 1 512 512]);
print(figHandles(2),'-dpng','-r300',fullfile(pivdir,'div_param_PA.png'))

ddivdt = diff(div)./diff(t);
ddivdt(end+1) = nan;
intertime = t;
param_out{19} = arr2str_rik(ddivdt);
[contraction_rate,relaxation_rate] = rate_computation(ddivdt,intertime,pivdir);
param_out{12} = num2str(contraction_rate);
param_out{13} = num2str(relaxation_rate);

fid = fopen(fullfile(pivdir,'contractility_params_PA.txt'),'w');
fprintf(fid,'T_peak: %7.9f \n',t_peak);
fprintf(fid,'T_rise: %7.9f \n',t_rise);
fprintf(fid,'T_fall: %7.9f \n',t_fall);
fprintf(fid,'TCT: %7.9f \n',total_contraction_time);
fprintf(fid,'D_peak: %7.9f \n',d_peak);
fprintf(fid,'D_valley: %7.9f \n',d_valley);
fprintf(fid,'D_high: %7.9f \n',d_high);
fprintf(fid,'D_low: %7.9f \n',d_low);
fprintf(fid,'AUC: %7.9f \n',auc);
fprintf(fid,'Power: %7.9f \n',pow);
fprintf(fid,'CR: %7.9f \n',contraction_rate);
fprintf(fid,'RR: %7.9f \n',relaxation_rate);
fclose(fid);

fid = fopen(fullfile(pivdir,'contraction_parameters_bin.bin'),'w');
fwrite(fid,[t_peak,t_rise,t_fall,total_contraction_time,d_peak,d_valley,d_high,d_low,...
    auc,pow,contraction_rate,relaxation_rate],'double');
fclose(fid);

outstr = [strjoin(param_out,','),'\n'];
fprintf(output_file_connection,outstr);
close(figHandles);
catch ME
    outstr = [strjoin(param_out,','),'\n'];
    fprintf(output_file_connection,outstr);
    close(figHandles);
    rethrow(ME)
end
end
function rg = range_rik(arr)
    rg = max(arr(:))-min(arr(:));
end
function [ddivdt_positive_peak,ddivdt_negative_peak] = rate_computation(ddivdt,time_interval,filpiv)
figHandles = createFigHandles(1);

% ddivdt_raw = ddivdt;
% filterSize = 4;
% b = (1/filterSize)*ones(1,filterSize);
% a = 1;
% ddivdt2 = filtfilt(b,a,ddivdt);
% ddivdt = ddivdt(ceil(filterSize/2+1):end-ceil(filterSize/2+1));

time_interval_positive = zeros(size(ddivdt));
time_interval_negative = zeros(size(ddivdt));

ddivdt_positive = zeros(size(ddivdt));
ddivdt_negative = zeros(size(ddivdt));

% Plot red-blue line if its positive-negative
time_interval_positive = time_interval;
time_interval_positive(ddivdt<0) = NaN;
ddivdt_positive = ddivdt;
ddivdt_positive(ddivdt<0) = NaN;

% In reality we plot all the whole line in blue (as if it were negative)
% but the red will be overimposed
time_interval_negative = time_interval;
ddivdt_negative = ddivdt;

set(0,'CurrentFigure',figHandles(1));clf
plot(time_interval,ddivdt,'Color',[.75 .75 .75],'LineWidth',0.5)
hold on
plot(time_interval_negative,ddivdt_negative,'b');
plot(time_interval_positive,ddivdt_positive,'r');

axis([time_interval(1) time_interval(end) , ylim*1.1]);
set(figHandles(1), 'PaperPositionMode', 'manual');
set(figHandles(1), 'PaperUnits', 'points');
set(figHandles(1), 'PaperPosition', [1 1 1024 512]);
print(figHandles(1),'-dpng','-r300',fullfile(filpiv,'ddivdt_PA.png'))
set(0,'CurrentFigure',figHandles(1));

[ddivdt_positive_peak,ddivdt_positive_loc_indx] = max(ddivdt);            %  contraction rate
line([xlim],[ddivdt_positive_peak, ddivdt_positive_peak],'Color','r','LineStyle','-.')


[ddivdt_negative_peak,ddivdt_negative_loc_indx] = min(ddivdt);
line([xlim],[ddivdt_negative_peak, ddivdt_negative_peak],'Color','b','LineStyle','-.')

line([xlim],[0, 0],'Color','k','LineStyle','-','LineWidth',1)

yr = max(abs([min(ddivdt_negative),max(ddivdt_positive)]));
text(time_interval(end)*0.03,ddivdt_positive_peak+0.05*yr,sprintf('CR = %1.4e',ddivdt_positive_peak));
text(time_interval(end)*0.03,ddivdt_negative_peak-0.1*yr,sprintf('RR = %1.4e',ddivdt_negative_peak));
xr = xlim;
axis([xr(1),xr(2),-1.15*yr,1.15*yr])
% line([time_interval(ddivdt_positive_loc_indx);time_interval(ddivdt_positive_loc_indx)],[ylim],'Color','r','Linestyle','--')
% line([time_interval(ddivdt_negative_loc_indx);time_interval(ddivdt_negative_loc_indx)],[ylim],'Color','b','Linestyle','--')

set(figHandles(1), 'PaperPositionMode', 'manual');
set(figHandles(1), 'PaperUnits', 'points');
set(figHandles(1), 'PaperPosition', [1 1 512 512]);
print(figHandles(1),'-dpng','-r300',fullfile(filpiv,'ddivdt_param_PA.png'))

close(figHandles)
end
function [t,div,phase_time,s_avg,t_peak] = phase_average(signal,time,raw_signal,raw_time,lag,lag_locs,smooth_signal,filpiv)
figHandles = createFigHandles(1);
mpd = lag(lag_locs(2)); % Min peak distance (in interpolated timepoints)

t_peak = time(mpd);      % Period of the signal

% [pks,locs] = findpeaks(smooth_signal,'minpeakheight',mean(signal),'minpeakdistance',floor(0.5*mpd));
[pks,locs] = findpeaks(smooth_signal,'minpeakheight',1.25*mean(signal),'minpeakdistance',floor(0.5*mpd));

pks(locs<10) = [];
locs(locs<10) = [];
pks(locs>(length(signal)-10)) = [];
locs(locs>(length(signal)-10)) = [];

Npks = length(pks);
if Npks > 1
    windowSize = mpd+(4-mod(mpd,4)); % Take windowSize divisible by 4
    cycle_time = linspace(0,t_peak,windowSize);
    cycle_clock = [];
    current_clock = [];
    peakID = [];
    s_avg = nan(Npks,windowSize);
    for ipk = 1:Npks
        iloc = locs(ipk);
        iini = max(1,iloc-windowSize/4+1);
        iend = min(length(signal),iloc+3*windowSize/4);
        if (iloc-windowSize/4<=0)
            a = 1-(iloc-windowSize/4);
        else
            a = 1;
        end
        if (iloc+3*windowSize/4>=length(signal))
            b = windowSize-((iloc+3*windowSize/4)-length(signal));
        else
            b = windowSize;
        end
        piece = nan (1,windowSize);
        piece(a:b) = signal(iini:iend);
        
        if ipk == 1
            model = piece;
            model(isnan(model))=0;
        else
            [cor,phase] = xcorr(piece(~isnan(piece)),model(~isnan(piece)));
            %         correlations(ipk,:) = cor;
            %         phases(ipk,:) = phase;
            [~,I] = max(cor);
            iini = iini + phase(I);
            iend = iend + phase(I);
            if iini<1
                iini = 1;
                a = 1;
            end
            if iend > length(signal)
                iend = length(signal);
                b = 1+(iend-iini);
            end
        end
        piece(a:b) = signal(iini:iend);
        cycle_clock(ipk,a:b) = cycle_time(a:b)';
        current_clock(ipk,a:b) = time(iini:iend)';
        peakID(ipk,a:b) = ipk*ones(size(1,iend-iini));
        s_avg(ipk,a:b) = piece(a:b)';
    end
    % to output the phase time
    [c,ic] = unique(current_clock(:));
    dum = cycle_clock(:);
    d = dum(ic);
    phase_time = interp1(c,d,raw_time,'linear');
    
    set(0,'CurrentFigure',figHandles(1));clf
    plot(current_clock',s_avg');
    
    set(figHandles(1), 'PaperPositionMode', 'manual');
    set(figHandles(1), 'PaperUnits', 'points');
    set(figHandles(1), 'PaperPosition', [1 1 1024 512]);
    print(figHandles(1),'-dpng','-r300',fullfile(filpiv,'divPA_cycles.png'))
    
    idx = isnan(s_avg);
    buf = s_avg;
    buf(idx) = 0;
    div = sum(buf,1)./sum(~idx,1);
    t = time(1:windowSize);
    
    fid = fopen(fullfile(filpiv,'div_PA.txt'),'w');
    fprintf(fid,'%7.9f %7.9f\n',[t(:),div(:)]');
    fclose(fid);
    
    current_clock = reshape(current_clock',[],1);
    cycle_clock   = reshape(cycle_clock',[],1);
    cycle_number  = reshape(peakID',[],1);
    current_div   =	reshape(s_avg',[],1);
    
    TMP = [current_clock,cycle_clock,cycle_number,current_div];
    TMP(isnan(current_div),:) = [];
    
    
    fid = fopen(fullfile(filpiv,'synch_PA.txt'),'w');
    fprintf(fid,'%7.9f %7.9f %2.2d %7.9f\n',TMP');
    fclose(fid);
    
else
    div = signal;
    t = time;
    phase_time = raw_time;
    s_avg = nan;
    t_peak = nan;
end
close(figHandles)
end

function figHandles = createFigHandles(nFigs)
% close all
set(0, 'DefaultFigureVisible', 'off');
for i = 1:nFigs
    figure;clf
end

set(0, 'DefaultFigureVisible', 'on');
allFigHandles = findall(0, 'Type', 'figure');
figHandles = allFigHandles(1:nFigs);
set(figHandles(:), 'visible', 'on')
end

function [str]=arr2str_rik(arr)
allOneString = sprintf('%f;' , arr);
allOneString = allOneString(1:end-1);% strip final comma
str = ['{',allOneString,'}'];
end