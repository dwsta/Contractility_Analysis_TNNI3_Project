[xvec,yvec,tvec,u,v] = readPIV_bin('G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_TopLeft\contractility_run_20221105_190855\output\Well__E_002_Tritc\piv\smooth_deformations_pass2.bin');
[X,Y,T] = meshgrid(xvec,yvec,tvec);
div = divergence_rik(X,Y,u,v);
div_trace = squeeze(sqrt(mean(div.^2,[1 2],"omitnan")));

fps = 25;
raw_time = tvec/fps;
min_sampling_time = min(diff(raw_time));
time = min(raw_time):min_sampling_time:max(raw_time);
signal = interp1(raw_time,div_trace,time,'spline');
%%
myVideo = VideoWriter('animated_trace','MPEG-4'); %open video file
myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
open(myVideo)
fig = figure;
fig.Color = [1 1 1];
ax = axes();

for itime = 1:length(time)
    plot(time(1:itime), signal(1:itime), 'LineWidth', 2)
    box off
    pbaspect([5 1 1])
    ax.XLim = [0, max(time)];
    ax.YLim = [0, 0.03];
    ax.YTick = [0.00, 0.01, 0.02, 0.03];
    ax.LineWidth = 1;
    ax.FontSize = 12;
    drawnow
%     pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)
close(fig)
