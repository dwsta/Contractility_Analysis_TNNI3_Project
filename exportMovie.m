function exportMovie(quiverdata,cdata,cellboundaries,outputres,framerate,outputdir,outputname,xl,yl)
% quiverdata = structure (can be []) with the fields
% .xvec, .yvec, .tvec = vectors with x,y,t gridpoints
% .U, .V = [length(yvec) by length(xvec) by length(tvec)] matrices with the
%           U and V vector data
% .qscale = scale for quiver
% .qcolor = color of the arrows
%
%
% cdata  = structure (can be empty) with color data for imagesc
% .XData .YData: to be used as imagesc's 'XData' 'YData'
% .T = vector of timepoints
% .caxis = color axis limits for cdata
% .cmap  = colormap for cdata
%
% cellboundaries = structure (can be empty) with cell boundary
%
% outputres = resolution of the output video
% framerate = framerate of the output video
% outformat = format of the output video
% outputdir  = directory to save the output video
% outputname = name of the output video

fig = figure('Visible','off');

fig.Units = 'pixels';
fig.Position = [1 1 outputres(1) outputres(2)];
fig.MenuBar = 'none';
fig.ToolBar = 'none';
fig.Color = [0 0 0];
if ~isempty(cdata)
    Nframes = cdata.tvec(end);
    ax1 = axes();
switch cdata.plottype
    case 'video'
    himagesc = imagesc(ax1,cdata.Data(:,:,1));
    ax1.Colormap = cdata.colormap;
    ax1.CLim = cdata.CLim;
    case 'cell_labels'
    himagesc = imagesc(ax1,cdata.Data);
end
    himagesc.XData = cdata.XData;
    himagesc.YData = cdata.YData;
    ax1.XLim = xl;
    ax1.YLim = yl;
    ax1.Units = 'pixels';
    ax1.Position = [1 1 outputres(1) outputres(2)];
    ax1.Color='none';
    ax1.YDir = 'reverse';
    ax1.XColor = 'none';
    ax1.YColor ='none';
    ax1.PlotBoxAspectRatio = [1 1 1];
    ax1.XTick = [];
    ax1.YTick = [];
    hold on
end
if ~isempty(cellboundaries)
    ax3 = axes();
    hbd = imagesc(ax3,cellboundaries.CData);
    hbd.AlphaData = cellboundaries.AlphaData;
    ax3.Colormap = cellboundaries.colormap;
    ax3.CLim = cellboundaries.CLim;
    ax3.XLim = xl;
    ax3.YLim = yl;
    ax3.Units = 'pixels';
    ax3.Position = [1 1 outputres(1) outputres(2)];
    ax3.Color='none';
    ax3.YDir = 'reverse';
    ax3.XColor = 'none';
    ax3.YColor ='none';
    ax3.PlotBoxAspectRatio = [1 1 1];
    ax3.XTick = [];
    ax3.YTick = [];
    hold on
end
if ~isempty(quiverdata)
    Nframes = quiverdata.tvec(end);
    ax2 = axes();
    [Xquiv,Yquiv]=meshgrid(quiverdata.xvec,quiverdata.yvec);
    hquiver = quiver(ax2,Xquiv,Yquiv,...
        quiverdata.U(:,:,1)*quiverdata.scale,...
        quiverdata.V(:,:,1)*quiverdata.scale,...
        0,'color',quiverdata.color,'LineWidth',quiverdata.linewidth);
    ax2.XLim = xl;
    ax2.YLim = yl;
    ax2.Units = 'pixels';
    ax2.Position = [1 1 outputres(1) outputres(2)];
    ax2.Color='none';
    ax2.YDir = 'reverse';
    ax2.XColor = 'none';
    ax2.YColor ='none';
    ax2.PlotBoxAspectRatio = [1 1 1];
    ax2.XTick = [];
    ax2.YTick = [];
end
% linkaxes([ax1, ax2])

%%

[filename,pathname] = uiputfile({'*.mp4'}, 'Save as',fullfile(outputdir, outputname));
v = VideoWriter(fullfile(pathname,filename),'MPEG-4');
v.FrameRate = framerate;
% v.Quality = 100;
open(v);
tic
f = waitbar(0,'Saving Video...');
for iframe = 1 : Nframes
    waitbar(iframe/Nframes,f,'Saving Video...');

    if ~isempty(cdata) && strcmp(cdata.plottype,'video')
        [~,cframe] =  min(abs(iframe-cdata.tvec));
        himagesc.CData = cdata.Data(:,:,cframe);
    end
    if ~isempty(quiverdata)
        [~,qframe] =  min(abs(iframe-quiverdata.tvec));
        hquiver.UData = quiverdata.U(:,:,qframe)*quiverdata.scale;
        hquiver.VData = quiverdata.V(:,:,qframe)*quiverdata.scale;
    end
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(f);
close(v);
toc

