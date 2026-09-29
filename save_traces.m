function save_traces(trace_matrix,peakIDs,tvec,fnameout,cfg_data)
% Saves the trace(s) plots. 
% 'trace_matrix' is a matrix of Ntraces rows and Mtimepoints columns (i.e. one trace per row). 
% 'peakIDs' will color the peaks. It can be empty (i.e. []), in which case the
%  resulting plot will onlu have one colorline
% 'frames' is a vector of 1xMtimepoints that contains the time (x-axis) data.
% The plots will be saved in a file at 'fnameout' (provide a full path)
% Provide also the 'cfg_data' to do conversion frame -> time and figure out
% title (Force or Divergence)

if isempty(peakIDs)
    peakIDs = ones(size(trace_matrix));
end

% Default options
nRows = 5; % Number of rows in the tiled layout
nCols = 5; % Number of columns in the tiled layout
line_opts = {'linewidth',1.5};


% How many traces to be plotted
nTraces = size(trace_matrix,1);
    
fig = figure();

%% Whole ROI plot (i.e. Ntraces == 1)
if nTraces == 1
    if cfg_data.TFM.DoTFM
        yname = {'Whole ROI Traction';' Stress RMS (Pa)'};
    else
        yname = {'Whole ROI Divergence';' RMS (dimensionless)'};
    end
    p=plot_multicolor_rik(tvec,trace_matrix,peakIDs,line_opts{:});
    axis([min(tvec) max(tvec) 0 max(trace_matrix(:))*1.05]);
    box off
    set(gca,'FontSize',10)
    set(gca,'linewidth',1.25)
    ylabel(yname)
    xlabel('Time (s)')
    set(fig, 'Units', 'inches');
    set(fig, 'Position', [0 0 4 2]);
    exportgraphics(fig,fnameout,'Resolution',150)
    close(fig)
    return
end

%% Single cell plots (i.e. Ntraces > 1)
if cfg_data.TFM.DoTFM
    yname = 'Single Cell Traction Stress RMS (Pa)';
else
    yname = 'Single Cell Divergence RMS (unitless)';
end


%DWS error handle
try
    xl = [min(tvec) max(tvec)];
    yr = range(trace_matrix(:));
    yl = [min(trace_matrix(:))-0.05*yr max(trace_matrix(:))+0.05*yr];
    nPages = ceil(nTraces/nCols/nRows);
    for iPage = 1:nPages
        %disp("DWS iPage = " + iPage);
        clf(fig)
        set(fig, 'Units', 'inches');
        set(fig, 'Position', [0 0 11 8.5]);
        t = tiledlayout(nRows,nCols);
        t.Units = 'inches';
        t.OuterPosition = [0.25 0.25 10.75 8.25];
        start_trace = 1 + (iPage-1)*nRows*nCols;
        end_trace = min(iPage*nRows*nCols,nTraces);
        for iTile = start_trace:end_trace
            nexttile(t)
            p=plot_multicolor_rik(tvec,trace_matrix(iTile,:),peakIDs(iTile,:),line_opts{:});
            title(['Cell ',num2str(iTile)]);
            axis([xl, yl]);
            box off
            set(gca,'FontSize',10)
            set(gca,'linewidth',1.25)
            ylabel('')
            xlabel('')
        end
        t.XLabel.String = 'Time (s)';
        t.YLabel.String = yname;
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
        set(fig, 'PaperUnits', 'inches');
        set(fig, 'PaperSize', [11 8.5]);
        set(fig, 'PaperPositionMode', 'manual');
        set(fig, 'PaperPosition', [0 0 11 8.5]);
        append_flag = iPage>1; % If it is the first page, this lets the pdf be overwritten, otherwise append the page
        %disp("about to export " + iPage)
        %DWS - this broke for unclear reasons. Outputting one page at a
        %time works
        %exportgraphics(t,fnameout,'Resolution',300,'append',append_flag) %DWS edit, capitalized?
        %make a new filename, insert page number before 3 character
        %filetype
        fnameout_page = insertBefore(fnameout,strlength(fnameout)-3,"_page_" + num2str(iPage));
        %export anew
        exportgraphics(t,fnameout_page,'Resolution',300)
        %disp("Finished DWS iPage = " + iPage);
    end

catch ME
    disp("something went wrong with saving")
    disp(ME)
end

close(fig)


