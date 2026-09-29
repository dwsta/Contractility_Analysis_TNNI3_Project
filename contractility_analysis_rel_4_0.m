% Copyright (c) 2022, Ricardo Serrano
% All rights reserved.

function contractility_analysis_rel_4_0(app)
% Create Time Log
logfile = fullfile(app.OutputDirectoryEditField.Value,['contractility_log_',datestr(now,30),'.log']);
app.Logfile = logfile;
fid = fopen(logfile,'w');
fclose(fid);

rootdir = app.OutputDirectoryEditField.Value;
app.UITableStatus.Data.Status(:) = {'Queued'};
drawnow;
figure(app.UIFigure);
cfg_data = app.Config;

% Make Output Folders for Traces
filoutRAW = fullfile(rootdir,'raw');
if ~exist(filoutRAW,'dir'); mkdir(filoutRAW); end

filoutVEL = fullfile(rootdir,'vel');
if (cfg_data.ReferenceSelection.ReferenceSelectionMode ~= 0)
    if ~exist(filoutVEL,'dir'); mkdir(filoutVEL); end
end

% Make results csv
output_file_connections = create_csv(cfg_data,rootdir);
csv_whole_FOV = output_file_connections(1);
if cfg_data.SingleCell.UseSingleCell
    csv_CBC = output_file_connections(2);
end

%
%
%                         Main Loop
%
%
for iexp = 1:height( app.UITableStatus.Data)
    try
        close all
        
        %if we tried to quit and there was an error that stopped the
        %interrupt, drop out here.
        if app.QuitSignal
            break
        end

        % Define some variables for shorthand
        alias =  app.UITableStatus.Data.Alias{iexp};
        imgpath =  app.UITableStatus.Data.Location{iexp};
        ext = 'tif'; % For now we only work with .tif files

        % Make Output Folder for each Well
        pivdir = fullfile(rootdir,'output',alias,'piv');
        if ~exist(pivdir,'dir'); mkdir(pivdir); end

        %
        %   Load images
        %
        updateStatus(app,iexp,{'Loading Images'});
        writeToLog(app.Logfile,[alias,' Loading Images']);
        [IMAGES, IMAGE_NAMES] = imageLoader(imgpath,ext);
        if app.QuitSignal
            break
        end

        %
        %   Reference Selection
        %
        updateStatus(app, iexp,{'Reference Selection'});
        writeToLog(app.Logfile,[alias,' Reference Selection']);
        [reference_frames, vecind,IMAGES] = reference_selection(IMAGES,imgpath,ext,cfg_data,pivdir,filoutVEL,iexp,alias);
        if app.QuitSignal
            break
        end

        %
        %   Drift correction. Drift is calculated for all frames in the
        %   movie, even if some of those won't be used for computation. The
        %   reason for it is that we want to be able to align the movie
        %   afterwards in the player and will need these drifts.
        %
        if cfg_data.OtherParameters.CorrectImageDrift
            updateStatus(app, iexp,{'Whole-Image Drift Correction'});
            writeToLog(app.Logfile,[alias,' Whole-Image Drift Correction']);
            [IMAGES, DRIFT] = drift_correction(IMAGES,cfg_data,reference_frames,pivdir);
            xdrift = round(squeeze(DRIFT.Xshift)');
            ydrift = round(squeeze(DRIFT.Yshift)');
        else
            xdrift = zeros(1,1,size(IMAGES,3));
            ydrift = zeros(1,1,size(IMAGES,3));            
        end

        %
        %   Multi-pass PIV
        %
        [output] = compute_deformations(app,IMAGES,cfg_data,pivdir,alias,iexp,reference_frames,vecind,xdrift,ydrift);
        if app.QuitSignal
            break
        end
        %
        %   Peak contraction frames
        %
        Ngrids = length(cfg_data.Deformation);
        ptime = squeeze(output(Ngrids).tvec);
        psignal = squeeze(nanmean(sqrt(output(Ngrids).U.^2+output(Ngrids).V.^2),[1 2]))';
        register_peak_contraction_frames(rootdir,alias,ptime,psignal,vecind)
        
        %
        %   Traction Force Microscopy + Elastography
        %
        if cfg_data.TFM.DoTFM
            updateStatus(app,iexp,{'Traction Force Microscopy'})
            writeToLog(app.Logfile,[alias,' Traction Force Microscopy']);
            TFM_step(rootdir,alias,cfg_data)
            elastography_step(rootdir,alias,cfg_data)
        end
        if app.QuitSignal
            break
        end

        %
        %   Data output (whole FOV)
        %
        updateStatus(app,iexp,{'Whole-FOV analysis'})
        writeToLog(app.Logfile,[alias,' Whole-FOV analysis']);
            
        whole_FOV_analysis(rootdir,alias,cfg_data,csv_whole_FOV)


        %
        %   Cell segmentation with Cellpose
        %
        if cfg_data.SingleCell.UseSingleCell
            updateStatus(app,iexp,{'Segmenting Cells'})
            writeToLog(app.Logfile,[alias,' Cell Segmentation']);
            % Move the image to a new folder called 'mask'. Each alias has its own
            imdir = fullfile(rootdir,'output',alias,'mask');
            mkdir(imdir);
            CP_imnames = dir(fullfile(imgpath,cfg_data.SingleCell.Directory,'*.tif'));
            src = fullfile(CP_imnames(1).folder,CP_imnames(1).name);
            if cfg_data.SingleCell.UseCellpose
                dst = fullfile(imdir,['cp_input_',CP_imnames(1).name]);
            else
                dst = fullfile(imdir,['mask_',CP_imnames(1).name]);
            end
            copyfile(src,dst);
            % Main function of the segmentation
            singleCellSegmentation(imdir,cfg_data)
            %
            %   Data output (cell by cell)
            %
            cell_by_cell_analysis(rootdir,alias,cfg_data,csv_CBC)
        end

        if app.QuitSignal
            break
        end

        %DWS add save aligned reference and minimum frame
        minframestring = fileread(fullfile(rootdir,'output',alias,'whole-ROI_analysis','min_frame.txt'));
        minframe = str2num(minframestring);
        mkdir(fullfile(rootdir,'output',alias,'aligned_frames'));
        imwrite(IMAGES(:,:,reference_frames(1)),fullfile(rootdir,'output',alias,'aligned_frames','ref_frame.tif'));
        imwrite(IMAGES(:,:,minframe),fullfile(rootdir,'output',alias,'aligned_frames','min_frame.tif'));
        

        
        %
        %   Finished analysis of this video
        %
        updateStatus(app,iexp,{sprintf('COMPLETED %s',datestr(now))});
        writeToLog(app.Logfile,[alias,' COMPLETED']);

    catch ME
        updateStatus(app,iexp, {sprintf('ERROR %s',datestr(now))});
        writeToLog(app.Logfile,[alias,' ERROR']);

        fid = fopen(fullfile(rootdir,'crashes.log'),'a');
        fprintf(fid,'%s %s \n',alias,datestr(now));

        msgText = getReport(ME,'extended','hyperlinks','off');

        fprintf(fid,'%s',msgText);

        fprintf(fid,'\n\n\n\n');

        fclose(fid);
        close all
    end
end
if app.QuitSignal
    app.UITableStatus.Data.Status(strcmp(app.UITableStatus.Data.Status,{'Queued'})) = {'Stopped'};
    return;
end
fclose all;
end

function [output] = compute_deformations(app,IMAGES,cfg_data,filpiv,alias,iexp,reference_frames,vecind,xdrift,ydrift)
% Multipass PIV with UQ
Ngrids = length(cfg_data.Deformation);
output = struct;
for igrid = 1:Ngrids
    dum_cfg = cfg_data;
    dum_cfg.Deformation = dum_cfg.Deformation(igrid);
    if igrid == 1
        [X0,Y0,T0,U0,V0] = deal([]);
    else
        [X0,Y0,T0] = meshgrid(xvec,yvec,tvec); U0 = nanmean(U,4); V0 = nanmean(V,4);
    end
    %     updateStatus(app,iexp,{sprintf('Image quality pass %s %s',num2str(igrid),datestr(now))});
    %     writeToLog(app.Logfile,[alias,' Image quality pass ',num2str(igrid)]);

    % Compute cross-correlation correction
    %     CIQ = image_quality(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0);

    updateStatus(app,iexp,{sprintf('PIV pass %s %s',num2str(igrid),datestr(now))});
    writeToLog(app.Logfile,[alias,' PIV pass ',num2str(igrid)]);
    [xvec,yvec,tvec,repvec,U,V] = PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecind,X0,Y0,T0,U0,V0,[]);
    
    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    fwrite(fid2,repvec,'single');
    fwrite(fid2,xdrift(vecind),'single');
    fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

    % PIV post-processing with smoothn
    updateStatus(app,iexp,{sprintf('post-processing PIV pass %s %s',num2str(igrid),datestr(now))});
    writeToLog(app.Logfile,[alias,' post-processing PIV pass ',num2str(igrid)]);
    %         [U,V] = filter_PIV(U,V);

    if cfg_data.Performance.UseGPU % CPU -> GPU
        U = gpuArray(U);
        V = gpuArray(V);
    end

    U = smoothn_rik(U,0.001,'robust');
    V = smoothn_rik(V,0.001,'robust');

    if cfg_data.Performance.UseGPU % GPU -> CPU
        U = gather(U);
        V = gather(V);
    end

    output(igrid).xvec = xvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).yvec = yvec + dum_cfg.Deformation.wdw_size/2;
    output(igrid).tvec = tvec;
    output(igrid).repvec = repvec;
    output(igrid).U = U;
    output(igrid).V = V;

    [NY,NX,NF,NR] = size(U);
    fid2 = fopen(fullfile(filpiv,['smooth_deformations_pass',num2str(igrid),'.bin']),'w');
    fwrite(fid2,'v3','uchar');
    fwrite(fid2,[NX,NY,NF,NR],'single');
    fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
    fwrite(fid2,tvec,'single');
    fwrite(fid2,repvec,'single');
    fwrite(fid2,xdrift(vecind),'single');
    fwrite(fid2,ydrift(vecind),'single');
    fwrite(fid2,U,'single');
    fwrite(fid2,V,'single');
    fclose(fid2);

end
end

function output_file_connections = create_csv(cfg_data,rootdir)
% Helper function to create the csv output files and remove clutter from the main loop
% For whole FOV
if cfg_data.TFM.DoTFM
    outcsv = fullfile(rootdir,['contractility_traction_',datestr(now,30),'.csv']);
    headers = {'alias','signal','time','peak_ID','mean_amplitude','mean_peak_duration','mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30', 'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'mean_E','mean_nu','mean_G', 'mean_K', ...
    'sd_amplitude', 'sd_peak_duration' ,'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', 'sd_peak_value', 'sd_baseline', 'sd_valley',...
    'sd_E','sd_nu','sd_G', 'sd_K'};
else
    outcsv = fullfile(rootdir,['contractility_motion_',datestr(now,30),'.csv']);
    headers = {'alias','signal','time','peak_ID','mean_amplitude','mean_peak_duration','mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30', 'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'sd_amplitude', 'sd_peak_duration' ,'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', 'sd_peak_value', 'sd_baseline', 'sd_valley'};

end
output_file_connections(1) = fopen(outcsv,'w');
fprintf(output_file_connections(1),'%s',strjoin(headers,','));
% For cell by cell
if cfg_data.SingleCell.UseSingleCell
    if cfg_data.TFM.DoTFM
        outcsv2 = fullfile(rootdir,['contractility_traction_cell-by-cell_',datestr(now,30),'.csv']);
    else
        outcsv2 = fullfile(rootdir,['contractility_motion_cell-by-cell_',datestr(now,30),'.csv']);
    end
    headers = {'alias','cellID','area','centFOVd_x','centFOVd_y','majoraxisLength','minoraxisLength','orientation','circularity','perimeter',...
    'signal','time','peak_ID','mean_amplitude','mean_peak_duration', 'mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30',...
    'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'sd_amplitude','sd_peak_duration', 'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', ...
    'sd_peak_value', 'sd_baseline', 'sd_valley'};
output_file_connections(2) = fopen(outcsv2,'w');
    fprintf(output_file_connections(2),'%s',strjoin(headers,','));
end
end


function CIQ = image_quality(IMAGES,cfg_data,reference_frames,vecind,X0,Y0,T0,U0,V0)
% Generate images for IQ study
IMAGESIQ = generate_IMAGESIQ(IMAGES,reference_frames);
vecindIQ = 2:size(IMAGESIQ,3);
C0 = [];
[xvec,yvec,tvec,repvec,U,V,~,~,~,~,~,C,~,~] = PIV_multipass_UQ(IMAGESIQ,cfg_data,1,vecindIQ,X0,Y0,T0,U0,V0,C0);
CIQ = repmat(nanmean(C,5),1,1,1,1,length(vecind));
CIQ = reshape(CIQ,cfg_data.Deformation.wdw_size,cfg_data.Deformation.wdw_size,[]);
end

function IMAGESIQ = generate_IMAGESIQ(IMAGES,reference_frames)
IMAGESIQ(:,:,1) = IMAGES(:,:,reference_frames(1));
ntries = 10;
d = 2*pi/ntries;
alp = linspace(0,2*pi-d,ntries);
R = 7;
for iframe = 1 :ntries
    ux = R*cos(alp(iframe));
    uy = R*sin(alp(iframe));
    IMAGESIQ(:,:,iframe+1) = imtranslate(IMAGES(:,:,reference_frames(end)),[ux uy]); % If there are more than 1 reference frames, compare first and last
end
end