function [IMAGES,DRIFT] = drift_correction(IMAGES,cfg_data,reference_frames,filpiv)
[NY,NX,NF] = size(IMAGES);
vecind = 1 : NF; % Calculate the drift for all frames
NR = 1; % Number of UQ repetitions
dum_cfg = cfg_data;
dum_cfg.Deformation = [];
dum_cfg.Deformation.wdw_size = min([NX,NY]);
dum_cfg.Deformation.wdw_spacing = min([NX,NY]);
dum_cfg.Deformation.Nreps = NR;
dum_cfg.Deformation.resPercent = 0;
dum_cfg.Deformation.SubPixelInterpolationMode = cfg_data.Deformation.SubPixelInterpolationMode;

if cfg_data.OtherParameters.UseFeatureDriftCorrection
    [xvec,yvec,tvec,repvec,Ushift,Vshift,angles,newIMAGE] = driftcorrection_Features(IMAGES,cfg_data,reference_frames,vecind);
    IMAGES = newIMAGE;
else
    [xvec,yvec,tvec,repvec,Ushift,Vshift] = PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecind,[],[],[],[],[],[]);
    
    % Align the images
    for indx = 1:length(vecind)
        iframe = vecind(indx); 
        IMAGES(:,:,iframe) = circshift(IMAGES(:,:,iframe),[-round(Vshift(:,:,indx)) -round(Ushift(:,:,indx))]);
    end
end
% Return the calculated drifts (needed to keep track of where the actual
% XY coordinates of the PIV field)
DRIFT.Xshift = round(Ushift);
DRIFT.Yshift = round(Vshift);

% Since the video (variable IMAGES) is aligned we can now crop it.
% Negative X drift -> images need to be cropped from the left
left_crops = squeeze(DRIFT.Xshift(DRIFT.Xshift <0)); 

% If nothing was negative, then no crop from the left
xcrop_left  = max([0 ; -round(left_crops)]); % Left_crops are negative, flip sign

% Analogous logic for the right
right_crops = squeeze(DRIFT.Xshift(DRIFT.Xshift >0)); 
xcrop_right = max([0; round(right_crops)]);

% Negative Y drift -> images need to be cropped from the top
top_crops = squeeze(DRIFT.Yshift(DRIFT.Yshift <0)); 
ycrop_top  = max([0 ; -round(top_crops)]);

% Positive Y drift -> images need to be cropped from the bottom
bottom_crops = squeeze(DRIFT.Yshift(DRIFT.Yshift >0)); 
ycrop_bottom = max([0 ; round(bottom_crops)]);
                                                   
% Crop the images
IMAGES = IMAGES(1 + ycrop_top : end - ycrop_bottom,...
                             1 + xcrop_left: end - xcrop_right ,:);



fid2 = fopen(fullfile(filpiv,'image_drift.bin'),'w');
fwrite(fid2,'v3','uchar');
fwrite(fid2,[length(xvec),length(yvec),length(vecind),length(repvec)],'single');
fwrite(fid2,xvec + dum_cfg.Deformation.wdw_size/2,'single');
fwrite(fid2,yvec + dum_cfg.Deformation.wdw_size/2,'single');
fwrite(fid2,tvec,'single');
fwrite(fid2,repvec,'single');
fwrite(fid2,Ushift,'single');
fwrite(fid2,Vshift,'single');
fclose(fid2);

if cfg_data.OtherParameters.UseFeatureDriftCorrection
    fid3 = fopen(fullfile(filpiv,'angles.bin'),'w');
    fwrite(fid3,angles,'single');
    fclose(fid3);
end

end