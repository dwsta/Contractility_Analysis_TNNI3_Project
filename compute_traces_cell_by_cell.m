function [trace_matrix,frames] = compute_traces_cell_by_cell(cfg_data, rootdir, alias, mask)

% Load TFM/PIV data
if cfg_data.TFM.DoTFM
    pivdir = fullfile(rootdir,'output',alias,'piv');
    tfmdir = fullfile(rootdir,'output',alias,'tfm');
    filmsm = 'monolayer_stresses.bin';
    ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
    %DWS change - reading xdrift and ydrift
    [xvec,yvec,tvec,Sxx,Syy,Sxy,xdrift,ydrift] = readMSM_bin(fullfile(tfmdir,filmsm));
    % Remove the measurements at reference frames, as it can introduce errors
    % later in the analysis
    vecind = load(fullfile(pivdir,'vecind.txt'));
    [tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
    idx(idx ==0)=[];
    if ~isempty(idx)
        tvec(idx)=[];
        Sxx(:,:,idx,:)=[];
        Syy(:,:,idx,:)=[];
    end
    frames = tvec;
    [X,Y] = meshgrid(xvec,yvec);
    raw_signal = Sxx+Syy;
else
    Ngrids = length(cfg_data.Deformation);
    pivdir = fullfile(rootdir,'output',alias,'piv');
    filpiv = ['smooth_deformations_pass',num2str(Ngrids),'.bin'];
    ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
    %DWS change - loaded xdrift and ydrift
    [xvec,yvec,tvec,U,V,xdrift,ydrift] = readPIV_bin(fullfile(pivdir,filpiv));
    % Remove the measurements at reference frames, as it can introduce errors
    % later in the analysis
    vecind = load(fullfile(pivdir,'vecind.txt'));
    [tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
    idx(idx ==0)=[];
    if ~isempty(idx)
        tvec(idx)=[];
        U(:,:,idx,:)=[];
        V(:,:,idx,:)=[];
    end
    frames = tvec;
    [X,Y] = meshgrid(xvec,yvec);
    div = divergence_rik(X,Y,U,V);
    raw_signal = div;
end

%DWS correct drift
%TODO - make conditional - is there a situation where you would not want to
%do this?

%save current mask for display later
displaymask = mask;

%2 - compute average drift (?) vs drift of first frame
% already rounded
xshift = median(xdrift); %median drift should be representative
yshift = median(ydrift); %median drift should be representative
fprintf('xshift = %d yshift = %d\n',xshift,yshift);
%3 - circular shift mask
mask = circshift(mask,[-yshift -xshift]);



%4 - crop mask (see procedure in drift_correction)

% If nothing was negative, then no crop from the left
xcrop_left  = max([0 ; -xshift]); % Left_crops are negative, flip sign

% Analogous logic for the right
xcrop_right = max([0; xshift]);

% Negative Y drift -> images need to be cropped from the top
ycrop_top  = max([0 ; -yshift]);

% Positive Y drift -> images need to be cropped from the bottom
ycrop_bottom = max([0 ; yshift]);
                                                   
% Crop the mask
mask = mask(1 + ycrop_top : end - ycrop_bottom,...
                             1 + xcrop_left: end - xcrop_right);

%creating a display mask that will be compatible with the visualization
%because we didn't shift this cell image, we flip left/right and top/bottom
%and set anything outside to 0
% displaymask(1:ycrop_bottom+1,1:end)=0;
% displaymask(ycrop_top+1:end,1:end)=0;
% displaymask(1:end,1:xcrop_right+1)=0;
% displaymask(1:end,1:xcrop_left+1)=0;

%5 - resize as below

% Resize the mask to match resolution of the vector field
IND = sub2ind(size(mask),round(Y(:)),round(X(:))); % TFM requires an even number of gridpoints, which may cause X and Y to be decimal numbers. 
% Rounding is a good approximation because the resolution of the mask is
% that of the frames of the video (i.e. 2048x2048) vs that of PIV is much
% smaller (63x63) TFM would be (64x64).
mask_resized = reshape(mask(IND),size(X));
% Compute traces cell-by-cell
Ncells = max(mask(:));
[NX,NY,NZ] = size(raw_signal);
raw_signal_r = reshape(raw_signal,[NX*NY],NZ);
for icell = 1:Ncells
    m = mask_resized==icell;
    trace_matrix(icell,:) = squeeze(mean(raw_signal_r(m(:),:),1,'omitnan'));
end
frames = frames';

%DWS - save masks
maskdir = fullfile(rootdir,'output',alias,'mask');
imwrite(mask,fullfile(maskdir,[alias,'_alignedmask.tif']),'tiff');
imwrite(mask_resized,fullfile(maskdir,[alias,'_resizedmask.tif']),'tiff');
imwrite(displaymask,fullfile(maskdir,[alias,'_displaymask.tif']),'tiff');


