function [aligned_images] = align_images(IMAGES,drift_file)
[Xdrift,Ydrift] = readDrift_bin(drift_file);
% Align the images
for indx = 1:size(IMAGES,3)
    IMAGES(:,:,indx) = circshift(IMAGES(:,:,indx),[-round(Ydrift(indx)) -round(Xdrift(indx))]);
end


% Since the video (variable IMAGES) is aligned we can now crop it.
% Negative X drift -> images need to be cropped from the left
left_crops = squeeze(Xdrift(Xdrift <0)); 

% If nothing was negative, then no crop from the left
xcrop_left  = max([0 ; -round(left_crops)]); % Left_crops are negative, flip sign

% Analogous logic for the right
right_crops = squeeze(Xdrift(Xdrift >0)); 
xcrop_right = max([0; round(right_crops)]);

% Negative Y drift -> images need to be cropped from the top
top_crops = squeeze(Ydrift(Ydrift<0)); 
ycrop_top  = max([0 ; -round(top_crops)]);

% Positive Y drift -> images need to be cropped from the bottom
bottom_crops = squeeze(Ydrift(Ydrift >0)); 
ycrop_bottom = max([0 ; round(bottom_crops)]);
                                                   
% Crop the images
aligned_images = IMAGES(1 + ycrop_top : end - ycrop_bottom,...
                             1 + xcrop_left: end - xcrop_right ,:);