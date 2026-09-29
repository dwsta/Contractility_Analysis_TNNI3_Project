[xvec,yvec,tvec,Sxx,Syy,Sxy,xdrift,ydrift] = readMSM_bin('F:\Playaround_0407data\contractility_run_20240121_170445\output\CyanEx\tfm\monolayer_stresses.bin');
xshift = median(xdrift);
yshift = median(ydrift);

mask = imread('F:\Playaround_0407data\contractility_run_20240121_170445\output\CyanEx\mask\mask_mask_manual.tif');
mask2 = circshift(mask,[yshift,xshift]);
mask3 = mask2;

% If nothing was negative, then no crop from the left
xcrop_left  = max([0 ; xshift]); % Left_crops are negative, flip sign

% Analogous logic for the right
xcrop_right = max([0; -xshift]);

% Negative Y drift -> images need to be cropped from the top
ycrop_top  = max([0 ; yshift]);

% Positive Y drift -> images need to be cropped from the bottom
ycrop_bottom = max([0 ; -yshift]);

mask3(1:ycrop_top+1,1:end)=0;
mask3(ycrop_bottom+1:end,1:end)=0;
mask3(1:end,1:xcrop_left+1)=0;
mask3(1:end,1:xcrop_right+1)=0;

%mask = mask(1 + ycrop_top : end - ycrop_bottom,...
%                             1 + xcrop_left: end - xcrop_right);

imshow(mask,[]);
figure;
imshow(mask2,[]);
figure;
imshow(mask3,[]);