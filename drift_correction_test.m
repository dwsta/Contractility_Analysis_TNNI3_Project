%%

imgpath = 'G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\scan\B003_r0004_c0005\beads_interleaved\';
cfg_data = loadJsonConfig('G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\scan\contractility_run_20230319_162349\config.json');
ext = 'tif';
[IMAGES, IMAGE_NAMES] = imageLoader(imgpath,ext,'noparallel');
pivdir = 'G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\scan\contractility_run_20230319_162349\output\B003_r0004_c0005_beads_interleaved\piv\';
iexp = 1;
alias = 'TEST';
filoutVEL = 'G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\scan\contractility_run_20230317_151532\output\B003_r0004_c0005_beads_interleaved\vel\';

%%
[reference_frames, vecind,IMAGES] = reference_selection(IMAGES,imgpath,ext,cfg_data,pivdir,filoutVEL,iexp,alias);


%% Drift calculation
[IMAGES_dedrift,DRIFT] = drift_correction(IMAGES,cfg_data,reference_frames,vecind,pivdir);

%%
figure(1)
axis equal
ax = imshow(IMAGES_dedrift(:,:,1));
caxis([0 1]*3000)

for iframe = 1:size(IMAGES_dedrift,3)
    set(ax,'CData',IMAGES_dedrift(:,:,iframe));
    title(num2str(iframe))
    pause(0.01)
end
%% 
drift_file = 'G:\Shared drives\David_Ricardo_TFM_codes\TFM_DEMO_TEST\scan\contractility_run_20230319_162349\output\B003_r0004_c0005_beads_interleaved\piv\image_drift.bin';
[IMAGES_aligned] = align_images(IMAGES,drift_file);

figure(2)
axis equal
ax = imshow(IMAGES_aligned(:,:,1));
caxis([0 1]*3000)

for iframe = 1:size(IMAGES_aligned,3)
    set(ax,'CData',IMAGES_aligned(:,:,iframe));
    title(num2str(iframe))
    pause(0.01)
end