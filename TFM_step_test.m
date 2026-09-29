PIV_file = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\output\Well__E_007_Tritc\piv\smooth_deformations_pass4.bin';
E = 8e3;
h = 100;
fcalx = 0.35;
TFM_dir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\contractility_run_20220723_193555_multipass_filter\output\Well__E_007_Tritc\tfm';
TFM_step(PIV_file,E,h,fcalx,TFM_dir);
[xvec,yvec,tvec,U,V] = readPIV_bin(PIV_file);

[xvec,yvec,tvec,TX,TY] = readPIV_bin(fullfile(TFM_dir,'traction_stresses.bin'));
