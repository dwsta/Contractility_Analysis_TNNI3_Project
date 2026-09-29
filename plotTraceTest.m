filpiv = 'D:\2022-08-31\plate1_30kcells_30min\contractility_run_20220912_093238\output\Well__E_017_MTDeepRed\piv\';
[xvec,yvec,tvec,U,V] = readPIV_bin(fullfile(filpiv,'smooth_deformations_pass3.bin'));
[X,Y] = meshgrid(xvec,yvec);
DIV = divergence_rik(X,Y,U,V);
mag_trace = squeeze(nanmean(abs(DIV),[1,2]));
raw_time = tvec;
raw_signal = mag_trace;
getParameters_phase_avg_4_0(filpiv,raw_time,raw_signal)