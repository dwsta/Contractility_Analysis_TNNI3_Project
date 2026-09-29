
%
% 
% temp = readtable("D:\Data\20240413_otrin_2185\dws_bluebead384_NOfurared_5s_20240413135121_test\contractility_run_20240413_141420\contractility_traction_20240413T141442.csv")


%temp = readtable("D:\Data\20240413_otrin_2185\dws_bluebead384_NOfurared_5s_20240413135121_test\contractility_run_20240413_141420\output\Well__C_012_BLUE_BEADS\piv\reference_frames.txt");
%t = temp{1,"Var1"};


temp = fileread("D:\Data\20240413_otrin_2185\dws_bluebead384_NOfurared_5s_20240413135121_test\scan\contractility_run_20240413_154034\output\Well__C_012_BLUE_BEADS\whole-ROI_analysis\min_frame.txt");
t = str2double(temp);