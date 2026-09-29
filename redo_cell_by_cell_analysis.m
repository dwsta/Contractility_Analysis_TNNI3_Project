rootdir = 'D:\20211213_Stiffness_and_Force\contractility_run_20221130_144728';
% jobname = 'jobfile.csv';
jobname = 'jobfile_ista.csv';
cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));
% cfg_data.TFM.DoTFM = 0;
output_filename= fullfile(rootdir,'results_msm_cbc.csv');
headers = {'alias','cellID','area','centroid_x','centroid_y','majoraxisLength','minoraxisLength','orientation','circularity','perimeter',...
    'signal','time','peak_ID','mean_amplitude','mean_peak_duration', 'mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30',...
    'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'sd_amplitude','sd_peak_duration', 'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', ...
    'sd_peak_value', 'sd_baseline', 'sd_valley'};
output_file_connection = fopen(output_filename,'w');
fprintf(output_file_connection,'%s',strjoin(headers,','));
t = readJobFile2(jobname, rootdir);
for iexp = 1:height(t)
    try
    location = t(iexp,:).Location{1};
    location(1) = 'D';
    alias = t(iexp,:).Alias{1};
    elastography_step(rootdir,alias,cfg_data)
    cell_by_cell_analysis(rootdir,alias,cfg_data,output_file_connection);
    catch ME
        display('Skipping: ')
        t(iexp,:)
    end
end
fclose(output_file_connection)