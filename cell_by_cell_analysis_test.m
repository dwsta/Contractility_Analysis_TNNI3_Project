rootdir = 'G:\Shared drives\Stanford UCSD_sharedDrive\contractility_analysis_rel_4_0\TaxolFullMovies\Pretreatment_Middle\contractility_run_20221122_211235';
alias = 'Fitc_2';
cfg_data = loadJsonConfig(fullfile(rootdir,'config.json'));

output_filename= 'results_cbc.csv';
headers = {'alias','cellID','area','centroid_x','centroid_y','orientation','circularity','majoraxisLength','minoraxisLength','perimeter','signal','time','peak_ID','mean_amplitude', 'mean_rise_time', 'mean_fall_time', 'mean_pw90', 'mean_pw50', 'mean_pw30', 'mean_peak_value', 'mean_baseline', 'mean_valley',...
    'sd_amplitude', 'sd_rise_time', 'sd_fall_time', 'sd_pw90', 'sd_pw50', 'sd_pw30', 'sd_peak_value', 'sd_baseline', 'sd_valley'};
output_file_connection = fopen(output_filename,'w');
fprintf(output_file_connection,'%s',strjoin(headers,','));
cell_by_cell_analysis(rootdir,alias,cfg_data,output_file_connection)
fclose(output_file_connection)