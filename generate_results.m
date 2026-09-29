jobpath = 'D:\20211213_Stiffness_and_Force\contractility_run_20221004_220636\';
jobname = 'jobfile.csv';
jobfile=readJobFile2(jobname,jobpath);
outcsv = fullfile(jobpath,'tractions.csv');

output_file_connection = fopen(outcsv,'w');
headers = {'Alias','T peak','T rise','T fall','TCT',...
    'D peak','D valley','D high','D low','AUC','Power', 'CR', 'RR',...
    'Trace','Time','Phase','Average Trace','dPA_dt','Time for Average'};
fprintf(output_file_connection,'%s\n',strjoin(headers,','));
for iexp = 1 : height(jobfile)
    alias = jobfile.Alias{iexp};
    pivdir = [jobpath,'output\',alias,'\piv'];
    tfmdir = [jobpath,'output\',alias,'\tfm'];
    filtfm = 'traction_stresses.bin';
    ref_frames = load(fullfile(pivdir,'reference_frames.txt'));
    [xvec,yvec,tvec,Tx,Ty] = readPIV_bin(fullfile(tfmdir,filtfm));
    % Remove the measurements at reference frames, as it can introduce errors
    % later in the analysis
    vecind = load(fullfile(pivdir,'vecind.txt'));
    [tf,idx] = ismember(ref_frames,vecind); % This finds the instances we need to remove from tvec, and slices in U and V
    idx(idx==0)=[];
    if ~ isempty(idx)
        tvec(idx)=[];
        Tx(:,:,idx,:)=[];
        Ty(:,:,idx,:)=[];
    end
    raw_time = tvec;
%     [X,Y] = meshgrid(xvec,yvec);
%     div = divergence_rik(X,Y,U,V);
%     raw_signal = squeeze(mean(sqrt(div.^2),[1 2],'omitnan'));
    raw_signal = squeeze(sqrt(mean(Tx.^2+Ty.^2,[1 2],'omitnan')));
    computeParameters_phase_avg(tfmdir,alias,raw_time,raw_signal,output_file_connection)
end
fclose(output_file_connection);
