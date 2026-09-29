function t = readJobFile2(jobname, jobpath)
[~,name,ext] = fileparts(jobname);
switch ext
    case '.txt'
        list  = readJobFile(jobname,jobpath);
        alias = readJobFile(['Alias_',jobname],jobpath);
        t = table;
        t.Location = list';
        t.Alias    = alias';
    case '.csv'
        t = readtable(fullfile(jobpath, jobname),'Delimiter',',');
    otherwise
        uialert('Jobfile format not recognized')
end
end