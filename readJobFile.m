function list = readJobFile(jobname, jobpath)
fid = fopen([jobpath,jobname]);
tline = fgetl(fid);
list{1} = tline;
iii=1;
while ischar(tline)
    iii=iii+1;
    tline = fgetl(fid);
    list{iii} = tline;
end
list(end)=[];
fclose(fid);
end