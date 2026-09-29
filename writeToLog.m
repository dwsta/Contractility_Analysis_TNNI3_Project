function writeToLog(logfile,msg)
fid = fopen(logfile,'a');
fprintf(fid,[msg,' ',datestr(now),'\n']);
fclose(fid);
end