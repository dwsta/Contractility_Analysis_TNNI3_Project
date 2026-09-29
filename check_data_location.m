function [videoFound,VideoDir] = check_data_location(VideoDir)
% The drive letter may change if the drive is mounted/re-mounted. First,
% check that the data exist in the attempted location, recorded in the
% jobfile. If it doesn't check other drives for possible matches.
if ~exist(VideoDir,'dir')
    drives = erase(getdrives,':\');
    dum = strsplit(VideoDir,':\');
    strippedDir = dum{2};
    for ind=1:length(drives)
        cur_drive = drives{ind};
        cur_dir   = [cur_drive,':\',strippedDir];
        if exist(cur_dir,'dir')
            videoFound = 2;
            VideoDir = cur_dir;
            return
        end
    end
    videoFound = 0;
    VideoDir = [];
else
    videoFound = 1;
end