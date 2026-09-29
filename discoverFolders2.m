function discovered_folders = discoverFolders2(folder,d)
% folder = '/media/extensiondrive/Wesley/20150609_TCT/';
% Discover images with extension .ext inside the folder.
% List contents inside the folder
if ispc
    sep = '\';
elseif ismac || isunix
    sep = '/';
end

a = dir(folder);
d.Message = ['Looking inside ', folder];
if d.CancelRequested
    discovered_folders = [];
    return
end
filter_out = cellfun('isempty',regexpi({a.name},['\w'],'match','once'));

% Get all the subfolders contained in folder
a(filter_out) = [];

discovered_folders = cellfun(@(c)[folder,c],{a([a.isdir]).name},'uni',false)';

for ifolder = 1:length(a)
    if a(ifolder).isdir
        subfolder = discoverFolders2([folder,a(ifolder).name,sep],d);    
        discovered_folders=[discovered_folders;subfolder];    
    end
end