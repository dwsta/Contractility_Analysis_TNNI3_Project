% folder = '/media/extensiondrive/Wesley/20150609_TCT/scan/';
function folder_images = discoverImages2(folder,hFig)
if ispc
    sep = '\';
elseif ismac || isunix
    sep = '/';
end

if folder(end)~= sep; folder = [folder,sep]; end

d = uiprogressdlg(hFig,'Title','Please Wait',...
    'Message','Finding Candidate Locations',...
    'Indeterminate','on',...
    'Cancelable','on');

discovered_folders = discoverFolders2(folder,d);
discovered_folders{end+1} = folder;
folder_images = {};
num_images=[];
ext = 'tif';

for ifolder = 1:length(discovered_folders)
    if d.CancelRequested
        folder_images = [];
        return
    end
    a = dir([discovered_folders{ifolder},sep,'*',ext]);
    if ~isempty(a)
        folder_images = [folder_images;discovered_folders{ifolder}];
        num_images    = [num_images;length(a)];
    end
end
