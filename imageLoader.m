function [IMAGES,IMAGE_NAMES] = imageLoader(imgpath,ext,varargin)

% For now we only accept videos as image sequences
buf = ['*',ext];
IMAGE_NAMES = dir(fullfile(imgpath,buf));
IMAGE_NAMES = {IMAGE_NAMES.name};

N_FRAMES = length(IMAGE_NAMES);

info = imfinfo(fullfile(imgpath,IMAGE_NAMES{1}));
IMAGES = zeros(info.Height,info.Width,N_FRAMES,['uint',num2str(info.BitDepth)]);

try 
    nopar = strcmp(varargin(1),'noparallel'); 
catch
    nopar = 0;
end

if nopar
    for iframe = 1 : N_FRAMES
        m=memmapfile(fullfile(imgpath,IMAGE_NAMES{iframe}),'Offset',info.StripOffsets,'Format',{['uint',num2str(info.BitDepth)],[info.Width,info.Height],'x'},'Writable',false);
        IMAGES(:,:,iframe) = m.Data.x';
    end

else
    % Load whole video into memory
    parfor iframe = 1 : N_FRAMES
        m=memmapfile(fullfile(imgpath,IMAGE_NAMES{iframe}),'Offset',info.StripOffsets,'Format',{['uint',num2str(info.BitDepth)],[info.Width,info.Height],'x'},'Writable',false);
        IMAGES(:,:,iframe) = m.Data.x';
    end
end
end