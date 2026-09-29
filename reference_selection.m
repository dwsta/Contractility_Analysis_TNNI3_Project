function [reference_frames,vecind,IMAGES] = reference_selection(IMAGES,imgpath,ext,cfg_data,filpiv,filoutVEL,iexp,alias)
switch cfg_data.ReferenceSelection.ReferenceSelectionMode
    case 0        % Manual mode
        A = readlines(cfg_data.ReferenceSelection.manualReferenceFile);
        iexpline = A(~cellfun(@isempty,regexp(A,alias)));
        iexpline = strsplit(iexpline,',');
        aliasCheck = iexpline{1};
        if strcmp(alias,aliasCheck)
            Nframes = size(IMAGES,3);
            vecind = [1:Nframes];
            buf = textscan(iexpline{2},'%s','delimiter',';');
            buf = buf{1};
            Nrefs = length(buf);
            reference_frames = [(Nframes+1):(Nframes+Nrefs)];
            for iref = 1:Nrefs
                IMAGES(:,:,Nframes+iref) = imread(buf{iref});
            end
        else
            error('Mismatch between Reference frame file and jobfile.');
        end
    case 3 % Different folder
        Nframes = size(IMAGES,3);
        vecind = [1:Nframes];
        ref_imgpath = [fullfile(imgpath,cfg_data.ReferenceSelection.ReferenceFolderLocation),filesep];
        [refIMAGES, refIMAGE_NAMES] = imageLoader(ref_imgpath,ext);
        Nrefs = size(refIMAGES,3);
        reference_frames = [(Nframes+1):(Nframes+Nrefs)];
        IMAGES(:,:,Nframes+1:Nframes+Nrefs) = refIMAGES;
    case 2   % Frame to frame
        try
            avg_mag_vel = double(squeeze(mean(mean(abs(diff(single(IMAGES),1,3)),1),2)))';
        catch ME
            %DWS added because of memory problems with the "single" call that duplicates the stack. If this fails I don't
            %know what to do
            avg_mag_vel = double(squeeze(mean(mean(abs(diff(IMAGES,1,3)),1),2)))';
        end
    
    case 1   % PIV
        dum_cfg = cfg_data;
        dum_cfg.Deformation = cfg_data.ReferenceSelection;
        NF = size(IMAGES,3);
        for iii = 1 : NF-1
            reference_frames = iii;
            vecind = iii+1;
            [xvec,yvec,tvec,repvec,U,V] = PIV_multipass_UQ(IMAGES,dum_cfg,reference_frames,vecind,[],[],[],[],[],[]);
            avg_mag_vel(iii)= mean(sqrt(U(:).^2+V(:).^2),'omitnan');
        end
        avg_mag_vel(NF)=avg_mag_vel(NF-1);
end
% Skip frames with low contraction
if (cfg_data.ReferenceSelection.ReferenceSelectionMode == 1 || ...
        cfg_data.ReferenceSelection.ReferenceSelectionMode == 2)
    % To find the reference frames, find points whose neighbors are also
    % low movement. For that, perform a moving maximum and then find the
    % minimum k points, where k will be the number of reference frames
    % requested by the user
    stencil = 7;
    velmovmax = movmax(avg_mag_vel,stencil);
    [~,reference_frames] =  mink(velmovmax,cfg_data.ReferenceSelection.Nref_frames);

    valley_vel = quantile(avg_mag_vel,0.05);
    range_vel  = range(avg_mag_vel);

    resting = velmovmax < (0.05*range_vel +valley_vel);
    skip = resting;
    N_FRAMES = length(avg_mag_vel);
    stencil2 = 5;
    for iii = 1 : N_FRAMES
        if iii < stencil2     % Do not skip first frames
            skip(iii) = 0;
            continue
        elseif iii > N_FRAMES-stencil2  % Do not skip last frames
            skip(iii) = 0;
            continue
        end
        instencil = iii - floor(stencil2/2) : iii + floor(stencil2/2);
        % Remove out of stencil points
        instencil(instencil<1) =[];
        instencil(instencil>N_FRAMES-1) =[];
        % If frames before and after are resting and all are being skipped,
        % then don't skip this one
        if all(resting(instencil)) && all(skip(instencil))
            skip(iii) = 0;
        end
    end
    vecind = find(~skip);
    % Print automatic reference selection signal
    printAutoselection(filpiv,reference_frames,vecind,avg_mag_vel);
    copyfile(fullfile(filpiv,'autoselection.png'),fullfile(filoutVEL,[alias,'.png']));

end

% Save reference frames and vector of frames that will be analyzed with PIV (vecind)
fid2 = fopen(fullfile(filpiv,'reference_frames.txt'),'w');
fprintf(fid2,repmat('%3d ',[1 length(reference_frames)]),reference_frames);
fclose(fid2);
fid = fopen(fullfile(filpiv,'vecind.txt'),'w');
fprintf(fid,'%3d ',vecind');
fclose(fid);

end
function printAutoselection(filpiv,reference_frames,vecind,avg_mag_vel)
% OUTPUT
hf = figure;clf
set(0, 'CurrentFigure', hf)
plot(1:length(avg_mag_vel),avg_mag_vel)
hold on
plot(vecind,median(avg_mag_vel)*ones(1,length(vecind)),'r.')
plot(reference_frames,mean(avg_mag_vel)*ones(1,length(reference_frames)),'k.')
h = legend('$\bar{|vel|}$','Subset for PIV','Reference frames');
xlabel('Time (frames)')
set(0, 'CurrentFigure', hf)
set(h,'Interpreter','latex')
set(hf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4096 1024]/200);
% print('-dpng','-r200',[filpiv,'autoselection.png']);
saveas(hf,fullfile(filpiv,'autoselection'),'png')
fid = fopen(fullfile(filpiv,'reference_frames.txt'),'w');
fprintf(fid,'%3d %3d %3d %3d',reference_frames);
fclose(fid);
fid = fopen(fullfile(filpiv,'vecind.txt'),'w');
fprintf(fid,'%3d ',vecind');
fclose(fid);
fid = fopen(fullfile(filpiv,'velocity.txt'),'w');
fprintf(fid,' %7.5f %7.9f \n',[1:length(avg_mag_vel)-1;avg_mag_vel(2:end)]);
fclose(fid);
close(hf)

end
