function alias_list=makeAliasFile(jobname,jobpath)
if ispc
    sep = '\';
elseif ismac || isunix
    sep = '/';
end

exp_list = readJobFile(jobname, jobpath);
for iexp = 1:length(exp_list)
    [C,matches] = strsplit(exp_list{iexp},sep);
    t{:,iexp} = C(:);
end
deepestLevel = max(cell2mat(cellfun(@size,t,'UniformOutput',false))); % Sometimes some folders contain more level of subfolders than others. To create the table we need uniform entries so we'll fill them with empty spaces.
for iexp = 1:length(t)
    if length(t{iexp}) < deepestLevel
        buf = t{iexp};
        buf(end+1:deepestLevel) = {''};
        t{iexp} = buf;
    end
end
f = figure('Position',[100 100 800 600])
dum    = false(size(t{1}));
% dum(:) = ;
data = horzcat(num2cell(dum),t{:});
columnformat   = cell(1,size(t{1},1));
columnformat(:) = {'char'};
columnformat(1) = {'logical'};
columneditable  = false(size(columnformat));
columneditable(1) = true;
u = uitable('Units','normalized','Position',[0.05 0.1 0.9 0.8],...
    'data',data,'ColumnFormat', columnformat,'ColumnEditable',columneditable);
txt_title = uicontrol('Style', 'text','Units','normalized', 'Position', [0.1 0.9 0.9 0.1], 'String', 'Select rows that will generate alias');

show_warning = 'on';
txt_warning = uicontrol('Style', 'text','Units','normalized', 'Position', [0.2 0.01 0.7 0.08],...
    'String', 'WARNING: Check that there is a unique alias for each experiment. Experiments with same alias will result in overwritten files.',...
    'Visible',show_warning);


h = uicontrol(f,'Units','normalized','Position', [0.05 0.01 0.1 0.08], 'String', 'Create Alias File', ...
    'Callback', 'uiresume(gcbf)');

uiwait(gcf);
tableData = get(u, 'data');

name_parts=tableData(cell2mat(tableData(:,1)),2:end)';
for ialias = 1:size(name_parts,1)
    alias_list(ialias) = {strjoin(name_parts(ialias,:),'_')};
end

% if length(unique(alias_list)) >= length(exp_list)
%    	show_warning = 'off';
% end



fid = fopen([jobpath,'Alias_',jobname],'w');
fprintf(fid,'%s\n',alias_list{:});
fclose(fid);
close(gcf)
end