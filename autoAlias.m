function [aliases,chosenLevels] = autoAlias(LocationList)

% Automatically find a combination of fileparts that generate unique
% aliases
data = LocationList2FilePartsTable(LocationList);

% % Example of Jobfile
% jobpath = 'D:\20211213_Stiffness_and_Force\contractility_output\';
% jobname = 'jobfile_ALL.txt';
% 
% exp_list = readJobFile(jobname, jobpath);

% if ispc
%     sep = '\';
% elseif ismac || isunix
%     sep = '/';
% end
% 
% 
% % Reshape exp_list to be in table compatible format
% for iexp = 1:length(exp_list)
%     [C,matches] = strsplit(exp_list{iexp},sep);
%     t{:,iexp} = C(:);
% end
% deepestLevel = max(cell2mat(cellfun(@size,t,'UniformOutput',false))); % Sometimes some folders contain more level of subfolders than others. To create the table we need uniform entries so we'll fill them with empty spaces.
% for iexp = 1:length(t)
%     if length(t{iexp}) < deepestLevel
%         buf = t{iexp};
%         buf(end+1:deepestLevel) = {''};
%         t{iexp} = buf;
%     end
% end

% if a level is the same for all entries, it is useless for
% making unique aliases
uselessLevels = [];
candidateLevels = [];
% data = horzcat(t{:});
for ilevel = 1:size(data,1)
    dum = unique(data(ilevel,:));
    if length(dum) == 1;
        uselessLevels = [uselessLevels, ilevel];
    else
        candidateLevels = [candidateLevels,ilevel];
    end
end

if isempty(candidateLevels)
    candidateLevels = size(data,1)-1;
end
% From the candidate levels, try combinations until each exp_list entry has 
% a unique aliases. Try to use a reduced amount of levels to achieve this.
flag = 0; % Use this flag to break out of nested loop
for nLevels = 1 : length(candidateLevels)
    combinations = nchoosek(candidateLevels,nLevels);
    for icombo = 1 : size(combinations,2)
        curCombo = combinations(icombo,:);
        parts =  data(curCombo,:);
        aliases = {};
        for icol = 1 : size(data,2)
            aliases{icol} = strjoin(parts(:,icol),'_');
        end
        if length(unique(aliases)) == size(data,2)
            chosenLevels = curCombo;
            flag = 1;
            break
        end
    end
    if flag
        break
    end
end
aliases = aliases';