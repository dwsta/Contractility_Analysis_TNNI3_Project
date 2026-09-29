function FilePartsTable = LocationList2FilePartsTable(LocationList)
if ispc
    sep = '\';
elseif ismac || isunix
    sep = '/';
end

for iexp = 1:length(LocationList)
    [C,matches] = strsplit(LocationList{iexp},sep);
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
FilePartsTable = horzcat(t{:});
