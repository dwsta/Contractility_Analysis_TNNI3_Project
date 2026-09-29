function p = plot_multicolor_rik(x,y,group,varargin)
group_keys = unique(group);
nGroups = length(group_keys);
colorgroup = lines(nGroups);
hold on
for iGroup = 1:nGroups
    group_elements = find(group==group_keys(iGroup));
    if iGroup<nGroups
        %DWS EDIT - got rid of the +1 for now
        group_elements(end) = group_elements(end)+0; % This adds overlap between a group and the next one
    end

    %DWSEdit
    if group_keys(iGroup) ~= 0
        p = plot(x(group_elements),y(group_elements),varargin{:});
    else
        p = plot(x(group_elements),y(group_elements),'.');
        %p = plot(x(group_elements),y(group_elements),varargin{:}); %DWS bug fix temp
    end
    p.Color = colorgroup(iGroup,:);
end
hold off
% drawnow
% 
% nGroups = length(unique(group));
% colorgroup = [uint8(lines(nGroups)*255) uint8(ones(nGroups,1))].';
% colorpoint = colorgroup(:,group);
% set(p.Edge,'ColorBinding','interpolated', 'ColorData',colorpoint)

