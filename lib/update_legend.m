function update_legend(ah, plotnameOrdered, symbols)
% reorder entries in legend and set visibility of line objects according to
% the given names 
% Xiaoyan, 2017


% get all child line objects
gc = get(ah, 'children');
allChildren = [];
for i=1:length(gc)
    if strcmp(gc(i).Type, 'line') 
        allChildren = [allChildren; i];
    end
end
plotname = {gc(allChildren).DisplayName};

% get the ones that are required
clear h
j = 1;
for i = 1:length(plotnameOrdered)
    l = find(strcmp(plotnameOrdered{i}, plotname));
    if ~isempty(l)
        h(j) = gc(l);
        h(j).Visible = 'On';
        if nargin > 2
            h(j).Marker = symbols{i}(2);
            h(j).Color = symbols{i}(1);
        end
        j=j+1;
    end
end

% set all others invisible
other_h = setdiff(gc(allChildren), h);
for i = 1:length(other_h)
    other_h(i).Visible = 'Off';
end

% new legend
legendnames = plotnameOrdered(ismember(plotnameOrdered, plotname));
legend(h(:), legendnames, 'color', [.6 .6 .6]);

end
