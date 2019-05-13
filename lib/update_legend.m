function update_legend(ah, plotnameOrdered, symbols, symsize, legendlocation)
% update_legend(ah, plotnameOrdered, symbols, symsize, legendlocation)
% reorder entries in legend and set visibility of line objects according to
% the given names 
% Xiaoyan, 2018

if nargin < 5
    legendlocation = 'NorthEastOutside';
end

if ischar(plotnameOrdered)
    plotnameOrdered = {plotnameOrdered};
end

if nargin>2 && ischar(symbols)
    symbols = {symbols};
end

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
if isempty(plotnameOrdered)
    plotnameOrdered = plotname;
end

for i = 1:length(plotnameOrdered)
    l = find(strcmp(plotnameOrdered{i}, plotname));
    if ~isempty(l)
        h(j) = gc(l);
        h(j).Visible = 'On';
        if nargin > 2
            h(j).Marker = symbols{i}(2);
            h(j).Color = symbols{i}(1);
        end
        if nargin > 3
            h(j).MarkerSize = symsize;
        end
        j=j+1;
    end
end

% new legend
legendnames = plotnameOrdered(ismember(plotnameOrdered, plotname));
legendnames = cellfun(@(v) strrep(v, '\_', '_'), legendnames, 'uni', 0);
legendnames = cellfun(@(v) strrep(v, '_', '\_'), legendnames, 'uni', 0);

% set all others invisible
try
    other_h = setdiff(gc(allChildren), h);
    legend(h(:), legendnames, 'color', [.6 .6 .6], 'location', legendlocation);
catch
    other_h = gc(allChildren);
    legend off;
end
for i = 1:length(other_h)
    other_h(i).Visible = 'Off';
end


end
