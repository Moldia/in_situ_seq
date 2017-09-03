function [threshold, N] = qualitybar...
    (lower, upper, category, quality, varargin)
% [threshold, N] = qualitybar...
%     (lower, upper, category, quality, varargin)
% draw a bar series to show proportions of 
% expected/unexpected/homomer/below QT at different QTs
% and give a table with read counts
%
% barQuality new
% Xiaoyan, 2014-11-30

threshold = linspace(lower, upper, 21);
inter = (upper-lower)/20;
N = zeros(length(threshold),4);

for i = 1:length(threshold)
    f = category(quality > threshold(i));
    % below Th
    N(i,4) = length(category)-length(f);
    % expected
    N(i,1) = nnz(f~=0 & f~=-1);
    % unexpected
    N(i,2) = nnz(f==0);
    % homomers
    N(i,3) = nnz(f==-1);
end


f = figure; subplot(1, 9, [4 9]);
h = bar(threshold, N, 0.7, 'stacked');
axis([lower-inter upper+inter 0 length(category)*1.02]);
colormap('pink');
set(gca, 'box', 'off','XTick',lower:inter*4:upper);
%set(h,'EdgeColor','none');
set(h, 'EdgeColor', rgb('goldenrod'), 'LineWidth', 1.2);
xlabel('threshold', 'FontSize', 10);
ylabel('frequency', 'FontSize', 10);
title('Reads at different thresholds')
legend({'expected'; 'unexpected'; 'homomer'; 'below QT'});

h = uitable(f, 'data', N, 'ColumnName',...
    {'expected' 'unexpected' 'homomer' 'below QT'},...
    'RowName',threshold);
set(h, 'units', 'normalized', 'Position', [0 0 0.3 1]);
set(f, 'name', 'quality bar', 'units', 'normalized',...
    'position', [.2 .15 .5 .7]);

end
