function make_table_barplot(name, count, figname)
% make_table_barplot(name, count, figname)
% make a bar plot with a ui table
% Xiaoyan, 2017


f = figure; subplot(1,9,[4 9]);
bar(1:length(name), count);
set(gca, 'XTick', 1:1:length(name),...
    'XTickLabelRotation', 90, 'box', 'off');
if length(name) <= 20
    set(gca,'XTick', 1:1:length(name),...
        'XTickLabel', name,...
        'XTickLabelRotation', 90, 'FontSize', 9, 'box', 'off');
end
axis([0 length(name)+1 0 max(count)*1.1]);
title(figname);
ylabel('count');
xlabel('read');

h = uitable(f, 'data', [name num2cell(count(:))],...
    'ColumnName', {'name' 'count'},...
    'ColumnFormat', {'char' 'numeric'},...
    'RowName', 1:length(name));
set(h,'units','normalized','Position',[0 0 0.3 1]);
set(f,'name','read distribution','units','normalized',...
    'position',[.2 .2 .5 .5]);

end
