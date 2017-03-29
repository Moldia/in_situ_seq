function N = tableGene_2(threshold,seq_quality_min,catagory,exp_tags,output_prefix)

% Make a table of read counts based on gene names given different
% thresholds
% Output a csv file
% 
% tableGene_2
% Xiaoyan, 2014-11-29

b = -1:max(catagory);
b = b';
N = zeros(length(b)+3,length(threshold));
for i = 1:length(threshold)
    f = seq_quality_min>threshold(i);
    % below Th
    empty = length(find(~f));
    f = catagory(f);
    [a,~] = hist(f,b);
    a = a';
    N(:,i) = [a(3:end);sum(a(3:end));a(2);a(1);empty;sum(a)+empty];
end

rown = [exp_tags';'total expected';'NNNN';'homomer';'below QT';'TOTAL'];

f = figure;
h = uitable(f,'data',N,'ColumnName',threshold,'RowName',rown);
set(h,'units','normalized','Position',[0 0 1 1]);
set(f,'name','read distribution table');

datawrite = [cellstr(['threshold';rown]),num2cell([threshold;N])]';
fid = fopen([output_prefix 'ReadTable' '.csv'], 'w');
fprintf(fid,'%s\n',['This is a table showing the read counts at different thresholds ',...
    '(defined by user)']);
temp = '%s';
for i = 1:length(threshold)
    temp = [temp ',%d'];
end
temp = [temp '\n'];
fprintf(fid,temp,datawrite{:});
fclose(fid);

