function tableCount = histread(name, varargin)
% plot the number of reads per gene name and ouput a table
% Xiaoyan, 2017

try
    idxHomo = varargin{1};
catch
    idxHomo = [];
end
name(idxHomo) = {'Homomer'};

try
    figname = varargin{2};
catch
    figname = '';
end

% group reads with the same name tag
[uniName, ~, idxName] = unique(name);
idxnameNNNN = find(strcmp(uniName, 'NNNN'));
idxnamehomo = find(strcmp(uniName, 'Homomer'));
orderName = 1:length(uniName);
orderName([idxnameNNNN, idxnamehomo]) = [];
orderName = [orderName, idxnameNNNN, idxnamehomo];

countName = hist(idxName, 1:length(uniName));
countName = countName(orderName);

tableCount = [uniName(orderName) num2cell(countName(:))];

make_table_barplot(uniName(orderName), countName, ['Gene count ' figname]);
 

end