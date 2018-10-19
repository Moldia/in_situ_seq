function tableCount = histread(name, varargin)
% tableCount = histread(name, varargin)
% plot the number of reads per gene name and ouput a table
% Xiaoyan, 2018

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
[uNames, ~, iName] = unique(name);
iNNNN = find(strcmp(uNames, 'NNNN'));
iHomo = find(strcmp(uNames, 'Homomer'));
orderName = 1:length(uNames);
orderName([iNNNN, iHomo]) = [];
orderName = [orderName, iNNNN, iHomo];

countName = hist(iName, 1:length(uNames));
countName = countName(orderName);

tableCount = [uNames(orderName) num2cell(countName(:))];

make_table_barplot(uNames(orderName), countName, ['Gene count ' figname]);
 

end
