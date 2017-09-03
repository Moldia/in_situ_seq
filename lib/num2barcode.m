function barcode = num2barcode(numlist, varargin)
% barcode = num2barcode(numlist, varargin)
% change list of numbers to list of barcode letters (all same length)
% default letter list: A, C, G, T, custom list can be used
% Xiaoyan, 2017


if isempty(varargin)
    letters = {'A','C','G','T'};
else
    letters = varargin{1};
end

nHybs = length(num2str(numlist(1)));

nummatrix = zeros(length(numlist), nHybs);
for i = 1:nHybs
    nummatrix(:,nHybs-i+1) = mod(numlist, 10);
    numlist = floor(numlist/10);
end

barcode = cell(length(numlist),1);
for i = 1:nHybs
    try
        barcode = strcat(barcode,...
            cellfun(@(v) letters{v}, num2cell(nummatrix(:,i)), 'uni', 0));
    catch
        error('Barcode number exceeds possible letters.');
    end
end

end

