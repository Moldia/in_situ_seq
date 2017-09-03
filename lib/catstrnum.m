function strnum = catstrnum(str, numlist, varargin)
% strnum = catstrnum(str, numlist, varargin)
% concatenate string and a list of numbers
% Xiaoyan, 2017


if isempty(varargin)
    numstr = cellfun(@num2str, num2cell(numlist), 'uni', 0);
else
    numstr = cellfun(@(v) paddigits(v, varargin{1}), num2cell(numlist), 'uni', 0);
end
strnum = strcat(str, numstr);

end

function padded = paddigits(oriNum, ndigits)

padded = num2str(oriNum);
while length(padded) < ndigits
    padded = ['0' padded];
end

end
