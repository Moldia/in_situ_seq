function iHomomer = findhomomer_num(numlist, expectednum, varargin)
% idxHomomer = findhomomer_num(numlist, expectednum, varargin)
% find homomer (unexpected reads) indices
% Xiaoyan, 2017

nHybs = length(num2str(numlist(1)));

if isempty(varargin)
    nDiffBases = 4;
else
    nDiffBases = varargin{1};
end

% possible homomers
homomatrix = repmat((1:4)', 1, nHybs);
homomer = zeros(nDiffBases,1);
for i = 1:nHybs
    homomer = homomer + homomatrix(:,i)*10^(nHybs-i);
end

% exclude when homomer reads are one of expected reads
homomer = homomer(~ismember(homomer, expectednum));

iHomomer = [];
for i = 1:length(homomer)
    idxhomo = find(numlist == homomer(i));
    iHomomer = [iHomomer; idxhomo];
end
iHomomer = sort(iHomomer);

end
