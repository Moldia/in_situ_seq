function iHomomer = findhomomer_read(readlist, expTags, code)
% idxHomomer = findhomomer_read(readlist, expTags, code)
% find homomer (unexpected reads) indices
% Xiaoyan, 2017


if nargin < 3
    code = {'A' 'C' 'G' 'T'};
end

% possible homoemrs
nHybs = length(readlist{1});
homomer = cellfun(@(v) repmat(v, 1, nHybs), code, 'uni', 0);
homomer = homomer(~ismember(homomer, expTags));

[uniRead, ~, idxRead] = unique(readlist);
homomer = cellfun(@(v) find(strcmp(v, uniRead)), homomer, 'uni', 0);
homomer = [homomer{:}];

iHomomer = [];
for i = 1:length(homomer)
    idxhomo = find(idxRead == homomer(i));
    iHomomer = [iHomomer; idxhomo];
end
iHomomer = sort(iHomomer);

end
