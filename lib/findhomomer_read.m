function idxHomomer = findhomomer_read(readlist, expTags, varargin)
% find homomer (unexpected reads) indeces
% Xiaoyan, 2017


if isempty(varargin)
    code = {'A' 'C' 'G' 'T'};
else
    code = varargin{1};
end

% possible homoemrs
nHybs = length(readlist{1});
homomer = cellfun(@(v) repmat(v, 1, nHybs), code, 'uni', 0);
homomer = homomer(~ismember(homomer, expTags));

[uniRead, ~, idxRead] = unique(readlist);
homomer = cellfun(@(v) find(strcmp(v, uniRead)), homomer, 'uni', 0);
homomer = [homomer{:}];

idxHomomer = [];
for i = 1:length(homomer)
    idxhomo = find(idxRead == homomer(i));
    idxHomomer = [idxHomomer; idxhomo];
end
idxHomomer = sort(idxHomomer);

end