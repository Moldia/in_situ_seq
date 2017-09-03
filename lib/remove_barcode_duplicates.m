function [taglist, map_TagName] = remove_barcode_duplicates(taglist)
% [taglist, map_TagName] = remove_barcode_duplicates(taglist)
% input: 2- or 3-column taglist, check barcode duplicates
% Xiaoyan, 2017


[uniTag, idxtag1, idxTag] = unique(taglist(:,1));
[~, ~, idxName] = unique(taglist(:,2), 'stable');
uniComb = fliplr(unique([idxName, idxTag], 'rows'));

[nNames, ~] = hist(uniComb(:,1), 1:length(uniTag));
if nnz(nNames>1)
    dupTags = find(nNames>1);
    for i = 1:length(dupTags)
        fprintf('Multiple usage of barcode %s found.\n', uniTag{dupTags(i)});
    end
    error('Every barcode should be uniquely mapped to one gene!')
end

taglist = taglist(idxtag1(uniComb(:,1)),:);
map_TagName = idxName(idxtag1(uniComb(:,1)));

end
