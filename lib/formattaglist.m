function [taglist, isRandSym] = formattaglist(taglist)
% [taglist, isRandSym] = formattaglist(taglist)
% input: old .m style taglist, output: 3-column taglist with symbols
% Xiaoyan, 2017

try
    splitedtag = cellfun(@strsplit, taglist(:,1), 'uni', 0 );
    try
        symbols = taglist(:,2);
        isRandSym = 0;
    catch
        symbols = repmat(symlist,...
            ceil(size(taglist,1)/length(symlist)), 1);
        symbols = symbols(1:length(taglist));
        isRandSym = 1;
    end

    taglist = cell(length(splitedtag), 3);
    for i = 1:length(splitedtag) 
        taglist(i,:) = [splitedtag{i}, symbols(i)];
    end
catch
    taglist = cell(0,3);
end
    