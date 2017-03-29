function out = CheckBarcode(expected_list,taglist,use_old)

% Check the multiple usage of the same barcode, and make a new taglist
% If multiple usage exists, gives an error message and terminate.
%
% Checkbarcode v2.1
% Apr 11, 2014, Xiaoyan

[m,n] = size(taglist);
if m==0 || n==0
    error('No information in taglist.');
end
if use_old
    if n==3
        n = 2;
    else
        error('Old styled taglist must have three columns.');
    end
end
taglist_re = cell(length(expected_list),n+1);
if n==2
    if use_old
        for i = 1:length(expected_list)
            taglist_re(i,:) = [taglist(i,1) {expected_list(i)} taglist(i,3)];
        end
    else
        for i = 1:length(expected_list)
            taglist_re(i,:) = [taglist(i,1) {expected_list(i)} taglist(i,2)];
        end
    end
elseif n==1
    for i = 1:length(expected_list)
        taglist_re(i,:) = [taglist(i,1) {expected_list(i)}];
    end
else
    error('Wrong column number in taglist.');
end

[n,m] = hist(expected_list,unique(expected_list));
dup = find(n>1);
del = [];
for i = 1:length(dup)
    dup_index = find(expected_list == m(dup(i)));
    tagerror = ['Multiple usage of barcode '...
        '(transformed based on sequencing order) '...
        num2letters(m(dup(i))) ' is found. '];
    error(tagerror);
    %del = [del;dup_index(2:end)];
end

%taglist_new = removerows(taglist_re,'ind',del);
%out = taglist_new;
out = taglist_re;
end

