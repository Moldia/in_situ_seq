function [taglist,expected_list,expected_digits,exp_letters,exp_tags,exp_groups] =  ...
    maketag(taglist,num_hybs,cycle5_empty_YN,use_old_style_taglist_YN,...
    grouped_YN,abnormal_sequencing_YN,sequencing_order)
% prepare taglist for decoding
% v 1.3
% Xiaoyan, 2015-5-3


if isempty(taglist)
    error('Please replace the empty taglist');
end
for n = 1:size(taglist,1)
    exp_tags(n) = taglist(n,1);
    [token, remain] = strtok(exp_tags(n));
    exp_letters(n) = token;
end

% convert letters to digits
%--------------------
if abnormal_sequencing_YN==0
    if cycle5_empty_YN==0
        expected_digits = letter2num(exp_letters,num_hybs);
    elseif num_hybs==5 && cycle5_empty_YN
        expected_digits = letter2num(exp_letters,4);
    else
        error('When cycle5 is empty, num_hybs must be 5.');
    end
else
    if num_hybs==5 && cycle5_empty_YN
        expected_digits = letter2num_abnormal(exp_letters,sequencing_order,4);
    else
        expected_digits = letter2num_abnormal(exp_letters,sequencing_order,num_hybs);
    end
end
expected_list = expected_digits(:,1);

% add digit column to taglist, check multiple usage of barcode
%--------------------
taglist = CheckBarcode(expected_list,taglist,use_old_style_taglist_YN);

exp_tags={}; exp_letters={}; exp_groups = cell(size(taglist,1),1);
for n = 1:size(taglist,1)
    exp_tags(n) = taglist(n,1);
    [token,remain] = strtok(exp_tags(n));
    exp_letters(n) = token;
    [token,remain] = strtok(remain);
    exp_tags(n) = token;
    if grouped_YN
        token = strtok(remain);
        exp_groups(n) = token;
    end
end
