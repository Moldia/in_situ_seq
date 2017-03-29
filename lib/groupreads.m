function [uni_seq_res,count_code,uni_letters,uni_name_tag,uni_group_tag,catagory_read] = ...
    groupreads(list,seq_res,letters,expected_list,exp_tags,exp_groups)
% group and classify reads based on expected tags and names
% Xiaoyan, 2014-11-28


% find unique reads
[uni_num,idx_first,~] = unique(list);
% histogram
if length(uni_num)==1
    count_code = length(list);
else
    [count_code,~] = hist(list,uni_num);
end


uni_seq_res = [];
uni_letters = [];
uni_name_tag = cell(length(uni_num),1);
uni_group_tag = cell(length(uni_num),1);

uni_seq_res = seq_res(idx_first,:);
uni_letters = letters(idx_first,:);

uni_name_tag(~ismember(uni_num,expected_list)) = {'NNNN'};
uni_group_tag(~ismember(uni_num,expected_list)) = {'NNNN'};

exp_reads = find(ismember(uni_num,expected_list));
for i = 1:length(exp_reads)
    uni_name_tag(exp_reads(i)) = ...
        exp_tags(uni_num(exp_reads(i))==expected_list);
    uni_group_tag(exp_reads(i)) = ...
        exp_groups(uni_num(exp_reads(i))==expected_list);    
end


catagory_read = zeros(length(list),1);
for i = 1:length(expected_list)
    catagory_read(list==expected_list(i)) = i;
end

