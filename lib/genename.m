function [name_tag,group_tag] = genename(grouped_YN,allbt,expected_list,exp_tags,exp_groups)
% assign gene names based on reads and expected tags
% Xiaoyan, 2014-11-29


name_tag = cell(length(allbt),1);
group_tag = cell(length(allbt),1);


name_tag(~ismember(allbt,expected_list)) = {'NNNN'};
for i = 1:length(expected_list)
    name_tag(allbt == expected_list(i)) = exp_tags(i);
end

if grouped_YN
    group_tag(ismember(allbt,expected_list)==0) = {'NNNN'};
    for i = 1:length(expected_list)
        group_tag(allbt == expected_list(i)) = exp_groups(i);
    end
    
end
