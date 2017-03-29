function [list,expected_list,unitag,sym] = ...
    groupreadsforplotting(list,expected_list,exp_tags,exp_groups,plot_based_on_group_YN,sym)

% group reads for plotting
% Xiaoyan, 2014-11-29

if plot_based_on_group_YN
    [unitag,~,idx_re] = unique(exp_groups,'stable');
else
    [unitag,~,idx_re] = unique(exp_tags,'stable');
end
for i = 1:max(idx_re)
    dup = find((idx_re==i));
    if length(dup)>1
        for j = 2:length(dup)
            list(list==expected_list(dup(j))) = expected_list(dup(1));
            expected_list(dup(j)) = expected_list(dup(1));
            sym(dup(j)) = sym(dup(1));
        end
    end
end
[expected_list,temp_first_idx,~] = unique(expected_list,'stable');
sym = sym(temp_first_idx);
