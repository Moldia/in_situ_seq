function out = histRead_group(list,exp_tags,exp_groups,...
    expected_list,name,homomer_index)

% Plot the number of reads per gene name, give a table, with gene group
% also plot the number of reads per gene group, with a table
%
% new histRead_group
% Xiaoyan, 2016-1-25

[unitag,~,re_idx] = unique(exp_tags,'stable');
% check name replicates
[a,~] = hist(re_idx,1:length(unitag));
if nnz(a>1)
    b = find(a>1);
    for i = 1:length(b)
        dup = find(re_idx==b(i));
        
        % change the seq num of name replicates
        for j = 2:length(dup)
            list(list==expected_list(dup(j))) = expected_list(dup(1));
        end
        expected_list(re_idx==b(i)) = expected_list(dup(1));
        exp_groups(re_idx==b(i)) = exp_groups(dup(1));
    end
    expected_list = unique(expected_list,'stable');
end

% give each read a "plot number" based on gene name/homomer/unexpected
plot_num = zeros(length(list),1);
% expected reads
plot_num(ismember(list,expected_list)) = list(ismember(list,expected_list));
% homomers
plot_num(homomer_index) = 666666;
% other unexpected
plot_num(plot_num==0) = 555555;
 
 plot_Xaxis = [expected_list;555555;666666];
 plot_Xaxis_sorted = sort(plot_Xaxis);
 
 for k = 1:length(plot_Xaxis)
     a = find(expected_list == plot_Xaxis(k));
     if a
         plot_label{k} = unitag{a};
     elseif plot_Xaxis(k)==555555
         plot_label{k} = 'NNNN';
     elseif plot_Xaxis(k)==666666
          plot_label{k} = 'homomer';
     end
 end

 for k = 1:length(plot_Xaxis)
     a = find(expected_list == plot_Xaxis(k));
     if a
         plot_group{k} = exp_groups{a};
     elseif plot_Xaxis(k)==555555
         plot_group{k} = 'NNNN';
     elseif plot_Xaxis(k)==666666
         plot_group{k} = 'homomer';
     end
 end
  
 [cc,dd] = hist(plot_num,plot_Xaxis_sorted);
 
 % create a figure and make a bar graph
 %---------------------------
 f=figure; subplot(1,9,[4 9]), bar(1:length(dd),cc);
 set(gca,'XTick',1:1:length(dd),'XTickLabelRotation',90,'box','off');
 if length(dd) <= 20
    set(gca,'XTick',1:1:length(dd),'XTickLabel',plot_label,...
        'XTickLabelRotation',90,'FontSize',9,'box','off'); 
 end
 axis([0 length(dd)+1 0 max(cc)*1.1]);
 title(['Distribution of reads ' name]);
 ylabel('frequency');
 xlabel('read');
 
 h = uitable(f,'data',[plot_label(:) num2cell(cc(:)) plot_group(:)],...
     'ColumnName',{'gene name' 'frequency' 'gene group'},...
     'ColumnFormat',{'char' 'numeric' 'char'},...
     'RowName',1:length(dd));
 set(h,'units','normalized','Position',[0 0 0.3 1]);
 out = [plot_label(:) num2cell(cc(:)) plot_group(:)];
 
 set(f,'name','read distribution','units','normalized',...
         'position',[.2 .2 .5 .5]);
end
