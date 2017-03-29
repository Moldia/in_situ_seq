function out=letter2num(letter_list,num_hybs)

% change the letters into digits
% for abnormal sequencing order
% with empty base ('N')
%
% letter2num_abnormal v2.1
% 2016-2-21, Xiaoyan


letters = {'A';'C';'G';'T';'N'};
tag_list = [];
for i = 1:length(letter_list)
    tag = []; 
    tag_temp = letter_list(i);
    tag_temp = tag_temp{1};
    
     for j = 1:num_hybs
        tag_lett = tag_temp(j);
        a = find(strcmp(letters,tag_lett));
        
        if a
            taglist_num(i,j) = a;
            a = num2str(a);
            tag = [tag a];
        else 
            error('Unexpected barcode letter found.');
        end
       
     end
     
     tag_list = [tag_list; tag]; 
end

exp_list = str2num(tag_list);

out = [exp_list taglist_num];