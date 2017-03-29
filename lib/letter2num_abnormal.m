function out=letter2num_abnormal(letter_list,order,num_hybs)

% change the letters into digits
% for abnormal sequencing order
% with empty base ('N')
%
% letter2num_abnormal v2.1
% last update, 2016-5-6, Xiaoyan


letters = {'A';'C';'G';'T';'N'};

if length(order)~=num_hybs
    error('Abnormal sequencing: sequencing length and num_hybs do not match.');
elseif isempty(order(order~='0'))
    error('Abnormal sequencing: sequencing length is zero.');
end

tag_list = [];
for i = 1:length(letter_list)
    tag = []; 
    tag_temp = letter_list(i);
    tag_temp = tag_temp{1};
    
    p = 1;
    
     for m =1:length(order(1,:))
         j = order(1,m);
         if str2num(j)
             tag_lett = tag_temp(str2num(j));
             a = find(strcmp(letters,tag_lett));
             
             if a
                 taglist_num(i,p) = a;
                 a = num2str(a);
                 tag = [tag a];
             else
                 error('Unexpected barcode letter found.');
             end
             p = p+1;
         end
         
     end
     
     tag_list = [tag_list; tag];
end

exp_list = str2num(tag_list);

out = [exp_list taglist_num];
end