function guessRead(expect_matrix,list,name_tag_list,exp_tags,output_prefix)

% Map the unexpected reads to the closed expected one(s)
% Allow only one mismatch
% Give a uitable as output
% Output a csv file
%
% guessRead v2.2
% May 12, Xiaoyan


% extract unexpected reads, and make a digit matrix of unique reads
%----------------------------
NNNN  = list(strcmp(name_tag_list,'NNNN'));
[numN,uni_N] = hist(NNNN,unique(NNNN));

uni_N_str = num2str(uni_N);
N_matrix = [];

for j = 1:length(uni_N_str)
    N_digits = [];
    for i = 1:length(expect_matrix(1,:))
        N_digits = [N_digits str2num(uni_N_str(j,i))];
    end
    N_matrix = [N_matrix; N_digits];
end

% "sequence alignment"
%----------------------------
for kk = 1: length(N_matrix(:,1))
    N_read = N_matrix(kk,:);
    for j = 1:length(expect_matrix(:,1))
        exp_read = expect_matrix(j,:);
        x = 0;
        for i = 1:length(N_matrix(1,:))
            if N_read(i) == exp_read(i)
                x = x+1;
            end
        end
        Ascore(kk,j) = x;
    end
end

% find one mismatch
%----------------------------
f = Ascore(:,:)==length(expect_matrix(1,:))-1;
show=cell(length(uni_N),length(expect_matrix(:,1)));
show(f)={'match'};
show(~f)={''};
fi=figure;
for i = 1:length(numN)
    rown{i}= num2letters(uni_N(i));
end
columnwide=[];
for i = 1:length(expect_matrix(:,1))
    columnwide = [columnwide, num2cell(80)];
end
column = ['count'; exp_tags];
h = uitable(fi,'data',[num2cell(numN(:)) show],'ColumnName',column,'RowName',rown);
set(h,'units','normalized','ColumnWidth',columnwide,'Position',[0 0 1 1]);
set(fi,'name','guess reads');


fid = fopen([output_prefix 'GuessRead' '.csv'], 'w');
fprintf(fid,'%s\n','This is a table showing the guessing of reads ');
fprintf(fid,'');
for j = 1:length(column)
    fprintf(fid,',%s',column{j});
end
fprintf(fid,'\n');
for i = 1:length(rown)
    fprintf(fid,'%s,%d',rown{i},numN(i));
    for j = 1:length(column)-1
        fprintf(fid,',%s',show{i,j});
    end
    fprintf(fid,'\n');
end
fclose(fid);

end

