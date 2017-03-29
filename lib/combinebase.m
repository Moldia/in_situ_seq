function [seq_num,inuni] = combinebase(abnormal_sequencing_YN,num_hybs,cycle5_empty_YN,...
    num_blobs,seq_res,sequencing_order)
% combine bases into reads
% Xiaoyan, 2014-11-28


seq_num = zeros(num_blobs,1); in=[];

if abnormal_sequencing_YN==0
    if cycle5_empty_YN
        p = 4;
    else
        p = num_hybs;
    end
    for b = 1:p
        seq_num = seq_res(:,b).*10^(p-1)+seq_num;
        p = p-1;
        in = [in,b];
    end
else
    if num_hybs==5 && cycle5_empty_YN
        s=4;
    else
        s=num_hybs;
    end
    p = nnz(sequencing_order~='0');
    for b = 1:s
        m = str2num(sequencing_order(b));
        if m
            seq_num = seq_res(:,b).*10^(p-1)+seq_num;
            p = p-1;
            in = [in,b];
        end
    end
end
inuni = unique(in);
