function out = findHomomer(expected_list,seq_num)

% Find the index for homomer reads (which are not expected) according to 
% the given taglist
%
% findHomomer v2.0
% Mar 25, 2014, Xiaoyan

% m ~ num_hybs
m = length(num2str(expected_list(1,:)));

% create homomer list
Nread = zeros(4,1);
for j = 1:4
    for i = 1:m
        Nread(j,:) = j*10^(i-1)+Nread(j,:);
    end
end

% exclude when homomer reads are one of expected reads
f = ismember(Nread,expected_list);
Nhomo = Nread(~f);

F = [];
if ~isempty(Nhomo)
    for i = 1:length(Nhomo)
        f = find(seq_num==Nhomo(i));
        F = [F;f];
    end
end

out = sort(F);
    
