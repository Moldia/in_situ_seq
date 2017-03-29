function s_letters = num2letters(numlist)

% Convert digits to letters
% New order of letters
% N and O for multiple anchor sequencing special cases
%
% num2letters v2.1
% May 16, 2014, Xiaoyan

letters = ['A';'C';'G';'T';'N';'O'];
b = num2str(numlist);
s_letters=[];
for i = 1:length(b)
    s_letters = [s_letters letters(str2num(b(i)))];
end

