function m_letters = num2letters_2(m_num)

% Convert digits to letters
% New order of letters
% N and O for multiple anchor sequencing special cases
% Input sequencing result matrix
%
% new num2letters
% Sep 4, 2014, Xiaoyan

letters = ['A';'C';'G';'T';'N';'O'];
m_letters = repmat(char(0),size(m_num,1),size(m_num,2));

for i = 1:size(m_num,1)
    for b=1:size(m_num,2)
        m_letters(i,b) = letters(m_num(i,b));
    end
end

