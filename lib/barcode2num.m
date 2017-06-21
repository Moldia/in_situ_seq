function [numlist, nummatrix] = barcode2num(letterlist, letters)
% [numlist, nummatrix] = barcode2num(letterlist, letters)
% change the barcode letters to numbers (all same length)
% Xiaoyan, 2017


if nargin < 2
    letters = {'A', 'C', 'G', 'T'};
end

if isempty(letterlist)
    numlist = [];
    nummatrix = [];
else    
    if ~iscell(letterlist)
        letterlist = {letterlist};
    end
    
    letterlist = letterlist(:);
    nHybs = length(letterlist{1});
    
    nummatrix = cell2mat(cellfun(@base2num, letterlist, 'uni', 0));
    numlist = zeros(size(nummatrix,1), 1);
    for i = 1:nHybs
        numlist = numlist + nummatrix(:,i)*10^(nHybs-i);
    end
    
end

    function num = base2num(letter)
        num = [];
        for i = 1:nHybs
            base = letter(i);
            a = find(strcmp(letters, base));
            
            try
                num(i) = a;
            catch
                error('Unexpected barcode letter found.');
            end
        end
    end

end
