function padded = paddigits(oriNum, ndigits)
% padded = paddigits(oriNum, ndigits)
% pad a number to given digits (string)
% Xiaoyan, 2017


padded = num2str(oriNum);
while length(padded) < ndigits
    padded = ['0' padded];
end

end
