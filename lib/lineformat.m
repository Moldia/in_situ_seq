function fmt = lineformat(vartype, repeat)
% line format for csv file writing
% Xiaoyan, 2017


fmt = repmat([vartype, ','], 1, repeat);
fmt = [fmt(1:end-1), '\n'];

end