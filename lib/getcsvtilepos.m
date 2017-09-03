function startpos = getcsvtilepos(csvfile)
% startpos = getcsvtilepos(csvfile)
% get tile starting positions from CP input file
% Xiaoyan, 2017

fid = fopen(csvfile);
tilepos = textscan(fid, '%s%s%s%*[^\n]', 'headerlines', 1, 'delimiter', ',');
tilepos = cellfun(@str2double, tilepos, 'uni', 0);
tilepos = reshape([tilepos{:}], [], 3);
startpos = unique(tilepos, 'rows');

fclose(fid);

end

