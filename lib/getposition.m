function tile_start_pos = getposition(csv_file_contain_tile_position)
% get tile starting positions from CP input file
% Xiaoyan 2014-11-29
% depricated

fID = fopen(csv_file_contain_tile_position);
tile_pos_file = textscan(fID, '%s%s%s%*[^\n]','headerlines',1,'delimiter',',');
python_tile = tile_pos_file{1};
python_x = tile_pos_file{2};
python_y = tile_pos_file{3};
for l = 1:length(python_tile)
    tile_start_pos(l,1) = str2double(python_tile{l});
    tile_start_pos(l,2) = str2double(python_x{l});
    tile_start_pos(l,3) = str2double(python_y{l});
end
fclose(fID);

end

