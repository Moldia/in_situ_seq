function qthreshold(decodedir, input)
% qthreshold(decodedir, input)
% quality thresholding for in situ sequencing
% input from InSituSequencing script
% Xiaoyan, 2018

drawnow

disp('loading files..')
decodedata = importdata(fullfile(decodedir, 'beforeQT_details.csv'));
header = decodedata.textdata(1,:);
textdata = decodedata.textdata(2:end, 1:2);
decodedata = decodedata.data;

% thresholding
headerq = find(~cellfun(@isempty, strfind(lower(header), 'quality')));
headera = find(...
    ~cellfun(@isempty, strfind(lower(header), 'general')) |...
    ~cellfun(@isempty, strfind(lower(header), 'anchor')));
toinclude = decodedata(:,headerq-2) > input.quality_threshold &...
    decodedata(:,headera-2) > input.general_stain_threshold;
decodedata = decodedata(toinclude,:);
if isempty(decodedata)
    error('No read above quality and general stain threshold.')
end
textdata = textdata(toinclude,:);

% expected barcode
[uTags, idxtag1, iTag] = unique(textdata(:,1));
cTag = hist(iTag, 1:length(uTags));
[uNames, ~, iName] = unique(textdata(:,2));
TagName = uNames(iName(idxtag1));
iNNNN = strcmp(uNames, 'NNNN');
try
    unexpected = iName(idxtag1)==find(iNNNN);
catch
    unexpected = false(numel(uTags), 1);    % no NNNN
end
iHomo = findhomomer_read(textdata(:,1), uTags(~unexpected));

% figures
countgene = histread(textdata(:,2), iHomo, ['QT-' num2str(input.quality_threshold)]);
[~, idx] = sort(cTag, 'descend');
idx = [idx(~ismember(idx, find(unexpected))), idx(ismember(idx, find(unexpected)))];
make_table_barplot(strcat(uTags(idx), ': ', TagName(idx)),...
    cTag(idx), ['Read count QT-' num2str(input.quality_threshold)]);
drawnow;

disp('writing files..');
% output
outprefix = ['QT_' num2str(input.quality_threshold) '_' num2str(input.general_stain_threshold)];

% code_n_count file
fid = fopen(fullfile(decodedir, [outprefix '_code_n_count.csv']), 'w');
fprintf(fid,'Code,Count,GeneName\n');
towrite = [uTags, num2cell(cTag'), TagName]';
fprintf(fid, '%s,%d,%s\n', towrite{:});
fclose(fid);

% gene_n_count
fid = fopen(fullfile(decodedir, [outprefix '_gene_n_count.csv']), 'w');
fprintf(fid, 'GeneName,Count\n');
countgene = countgene';
fprintf(fid, '%s,%d\n', countgene{:});
fclose(fid);

% details
fid = fopen(fullfile(decodedir, [outprefix '_details.csv']), 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
towrite = [textdata, num2cell(decodedata)]';
fprintf(fid, ['%s,%s,', lineformat('%d', length(header)-2)], towrite{:});
fclose(fid);

% details, no unexpected reads
fid = fopen(fullfile(decodedir, ['QT_' num2str(input.quality_threshold) '_details_noNNNN.csv']), 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
towrite = [textdata, num2cell(decodedata)];
towrite = (towrite(~strcmp(towrite(:,2), 'NNNN'),:))';
fprintf(fid, ['%s,%s,', lineformat('%d', length(header)-2)], towrite{:});
fclose(fid);

end
