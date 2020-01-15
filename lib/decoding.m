function decoding(input)
% decoding(input)
% decoding for in situ sequencing analysis
% input from InSituSequencing script
% Xiaoyan, 2018


drawnow
% taglist, expected barcode and tile position
% expected barcode taglist
if ischar(input.taglist)
    if isempty(input.taglist)
        taglist = cell(0, 2);
    else
        taglist = importdata(input.taglist);
        taglist = cellfun(@(v) strsplit(v, ','), taglist, 'uni', 0);
        taglist = reshape([taglist{:}], length(taglist{1}), [])';
    end
else
    taglist = formattaglist(input.taglist);
end
taglist = taglist(:,1:2);

% abnormal sequencing order
nHybs = input.num_hybs;
if input.abnormal_sequencing_YN
    seqorder = input.sequencing_order;
    seqorder(strfind(seqorder, '0')) = [];
    numseqorder = [];
    for iCell = 1:length(seqorder)
        numseqorder(iCell) = str2double(seqorder(iCell));
    end
    taglist(:,1) = cellfun(@(v) v(numseqorder), taglist(:,1), 'uni', 0);
    readorder = 1:nHybs;
    readorder(strfind(input.sequencing_order, '0')) = [];
else
    readorder = 1:nHybs;
end

taglist = remove_barcode_duplicates(taglist);
[expectedNum, ~] = barcode2num(taglist(:,1));
genes = taglist(:,2);

% tile position
if ~isempty(input.csv_file_contain_tile_position)
    tilepos = getcsvtilepos(input.csv_file_contain_tile_position);
end

%% load and procees intensity data 
disp('loading Cellprofiler results..')
seqdata = decodinginput(input.input_file,...
    input.General_Alignment_ACGT_column_number,...
    input.XYPosition_ParentCell_column_number);
nBlobs = size(seqdata, 1)/nHybs;
if floor(nBlobs) ~= nBlobs
    error('Oops! The number of blobs is not integer. Check num_hybs again.')
end

% preallocate
iBlob = zeros(nBlobs, 2);
basecall = zeros(nBlobs, nHybs);
maxChannel = zeros(nBlobs, nHybs);
sumChannel = zeros(nBlobs, nHybs);
originalChannel = zeros(nHybs, 4, nBlobs);
anchor = zeros(nBlobs, nHybs);
alignment = zeros(nBlobs, nHybs);
iCell = zeros(nBlobs, 1);
globalpos = zeros(nBlobs, 2);

% base-calling
startblob = 0;
startcell = 0;
fprintf('number of tiles\t number of blobs\n');
fprintf('\t%6d\t\t%6d\n', max(seqdata(:,3)), nBlobs);
fprintf('# tile\tnumber of blobs\n');

for t = 1:max(seqdata(:,3))
    tiledata = seqdata(seqdata(:,3) == t,:);
    nBlobsInTile = size(tiledata,1)/nHybs;
    nCellsInTile = max(tiledata(:,12));
    
    fprintf('%6d\t\t%6d\n', t, nBlobsInTile);
    
    if nBlobsInTile
        [~, idxHyb1, idxBlob] = unique(tiledata(:,2));
        [~, order] = sort(idxBlob);
        
        % blob index 
        iBlob(startblob+1:startblob+nBlobsInTile, 1) =...
            (1:nBlobsInTile)' + startblob;
        
        % tile ID
        iBlob(startblob+1:startblob+nBlobsInTile, 2) = t;
        
        % max intensity and base-calling
        [maxint, maxidx] = max(tiledata(order,6:9), [], 2);
        maxChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxint, nHybs, nBlobsInTile))';
        basecall(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(maxidx, nHybs, nBlobsInTile))';
        
        % original intensities
        sumChannel(startblob+1:startblob+nBlobsInTile,:) = ...
            reshape(sum(tiledata(order,6:9), 2), nHybs, nBlobsInTile)';
        temp = reshape(tiledata(order,6:9)', 4, nHybs, []);
        temp = permute(temp, [2, 1, 3]);
        originalChannel(:,:,startblob+1:startblob+nBlobsInTile) = temp;
        
        % alignment and general stain
        alignment(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,5), nHybs, nBlobsInTile))';
        anchor(startblob+1:startblob+nBlobsInTile,:) = ...
            (reshape(tiledata(order,4),nHybs,nBlobsInTile))';
        
        % global coordinates
        try
            globalpos(startblob+1:startblob+nBlobsInTile,:) =...
                bsxfun(@plus,...
                tilepos(tilepos(:,1)==t,2:3),...
                tiledata(idxHyb1,10:11));
        catch
            globalpos(startblob+1:startblob+nBlobsInTile,:) =...
                tiledata(idxHyb1,10:11);
        end
        
        % parent cell
        iCell(startblob+1:startblob+nBlobsInTile) = ...
            tiledata(idxHyb1,12) + tiledata(idxHyb1,12)*startcell;        
        
        startblob = startblob + nBlobsInTile;
        startcell = startcell + nCellsInTile;
    end
    
end
clear temp

% remove '0' hybs as specified in sequencing order
basecall = basecall(:,readorder);
propsPerBase = cat(3, maxChannel, sumChannel, anchor, alignment);
propsPerBase = propsPerBase(:,readorder,:);

% bases into reads
readsNum = zeros(size(basecall,1),1);
for i = 1:length(readorder)
    readsNum = readsNum + basecall(:,i)*10^(length(readorder)-i);
end

%% sequencing quality
quality = propsPerBase(:,:,1)./propsPerBase(:,:,2);
[toinclude, propsReads, hybMinProps] = qualitycheck...
    (input, quality, propsPerBase(:,:,3), propsPerBase(:,:,4), iCell, iBlob, globalpos);

% mat file
mkdir(input.output_directory);
try
    if input.savemat
        save(fullfile(input.output_directory, 'beforeQT.mat'), 'toinclude', 'readsNum', 'originalChannel',...
            'propsReads', 'hybMinProps', 'propsPerBase');
    end
end

% all reads before QT
readsNum = readsNum(toinclude);
propsReads = propsReads(toinclude,:);
originalChannel = originalChannel(:,:,toinclude);
hybMinProps = hybMinProps(toinclude,:);
propsPerBase = propsPerBase(toinclude,:,:);

% reads and gene name
[uReads, ~, iRead] = unique(readsNum);
cRead = hist(iRead, length(uReads));
iExpected = cellfun(@(v) find(v==expectedNum), num2cell(uReads), 'uni', 0);
unexpected = cellfun(@isempty, iExpected);
uNames = cell(length(uReads), 1);
uNames(unexpected) = {'NNNN'};
uNames(~unexpected) = genes([iExpected{~unexpected}]);
uTags = num2barcode(uReads);

% find unexpected homomer reads
idxHomo = findhomomer_num(readsNum, expectedNum);

disp('assigning gene names..');
Name = cell(length(readsNum), 1);
for i = 1:length(uTags)
    Name(iRead==i) = uNames(i);
end

%% output figures and files
disp('drawing figures..');
tic
countgene = histread(Name, idxHomo, 'beforeQT');

% bar plot
[~, idx] = sort(cRead, 'descend');
idx = [idx(~ismember(idx,find(unexpected))), idx(ismember(idx,find(unexpected)))];
make_table_barplot(strcat(uTags(idx), ': ', uNames(idx)),...
    cRead(idx), 'Read count before QT');

% quality bar
category = ones(length(readsNum), 1);
category(ismember(iRead, find(unexpected))) = 0;
category(idxHomo) = -1;
[thres, count] = qualitybar(0, 1, category, propsReads(:,6));
drawnow;
toc

disp('writing files..');

% code_n_count file
fid = fopen(fullfile(input.output_directory, 'beforeQT_code_n_count.csv'), 'w');
fprintf(fid,'Code,Count,GeneName\n');
towrite = [uTags, num2cell(cRead'), uNames]';
fprintf(fid, '%s,%d,%s\n', towrite{:});
fclose(fid);

% gene_n_count
fid = fopen(fullfile(input.output_directory, 'beforeQT_gene_n_count.csv'), 'w');
fprintf(fid, 'GeneName,Count\n');
countgene = countgene';
fprintf(fid, '%s,%d\n', countgene{:});
fclose(fid);

% quality bar
fid = fopen(fullfile(input.output_directory, 'beforeQT_qualitybar.csv'), 'w');
fprintf(fid, 'threshold,expected,unexpected,homomer,belowQT\n');
towrite = [thres',count];
fprintf(fid,'%.2f,%d,%d,%d,%d\n',towrite');
fclose(fid);

% details
fid = fopen(fullfile(input.output_directory, 'beforeQT_details.csv'), 'w');
fprintf(fid, 'Read,Gene,PosX,PosY,ParentCell,Tile,MinAnchor,MinQuality,MinAlign\n');
towrite = [uTags(iRead), Name, num2cell(propsReads)]';
fprintf(fid, ['%s,%s,', lineformat('%d', size(propsReads,2))], towrite{:});
fclose(fid);

        
fprintf('Decoding beforeQT finished.\n\n');

%% functions
   function [toinclude, propsReads, hybMinProps] = qualitycheck...
           (input, quality, anchor, alignment, idCell, idBlob, globalpos)
        quality(isnan(quality)) = 0;
        
        [minQ, hybMinQ] = min(quality, [], 2);
        [minAnchor, hybMinAnchor] = min(anchor, [], 2);
        [minAlign, hybMinAlign] = min(alignment, [], 2);
        propsReads = [globalpos, idCell, idBlob(:,2), minAnchor, minQ, minAlign];
        hybMinProps = [hybMinAnchor, hybMinQ, hybMinAlign];

        toinclude = minQ > 0 & minAnchor > 0;
        if input.check_parent_cell_YN
            toinclude = toinclude & idCell > 0;
        end
        if input.check_alignment_YN
            toinclude = toinclude & minAlign > input.alignment_min_threshold;
        end
        
        if nnz(toinclude)
        else
            error('No reads available. Try without parent cell check or alignment score.');
        end
   end

end
