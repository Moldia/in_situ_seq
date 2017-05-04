function decoding(input)
% decoding for in situ sequencing analysis
% input from Sequencing script
% Xiaoyan, 2017


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
    for i = 1:length(seqorder)
        numseqorder(i) = str2double(seqorder(i));
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
idBlob = zeros(nBlobs, 2);
basecalling = zeros(nBlobs, nHybs);
maxChannel = zeros(nBlobs, nHybs);
sumChannel = zeros(nBlobs, nHybs);
oriChannel = zeros(nHybs, 4, nBlobs);
anchor = zeros(nBlobs, nHybs);
alignment = zeros(nBlobs, nHybs);
idCell = zeros(nBlobs, 1);
globalpos = zeros(nBlobs, 2);

% base-calling
startblob = 0;
startcell = 0;
fprintf('number of tiles\t number of blobs\n');
fprintf('\t%6d\t\t%6d\n', max(seqdata(:,3)), nBlobs);
fprintf('# tile\tnumber of blobs\n');

for t = 1:max(seqdata(:,3))
    tiledata = seqdata(seqdata(:,3) == t,:);
    nTileblobs = size(tiledata,1)/nHybs;
    nCellsTile = max(tiledata(:,12));
    
    fprintf('%6d\t\t%6d\n', t, nTileblobs);
    
    if nTileblobs
        [~, idxhyb1, idxblob] = unique(tiledata(:,2));
        [~, order] = sort(idxblob);
        
        % blob index and tile ID
        idBlob(startblob+1:startblob+nTileblobs, 1) =...
            (1:nTileblobs)' + startblob;
        idBlob(startblob+1:startblob+nTileblobs, 2) = t;
        
        % max intensity and base-calling
        [maxint, maxidx] = max(tiledata(order,6:9), [], 2);
        maxChannel(startblob+1:startblob+nTileblobs,:) = ...
            (reshape(maxint, nHybs, nTileblobs))';
        basecalling(startblob+1:startblob+nTileblobs,:) = ...
            (reshape(maxidx, nHybs, nTileblobs))';
        
        % original intensities
        sumChannel(startblob+1:startblob+nTileblobs,:) = ...
            reshape(sum(tiledata(order,6:9), 2), nHybs, nTileblobs)';
        temp = reshape(tiledata(order,6:9)', 4, nHybs, []);
        temp = permute(temp, [2, 1, 3]);
        oriChannel(:,:,startblob+1:startblob+nTileblobs) = temp;
        
        % alignment and general stain
        alignment(startblob+1:startblob+nTileblobs,:) = ...
            (reshape(tiledata(order,5), nHybs, nTileblobs))';
        anchor(startblob+1:startblob+nTileblobs,:) = ...
            (reshape(tiledata(order,4),nHybs,nTileblobs))';
        
        % global coordinates
        try
            globalpos(startblob+1:startblob+nTileblobs,:) =...
                bsxfun(@plus,...
                tilepos(tilepos(:,1)==t,2:3),...
                tiledata(idxhyb1,10:11));
        catch
            globalpos(startblob+1:startblob+nTileblobs,:) =...
                tiledata(idxhyb1,10:11);
        end
        
        % parent cell
        idCell(startblob+1:startblob+nTileblobs) = ...
            tiledata(idxhyb1,12) + tiledata(idxhyb1,12)*startcell;        
        
        startblob = startblob+nTileblobs;
        startcell = startcell+nCellsTile;
    end
    
end
clear temp

% remove '0' hybs as specified in sequencing order
basecalling = basecalling(:,readorder);
propsPerBase = cat(3, maxChannel, sumChannel, anchor, alignment);
propsPerBase = propsPerBase(:,readorder,:);

% bases into reads
seqnum = zeros(size(basecalling,1),1);
for i = 1:length(readorder)
    seqnum = seqnum + basecalling(:,i)*10^(length(readorder)-i);
end

%% sequencing quality
quality = maxChannel./sumChannel;
qualitycheck;

% all reads before QT
seqnum = seqnum(toinclude);
propsReads = propsReads(toinclude,:);
oriChannel = oriChannel(:,:,toinclude);
hybMinProps = hybMinProps(toinclude,:);
propsPerBase = propsPerBase(toinclude,:,:);

% gene name
[uNums, ~, idxNum] = unique(seqnum);
cNums = hist(idxNum, length(uNums));
uReads = cellfun(@(v) find(v==expectedNum), num2cell(uNums), 'uni', 0);
unexpected = cellfun(@isempty, uReads);
uNames = cell(length(uNums),1);
uNames(unexpected) = {'NNNN'};
uNames(~unexpected) = genes([uReads{~unexpected}]);
uReads = num2barcode(uNums);

% find unexpected homomer reads
idxHomo = findhomomer_num(seqnum, expectedNum);

tic
Name = cell(length(seqnum),1);
for i = 1:length(uNums)
    Name(idxNum==i) = uNames(i);
end
toc

%% output figures and files
disp('drawing figures..');
tic
countgene = histread(Name, idxHomo, 'beforeQT');

% bar plot
[~, idx] = sort(cNums, 'descend');
idx = [idx(~ismember(idx,find(unexpected))), idx(ismember(idx,find(unexpected)))];
make_table_barplot(strcat(uReads(idx), ': ', uNames(idx)),...
    cNums(idx), 'Read count before QT');

% quality bar
category = ones(length(seqnum), 1);
category(ismember(idxNum, find(unexpected))) = 0;
category(idxHomo) = -1;
[thres, count] = qualitybar(0, 1, category, propsReads(:,6));
drawnow;
toc

disp('writing files..');
outdir = input.output_directory;
if ~strcmp(outdir(end), filesep)
    outdir = [outdir filesep];
end
mkdir(outdir);

% code_n_count file
fid = fopen([outdir 'beforeQT_code_n_count.csv'], 'w');
fprintf(fid,'Code,Count,GeneName\n');
towrite = [uReads, num2cell(cNums'), uNames]';
fprintf(fid, '%s,%d,%s\n', towrite{:});
fclose(fid);

% gene_n_count
fid = fopen([outdir 'beforeQT_gene_n_count.csv'], 'w');
fprintf(fid, 'GeneName,Count\n');
countgene = countgene';
fprintf(fid, '%s,%d\n', countgene{:});
fclose(fid);

% quality bar
fid = fopen([outdir 'beforeQT_qualitybar.csv'], 'w');
fprintf(fid, 'threshold,expected,unexpected,homomer,belowQT\n');
towrite = [thres',count];
fprintf(fid,'%.2f,%d,%d,%d,%d\n',towrite');
fclose(fid);

% details
fid = fopen([outdir 'beforeQT_details.csv'], 'w');
fprintf(fid, 'Read,Gene,PosX,PosY,ParentCell,Tile,MinAnchor,MinQuality,MinAlign\n');
towrite = [uReads(idxNum), Name, num2cell(propsReads)]';
fprintf(fid, ['%s,%s,', lineformat('%d', size(propsReads,2))], towrite{:});
fclose(fid);

% mat file
try
    if input.savemat
        save([outdir 'beforeQT.mat'], 'toinclude', seqnum, 'oriChannel',...
            'propsReads', 'hybMinProps', 'propsPerBase');
    end
end
        
fprintf('Decoding beforeQT finished.\n\n');

%% functions
   function qualitycheck
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
