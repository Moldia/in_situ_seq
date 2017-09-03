function seqtiling(input)
% seqtiling(input)
% re-tile images for in situ sequencing image analysis
% input from Sequencing script
% Xiaoyan, 2017


nChannels = input.channel_end - input.channel_start + 1;
nBases = input.base_end - input.base_start + 1;

checkinput;
checkimages;

% pad and tile images
nImages = 0;
imdirs = cell(ntilesX*ntilesY, input.base_end, input.channel_end);
for b = input.base_start:input.base_end
    for c = input.channel_start:input.channel_end
        if input.in_subfolder_YN
            imfile = fullfile(input.folder_image,...
                [input.filename_base_prefix num2str(b)],...
                [input.filename_base_prefix num2str(b) ...
                input.filename_channel_prefix num2str(c) input.filename_suffix]);
        else
            imfile = fullfile(input.folder_image,...
                [input.filename_base_prefix  num2str(b) ...
                input.filename_channel_prefix num2str(c) input.filename_suffix]);
        end
        
        I = imread(imfile);
        nImages = nImages+1;
        I = padimg(I, padX(b,c), padY(b,c));
        suffix = strsplit(input.filename_suffix, '.');
        suffix = suffix{1};
        outdir = fullfile(input.folder_image,...
            ['Tiled_' input.filename_base_prefix num2str(b)...
            input.filename_channel_prefix num2str(c) suffix]);
        
        disp(['tiling image ' num2str(nImages) '..']);
        tileimage(I, input.tile_size, outdir);
        [imdir, imfiles] = imfiles_in_grid(outdir, 'tile', '.tif', 0, [ntilesX, ntilesY]);
        imdirs(:,b,c) = reshape(imdir', [], 1);
    end
end
imfiles = reshape(imfiles', [], 1);

tileposX = repmat(0:input.tile_size:input.tile_size*(ntilesX-1), ntilesY, 1);
tileposY = repmat((0:input.tile_size:input.tile_size*(ntilesY-1))', 1, ntilesX);

% prepare to write position csv file
towrite = [];
header = input.channel_order;
metadata = num2cell([(1:ntilesX*ntilesY)', reshape(tileposX',[],1), reshape(tileposY',[],1)]);
for b = input.base_start:input.base_end
    perbase = [];
    for c = input.channel_start:input.channel_end
        perbase = [perbase, imdirs(:,b,c), imfiles(:)];
    end
    towrite = [towrite; [metadata, num2cell(repmat(b, ntilesX*ntilesY, 1)), perbase]];
end

% general stain
cGeneral = find(ismember(input.channel_order,...
    {'General_stain', 'general_stain', 'General_blob', 'general_blob',...
    'General stain', 'general stain', 'General blob', 'general blob',...
    'AF750', 'Cy7'}));
if length(cGeneral) == 1
    towrite = [towrite,...
        repmat([imdirs(:,input.base_start,cGeneral), imfiles], nBases, 1)];
    header{cGeneral} = 'Spec_blob';
    header = [header, 'General_blob'];
end

% DAPI (from base1 only)
cNuclei = find(ismember(input.channel_order,...
    {'Nuclei','nuclei','DAPI','Hoechst'}));
if length(cNuclei) == 1
    towrite(:,4+(cNuclei-1)*2+1,:) =...
        repmat(towrite(1:ntilesX*ntilesY,4+(cNuclei-1)*2+1), nBases, 1);
end
[~, idx] = sort(cell2mat(towrite(:,1)));
towrite = towrite(idx,:)';

% organize file header
header = [strcat('Image_PathName_', header); strcat('Image_FileName_', header)];
header = [{'Metadata_position','Tile_xPos','Tile_yPos','Hyb_step'}, reshape(header, 1, [])];

% write file
fid = fopen(fullfile(input.folder_image, [input.CSV_filename_prefix '.csv']), 'w');
fprintf(fid, lineformat('%s', length(header)), header{:});
fprintf(fid, ['%d,%d,%d,%d,' lineformat('%s', length(header)-4)], towrite{:});
fclose(fid);

        
    function checkinput
        if nChannels<0 || nBases<0 || nChannels~=uint8(nChannels) || nBases~=uint8(nBases)
            error('Number of channels and bases must be positive integers.');
        end
        
        if length(input.channel_order) ~= nChannels
            error(['Channel number and the number of elements '...
                'in channel_order do not match.']);
        end
    end

    function checkimages
        imdimensions = zeros(input.base_end, input.channel_end, 2);
        disp('checking images..');
        for b = input.base_start:input.base_end
            for c = input.channel_start:input.channel_end
                if input.in_subfolder_YN
                    imfile = fullfile(input.folder_image,...
                        [input.filename_base_prefix num2str(b)],...
                        [input.filename_base_prefix num2str(b)...
                        input.filename_channel_prefix num2str(c)...
                        input.filename_suffix]);
                else
                    imfile = fullfile(input.folder_image,...
                        [input.filename_base_prefix num2str(b)...
                        input.filename_channel_prefix num2str(c)...
                        input.filename_suffix]);
                end
                
                f = imfinfo(imfile);
                imdimensions(b,c,:) = [f.Width, f.Height];
            end
        end
        maxX = max(reshape(imdimensions(:,:,1), [], 1));
        maxY = max(reshape(imdimensions(:,:,2), [], 1));
        ntilesX = ceil(maxX/input.tile_size);
        ntilesY = ceil(maxY/input.tile_size);
        
        padX = input.tile_size*ntilesX - imdimensions(:,:,1);
        padY = input.tile_size*ntilesY - imdimensions(:,:,2);
    end

disp('Tiling finished.');
end

