function  Decoding_Sequencing_2(input_file,num_hybs,taglist,grouped_YN,...
    calculate_global_position_YN,csv_file_contain_tile_position,...
    output_directory_decode,output_filename_decode_prefix,...
    check_parent_cell_YN,check_alignment_YN,alignment_min_threshold,...
    abnormal_sequencing_YN,sequencing_order,...
    General_Alignment_ACGT_column_number,XYPosition_ParentCell_column_number,varargin)

% Decode the CellProfile (CP) result to reads
%
% Taglist settings:
%   Automatically converts barcode letters to digits.
%   Optional: Group the tags and give group tag.
%       Eg. 'ATCG STK15 mutation','go';
%       The same group tag can be shared by several barcodes.
%
% Calculate global position
%    Specify the CSV file used as CellProfiler input. It contains the x-
%    and y- starting position of tiles.
%
% Ouput directoy: a sub-folder will be created within the current working
% directory.
%
% Option1: Check parent cell.
%   All blobs without a parent cell will be discarded.
%   Set to 1 if the blob-cell relation is fairly good in CellProfiler.
%
% Option2: Check alignment.
%   The recommended alignemnt score threshold is 1.8-2
%   Only blobs with good alignemnt in all hybridization steps will be
%   saved. Others will be discarded.
%   Not always necessary, since the general stain for each hybridization
%   step will be checked later in thresholding.
%
% Option3: Remove T in analysis.
%   The T channel will be directly excluded in decoding.
%   But an additional read guess regarding T will be done and give a ui
%   table.
%   For now, these reads are still considered as unexpected.
%
% Option4: Abnormal sequencing.
%   If you use some special sequencing order,
%   or if something has gone wrong in CellProfiler analysis...
%   The length of sequencing order has to match with number of hybs (except
%   when there is 5th empty cycle). The number at each position should be
%   bettween 0 to 4 (base position in barcode).
%   Example1: reverse sequencing.
%       sequencing_order = '4321';
%   Example2: sequenced base3 twice and you want reads from both
%       num_hybs = 5; sequencing_order = '12334';
%   Example3: sequenced base3 twice and you want read only from the second
%       num_hybs = 5; sequencing_order = '12034';
%   Example4: base2 has very poor alignment and you want exclude it
%       num_hybs = 4; sequencing_order = '1034';
%
%
% Decoding_Sequencing v2.5.1, lite
% Xiaoyan, 2015-12-7

%% initiate
drawnow;
disp('Initiating Decoding_Sequencing.');
total_t1 = clock;

%% import and structure data, reform taglist
checkinputandimport;

[taglist,expected_list,expected_digits,~,exp_tags,exp_groups] =  ...
    maketag(taglist,num_hybs,0,0,grouped_YN,abnormal_sequencing_YN,sequencing_order);

%% preallocate arrays
seq_res = zeros(num_blobs,num_hybs);
channel_strength_max = zeros(num_blobs,num_hybs);
channel_strength_sum = zeros(num_blobs,num_hybs);
num_channels = 4;
channel_strength_original = zeros(num_hybs,num_channels,num_blobs);
seq_quality = zeros(num_blobs,num_hybs);
general_strength = zeros(num_blobs,num_hybs);
alignment_score = zeros(num_blobs,num_hybs);
tile_ID = zeros(num_blobs,1);
tile_cell_ID = zeros(num_blobs,1);
cell_ID = zeros(num_blobs,1);
tile_blob_ID = zeros(num_blobs,1);
blob_ID = zeros(num_blobs,1);
tile_x_pos = zeros(num_blobs,1);
tile_y_pos = zeros(num_blobs,1);

if calculate_global_position_YN
    global_x_pos = zeros(num_blobs,1);
    global_y_pos = zeros(num_blobs,1);
end


%% main function: extract data
start_blob_ID=0;
start_cell_ID=0;
fprintf('number of tiles\t number of blobs\n');
fprintf('\t%6d\t\t%6d\n',num_tiles,num_blobs);
fprintf('# tile\tnumber of blobs\n');

for t=1:num_tiles
    temp_tile_seq_data = seqdata(seqdata(:,3) == t,:);
    temp_tile_num_blobs = length(temp_tile_seq_data(:,2))/num_hybs;
    temp_tile_num_cells = max(temp_tile_seq_data(:,12));
    
    %------------------------
    if calculate_global_position_YN
        temp_tile_x_start = mean(tile_start_pos(tile_start_pos(:,1) == t,2));
        temp_tile_y_start = mean(tile_start_pos(tile_start_pos(:,1) == t,3));
    end
    fprintf('%6d\t\t%6d\n',t,temp_tile_num_blobs);
    
    %------------------------
    if temp_tile_num_blobs
        [~,order] = sort(temp_tile_seq_data(:,2));
        
        % extract info about blobs within the tile
        temp_tile_intensity = temp_tile_seq_data(order,6:6+num_channels-1);
        temp_tile_general = temp_tile_seq_data(order,4);
        temp_tile_align = temp_tile_seq_data(order,5);
        temp_tile_x_pos = temp_tile_seq_data(order,10);
        temp_tile_y_pos = temp_tile_seq_data(order,11);
        temp_tile_cell_id = temp_tile_seq_data(order,12);
        
        % max intensities
        [maxint,maxidx] = max(temp_tile_intensity,[],2);
        channel_strength_max(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxint,num_hybs,temp_tile_num_blobs))';
        
        % seq result matrix
        seq_res(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(maxidx,num_hybs,temp_tile_num_blobs))';
        
        % blob info
        tile_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = t;
        tile_blob_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            (1:temp_tile_num_blobs)';
        blob_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            (1:temp_tile_num_blobs)'+repmat(start_blob_ID,temp_tile_num_blobs,1);
        alignment_score(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(temp_tile_align,num_hybs,temp_tile_num_blobs))';
        
        
        temp_list = 1:temp_tile_num_blobs*num_hybs;
        if num_hybs>1
            temp_list = mod(temp_list,num_hybs)==1;
        end
        
        tile_cell_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_cell_id(temp_list);
        cell_ID(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_cell_id(temp_list)+(temp_tile_cell_id(temp_list)~=0)*start_cell_ID;
        tile_x_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_x_pos(temp_list);
        tile_y_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
            temp_tile_y_pos(temp_list);
        % global coordinates
        if calculate_global_position_YN
            global_x_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                repmat(temp_tile_x_start,temp_tile_num_blobs,1) + temp_tile_x_pos(temp_list);
            global_y_pos(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                repmat(temp_tile_y_start,temp_tile_num_blobs,1) + temp_tile_y_pos(temp_list);
        end
        
        % original intensities
        for c = 1:num_channels
            channel_strength_original(:,c,start_blob_ID+1:start_blob_ID+temp_tile_num_blobs) = ...
                reshape(temp_tile_intensity(:,c),num_hybs,1,temp_tile_num_blobs);
        end
        general_strength(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            (reshape(temp_tile_general,num_hybs,temp_tile_num_blobs))';
        channel_strength_sum(start_blob_ID+1:start_blob_ID+temp_tile_num_blobs,:) = ...
            reshape(sum(temp_tile_intensity,2),num_hybs,temp_tile_num_blobs)';
        
        start_blob_ID = start_blob_ID+temp_tile_num_blobs;
        start_cell_ID = start_cell_ID+temp_tile_num_cells;
    end
    
end
if ~calculate_global_position_YN
    global_x_pos = tile_x_pos;
    global_y_pos = tile_y_pos;
end
clear -regexp ^temp
clear t c

%% bases into reads
[seq_num,inuni] = combinebase(abnormal_sequencing_YN,num_hybs,0,...
    num_blobs,seq_res,sequencing_order);

%% sequencing quality
seq_quality = channel_strength_max./channel_strength_sum;
qualityitems;
qualitycheck;

%% all reads before QT
allbt = seq_num(ffbt);

if abnormal_sequencing_YN
    seq_res_allbt = seq_res(ffbt,inuni);
    seq_quality_allbt = seq_quality(ffbt,inuni);
    general_strength_allbt = general_strength(ffbt,inuni);
    alignment_score_allbt = alignment_score(ffbt,inuni);
    channel_strength_max_allbt = channel_strength_max(ffbt,inuni);
    channel_strength_original_allbt = channel_strength_original(inuni,:,ffbt);
else
    seq_res_allbt = seq_res(ffbt,:);
    seq_quality_allbt = seq_quality(ffbt,:);
    general_strength_allbt = general_strength(ffbt,:);
    alignment_score_allbt = alignment_score(ffbt,:);
    channel_strength_max_allbt = channel_strength_max(ffbt,:);
    channel_strength_original_allbt = channel_strength_original(:,:,ffbt);
    
end
letters_allbt = num2letters_2(seq_res_allbt);

tile_allbt = tile_ID(ffbt);
tile_cell_allbt = tile_cell_ID(ffbt);
cell_allbt = cell_ID(ffbt);
tile_blob_allbt = tile_blob_ID(ffbt);
blob_allbt = blob_ID(ffbt);
tile_x_allbt = tile_x_pos(ffbt);
tile_y_allbt = tile_y_pos(ffbt);
global_x_allbt = global_x_pos(ffbt);
global_y_allbt = global_y_pos(ffbt);

seq_quality_min_allbt = seq_quality_min(ffbt,:);
general_strength_min_allbt = general_strength_min(ffbt,:);
alignment_score_min_allbt = alignment_score_min(ffbt,:);
seq_quality_min_ind_allbt = seq_quality_min_ind(ffbt,:);
general_strength_min_ind_allbt = general_strength_min_ind(ffbt,:);
alignment_score_min_ind_allbt = alignment_score_min_ind(ffbt,:);

%% assign gene name and group tag to each read
[name_tag_allbt,group_tag_allbt] = ...
    genename(grouped_YN,allbt,expected_list,exp_tags,exp_groups);

%% group reads
[uni_seq_res_allbt,count_code_allbt,uni_letters_allbt,uni_reads_name_tag_allbt,uni_reads_group_tag_allbt,catagory_read_allbt] = ...
    groupreads(allbt,seq_res_allbt,letters_allbt,expected_list,exp_tags,exp_groups);

% find unexpected homomer reads
homomer_index_allbt = findHomomer(expected_list,allbt);
catagory_read_allbt(homomer_index_allbt) = -1;

%% output figures and files
disp('busy drawing figures..');
tic
if grouped_YN
    count_gene= histRead_group(allbt,exp_tags,exp_groups,...
        expected_list,'before QT',homomer_index_allbt);
    count_group = histRead(allbt,exp_groups,...
        expected_list,'before QT',homomer_index_allbt,1);
else
    count_gene = histRead(allbt,exp_tags,...
        expected_list,'before QT',homomer_index_allbt,1);
end
[thres,count] = barQuality_2(0,1,catagory_read_allbt,seq_quality_min_allbt,1);
drawnow;
toc

writefiles;

%% functions
% function to check input arguments, and import data
    function checkinputandimport
        if exist(input_file,'file')
        else error('Could not find the input file.');
        end
        
        if calculate_global_position_YN
            if exist(csv_file_contain_tile_position,'file')
            else error('Could not find the position file.');
            end
            
        end
        
        % import data from files
        %--------------------
        seqdata = decodinginput(input_file,General_Alignment_ACGT_column_number,...
            XYPosition_ParentCell_column_number);
        
        if calculate_global_position_YN
            tile_start_pos = getposition(csv_file_contain_tile_position);
        end
        
        num_blobs = size(seqdata,1)/num_hybs;
        tile_ID_all = seqdata(:,3);
        num_tiles = max(tile_ID_all);
        if floor(num_blobs)==num_blobs
        else
            error('Oops! The number of blobs is not integer. Check num_hybs again.');
        end
    end

% function to get quality items
    function qualityitems
        seq_quality(isnan(seq_quality)) = 0;
        
        if abnormal_sequencing_YN
            [seq_quality_min,seq_quality_min_ind] = min(seq_quality(:,inuni),[],2);
            [general_strength_min,general_strength_min_ind] = min(general_strength(:,inuni),[],2);
            [alignment_score_min,alignment_score_min_ind] = min(alignment_score(:,inuni),[],2);
        else
            [seq_quality_min,seq_quality_min_ind] = min(seq_quality,[],2);
            [general_strength_min,general_strength_min_ind] = min(general_strength,[],2);
            [alignment_score_min,alignment_score_min_ind] = min(alignment_score,[],2);
        end
    end

% function to remove nonsensical reads
    function qualitycheck
        if check_parent_cell_YN && check_alignment_YN
            ffbt = (seq_quality_min>0 & general_strength_min>0 & cell_ID>0 & alignment_score_min>=alignment_min_threshold);
        elseif check_parent_cell_YN && check_alignment_YN==0
            ffbt = (seq_quality_min>0 & general_strength_min>0 & cell_ID>0);
        elseif check_parent_cell_YN==0 && check_alignment_YN
            ffbt = (seq_quality_min>0 & general_strength_min>0 & alignment_score_min>=alignment_min_threshold);
        else
            ffbt = (seq_quality_min>0 & general_strength_min>0);
        end
        
        if nnz(ffbt)
        else
            error('No reads available. Try without parent cell check or alignment score.');
        end
    end

% function to write output files
    function writefiles
        disp('writing files..');
        tic
        output_directory_decode = [output_directory_decode '\'];
        if exist(output_directory_decode, 'dir')
        else mkdir (output_directory_decode);
        end
        
        % write code_n_count file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_code_n_count' '.csv'], 'w');
        if grouped_YN
            fprintf(fid,'Code,Count,GeneName,Group\n');
            temp_write = [cellstr(uni_letters_allbt),...
                num2cell(count_code_allbt'),...
                uni_reads_name_tag_allbt,...
                uni_reads_group_tag_allbt];
            for row = 1:size(temp_write,1)
                fprintf(fid,'%s,%d,%s,%s\n',temp_write{row,:});
            end
            
        else
            fprintf(fid,'Code,Count,GeneName\n');
            temp_write = [cellstr(uni_letters_allbt),...
                num2cell(count_code_allbt'),...
                uni_reads_name_tag_allbt];
            for row = 1:size(temp_write,1)
                fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
            end
        end
        fclose(fid);
        
        % write gene_n_count file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_gene_n_count' '.csv'], 'w');
        if grouped_YN
            fprintf(fid,'GeneName,Count,Group\n');
            for row=1:length(count_gene(:,1))
                fprintf(fid,'%s,%d,%s\n',count_gene{row,:});
            end
        else
            fprintf(fid,'GeneName,Count\n');
            for row=1:length(count_gene(:,1))
                fprintf(fid,'%s,%d\n',count_gene{row,:});
            end
        end
        fclose(fid);
        
        % write group_n_count file
        %--------------------
        if grouped_YN
            fid = fopen([output_directory_decode output_filename_decode_prefix...
                '_group_n_count' '.csv'], 'w');
            fprintf(fid,'GeneGroup,Count\n');
            for row=1:length(count_group(:,1))
                fprintf(fid,'%s,%d\n',count_group{row,:});
            end
            fclose(fid);
        end
        
        % write qualitybar file
        %--------------------
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_qualitybar' '.csv'], 'w');
        fprintf(fid,'threshold,expected,unexpected,homomer,belowQT\n');
        temp_write = [thres',count];
        fprintf(fid,'%.2f,%d,%d,%d,%d\n',temp_write');
        fclose(fid);
        
        % write details file
        %--------------------
        tempheader = {'letters','name','group',...
            'tile_X_pos','tile_Y_pos','global_X_pos','global_Y_pos',...
            'parent_cell','tile_ID',...
            'general_stain_min','seq_quality_min','alignment_score_min'};
        tempwrite = [cellstr(letters_allbt),name_tag_allbt,group_tag_allbt,...
            num2cell([tile_x_allbt,tile_y_allbt,global_x_allbt,global_y_allbt,...
            cell_allbt,tile_allbt,...
            general_strength_min_allbt,seq_quality_min_allbt,alignment_score_min_allbt])];
        
        tempidx = [];
        if calculate_global_position_YN
            tempidx = [tempidx,4,5];
        else
            tempidx = [tempidx,6,7];
        end
        if ~check_alignment_YN
            tempidx = [tempidx,12];
        end
        if ~check_parent_cell_YN
            tempidx = [tempidx,8];
        end
        if grouped_YN
            tempidx = [tempidx,1];
        else
            tempidx = [tempidx,3];
        end
        
        tempheader(tempidx) = [];
        tempwrite(:,tempidx) = [];
        
        
        fid = fopen([output_directory_decode output_filename_decode_prefix...
            '_details' '.csv'], 'w');
        % header
        tempformat = [];
        for i = 1:length(tempheader)
            tempformat = [tempformat,'%s,'];
        end
        tempformat(end) = [];
        tempformat = [tempformat,'\n'];
        fprintf(fid,tempformat,tempheader{:});
        
        % main table
        tempformat = '%s,%s';
        for i = 1:length(tempheader)-2
            tempformat = [tempformat,',%d'];
        end
        tempformat = [tempformat,'\n'];
        for i = 1:size(tempwrite,1)
            fprintf(fid,tempformat,tempwrite{i,:});
        end
        
        fclose(fid);
    end

toc
clear -regexp ^temp


%% save mat file
total_t2 = clock;
total_t = etime(total_t2,total_t1);
fprintf('%s%6f%s\n','Total elapsed time: ',total_t,'s');
clear -regexp ^temp ^total
disp('saving workspace variables..');
save([output_directory_decode output_filename_decode_prefix '.mat']);

fprintf('Decoding beforeQT finished.\n\n');
end