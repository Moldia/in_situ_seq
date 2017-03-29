function Threshold_Sequencing_2(output_directory_decode,...
    output_filename_decode_prefix,output_filename_afterQT_prefix,...
    quality_threshold,general_strain_threshold,varargin)

% Take reads above threshold.
% Two thresholds: general stain inteisty and sequencing quality.
% Can remove homomer reads.
%
% Thershold_Sequencing v2.5, lite
% Xiaoyan, 2015-6-30

%% initiate
drawnow;
disp('Initiating Threshold_Sequencing.');
output_directory_decode = [output_directory_decode '\'];
if exist([output_directory_decode output_filename_decode_prefix '.mat'],'file')~=2
    error(['Could not find the beforeQT file. '...
        'Make sure the working directory is correct.']);
else
    load([output_directory_decode output_filename_decode_prefix '.mat']);
end

%% quality threshold
ff = (seq_quality_min_allbt>quality_threshold & ...
    general_strength_min_allbt>general_strain_threshold);
count_with_homomer = nnz(ff);

if ~nnz(ff)
    error('No reads available. Try a lower threshold.');
end

%% all reads after QT
allqt = allbt(ff);
seq_res_allqt = seq_res_allbt(ff,:);
letters_allqt = letters_allbt(ff,:);
name_tag_allqt = name_tag_allbt(ff);
group_tag_allqt = group_tag_allbt(ff);

seq_quality_allqt = seq_quality_allbt(ff,:);
general_strength_allqt = general_strength_allbt(ff,:);
alignment_score_allqt = alignment_score_allbt(ff,:);

channel_strength_max_allqt = channel_strength_max_allbt(ff,:);
channel_strength_original_allqt = channel_strength_original_allbt(:,:,ff);

tile_allqt = tile_allbt(ff);
tile_cell_allqt = tile_cell_allbt(ff);
cell_allqt = cell_allbt(ff);
tile_blob_allqt = tile_blob_allbt(ff);
blob_allqt = blob_allbt(ff);
tile_x_allqt = tile_x_allbt(ff);
tile_y_allqt = tile_y_allbt(ff);
global_x_allqt = global_x_allbt(ff);
global_y_allqt = global_y_allbt(ff);

seq_quality_min_allqt = seq_quality_min_allbt(ff,:);
general_strength_min_allqt = general_strength_min_allbt(ff,:);
alignment_score_min_allqt = alignment_score_min_allbt(ff,:);

if exist('seq_quality_min_ind_allbt','var')
    seq_quality_min_ind_allqt = seq_quality_min_ind_allbt(ff,:);
end
if exist('general_strength_min_ind_allbt','var')
    general_strength_min_ind_allqt = general_strength_min_ind_allbt(ff,:);
end
if exist('alignment_score_min_ind_allbt','var')
    alignment_score_min_ind_allqt = alignment_score_min_ind_allbt(ff,:);
end

%% group reads
[uni_seq_res_allqt,count_code_allqt,uni_reads_letters_allqt,uni_reads_name_tag_allqt,uni_reads_group_tag_allqt,catagory_read_allqt] = ...
    groupreads(allqt,seq_res_allqt,letters_allqt,expected_list,exp_tags,exp_groups);

% find unexpected homomer reads
homomer_index_allqt = findHomomer(expected_list,allqt);
catagory_read_allqt(homomer_index_allqt) = -1;

%% output figures
disp('busy drawing figures..');
tic
if grouped_YN
    count_gene = histRead_group(allqt,exp_tags,exp_groups,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt);
    count_group = histRead(allqt,exp_groups,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt,1);
else
    count_gene = histRead(allqt,exp_tags,expected_list,...
        ['QT-' num2str(quality_threshold)],homomer_index_allqt,1);
end
drawnow;
toc

%% output files
disp('writing files..');
tic
% write code_n_count file
%--------------------
fid = fopen([output_directory_decode output_filename_afterQT_prefix...
    '_code_n_count' '.csv'], 'w');
if grouped_YN
    fprintf(fid,'Code,Count,GeneName,Group\n');
    temp_write = [cellstr(uni_reads_letters_allqt),...
        num2cell(count_code_allqt'),...
        uni_reads_name_tag_allqt,...
        uni_reads_group_tag_allqt];
    for row = 1:size(temp_write,1)
        fprintf(fid,'%s,%d,%s,%s\n',temp_write{row,:});
    end
else
    fprintf(fid,'Code,Count,GeneName\n');
    temp_write = [cellstr(uni_reads_letters_allqt),...
        num2cell(count_code_allqt'),...
        uni_reads_name_tag_allqt];
    for row = 1:size(temp_write,1)
        fprintf(fid,'%s,%d,%s\n',temp_write{row,:});
    end
end
fclose(fid);

% write gene_n_count file
%--------------------
fid = fopen([output_directory_decode output_filename_afterQT_prefix...
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
    fid = fopen([output_directory_decode...
        output_filename_afterQT_prefix '_group_n_count' '.csv'], 'w');
    fprintf(fid,'GeneGroup,Count\n');
    for row=1:length(count_group(:,1))
        fprintf(fid,'%s,%d\n',count_group{row,:});
    end
    fclose(fid);
end

% write details file
%--------------------
tempheader = {'letters','name','group',...
    'tile_X_pos','tile_Y_pos','global_X_pos','global_Y_pos',...
    'parent_cell','tile_ID',...
    'general_stain_min','seq_quality_min','alignment_score_min'};
tempwrite = [cellstr(letters_allqt),name_tag_allqt,group_tag_allqt,...
    num2cell([tile_x_allqt,tile_y_allqt,global_x_allqt,global_y_allqt,...
    cell_allqt,tile_allqt,...
    general_strength_min_allqt,seq_quality_min_allqt,alignment_score_min_allqt])];

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

fid = fopen([output_directory_decode output_filename_afterQT_prefix...
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

toc

% save mat file
%--------------------
clear -regexp ^temp
clear fid i
disp('saving workspace variables..');
save([output_directory_decode output_filename_afterQT_prefix '.mat']);

fprintf('Thresholding finished.\n\n');
end

