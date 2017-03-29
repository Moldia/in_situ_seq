% Integrated analysis for sequencing.
% For a specific module, if you don't know how to set the parameteres,
% eg. Tiling_Sequencing, type "help Tiling_Sequencing" in the command
% window.
% Ideally keep one file for one experiment
% Five major functions: 1)Tiling, 2)Decoding, 3)Preliminary analysis before QT,
% 4)Thresholding, 5)Plotting 
% Carefully set the paramets.
% All variables ending with _YN are yes or no (1 or 0) questions.
%
% Sequencing v3
% Xiaoyan, 2015-6-30


clear, clc; close all; drawnow;

% Choose functions to run
run_Tiling_YN = 0;
run_Decode_YN = 1;
run_Analysis_YN = 0;
run_Threshold_YN = 1;
run_Plotting_Global_YN = 1;
%================================================
% set parameters
%----------------------------
% Tiling_Sequencing
    folder_image = 'D:\Homer\Duff'; % preferably full path name
    filename_base_prefix = 'base';  % keep single quote marks
        in_subfolder_YN = 0;
    filename_channel_prefix = '_s3c';
    filename_suffix = '_ORG.tif';
    base_max = 4;       channel_max = 6;
    x_size = 1000;      y_size = 1000;
    % options
    create_CSV_file_YN = 1; 
        CSV_filename_prefix = 'Tiled';
        channel_order = {'Nuclei' 'T' 'G' 'C' 'A' 'General_stain'};
%----------------------------
% Decoding_Sequencing
    input_file = 'blobs.csv';
    % don't change this unless your input file has some weird form
    General_Alignment_ACGT_column_number = [4,5,6,7,8,9];    % use 0 if any of them is MISSING in the file
    XYPosition_ParentCell_column_number = [10,11,12];
    
    num_hybs = 4;  % if 5th empty cycle exists, num_hybs=5, no upper limit
    taglist = ID_empty;
        grouped_YN = 0;
    calculate_global_position_YN = 1;
        csv_file_contain_tile_position = 'TilePosition.csv'; %full path if not in the current directory
    output_directory_decode = 'Decoding';   
    output_filename_decode_prefix = 'beforeQT';
    % options
    check_parent_cell_YN = 0;       
    check_alignment_YN = 0;
        alignment_min_threshold = 1.8;
    abnormal_sequencing_YN = 0;
        sequencing_order = '1234';  % keep the quote marks, same length as (num_hybs - cycle5_empty_YN)
%----------------------------
% Analysis_Sequencing_beforeQT
    narrow_down_quality_range_YN = 0;
        lower_limit = 0.3;        upper_limit = 0.6;
    table_gene_counts_at_different_thresholds_YN = 0;
        threshold = [0.1, 0.3, 0.4, 0.8];
    plot_general_quality_items_YN = 0;
    plot_seq_spec_channel_histogram_YN = 0;
    guess_T_YN = 1;
    guess_closest_expected_reads_YN = 0;
    plot_quality_vs_general_stain_YN = 0;
    plot_spectrum_YN = 0;
%----------------------------
% Threshold_Sequencing
    quality_threshold = 0.4;        general_strain_threshold = 0.001;
    output_filename_afterQT_prefix = ['QT_' num2str(quality_threshold) '_' num2str(general_strain_threshold)];
%----------------------------
% Plotting_global_Sequencing
    background_image = 'D:\Homer\Duff\Image_20%.jpg'; 
        scale = 0.2;
    I_want_to_plot_on_white_backgound = 0;
    taglist_plot = taglist; % default: use the same one as in Decoding
        use_default_symbol_list_YN = 0;
    symbol_size = 8; % default: 6
    % options
    exclude_NNNN_YN = 0;
    plot_reads_beforeQT_YN = 0;
    plot_based_on_group_YN = 0;
    plot_base1_general_stain = 0; 
%================================================
if run_Tiling_YN || run_Decode_YN || run_Analysis_YN || run_Threshold_YN || run_Plotting_Global_YN
else
    error('Choose at least one function.');
end

if run_Tiling_YN
    Tiling_Sequencing_2(folder_image,filename_base_prefix,in_subfolder_YN,...
        filename_channel_prefix,filename_suffix,base_max,channel_max,...
        x_size,y_size,create_CSV_file_YN,CSV_filename_prefix,channel_order);
end
if run_Decode_YN
    Decoding_Sequencing_2(input_file,num_hybs,taglist,grouped_YN,...
        calculate_global_position_YN,csv_file_contain_tile_position,...
        output_directory_decode,output_filename_decode_prefix,...
        check_parent_cell_YN,check_alignment_YN,alignment_min_threshold,...
        abnormal_sequencing_YN,sequencing_order,...
        General_Alignment_ACGT_column_number,XYPosition_ParentCell_column_number);
end
if run_Analysis_YN
    Analysis_Sequencing_beforeQT(output_directory_decode,...
        output_filename_decode_prefix,...
        narrow_down_quality_range_YN,lower_limit,upper_limit,...
        table_gene_counts_at_different_thresholds_YN,threshold,...
        plot_general_quality_items_YN,plot_seq_spec_channel_histogram_YN,...
        guess_closest_expected_reads_YN,plot_quality_vs_general_stain_YN,...
        plot_spectrum_YN,guess_T_YN);
end
if run_Threshold_YN
    Threshold_Sequencing_2(output_directory_decode,...
        output_filename_decode_prefix,output_filename_afterQT_prefix,...
        quality_threshold,general_strain_threshold);
end 
if run_Plotting_Global_YN
    Plotting_global_Sequencing_2(output_directory_decode,...
        output_filename_afterQT_prefix,background_image,scale,...
        taglist_plot,use_default_symbol_list_YN,symbol_size,...
        plot_reads_beforeQT_YN,exclude_NNNN_YN,...
        plot_based_on_group_YN,...
        plot_base1_general_stain,I_want_to_plot_on_white_backgound)
end

clear;   
    