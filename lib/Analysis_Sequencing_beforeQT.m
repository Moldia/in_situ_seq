function Analysis_Sequencing_beforeQT(output_directory_decode,...
        output_filename_decode_prefix,...
        narrow_down_quality_range_YN,lower_limit,upper_limit,...
        table_gene_counts_at_different_thresholds_YN,threshold,...
        plot_general_quality_items_YN,plot_seq_spec_channel_histogram_YN,...
        guess_closest_expected_reads_YN,plot_quality_vs_general_stain_YN,...
        plot_spectrum_YN,varargin)

% Preliminery analysis on data before QT.
% Help to set thresholds properly in thresholding.
%
% Analysis_Sequencing_beforeQT v2.2
% Xiaoyan, 2015-6-30

%% initiate
drawnow;
disp('Initiating Analysis_Sequencing_beforeQT.');

switch length(varargin)
    case 1
        guessTbase = varargin{1};
    otherwise
        guessTbase = 0;
end
output_directory_decode = [output_directory_decode '\'];
load([output_directory_decode output_filename_decode_prefix '.mat']);
output_prefix = [output_directory_decode,output_filename_decode_prefix,...
    '_analysis\'];
if exist(output_prefix, 'dir')  
else mkdir (output_prefix);
end 


%%
if narrow_down_quality_range_YN
    [thres,count] = ...
        barQuality_2(lower_limit,upper_limit,catagory_read_allbt,seq_quality_min_allbt,1);
    fid = fopen([output_prefix 'QualityBar' '.csv'], 'w');
    fprintf(fid,'threshold,expected,unexpected,homomer,belowQT\n');
    temp_write = [thres',count];
    fprintf(fid,'%.3f,%d,%d,%d,%d\n',temp_write');
    fclose(fid);
end

%%
if table_gene_counts_at_different_thresholds_YN
    tableGene_2(threshold,seq_quality_min_allbt,catagory_read_allbt,exp_tags,output_prefix);
end

%%
if plot_general_quality_items_YN
    fh = figure;
    histRCP(general_strength_allbt,fh,output_prefix);
    histQuality(seq_quality_allbt,fh,output_prefix);
    set(fh,'name','general stain and quality score histogram');
    fh2 = figure;
    [channelhist,b] = histChannel(channel_strength_max_allbt,...
        seq_res_allbt,fh2,1,output_prefix);
    set(fh2,'name','sequencing channel strength histogram - base called');
end

%%
if plot_seq_spec_channel_histogram_YN
    fh3 = figure;
    if exist('channelhist','var')==1
    else
        [channelhist,b] = histChannel(channel_strength_max_allbt,seq_res_allbt,fh3,0,output_prefix);
    end
    seqchannelhist = histSeqChannel(channel_strength_original_allbt,...
        channelhist,b,fh3,1);
    set(fh3,'name','sequencing channel strength histogram - all raw signals');
end

%%
if guess_closest_expected_reads_YN
    guessRead(expected_digits(:,2:end),allbt,name_tag_allbt,taglist(:,1),output_prefix);
end

%%
if plot_quality_vs_general_stain_YN
    fh4 = figure;
    heatThreshold(seq_quality_allbt,general_strength_allbt,fh4,0);
end

%%
if plot_spectrum_YN
    fh5 = figure;
    scatterSpect(seq_res_allbt,channel_strength_original_allbt,fh5,1,0,0);
end

%%

if guessTbase
    [Tdata,Tcolumnn,Trown] = ...
        GuessT(expected_digits(:,2:end),allbt,name_tag_allbt,taglist(:,1));
    fid = fopen([output_prefix 'GuessT' '.csv'], 'w');
    fprintf(fid,'%s\n','This is a table showing the guessing of T');
    fprintf(fid,'');
    for j = 1:length(Tcolumnn)
        fprintf(fid,',%s',Tcolumnn{j});
    end
    fprintf(fid,'\n');
    for i = 1:length(Trown)
        fprintf(fid,'%s,%d',Trown{i},Tdata{i,1});
        for j = 2:length(Tcolumnn)
            fprintf(fid,',%s',Tdata{i,j});
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

fprintf('Analysis beforeQT finished.\n\n'); 
end
