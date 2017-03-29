function    Plotting_global_Sequencing_2(output_directory_decode,...
    output_filename_afterQT_prefix,background_image,resize_factor,...
    taglist_plot,use_default_symbol_list_YN,symbol_size,...
    plot_reads_beforeQT_YN,exclude_NNNN_YN,...
    plot_based_on_group_YN,plot_base1_general_stain,plot_on_blank,varargin)

% Plot reads based on the given taglist_plot.
%
% The background image specified has to be full size, but can be of lower
% quality.
%   Resize factor: the image will be resized based on the factor. It should
%   be between 0 and 1.
%
% Taglist_plot: the reads will be regrouped and curated again based on this
% taglist. It will only affect plotting, not the original decoding and
% thresholding results.
% If not specified, the same taglist as in decoding will be used.
%   Use old style: set to 1 if using the digit column exits in the
%   taglist. It will be ignored.
%
%   If you are too lazy to set the symbols for the taglist, set
%   use_default_symbol_list to 1. A premade symbol list will be used.
%
% Symbol size: default 6. Adjust based on your own need.
%
% Output directory plot: the plotted images will be saved in as a subfolder
% of current folder, with the specified name.
% Output filename plot prefix: filename prefix.
%
% Option1: exclude NNNN.
%   Not to plot unexpected reads.
%
% Option2: plot reads before QT
%
% Option3: based on group.
%   Plot and give the legend name as in the group tag.
%
% Option4: plot base1 general stain.
%   Might be helpful when manually correcting cell-blob relation or see the
%   image position shift.
%
% Plotting_global_Sequencing v2.5.1.1, lite
% Xiaoyan, 2016-2-21

%% initiate
drawnow;
disp('Initiating Plotting_global_Sequencing.');

output_directory_decode = [output_directory_decode '\'];
load([output_directory_decode output_filename_afterQT_prefix '.mat']);

if ~calculate_global_position_YN
    error('Global coordinates are not calculated. Check settings in Decoding.');
end

%% check and reform taglist
[taglist,expected_list,~,~,exp_tags,exp_groups] =  ...
    maketag(taglist_plot,num_hybs,0,0,...
    plot_based_on_group_YN,abnormal_sequencing_YN,sequencing_order);

%% prepare symbol list
if use_default_symbol_list_YN
    sym = default_symbolist;
else
    n = size(taglist,2);
    if n==3
        sym = taglist(:,3);
    else
        disp('No symbol column detected in taglist_plot. Continue with default symbols.');
        sym = default_symbolist;
    end
end
n = 1;
while length(sym) < length(taglist)
    sym = [sym; sym(n)];
    n = n+1;
end

%% prepare the read list
if plot_reads_beforeQT_YN
    list = allbt; global_x_plot = (global_x_allbt-1/resize_factor/2+.5)*resize_factor+1; global_y_plot = (global_y_allbt-1/resize_factor/2+.5)*resize_factor+1;
else
    list = allqt; global_x_plot = (global_x_allqt-1/resize_factor/2+.5)*resize_factor+1; global_y_plot = (global_y_allqt-1/resize_factor/2+.5)*resize_factor+1;
end

%% check image size
if plot_on_blank
else
    ii = imfinfo(background_image);
    if ii.Width<max(global_x_plot) || ii.Height<max(global_y_plot)
        error('Image size smaller than maximum coordinates. Check image again.');
    end
end

%% merge reads with same gene name/gene group
[list,expected_list,unitag,sym] = ...
    groupreadsforplotting(list,expected_list,exp_tags,exp_groups,plot_based_on_group_YN,sym);

%% import background image
if plot_on_blank
else
    if exist(background_image,'file')
        disp('working hard to load image..');
        I = imread(background_image);
    else
        error('Could not find the background image.');
    end
end

%% start plotting
disp('start plotting..');
f = figure;
axis image;
set(f,'units','normalized','position',[.2 .05 .75 0.85]);
fa = get(f,'Children');
set(fa,'position',[0.025 0.05 0.95 0.9]);

if plot_on_blank
    set(fa,'YDir','reverse');
else
    imshow(I,[]);
end
hold on;

n = 0;
if plot_base1_general_stain
    plot(global_x_plot,global_y_plot,'c*','markersize',symbol_size);
    data = length(list);
    n = n+length(list);
    plot_tag = 'general stain';
    rowname = 'general_stain';
else
    plot_tag = []; plot_tag_empty = [];
    A=[]; B=[]; count=[];
    for i = 1:length(expected_list)
        a = find(list == expected_list(i));
        if a
            plot(fa,global_x_plot(a),global_y_plot(a),...
                sym{i},'markersize',symbol_size);
            n = n+length(a);
            plot_tag = [plot_tag unitag(i)];
            A = [A; a]; count = [count; length(a)];
        else
            plot_tag_empty = [plot_tag_empty unitag(i)];
            B = [B; 0];
        end
    end
    
    % plot the rest of the reads as NNNN
    if exclude_NNNN_YN
    else
        if length(A) ~= length(list)
            global_x_N = removerows(global_x_plot,'ind',A);
            global_y_N = removerows(global_y_plot,'ind',A);
            plot(global_x_N,global_y_N,'ch');
            plot_tag = [plot_tag 'NNNN'];
            count = [count; length(global_x_N)];
            n = n+length(global_x_N);
        else
            plot_tag_empty = [plot_tag_empty 'NNNN'];
            B = [B; 0];
        end
    end
    if isempty(B)
        data = count; rowname = [plot_tag];
    else
        data = [count;B]; rowname = [plot_tag, plot_tag_empty];
    end
end

% axis and legend
axis(fa,'off');
legend(fa,plot_tag,'location','NorthEastOutside','color',[.6 .6 .6]);

% table to show counts
f2 = figure;
set(f2,'units','normalized','position',[0.05 0.2 0.14 0.4]);
h = uitable(f2,'data',data,'ColumnName',{'count'},'RowName',rowname);
set(h,'units','normalized','position',[0 0 1 1]);

show = ['Number of blobs plotted: ' num2str(n) '.'];
disp(show);
fprintf('Plotting finished.\n\n');

end




