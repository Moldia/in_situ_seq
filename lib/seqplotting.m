function seqplotting(decodedir, taglist, qinput, pinput)
% seqplotting(decodedir, taglist, qinput, pinput)
% plotting for in situ sequencing
% input from InSituSequencing script
% Xiaoyan, 2018


drawnow

if pinput.plot_reads_beforeQT_YN || pinput.plot_ref_general_stain
    decodefile = fullfile(decodedir, 'beforeQT_details.csv');
else
    decodefile = fullfile(decodedir,...
        ['QT_' num2str(qinput.quality_threshold) '_' num2str(qinput.general_stain_threshold) '_details.csv']);
end

[name, pos] = getinsitudata(decodefile);
[uNames, ~, iName] = unique(name);

% taglist and symbol
if ischar(taglist)
    if isempty(taglist)
        plotname = uNames;
        plotname = plotname(~strcmp(plotname, 'NNNN'));
    else
        taglist = importdata(taglist);
        taglist = cellfun(@(v) strsplit(v, ','), taglist, 'uni', 0);
        taglist = reshape([taglist{:}], length(taglist{1}), [])';
        plotname = taglist(:,2);
    end        
else
    [taglist, isRandomSym] = formattaglist(taglist);
    taglist = taglist(:,[true true ~isRandomSym]);
    plotname = taglist(:,2);
end

% take only genes that are actually detected
tempidx = find(ismember(plotname, uNames));
plotname = plotname(tempidx);

% remove name duplicates in taglist
plotnameOriginalOrder = unique(plotname, 'stable');
[plotname, idx] = unique(plotname);
idx = tempidx(idx);

% prepare symbols if not provided
try
    sym = taglist(idx, 3);
catch
    sym = repmat(symlist, ceil(length(plotname)/length(symlist)), 1);
    sym = sym(1:length(plotname));
end

% add NNNN to legend as well as to symbols
if ~pinput.exclude_NNNN_YN && ismember('NNNN', uNames)
    plotname = [plotname; {'NNNN'}];
    plotnameOriginalOrder = [plotnameOriginalOrder; {'NNNN'}];
    sym = [sym; {'ch'}];
end

% start plotting
figure; 
if pinput.I_want_to_plot_on_white_backgound
    plotonblank;
else
    I = imread(pinput.background_image);
    imshow(I, [])
    pos = correctcoord(pos, pinput.scale);
end

hold on;
if pinput.plot_ref_general_stain
    plot(pos(:,1), pos(:,2), '.');
    legend({'all blobs'}, 'color', [.6 .6 .6]);
else
    for i = 1:length(plotname)
        subreads = iName == find(strcmp(plotname{i}, uNames));
        plot(pos(subreads,1), pos(subreads,2), sym{i});
    end
    legend(plotname, 'color', [.6 .6 .6]);
    update_legend(gca, plotnameOriginalOrder);
end

end
