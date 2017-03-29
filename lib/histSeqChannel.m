function channel_hist = histSeqChannel(channel_strength_list,...
    channelhist,bin,fh,createfigure)

% Plot histogram of sequencing specific channels
% Raw data in solid color (called and non-called)
% Called signals in transparent front layer
%
% histSeqChannel v2.2
% 2015-6-29, Xiaoyan

bmax = length(channel_strength_list(:,1,1));

x_max = max(bin);
x_min = min(bin);
y_max = zeros(4,1);
yinter = (1-0.15)/4;
xinter = (1-0.15)/bmax;
bh=[];
bh2 = [];

figure(fh);
channel_hist=struct;
% i for base position, j for nucleotide
for i = 1:bmax
    A=[];
    field = ['base' num2str(i)];
    channelcalled = getfield(channelhist,field);
    for j = 1:4
        if size(channel_strength_list,2)>=j
            channel = channel_strength_list(i,j,:);
        else
            channel = zeros(1);
        end
        [a,b] = hist(channel(:),bin);
%         [a,b] = hist(channel(:),50);
        A = [A,a'];
        if createfigure
            % create subplot handle
            sh((j-1)*bmax+i) = subplot('position',...
                [0.1+(0.01+xinter)*(i-1)...
                1-0.01-yinter-(0.03+yinter)*(j-1)...
                xinter yinter]);
            bh((j-1)*bmax+i) = bar(sh((j-1)*bmax+i),b,a); hold on;
            bh2((j-1)*bmax+i) = bar(sh((j-1)*bmax+i),b,channelcalled(:,j));
            hold off;
            if max(a)>y_max(j)
                y_max(j) = max(a);
            end
        end
    end
    channel_hist = setfield(channel_hist,field,A);
end


% set xTick,yTick,xlable,ylabel      
if createfigure
%     set(bh2(:),'facecolor',rgb('white'))
    % doesn't work after R2014b
%     childH = get(bh(:),'child');
%     childH2 = get(bh2(:),'child');
%     set(childH,'facea',0.8);
%     set(cell2mat(childH2),'facealpha',0.3);
    opengl software;
    set(sh(:),'xTick',[0.01,0.1,0.2:0.2:x_max],'yTick',[50,200,500:500:max(y_max)],...
        'box','off');
    for i = 1:length(sh)
        switch mod(i,bmax)
            case 1
                ylabel(sh(i),num2letters(floor((i-1)/bmax)+1),...
                    'fontsize',10,'fontweight','bold');
            otherwise
                set(sh(i),'yTickLabel',[],'YScale','log');
        end
    end
    for i = 1:length(sh)
        switch floor((i-1)/bmax)
            case 0
                set(bh(i),'facecolor',rgb('indigo'));
                set(sh(i),'xTickLabel',[],'YScale','log');
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(1)*1.1]);
            case 1
                set(bh(i),'facecolor',rgb('royalblue'));
                set(sh(i),'xTickLabel',[],'YScale','log');
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(2)*1.1]);
            case 2
                set(bh(i),'facecolor',rgb('forestgreen'));
                set(sh(i),'xTickLabel',[],'YScale','log');
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(3)*1.1]);
            case 3
                set(bh(i),'facecolor',rgb('firebrick'));
                xlabel(sh(i),['base' num2str(i-bmax*3)],...
                    'fontsize',10,'fontweight','bold');
                set(sh(i),'YScale','log');
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(4)*1.1]);
        end
    end
end

end
