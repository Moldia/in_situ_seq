function scatterSpect(seq_res,channel_strength_list,fh1,...
    createfigure,createheat,mergehybs)

% Plot "spectrum" of blobs in each cycle, per base
%
% scatterSpect v2.2
% May 17, 2014, Xiaoyan


bmax = length(seq_res(1,:));
blobn = length(seq_res(:,1));
y_max = zeros(4,1)+0.001;
yinter = (1-0.15)/4;
xinter = (1-0.15)/bmax;
cmax = 20;

color = {'b',rgb('forestgreen'),'r',rgb('darkturquoise')};
if createheat
    fh2 = figure;
end
channel_spect=struct;
% i for base position, j for nucleotide
for i = 1:bmax
    A=[];
    for j = 1:4
        M = find(seq_res(:,i)==j);
        if M
            channel = zeros(length(M),4);
            x = zeros(length(M),4);
            for k = 1:4
                if size(channel_strength_list,2)>=k
                    channel(:,k) = channel_strength_list(i,k,M);
                else
                    channel(:,k) = zeros(length(M),1)-1;
                end
                x(:,k) = k-0.2 + 0.4.*rand(length(M),1);
            end
        else
            channel = zeros(1,4)-1;
            x = zeros(1,4);
        end
        if createfigure
            % create subplot handle
            figure(fh1); sh((j-1)*bmax+i) = subplot('position',...
                [0.1+(0.01+xinter)*(i-1)...
                1-0.01-yinter-(0.03+yinter)*(j-1)...
                xinter yinter]);
            plotSpread(channel,'distributionColors',color);
            
            if createheat
                figure(fh2);
                sh2((j-1)*bmax+i) = subplot('position',...
                    [0.1+(0.01+xinter)*(i-1)...
                    1-0.01-yinter-(0.03+yinter)*(j-1)...
                    xinter yinter]);
                x_heat = [x(:);-1;6];channel_heat = [channel(:);-1;2];
                hist3([x_heat,channel_heat],[70 70]);
                axis([-0.5 5.5 0 1]);
                set(get(sh2((j-1)*bmax+i),'child'),'FaceColor','interp','CDataMode','auto');
                view(2); colormap(jet); grid off;
            end
            
            if (max(channel(:))>y_max(j))
                y_max(j) = max(channel(:));
            end
            
        end
    end       
end

if mergehybs
    fh3 = figure;
    for j = 1:4
        [M,N] = find(seq_res==j);
        if M
            channelmerge = zeros(length(M),4);
            x = zeros(length(M),4);
            for k = 1:4
                if size(channel_strength_list,2)>=k
                    for l = 1:length(M)
                        channelmerge(l,k) = channel_strength_list(N(l),k,M(l));
                    end
                else
                    channelmerge(:,k) = zeros(length(M),1)-1;
                end
                x(:,k) = k-0.2 + 0.4.*rand(length(M),1);
            end
        else
            channelmerge = zeros(1,4)-1;
            x = zeros(1,4);
        end
        
        figure(fh3),sh3(j)=subplot(2,2,j); hold on;
        plotSpread(channelmerge,'distributionColors',color);
        xlim([-0.5 5.5]);
        xlabel(num2letters(j));
        set(gca,'xTick',1:4,'xTickLabel',{'A' 'C' 'G' 'T'},...
            'yTick',[0.01,0.1,0.2:0.2:max(y_max)]);
        hold off;
    end
    
end

% set xTick,yTick,xlable,ylabel
if createfigure
    set(sh(:),'xTick',1:4,'yTick',[0.01,0.1,0.2:0.2:max(y_max)],...
        'box','off','xTickLabel',{'A' 'C' 'G' 'T'},...
        'fontsize',8,'fontweight','bold');
    for i = 1:length(sh)
        switch mod(i,bmax)
            case 1
                ylabel(sh(i),num2letters(floor((i-1)/bmax)+1),...
                    'fontsize',10,'fontweight','bold');
            otherwise
                set(sh(i),'yTickLabel',[]);
        end
    end
    for i = 1:length(sh)
        switch floor((i-1)/bmax)
            case 0
                axis(sh(i),[-0.5 5.5 0 y_max(1)*1.1]);
            case 1
                axis(sh(i),[-0.5 5.5 0 y_max(2)*1.1]);
            case 2
                axis(sh(i),[-0.5 5.5 0 y_max(3)*1.1]);
            case 3
                xlabel(sh(i),['base' num2str(i-bmax*3)],'fontsize',10);
                axis(sh(i),[-0.5 5.5 0 y_max(4)*1.1]);
        end
        
    end
    set(fh1,'name','spectrum characteristics - scatter plot','visible','on');

    if createheat
        set(sh2(:),'xTick',1:4,'yTick',[0.01,0.1,0.2:0.2:max(y_max)],...
            'box','off','xTickLabel',{'A' 'C' 'G' 'T'},...
            'fontsize',8,'fontweight','bold');

        for i = 1:length(sh2)
            switch mod(i,bmax)
                case 1
                    ylabel(sh2(i),num2letters(floor((i-1)/bmax)+1),...
                        'fontsize',10,'fontweight','bold');
                otherwise
                    set(sh2(i),'yTickLabel',[]);
            end
        end
        for i = 1:length(sh2)
            switch floor((i-1)/bmax)
                case 0
                    axis(sh2(i),[-0.5 5.5 0 1]);
                case 1
                    axis(sh2(i),[-0.5 5.5 0 1]);
                case 2
                   axis(sh2(i),[-0.5 5.5 0 1]);
                case 3
                    xlabel(sh2(i),['base' num2str(i-bmax*3)],'fontsize',10);
                    axis(sh2(i),[-0.5 5.5 0 1]);
            end
        end
        set(fh2,'name','spectrum characteristics - heat map','visible','on');
    end

end
    
end