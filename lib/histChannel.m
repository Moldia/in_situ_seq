function [channel_hist,b] = histChannel(channel_strength_list,seq_res,...
    fh,createfigure,output_prefix)

% Plot separate histograms for each nucleotide for each base (called
% bases)w
% Ouput a csv file
% 
% histChannel v2.2
% May 3, 2014, Xiaoyan

bmax = length(seq_res(1,:));

x_max = max(channel_strength_list(:));
x_min = min(channel_strength_list(:));
y_max = zeros(4,1);
yinter = (1-0.15)/4;
xinter = (1-0.15)/bmax;
bh=[];

figure(fh);
channel_hist=struct;
% i for base position, j for nucleotide
for i = 1:bmax
    A=[];
    for j = 1:4
        M = find(seq_res(:,i)==j);
        if M
            channel = channel_strength_list(M,i);
        else
            channel = zeros(1);
        end    
        [a,b] = hist(channel,linspace(x_min,x_max,101));
%         [a,b] = hist(channel,50);
        A = [A,a'];
        if createfigure
            % create subplot handle
            sh((j-1)*bmax+i) = subplot('position',...
                [0.1+(0.01+xinter)*(i-1)...
                1-0.01-yinter-(0.03+yinter)*(j-1)...
                xinter yinter]);
            bh((j-1)*bmax+i) = bar(b,a);

            if max(a)>y_max(j)
                y_max(j) = max(a);
            end
        end
    end
    field = ['base' num2str(i)];
    channel_hist = setfield(channel_hist,field,A);
end

% set xTick,yTick,xlable,ylabel 
if createfigure
    set(sh(:),'xTick',[0.01,0.1,0.2:0.2:x_max],'yTick',[10,50,100,200:200:max(y_max)],'box','off');
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
                set(bh(i),'facecolor',rgb('indigo'));
                set(sh(i),'xTickLabel',[]);
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(1)*1.1]);
            case 1
                set(bh(i),'facecolor',rgb('royalblue'));
                set(sh(i),'xTickLabel',[]);
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(2)*1.1]);
            case 2
                set(bh(i),'facecolor',rgb('forestgreen'));
                set(sh(i),'xTickLabel',[]);
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(3)*1.1]);
            case 3
                set(bh(i),'facecolor',rgb('firebrick'));
                xlabel(sh(i),['base' num2str(i-bmax*3)],...
                    'fontsize',10,'fontweight','bold');
                axis(sh(i),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 ...
                    1 y_max(4)*1.1]);
        end
    end
end



fid = fopen([output_prefix 'SeqChannelHistogram' '.csv'], 'w');
fprintf(fid,'%s\n',['This is a table showing the distribution of '...
    'sequencing channel intensities'...
    ' (only the subset of data when the corresponding base is called']);
for i = 1:bmax
    fprintf(fid,'\n%s\n',['hyb' num2str(i)]);
    fprintf(fid,'bin,A,C,G,T\n');
    field = ['base' num2str(i)];
    printhis = getfield(channel_hist,field);
    for j = 1:length(b)
        fprintf(fid,'%.3f,%d,%d,%d,%d\n',b(j),printhis(j,:));
    end
end
fclose(fid);

end
