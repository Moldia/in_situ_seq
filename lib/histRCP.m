function out = histRCP(general_strength_list,fh,output_prefix)

% Plot histogram of RCP intensities of each base
% Output a csv file
% 
% histRCP v2.1
% May 13, 2014, Xiaoyan

bmax = length(general_strength_list(1,:));
x_max = max(general_strength_list(:));
x_min = min(general_strength_list(:));
y_max = 0;
A = [];

figure(fh);
for i = 1:bmax
    if x_min==x_max
        a = length(general_strength_list(:,i));
        b = x_min;
    else
        [a,b] = hist(general_strength_list(:,i),linspace(x_min,x_max,31));
    end
    inter = (1-0.15)/bmax;
    
    % create subplot handle
    sh(i) = subplot('position',[0.1+(0.01+inter)*(i-1) 0.56 inter 0.4]);
    bh(i) = bar(sh(i),b,a);
    xlabel(['base' num2str(i)],'fontsize',8,'fontweight','bold');
    
    if max(a)>y_max
        y_max = max(a);
    end
    A = [A a'];
end

set(bh(:),'facecolor',rgb('royalblue'));

set(sh(2:end),'yTickLabel',[]);
if x_min==x_max
    axis(sh(:),[-1 2 1 y_max*1.1]);
    set(sh(:),'xTick',x_min,...
        'yTick',[20,100,200,500:500:y_max],'box','off');
else
    axis(sh(:),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 1 y_max*1.1]);
    set(sh(:),'xTick',[0.01,0.1,0.2:0.2:x_max],...
        'yTick',[20,100,200,500:500:y_max],'box','off');
end
ylabel(sh(1),'general stain','fontsize',10,'fontweight','bold');


fid = fopen([output_prefix 'GeneralStainHistogram' '.csv'], 'w');
fprintf(fid,'%s\n','This is a table showing the distribution of general stain intensities');
fprintf(fid,'bin');
for i = 1:bmax
    fprintf(fid,',%s',['hyb' num2str(i)]);
end
fprintf(fid,'\n');
for i = 1:length(b)
    fprintf(fid,'%.3f',b(i));
    for j = 1:bmax
        fprintf(fid,',%d',A(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

out = [b' A];

end