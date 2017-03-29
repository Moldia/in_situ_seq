function out = histQuality(seq_quality_list,fh,output_prefix)

% Plot histogram of quality scores of each base
% Output a csv file
%
% histQuality v2.2
% May 4, 2014, Xiaoyan

bmax = length(seq_quality_list(1,:));
x_max = max(seq_quality_list(:));
x_min = min(seq_quality_list(:));
y_max = 0;
A = [];

figure(fh);
for i = 1:bmax
    [a,b] = hist(seq_quality_list(:,i),linspace(x_min,x_max,31));
    inter = (1-0.15)/bmax;
    
    % create subplot handle
    sh(i) = subplot('position',[0.1+(0.01+inter)*(i-1) 0.1 inter 0.4]);
    bh(i) = bar(sh(i),b,a);
    xlabel(['base' num2str(i)],'fontsize',8,'fontweight','bold');
    
    if max(a)>y_max
        y_max = max(a);
    end
    A = [A a'];
end

set(bh(:),'facecolor',rgb('olivedrab'));    
set(sh(:),'xTick',0:0.2:x_max,'yTick',[20,100,200:200:y_max],'box','off');
set(sh(2:end),'yTickLabel',[]);
axis(sh(:),[x_min-(x_max-x_min)/31 x_max+(x_max-x_min)/31 1 y_max*1.1]);
ylabel(sh(1),'quality score','fontsize',10,'fontweight','bold');


fid = fopen([output_prefix 'QualityScoreHistogram' '.csv'], 'w');
fprintf(fid,'%s\n','This is a table showing the distribution of quality scores');
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