function heatThreshold(seq_quality_list,general_strength_list,fh1,createheat)

% Plot general stain vs. seq quality (two items used in thresholding).
% Give two plots. One scatter plot, and the other heat map.
%
% dotThreshold v2.1
% Apr 28, 2014, Xiaoyan

if createheat
    fh2 = figure;
end

bmax = length(general_strength_list(1,:));
x_min = min(seq_quality_list(:));
x_max = max(seq_quality_list(:));
y_min = min(general_strength_list(:));
y_max = max(general_strength_list(:));
cmax = length(general_strength_list(:,1))/100;


for i = 1:bmax
    figure(fh1); sh(i) = subplot(2,ceil(bmax/2),i);
    ph(i) = scatter(sh(i),seq_quality_list(:,i),general_strength_list(:,i),7,'ks');
    if createheat
        figure(fh2); sh2(i) = subplot(2,ceil(bmax/2),i);
        hist3([seq_quality_list(:,i),general_strength_list(:,i)],[50 50]);
        view(2);
        set(get(sh2(i),'child'),'FaceColor','interp','CDataMode','auto');
        colormap(jet), caxis([0 cmax]); grid off; colorbar;
        xlabel(sh2(i),['seq quality - base' num2str(i)]);
        ylabel(sh2(i),'general stain');
    end
    xlabel(sh(i),['seq quality - base' num2str(i)]);
    ylabel(sh(i),'general stain');
end


%axis(sh(:),[x_min-0.05 x_max+0.05 y_min-0.01 y_max+0.01]);
if createheat
    set(sh2(:),'xTick',0:0.2:1,'yTick',0:0.2:1);
    set(ph2(:),'markeredgecolor',rgb('royalblue'));
    set(fh2,'name','general stain vs seq quality - heat map');
end
set(sh(:),'xTick',0:0.2:1,'yTick',[0.01,0.05,0.1,0.2:0.2:1]);
set(fh1,'name','general stain vs seq quality - scatter plot');

end