function posNew = correctcoord(pos,resize_factor)
% posNew = correctcoord(pos,resize_factor)
% correct coordinates for plotting and image processing
% Xiaoyan, 2017


posNew = (pos-1/resize_factor/2+.5)*resize_factor + 1;

end
