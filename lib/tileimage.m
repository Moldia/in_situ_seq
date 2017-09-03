function tileimage(im, tilesize, outdir)
% tileimage(im, tilesize, outdir)
% create non-overlapping image tiles
% Xiaoyan, 2017

imsize = size(im);

ntileX = ceil(imsize(2)/tilesize);
ntileY = ceil(imsize(1)/tilesize);

% pad image
im = padimg(im, tilesize*ntileX-imsize(2), tilesize*ntileY-imsize(1), 'SE');

% tile and save
mkdir(outdir);
for i = 1:ntileY
    for j = 1:ntileX
        xstart = tilesize*(j-1);
        ystart = tilesize*(i-1);
        tile = im(ystart+1:ystart+tilesize, xstart+1:xstart+tilesize, :);
        imwrite(tile,...
            [outdir, '\tile' num2str(ntileX*(i-1)+j), '.tif']);
    end
end

end
