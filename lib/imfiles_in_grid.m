function [imdir, imfiles, imOrder] = imfiles_in_grid(imfolder, imprefix, imsuffix,...
    ndigits, ntiles, varargin)
% [imdir, imfiles] = imfiles_in_grid(imfolder, imprefix, imsuffix,...
%     ndigits, ntiles, varargin)
% return a cell matrix with each cell containing tiled image file name
% Xiaoyan, 2017


% tile directions, default: start from upper left, parallel
if isempty(varargin)
    order = 'EES';
else
    order = varargin{1};
end

if ntiles(1)*ntiles(1) > 1
    % first dimension
    imOrder = zeros(ntiles(2), ntiles(1));
    switch order(1)
        case 'E'
            imOrder(1,:) = 1:ntiles(1);
            trans = 0;
        case 'W' 
            imOrder(1,:) = ntiles(1):-1:1;
            trans = 0;
        case 'S'
            imOrder(:,1) = 1:ntiles(2);
            imOrder = imOrder';
            trans = 1;
        case 'N'
            imOrder(:,1) = ntiles(2):-1:1;
            imOrder = imOrder';
            trans = 1;
    end

    if strcmp(order(1), order(2))
        imOrder(2,:) = imOrder(1,:);            % parallel
    else
        imOrder(2,:) = fliplr(imOrder(1,:));    % snake
    end

    try imOrder = repmat(imOrder(1:2,:), size(imOrder,1)/2, 1);
    catch
        imOrder = [repmat(imOrder(1:2,:), floor(size(imOrder,1)/2), 1); imOrder(1,:)];
    end

    % second dimension
    if strcmp(order(end), 'S') || strcmp(order(end), 'E')
        rowstart = 1:size(imOrder,1);
    else
        rowstart = size(imOrder,1):-1:1;
    end
    rowstart = (rowstart(:)-1)*size(imOrder,2);
    imOrder = imOrder + repmat(rowstart, 1, size(imOrder,2));

    if trans
        imOrder = imOrder';
    end
else
    % one and only one tile in the file
    imOrder = 1;
end

imfiles = catstrnum(imprefix, imOrder, ndigits);
imfiles = strcat(imfiles, imsuffix);
imfiles = reshape(imfiles, ntiles(2), ntiles(1));

try
    if ~strcmp(imfolder(end), filesep)
        imfolder = [imfolder filesep];
    end
end
imdir = repmat({imfolder}, ntiles(2), ntiles(1));

end

