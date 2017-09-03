function im = padimg(im, deltaX, deltaY, varargin)
% im = padimg(im, deltaX, deltaY, varargin)
% pad image with zeros
% if varargin not specified, pad right and lower
% Xiaoyan, 2017


imsize = size(im);
try 
    imsize(3);
catch
    imsize(3) = 1;
end

if deltaX > 0
    padX = zeros(imsize(1), deltaX, imsize(3));
end

if deltaY > 0
    padY = zeros(deltaY, imsize(2)+deltaX, imsize(3));
end

try
    side = varargin{1};
catch
    side = 'ES';
end

% pad/crop different side
if strfind(side, 'E')
    try im = [im, padX];
    catch
        im = im(:,1:imsize(2)+deltaX,:);
    end
end

if strfind(side, 'W')
    try im = [padX, im];
    catch
        im = im(:,abs(deltaX)+1:end,:);
    end
end

if strfind(side, 'S')
    try im = [im; padY];
    catch errorm
        if strfind(errorm.identifier, 'dimensionMismatch')
            padY = zeros(deltaY, imsize(2)+2*deltaX, imsize(3));
            im = [im; padY];
        else
            im = im(1:imsize(1)+deltaY,:,:);
        end
    end
end

if strfind(side, 'N')
    try im = [padY; im];
    catch errorm
        if strfind(errorm.identifier, 'dimensionMismatch')
            padY = zeros(deltaY, imsize(2)+2*deltaX, imsize(3));
            im = [im; padY];
        else
            im = im(abs(deltaY)+1:end,:,:);
        end
    end
end

end
