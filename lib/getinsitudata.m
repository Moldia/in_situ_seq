function [name, pos, varargout] = getinsitudata(filename, varargin)
% [name, pos, varargout] = getinsitudata(filename, varargin)
% get transcript info from decoded file
% Xiaoyan, 2107

switch length(varargin)
    case 0
        colName = 2;
        colPos = 1;
    case 1
        colName = varargin{1};
        colPos = 1;
    otherwise
        colName = varargin{1};
        colPos = varargin{2};
end

data = importdata(filename, ',', 1);
name = data.textdata(2:end,colName);
pos = data.data(:,colPos:colPos+1);

for i = 1:length(varargin)-2
    varargout{i} = data.data(:,varargin{i+2});
end
    
end