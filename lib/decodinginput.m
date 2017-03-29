function seqdata = decodinginput(inputfile, colChannels, colOthers)
% return sequencing data matrix from CP output
% Xiaoyan 2014-11-28


% input check 
temp = 3 + length(colChannels) + length(colOthers);
if temp ~= 12
    error('Wrong column number detected.');
end

% import
seq = importdata(inputfile, ',', 1);
blobs = size(seq.data, 1);

% preallocate
seqdata = zeros(blobs, 12);

% check column numbers
column1_exist = find(colChannels);
column1_missing =  find(colChannels==0);
column2_exist = find(colOthers);
column2_missing = find(colOthers==0);

if colOthers(1)==0 || colOthers(2)==0
    error('X and Y position mush exist.')
end

% reform into output matrix
seqdata(:,1:3) = seq.data(:,1:3);
seqdata(:,column1_exist+3) = ...
    seq.data(:,colChannels(column1_exist));
seqdata(:,column1_missing+3) = seqdata(:,column1_missing+3) + 1E-6;
seqdata(:,10:11) = ...
    seq.data(:,...
    colOthers(column2_exist(1)):...
    colOthers(column2_exist(2)));
if ~isempty(column2_missing)
    seqdata(:,12) = zeros(blobs,1);
else
    seqdata(:,12) = seq.data(:,colOthers(column2_exist(3)));
end
end
            
        
  
  

