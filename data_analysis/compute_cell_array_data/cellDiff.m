% function diffArray = cellDiff(cellArr,diffDim,minuendPage,subtrahendPage,catFlag)
%
% Computes cell-by-cell differences between the contents of two pages of a 
% cell array of arbitrary dimensionality, analogous to matrix operations
% like a(:,3,:) - a(:,5,:). The content of a component in the subtrahend page 
% are simply subtracted from that of the corresponding cell (same subscript
% index except for dimension diffDim) in the minuend page. Page numbers are
% passed by minuendPage and subtrahendPage. The array dimension along which
% subtraction should take place is given by diffDim. The output cell array
% containing the results of the subtractions has the same dimensionality as
% the input array and the same size, except for dimension diffDim, along
% which the size is one (if catFlag == 0) or size(cellArr,diffDim)+1 (if
% catFlag == 1). Note: The nature of the contents of involved cells is
% irrelevant as long as they provide legal minuends and subtrahends.
%
% catFlag (default 0) is an optional argument. If 0, function output is
% simply a cell array page of differences, if 1 function output is the page
% of differences concatenated to the input cell array as the last page
% along diffDim.

function diffArray = cellDiff(cellArr,diffDim,minuendPage,subtrahendPage,catFlag)

if nargin < 5
    catFlag = 0;
end   

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArr)
    ssndcs(lndx,:) = ind2subAll(size(cellArr),lndx);
end

% find all linear indices of cellArr where diffDim is equal to the minuend page
minuendSet_lndcs = find(ssndcs(:,diffDim) == minuendPage);
% find all linear indices of cellArr where diffDim is equal to the subtrahend page
subtrahendSet_lndcs = find(ssndcs(:,diffDim) == subtrahendPage);


% Compute difference
difference = cellfun(@(m,s) m-s,cellArr(minuendSet_lndcs), cellArr(subtrahendSet_lndcs),'uniformoutput',0);

% Reshape into the shape of a page of diffDim
szPage = size(cellArr);
szPage(diffDim) = 1;
difference = reshape(difference,szPage);

if catFlag == 1
    % page with differences is concatenated with original cell array, thus
    % increasing its size by one along diffDim.
    diffArray = cat(diffDim,cellArr,difference);
else
    diffArray = difference;    
end




 