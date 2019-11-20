% function orArray = cellLogicalOr(cellArr,alongDim,pageIndices,catFlag)
%
% Computes cell-by-cell logical OR between the contents of two or more pages
% of a cell array of arbitrary dimensionality, analogous to matrix operations
% like a(:,3,:) | a(:,5,:). The contents of a component in one page 
% are simply logically combined with that of the corresponding cell (same
% subscript index except for dimension alongDim) in the other page. Page
% numbers are passed by pageIndices as a row vector. The array dimension
% along which the logic operation should take place is given by alongDim. The
% output cell array containing the results of the logic operation has the same
% dimensionality as the input array and the same size, except for dimension
% alongDim, along which the size is 1 (if catFlag == 0) or
% size(cellArr,alongDim)+1 (if catFlag == 1). Note: The type (matrix,
% scalar,...) of the contents of involved cells is irrelevant as long as
% they provide legal operands to the logical operation.
%
% catFlag (default 0) is an optional argument. If 0, function output is
% simply a cell array page of logicals, if 1 function output is the page
% of logicals concatenated to the input cell array as the last page
% along alongDim.

function orArray = cellLogicalOr(cellArr,alongDim,pageIndices,catFlag)

if nargin < 4
    catFlag = 0;
end   

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArr)
    ssndcs(lndx,:) = ind2subAll(size(cellArr),lndx);
end

% for each page give in pageIndices, find linear indices where alongDim is
% equal to that page. (in resulting array, rows are pages and cols are
% elements of that page)
for curPage = 1:numel(pageIndices)   
    curPageNdx = pageIndices(curPage);
    page_lndcs(curPage,:) = find(ssndcs(:,alongDim) == curPageNdx);
end

% Compute sum of the pages
for curPage = 1:2:size(page_lndcs,1)
    if (curPage+1) <= size(page_lndcs,1)
        curSum = cellfun(@(p1,p2) p1 | p2,cellArr(page_lndcs(curPage,:)), cellArr(page_lndcs(curPage+1,:)),'uniformoutput',0);
        if curPage == 1
            sumOfPages = curSum;
        else
            sumOfPages = cellfun(@(p1,p2) p1 | p2, curSum, sumOfPages, 'uniformoutput', 0);
        end
    elseif curPage == size(page_lndcs,1)
        sumOfPages = cellfun(@(p1,p2) p1 | p2, cellArr(page_lndcs(curPage,:)), sumOfPages, 'uniformoutput', 0);
    end
    
end


% Reshape into the shape of a page of alongDim
szPage = size(cellArr);
szPage(alongDim) = 1;
sumOfPages = reshape(sumOfPages,szPage);

if catFlag == 1
    % page with sumOfPages is concatenated with original cell array, thus
    % increasing its size by one along alongDim.
    orArray = cat(alongDim,cellArr,sumOfPages);
else
    orArray = sumOfPages;    
end




 