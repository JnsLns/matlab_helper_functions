% function meanArray = cellMean(cellArr,meanDim,pageIndices,catFlag)
%
% Computes cell-by-cell means over the contents of two or more pages of a 
% cell array of arbitrary dimensionality. The contents of a component in
% one page are simply added to that of the corresponding cell (same subscript
% index except for dimension meanDim) in the other page(s) and the overall
% sum is divided by the number of values thus included in the sum.
% Note that empty array cells are not included in the mean, that is, they
% do not count against the number of values by which the sum is divided
% (and obviously they do not contribute to the sum). Cells holding matrices
% that contain only NaNs are treated in the same way as empty cells.
% The contents of a cell in the output array are determined in the following
% way, depending on the contents of the averaged set of cells from the
% input array:
%
% - Only empty cells -> output empty cell
% - Only NaNs -> output NaN
% - Only empty cells mixed with NaNs -> output NaN
% - NaNs and/or empty cells mixed with at least one double matrix, vector,
%   or scalar -> output is mean over all doubles in the set (sum divided
%   by number of doubles).
% - Only doubles -> mean over these 
% 
% Page numbers are passed by pageIndices as a row vector. The array dimension
% along which averaging should take place is given by meanDim. The output
% cell array containing the results of the averaging has the same dimensio-
% nality as the input array and the same size, except for dimension meanDim,
% along which the size is one (if catFlag == 0) or size(cellArr,meanDim)+1
% (if catFlag == 1). Note: The nature of the contents of involved cells is
% irrelevant as long as they provide legal operands to addition and division
% operations.
%
% IMPORTANT NOTE: Empty cells included in an average are ignored (see above)
% if there is at least one non-empty cell in the average at hand. An average
% over only empty cells returns an empty cell as well.
%
% pageIndices is optional. If not specified or set to [0], all pages are
% averaged.
%
% catFlag (default 0) is an optional argument. If 0, function output is
% simply a cell array page of means, if 1, function output is the page
% of means concatenated to the input cell array as the last page
% along meanDim.



function meanArray = cellMean(cellArr,meanDim,pageIndices,catFlag)

cellArr_backup = cellArr;

if nargin < 4
    catFlag = 0;
end   

if nargin < 3 || all(pageIndices == 0)
    pageIndices = 1:size(cellArr,meanDim);
end

if numel(pageIndices)<2
    error('You cannot take a mean over less than two values/array pages!')
end

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArr)
    ssndcs(lndx,:) = ind2subAll(size(cellArr),lndx);
end

% for each page given in pageIndices, find linear indices where meanDim is
% equal to that page. (in resulting array, rows are pages and cols are
% elements of that page)
for curPage = 1:numel(pageIndices)   
    curPageNdx = pageIndices(curPage);
    page_lndcs(curPage,:) = find(ssndcs(:,meanDim) == curPageNdx);
end

% Check for each colum of cells along meanDim whether it includes empty
% cells, if so, fill the empty cells with a zero matrix of the same size as the
% matrix in the cell.
% Store for each column along meanDim the number of non-empty cells in order
% to later divide by that number instead of dividing by the size of the
% cell array (in other words: although a non-zero matrix is entered there
% to allow summing, that cell is not included in the mean in terms of not
% counting against the number of values by which the sum is divided)
% Also check whether all matrices along the column are of the same size
% (i.e., whether they are legal operands to the addition operation); error
% if not.
%
% Make lookup array that has same size and dimensionality as cellArr,
% except for meanDim, along which it has size one; fill each cell with its
% subscript index
warningEmptiesIssued = 0;
sz = size(cellArr);
sz(meanDim) = 1;
subsArr = cell(sz);
nNonEmptyCells = zeros(sz);
for curLndx = 1:numel(subsArr)
    subsArr{curLndx} = ind2subAll(size(subsArr),curLndx);
end
% Go through subsArr's cells and use subscripts as indices into columns along meanDim
for curLndx = 1:numel(subsArr)    
    % get subscript index of current cell
    curSubs = subsArr{curLndx};        
    % first, check largest numel of any array in cells along current column of meanDim
    % (that array's size is used as size against which that of other arrays is compared below)
    fro = curSubs; fro(meanDim) = 1; to = curSubs; to(meanDim) = size(cellArr,meanDim);
    curEntireColumn = subArray(cellArr,fro,to);
    curSubArr = curEntireColumn(pageIndices); % consider only pages included in sum
    numelsInCurSubArray = cellfun(@(curEl) numel(curEl),curSubArr);    
    [~,maxInd] = max(numelsInCurSubArray); % largest number of elements 
    baseSize = size(curSubArr{maxInd});     
    nNonEmptyCells(curLndx) = sum(numelsInCurSubArray>0);
    % Check for cells in current subarray that contain matrices which only
    % contain nans (these are treated exactly like empty cells)    
    onlyNanPagesLogndx_curSubArray = cell2mat(cellfun(@(curEl) all(isnan(curEl(:)))&~isempty(curEl),squeeze(curSubArr),'uniformoutput',0));            
    onlyNanPagesLogndx_curEntireColumn = cell2mat(cellfun(@(curEl) all(isnan(curEl(:)))&~isempty(curEl),squeeze(curEntireColumn),'uniformoutput',0));        
    nNonEmptyCells(curLndx) = nNonEmptyCells(curLndx)-sum(onlyNanPagesLogndx_curSubArray);
    % go through cellArr's cells having that subscript index along meanDim
    for curPgNdx = 1:numel(pageIndices)        
        curPage = pageIndices(curPgNdx); % Do only pages included in sum        
        curSubs(meanDim) = curPage;
        curCellLindex = sub2indAll(size(cellArr),curSubs);
        if ~isempty(cellArr{curCellLindex}) && ~onlyNanPagesLogndx_curEntireColumn(curPage)                     
            % compare size of contents of current cell to baseline size
            if ~isequal(size(cellArr{curCellLindex}), baseSize)
                error('You cannot take means over matrices of differing sizes.')
            end            
        else
           % replace empty cells or those which contain matrices with only nans
           % with zero matrix of same size as baseSize (to enable summing in code below)
           cellArr{curCellLindex} = zeros(baseSize);
           if ~warningEmptiesIssued 
                warning('cellMeanFun:EmptyCells',['Found empty cell or cell with NaN-matrix during usage of cellMean; first found is at ', num2str(curSubs), '; such cells are ignored, that is, they are not included in the mean in any way.']);
                warningEmptiesIssued = 1;
           end
        end                
    end    
end

% Do the actual summing
for curPage = 1:2:size(page_lndcs,1)   
    if (curPage+1) <= size(page_lndcs,1)                                        
        curSum = cellfun(@(p1,p2) p1+p2,cellArr(page_lndcs(curPage,:)), cellArr(page_lndcs(curPage+1,:)),'uniformoutput',0);
        if curPage == 1
            sumOfPages = curSum;
        else
            sumOfPages = cellfun(@(p1,p2) p1+p2, curSum, sumOfPages, 'uniformoutput', 0);
        end
    elseif curPage == size(page_lndcs,1)
       sumOfPages = cellfun(@(p1,p2) p1+p2, cellArr(page_lndcs(curPage,:)), sumOfPages, 'uniformoutput', 0);             
    end
end

% Reshape into the shape of a page of meanDim
szPage = size(cellArr);
szPage(meanDim) = 1;
sumOfPages = reshape(sumOfPages,szPage);

% convert sum to mean (divide each sum by the number of non-empty cells for the column that this
% sum is based on)
meanOfPages = cellfun(@(sop,nec) sop/nec,sumOfPages,num2cell(nNonEmptyCells),'uniformoutput',0);

if catFlag == 1
    % page with meanOfPages is concatenated with original cell array, thus
    % increasing its size by one along meanDim.
    meanArray = cat(meanDim,cellArr_backup,meanOfPages);
else
    meanArray = meanOfPages;    
end




 