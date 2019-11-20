% function sumArray = cellSum(cellArr,sumDim,pageIndices,catFlag)
%
% Computes cell-by-cell sums between the contents of two or more pages of a 
% cell array of arbitrary dimensionality, analogous to matrix operations
% like a(:,3,:) + a(:,5,:). The contents of a component in one page 
% are simply added to that of the corresponding cell (same subscript
% index except for dimension sumDim) in the other page. Page numbers are
% passed by pageIndices as a row vector. The array dimension along which
% addition should take place is given by sumDim. The output cell array
% containing the results of the addition has the same dimensionality as
% the input array and the same size, except for dimension sumDim, along
% which the size is one (if catFlag == 0) or size(cellArr,sumDim)+1 (if
% catFlag == 1). Note: The nature of the contents of involved cells is
% irrelevant as long as they provide legal operands to the addition operation.
%
% IMPORTANT NOTE: Empty cells in one or more of the summed cells are treated
% as zero matrices of the same size as the other matrices in the sum (i.e.,
% they are basically ignored). A sum of only empty cells, however, returns
% an empty cell as well.
%
% pageIndices is optional. If not specified or set to [0], all pages are summed.
%
% catFlag (default 0) is an optional argument. If 0, function output is
% simply a cell array page of sums, if 1, function output is the page
% of sums concatenated to the input cell array as the last page
% along sumDim.



function sumArray = cellSum(cellArr,sumDim,pageIndices,catFlag)

cellArr_backup = cellArr;

if nargin < 4
    catFlag = 0;
end   

if nargin < 3 || all(pageIndices == 0)
    pageIndices = 1:size(cellArr,sumDim);
end

if numel(pageIndices)<2
    error('You cannot sum less than two values/array pages!')
end

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArr)
    ssndcs(lndx,:) = ind2subAll(size(cellArr),lndx);
end

% for each page given in pageIndices, find linear indices where sumDim is
% equal to that page. (in resulting array, rows are pages and cols are
% elements of that page)
for curPage = 1:numel(pageIndices)   
    curPageNdx = pageIndices(curPage);
    page_lndcs(curPage,:) = find(ssndcs(:,sumDim) == curPageNdx);
end

% Check for each colum of cells along sumDim whether it includes empty
% cells, if so, fill the empty cells with a zero matrix of the same size as the
% matrix in the cell.
% Also check whether all matrices along the column are of the same size
% (i.e., whether they are legal operands to the addition operation); error
% if not.
%
% Make lookup array that has same size and dimensionality as cellArr,
% except for sumDim, along which it has size one; fill each cell with its
% subscript index
warningIssued = 0;
sz = size(cellArr);
sz(sumDim) = 1;
subsArr = cell(sz);
for curLndx = 1:numel(subsArr)
    subsArr{curLndx} = ind2subAll(size(subsArr),curLndx);
end
% Go through subsArr's cells and use subscripts as indices into columns along sumDim
for curLndx = 1:numel(subsArr)    
    % get subscript index of current cell
    curSubs = subsArr{curLndx};        
    % first, check largest numel of any array in cells along current column of sumDim
    % (that array's size is used as size against which that of other arrays is compared below)
    fro = curSubs; fro(sumDim) = 1; to = curSubs; to(sumDim) = size(cellArr,sumDim);
    curSubArr = subArray(cellArr,fro,to);
    curSubArr = curSubArr(pageIndices); % consider only pages included in sum
    [~,maxInd] = max(cellfun(@(column) numel(column),curSubArr));
    baseSize = size(curSubArr{maxInd});     
    % go through cellArr's cells having that subscript index along sumDim
    for curPgNdx = 1:numel(pageIndices)        
        curPage = pageIndices(curPgNdx); % Do only pages included in sum        
        curSubs(sumDim) = curPage;
        curCellLindex = sub2indAll(size(cellArr),curSubs);
        if ~isempty(cellArr{curCellLindex})                        
            % compare size of contents of current cell to baseline size
            if ~isequal(size(cellArr{curCellLindex}), baseSize)
                error('You cannot sum matrices of differing sizes.')
            end
        else
           % replace empty cells with zero matrix same size as baseSize
           % (to enable summing in code below)
           cellArr{curCellLindex} = zeros(baseSize);
           if ~warningIssued 
                warning(['Found empty cell for a summand during usage of cellSum (first found is at ', num2str(curSubs), '; treating as zero).']);
                warningIssued = 1;
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

% Reshape into the shape of a page of sumDim
szPage = size(cellArr);
szPage(sumDim) = 1;
sumOfPages = reshape(sumOfPages,szPage);

if catFlag == 1
    % page with sumOfPages is concatenated with original cell array, thus
    % increasing its size by one along sumDim.
    sumArray = cat(sumDim,cellArr_backup,sumOfPages);
else
    sumArray = sumOfPages;    
end




 