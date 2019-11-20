function collapsedArray = cellCollapse(cellArr,collapseDim,catDim,catFlag)
%
% "Collapses" a cell array along one of its dimension (collapseDim). Any
% data in the cells along that dimension is concatenated with
% cat(catDim, data_1, ..., data_n) and stored in the first page of cells
% of collapseDim. If catFlag==0 (default, argument is optional), all other
% pages are removed, so that the returned array has size 1 along the di-
% mension collapseDim. Other than that, dimension sizes and overall
% dimensionality as well as the order of dimensions is retained.
% If catFlag==1, then no pages are removed and the new array page holding the
% concatenated arrays is concatenated to the original cell array as last page
% along collapseDim.
%
% Function works for both matrices and cell arrays as data within the
% input array. Note that the contents of cellArr must be valid inputs
% to cat(), that is, they must have the same size (except for
% size(data,catDim) itself). If not specified, catDim is by default set
% to the same value as collapseDim.
%
%  Examples:
%
%  a = {[1 2],[3 4],[5 6];[7 8],[9 10],[11 12]}
%  a =
%       [1 2] [3 4]  [5 6]
%       [7 8] [9 10] [11 12]
%
%  b = cellCollapse(a,2,2)
%  b =
%       [1 2 3 4 5 6]
%       [7 8 9 10 11 12]
%
%  c = cellCollapse(a,2,1)
%  c =
%       [1 2
%        3 4
%        5 6]
%
%       [7 8
%        9 10
%       11 12]

if nargin < 4
   catFlag = 0; 
end

if nargin < 3
   catDim = collapseDim; 
end

if catFlag
    cellArr_backup = cellArr;
end

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArr)
    ssndcs(lndx,:) = ind2subAll(size(cellArr),lndx);
end

% find all linear indices of cellArr where collapseDim is 1 (array will be
% collapsed onto that subarray)
collapseOntoSet_lndcs = find(ssndcs(:,collapseDim) == 1);

% loop over components of subarray onto which collapsing is done
for curTgt_num = 1:numel(collapseOntoSet_lndcs)
    
    % loop over components along collapse dim and collapse onto first one
    % (which is part of the subarray from above)
    for curSource_num = 2:size(cellArr,collapseDim)
        
        curTgt_lndx = collapseOntoSet_lndcs(curTgt_num);
        curTgt_ssndx = ssndcs(curTgt_lndx,:);
        
        curSource_ssndx = curTgt_ssndx;
        curSource_ssndx(collapseDim) = curSource_num;
        
        curSource_lndx = find(ismember(ssndcs,curSource_ssndx,'rows'),1);
                
        cellArr{curTgt_lndx} = cat(catDim, cellArr{curTgt_lndx}, cellArr{curSource_lndx});              
                        
    end
    
end

% Remove all cells but those that hold the collapsed data (dimensionality
% is retained).
ssStr = '';
for i = 1:collapseDim-1
    ssStr = [ssStr, ' :,'];
end
ssStr = [ssStr, ' 2:end'];
if collapseDim ~= numel(size(cellArr))
    for i = collapseDim+1:numel(size(cellArr))
        ssStr = [ssStr, ', :'];
    end
end
eval(['cellArr(', ssStr, ') = [];']);

if catFlag
    collapsedArray = cat(collapseDim,cellArr_backup,cellArr);
else
    collapsedArray = cellArr;
end
