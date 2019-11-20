function stdArray = cellStdOfDiff(cellArrStds,cellArrNumCases,alongDim,pageA,pageB,catFlag)
%
%
% Similar as cellDiffs() but computes the pooled standard deviation for a
% difference of means from the two known separate standard deviations s1 and
% s2, according to the formula:
% 
%    s_pooled = sqrt( ((s1.^2)*(n1-1) + (s2.^2)*(n2-1)) / ((n1-1)+(n2-1)) )
%
% This formula effectively reverts sds to sum of squares and, by dividing
% by the overall number of cases, arrives at pooled variance, while the 
% square root leads to the pooled standard deviation. Note that the standard 
% error in the denominator of a two-sample t-test can be computed from this 
% value by:
%
%    se_pooled = s_pooled * sqrt(1/n1 + 1/n2)    

if nargin < 6
    catFlag = 0;
end   

if isa(cellArrNumCases,'double')
    cellArrNumCases = num2cell(cellArrNumCases);
end

% make "lookup table" for conversion of linear indices to subscript indices
% in the current array (in this table, rows are linear indices, columns are
% dimensions of the input array)
for lndx = 1:numel(cellArrStds)
    ssndcs(lndx,:) = ind2subAll(size(cellArrStds),lndx);
end

% find all linear indices of cellArrStds where alongDim is equal to page a
pageASet_lndcs = find(ssndcs(:,alongDim) == pageA);
% find all linear indices of cellArrStds where alongDim is equal to page b
pageBSet_lndcs = find(ssndcs(:,alongDim) == pageB);

% Compute difference
sdev = cellfun(@(s1,s2,n1,n2) sqrt( ((s1.^2)*(n1-1) + (s2.^2)*(n2-1)) / ((n1-1)+(n2-1)) ), ...
    cellArrStds(pageASet_lndcs), cellArrStds(pageBSet_lndcs),... 
    cellArrNumCases(pageASet_lndcs), cellArrNumCases(pageBSet_lndcs), ...
    'uniformoutput',0);

% Reshape into the shape of a page of alongDim
szPage = size(cellArrStds);
szPage(alongDim) = 1;
sdev = reshape(sdev,szPage);

if catFlag == 1
    % page with sdevs is concatenated with original cell array, thus
    % increasing its size by one along alongDim.
    stdArray = cat(alongDim,cellArrStds,sdev);
else
    stdArray = sdev;    
end




 