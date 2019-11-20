function SEArray = cellSEOfSum(cellArrStds,cellArrNumCases,alongDim,pageA,pageB,catFlag)
%
%
% Similar as cellDiffs() but computes the standard deviation of the sampling 
% distribution of the sum or difference of means, that is, the standard error,
% according to the formula:
% 
% SE = sqrt( s1^2 / n1 + s2^2 / n2 )
%
% where s1 and s2 are the standard deviations of the two means and n1 and
% n2 are the number of values upon which the means are based.

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
sdev = cellfun(@(a,b,ancs,bncs) sqrt((a.^2)./ancs + (b.^2)./bncs), ...
    cellArrStds(pageASet_lndcs), cellArrStds(pageBSet_lndcs),...
    cellArrNumCases(pageASet_lndcs), cellArrNumCases(pageBSet_lndcs),...
    'uniformoutput',0);
% Reshape into the shape of a page of alongDim
szPage = size(cellArrStds);
szPage(alongDim) = 1;
sdev = reshape(sdev,szPage);

if catFlag == 1
    % page with sdevs is concatenated with original cell array, thus
    % increasing its size by one along alongDim.
    SEArray = cat(alongDim,cellArrStds,sdev);
else
    SEArray = sdev;    
end




 