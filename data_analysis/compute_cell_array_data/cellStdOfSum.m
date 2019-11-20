function stdArray = cellStdOfSum(cellArrStds,cellArrNumCases,alongDim,pageA,pageB,catFlag)
%
%
% Similar as cellDiff() and cellSum() but computes the pooled standard
% deviation for a sum/difference of means from the two known separate
% standard deviations s1 and s2 (where the two samples must be statistically
% independent).
%
% cellArrStds is the cell array of input standard deviations
%
% cellArrNumCases can EITHER be a cell array or matrix that has same size as
% cellArrStds and holds the case numbers based on which the standard
% deviation in the latter were computed. In that case, the pooled standard
% deviation is computed using the formula (*)
% 
%    s_pooled = sqrt( ((s1.^2)*(n1-1) + (s2.^2)*(n2-1)) / ((n1-1)+(n2-1)) )
%
% which means that the variances are weighted by case numbers (actually,
% by degress of freedom) before summing and square rooting (i.e., the result 
% will be the same as if computing the "grand" standard deviation across all
% original values at once). OR it can be [0], in which case the used formula is  
%
%    s_pooled = sqrt((s1.^2 + s2.^2)/2)
% 
% which is the more common way of computing the standard deviation of a sum
% of random variables and is equivalent to the above formula when both 
% input standard deviations are based in equal case numbers.
%
% alongDim gives the dimension of the input array along which pages will be
% joined.
%
% pageA and pageB give the array pages whose standard deviation will be  
% combined.
%
% catFlag determines whether (1) or not (0) the resulting page of standard 
% deviations is concatenated with the original cell array for output (thus
% increasing its size by one along alongDim). If not, output is just the
% page with the combined standard deviations.
%
%
% (*) Note that the standard  error in the denominator of a two-sample
% t-test can be computed from this  value by:
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
if numel(cellArrNumCases)==1 && cellArrNumCases{1}==0
    sdev = cellfun(@(s1,s2) sqrt((s1.^2+s2.^2)/2), ...
        cellArrStds(pageASet_lndcs), cellArrStds(pageBSet_lndcs),...
        'uniformoutput',0);    
else    
    sdev = cellfun(@(s1,s2,n1,n2) sqrt( ((s1.^2)*(n1-1) + (s2.^2)*(n2-1)) / ((n1-1)+(n2-1)) ), ...
        cellArrStds(pageASet_lndcs), cellArrStds(pageBSet_lndcs),...
        cellArrNumCases(pageASet_lndcs), cellArrNumCases(pageBSet_lndcs), ...
        'uniformoutput',0);
end

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




 