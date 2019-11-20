% function outArr = cellPageDelete(arr,dim,pageInd)
%
% This function allows to remove entire pages from a multidimensional 
% array without knowing the array's dimensionality (ndims). (This is
% normally achieved by the subscript index-based syntax, making it difficult
% to work dynamically with arrays whose dimensionality isn't known at
% coding time).
%
% arr: input array.
%
% dim: dimension along which page(s) should be removed.
%
% pageInd: page(s) along dim that should be deleted. Scalar for single
% page, row vector for multiple pages.
%
%
% Examples:
%   
%   For ndims(arr) = 3, arr = cellPageDelete(arr,2,3) is equivalent to
%   arr(:,3,:) = []
%
%   For ndims(arr) = 5, arr = cellPageDelete(arr,4,[3 8]) is equivalent to
%   arr(:,:,:,[3 8],:) = [];


function outArr = cellPageDelete(arr,dim,pageInd)

ssStr = '';
for i = 1:dim-1
    ssStr = [ssStr, ' :,'];
end
ssStr = [ssStr, '[', num2str(pageInd(:)') ']' ];

if dim ~= numel(size(arr))
    for i = dim+1:numel(size(arr))
        ssStr = [ssStr, ', :'];
    end
end

eval(['arr(', ssStr, ') = [];']);

outArr = arr;