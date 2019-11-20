function sub = ind2subAll(siz, ndx)
% function sub = ind2subAll(siz,ndx)
%
% Wraps MATLAB function ind2sub (which converts a linear index into a set
% of subscript indices based on the target array's size, 'siz') but
% returns the subscript indices as a row vector instead of requiring
% multiple output arguments to be specified.
%
% See also IND2SUB

str = '';
for i = 1:length(siz)
str = strcat(str, ['sub(', num2str(i), '), ']);
end
str(end) = [];

eval(['[',str,'] = ind2sub(siz, ndx);']);

end
