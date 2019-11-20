function ndx = sub2indAll(siz, ssInput)
% function ndx = sub2indAll(siz,ssInput)
%
% Wraps MATLAB function sub2ind (which converts set of subscript indices to
% a single linear index based on the target array's size, 'siz') but
% allows passing the subscript input as a row vector 'ssInput' instead of
% requiring to specify multiple input arguments.
%
% See also SUB2IND

ssInput = sprintf('%.0f,' , ssInput);
ssInput = ssInput(1:end-1);
eval(['ndx = sub2ind(siz, ', ssInput, ');']);

end