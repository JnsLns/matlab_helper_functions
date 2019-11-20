% foundInds = strfindn(str,pat,num,backwards)
%
% Searches string str for integer num occurrences of string pat. Elements
% of the string str already counted as part of one occurrence will not be
% considered as part of another occurrence (e.g., 'ppp' contains pattern 'pp'
% only once, not twice). If backwards ==  1, str is searched backwards,
% else (backwards == 0, default) search starts from the first element.
%

% Jonas Lins, Jan 26 2016

function foundInds = strfindn(str,pat,num,backwards)

strlen = length(str);
patlen = length(pat);

if strlen < patlen
    error('String str must have at least as many elements as string pat');
end

foundInds = [];
foundNo = 0;

if nargin < 4 || backwards ~= 1
    backwards = 0;
end

i = 1;
% loop over string elements
while i <= strlen-patlen+1 && foundNo < num
    
    if backwards % from last to first element
        currentInd = strlen-patlen+2-i;
    else
        currentInd = i;
    end
    
    % compare
    if str(currentInd:currentInd+patlen-1) == pat
        
        % store starting index of found pattern
        foundInds(end+1) = currentInd;
        foundNo = foundNo + 1;
        % make sure once found, parts of a pattern are not counted again
        i = i + patlen;
        
    else
        
        i = i + 1;
        
    end
    
end


end