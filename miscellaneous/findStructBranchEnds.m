function branchEndPaths = findStructBranchEnds(inputStruct)

% !!! This works in principle but may not work in some cases since I did
% this quickly and more or less dirty... should be tested more thoroughly
% before using in new code !!!

% Returns a one-column cell array containing strings that are the "paths"
% to the "branch ends" of a nested structure (i.e., to those fields that do
% not contain further structs themselves, but other data types). 
%
% Note that this function can be called only with a plain structname as
% argument; expressions are not allowed (due to usage of inputname(), which
% cannot deal with expressions).

% Give input data the same variable name as used during the function call
eval([inputname(1) '= inputStruct;']);
addresses = {inputname(1)};
curLayer = 1;
branchEndPaths = {};
while curLayer > 0
    
    for curAddressNum = 1:numel(addresses(curLayer,:))
        
        % get starting address for this loop
        curAddress = addresses(curLayer,curAddressNum);
        % delete current address so it's not done again
        addresses{curLayer,curAddressNum} = [];
        
        if isempty(curAddress{:})
            if curAddressNum == numel(addresses(curLayer,:))
                curLayer = curLayer - 1; % back one layer
                if curLayer == 0 % end criterion
                    break;
                end
            end
            continue;
        end
        
        % in case current field contains another struct, store the fields
        % address to later examine it; else, store its address in
        % dataAdresses
        if isstruct(eval(cell2mat(curAddress)))
            curLayerFields = fieldnames(eval(cell2mat(curAddress)))';
            for curFieldNum = 1:numel(curLayerFields)
                curFieldName = curLayerFields(curFieldNum);
                addresses{curLayer+1,curFieldNum} = cell2mat([cell2mat(curAddress),'.',curFieldName]);
            end
        else
            branchEndPaths{end+1,1} = cell2mat(curAddress);
            continue;
        end
        
        % go down one layer and restart the for loop
        curLayer = curLayer + 1;
        break;
        
    end
    
end




end



