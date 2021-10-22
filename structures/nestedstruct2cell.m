function [values,names,implicit] = nestedstruct2cell(s)
% [VALUES,NAMES] = NESTEDSTRUCT2CELL(S) - Get a list of filed values [and names] for a structure, 
%   including the field values and names of all nested structures.

    [names,implicit] = nestedfieldnames(s);
    values = cell(size(names));
    n = numel(names);
    for j = 1:n
        f = strsplit(names{j},'.');
        values{j} = getfield(s,{1},f{:});
    end
end
