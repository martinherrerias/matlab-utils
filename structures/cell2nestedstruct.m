function S = cell2nestedstruct(values,fieldnames)
% S = CELL2NESTEDSTRUCT(VALUES, FIELDNAMES) - create a nested structure with nested-field-names 
%   given in array-cell 'fieldnames', and values given in nested array-cell 'values'.
%   e.g. S = CELL2NESTEDSTRUCT({1,2,3},{'a','b.c','b.d'}) -> S.a = 1, S.b.c = 2, S.b.d = 3
%
% See also: GETNESTEDFIELD, SETNESTEDFIELD, NESTEDSTRUCT

    assert(iscell(fieldnames) && iscell(values) && isequal(size(fieldnames),size(values)),...
        'nestedstruct:args','Expecting equal-sized cell-arrays');

    S = struct();
    for j = 1:numel(fieldnames)
        f = strsplit(fieldnames{j},'.');
        S = setfield(S,{1},f{:},values{j});
    end
end

