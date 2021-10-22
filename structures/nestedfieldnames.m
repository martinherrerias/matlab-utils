function [names,implicit] = nestedfieldnames(s)
% [NAMES,IMPLICIT] = NESTEDFIELDNAMES(S) - Get a list of field names for a structure, including   
%   the subfields of all nested structures. 
%
%   NAMES - cell array of field-names of the form 'A.B.C,..' for any (sub) fields that are NOT
%       structures (i.e. have no subfields of their own). Useful for reconstruction of the
%       hierarchical nested structure using e.g. CELL2NESTEDSTRUCT.
%   IMPLICIT - complement of NAMES: cell array of field-names for structure (sub) fields that are
%       themselves structures. Useful for comparing hierarchies, e.g. COMPLETESTRUCT. 
%
% See also: CELL2NESTEDSTRUCT, GETNESTEDFIELD, SETNESTEDFIELD

    names = fieldnames(s);
    subfields = {};
    implicit = {};
    n = numel(names);
    havesubfields = false(n,1);
    for i = 1:numel(s)
        for j = 1:n
            if isstruct(s(i).(names{j}))
                havesubfields(j) = true;
                [sf,si] = nestedfieldnames(s(i).(names{j}));
                sf = cellfun(@(s) [names{j} '.' s],sf,'unif',0);
                si = cellfun(@(s) [names{j} '.' s],si,'unif',0);
                subfields = [subfields(:); sf(:)];
                implicit = [implicit(:); si(:)];
            end
        end
    end
    
    implicit = [names(havesubfields); implicit(:)];
    names = [names(~havesubfields); subfields(:)];
    if i > 1
        names = unique(names,'stable'); 
        implicit = unique(implicit,'stable');
    end
end
