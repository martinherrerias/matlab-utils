function isnf = isnestedfield(S,f)
% ISNESTEDFIELD(S,f) - Generalized version of ISFIELD, that takes objects as well as structures,
%   and nested field names. True if S has a [nested] field/prop f, e.g. S.a.b.c for f = 'a.b.c'
%
% See also: GETNESTEDFIELD, SETNESTEDFIELD, NESTEDFIELDNAMES, ISFIELD

    if iscell(f), isnf = cellfun(@(x) isnestedfield(S,x),f); return; end
    
    assert((isstruct(S) || isobject(S)) && ischar(f),'Invalid arguments');
    if isempty(f), isnf = false(size(f)); return; end
    
    ff = strsplit(f,'.');
    if isfield(S,ff{1}) || (~isstruct(S) && isprop(S,ff{1}))
        if numel(ff) > 1
            if isstruct(S.(ff{1})) || isobject(S.(ff{1}))
                isnf = isnestedfield(S.(ff{1}),strjoin(ff(2:end),'.'));
            else
                isnf = false;
            end
        else
            isnf = true;
        end
    else
        isnf = false;
    end
end
