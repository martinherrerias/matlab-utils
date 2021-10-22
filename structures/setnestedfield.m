function S = setnestedfield(S,field,val)
% S = SETNESTEDFIELD(S,field,val) - S the value of S.(field) = val, where field is a [nested] 
%   field name of the form 'A.B.C.D..'
%
% See also: NESTEDSTRUCT, GETNESTEDFIELD, NESTEDFIELDNAMES

    assert((isstruct(S) || isobject(S)) && ischar(field),'Invalid arguments');
    if ~isscalar(S)
       warning('setnestedfield:nonscalar','SETNESTEDFIELD should be used with scalar structures only!'); 
    end

    if nargout < 1 && ~isa(S,'handle')
        warning('setnestedfield:waste','Use S = SETNESTEDFIELD(S,..) for non-handle classes')
    end
    
    f = strsplit(field,'.');

    for j = 1:numel(f)-1
        if j == 1, s = S; end
        if ~isfield(s,f{j}), break; end
        s = s.(f{j});
        if ~(isstruct(s) || isobject(s))
            warning('setnestedfield:nostruct',...
                'Replacing non-structure field %s to add %s',strjoin(f(1:j),'.'),field);
            S = setfield(S,{1},f{1:j},struct());
        end
    end
    S = setfield(S,{1},f{:},val);
end
