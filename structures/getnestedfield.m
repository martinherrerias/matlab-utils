function v = getnestedfield(S,field)
% v = GETNESTEDFIELD(S,field) - Return the value of S.(field), where field is a [nested] field name
%   of the form 'A.B.C.D..'
%
% See also: NESTEDSTRUCT, SETNESTEDFIELD, NESTEDFIELDNAMES

    f = strsplit(field,'.');
    v = getfield(S,{1},f{:});
end
