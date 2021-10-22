function S = nestedstruct(varargin)
% S = NESTEDSTRUCT(FIELD,VAL,FIELD2,VAL2,...) - Create a nested structure with nested-field-names 
%   and values provided as argument pairs FIELD,VAL,FIELD2,VAL2,...
%   e.g. S = nestedstruct('a',1,'b.c',2,'b.d',3) -> S.a = 1, S.b.c = 2, S.b.d = 3
%
% See also: GETNESTEDFIELD, SETNESTEDFIELD, CELL2NESTEDSTRUCT

    assert(mod(nargin,2) == 0,'Arguments must come in pairs');
    varargin = reshape(varargin,2,[]);
    S = cell2nestedstruct(varargin(2,:),varargin(1,:));
end

