function varargout = compatiblesize(varargin)
% COMPATIBLESIZE(A,B,C,..) - Check if arrays A,B,.. are (singleton-expansion) size-compatible,
%   i.e. that they can work together under element-wise operations, e.g. A + B.*C..
%   Throw an error (as caller) otherwise.
%
%   More precisely, check that for all X in {A,B,C,..}, and along any dimension j:
%       size(X,j) == S(j) || size(X,j) == 1, where S(j) is the maximum size along dimension j.
%
%   In agreement with allowed expressions like [].*1 and repmat(1,0,..), the combination of 
%   scalars and empty dimensions is also allowed, i.e. size(X,j) <= 1 for all X.
%
% [AA,BB,CC,..] = COMPATIBLESIZE(A,B,C,..) - If everything goes right, return S-sized expanded 
%   arrays AA,BB,CC,..
%
% S = COMPATIBLESIZE(A,B,C,..'-size') - Return (common) expanded size S as a single output.

    narginchk(2,Inf);
    nargoutchk(0,nargin);
    
    outputsize = nargout <= 1 && ischar(varargin{end}) && strcmpi(varargin{end},'-size');
    if outputsize, varargin(end) = []; end
    
    s = cellfun(@size,varargin,'unif',0);
    d = max(cellfun(@numel,s));
    for j = 1:numel(s), s{j}(end+1:d) = 1; end
    
    S = cat(1,s{:});
    s = max(S,[],1);
    ok = all(all(S == s | S == 1,1) | all(S <= 1,1));
        
    if ~ok
        names = arrayfun(@inputname,1:numel(varargin),'unif',0);
        if any(cellfun(@isempty,names)), names = 'Arguments';
        else, names = shortliststr(names);
        end
        ERR = MException('args:compatiblesize',[names ' have inconsistent size.']);
        throwAsCaller(ERR);
    end
        
    if outputsize, varargout{1} = s; return; end
    
    % Collapse empty dimensions, i.e. []+1 = []
    s(any(S == 0,1)) = 0;
    
    varargout = varargin(1:nargout);
    for j = 1:nargout
        if isequal(S(j,:),s), continue; end
        mult = s./S(j,:);
        mult(s == 0) = 0;
        varargout{j} = repmat(varargout{j},mult);
    end
end