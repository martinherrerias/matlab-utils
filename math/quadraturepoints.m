function [X,W,n] = quadraturepoints(n,a,b,u)
% QUADRATUREPOINTS - Wrapper for LGWT (Legendre-Gauss Nodes and Weights) for M-dimensional 
%   and/or per-unit input.
%
% [X,W] = QUADRATUREPOINTS(N,A,B) - (scalar A,B,N) Identical to LGWT, Returns the Legendre-Gauss 
%   nodes X and weights W on an interval [A,B] with truncation order N.
%
% [X,W] = QUADRATUREPOINTS(N,A,B) - (M-vector A,B,N) Returns an M-dimensional grid of nodes X and 
%   weights W along intervals [A(j),B(j)] with truncation order N(j) for each dimension j = 1:M.
%
% [X,W] = QUADRATUREPOINTS(-D,A,B,[U]) - (Note D < 0) Use -D points per unit |B-A|, i.e. set
%   N = ceil|(B-A)·D/U |, where optional U defines the scale (unit) of [A,B].
%
%   Note that N and/or U can be M-vectors, allowing for different spacing along dimensions.
%   Mixing of positive-integers (absolute number of points) and negative floats (points-per-unit) 
%   is also possible.
%
% See also: LGWT, SAMPLESYSTEM, IMPORTPOLYGONSFILE

    narginchk(3,4)
    if nargin < 4 || isempty(u), u = 1; end
        
    m = max([numel(n),numel(a),numel(b),numel(u)]);

    a = check(a,m,'A'); b = check(b,m,'B'); u = check(u,m,'U'); n = check(n,m,'N');
        
    % turn per-unit values into absolute number of points
    pu = n < 0;
    n(pu) = ceil(abs((b(pu)-a(pu)).*n(pu)./u(pu)));
    
    assert(all(n > 0 & mod(n,1)==0),'N must contain either positive integers or negative floats');
    
	[x,w] = arrayfun(@lgwt,n,a,b,'unif',0);
	[x{:}] = ndgrid(x{:});
    [w{:}] = ndgrid(w{:});
	
    W = ones(size(w{1}));
    for j = 1:m
        x{j} = x{j}(:);
        W = W.*w{j};
    end
    X = [x{:}];
    W = W(:);
        
end

function x = check(x,m,tag)
    if ~(isnumeric(x) && isreal(x) && all(isfinite(x)) && any(numel(x) == [1,m]))
        ERR = MException('quadpts:check','%s is not a real numeric scalar/%d-vector',tag,m);
        throwAsCaller(ERR);
    end   
    if m > 1 
        if isscalar(x), x = repmat(x,m,1); else, x = x(:); end
    end
end