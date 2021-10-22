function [p,bounded,fevals] = bisection(f,a,b,varargin)
% P = BISECTION(F,A,B,TOL) - Find a bounded root of function F(x) in interval (A,B) within 
%   tolerance TOL by the interval bisection method.
%
%   F - Function handle (scalar or element-wise function).
%   A, B - Interval bounds (scalars or equal-sized arrays). Should include at least one zero, i.e.
%       F(A)·F(B) <= 0 for all A, B. Closed intervals (a == b) and no-root intervals F(A)·F(B) > 0
%       do not cause the method to fail, but just to issue warnings.
%   TOL - Search tolerance Tr/[Tx,Ty]/[Tx,Ty,Tr], parsed as [~,fx,fy] = PARSETOLERANCE(TOL).
%       Convergence criterion is |A-B| <= fx(P) | |f(P)| <= ty(f(P))
%
% P = BISECTION(F,A,B,TOLX,TOLY) - Use pre-parsed tolerance functions TOLX(X) and TOLY(Y).
%
% See also: PARSETOLERANCE, FZERO
    
    narginchk(3,5);
    if nargin < 4, varargin{1} = []; end
    if numel(varargin) < 2
        [~,tolx,toly] = parsetolerance(varargin{1},'minreltol',eps(1),'minabstol',eps(0));
    else
        [tolx,toly] = deal(varargin{:});
    end

    fa = f(a); 
    fevals = 1;
    
    [a,b] = compatiblesize(a,b,fa);
    % assert(isequal(size(a),size(b)),'a and b must have equal size');
    assert(numel(fa) == numel(a),'element-wise function required');
    
    if isempty(a), p = a; return; end
    
    fb = f(b);
    p = (a + b)/2;
    fp = f(p);
    fevals = fevals+2;
    
    % Handle roots at the extremes: close interval on root and set fa = fb = fp = 0
    d = ndims(a)+1;
    [F,idx] = min(abs(cat(d,fa,fb,fp)),[],d);  % idx = j > 0 means a/b/p is already a solution
    unsolved = F > toly(F);

    closed = (b == a) & unsolved;
    if any(closed,'all')
       warning('bisection:closedintervals','Closed (zero-distance) intervals included');
    end

    bounded = fa.*fb <= 0;
    notbounded = ~bounded & unsolved;
    % if any(notbounded)
    % % see if, by any chance, a -> (a+b)/2, or b -> (a+b)/2 are propper intervals
    %     lucky = notbounded & (fa.*fp < 0);
    %     needsupdate = false;
    %     if any(lucky)
    %         b(lucky) = p(lucky); 
    %         fb(lucky) = fp(lucky);
    %         bounded(lucky) = true;
    %         needsupdate = true;
    %     end
    %     lucky = notbounded & (fb.*fp < 0);
    %     if any(lucky)
    %         a(lucky) = p(lucky); 
    %         fa(lucky) = fp(lucky);
    %         bounded(lucky) = true;
    %         needsupdate = true;
    %     end
    %     if needsupdate
    %         p = (a + b)/2;
    %         fp = f(p);
    %         fevals = fevals + 1; % debug
    %     end
    % end

    if any(notbounded)
        warning('bisection:notbounded','At least one root is not included in provided interval');
    end
    
    % save the best evaluated point so far as 'root', for ~unsolved and notbounded points...
    idx(unsolved & ~notbounded) = 0; 
    for j = 1:3
        s = (idx == j);
        if ~any(s), continue; end
        switch j
            case 1, p(s) = a(s); b(s) = a(s); fa(s) = 0; fb(s) = 0; fp(s) = 0;
            case 2, p(s) = b(s); a(s) = b(s); fa(s) = 0; fb(s) = 0; fp(s) = 0; 
            case 3, a(s) = p(s); b(s) = p(s); fa(s) = 0; fb(s) = 0; fp(s) = 0; 
        end
    end
    clear F idx

    convergedonx = false(size(p)); convergedony = false(size(p));
    while ~all(convergedony(bounded) | convergedonx(bounded))
       neg = (fa.*fp < 0);
       b(neg) = p(neg); 
       fb(neg) = fp(neg);
       a(~neg) = p(~neg); 
       fa(~neg) = fp(~neg);
       p = (a + b)/2; 
       fp = f(p);
       fevals = fevals + 1; % debug
       convergedonx = abs(a-b) <= tolx(p);
       convergedony = abs(fp) <= toly(fp);
    end
end
