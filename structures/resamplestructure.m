function [B,F] = resamplestructure(A,m,varargin)
% B = RESAMPLESTRUCTURE(A,m) - Take an array A, and resample its N rows at m equal intervals, 
%   to get an m(N-1)+1 row array.
%
% [B,F] = RESAMPLESTRUCTURE(..) - Return a sparse interpolation matrix F, such that B = F·A.
%
% [R,F] = RESAMPLESTRUCTURE(S,m,N) - where S is a structure, applies RESAMPLESTRUCTURE(A,m) to 
%   any fields and sub-fields A of S which are N-row numeric arrays. In that case R will be a 
%   copy of the structure S, with (some) fields expanded to m(N-1)+1 rows.
%
%   If N is not provided, is will be set to N = mode(n(n>1)), where n is the number of rows of
%   all fields and subfields of S.
%
% .. RESAMPLESTRUCTURE(A,m,..,'-centered') - Use sampling points t(j)-1/2:1/m:t(j)+1/2 instead of
%   t(j):1/m:t(j+1) so that the result is an N·m row array, and A = AVGDOWNSAMPLE(B,m). Points at
%   the edges are extrapolated linearly.
% 
% R = RESAMPLESTRUCTURE(S,m,..) - Where m is a rational number m = P/Q or a [P,Q] integer pair,
%   will resample by P and then downsample by Q, to return floor(PN/Q) points, if centered, or
%   floor(P(N-1)/Q + 1/Q) otherwise. 
%
% B = RESAMPLESTRUCTURE(..,'-full') - When downsampling (Q > 1), extrapolate the required samples
%   at the end to return ceil(PN/Q) or ceil(P(N-1)/Q + 1/Q)
%
% .. RESAMPLESTRUCTURE(..,'offset',X) - Apply an offset X/m to the sampling points obtained by 
%   any of the methods above. That is, make sure that the first resampled point is at X/m from
%   A(1,:). A negative offset results in extrapolation. E.g.:
%
%       resamplestructure([0;1],2,'offset',0) = [0 0.5 1]'
%       resamplestructure([0;1],2,'offset',0.5) = [0.25 0.75 1.25]'
%       resamplestructure([0;1],2,'offset',-1) = [-0.5 0 0.5]'
%
% .. RESAMPLESTRUCTURE(..,'offset',[a,b]) - Uses offset a/m for the first point, and b*/m for the
%   last point. b* will exactly match b as long as the length L = b/m+(N-1)-a/m of the resulting
%   interval is exactly divisible by 1/m. Otherwise b will be rounded up or down depending on the
%   '-full' flag (see below). E.g.
%       resamplestructure([0;1],5,'offset',[-1,1]) = -0.2:0.2:1.2'
%       resamplestructure(..,'offset',[-m/2+1/2,m/2-1/2]) is equivalent to using '-centered'
%
% See also: AVGDOWNSAMPLE, FILTERSTRUCTURE, INTERPMATRIX

    if nargin == 0 && debugging(), test(); return; end % DEBUG
    narginchk(2,6);
    
    assert(isnumeric(m) && isreal(m) && all(m > 0),'Expecting real, positive M');

    [opt,varargin] = getflagoptions(varargin,{'-centered','-full','-omitnan'});
    opt.offset = 0;
    opt.tol = 1e-6;
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    if isscalar(m)
        if mod(m,1) == 0, q = 1; 
        else
            [p,q] = rat(m,opt.tol*m);
            if p/q ~= m
               warning('Using rational approximation %g ~ %d/%d (%g error)',m,p,q,m-p/q);
            end
            m = p;
        end
    else
        assert(numel(m) == 2 && all(mod(m,1) == 0),'Expecting scalar M or rational [N,D]');
        m = m/gcd(m(1),m(2));
        q = m(2); 
        m = m(1);
    end

    if ~isempty(varargin)
        assert(numel(varargin) == 1,'Unrecognized arguments');
        N = varargin{1};
    else
        if isnumeric(A) || isdatetime(A) || isa(A,'tabular') || ...
                (isobject(A) && ismethod(A,'filterstructure'))
            N = size(A,1);
        else
            N = cellfun(@(f) size(getnestedfield(A,f),1),nestedfieldnames(A));
            N(N<=1) = [];
            assert(~isempty(N),'Failed to find N > 1');
            N = mode(N);
        end
    end
    assert(N > 0,'Cannot resample single-row');

    parsestruct(opt,{'tol','offset'},'-n','-r','-f');
    assert(opt.tol >= 0 && opt.tol <= 0.5,'Invalid tolerance');
    
    switch numel(opt.offset)
        case 1, opt.offset(2) = opt.offset(1);
        case 2, opt.offset = opt.offset(:)';
        otherwise, error('Expecting scalar or 2-vector offset')
    end
    
    if q > 1
        opt.offset = opt.offset*q;
        opt.tol = opt.tol*q;
    end
    
    if opt.centered
        opt.offset = opt.offset + [- m/2 + 1/2, m/2 - 1/2];
    else
        opt.offset = opt.offset + [-1,1]*(q-1)/2;
    end
    
    % ensure interpolation interval is an integer multiple of m/p
    L = opt.offset(2)+(N-1)*m - opt.offset(1) + 1; % length of interval, in units of m
    dL = mod(L,q);
    if dL == 0
        % already a perfect fit
    elseif opt.full || L < q || (q-dL) < opt.tol
        opt.offset(2) = opt.offset(2) + (q-dL);
        L = L + (q-dL);
    else
        opt.offset(2) = opt.offset(2) - dL;
        L = L - dL;
    end
    
    % make sure last point is not skipped due to numerical precission issues
    opt.offset(2) = opt.offset(2) + 1/(2*m); 
    
    if m == q && all(opt.offset == 0)
        B = A; 
        if nargout > 1, F = speye(N); end
        return; 
    end
    
    % Original N points can be imagined as a grid of units m, that is 0:m:(N-1)m
    % query points then range from offset(1):1:(N-1)m + offset(2)
    
    % irregular query points at the start -- that is, up to m (second point) for offset <= m, 
    % or up to ceil(offset)·m, for offsets > m -- will produce the first b+1 columns of F
    b = max(1,ceil(opt.offset(1)/m)); 
    b = min(b,N-1); 
    xb = opt.offset(1):min(b*m,opt.offset(2)+(N-1)*m);
    B = interpmatrix(xb',(0:m:b*m)');  % start kernel
    
    % irregular sample points at the end -- that is, from floor(offset+N)m for offset < -1, 
    % or from (N-2)m (second-last point) for offsets >= -1 -- will produce the last e+1
    % columns of F
    % NOTE: e might be equal to b, or even b-1, either because N = 2, or offset(1) > (N-1)m,
    % the use of max(e,b) for xe avoids redundancy with xb...
    e = min(N-2,floor(opt.offset(2)/m+N-1));
    n = (e-b);
    xe = xb(end)+max(0,n)*m+1:opt.offset(2)+(N-1)*m;    
    % xe = fliplr(opt.offset(2)+(N-1)*m:-1:max(e,b)*m);
    E = interpmatrix(xe',(e*m:m:(N-1)*m)');   % end kernel
    
    n = (e-b);
    if ~isempty(xe) && xb(end) + max(0,n)*m == xe(1)
    % Remove the first row of E if redundant
        E(1,:) = [];
    end
    
    if n > 0
    % regular query points from b to e -- n = (e-b) repetitions of the same kernel matrix -- 
    % will make out the central mn × n+1 block of F
        
        x = xb(end)+1:(b+1)*m;
        k = interpmatrix(x',[b*m;(b+1)*m]);   % center kernel [m,2]
        
        assert(nnz(B)+nnz(E)+n*numel(k) < maxarraysize()/3,'Not enough memory');

        ii = repmat(1:n*m,2,1)' + size(B,1);
        jj = repelem(1:n,m)'+[0,1] + size(B,2)-1;
        s = size(B)+size(E)+[m*n,n-1];
        F = sparse(ii,jj,repmat(k,n,1),s(1),s(2),nnz(B)+nnz(E)+n*numel(k));
    else
    % n = -1 means that the two first columns of E are the two last columns of B
    % n = 0 means that the first column of E is the last column of B
        s = size(B)+size(E)+[0,n-1];
        F = spalloc(s(1),s(2),nnz(B)+nnz(E));
    end

    F(1:size(B,1),1:size(B,2)) = B;
    if ~isempty(E), F(end-size(E,1)+1:end,end-size(E,2)+1:end) = E; end

    if q > 1
    % Downsample by q
        K = repelem(speye(L/q)./q,1,q);
        F = K*F;
    end
    
    if isempty(A), B = A; return; end
    
    if isnumeric(A)
        if ~ismatrix(A)
            sz = size(A); 
            A = reshape(A,N,[]); 
        else, sz = []; 
        end
        
        if opt.omitnan
            W = isnan(A);
            A(W) = 0;
            W = double(~W);
        end
        
        % B = F*A (MTIMES (*) is not supported for one sparse argument and one single argument)
        B = zeros(size(F,1),size(A,2),'like',A);
        B(:) = F*double(A);
        
        if opt.omitnan
           B(:) = B./(F*W); 
        end
        
        if ~isempty(sz)
            B = reshape(B,[size(B,1),sz(2:end)]); 
        end
    else
        B = filterstructure(A,F,N);
    end
end

function test()
    
    N = 50;
    M = 20; % make sure M < N, for downsampling up to 1/M to make sense

    fprintf('resamplestructure test: %d points, M in [1/%d,%d]\n',N,M,M);
    maxerror = 0;
    for j = 1:1000
        m = randi(M,1,2);   % rational M(1)/M(2)
        b = round(1-2*rand(1,2),2);
        x = resamplestructure((0:N-1)',m,'offset',b)';
        x0 = b(1)*m(2)/m(1):m(2)/m(1):(N-1)+b(2)*m(2)/m(1);
        maxerror = max(maxerror,max(abs(x-x0)));
        if maxerror > 1e-12
            keyboard();
        end
    end
    fprintf('max. error = %g\n',maxerror);
end