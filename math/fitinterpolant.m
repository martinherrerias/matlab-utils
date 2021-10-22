function [x,V,g] = fitinterpolant(F,x,t,G,E)
% [x,V,g] = FITINTERPOLANT(F,{x1g,x2g,..xng},TOL,[METHOD]) - generate a non-uniform grid that 
%   approximates function F(X1,X2,..Xn) with tolerance TOL, starting with the grid of evaluation 
%   points NDGRID(x1g,x2g,...,xng) and refining it by progressive bisection along individual 
%   intervals in all dimensions, until for any new set of bisecting points Xq = {X1q,X2q,..Xnq} 
%   the interpolated values g(Xq{:}) = interpn(x{:},V,Xq{:},METHOD) are close enough to the 
%   function values F(Xq{:}). That is, |g(Xq{:}) - F(Xq{:})| <= TOL. 
%
% [x,V,g] = FITINTERPOLANT(F,{x1g,x2g,..xng},TOL,G,E) - Uses a custom interpolation function 
%   G(x{:},V,Xq{:}) instead of default @INTERPN, and/or a custom error function E(Vf,Vg), where 
%   Vf = F(Xq{:}) and Vg = g(Xq{:}), instead of the default @(vf,vg) abs(vf-vg). This allows F  
%   and G to output objects of any class (e.g. vectors, structures, classes,...) as long as a 
%   meaningful metric of error E can be defined between those objects.
%
% NOTES:
%   (*) The function is expected to be convex within each of the intervals defined by the original
%     grid points. If this is not the case, it is possible that grid refinemet stops prematurely 
%     at a local error-function zero. Include any known inflection points in the grid vectors,
%     and/or start with a sufficiently fine grid.
%
%   CAUTION: F is expected to be computationally expensive (the reason why you might want an  
%     interpolant on the first place). Repeated evaluations are avoided by keeping a list of all 
%     previously evaluated points, yet things can get quickly out of hand for problems with many 
%     dimensions, very convoluted functions, or too tight tolerances. Start small, and test your
%     way forward.
%
%   FITINTERPOLANT uses a simple greedy algorithm, testing all possible new 'slices' of the grid,
%   and keeping those with the highest mean error. The result is a non-optimal, non-uniform full
%   grid, which limits the choice of interpolating METHOD, but reduces the grid size and number of
%   function evaluations while preserving a simple implementation. 
%
%   It works well for smooth functions, in cases where a GRIDDEDINTERPOLANT would be used in any
%   case, but the required grid resolution is not know in advance. For functions with localized
%   strong non-linearities, or with many dimensions, consider using sparse-grid-interpolation,
%   e.g. A. Klimke and B. Wohlmuth, 2005. DOI: 10.1145/1114268.1114275
%
% INPUT:
%   F(X1q,X2q,..Xnq) - function handle, for N-variable function to be approximated. Must take N
%       same-size N-dimensional numeric arrays X1q,X2q,..Xnq and output also an N-dimensional
%       array of the same size as its inputs, although NOT NECESSARILY NUMERIC ($).
%   
%   {x1g,x2g,...,xng} - interpolation grid vectors in compact form (i.e. cell array of vectors).
%       Usually an N-cell array of 2-vectors, with lower/upper limits along each dimension.
%       For non-convex functions it might be necessary to insert intermediate points (*).
%
%   TOL - scalar, tolerance threshold value. Iteration stop criterion is  max(E(Vfq,Vgq)) <= TOL
%       where Vfq, Vgq are the function and interpolant values, respectively, evaluated on all new
%       points of a candidate splitting plane.
%
%   G(x1g,x2g,...,xng,V,X1q,X2q,..Xnq) - interpolation function handle, that that takes a set of 
%       grid vectors x1g,x2g,...,xng, the N-dimensional array of function values V, and a set of 
%       N-dimensional query point arrays X1q,X2q,..Xnq, to return an array of interpolated values.
%       The DEFAULT is G = @interpn. See INTERPN for details.
%
%   E(Vf,Vg) - error function handle, that takes same-size N-dimensional arrays Vf = F(Xq{:}) and 
%       Vg = G(..,Xq{:}) and returns a same-size array of positive values comparable to TOL.
%       DEFAULT is E = @(vf,vg) abs(vf-vg).
%
%   ($) NOTE: Vf,Vg need not be numeric, or even of the same type. The only restriction is that
%       Xq{:}, Vf = F(Xq{:}), Vg = G(..,Xq{:}) and E(Vf,Vg) all have the same size, and are
%       indexable. E(Vf,Vg) must be numeric/logical, as will be used to check the stop criterion 
%       max(E(Vf(idx),Vg(idx))) <= TOL.
%       
% OUTPUT:
%   {x1g,x2g,...,xng} - refined cell-array of interpolation grid vectors.
%
%   V - function values V = F(X1,X2,...,Xn) where [X1,X2,...,Xn] = NDGRID(x1g,x2g,...,xng)
%
%   g - ready to use function handle to replace F, g(Xq{:}) = G(x{:},V,Xq{:});
%
% EXAMPLES:
%   1. Approximate the function F(y,z) = y² + 2z in the interval y = [0,1], z = [0,1]:
%
%     F = @(y,z) y.^2 + 2*z;            % (!) NOTE: function must take/return same-size ND arrays!
%     x = {[0,1],[0,1]};                % start with interval boundaries
%     t = 0.01;                         % Use a 1% approximation tolerance
%     [x,V,g] = fitinterpolant(F,x,t);  % FITINERPOLANT call. Note that x{2} remains [0,1] 
%                                       % due to linear z, whereas x{1} (y) has now 9 grid points.
%     [X{1:2}] = ndgrid(x{:}); 
%     surf(X{:},V);                     % Plot approximating surface
%
%   2. Approximate the complex function F(a,b) = cos(a) + i sin(b) for a,b in [0,pi]
%      NOTE that defaults G = @INTERPN and E = @(vf,vg) abs(vf-vg) work fine with complex numbers.
%
%     F = @(a,b) complex(cos(a),sin(b));
%     x = {[0,pi/2,pi],[0,pi]};         % (!) NOTE: since cos(pi/2) = (cos(0) + cos(pi))/2
%                                       % the starting grid vector for A must be [0,pi/2,pi], and
%                                       % not [0,pi], as the first bisection (pi/2) would yield 
%                                       % F(pi/2,..) == G(pi/2,..) and splitting along dimension
%                                       % A would stop. In general, NON-CONVEX FUNCTIONS SHOULD 
%                                       % BE PRE-SPLIT AT INFLECTION POINTS.
%     t = 0.001;
%     [x,V,g] = fitinterpolant(F,x,t);
%
%     [X{1:2}] = ndgrid(x{:});              % Plot real/imaginary components
%     surf(X{:},real(V),'facecolor','b'); 
%     hold on; surf(X{:},imag(V),'facecolor','r');
%
% See also: INTERPN, NDGRID, GRIDDEDINTERPOLANT

    XX = []; VV = []; % locally persistent variables for EVALFUN
    MAXITER = 1000;

    narginchk(3,5);
    isnice = @(x) isvector(x) && numel(x) > 1 && issorted(x) && ~any(isnan(x) | isinf(x));
    assert(iscell(x) && all(cellfun(isnice,x)),'Bad grid vectors');
    assert(isa(F,'function_handle'),'Expecting function handle F');
    assert(isscalar(t) && t > 0,'Expecting scalar > 0 tolerance t');

    if nargin < 4 || isempty(G), G = @interpn; end
    if nargin < 5 || isempty(E), E = @(vf,vg) vf-vg; end
    if ischar(G)
        METHOD = G;
        G = @(varargin) interpn(varargin{:},METHOD);
    end
    assert(isa(G,'function_handle'),'Expecting function handle G');
    assert(isa(E,'function_handle'),'Expecting function handle E');
        
    d = numel(x); % number of dimensions
    X = cell(1,d);
    [X{:}] = ndgrid(x{:});
    try V = evalfun(F,X); 
    catch ERR
        struc2str = @(x) strjoin(strtrim(strsplit(evalc('disp(x)'),newline())),', ');
        error('Function handle F crashed with provided limits, ERR: %s',struc2str(ERR)); 
    end
    
    Xq = cell(1,d);
    for iter = 1:MAXITER
        N = cellfun(@numel,x);

        % Get function values at grid mid-points
        xq = cellfun(@interpn,x,'unif',0);     % refine grid inserting middle points
        [Xq{:}] = ndgrid(xq{:});
        Vq = evalfun(F,Xq);                    % get function values (avoiding repeated calls)
        Vg = G(x{:},V,Xq{:});                  % get interpolant values
        E2q = abs(E(Vq,Vg));                   % ... and interpolation error

        % indices of original planes: old{i}(j) == slice j on dimension i is not new
        old = arrayfun(@(n) mod(1:n,2)>0,2*N-1,'unif',0);

        % Decide which new query-planes (slices) are required, and which aren't:
        eq = slicestats(E2q,t);                         % get mean and max errors for each slice
        keep = old;                                     % keep{i}(j) == keep slice j on dimension i
        idx = NDmap.idx2sub(2*N-1,1:prod(2*N-1));       % for quick indexing of E2q
        anynew = false;
        while ~isempty(eq) && any(eq.max > t)
            [~,i] = max(eq.mean);                       % find the slice with max MAE
            keep{eq.dim(i)}(eq.slice(i)) = true;        % mark as required
            E2q(idx(:,eq.dim(i)) == eq.slice(i)) = 0;   % set errors on that slice to zero
            eq = slicestats(E2q,t);                     % recalculate stats
            anynew = true;
        end

        if ~anynew, break; end

        % Update base grids to include new (kept) slices
        x = cellfun(@interpn,x,'unif',0);
        x = arrayfun(@(j) x{j}(keep{j}),1:d,'unif',0);
        V = Vq(keep{:});
    end
    
    info = sprintf(['Interpolant with tol = %0.2g fitted to %s. \n',...
                '%d inputs required.'],t,func2str(F),numel(x));
    
    % Return ready-to-use handle to interpolant
    g = @interpwrapper;
    % g = @(varargin) G(x{:},V,varargin{:});

    function g = interpwrapper(varargin)
    % 
        if nargin ~= numel(x), error(info); end
        [varargin{:}] = compatiblesize(varargin{:});
        sz = size(varargin{1});
        type = class(V);
        for j = 1:numel(varargin)
            varargin{j} = typecast(varargin{j},type);
        end
        g = G(x{:},V,varargin{:});
        g = reshape(g,sz);
    end
    
    function V = evalfun(F,Xq)
    % EVALFUN(F,XQ) - return F(Xq{:}) for function handle F at grid points Xq, avoiding repeated
    %   evaluations of (costly) F by keeping a record of any previous evaluations:

        s = size(Xq{1});

        % Turn Xq into vectors, and pile them into a matrix Q
        Xq = cellfun(@(x) reshape(x,[],1),Xq,'unif',0);
        Q = cat(2,Xq{:});

        % If no records of previous evaluations exist, evaluate all
        if isempty(XX)
            XX = Q; 
            VV = F(Xq{:});
            V = reshape(VV,s);
            return;
        end

        % Find query points in list of previous evaluations
        [pasteval,pastidx] = ismember(Q,XX,'rows');
        %idx = idx(old);
        Xq = cellfun(@(x) x(~pasteval),Xq,'unif',0);
        vq = F(Xq{:});

        pastidx(~pasteval) = numel(VV)+(1:numel(vq))';
        XX = cat(1,XX,Q(~pasteval,:));
        VV = cat(1,VV,vq);

        V = reshape(VV(pastidx),s);
    end
end

function s = slicestats(E,t)
% Take d-dimensional abs. error array E, and calculate its mean and max along any possible slice.
% Return a table {mean,max,dim,slice}, leaving out rows/slices within tolerance, i.e. max < t

    d = ndims(E);
    for j = d:-1:1
        % 'flatten' m, M accross all dimensions but j
        m = mean(E,[1:j-1,j+1:d]);
        M = max(E,[],[1:j-1,j+1:d]);
        % m = E; M = E;
        % for l = [1:j-1,j+1:d]     
        %     m = mean(m,l);
        %     M = max(M,[],l);
        % end
        s{j}(:,1) = squeeze(m);   % mean of errors along all dimensions but j
        s{j}(:,2) = squeeze(M);   % max of ...  
        s{j}(:,3) = j;            % dimension index
        s{j}(:,4) = 1:numel(m);   % query-plane index
    end
    s = cat(1,s{:});
    s(s(:,2) < t,:) = [];         % remove query-planes already within tolerance
    
    s = array2table(s,'VariableNames',{'mean','max','dim','slice'});
end