function [W,e] = interpmatrix(Q,P,varargin)
% W = INTERPMATRIX(Q,P) - Given a set Q of d-dimensional n query-points (n x d array), and a set P 
%   of m nodal-points (m x d array); return an [n,m] sparse weight-matrix W such that:
%
%           W·z' ~= F(Q) : F = scatteredInterpolant(P,z,METHOD*);
%
%   Where z is an m-vector of values (for each nodal-point). METHOD is 'linear', by default, 
%   unless the '-smooth' flag or 'method1d',..,'methodNd',.. option pairs are used.
%
% [W,E] = INTERPMATRIX(..,'-extrap') : unless this flag is used, values outside the convex hull/ 
%   span of P will not be calculated. The n-vector of flags E will contain TRUE for points that
%   are extrapolated. NOTE that without the '-extrap' flag, all rows W(E,:) == 0
%
% [W,E] = INTERPMATRIX(..,'extrap',METHOD) : extrapolate, but using a different METHOD
%
% W = INTERPMATRIX(..,'-smooth') : [Currently slow!] Sets default METHOD to:
%       - makima/cubic (version dependent): for 1-dimensional / colinear points P
%       - natural-neighbor: for 2+ dimensional points P
%
% W = INTERPMATRIX(..,'method1d',M2,'methodNd',M3) : Directly provide interpolation method(s)
%   for 1D and/or ND scatettered data.
%
% [W,E] = INTERPMATRIX(..,'-sph') : for 3D point sets Q, P on the surface of the unit sphere, use
%   geodesic-arcs instead of euclidean distance for all algorithms. 
%
%   NOTE: currently only works for 'linear', returning spherical barycentric coordinates [1]
%       Unlike barycentric coordinates for flat triangles, these DO NOT ADD UP TO 1. 
%
% [W,E] = INTERPMATRIX(..,'maxarc',S) - When using the '-sph' flag, every point can be interpolated
%   (from one side or the other). If an explicit maximum arc distance S (radians) is provided,
%   points whose distance to the nearest nodal point > S will be flagged as extrapolated.
%
% TODO: Vectorize extrapolation for all cases
%       Vectorize smooth methods (e.g. natural neighbor,inverse distance, RBF/Kernel,...)
%       Smooth methods for spherical interpolation
%
% [1] P. Alfeld, M. Neamtu, and L. L. Schumaker, “Bernstein-Bézier polynomials on spheres and 
% sphere-like surfaces,” Computer Aided Geometric Design, vol. 13, no 4, pp. 333–349. 1996.
%
% See also: RESAMPLESTRUCTURE, AVGDOWNSAMPLE, FILTERSTRUCTURE

    [opt,varargin] = getflagoptions(varargin,{'-smooth','-sph','-extrap'});
    opt.maxarc = [];
    opt.tol = 0;
    
    if opt.smooth
        opt.methodNd = 'natural';
        if str2double(regexprep(version(),'(\d+\.\d+).*','$1')) <= 9.2
            opt.method1d = 'pchip';
        else
            opt.method1d = 'makima';
        end
    else
        opt.methodNd = 'linear';
        opt.method1d = 'linear'; 
    end
    opt = getpairedoptions(varargin,opt,'restchk');
    if islogical(opt.extrap) && opt.extrap
        extrap_method = 'linear';
    elseif ischar(opt.extrap)
        extrap_method = opt.extrap;
        opt.extrap = true;
    else
        extrap_method = '';
        opt.extrap = false;
    end
    if ~isempty(extrap_method)
        methods_Nd = {opt.methodNd,extrap_method};
        methods_1d = {opt.method1d,extrap_method};
    else
        methods_Nd = {opt.methodNd};
        methods_1d = {opt.method1d};
    end

    d = size(Q,2);
    nq = size(Q,1);
    np = size(P,1);
    
    assert(isnumeric(P) && isnumeric(Q) && size(P,2) == d,'Invalid P and/or Q');
    assert(size(unique(P,'rows'),1) == np,'Nodal-points must be unique');
    
    if opt.sph
        assert(np > 2 && d == 3,'Spherical interpolation requires N > 2 points in 3D');
    end

    P = double(P);
    Q = double(Q);
    
    switch np
    case 1
    % Set every Q(j) = P(1)
            W = ones(nq,np);
            e = true(nq,1);
    case 2
    % Interpolate linearly in direction P(1) -> P(2)
        u = P(2,:) - P(1,:);
            W = zeros(nq,np);
        W(:,2) = (Q-P(1,:))*u'/(u*u');
        W(:,1) = 1-W(:,2);
        W = sparse(W);
        
    otherwise
        if rank(P(2:end,:)-P(1,:)) < 2
        % Colinear points
        
            % TODO: simplification for linear case!
        
            u = P(2:end,:)-P(1,:);
            u = mean(u./rssq(u,2),1);
            p = (P-P(1,:))*u';
            q = (Q-P(1,:))*u';
            
            w = cell(1,np);
            for j = 1:np
                v = full(sparse(j,1,1,np,1,1));
                F = griddedInterpolant(p,v,methods_1d{:});                
                w{j} = sparse(F(q));
            end
            W = cat(2,w{:});
            return;
        end
        
        if ~strcmpi(opt.methodNd,'linear')
        % Slow & Smooth inerpolation: calculate the influence of each j in P over points Q
            w = cell(1,np);
            F = scatteredInterpolant(P,rand(np,1),methods_Nd{:});
            for j = 1:np
                    F.Values = full(sparse(j,1,1,np,1,1));
                w{j} = sparse(F(Q));
            end
            W = cat(2,w{:});
            return;
        end
        
        if ~opt.sph
        % General case: use Delaunay-Triangulation and barycentric coordinates

            DT = delaunayTriangulation(P);
            [tid,B] = pointLocation(DT,Q);

            e = isnan(tid); % query points that are outside convex hull

            % For interior points: barycentric coordinates B already provide weights...
            i = repmat(find(~e),1,size(B,2));
            j = DT(tid(~e),:);
            W = sparse(i,j,B(~e,:),nq,np);
            
            % For extrapolation....
            % TODO: seach closest edge/vertex only once, use barycentric coords. as above 
            if any(e) > 0 && opt.extrap

                if opt.sph
                   warning('Extrapolation of spherical data can cause unpredictable results!'); 
                end

                % get a list b of vertices at- or attached-to the convex hull...
                b = DT.convexHull;
                if ismatrix(b), b = unique(b(:)); end
                b = vertexAttachments(DT,b);
                b = DT([b{:}]',:);
                b = unique(b(:));
                nb = numel(b);

                w = cell(1,nb);
                F = scatteredInterpolant(P(b,:),rand(nb,1),methods_Nd{:});
                for j = 1:nb
                % calculate the influence of each b(j) over points e
                    F.Values = full(sparse(j,1,1,nb,1,1));
                    w{j} = sparse(F(Q(e,:)));
                end
                W(e,b) = cat(2,w{:});
            end
        else
        % Spherical variation: use alpha-shape (alpha~1) instead of convex hull, include origin 
        % to get (mostly) "radial" facets (a triangle on the surface + the origin).
        
            P = P./rssq(P,2);
            Q = Q./rssq(Q,2);
            [W,e] = spherical(P,Q,opt.maxarc,opt.tol);
        end
    end
end

function [W,e] = spherical(P,Q,maxarc,tol)
% Interpolate using alpha-shape (3D) of grid points (plus origin), with alpha ~ 1.0
       
    % The volume formed by the grid points must be expanded by R > 2·ALPHA, so that all query
    % points Q contained in the solid angle covered by alpha-shape(P,ALPHA) fall inside some
    % tetrahedra.
    R = 4;
    
    if ~isempty(maxarc)
        assert(isfinite(maxarc) && isreal(maxarc) && maxarc > 0 && maxarc < pi,...
            'Expecting maxarc interpolation distance in range (0,pi)');
        R = max(R,1/cos(maxarc/2));
    end
    ALPHA = R/2;

    Nv = size(P,1); % number of (known) vertices
    Nq = size(Q,1); % number of query points
    P(Nv+1,:) = 0;  % add origin as last point

    shp = alphaShape(R*P,ALPHA*R);
    shp.HoleThreshold = shp.volume();
    T = shp.alphaTriangulation;

    T(~any(T > Nv,2),:) = [];        % remove tetrahedra without origin
    T = sort(T,2);                   % sort so that last vertex for each facet is origin
    T = triangulation(T,R*P);
    
    % [tid,B] = pointLocation(T,Q);    % get interpolation weights as barycentric coords.
    % B(:,4) = [];                     % remove weight of last point (origin)
    % B = B*R;                         % ... and renormalize

    % PROVISIONAL?: pointLocation scales badly for large Nq: For a dense triangulation, chances
    %  are most points will lie in a triangle enclosed by their closest neighbors...
    [tid,B] = in_nearest_triangle(Q,P(1:end-1,:),T.ConnectivityList(:,1:3));

    valid = ~isnan(tid);
    if any(~valid)
    % ... call pointLocation only for the exceptions
        [tid(~valid),b] = pointLocation(T,Q(~valid,:));
        B(~valid,:) = b(:,1:3)*R;  
    end

    % B = B./sum(B,2);                 
    valid = ~isnan(tid);             % query points inside the alpha-shape

    if ~isempty(maxarc)
    % Check that the weighted interpolation distance < tol
        distk = @(k) rssq(P(T(tid(valid),k),:)-Q(valid,:),2); % distance to vertex k
        L = sum([distk(1),distk(2),distk(3)].*B(valid,:),2);
        valid(valid) = L <= maxarc;
    elseif any(~valid)
        x = acos(1/R)*2;
        warning('Points beyond default max. interpolation arc %0.2f rad (%0.1f°)',x,x*180/pi);
    end
    
    % Put together interpolation matrix
    i = repmat(find(valid),1,3);
    j = T(tid(valid),1:3);
    w = B(valid,:);
    
    if tol > 0
        % Remove vertices whose weight < tol 
        f = abs(w) >= tol;
        i = i(f);
        j = j(f);
        w = w(f);

        if any(~valid)
        % Allow query-points near grid points (within tol) even outside the alpha-shape
        
            [idx,d] = knnsearch(P(1:Nv,:),Q(~valid,:),'K',1,'distance','cosine');
            matches = acos(1-d) <= tol;
            if any(matches)
                j = [j(:);idx(matches)];
                w = [w(:);ones(nnz(matches),1)];
                idx = find(~valid);
                i = [i(:);idx(matches)];
                valid(~valid) = matches;
            end
        end
    end

    W = sparse(i,j,w,Nq,Nv);
    e = ~valid;
end

function [tid,B] = in_nearest_triangle(Q,P,T)
% See if the triangle formed by the three nearest-neighbors of Q is part of the spherical 
%   triangulation T, if so, get Q's barycentric coordinates in terms of that triangle, and 
%   check that Q is inside.
%
% FUTURE: check other, neighboring triangles?

    idx = knnsearch(P,Q,'K',3,'distance','cosine');
    [easy,tid] = ismember(sort(idx,2),T,'rows');
    T = T(tid(easy,:),1:3);
    B = zeros(size(Q,1),3);
    B(easy,:) = sphbarycoords(Q(easy,:),P(T(:,1),:),P(T(:,2),:),P(T(:,3),:));
    really_easy = all(B(easy,:) >= 0,2);
    easy(easy) = really_easy;
    tid(~easy) = NaN;
    B(~easy,:) = NaN;
end

function B = sphbarycoords(q,p1,p2,p3)
% Spherical barycentric coordinates for point(s) q in triangle(s) p1-p2-p3.
% From:
% Lei, K., Qi, D., Tian, X., 2020. A New Coordinate System for Constructing Spherical Grid Systems. 
% Applied Sciences 10, 655. https://doi.org/10.3390/app10020655

    dotq1 = dot(q,p1,2);
    dotq2 = dot(q,p2,2);
    dotq3 = dot(q,p3,2);

    dot12 = dot(p1,p2,2);
    dot23 = dot(p2,p3,2);
    dot13 = dot(p1,p3,2);
    x23 = cross(p2,p3,2);
    x31 = cross(p3,p1,2);
    x12 = cross(p1,p2,2);

    S = atan2(dot(p1,x23,2),1 + dot12 + dot23 + dot13);
    B(:,1) = atan2(dot(q,x23,2),1 + dotq2 + dot23 + dotq3)./S;
    B(:,2) = atan2(dot(q,x31,2),1 + dotq3 + dot13 + dotq1)./S;
    B(:,3) = atan2(dot(q,x12,2),1 + dotq1 + dot12 + dotq2)./S;
end
