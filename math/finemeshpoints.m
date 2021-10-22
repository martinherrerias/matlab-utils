function [V,P,w] = finemeshpoints(T,n)
% [V,P,W] = FINEMESHPOINTS(T,n) - return O(N·n²/2) vertices V given triangulation T, by
%   dividing the N triangles into (n-1)² smaller facets (inserting vertices at barycentric 
%   coordinates {i/n, j/n, (1-i-j)/n}). 
%   Index matrix P will contain indices of the n(n+5)/2 vertices that "belong" to each original
%   triangle. That is, vertices P(k,:) are within (or at the outline of) triangle T(k,:). 
%   W is an n(n+5)/2 weight vector. Weights are proportional to 1 for original vertices (corners)
%   3 for new vertices on the edges, and 6 for new vertices inside the triangle.

    if nargin == 0, test(); return; end

    V = single(T.Points);
    E = T.edges;
    P = T.ConnectivityList;
    
    nV = size(V,1);
    nE = size(E,1);
    nT = size(P,1);
    
    nmax = floor(sqrt(maxarraysize(class(V))/(3*nT))-2)/2;
    validateattributes(n,{'numeric'},{'integer','scalar','nonnegative'});
    if n > nmax
       warning('Not enough memmory, using n = %d',nmax);
       n = nmax;
    end
    if n == 0, w = [1,1,1]/3; return; end

    e0 = linspace(0,1,n+2);
    e0 = e0(2:end-1);
    
    % add nxe new vertices on every edge (barycentric coords e0, 1-e0, 0
    nxe = numel(e0);
    M = V(E(:,1),:).*shiftdim(e0,-1)+ V(E(:,2),:).*shiftdim(1-e0,-1);
    M = reshape(permute(M,[1,3,2]),nE*nxe,[]);
    iM = reshape(nV+(1:nE*nxe)',nE,nxe);
    
    % get barycentric coordinates of "middle" points
    [r,c] = find(e0+e0' < 1);
    c0 = [e0(r);e0(c);1-e0(r)-e0(c)]';
    nxc = size(c0,1);

    % ... and nxc new vertices inside each triangle
    C = T.barycentricToCartesian(repmat((1:nT)',nxc,1),repelem(c0,nT,1));
    iC = reshape(iM(end)+(1:nT*nxc)',nT,nxc);
    
    V = [V;M;C];
    P(:,3+3*nxe+(1:nxc)) = iC;

    edgemap = sparse(E,fliplr(E),repmat((1:nE)',2,1)',nV,nV);
    for j = 1:3
       e = edgemap(sub2ind([nV,nV],P(:,j),P(:,mod(j,3)+1)));
       P(:,3+(j-1)*nxe+(1:nxe)) = iM(e,:);
    end
    w = [1,1,1,repmat(3,1,3*nxe),repmat(6,1,nxc)];
    w = w/sum(w);
end

function test()
    N = 42;
    M = 9;

    V = rand(N,2);
    T = triangulation(delaunay(V),V);

    clf(); hold on;
    triplot(T);

    [V,P] = finemeshpoints(T,M);
    plot(V(:,1),V(:,2),'g.')
    
    k = randi(size(P,1));
    plot(V(P(k,:),1),V(P(k,:),2),'m.')
end
