function [X,W,T] = sphquadraturepoints(N)
% [X,W,T] = SPHQUADRATUREPOINTS(N) - Returns 10·4^Z ~ N points X quasi-uniformly distributed on the
%   hemispherical unit dome, along with quadrature weights W for 1st order integration. X are the
%   centroids of the triangles of a quasi-regular polyhedral convex hull T (see below), and W their
%   solid angles.
%
%   The particular set 10·4^Z = {40, 160, 640, 2560, 10240,..} arises from starting with half an
%   icosahedron, and recursively inserting vertices at edge midpoints, to achieve a quasi-regular
%   triangulation T that covers the half dome.
%
%   If N is omitted, 1/SimOptions.RelTol will be used. Z is set as max{40,ceil[ log4((N-2)/10) ]}
%   NOTE that whether the actual integration error falls within RelTol depends on the function to
%   be integrated. For standard irradiance distrubutions, numerical tests suggest that setting 
%   N = 1/tol does ensure max. errors < 1/N, with RMS ~ 0.3/N. 
%
% FUTURE: should be replaced by integration with quadrature values, ideally at a higher-order
%   e.g. J. A. Reeger and B. Fornberg, “Numerical Quadrature over the Surface of a Sphere: 
%   Numerical Quadrature over the Surface of a Sphere,” Studies in Applied Mathematics,
%   vol. 137, no. 2, pp. 174–188, Aug. 2016.
%
% See also: SPHEREPOINTS, QUADRATUREPOINTS

    if nargin < 1 || isempty(N), N = 1./getSimOption('RelTol'); end
    Z = ceil(log(max(40,N)/10)/log(4));
    N = 10*4^Z;
    
    % Save (Load) precalculated results for faster execution at runtime
    file = sprintf('sphere_qpts_%d.mat~',N);
    file = fullfile(fileparts(mfilename('fullpath')),file);
    if ~isempty(dir(file))
        load(file,'-mat','X','W','T');
        return;
    end

    % Start with 2+10*4^Z vertices (with icosahedral symmetry)
    P = spherepoints(N+2,'regular',true,'symmetric',true);
    
    % solve numerical precission issue near XY plane
    TINY = 1e-6*N.^(-0.5);
    P(abs(P(:,3)) < TINY,3) = 0;
     
    % Triangulate, and keep the upper half
    P = sortrows(P,3,'descend');
    T = convhull(P);
    T(any(reshape(P(T,3) < 0,[],3),2),:) = [];
    P(max(T(:))+1:end,:) = [];
    T = triangulation(T,P);
    
    % Set quadrature points as the centroids of the triangles
    P = permute(reshape(P(T.ConnectivityList(:),:),[],3,3),[1 3 2]);
    X = single(mean(P,3));
    
    % ... and set weights as the solid angle of each triangle
    W = 2*atan(dot(P(:,:,1),cross(P(:,:,2),P(:,:,3),2),2)./...
        (1+dot(P(:,:,1),P(:,:,2),2)+dot(P(:,:,2),P(:,:,3),2)+dot(P(:,:,1),P(:,:,3),2)));
    W = single(W);
    
    % a = vectorangles(P(:,:,1),P(:,:,2));
    % b = vectorangles(P(:,:,2),P(:,:,3));
    % c = vectorangles(P(:,:,1),P(:,:,3));
    % W = huilier(a,b,c);
    
    save(file,'X','W','T');
end

% function a = vectorangles(x,y)
% % W. Kahan, 2006. "How Futile are Mindless Assessments of Roundoff in Floating-Point Computation?"
% % https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf    
% 
%     x = x./vecnorm(x,2,2);
%     y = y./vecnorm(y,2,2);
%     a = 2 * atan2(vecnorm(x - y,2,2) , vecnorm(x + y,2,2));
% end
% 
% function A = huilier(a,b,c)
% % https://en.wikipedia.org/wiki/Spherical_trigonometry#Area_and_spherical_excess
%     a = 0.5*a;
%     b = 0.5*b;
%     c = 0.5*c;
%     s = 0.5*(a + b + c);
%     A = 4*atan(sqrt(tan(s).*tan(s-a).*tan(s-b).*tan(s-c)));
% end
