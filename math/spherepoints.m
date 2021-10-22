function [X,d] = spherepoints(N,varargin)
% X = SPHEREPOINTS(N) - Distribute N points evenly* on a unit sphere
%
% (*) Points are distributed starting with a parametric spiral distribution,
%   (elevation angle phi, azimuth angle alpha):
%
%     phi = asin(1-2·(0.5:N)/N) - from forcing spherical cap 2·pi·(1-sin(phi)) == n/N·4·pi
%     t² = 8/sqrt(3)·pi/N - from division of 4·pi sphere into equilateral triangles sqrt(3)/2·t²
%     alpha = cumsum(t/cos(phi)) - t ~ cos(phi)·d(alpha)/dn
%
%   Uniform random noise of +/- opt.NOISE·t is added on elevation and azimuth. Noise values of up to
%   ~1/4 tend to increase uniformity of distribution, but require more iterations to converge.
%
%   The distribution is then iteratively corrected applying a 'force' proportional to 1/r² from each
%   point's K nearest neighbors, with K = 6·sum(1:m) starting with m = opt.M_START and increasing 
%   the 'reach' (m++) up to K >= opt.N_MAX neighbors. For every given number of neighbors K, several
%   iterations can be performed, the stop criteria being that the standard deviation of the minimum
%   distance between neighboring points no longer decreases.
%
% X = SPHEREPOINTS(N,'symmetric',true) - return points over a half dome only
%       FUTURE: it'd be really fancy if you could enforce ANY symmetry group passed by the user:
%       replicate a random point > build voronoi cells > filter points in a single cell > copy to
%       all other > redistribute (only for points in control cell) > copy, redistribute, copy, ...
%
% X = SPHEREPOINTS(N,'regular',true) - for particular choices of N, returns points distributed
%   analytically:
%
%       N = 2 returns 2 poles
%       N = {4,6,8,12,20} returns the vertices of platonic solids
%       N = 32 returns the vertices of a Pentakis dodecahedron
%       N = {42, 162, 642, 2562, 10242,..} = 2 + 10·4^Z, for Z in {1,2,3,..} starts with an 
%           icosahedron, and recursively inserts vertices at each edge midpoint, renormalizing
%           at each step to bring them to the surface of the sphere.
%
% SEE ALSO: VORONOISPHERE
        
    opt.N_MAX = 330;     % max. number of neighbors to consider
    opt.NOISE = 1/16;    % random noise: improves uniformity, takes longer
    opt.M_START = 4;     % K = 6·sum(1:m) = 3·(m+1)·m for m = M_START
    opt.SYMMETRIC = false;
    opt.PLOT = false;
    opt.REGULAR = false;
    
    [opt,remargs] = getpairedoptions(varargin,opt);
    assert(isempty(remargs),'Bad arguments');
    
    if opt.REGULAR, X = regular(N); d = NaN; return; end
    
    if opt.SYMMETRIC
        t = sqrt((4*pi/sqrt(3))/N);         % fill half sphere only
        phi = asin(1-(0.5:N)/N);  
        h = sqrt(3)/2*t;
        phi = round((phi-h/2)/h)*h + h/2;   % break spiral into planes h/2, 3h/2, 5h/2, ...
    else
        t = sqrt((8*pi/sqrt(3))/N);
        phi = asin(1-2*(0.5:N)/N);
    end
    alpha = cumsum(t./cos(phi));

    phi = phi + t*opt.NOISE*(2*rand(1,N)-1);      % add random perturbation
    alpha = alpha + t*opt.NOISE*(2*rand(1,N)-1);
        
    X = sph2cartV(alpha(:)*180/pi,phi(:)*180/pi);

    % X = randn(N,3);
    % X = X./rssq(X,2);

    SSD = [];

    for m = opt.M_START:N/6  % will break sooner, at floor(sqrt(1+4/3·N_MAX)-1)/2
        K = 3*(m+1)*m;
        if K > 6*N || K > opt.N_MAX, break;
        elseif K > N-1, K = N-1;
        end

        if ~opt.SYMMETRIC
            Y = X;
        else
        % Consider XY-plane-mirrored points as well
            X(:,3) = abs(X(:,3));
            Y = [X;X];
            Y(N+1:end,3) = -Y(N+1:end,3);
        end
        
        cIdx = knnsearch(Y,X,'K',K+1); 

        D = zeros(N,3,K);
        sD = Inf;
        e = Inf;

        while any(e > t/1000)

            for j = 1:K, D(:,:,j) = X - Y(cIdx(:,j+1),:); end
            r = rssq(D,2);
            f = sum(D./r.^3,3);
            d = min(r,[],3);
            md = mean(d);
            last_sD = sD;
            sD = std(d);

            f = f - dot(X,f,2).*X;      % remove radial component
            e = abs(mean(d)-d);         % error metric is difference with mean min-distance
            f = 0.5*e.*f./rssq(f,2);    % restrict movement to half the error

            SSD = cat(1,SSD,mean(e)/md);
            
            if sD >= last_sD, break; end

            X = X + f;
            if opt.SYMMETRIC, X(:,3) = abs(X(:,3)); end 
            X = X./rssq(X,2); 
        end
    end
    
    if opt.PLOT
        msg{1} = sprintf('exiting with K = %d',K);
        msg{2} = sprintf('MAE = %0.2f%%, STD = %0.2f%%',SSD(end)*100,sD/mean(d)*100);

        figure(gcf());
        clf(); set(gcf(),'position',get(groot,'Screensize').*[0.1,0.1,0.8,0.5]);

        subplot(1,3,2); cla(); plot(SSD*100);
        set(gca,'YScale','log')
        yticks(1:9);
        axis([1 numel(SSD) 1 10]);
        xlabel('iterations');
        ylabel('min. point-distance MAE (%)');
        text(numel(SSD)-1,9,msg,'HorizontalAlignment','Right','VerticalAlignment','Top');

        subplot(1,3,3); cla(); 
        hist(100*e/mean(d)); 
        xlabel('min. point-distance variation (%)');

        subplot(1,3,1);
        cla(); hold on; axis equal;
        [Xs,Ys,Zs] = sphere(40);
        surf(Xs,Ys,Zs,'FaceColor',[1 1 1],'FaceAlpha',0.6,'EdgeAlpha',0.0);
        scatter3(X(:,1),X(:,2),X(:,3),1+min(9,2000/size(X,1)),'filled');
        view(60,30);
    end
end

function n = regular(N)
    switch N
        case 2, n = [0,0,1;0,0,-1];                             % 2 poles
        case 4, n = multisensor(0,0,120,[0 120 240])';          % tetrahedral
        case 6, n = [0,0,1;0,1,0;1,0,0;-1,0,0;0,-1,0;0,0,-1];   % octahedral vertices
        case 8, n = multisensor(45,0:90:270,135,0:90:270)';     % cube vertices
        case 12                                                 % icosahedral vertices
            a = atand(2);
            w = (0:4)*72;
            n = multisensor(0,0,a,w,180-a,w+36,180,0)';         
        case 20                                                 % dodecahedral vertices
            n = regular(12);
            TR = triangulation(convhull(n),n);
            n = TR.circumcenter;
            n = n./rssq(n,2);
        case 32
            n = [regular(12);regular(20)];
        otherwise
            M = log((N-2)/10)/log(4);
            if N >= 42 && mod(M,1) == 0                         % recursively divided icosahedron
                n = regular(12);
                for j = 1:M
                    T = convhull(n);
                    TR = triangulation(T,n);
                    E = TR.edges;
                    m = (n(E(:,1),:) + n(E(:,2),:))/2;  % add edge middle-points
                    n = [n;m./rssq(m,2)]; %#ok 
                end
                return
            end
        error('Unknown regular distribution');
    end
end