function varargout = densityplot(x,y,varargin)
% DENSITYPLOT(X,Y,S) - scatter plot SCATTER(X,Y,S,C) where colors C are set proportional to the 
%   kernel density of the distribution of points X,Y. Specifically:
% 
%   C = CMAP(d), where CMAP is a custom mapping function (see below), and:
%   d = INTERPN(Gx,Gy,D,x,y), D = KSDENSITY([x,y],[Gx,Gy],'bandwidth',BANDWIDTH)
%
% [H,D,OPT] = DENSITYPLOT(X,Y,S,..) - returns the graphics handle H, a vector of estimated kernel
%   densities D, and the structure of options OPT.
% 
%   The grid of density estimation points Gx,Gy is controlled by argument-value pairs 
%   ..,'resolution',N,.. and ..,'limits',[XMIN,XMAX,YMIN,YMAX]. 
%
%   ..,'cmap',F.. Use color mapping function F = @(x), default is to use ordered density rank,
%   i.e: @tiedrank. Alternatives could be linear or log. maps, e.g. @(x) x, or @(x) log(1+x).
%
%   ..,'maxpts',MAXPTS (default 10000) - If the number of points exceeds MAXPTS, only a random
%   subset of approx MAXPTS points will be used to estimate the set's KSDENSITY.
%
%   ..,'bandwidth',BANDWIDTH - overrides the default BANDWIDTH = 0.5/N
%
%   ..,'-notfilled' flag overrides SCATTER's 'filled', i.e. 'filled' becomes the default.
%
%   ..,'ax',H plot on axis H instead of the default gca()
%
%   ..,'Name',Value - any remaining arguments will be passed to SCATTER function
%
% See also: SCATTER, KSDENSITY

    % Parsing
    narginchk(2,Inf);
    assert(isequal(size(x),size(y)),'X and Y must have the same size');
    
    if nargin < 3 || ~isnumeric(varargin{1}), s = 10;
    else
        s = varargin{1};
        varargin = varargin(2:end);
    end
    if ~isempty(varargin) && isnumeric(varargin{1})
       warning('Overriding color information');
       varargin = varargin(2:end);
    end
    [opt,varargin] = getflagoptions(varargin,{'filled','notfilled'});
    opt.resolution = 100;
    opt.maxpts = 10000;
    opt.bandwidth = [];
    opt.limits = [];
    opt.cmap = @tiedrank;
    opt.ax = [];
    [opt,varargin] = getpairedoptions(varargin,opt);
    if ~opt.notfilled, varargin{end+1} = 'filled'; end
    if isempty(opt.ax), opt.ax = gca(); end
    
    if nargout > 1, varargout{2} = zeros(size(x),'like',x); end
    
    idx = isfinite(x) & isfinite(y);
    x = double(x(idx));
    y = double(y(idx));
    
    if isscalar(opt.resolution), opt.resolution(2) = opt.resolution; end
    if isempty(opt.limits)
        opt.limits = [min(x,[],'omitnan'),max(x,[],'omitnan'),...
                      min(y,[],'omitnan'),max(y,[],'omitnan')];
    end
    if isempty(opt.bandwidth)
        opt.bandwidth(1) = 0.5*diff(opt.limits(1:2))/opt.resolution(1); 
        opt.bandwidth(2) = 0.5*diff(opt.limits(3:4))/opt.resolution(2); 
    end
    
    if isempty(opt.cmap), opt.cmap = @tiedrank; end
    
    if numel(x) < 3 || ~any(opt.limits([2,4]) - opt.limits([1,3]) > eps(1))
        warning('Cannot estimate point density');
        varargout{2}(:) = NaN;
        h = scatter(opt.ax,x,y,s,varargin{:});
    else
    % Actual work
        gx = linspace(opt.limits(1),opt.limits(2),opt.resolution(1)+1);
        gy = linspace(opt.limits(3),opt.limits(4),opt.resolution(2)+1);
        [gx,gy] = ndgrid(gx,gy);

        if numel(x) > opt.maxpts
            i = rand(numel(x),1) < opt.maxpts/numel(x);
        else
            i = true(numel(x),1);
        end
        gd = ksdensity([x(i),y(i)],[gx(:),gy(:)],'bandwidth',opt.bandwidth);
        d = interpn(gx,gy,reshape(gd,size(gx)),x,y);

        h = scatter(opt.ax,x,y,s,opt.cmap(d),varargin{:});
        if nargout > 1, varargout{2}(idx) = d; end
    end
    
    if nargout > 0, varargout{1} = h; end
    if nargout > 2, varargout{3} = opt; end
end