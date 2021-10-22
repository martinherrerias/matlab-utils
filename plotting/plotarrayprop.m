function varargout = plotarrayprop(Trck,varargin)
% PLOTARRAYPROP(TRCK,P) - create a 2D colormap representation of physical property P over the 
%   tracker array defined by structure TRCK. P in this case must be a [Ntr/Nat,Nm/1] array, i.e.
%   have one row for each [analysed] mount, and optionally one column for each mount module.
%
% PLOTARRAYPROP(TRCK,ARRDEF,E) - create a 2D colormap representation of electrical property E over 
%   the tracker array defined by structure TRCK and connected according to NDMAP ARRDEF. 
%   E in this case must be a [Ni,Ns/1,Nm/1] array, i.e. have one row for each MPPT, and optionally
%   one column for each string, and one plane for each module.
%
% NOTES:
%   PLOTARRAYPROP will bin the  values of E/P, and plot all modules in the array as colored patches
%   in the current figure, with color according to their corresponding value.
%
%   NaN and -Inf/Inf values are allowed, and will not be considered for setting colorscale-limits.
%   NaN-valued modules will be plotted as 60% gray, -Inf and Inf will be included in the first and
%   last bins, unless explicit EDGES are used, in that case they might be treated as NaNs.
%
% INPUT:
%   TRCK - structure with fields 'type','centers', 'geom', and depending on TRCK.type: 'tilt',
%       'azimuth','slope','tracklimits', etc. (see MOUNTROTATIONS). Field 'analysedtrackers'  
%       will be used if size(X,1) = Na < Ntr (see below).
%       TRCK.geom: pvArea object representing a single mount/tracker. TRCKAREAS.elements(k).border 
%       must return a POLYGON object representing the outline of module k (1:Nm) in every mount.
%
%   P: [Ntr/Nat,Nm/1] numerical array of physical propperty values for each mount [and module]
%
%   ARRDEF: NDMAP object containing an array-definition with ARRDEF.psize == [Ntr/Nat,Np] and
%       ARRDEF.esize == [Ni,Ns,Nm].
%
%   E: [Ni,Ns/1,Nm/1] numerical array of electrical propperty values for each inverter [and string]
%           [and module]. Requires ARRDEF.
%
% PLOTARRAYPROP(...,EDGES,CMAP,HF)
%   'edges', EDGES : specify a number of bins or a vector of bin-edges for E/P.
%   'colormap', CMAP : use a non-standard colormap name (e.g. 'spring','hot',...) or provide an 
%       Nb·3 list of RGB colors for the respective bins; and finally plot the patches in figure HF instead of the
%   default clf().
%
% EXAMPLES:
%   [M,A] = samplesystem('0a','landscape',true);
%   plotarrayprop(M,M.centers(3,:)');   % plot mount height
%   plotarrayprop(M,A,(1:A.esize(1))'); % identify mounts by inverter
%
% See also: MOUNTROTATIONS, PVAREA, PLOTTRACKERARRAY, PLOTARRAYDEF

    opt.colormap = 'jet';
    opt.nocolor = [1,1,1]*0.6;
    opt.edges = 20;
    opt.axes = [];
    [opt,remargs] = getpairedoptions(varargin,opt);
    if isempty(opt.axes), AX = gca(); else, AX = opt.axes; end
    
    CMAPS = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'};
    
    Ntr = size(Trck.centers,2);
    if isfield(Trck,'analysedtrackers'), Na = numel(Trck.analysedtrackers); else, Na = Ntr; end
    Np = numel(Trck.geom.elements);

    switch numel(remargs)
      case 1
        x = varargin{1};
      case 2
      % Convert electrical to physical coordinates
        AD = varargin{1};
        x = varargin{2};         
        assert(isa(AD,'NDmap'),'Expecting NDmap object');
        assert(AD.psize(2) == Np && any(AD.psize(1) == [Na,Ntr]),'Invalid ARRDEF');

        [Nm,Ns,Ni] = size(AD.PIDX);
        if size(x,2) == 1 && Ns > 1, x = repmat(x,1,Ns,1); end
        if size(x,3) == 1 && Nm > 1, x = repmat(x,1,1,Nm); end
        assert(isequal(size(x),[Ni,Ns,Nm]),'plotarrayprop:edim','Expecting [Ni,Ns/1,Nm/1] array');
        
        y = NaN(AD.psize);
        x = permute(x,[3,2,1]);
        y(AD.pidxT(:)) = x(AD.PIDX > 0);
        x = y;

      otherwise, error('Unrecognized arguments')
    end
    
    if size(x,1) == Na && Na < Ntr
        x(Trck.analysedtrackers,:) = x;
        x(setdiff(1:Ntr,Trck.analysedtrackers),:) = NaN;
    end
    assert(size(x,1) == Ntr,'plotarrayprop:pdim','x must have Na or Ntr rows');
    
    if size(x,2) == 1 && Np > 1, x = repmat(x,1,Np); end
    assert(size(x,2) == Np,'plotarrayprop:pdim','x must have 1 or Nm columns');
    
    inbounds = ~(isnan(x) | isinf(x));
    assert(any(inbounds,'all'),'Nothing to plot');
    
    if isscalar(opt.edges)
        Nb = opt.edges;
        assert(Nb > 2,'plotarrayprop:nedges','Number of bins must be greater than 1');
        
        edges = linspace(min(x(inbounds)),max(x(inbounds)),Nb+1);
        edges(end) = edges(end) + eps(edges(end));
        if any(x == -Inf), edges(1) = -Inf; end
        if any(x == Inf), edges(end) = Inf; end
    else
        edges = opt.edges; 
    end
    assert(isvector(edges) && issorted(edges),'plotarrayprop:edges','EDGES is not a sorted vector');
    
    Nb = numel(edges) - 1;
    [~,x] = histc(x(:),edges);
    
    if ischar(opt.colormap)
        cmap = lower(opt.colormap);
        assert(any(strcmp(cmap,CMAPS)),'plotarrayprop:charcmap','Unrecognized color-map name');
        cmap = eval(sprintf('%s(%d);',cmap,Nb));
    else
        cmap = opt.colormap; 
    end
    assert(size(cmap,2) == 3 && size(cmap,1) == Nb,'plotarrayprop:sizecmap',...
        'CMAP must be a string (color-map name) or a (numel(EDGES) - 1) x 3 color list');
    
    cf = [opt.nocolor;cmap];
    cf = cf(x+1,:);
    
    R = mountrotations(Trck,90,90);
    if size(R,3) == 1, R = repmat(R,[1,1,Ntr]); end
    
    p0 = [Trck.geom.elements.border]; 
    p0 = repmat(p0,Ntr,1); % Ntr·Nm array of module-polygons
    
    for j = 1:Ntr
        p0(j,:) = polyrotate(p0(j,:),flatrotation(R(:,:,j)));
        p0(j,:) = polytranslate(p0(j,:),Trck.centers(1:2,j));
    end
			
    [V,~,F] = poly2vef(p0(:),true);

    H = patch(AX,'Faces',F,'Vertices',V,'FaceColor','flat','FaceVertexCData',cf,'EdgeColor','none');
    axis(AX,'equal');
    box(AX,'off');
    
    fmt = sprintf('%%0.%df',round(-log10(min(diff(edges))/10)));
    
    %set(gcf,'Name','Array Definition'); set(gcf,'Numbertitle','off');
    set(AX,'Xdir','reverse') % Equator (+y) facing down
    set(AX,'Ydir','reverse')
    colormap(AX,cmap);
    caxis(AX,[edges(1),edges(end)]);
    colorbar(AX,'ticks',edges,'ticklabels',cellstr(num2str(edges',fmt)));
    % colorbar();
    
    if nargout > 0, varargout{1} = H; end
    % % Print
    % view(180,90)
    % rez = 600; 
    % figpos = getpixelposition(h);
    % resolution = get(0,'ScreenPixelsPerInch');
    % set(h,'paperunits','inches','papersize',figpos(3:4)/resolution);...
    % set(h,'paperposition',[0 0 figpos(3:4)/resolution]);
    % print(h,'arrdef.png','-dpng',['-r',num2str(rez)],'-opengl')
end

function Q = flatrotation(R)
    x = R(1:2,1); x = x/hypot(x(1),x(2));
    y = R(1:2,2);
    y = y - y'*x*x; y = y/hypot(y(1),y(2));
    Q = [x,y];
end
