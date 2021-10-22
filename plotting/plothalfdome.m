function varargout = plothalfdome(varargin)
% PLOTHALFDOME() - Plot a unitary spherical-axis dome on current axes.
% PLOTHALFDOME(R,LAT) - use radius R, and optionally label the axes according to convention
%   +Y = equator, i.e. +Y = N if LAT < 0, S otherwise.
%
% PLOTHALFDOME(..,'keepaxis') - don't delete cartesian axis (axis image off).
% PLOTHALFDOME(..,'keepview') - don't rotate to default view(210,20)
%
% H = PLOTHALFDOME(..) - return an array of graphic objects, e.g. for later delete(H).

    DEFVIEW = [210,20];
    ELEV = [5,10,15,30,45,60]; % 2.^(0:6); % 1 2 4 8 .. 64
    
    if nargin > 0 && ishandle(varargin{1}) && isa(varargin{1},'matlab.graphics.axis.Axes')
       ax = varargin{1};
       varargin(1) = [];
    else
       ax = gca(); 
    end
    hold(ax,'on');

    [opt,args] = getflagoptions(varargin,{'-keepaxis','-keepview'});
    assert(numel(args) < 3,'Unrecognized arguments');
    if numel(args) < 1, R = 1; else, R = args{1}; end
    if numel(args) < 2, rotation = []; else, rotation = args{2}; end
    
    if isempty(rotation)
        axlbl = {'x','y'};
        axpos = [1,0 ; 0 1];
    elseif mod(rotation,90) == 0
        axlbl = circshift({'E','N','W','S'},-rotation/90);
        axlbl{1} = ['x, ' axlbl{1}];
        axlbl{2} = ['y, ' axlbl{2}];
        axpos = [1,0; 0 1; -1,0; 0,-1] ;
    else
        axlbl = {'x','y','E','N','W','S'};
        c = cosd(-rotation);
        s = sind(-rotation);
        axpos = [1,0; 0 1; c,s; -s,c; -c,-s; s,-c];
    end
    axpos(:,3) = 0;
    axlbl{end+1} = 'z';
    axpos(end+1,:) = [0,0,1];
    axpos = axpos.*1.2*R;

    % Draw abs axis
    h{1} = quiver3(ax,zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1],0,...
        'k-');
    h{2} = text(ax,axpos(:,1),axpos(:,2),axpos(:,3),axlbl,...
                'VerticalAlignment','middle','HorizontalAlignment','center'); 
    
    circle = pi/180*[0:5:360,NaN]';
    arc = pi/180*[0:5:180,NaN]';
  
    % Draw "meridians"
    [x,y,z] = sph2cart((0:30:150)*pi/180,arc,R);
    h{3} = plot3(ax,x,y,z,':','color',[1 1 1]*0.8);
    
    % ... "parallels"
    [x,y,z] = sph2cart(circle,ELEV*pi/180,R);
    [x,y,z] = compatiblesize(x,y,z);
    h{4} = plot3(ax,x,y,z,':','color',[1 1 1]*0.8);
    
    % and XY plane
    [x,y,z] = sph2cart(circle,0,R);
    h{5} = patch(ax,x,y,z,'facecolor',[1 1 1]*0.8,'facealpha',0.2,'edgecolor',[1 1 1]*0.8);

    axis square; axis equal;
    if ~opt.keepaxis, axis(ax,'image','off'); end
    if ~opt.keepview, view(ax,DEFVIEW); end
    
    if nargout > 0, varargout{1} = h; end
end