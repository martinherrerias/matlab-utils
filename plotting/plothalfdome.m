function varargout = plothalfdome(varargin)
% PLOTHALFDOME() - Plot a unitary spherical-axis half dome on current axes.
% PLOTHALFDOME(AX,..) - draw on axes AX
%
% ..,'R',1 - set radius
% ..,'rotation',[] - if specified, it will plot cardinal directions, and rotate XY axes accordingly
% ..,'view',[210,20] - sets view(AZ,EL), unless flag '-keepview' is set
% ..,'parallels',[5 10 15 30 45 60] - change location of reference parallels (elevation angles)
% ..,'meridians',0:30:150 - change location of reference meridians (azimuth angles)
%
% ..,'-full' - plot also lower hemisphere (changes default parallels to -60:30:60)
% ..,'-keepaxis' - don't delete cartesian axis (axis image off).
% ..,'-keepview' - don't rotate to default view(210,20)
%
% H = PLOTHALFDOME(..) - return an array of graphic objects, e.g. for later delete(H).
    
    opt.R = 1;
    opt.rotation = [];
    opt.view = [210,20];
    opt.parallels = [];
    opt.meridians = 0:30:150;
    
    if nargin > 0 && isa(varargin{1},'matlab.graphics.axis.Axes') && ishandle(varargin{1})
       opt.ax = varargin{1};
       varargin(1) = [];
    else
       opt.ax = gca(); 
    end
    [opt,~,def] = parseoptions(varargin,{'-keepaxis','-keepview','-full'},opt,'dealrest',2);

    if isempty(opt.parallels) && def.parallels
        if opt.full, opt.parallels = -60:30:60;
        else, opt.parallels = [5,10,15,30,45,60]; % 2.^(0:6); % 1 2 4 8 .. 64
        end
    end
    parsestruct(opt,{'parallels','meridians'},'numeric','size',[1,NaN],'finite','real');
    
    if def.full, opt.full = any(opt.parallels < 0); end
    
    hold(opt.ax,'on');
    
    if isempty(opt.rotation)
        axlbl = {'x','y'};
        axpos = [1,0 ; 0 1];
    elseif mod(opt.rotation,90) == 0
        axlbl = circshift({'E','N','W','S'},-opt.rotation/90);
        axlbl{1} = ['x, ' axlbl{1}];
        axlbl{2} = ['y, ' axlbl{2}];
        axpos = [1,0; 0 1; -1,0; 0,-1] ;
    else
        axlbl = {'x','y','E','N','W','S'};
        c = cosd(-opt.rotation);
        s = sind(-opt.rotation);
        axpos = [1,0; 0 1; c,s; -s,c; -c,-s; s,-c];
    end
    axpos(:,3) = 0;
    axlbl{end+1} = 'z';
    axpos(end+1,:) = [0,0,1];
    axpos = axpos.*1.2*opt.R;

    % Draw abs axis
    h{1} = quiver3(opt.ax,zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1],0,...
        'k-','HandleVisibility','off');
    h{2} = text(opt.ax,axpos(:,1),axpos(:,2),axpos(:,3),axlbl,...
                'VerticalAlignment','middle','HorizontalAlignment','center'); 
    
    circle = pi/180*[0:5:360,NaN]';
    if opt.full
        arc = pi/180*[-180:5:180,NaN]';
    else
        arc = pi/180*[0:5:180,NaN]';
    end
  
    % Draw "meridians"
    [x,y,z] = sph2cart((opt.meridians)*pi/180,arc,opt.R);
    h{3} = plot3(opt.ax,x,y,z,':','color',[1 1 1]*0.8,'HandleVisibility','off');
    
    % ... "parallels"
    [x,y,z] = sph2cart(circle,opt.parallels*pi/180,opt.R);
    [x,y,z] = compatiblesize(x,y,z);
    h{4} = plot3(opt.ax,x,y,z,':','color',[1 1 1]*0.8,'HandleVisibility','off');
    
    % and XY plane
    [x,y,z] = sph2cart(circle(1:end-1),0,opt.R);
    h{5} = patch(opt.ax,x,y,z,'facecolor',[1 1 1]*0.8,'facealpha',0.2,'edgecolor',[1 1 1]*0.8,...
        'HandleVisibility','off');
        
    axis square; axis equal;
    if ~opt.keepaxis, axis(opt.ax,'image','off'); end
    if ~opt.keepview, view(opt.ax,opt.view); end
    
    if nargout > 0, varargout{1} = h; end
end