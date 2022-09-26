function varargout = plotcontrols(types,labels,ranges,values,updatefcn,varargin)
% H = PLOTCONTROLS(TYPE,LABELS,RANGES,VALUES,UPDATEFCN,..)
%   Draw a series of UICONTROLS (currently only sliders or popup-menus), with given LABELS (*),
%   RANGES and starting VALUES. Their UI values are then linked to a persistent cell-array of
%   values X, whose change triggers a function call of the form UPDATEFCN(X{:},SRC,EVENT). The
%   idea is to collect the callbacks of several (related) controls into a single function, that
%   sees the updated values of all parameters on each call.
%
%   Run PLOTCONTROLS() and examine the code for PLOTCONTROLS.TEST for an example.
%
%   TYPES - cellstring of UICONTROL styles (or short-keys: 's' == 'slider','m' = 'popupmenu')
%
%   LABELS - cellstring of control labels. They are placed as an ANNOTATION objects on top of
%       each control, and can be addressed from each control H(j) as H(j).UserData.txtbox.
%       For sliders, labels are updated automatically to reflect the runtime value, i.e.
%       "label: VALUE".
%
%   RANGES - cell array, with contents matching to each control type:
%       + sliders: 2-4 vectors [min,max,[minorstep majorstep]] 
%       + popupmenus: cellstring of menu items
%
%   VALUES - cell array of starting values, matching each control type.
%
%   ..,'parent',gcf() - Where to place UICONTROLS
%   ..,'position',[0,0,1,1] - Extents (units of parent.Position) of the set of controls
%   ..,'-skipupdate' - skip first call to UPDATEFCN
%   ..,'-smooth' - call UPDATEFCN for real-time change in slider controls, not just when released.
%   ..,'-listvalues' - Pass popup-menu index (1,2,3,..) instead of actual pick {'A','B',..}
%                      to callback UPDATEFCN.
%   ..,'fontsize',10 - Label font size
%   ..,'format','%0.3g' - Number format to display on slider label 
    
    if nargin == 0, test(); return; end
    narginchk(6,14);
    
    [opt,varargin] = getflagoptions(varargin,{'-skipupdate','-listvalues','-smooth'});
    opt.fontsize = 10;
    opt.format = '%0.3g';
    opt.position = [0,0,1,1];
    opt.parent = gcf();
    opt = getpairedoptions(varargin,opt,'restchk');
    
    assert(ishandle(opt.parent),'Expecting figure handle');
    assert(isa(updatefcn,'function_handle'),'Expecting function handle');
    
    assert(all(cellfun(@iscell,{types,labels,ranges,values})),'Expecting 4 cell arrays as input');
    [types,labels,ranges,values] = compatiblesize(types,labels,ranges,values);

    X = values;
    m = numel(labels);

    % set(ax,'visible',false);
    Pos = getpositions(m,opt);
    
    H = arrayfun(@(j) drawcontrol(types{j},ranges{j},labels{j},values{j},Pos{j},j),1:m);

    if ~opt.skipupdate,updatefcn(X{:}); end
    
    if nargout > 0, varargout{1} = H; end
    


    function L = drawcontrol(type,range,label,value,pos,j)
    % Actual UICONTROL calls and callback function assignment. Useful parameters are stored
    % as fields in L.UserData.
    
        original_lbl = label;
        switch lower(type)
        case {'slider','s'}
            
            assert(isnumeric(range) && all(isfinite(range)) && isreal(range) && ...
                any(numel(range) == [2,4]),'Range must be a 2 or 4-vector of real values');
            assert(range(2) >= range(1),'Empty range');
            assert(value >= range(1) && value <= range(2),'Value outside provided range');
            if numel(range) > 2
                assert( range(3) <= range(4) && range(4) < range(2)-range(1),...
                    'Range [m,M,s,S] must comply with s < S < (M-m)');
            end
            
            L = uicontrol(opt.parent,'style','slider','min',range(1),'max',range(2),'value',value,...
                      'units','normalized','position',pos(1,:));
            if ~isempty(opt.format)
                label = strjoin({label,num2str(value,opt.format)},': ');
            end
            if numel(range) == 4
               L.SliderStep = range(3:4)./(range(2)-range(1)); 
            end
            
            L.UserData.smooth = opt.smooth;
            addlistener(L,'Value','PostSet',@callback); % for smooth value updates
            
        case {'popupmenu','dropdown','m'}
            
            assert(iscellstr(range) || isstring(range),'Expecting cellstr list for pop-menu');
            
            if ischar(value)
                [~,value] = ismember(value,range);
                if opt.listvalues, X{j} = value; end
            else
                if ~opt.listvalues, X{j} = range(value); end 
            end
            L = uicontrol(opt.parent,'style','popupmenu','string',range,'value',value,...
                'units','normalized','position',pos(1,:),'fontsize',opt.fontsize);

        case {'text'}
                        
            if ischar(value)
                [~,value] = ismember(value,range);
                if opt.listvalues, X{j} = value; end
            else
                if ~opt.listvalues, X{j} = range(value); end 
            end
            L = uicontrol(opt.parent,'style','text','string',value,...
                'units','normalized','position',pos(1,:),'fontsize',opt.fontsize,...
                'HorizontalAlignment','left');
                
        % case {'togglebutton','checkbox','radiobutton','edit','listbox'}
        otherwise
            
            error('plotcontrol currently only works with sliders and pop-up menus');
        end
        L.UserData.element = j;
        L.UserData.label = original_lbl;
        L.UserData.range = range;
        L.Callback = @callback;
        
        if ~strcmp(L.Style,'text')
            L.UserData.txtbox = annotation(opt.parent,'textbox','units','normalized','position',pos(2,:),...
                'margin',0,'EdgeColor','none','String',label,'verticalalignment','bottom',...
                'fontsize',opt.fontsize);
        end
    end

    function callback(L,event)

        switch event.EventName
            case 'Action', dostuff = true;   % menu change, slider release          
            case 'PostSet' 
                L = event.AffectedObject;
                dostuff = L.UserData.smooth; % smooth slider change
            otherwise
                error('well, this was unexpected');
        end
        k = L.UserData.element;
        v = L.Value;
        
        switch L.Style
        case {'popupmenu'}
            if opt.listvalues, X{k} = v;
            else
                X{k} = L.String{v};
            end
        case {'slider'}   
            if ~isempty(opt.format)
                L.UserData.txtbox.String = strjoin({L.UserData.label,num2str(v,opt.format)},': ');
            end
            if ~dostuff, return; end
            X{k} = v;
        case {'text'}
        otherwise
            error('you should not be here')
        end

        updatefcn(X{:},L,event);
    end
end

function Pos = getpositions(m,opt)

    if iscell(opt.position)
       Pos = opt.position;
       assert(numel(Pos) == m && all(cellfun(@(p) isequal(size(p),[4,2]),Pos)),...
           'Invalid position(s)');
       return;
    end
    validateattributes(opt.position,{'numeric'},{'real','vector','finite','numel',4});
        
    % p = frame.Position(3:4); 
    % x0 = frame.Position(1:2) + opt.position(1:2).*p;
    x0 = opt.position(1:2);
    p = opt.position(3:4)./[1,m];  % size of combined control+label
    
    opt.parent.Units = 'points';
    rows = p(2)*opt.parent.Position(4)/opt.fontsize; % text rows (height) for combined control+label
    opt.parent.Units = 'normalized';
    
    if rows > 3
        c = 1.5/rows;         % control height = 1 text row, units of p(2)
        t = (1 - 1.5/rows);   % tag height
    else
        c = 1/2; t = c;
    end
    
    Pos = arrayfun(@(j) [ x0 + p.*[0,j],p.*[1,c] ; x0 + p.*[0,j+c], p.*[1,t] ], m-1:-1:0,'unif',0);
end

function test()

    GUIfigure('plotcontrols','plotcontrols test','2:1'); clf();
    set(subplot(1,2,1),'visible','off');
    set(subplot(1,2,2),'visible','off'); axis equal;
    
    h = []; % stores graphics handle

    plotcontrols({'popupmenu','slider','s','s'},...
                 {'shape','red','green','blue'},...
                 {{'triangle','square','pentagon'},[0,100],[0,100],[0,100]},...
                 {'triangle',50,50,50},@update,'fontsize',14,'format','%0.0f%%',...
                 'position',subplot(1,2,1).Position,'-listvalues','-smooth');

    function update(shapeidx,r,g,b,~,~)
        n = shapeidx + 2;
        delete(h);
        h = polyplot(polygon(n),[r,g,b,100]/100,'k',subplot(1,2,2));
    end
end