function varargout = plotarrows(varargin)
% PLOTARROWS() - Draw left- and right-arrow buttons, linked to keyboard keys, to step the current
%  axis left or right.
%
% ..,AX ,..  link to a specific axis
% ..,POS ,.. for POS in {'northeast','southwest',.. } draw on a particular corner
% ..,S ,..   change size (default is 0.05 of AX.Parent.Position)
    
    narginchk(0,3);
    ax = gca(); s = 0.05; pos = 'southeast';
    
    used = false(size(varargin));
    for j = 1:numel(varargin)
       if ~used(1) && ishandle(varargin{j}) && isa(varargin{j},'matlab.graphics.axis.Axes')
           ax = varargin{j}; used(1) = true;
       elseif ~used(2) && ischar(varargin{j}), pos = varargin{j}; used(2) = true;
       elseif ~used(3) && isnumeric(s) && isscalar(s), s = varargin{j}; used(3) = true;
       else
           error('Unrecognized arguments');
       end
    end
    hfig = ax.Parent;
    
    switch lower(pos)
        case 'southeast', pos = {[1-2*s 0 s s], [1-s 0 s s]};
        case 'southwest', pos = {[0 0 s s], [s 0 s s]};
        case 'northeast', pos = {[1-2*s 1-s s s], [1-s 1-s s s]};
        case 'northwest', pos = {[0 1-s s s], [s 1-s s s]}; 
    end
    
    h = uicontrol('Parent',hfig,'Style','pushbutton','String','<','UserData',ax,...
            'Units','normalized','Position',pos{1},'Visible','on','Callback',@goleft);
    h(2) = uicontrol('Parent',hfig,'Style','pushbutton','String','>','UserData',ax,...
            'Units','normalized','Position',pos{2},'Visible','on','Callback',@goright);
        
    if nargout > 0, varargout{1} = h; end
end

function goleft(src,~)
    ax = src.UserData;
    L = xlim(ax);
    xlim(ax,L - diff(L));
end

function goright(src,~)
    ax = src.UserData;
    L = xlim(ax);
    xlim(ax,L + diff(L))
end