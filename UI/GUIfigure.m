function fh = GUIfigure(varargin)
% H = GUIFIGURE(KEY,TITLE,POS) - to avoid duplicate figures poping-up all over the place, a list 
%   of figure figure-handles attached to character KEYS is kept in a persistent variable.
%   If the provided KEY is not in use (it never existed, or its handle has been deleted), 
%   H = figure() is used, and H stored. If there is a valid handle H associated to that KEY,
%   H = figure(H) will return it.
%
%   If TITLE is provided, the figure's title is set to TITLE, otherwise KEY is used.
%
%   By default, 'GUImain' figure's position [y0 y0 w h] takes the left side of the screen, minus
%   margins, and all others are squares on the right side.
%   If POS is provided AND the figure is not already being displayed, the new figure's position 
%   will be adjsuted accordingly:
%
%     + POS can be a 2×4 matrix, such that [x0 y0 w h] = [W,H]·POS, for the screen's [W,H]
%     + POS can be a 4-vector, [x0 y0 w h] or relative vector [x0/W y0/H w/W h/H]
%     + POS = [] skips any position setting (use MATLAB internal defaults)
%     + POS = 'A:B' adjust the aspect-ratio of the standard-sized figure.
%     + POS = 'key' for key in {main,fig,big,max} uses default named sizes. Just try them out.
%
% IDX = GUIFIGURE(..,'scr',[W,H]) - override the default SCR = get(groot,'Screensize')(3:4)
%
% IDX = GUIFIGURE(KEY,'-silent') - returns the numerical index assigned to KEY, without displaying
%   any figure window.
%
% S = GUIFIGURE('-all') - returns the internal (structure) index, S.KEY = IDX.
%
% See also: SPLITGUI, FIGURE, GCF, FINDOBJ

    MARGIN = 0.05;  % used for default position offsets
    HEAD = 0.059;   % difference between fig.position and outerposition / H
    N = 10;
    SCR = get(groot,'Screensize'); 
    SCR = min(SCR(3:4),[1920,1080]); % [W,H]
    
    POS.main = [MARGIN,0;0,MARGIN;0.5-2*MARGIN,0;0,1-2*MARGIN-HEAD]';
    % POS.fig = [0.5,0;0,2*MARGIN;0,(1-4*MARGIN-HEAD);0,(1-4*MARGIN-HEAD)]';
    POS.fig = [0.5,0;0,2*MARGIN;0,(0.8-4*MARGIN-HEAD);0,(0.8-4*MARGIN-HEAD)]';
    POS.big = [0.5,0;0,MARGIN;0,(1-2*MARGIN-HEAD);0,(1-2*MARGIN-HEAD)]';
    POS.max = [0,0;0,0;1,0;0,1]';
    
    persistent fighandles;
    if isempty(fighandles), fighandles = struct(); end
    
    % fighandles = getSimOption('fighandles');
    [opt,varargin] = getflagoptions(varargin,{'-all','-silent'});
    opt.scr = SCR;
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    if opt.all, fh = fighandles; return; end
    
    assert(~isempty(varargin) && ischar(varargin{1}),'Expecting char-array KEY');
    key = matlab.lang.makeValidName(varargin{1});

    if numel(varargin) < 2, figtitle = key; else, figtitle = varargin{2}; end
    if numel(varargin) < 3
        pos = opt.scr*POS.fig;
        defaultpos = true;
    else
        pos = varargin{3};
        defaultpos = false;
        
        if ischar(pos)
            pos = lower(pos);
            if isfield(POS,pos)
            % fig, main, big
                pos = opt.scr*POS.(pos);
            else
            % Read ratio (e.g. '3:4')
                r = str2double(strsplit(pos,':'));
                assert(numel(r) == 2 && all(isfinite(r) & r > 0),'');
                
                pos = opt.scr*POS.fig;
                A = prod(pos(3:4));
                r = r(1)/r(2);
                pos(3:4) = sqrt([A*r,A/r]);
                r = max([1,pos(3:4)./(opt.scr*(1-MARGIN-HEAD))]); % scale down if it doesn't fit
                pos(3:4) = pos(3:4)/r;
                pos(2) = (opt.scr(2)-pos(4))/2; % centered on y
                pos(1) = min(opt.scr(1)*0.75,opt.scr(1)*(1-MARGIN)-pos(3)); % ~center on right side
            end
        elseif isvector(pos) && numel(pos) == 4
            pos = pos(:)';
            if all(pos <= 1)
                pos = pos.*[opt.scr,opt.scr];
            end
        elseif ~isempty(pos)
            assert(isequal(size(pos),[2,4]),'Expecting 4-vector or 2×4 matrix')
            pos = opt.scr*pos;
        end
    end
    
    if isfield(fighandles,key), fh = fighandles.(key); else, fh = NaN; end
    
    if isgraphics(fh) && isa(fh,'matlab.ui.Figure')
        figureexisted = true;
        if ~opt.silent, fh = figure(fh); end
    else
        if ~opt.silent, fh = figure(); end
        fighandles.(key) = fh;
        figureexisted = false;
    end
    if opt.silent, return; end
    
    if ~isempty(pos) && (~figureexisted || ~defaultpos)
    % Reset position only with explicit argument, or on new figures
    
        if defaultpos
        % Apply a small offset to avoid complete overlap

            f = fieldnames(fighandles);
            n = find(strcmp(key,f));      % index of current figure within all registered
            x = mod(N-n,N)/(N-1);         % a sort of offset fraction 0..1, based on n

            % margin available from top-right corner of figure [x,y]
            offset = opt.scr*(1-MARGIN) - [pos(1)+pos(3),pos(2)+pos(4)];
            
            offset = 2*(x-0.5)*max([offset,0]);
            pos(1:2) = pos(1:2) + offset;
        end

        set(fh,'position',pos);
    end
    
    if ~isempty(figtitle)
        set(fh,'Name',figtitle); 
        set(fh,'Numbertitle','off');
    end
end