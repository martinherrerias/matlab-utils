classdef plotter < handle
% PLOTTER - handle class that allows UI switching between a figure window and a waitbar.
%   The figure has control buttons  {Pause, Next, Fwd, Stop Plotting}, and the waitbar (if it
%   is actually running on UI mode), has a {Plot} button.
        
    properties
        name
        pausing      % if plotting, stop at every step - flag only, needs external control
    end 
    properties (Transient = true, GetAccess = public, SetAccess = protected)
        isUI        % Disables plotting mode, sets OPTWAITBAR's isUI
        plotting    % false = optwaitbar mode, true = figure mode
        fighandle   % handle of waitbar or plot
    end
    properties (Hidden = true)   
        btns        % structure of button handles {'sw2plot','sw2wait','pause','next','fwd'}
        queue       % cell-array of (string) instructions to be executed upon refresh(obj)
    end
    methods
        function obj = plotter(plotting,name,varargin)
        % OBJ = PLOTTER(PLOTING,VARARGIN) -Create [and draw] one of the two instances of PLOTTER:
        %
        %   a. If ~PLOTTING, create an OPTWAITBAR object, and if OBJ.isUI, adds an option button 
        %      that allows switching to plotting mode (see SWITCHTOPLOTTER).
        %      Makes a call to DRAWWAITBAR(VARARGIN{:}).
        %   b. If PLOTTING, draw GUIFIGURE('plotter',..) and add control buttons.
        %      Makes a call to DRAWPLOTTER(VARARGIN{:}).
        
            narginchk(2,Inf);
            plotter.clearexisting();
        
            UI = runningfromUI();
            if nargin < 1 || isempty(plotting), plotting = UI; end
            if (plotting == 1) && ~UI
                warning('plotter:UI','Plotting is disabled on non-UI mode');
                plotting = false; 
            end
            
            % Initialize object
            obj.name = name;
            obj.isUI = UI;
            obj.fighandle = NaN;
            obj.plotting = plotting;
            obj.queue = '';
            obj.pausing = false;           
            obj.btns = struct('sw2plot',NaN,'sw2wait',NaN,'pause',NaN,'next',NaN,'fwd',NaN);
 
            if (plotting == 1)
                drawplotter(obj,varargin{:}), 
            elseif (plotting == 0)
                drawwaitbar(obj,0,'Initializing...','name',obj.name);
            else
                % backdoor for subclass constructors to delay drawing
            end
        end  
        
        function drawwaitbar(obj,varargin)
            assert(~obj.plotting,'obj.drawwaitbar requires that ~obj.plotting');
            
            wb = optwaitbar(varargin{:});
            obj.fighandle = wb;
            obj.plotting = false;
            obj.pausing = false;

            if obj.isUI
                set(wb.ID,'CloseRequestFcn',@(~,~) obj.close);

                % Add switch-to-plotter button
                set(wb.ID,'position',get(wb.ID,'position').*[1 1 0 1] + [0 0 350 0]);
                obj.btns.sw2plot = uicontrol(wb.ID,'Style','pushbutton','String','Plot',...
                    'Units','normalized','Position',[0.8 0.25 0.15 0.5],'Visible','on',...
                    'Callback',@(~,~) obj.switchtoplotter());
            end
        end

        function drawplotter(obj,varargin)
            assert(obj.plotting,'obj.drawplotter requires that obj.plotting');
            obj.pausing = true; % Stop at first point, let the user click Fwd to disable

            plotter.clearexisting();
            
            h = GUIfigure('plotter',varargin{:}); clf(h);
            obj.fighandle = h;
            set(h,'CloseRequestFcn',@(~,~) obj.close);

            % Add 'Pause' button (set obj.pausing = true)
            obj.btns.pause = uicontrol(h,'Style','pushbutton','String','Pause','Units','normalized',...
                'Position',[0.64 0.95 0.1 0.04],'Visible','on','Callback',@(~,~) obj.stp);
            % Add 'Next' button (call to dbcont)
            obj.btns.next = uicontrol(h,'Style','pushbutton','String','Next','Units','normalized',...
                'Position',[0.74 0.95 0.1 0.04],'Visible','on','Callback',@(~,~) obj.nxt);
            % Add 'Fwd' button (set obj.pausing = false, and call dbcont)
            obj.btns.fwd = uicontrol(h,'Style','pushbutton','String','Fwd','Units','normalized',...
                'Position',[0.84 0.95 0.1 0.04],'Visible','on','Callback',@(~,~) obj.fwd);
            % Add 'Fwd' button (set obj.pausing = false, and call dbcont)
            obj.btns.sw2wait = uicontrol(h,'Style','pushbutton','String','Stop Plotting',...
                'Units','normalized','Position',[0.5 0.95 0.14 0.04],'Visible','on',...
                'Callback',@(~,~) obj.switchtowaitbar());
        end
        
        function updatewaitbar(obj,varargin)
        % UPDATEWAITBAR(obj,..) - shortcut for obj.fighandle.update(..)
            assert(~obj.plotting,'obj.updatewaitbar requires that ~obj.plotting');
			obj.fighandle.update(varargin{:});
        end

        function y = worthupdating(obj)
            if obj.plotting, y = false; return; end
            y = obj.fighandle.worthupdating;
        end
        
        function refresh(obj)
        % Execute queued instructions, empty queue, and update figure
            switch obj.queue
                case 'plot', obj.switchtoplotter();
                case 'wait', obj.switchtowaitbar();
                case '' % ...
                otherwise, error('Unrecognized queue');
            end
            obj.queue = '';
            drawnow()
        end
        
        function switchtowaitbar(obj,later)
        % SWITCHTOWAITBAR(OBJ) - Close plotter immediately and draw waitbar, stop pausing
        % SEITCHTOWAITBAR(OBJ,TRUE) - Add 'wait' signal to queue & wait for REFRESH
        
            if nargin > 1 && later, obj.queue = 'wait'; return; end
            close(obj,true);
            obj.plotting = false;
            drawwaitbar(obj,0,'Waiting for update...','name',obj.name); 
        end
        
        function switchtoplotter(obj,later)
        % SWITCHTOPLOTTER(OBJ) - Close waitbar immediately and draw plotter, pause
        % SEITCHTOPLOTTER(OBJ,TRUE) - Add 'plot' signal to queue & wait for REFRESH
        
            if ~obj.isUI
               warning('Cannot switch to plotter on non-UI mode');
               return;
            end
            
            if nargin > 1 && later, obj.queue = 'plot'; return; end
            close(obj,true);
            obj.plotting = true;
            drawplotter(obj);
            obj.pausing = false;
        end 

        function nxt(obj,~,~)
        % Callback for Next button
            obj.pausing = true;
            uiresume(obj.fighandle);
        end
        
        function fwd(obj,~,~)
        % Callback for Fwd button
            obj.pausing = false;
            uiresume(obj.fighandle);
        end
        
        function stp(obj,~,~)
            % Callback for Pause button
            obj.pausing = true;
        end 

        function btnh = getbtn(obj,btnid)
        % Get (available) button handles, based on (cell-array-of) string(s) btnid.
        % Valid ID's are:
        %   when obj.plotting: 'sw2wait','pause','next','fwd'
        %   otherwise: 'sw2plot'
        % If btnid is omitted or empty, a vector with all available btn handles will be returned.
        
            if obj.plotting, availablebtns = {'sw2wait','pause','next','fwd'};
            else, availablebtns = {'sw2plot'}; end
            
            if nargin < 2 || isempty(btnid), btnid = availablebtns; end
            
            if ~iscell(btnid), btnid = cellstr(btnid); end
            btnh = zeros(size(btnid));
            for j = 1:numel(btnh)
                assert(any(strcmpi(btnid{j},availablebtns)),'plotter:getbtn',...
                    'Button ID not recognized, or not available');
                btnh(j) = obj.btns.(btnid{j});
            end
        end
        
        function setbtn(obj,btnid,varargin)
        % Shortcut for set(obj.btns.(btnid),varargin{:})
        % btnid can be a cell-array of string ID's (for several buttons) or be left empty, for all 
        % available buttons. See getbtn for valid ID's.
            btnh = getbtn(obj,btnid);
            for j = 1:numel(btnh)
                set(btnh(j),varargin{:});
            end
        end
        
        function close(obj,fromcode)
        % Custom callback for CloseReqFcn, provide options to Continue/Debug/Quit
            if nargin < 2, fromcode = false; end

            if ~fromcode
                switch optquestdlg('What do you want to do?','Cancel Simulation',...
                        'Continue','Quit','Debug','Quit')
                case 'Continue' % do nothing
                case 'Debug'
                    uiresume(obj.fighandle);
                    keyboard
                otherwise
                    delete(obj.fighandle);
                    if ~feature('IsDebugMode')
                        error('plotter:noobj','Closed by user');
                    end
                end
            else
                if ishandle(obj.fighandle), delete(obj.fighandle); end
            end
        end
            
        function delete(obj)
        % Close figure/handle upon obj deletion
            if ishandle(obj.fighandle), close(obj.fighandle); end
        end
        
        function pause(obj)
        % stop execution until uiresume(h) is called by nxt() or fwd()
            uiwait(obj.fighandle);
        end
    end
    methods (Static = true)
        function clearexisting()
        % Check if GUIfigure('plotter') is open, and if so, close it. 
        %
        % An implicit reference to a previous PLOTTER object might be attached to the 
        % CloseRequestFcn of the figure, which causes MATLAB to crash in any attempt to reset it.
        % The problem seems to be that the reference prevents MATLAB from calling the destructor
        % which should close the figure when the object is no longer used.
        % See 'Java Objects Referencing MATLAB Objects' on:
        % https://www.mathworks.com/help/matlab/matlab_oop/handle-class-destructors.html#btavjj2-1
        
            h = GUIfigure('plotter','-silent');
            if ishandle(h), delete(h); end
        end
    end
end
