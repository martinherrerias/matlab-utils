classdef optwaitbar < handle
% OPTWAITBAR - class designed to wrap MATLAB's WAITBAR with added functionality:
%
%   + Command-line updates, of the form '#WAITBAR:(0.00%) MSG' for non-UI environments
%   + Integrated clock (STOPWATCH object) to provide automatic ETA estimates & messages
%   + Automatic closing of WAITBAR figure upon deletion of the wrapper OPTWAITBAR object
%
% The constructor is designed to take the same syntax as WAITBAR, but returns an OPTWAITBAR 
% object handle with the actual WAITBAR's handle (if running on UI mode) in OBJ.ID. Updates
% to an existing waitbar of the form OPTWAITBAR(X,H,..) where H is an OPTWAITBAR object
% work too. The only detail is that the constructor's output (H) cannot be optional, and tends
% to stay in the working space as ANS, potentially messing up automatic figure closure. Use of
% H.UPDATE(X,MSG,..) is therefore recommended.
%
% EXAMPLE:
%     h = optwaitbar(0,'starting...');   
%     N = 20;
%     pause(5)
%     figure('Name','optwaitbar history'); hold on
%     for j = 1:N
%         plot(h.elapsedtime,j/N,'o')
%         h.update(j/N,'waiting...','-addtime');
%         errorbar(h.elapsedtime+h.eta,1,h.eta_uncertainty,'horizontal');
%         pause(rand()+0.5);
%     end
%     clear h;
%
% See also: WAITBAR, RUNNINGFROMUI, STOPWATCH

properties (Transient = true, GetAccess = public, SetAccess = private) %(*)
    ID      % WAITBAR figure handle if OBJ.isUI, NaN otherwise.
    clock   % STOPWATCH object, created 
    isUI    % RUNNINGFROMUI() at the moment of OPTWAITBAR construction
end
properties (GetAccess = public, SetAccess = protected)
    msg              % last updated message cellstring, not including -addtime (see UPDATE)
    progress         % last updated progress [0,1]
    eta              % estimated time to reach progress = 1, based on past updates (seconds)
    eta_uncertainty  % uncertainty in OBJ.eta, based on past updates (seconds)
    memory           % records used to estimate OBJ.eta / OBJ.eta_uncertainty
end
properties (SetAccess = public, GetAccess = public)
    TAG         % when ~OBJ.isUI, print waitbar messages as: 'TAG: OBJ.msg', default is '#WAITBAR'
    LEASTWAIT   % function dt = @(t,x) such that OBJ.LEASTWAIT(OBJ.elapsedtime,OBJ.progress)
                % returns the minimum time (in seconds) between successive waitbar updates.
                % Default is @(t,x) min(60,max(1,t*5/60)), i.e. ramp from 1" to 60" in 12'.
end
properties (Dependent = true)
   elapsedtime      % for OBJ.clock.globaltime
   worthupdating    % minimum time between successive waitbar updates has passed
end
methods (Static = true)
    function obj = loadobj(obj,varargin)
    % Create transient propperties

        obj.isUI = runningfromUI();
        obj.clock = stopwatch();

        if obj.isUI
        % On interactive mode, draw standard waitbar, and create a figure-closer
            obj.ID = waitbar(obj.progress,varargin{:});
        else
            obj.ID = NaN;
        end
    end
end
methods
    function obj = optwaitbar(x,varargin)
    % H = OPTWAITBAR(X,MSG,..) displays/updates a waitbar figure if running MATLAB in interactive
    %    mode, otherwise it prints the message and progress-fraction to the console.
    %
    %   Interactive mode H.isUI is determined by RUNNINGFROMUI(). If H.isUI, OPTWAITBAR(...) 
    %   calls WAITBAR(...) and keeps the resulting figure handle in H.ID; otherwise no figure
    %   is created, and all updates will be printed to the console:
    %
    %       #WAITBAR(XXX.X%) MSG{1}. MSG{2}. ..
    %
    %   OPTWAITBAR is designed to take the same syntax as WAITBAR, but returns an OPTWAITBAR 
    %   object handle with an actual WAITBAR's handle (if running on UI mode) in OBJ.ID. Updates
    %   to an existing waitbar, of the form OPTWAITBAR(X,H,..) where H is an OPTWAITBAR object
    %   work too. The only detail is that this constructor's output (H) is not optional, so it 
    %   may stay in the working space as ANS, potentially messing up automatic figure closure. 
    %   Use of H.UPDATE(X,MSG,..) is therefore recommended for existing objects.
    %
    %   OPTWAITBAR(..,name,value) any name-value pairs will be passed to WAITBAR (when running in
    %   interactive mode), and ignored otherwise.
    %
    % See also: OPTWAITBAR, OPTWAITBAR.UPDATE, WAITBAR, RUNNINGFROMUI, STOPWATCH

        % Parse input
        if nargin < 1 || isempty(x), x = 0; end

        if ~isempty(varargin) && isa(varargin{1},'optwaitbar')
        % Update existing OPTWAITBAR object
            obj = varargin{1}; 
            if obj.isUI, assert(ishandle(obj.ID),'Waitbar has been closed'); end
            
            varargin(1) = [];
            obj.update(x,varargin{:});
            return;
        else
        % Create a new object...
            obj.msg = '';
            obj.progress = x;
            obj.eta = NaN;
            obj.eta_uncertainty = NaN;
            obj.memory = struct('Sx',0,'Sxx',0,'n',0);
            obj.TAG = '#WAITBAR';
            obj.LEASTWAIT = @(t,x) min(60,max(1,t*5/60)); % seconds
        
            obj = optwaitbar.loadobj(obj,varargin{:});
            obj.update(x,varargin{:});
            return;
        end
    end

    function delete(obj)
    % Close figure upon object destruction
        if ~isempty(obj.ID) && ishandle(obj.ID) && isvalid(obj.ID), delete(obj.ID); end
    end  

    function reset(obj)
    % RESET(OBJ) - Set OBJ.progress to zero, and reset all propperties used to get OBJ.eta
        obj.progress = 0;
        obj.eta = NaN;
        obj.eta_uncertainty = NaN;
        obj.memory = struct('Sx',0,'Sxx',0,'n',0);
        obj.clock.resetlocal();
    end
      
    function update(obj,x,varargin)
    % UPDATE(OBJ,X,MSG,...) - Update progress & message on UI waitbar figure or print a line to 
    %   the log, as long as OBJ.worthupdating (a minimum predefined time has past since the last
    %   update). If X > OBJ.progress (i.e. if there has been any progress since the last update)
    %   OBJ.progress, OBJ.eta, and OBJ.eta_uncertainty will be updated accordingly.
    %
    % UPDATE(..,'-addtime') - Append an additional line to the message, of the form:
    %   X elapsed, Y remaining. Where X, Y are the formatted OBJ.elapsedtime and OBJ.eta.
    %
    % UPDATE(..,'-worthy') - Override check of OBJ.worthupdating, and update no matter what.
    %
    % See also: OPTWAITBAR
    
        [opt,varargin] = getflagoptions(varargin,{'-addtime','-worthy'});
        opt.TAG = [];
        opt.LEASTWAIT = [];
        [opt,varargin] = getpairedoptions(varargin,opt);
        if ~isempty(opt.TAG), obj.TAG = opt.TAG; end
        if ~isempty(opt.LEASTWAIT), obj.LEASTWAIT = opt.LEASTWAIT; end
        
        x = max(0,min(x,1)); % Map any value of x to a value between 0 and 1
        
        dt = obj.clock.localtime();
        dy = x - obj.progress;
        m = dy/dt;
        
        % Don't do anything if trying to update within LEASTWAIT...
        if ~opt.worthy && ~obj.worthupdating && ~(m > 0 && obj.memory.n < 2)
            return;
        end

        if dy > 0
        % Estimate remaining time based on the current slope m, and all its past records

            % Update sums and moments Sx = sum(m), and Sxx = sum(m²), 
            obj.memory.Sx = obj.memory.Sx + m;
            obj.memory.Sxx = obj.memory.Sxx + m^2;
            obj.memory.n = obj.memory.n + 1;

            m = obj.memory.Sx/obj.memory.n;                                     % mean slope
            sm = sqrt((obj.memory.Sxx/obj.memory.n - m^2)/(obj.memory.n - 1));  % std. of m

            obj.eta = (1-x)/m;
            obj.eta_uncertainty = 1.96*obj.eta*sm/m;   % 95% expanded uncertainty on eta

            obj.progress = x;
            obj.clock.resetlocal();
        end
        
        if ~isempty(varargin) && (ischar(varargin{1}) || iscellstr(varargin{1}))
           obj.msg = varargin{1};
           varargin(1) = [];
        end
        if ~iscell(obj.msg), obj.msg = {obj.msg}; end

        message = obj.msg;
        if opt.addtime    
            message{end+1} = sprintf('%s elapsed, %s remaining',...
                            prettyinterval(obj.elapsedtime/(24*3600)),...
                            prettyinterval(obj.eta/(24*3600)));
            if isfinite(obj.eta_uncertainty)
                message{end} = ...
                    sprintf('%s±%0.0f%%',message{end},100*obj.eta_uncertainty/obj.eta);
            end
        end

        if obj.isUI
        % On interactive mode, update standard waitbar
            waitbar(x,obj.ID,strjoin(message,newline()),varargin{:});
        else
        % Print the message to command line

            fprintf('%s(%05.1f%%) %s\n',obj.TAG,x*100,singleline(message));
            % drawnow('update')
            % diary off; diary on; % DEBUG: force print to log
            % evalin('base','diary off; diary on; drawnow(''update'');');
        end

        function s = singleline(c)
        % S = SINGLELINE(C,SEP) - strjoin cell-string C, replacing any new-lines/jumps with '.'
            SEP = '. ';
            c = cellfun(@strtrim,c,'unif',0)';
            s = strjoin(c,SEP);
            s = regexprep(s,'[\r\n]+',SEP);
        end
    end
    
    function x = get.elapsedtime(obj), x = obj.clock.globaltime(); end
    function y = get.worthupdating(obj)
        y = obj.clock.localtime() >= obj.LEASTWAIT(obj.clock.globaltime(),obj.progress);
    end
    
    function y = ishandle(obj), y = isvalid(obj) && ishandle(obj.ID); end
    function uiresume(obj), uiresume(obj.ID); end
end
end

