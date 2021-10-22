classdef stopwatch < handle
% STOPWATCH Class
%
% Keep track of the computation time assigned to different tasks (LAPS) within a simulation, 
% and/or remember the time of particular events (SPLITS).
%
% EXAMPLE:
%
%   foo = @() pause(1); 
%   bar = @() pause(2);
%   
%   s = stopwatch();    % initialize object
%   foo();
%   s.add2lap('foo');   % add the time spent since s.lasttoc (STOPWATCH creation) to new lap 'foo' 
%   bar();
%   s.add2lap('bar');   % add the time spent since s.lasttoc (last call to add2lap) to new lap 'bar'
%   s.addsplit('bar');  % remeber the absolute (accumulated) current time, store it in split 'bar'
%   s.addsplit('baz');  % store it also in new split 'baz' (*)
%   foo();
%   s.add2lap('foo');   % add the time spent since s.lasttoc to existing lap 'bar'
%   s.addsplit('baz');  % rewrite existing split 'baz' - notice prior data is lost (*)
%   s.printrecord();
% 
%   foo();
%   s.resetlocal();     % set s.lasttoc to current time
%   s.add2lap('qux');   % the time since s.lasttoc is now ~0 since resetlocal
%   s.addsplit('qux');  % resetlocal doesn't affect absolute time (i.e. splits)...
%   s.resetglobal();
%   s.addsplit('quux'); % ... but resetglobal does
%   s.printrecord();
%
% See also: TIC, TOC, PROFILE

	properties
		clockID  % output of tic();
		lasttoc  % output of last toc();
        laps     % laps.(x) time spent on task x (seconds) 
        splits   % splits.(x) (absolute) time of event x (seconds)
    end
    properties (Hidden = true)  
        n   % persistent variable for REMAININGTIME estimation
    end
    methods (Static)
        function obj = loadobj(a)
            if isstruct(a)
            % Backwards compatibility: tags/intervals fields
                obj = stopwatch();
                obj.laps = cell2struct(num2cell(a.intervals(:)),a.tags(:));
            elseif isa(a,'stopwatch')
            % do nothing, really
                obj = a; 
            end
        end
    end
	methods
		function obj = stopwatch()
			obj.clockID = tic;
            obj.lasttoc = toc(obj.clockID);
            obj.laps = struct();
            obj.splits = struct();
            obj.n = 0;
		end
		
		function t = localtime(obj)
			t = toc(obj.clockID)-obj.lasttoc;
		end
		
		function t = globaltime(obj)
			t = toc(obj.clockID);
		end
		
		function resetlocal(obj)
			obj.lasttoc = toc(obj.clockID);
		end
		
		function resetglobal(obj)
			obj.clockID = tic;
			obj.lasttoc = toc(obj.clockID);
            obj.n = 0;
		end
		
		function add2lap(obj,field)
		% ADD2LAP(OBJ,NAME) - Adds the current OBJ.localtime to OBJ.laps.NAME
		%                     NOTE: If no such field exists, a new one will be created.
            
            if isfield(obj.laps,field), t0 = obj.laps.(field); else, t0 = 0; end
            obj.laps.(field) = t0 + obj.localtime;
			obj.resetlocal();
        end
        
        function addsplit(obj,field)
        % ADDSPLIT(OBJ,NAME) - shortcut for OBJ.splits.NAME = globaltime()
            
            obj.splits.(field) = obj.globaltime;
        end

		function varargout = printrecord(obj)
        % PRINTRECORD(OBJ) - Display all laps and splits with their respective time (HH:MM:SS).
                    
            f{1} = fieldnames(obj.laps);
            v{1} = cellfun(@(f) timestr(obj.laps.(f)),f{1},'unif',0);
            
            f{2} = fieldnames(obj.splits);
            v{2} = cellfun(@(f) timestr(obj.splits.(f)),f{2},'unif',0);
            
            nf = max(cellfun(@numel,cat(1,f{:})));
            nv = max(cellfun(@numel,v{1}));
            nf = max(nf,8); nv = max(nv,8); % make room for title 'Splits:'
            
            t{1} = {'Laps:'};
            t{2} = cellfun(@(f,v) [blanks(nf-numel(f)),f,': ',blanks(nv-numel(v)),v],f{1},v{1},'unif',0);
            t{3} = {'Splits:'};
            t{4} = cellfun(@(f,v) [blanks(nf-numel(f)),f,': ',blanks(nv-numel(v)),v],f{2},v{2},'unif',0);
            
            if isempty(t{2}), t{2} =  {[blanks(nf+nv-4) '(none)']}; end
            if isempty(t{4}), t{4} =  {[blanks(nf+nv-4) '(none)']}; end
            
            t = cat(1,t{:});
            if nargout > 0, varargout{1} = strjoin(t,newline());
            else
                fprintf('\n');
                fprintf('\t%s\n',t{:});
                fprintf('\n');
            end
            
            function s = timestr(x)
                ds = @(x) [num2str(floor(x/3600)) ':' datestr(x/(24*3600),'MM:SS')];
                
                if isscalar(x), s = ds(x);
                else
                   s = [ds(min(x)) ' - ' ds(max(x))];
                end
            end
        end
        
        function plot(obj)
        % PLOT(OBJ) - Bar plot for laps, line plot for splits.
        
            f = fieldnames(obj.laps);
            v = cellfun(@(f) obj.laps.(f),f);
            v = v/sum(v);
            
            GUIfigure('stopwatch','Timing',[0.1 0.1 0.6 0.4]);
            subplot(1,2,1)
            bar(v);
            ylim([0,1]);
            xticks(1:numel(f));
            xticklabels(f); grid on;
            ylabel('fraction of time')
            set(gca,{'XTickLabelRotation','TickLabelInterpreter'},{90,'none'});
            
            f = fieldnames(obj.splits);
            v = cellfun(@(f) obj.splits.(f)/60,f,'unif',0);
            nv = max(cellfun(@numel,v));
            v = cellfun(@(x) [x nan(1,nv-numel(x))],v,'unif',0);
            
            subplot(1,2,2)
            plot(cat(1,v{:}))
            xticks(1:numel(f));
            xticklabels(f); grid on;
            ylabel('minutes');
            set(gca,{'XTickLabelRotation','TickLabelInterpreter'},{90,'none'});
        end
        
        function M = merge(varargin)
        % M = merge(A,B,C..) - return a STOPWACH object M that has all the tags of objects A,B,C..
        % with their respective intervals added (whenever tags match).
        % clockID and lasttoc will be set to those of the object with the greatest GLOBALTIME.
       
            clocks = [varargin{:}];

            lapnames = arrayfun(@(c) fieldnames(c.laps),clocks,'unif',0);
            splitnames = arrayfun(@(c) fieldnames(c.splits),clocks,'unif',0);
            
            laplist = unique(cat(1,lapnames{:}),'stable');
            lapval = zeros(numel(laplist),1);
            
            splitlist = unique(cat(1,splitnames{:}),'stable');
            splitval = repmat({zeros(1,0)},numel(splitlist),1);
            
            if ~isempty(laplist)
                for j = 1:numel(clocks)
                    [~,iM] = ismember(lapnames{j},laplist);
                    lapval(iM) = lapval(iM) + cellfun(@(f) clocks(j).laps.(f),lapnames{j});
                end
            end
            if ~isempty(splitlist)
                for j = 1:numel(clocks)
                    [~,iM] = ismember(splitnames{j},splitlist);
                    splitval(iM) = cellfun(@(v,f) {[v,clocks(j).splits.(f)]},splitval(iM),splitnames{j});
                end
            end
            
            M = stopwatch();
            M.laps = cell2struct(num2cell(lapval(:)),laplist);
            M.splits = cell2struct(splitval(:),splitlist);
            M.lasttoc = max([clocks.lasttoc]);
            M.splits.merge = M.lasttoc;
        end
        function M = mergestructures(C,varargin)
        % Overloaded MERGESTRUCTURES for STOPWATCH objects == alias to MERGE
            M = merge(C);
        end
	end
end
