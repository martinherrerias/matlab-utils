function varargout = synchronizestruct(varargin)
% [A,B,C,..] = SYNCHRONIZESTRUCT(A,B,C,..) - TIMETABLE/SYNCHRONIZE how it should work.
%   Timetables (A,B,C..) are are handled separately, downsampled to 'step',F (F = @max by default,
%   i.e. set to the minimum resolution. Using @min, @mean,.. or just 'step',DURATION also works.
%
%   Time stamps must be START-OF-INTERVALS. They may hava a non zero offset with respect to hour 
%   starts, e.g. one might go 05:00,10:00,15:00 and the next one 01:00,06:00,11:00... this will 
%   be corrected by resampling/interpolation to match a common 'offset',X (0 by default)
%
%   TODO: resampling of non-doubles (logicals, datetime, integers,..)
%         move internal resample to resamplestructure?
%
% See also: TIMETABLE/SYNCHRONIZE, RETIMESTRUCT

    opt.offset = 0;
    opt.step = @max;
    [opt,varargin] = getpairedoptions(varargin,opt);
    assert(numel(varargin) > 1);

    [T,t,dt,idx,c] = cellfun(@parsetable,varargin,'unif',0);
    dt = cat(2,dt{:});
    assert(all(isfinite(dt)),'Timesteps are not regular');

    % Get resampling (m > 1) / downsampling (m < 1) rate 
    if isa(opt.step,'function_handle')
        m = dt./opt.step(dt);
    else
        assert(isduration(opt.step))
        m = dt./opt.step; 
    end
    dt = unique(dt./m);
    assert(isscalar(dt));
    x = -arrayfun(@(t) offset(t{1},dt,opt.offset),t); % units of new dt
    [T,t] = arrayfun(@resample,T,idx,t,m,x,'unif',0);
    
    [tc,idx{1},idx{2}] = intersect(t{1},t{2});
    for j = 3:numel(T)
        if isempty(tc)
            idx{j} = [];
            continue;
        end
        [tc,iA,idx{j}] = intersect(tc,t{j});
        idx(1:j-1) = cellfun(@(x) x(iA),idx(1:j-1),'unif',0);
    end

    varargout = cell(1,nargout);
    for j = 1:nargout
        if isa(T{j},'MeteoData')
            varargout{j} = filterstructure(T{j},idx{j},T{j}.Nt);
            varargout{j}.data.Properties.RowTimes = tc;
        else
            varargout{j} = T{j}(idx{j},:);
            varargout{j}.Properties.RowTimes = tc;
        end
        if isempty(c{j}), continue; end
        varargout{j} = c{j}(varargout{j}); % convert back non-table outputs
    end
end

function x = offset(t,dt,xC)
    x = mod(t(1) - dateshift(t(1),'start','hour') - xC,dt)/dt;
    if x > 0.5, x = x - 1; end
end

function [X,t,dt,idx,C] = parsetable(X)
    C = [];
    if iscell(X)
        [X,c] = struct2timetable(X{:});
        C = @(x) completestruct(c,timetable2struct(x));
    elseif isdatetime(X)
        X = timetable(X);
        C = @(x) x.Properties.RowTimes;
    elseif isa(X,'MeteoData')
        if ~isregular(X.data), X = checktimestamps(X); end
        X.interval = 'b';
        t = X.t;
        dt = X.timestep;
        idx = (1:X.Nt)';
        return;
    else
        assert(istimetable(X));
    end

    dt = X.Properties.TimeStep;
    if ~isfinite(dt), dt = []; end
    [t,dt,idx] = parsetime(X.Properties.RowTimes,'-regular','-unique','step',dt,'interval','b');
    if ~isfinite(dt)
        X.Properties.RowTimes = t(idx);
        X.Properties.TimeStep = dt;
    end
end

function [T,t] = resample(T,idx,t,m,offset)

    if iscell(T), T = T{1}; end
    if iscell(idx), idx = idx{1}; end
    if iscell(t), t = t{1}; end
    
    n = numel(idx);  
    assert(size(T,1) == n);

    if isequal(m,1) && offset == 0
        if isequal(idx(:)',1:n)
            t = parsetime(t,'-regular','interval','b');
            return; 
        end
        T = retime(T,t,'fillwithmissing');
        % T{idx,:} = T{:,:};
        % F = sparse(idx,1:n,1);
    else
        if ~isequal(m,1)
            [p,q] = rat(m);
            if p > q, msg = 'Resampling'; else, msg = 'Downsampling'; end
            msg = sprintf('%s %s by %d:%d',msg,inputname(1),p,q);
            m = [p,q];
        else
            msg = sprintf('Offseting %s (1/1 resample)',inputname(1));
        end
        if offset ~= 0
            msg = sprintf('%s, %0.2g step offset',msg,offset);
        end
        [T,F] = resamplestructure(T,m,'offset',offset);
        
        
        t0 = t(1);
        % t = t0 + F*days(t-t0);
        t = resamplestructure(days(t-t0),m,'offset',offset) + t0;
        t = parsetime(t,'-regular','interval','b'); % resolve numerical precision errors
        
        if isempty(T)
            T = timetable(t);
            return;
        end
%         warning('sync:resample','%s. All values will be converted to doubles',msg);
%             
%         % TODO: jump sparse*double propperly, parse-reparse other numeric types, logicals, etc. 
%         N = numel(t);
%         F = F*sparse(idx,1:n,1);
%         w = F*ones(n,1);
%         T{1:N,:} = F*double(T{1:n,:})./w;
%         if size(T,1) > N, T = T(1:N,:); end % downsampling
%         T.Properties.RowTimes = t;
    end
end

