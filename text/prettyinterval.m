function msg = prettyinterval(a,b,tol)
% MSG = PRETTYINTERVAL(DT) - format a time-duration interval DT
% MSG = PRETTYINTERVAL(A,B) - format an interval A-B (round to days!)
%
% MSG = PRETTYINTERVAL(DT) - formats a time interval DT (DURATION, or numeric scalar days) into 
%   a string with "reasonable" units, namely:
%
%       'Y.Y years'  for DT > 366
%       'M.M months' for 30 < DT <= 366
%       'D.D days'   for 3 < DT <= 30
%       'HH:MM:00'   for 1/24 < DT <= 3
%       'HH:MM:SS'   for DT up to 72 hours, 0.00 days for dt, X.XX weeks, ...
%
% MSG = PRETTYINTERVAL(A,B) - format the interval A-B (any A, B accepted by PARSETIME) as 
%   YYYYMMDD-YYYYMMDD, reducing the field(s) for year, month, and day if they are common between
%   A and B, e.g.:
%
%        20190101-0311   instead of 20190101-20190311
%        20190301-11     instead of 20190301-20190311  
%
%   Also, if the interval A-B covers a complete* year/month, all ensuing units are removed, e.g.:
%
%        201908         instead of 20190801-20190831
%        2019           instead of 20190101-20191231
%
% MSG = PRETTYINTERVAL(A,B,TOL) - (*) by default, B is expected to represent the end of the last 
%   day, and A the start of the first day of a non-empty interval. If this is not the case, a 
%   warning will be issued, and floor(A), ceil(B) will be used.
%
%   A tolerance TOL (DURATION or scalar hours) can be used to get a more readable, approximate
%   description of an interval. All combinations of integer days within [A - TOL, A + TOL] and
%   [B - TOL, B + TOL] will then be tested, and the shortest description selected.
%   e.g. prettyinterval('2019-12-31 22:00','2020-12-31 23:00',hours(2)) -> '2020'
%
% TODO: deal with sub-day intervals, change to ISO representation, at least for whole intervals 
%   (year/months) shifted by one day, e.g. 20191231P1Y or 20190802P1M.
%   Overload reverse syntax: [A,B] = prettyinterval(msg)
%
% See also: PARSETIME, ISODURATION

    narginchk(1,3);
    if nargin < 2 || isempty(b), b = NaN; end

    if nargin < 2
        
        if isduration(a), a = days(a); end
        assert(isnumeric(a) && isscalar(a) && isreal(a),'Unrecognized duration DT');
        
        if ~isfinite(a) || a > 128e+9  % (&) 128,849,018,881 is the limit for datevec()
            msg = 'NA'; 
            return; 
        end

        if a > 366, msg = sprintf('%0.4g years',a/365.25);
        elseif a > 30, msg = sprintf('%0.3g months',a/365.25*12);
        elseif a > 3, msg = sprintf('%0.3g days',a); 
        else
            if a > 1/24
                a = round(a*24*60)/(24*60); % round to minutes
            else
                a = round(a*24*3600)/(24*3600); % round to seconds
            end
            msg = datestr(a,'HH:MM:SS');
        end
    else       
        if nargin < 3 || isempty(tol), tol = 6/24; end
        
        a = parsetime(a,'TimeZone','keep');
        b = parsetime(b,'TimeZone','keep');
        if isnat(a) || isnat(b), msg = 'NA'; return; end  
        
        assert(b >= a && isequal(a.TimeZone,b.TimeZone),'Expecting b > a with equal time zones');

        if isduration(tol), tol = days(tol); end
        assert(isnumeric(tol) && isscalar(tol) && isfinite(tol) && isreal(tol) && ...
            tol >= 0 && tol < 31,'Expecting tolerance in the range [0,31] days'); 
        
        daystart = @(x) dateshift(x,'start','day','nearest');
        dayend = @(x) dateshift(x-eps(1),'end','day'); % keeps day starts/ends untouched
        
        if tol == 0
            A = daystart(a); B = dayend(b);
        else
        % Play around with tolerances...
            A = dayend(a-tol):daystart(a+tol);
            B = dayend(b-tol):daystart(b+tol);
            [A,B] = meshgrid(A,B);
        end
        
        dA = abs(A - a);
        dB = abs(B - b);
        valid = B > A & dB <= tol & dA <= tol;
        if ~any(valid)
            A = daystart(a); B = dayend(b); valid = true;
            warning('prettyinteval:notint','Interval is not an integer number of days within TOL');
        end
        
        if nnz(valid) == 1
            msg = shortdatestr(A(valid),B(valid));
        else
        % ... and pick the shortest description
            AB = unique([A(valid),B(valid)],'rows','stable');   
            [msgs,len] = arrayfun(@shortdatestr,AB(:,1),AB(:,2),'unif',false);
            len = [len{:}];
            n = min(len);
            idx = (len == n); % (possibly with ties)
            if nnz(idx) == 1
                msg = msgs{idx};
            else
            % ... that most closely describes the interval
                valid(valid) = idx;
                idx = find(idx);
                [~,j] = min(dA(valid) + dB(valid));
                msg = msgs{idx(j)};
            end
        end 
    end
end

function [msg,len] = shortdatestr(a,b)

    b = b-1; % start of last day
    A = datevec(a); 
    B = datevec(b);

    if year(a-1) < A(1) && year(b+1) > B(1)
    % whole-year coverage, e.g. 2019-2021
        msg = num2str(A(1));
        if B(1) ~= A(1)
           msg = [msg,'-',num2str(B(1))]; 
        end
    else
        if A(3) <= 1 && month(b+1) ~= B(2)
        % whole-month coverage, e.g. 201001-2021-05, 201901-03
            msg{1} = datestr(A,'yyyymm');
            msg{2} = datestr(B,'yyyymm');
            hasday = false;
        else
        % custom dates, e.g. 20190301-20210311, 20190101-0311, 20190301-11
            msg{1} = datestr(A,'yyyymmdd');
            msg{2} = datestr(B,'yyyymmdd');
            hasday = true;
        end
        if strcmp(msg{1}(1:4),msg{2}(1:4))
            msg{2}(1:4)=[];                     % remove common year
            if strcmp(msg{1}(5:6),msg{2}(1:2))
                msg{2}(1:2) = [];               % remove common month
                if hasday && strcmp(msg{1}(7:8),msg{2}(1:2))
                    msg{2}(1:2) = [];           % remove common day
                end
            end
            if isempty(msg{2}), msg(2) = []; end
        end
        msg = strjoin(msg,'-');
    end    
    
    len = numel(msg);
end