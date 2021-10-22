function s = deg2dms(ang, varargin)
% S = DEG2DMS(ANG, [DIRECTION])
% Takes a decimal degree angle and returns a string of the form '+00°00'00' or '00°00'00.00" E'
%   ANG provides the angle to convert, in degrees. When ANG is an array, the result will 
%       be a numel(ANG)-row char-array, with spaces used for padding, if necessary.
%   DIRECTION {'EW','NS'}, if provided, the result will be unsigned and followed by the
%       letters N, S, E, W, depending on sign (N and E are positive). Otherwise the result
%       will be explicitly-signed (+0°00'.. / -0°00'..).
% S = DEG2DMS(..,'precision',K) - changes the rounding precision of the result. K = 'nU' uses 
%       n decimals for unit U {d,m,s}. Default is '0s'. A numeric K is equivalent to 'Ks'.
% S = DEG2DMS(..,'degsym',c) uses custom character c instead of '°' for the degree symbol
%            (..,'minsym',c) and (..,'secsym',c) do the same for minutes, and seconds.
% EXAMPLES:
%     DEG2DMS(0) = +0°00'00"
%     DEG2DMS(-1,'precision',1) = -1°00'0.0"
%     DEG2DMS(2,'precision','2s') = +2°00'0.00"
%     DEG2DMS(-3,'precision','3m') = -3°0.000'
%     DEG2DMS(4.5,'EW') = 04°30'00" E
%     DEG2DMS(-6.7,'NS') = 06°42'00" S
%     DEG2DMS(1.5,'precision','0m','degsym',':','minsym','') = +1:30

% See also: DMS2DEG

    % if nargin == 0, test(); return; end

    [opt,remargs] = getpairedoptions(varargin,{'precision','degsym','minsym','secsym'},...
                                              {0          ,'°'     ,''''    ,'"'});
    assert(numel(remargs) < 2, 'Too many non-paired arguments');
    if isempty(remargs), direction = ''; else, direction = remargs{1}; end
    assert(isnumeric(ang),'Expecting numeric input');

    n = numel(ang);
    degformat = '%02';
    switch upper(direction)
        case 'EW'
            sufix = repmat(' E',n,1);
            sufix(ang < 0,2) = 'W';
        case 'NS' 
            sufix = repmat(' N',n,1);
            sufix(ang < 0,2) = 'S';
        otherwise
            assert(isempty(direction),'Unrecognized direction');
            sufix = char.empty(n,0);
            degformat = '%+02'; % use explicit signs when not printing direction
    end

    ang = ang(:);
    if ~isempty(sufix), ang = abs(ang); end
    
    lastorder = 'e'; % for Error
    if ischar(opt.precision)
        p = regexpi(opt.precision,'^(\d*)([dms])$','tokens');
        if ~isempty(p)
            lastorder = p{1}{2};
            precision = str2double(p{1}{1});
        end
    elseif isnumeric(opt.precision)
        lastorder = 's';
        precision = opt.precision;
    end
    assert(any(lastorder == 'dms') && mod(precision,1) == 0 && precision >= 0,...
        'Unrecognized precision setting');

    if lastorder == 'd'
        deg = ang;
        fmt = [degformat '.' num2str(precision) 'f' opt.degsym];
        s = arrayfun(@(d) sprintf(fmt,d),deg,'unif',0);
    else
        deg = fix(ang);
        ang = abs(ang - deg);
        min = ang * 60;
        fmt = [degformat 'd' opt.degsym];
        if lastorder == 'm'
            fmt = [fmt,'%02.',num2str(precision),'f',opt.minsym];
            s = arrayfun(@(d,m) sprintf(fmt,d,m),deg,min,'unif',0);
        else
            min = fix(min);
            ang = ang - min / 60;
            sec = ang * 3600;
            fmt = [fmt '%02d' opt.minsym '%02.' num2str(precision) 'f' opt.secsym '%c'];
            s = arrayfun(@(d,m,s) sprintf(fmt,d,m,s),deg,min,sec,'unif',0);
        end
    end
    s = [char(s),sufix];
end

% function test()
%     CALL = {'deg2dms(0)';
%             'deg2dms(-1,''precision'',1)';
%             'deg2dms(2,''precision'',''2s'')';
%             'deg2dms(-3,''precision'',''3m'')';
%             'deg2dms(4.5,''EW'')';
%             'deg2dms(-6.7,''NS'')';
%             'deg2dms(1.5,''precision'',''0m'',''degsym'','':'',''minsym'','''')'};
%     
%     % RES = cellfun(@eval,CALL,'unif',0);
%     RES = { '+0°00''00"';
%             '-1°00''0.0"';
%             '+2°00''0.00"';
%             '-3°0.000''';
%             '04°30''00" E';
%             '06°42''00" S';
%             '+1:30'};
%     assert(isequal(cellfun(@eval,CALL,'unif',0),RES));
%         
%     cellfun(@(x) fprintf('%s = %s\n',x,eval(x)),CALL);
% end