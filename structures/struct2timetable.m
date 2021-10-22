function [TT,R] = struct2timetable(S,t,varargin)
% [TT,REST] = STRUCT2TIMETABLE(S,T) - Convert (nested) structure S into a TIMETABLE object with 
%   timestamps T. Fields of S that do not have numel(T) rows will be returned in complementary
%   structure REST (or, if nargout < 2, discarded with a warning).
%
% [TT,REST] = STRUCT2TIMETABLE(S,F) - Specify a field name F instead of T. In that case S.(F) is 
%   used as the time vector, and the dimension name TT.Properties.DimensionNames{1} is set to F
%   i.e. T = TT.Properties.RowTimes = S.(F) = TT.(F)
%
% [TT,REST] = STRUCT2TIMETABLE(S,T,...) - Pass additional arguments to TIMETABLE. 
%
%  'VariableDescriptions',C - Set by default to the original (nested) structure names.
%  'VariableNames',C - By default valid, non-conflicting names extracted from struct. field names.
%  'DimensionNames',C - TIMETABLE defaults {'Time','Variables'}, with the exception above.
%
% See also: TIMETABLE2STRUCT, PARSETIME

    narginchk(2,Inf);
    if isempty(S), S = struct(); end
    
    DEFNAMES = defdimnames(); % {'Time','Variables'} as of R2020a
    
    if ischar(t) && strcmp(t,matlab.lang.makeValidName(t))
        assert(isfield(S,t),'Expecting time vector or structure field name');
        tname = t;
        t = S.(tname);
        S = rmfield(S,tname);
        OPT.DimensionNames = [{tname}, DEFNAMES(2)];
    else
        OPT.DimensionNames = DEFNAMES;
    end
    t = parsetime(t,'TimeZone','keep');
    
    [values,names] = nestedstruct2cell(S);
    valid = cellfun(@(c) size(c,1) == numel(t),values);
    R = struct();
    if ~all(valid)
        if nargout < 2
            warning('Discarding %s',shortliststr(names(~valid),'field'));
        else
            R = cell2nestedstruct(values(~valid),names(~valid));
        end
        values = values(valid);
        names = names(valid);
    end
    
    OPT.VariableNames = [];
    OPT.VariableDescriptions = names;
    [OPT,varargin] = getpairedoptions(varargin,OPT);

    if isempty(OPT.VariableNames)
        OPT.VariableNames = matlab.lang.makeValidName(names);
        OPT.VariableNames = matlab.lang.makeUniqueStrings(OPT.VariableNames,OPT.DimensionNames);
    end

    % (§) According to the documentation, TIMETABLE should take 'DimensionNames' as parameter, 
    % but it doesn't:
    %   TIMETABLE(T,..) takes INPUTNAME for the time vector (in this case 't')
    %   TIMETABLE(..,'RowTimes',t) takes the default 'Time' (DfltRowDimName)
    %
    % In order to avoid a (possibly non-existing) "conflicting name" crash, e.g.
    %   struct2timetable(..,'VariableNames',{'Time',..},'DimensionNames',{'not_Time',..})
    %
    %  ... use a dummy, conflict-free set of variable names (DimensionNames = default)
    varnames = matlab.lang.makeUniqueStrings(OPT.VariableNames,[DEFNAMES,OPT.DimensionNames]);
    TT = timetable(t,values{:},'VariableNames',varnames,varargin{:});
    
    TT.Properties.DimensionNames = OPT.DimensionNames; % ... set the actual DimensionNames
    
    if ~isequal(varnames,OPT.VariableNames)
    % ... and finally (try to) set the original VariableNames
        TT.Properties.VariableNames = OPT.VariableNames;
    end
    
    TT.Properties.VariableDescriptions = OPT.VariableDescriptions;
    
    function names = defdimnames()
    % Taken from TIMETABLE (returns {'Time','Variables'} as of R2020a)
        names = { getString(message('MATLAB:timetable:uistrings:DfltRowDimName')), ...
                  getString(message('MATLAB:timetable:uistrings:DfltVarDimName')) };
    end
end