function [S,t,P] = timetable2struct(TT,names,timename)
% [S,T] = TIMETABLE2STRUCT(TT,NAMES) - Convert TIMETABLE object TT into a (nested) structure 
%   with fieldnames NAMES, returning RowTimes in T.
%
% [S,T] = TIMETABLE2STRUCT(TT) - Use TT.Properties.VariableDescriptions as NAMES (e.g. as a
%   result of STRUCT2TIMETABLE), or TT.Properties.VariableNames when this fails.
%
% S = TIMETABLE2STRUCT(TT,[..],TNAME) - return RowTimes as an additional field S.(TNAME). 
%   If omitted, TNAME = TT.Properties.DimensionNames(1)
%
% NOTE: time-stamps are starts-of-interval, as is always the case with TIMETABLE objects.
%
% See also: STRUCT2TIMETABLE

    validnestedname = @(n) ~isempty(n) && all(cellfun(@isvarname,strsplit(n,'.')));
    if nargin < 2 || isempty(names)
        names = TT.Properties.VariableDescriptions';
        if numel(names) ~= size(TT,2) || ~all(cellfun(validnestedname,names))
            names = TT.Properties.VariableNames;
        end
    else
        names = cellstr(names);
        assert(size(TT,2) == numel(names),'Non matching TT & NAMES');
    end
    if nargin < 3 || isempty(timename), timename = TT.Properties.DimensionNames(1); end
    
    t = TT.Properties.RowTimes;
    P = TT.Properties;
    
    validnames = cellfun(validnestedname,names);
    if all(validnames)
        warning_disabler = naptime('MATLAB:table:ModifiedVarnames'); %#ok<NASGU>
    end
    TT = table2struct(TT,'ToScalar',true); % include warning for invalid names
    
    if ~all(validnames)
        fld = fieldnames(TT);
        names(~validnames) = fld(~validnames);
    end

    C = struct2cell(TT);
    if nargout < 2
        names(end+1) = timename;
        C{end+1} = t;
    end
    S = cell2nestedstruct(C(:),names(:) );
end