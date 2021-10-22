function [C,t2] = retimestruct(S,t,continuity,varargin)
% [C,tc] = RETIMESTRUCT(S,t,continuity,varargin) - provide timetable-like functionality
%   for nested structure S and time-step vector t. 
% EXAMPLE:
%   [C,t2] = retimestruct(S,t,'continuous','regular',method,'TimeStep',dt)
%
% See also: TIMETABLE/RETIME, SYNCHRONIZESTRUCT

    [TT,R] = struct2timetable(S,t);

    switch lower(continuity)
    case {'unset',[]}
    case {'continuous','step','event'}
        TT.Properties.VariableContinuity = repmat({continuity},1,size(TT,2));
        otherwise
        error('Expected {continuous, step, event, or unset}');
    end

    TT = retime(TT,varargin{:});
    [C,t2] = timetable2struct(TT);
    C = completestruct(C,R);
end