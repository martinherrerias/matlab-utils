function obj = naptime(warnIDs,state)
% OBJ = NAPTIME(WARNID) - turn off warning(s) WARNID until ONCLEANUP OBJ is destroyed.
% OBJ = NAPTIME(WARNID,STATE) - set warning(STATUS,WARNID) until OBJ is destroyed.
% OBJ = NAPTIME() - don't change anything (yet) just reset warnings when OBJ is destroyed

    if nargin < 2, state = 'off'; end
    
    assert(nargout > 0,'Cleanup object must be returned as output');
    
    warning_status = warning();
    
    if nargin > 0 && ~isempty(warnIDs)
        validateattributes(warnIDs,{'char','cell','string'},{'nonempty'},'','Warning ID(s)');
        cellfun(@(ID) warning(state,ID),cellstr(warnIDs));
    end
    
    obj = onCleanup(@() warning(warning_status));
end
