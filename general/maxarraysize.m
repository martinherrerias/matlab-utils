function n = maxarraysize(type)
% N = MAXARRAYSIZE(TYPE) - return the maximum number of elements of an array of class TYPE based 
%   on available memory.

    if nargin == 0 || isempty(type), type = 'double'; end
    if isa(type,'class'), type = type.Name; end
    assert(ischar(type),'Expecting MetaClass Object or class-name (char)');
    switch lower(type)
        case {'double','int64','uint64'}, s = 64;
        case {'single','int32','uint32'}, s = 32;
        case {'int16','uint16'}, s = 16;
        case {'int8','unit8','char'}, s = 8;
        case {'logical'}, s = 1;
        otherwise
            s = 1;
            warning('Unrecognized primitive, returning free memory in bytes');
    end

    try
        [userview, ~] = memory();
        n = floor(userview.MaxPossibleArrayBytes/s);
        return;
    catch ERR
        if ~strcmp(ERR.identifier,'MATLAB:memory:unsupported'), rethrow(ERR); end
    end

    [~,w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    % stats = stats(1)/1e6;
    n = floor((stats(3))*8e3/s);
end