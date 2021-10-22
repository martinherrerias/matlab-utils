function x = dms2deg(S)
% X = DMS2DEG(S) - Takes string(s) of the form '+00°00'00.00" N' and reads them as a 
% degree-minute-second series, returning their decimal equivalent. S can be a string,
% a char-array (where each row is a series) or a cell-array of strings.
% In the last cases the result will be a size(dms,1) or numel(dms) numeric vector.
% The function also works fine with decimal degrees, or degrees with decimal minutes,
% i.e. expresions like '0.00°N', or '0°00.5''' are not an issue.
% The separators between blocks of numbers can be any (group of) non-numeric character(s) 
% except for dots, commas, and the letters N,S,E,W (reserved for cardinal points).
% e.g. a string like ' -123 #30.0 $@  0*' will be read the same as -123°30'00".
% If the last letter matches W, or S the sign of the result will be inverted.
% The function is case insensitive, so DMS2DEG('52.1°w') is the same as DMS2DEG('52.1°W').
% Examples: DMS2DEG(' -123 #30.0 $@  0*') returns -123.5
%           DMS2DEG({'48°15´','48°15´W'}) returns [48.25;-48.25]
%           DMS2DEG('48.5°30´3600') returns 50.0 even though its just weird
% See also: DEG2DMS

    assert(ischar(S) || (iscell(S) && all(cellfun(@ischar,S))),...
            'Expecting string, character array, or cell-array of strings');
    if ~iscell(S), S = cellstr(S); end

    % Regular expresion patterns:
    num = '([\d.,]*)';      % floating point number
    sep = '[^\d.,NSEW]*';   % separator between blocks
    pat = ['\s*([+-]?)\s*',num,sep,num,sep,num,sep,'([NSEW]?)'];

    % Match strings to regexp pattern
    c = regexpi(S,pat,'tokens');

    c = cat(1,c{:}); c = cat(1,c{:}); % arrange into n·5 array of strings
    x = str2double(c(:,2:4)); % convert d,m,s to numbers, empties yield NaN

    % it's ok if minutes or seconds are empty, if degrees are empty result will be NaN
    x(isnan(x(:,2)),2) = 0;
    x(isnan(x(:,3)),3) = 0;

    signs = char(c(:,1)); % get signs from first token
    if isempty(signs), signs = ones(size(x,1),1); 
    else, signs = 1-(signs == '-')*2; end

    flip = char(c(:,5)); % revert signs for N, W
    if isempty(flip), flip = ones(size(x,1),1); 
    else, flip = 1-any(bsxfun(@eq,flip,'WSsw'),2)*2; end

    x = (x(:,1) + x(:,2)/60 + x(:,3)/3600).*signs.*flip;
end
