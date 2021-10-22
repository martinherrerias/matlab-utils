function [par,args] = getflagoptions(args,names,restchk)
% [PAR,REM] = GETFLAGOPTIONS(ARGS,NAMES) - Parse an input-argument list ARGS in search of flag
%   keys included in NAMES. Return a boolean structure PAR indicating whether each possible flag
%   was found, along with the reduced (i.e. not matching) set of arguments REM.
%
% PAR = GETFLAGOPTIONS(..,'restchk') - assert that all ARGS have been cast into PAR, i.e. that
%   all arguments have been recognized as known flags.
%
% INPUT:
%   ARGS: cell array of arguments to be parsed, any type. Only char- and string-elements will be
%       compared with NAMES.
%   NAMES: A) Single flag-name, e.g. '-opt' (must match exactly one of args{j})
%          B) Cell-string of different flag-names, e.g. {'-foo','-bar','-baz',...}
%          C) Cell of cell-strings of flag-aliases, where the first element is the desired field
%             name, e.g. {{'--foo','-f'},{'--bar','-b'},{'--baz','-z','-B'}}
%
% OUTPUT:
%   PAR: a structure with fields PAR.x = true if flag x was found among ARGS, PAR.x = false 
%       otherwise, for all x in NAMES*.
%       (*) Preceeding dashes will be removed, and names passed to matlab.lang.makeValidName.
%           The first name on every list of aliases will be used.
%       
%   REM: a cell-array with any arguments that were not parsed into the PAR structure.
%
% EXAMPLE: [opt,rem] = GETFLAGOPTIONS(({3, nan,'-foo','-a'},{'-a','-b'})
%          Returns opt.a = 1, opt.b = 0; rem = {3, NaN, '-foo'}
%
% See also: GETPAIREDOPTIONS

    narginchk(2,3);
    if nargin > 2
        assert(ischar(restchk) && strcmpi(restchk,'restchk'),'Unknown argument/flag');
        checkrestchk = true;
    else
        checkrestchk = false;
    end
    
    candidates = cellfun(@(x) ischar(x) || isstring(x),args);
    used = false(size(candidates));
    
    if isempty(names), par = struct(); return; end
    
    if ~iscell(names), names = cellstr(names); end
    
    matched = false(numel(names),1);
    keys = cell(numel(names),1);
    for j = 1:numel(names)    
        if ischar(names{j})
            keys(j) = names(j);
            isflagj = strcmpi(args(candidates),names{j});
        else
            keys(j) = names{j}(1);
            isflagj = cellfun(@(x) any(strcmpi(names{j},x)),args(candidates));
        end
        if any(isflagj)
            matched(j) = true;
            used(candidates) = used(candidates) | isflagj;
            candidates(candidates) = ~isflagj;
        end
    end
        
    % ... then reduce names to reduced* first element, e.g. {'--some flag','-s'} -> 'some_flag'
    keys = regexprep(keys,'^-*(.*)','$1');
    keys = matlab.lang.makeValidName(keys(:));

    args = args(~used); % remove found flags from argument list
    par = cell2struct(num2cell(matched),keys);
    
    if checkrestchk && ~isempty(args)
        evalin('caller','error(''Unrecognized argument-value pair(s)'')');
    end
end
