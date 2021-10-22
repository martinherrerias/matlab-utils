function [par,args,isdef] = getpairedoptions(args,varargin)
% [PAR,REST] = GETPAIREDOPTIONS(ARGS,DEFAULTS)
% [PAR,REST] = GETPAIREDOPTIONS(ARGS,NAMES,[DEFVAL])
% PAR = GETPAIREDOPTIONS(..,'restchk')
% PAR = GETPAIREDOPTIONS(..,'dealrest')
%
% Parse an input-argument list ARGS with parameter-value pairs.
%
% [PAR,REST] = GETPAIREDOPTIONS(ARGS,NAMES,DEFVAL) - Take the cell-array of arguments ARGS, and 
%   find all existing argument pairs of the form: {...,'par_i',val_i,...} for the parameter names
%   in NAMES = {'par_a','par_b'...}. If parameters are missing, and a cell-array of default values
%   DEFVAL = {def_a,def_b,...} is provided, assign to each its corresponding default value.
%
%   If DEFVAL is empty or missing, fields not found explicitly in ARGS will NOT be included!
%   A scalar DEFVAL will be expanded, e.g. {NaN}, or {[]} will place NaNs and empty values (resp.)
%   into any missing fields.
%
% [PAR,REST] = GETPAIREDOPTIONS(ARGS,DEFAULTS) - Do the same, but starting with structure of 
%   default options DEFAULTS. i.e. NAMES = fieldnames(DEFAULTS), DEFVAL = struct2cell(DEFAULTS).
%
% [PAR,REST,ISDEF] = GETPAIREDOPTIONS(ARGS,..) - Return a structure ISDEF similar to PAR, but with 
%   boolean values specifying whether each value comes from the defaults. NOTE that when ISDEF is
%   returned, PAR will be a full structure of parameters (even if no DEFVAL is provided) with 
%   empty values ([]) on otherwise missing fields.
%
% PAR = GETPAIREDOPTIONS(..,'restchk', [[MIN] MAX]) - assert that all but MIN to MAX ARGS have 
%   been cast into PAR, i.e. that after accounting for argument-value pairs, there are at least
%   MIN and no more than MAX arguments remaining. 'restchk' alone is equivalent to 'restchk', 0.
%   'restchk',N is equivalent to 'restchk',[0,N].
%
% PAR = GETPAIREDOPTIONS(..,'dealrest',M) - allow the first M NAMES to be provided as
%   positional arguments, allowing mixed syntaxes: ..,A,B,.. = ..,'a',A,'b',B,.. = ..,A,'b',B,..
%   An error will be raised if a non-empty argument is assigned to the same NAME both by a
%   name-value pair and by position. 'restchk' alone is equivalent to 'restchk',numel(NAMES).
%
% OUTPUT:
%   PAR: a structure with fields PAR.par_a = val_a, PAR.par_b = val_b ... for all found parameters,
%      and/or for all parameters with their corresponding defaults.
%   REST: a cell-array with any arguments that were not parsed into the PAR structure.
%
% TODO: check performance vs ARGUMENT block or INPUTPARSER!
%
% See also: CELL2STRUCT, GETFLAGOPTIONS, ADDPARAMETER

    narginchk(2,5);
    assert(iscell(args),'First argument must be a cell-array of arguments');
    
    checkrestchk = false;
    checkpositional = false;
    if isnumeric(varargin{end})
        nposargs = varargin{end}; 
        varargin(end) = []; 
    else
        nposargs = []; 
    end
    if ischar(varargin{end})
        switch varargin{end}
        case 'restchk'
            checkrestchk = true;
            if isempty(nposargs), nposargs = [0,0]; end
            if isscalar(nposargs), nposargs = [0,nposargs]; end
        case 'dealrest'
            checkpositional = true;
            if isempty(nposargs), nposargs = Inf; end
        otherwise
            error('Unknown option or flag: ''%s''',varargin{end}); 
        end
        varargin(end) = [];
    end
    
    nargin = numel(varargin)+1;
    assert(nargin >= 2,'Not enough input arguments');

    if isstruct(varargin{1})
    % Start with structure of defaults
        assert(nargin == 2,'Unrecognized arguments');
        par = varargin{1};
        names = fieldnames(par);
        n = numel(names);
        defval = true;
    else
        if nargin < 3, varargin{2} = {}; end % names, no default values
        [names,defval] = deal(varargin{:});
        assert(iscell(names) && iscell(defval),'Expecting a structure or two cell-arrays');
        n = numel(names);
        
        if isempty(defval)
        % No defaults provided, start with empty structure
            par = cell2struct(cell(n,1),names);
        else
        % Convert to structure of defaults
            if isscalar(defval) && n > 1, defval = repmat(defval,n,1); end
            assert(numel(defval) == n,'Expecting empty/N-cell-array of default values');
            par = cell2struct(defval(:),names(:));
        end
    end
        
    found = false(n,1);
    remargidx = 1:numel(args);
    for j = 1:n
    % Search each name in argument list...
        if isempty(args), break; end
        k = find(strcmpi(names{j},args));
        if isempty(k), continue; end
        
        found(j) = true;
        
        k(k == [0,k(1:end-1)+1]) = []; % resolve case ..'foo','foo'.. meaning par.foo = 'foo'
        if isscalar(k)
        % set par.name = val, remove ..,['name',val],.. from argument list
            if numel(args) < k+1
                ERR = MException('pairedopt:paired',...
                    'Named property "%s" has no matching value',names{j});
                throwAsCaller(ERR);
            end
            par.(names{j}) = args{k+1};
            args(k:k+1) = [];
            remargidx(k:k+1) = [];
        else
        % If the given option appears more than once ...
            c = args(k+1); % keep as cell
            args(k:k+1) = [];
            remargidx(k:k+1) = [];
            c = uniquecell(c);
            if isscalar(c)
            % ... and the option is just repeated, return single value
                warning('pairedopt:redundant','Redundant option: %s',names{j});
                par.(names{j}) = c{1};
            else
            % ... otherwise return cell-array of values
                warning('pairedopt:conflict','Conflicting option: %s',names{j});
                par.(names{j}) = c;
            end
        end
    end
    if nargout > 2
        isdef = cell2struct(num2cell(~found),names);
    elseif isempty(defval) && any(~found)
        par = rmfield(par,names(~found));
    end
    
    r = numel(args);
    if checkpositional && r > 0
        nposargs = min(n,nposargs);
        if r > nposargs
            ERR = MException('pairedopt:nposargs','Too many remaining positional arguments');
            throwAsCaller(ERR);
        end
        redundant = found(1:r)' & ~cellfun(@isempty,args);
        if any(redundant)
            ERR = MException('pairedopt:redundant',...
                shortliststr(names(redundant),'Redundant name-value/positional argument'));
            throwAsCaller(ERR);
        end
        if ~isequal(remargidx,1:r)
            warning('pairedopt:remargidx','Positional arguments appear after name-value pairs');
        end
        posargnames = names(~found);
        for j = 1:r
            par.(posargnames{j}) = args{j};
            isdef.(posargnames{j}) = false;
        end
    end

    if checkrestchk && (r < nposargs(1) || r > nposargs(2))
        msg = strtrim(regexprep(evalc('disp(args)'),'}\s*{',', '));
        ERR = MException('pairedopt:restchk','Unrecognized argument(s): %s',msg);
        throwAsCaller(ERR);
    end
end