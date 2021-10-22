function [x,ix] = parselist(x,LIST,varargin)
% [L,IDX] = PARSELIST(x,LIST) - Make sure x is a (case insensitive) subset of cellstr LIST
%   Read x = 'all' as LIST, and x = {'^',..} as setdiff(LIST,x(2:end))
%   Return position indices idx such that x ~ L = LIST(IDX).
%   If any x are not found in LIST, throw an error: "Unknown: 'a','b',.."
%
%   LIST can be a "dictionary" with multiple aliases for each list element, listed as a 2-column 
%   cell of the form { label1,{alias11,alias12,..} ;label2, ...}, or as a structure with fields
%   LIST.label1 = {alias11,alias12,..}.
%
% ..,'-regexp' - read LIST aliases as regular expressions, e.g. {length,'^len.*';...}. If LIST
%   is a simple cellstr, -regexp will actually behave like -regexpx (below) assuming that regular
%   expressions are on x, and taking elements L from LIST.
% ..,'-regexpx' - read X as regular expressions, to be tested on LIST / LIST aliases.
% ..,'-matchcase' - case-sensitive match (ignorecase is the default!)
% ..,'-soft' - by default, an error is generated if any x is not found in LIST. With '-soft'
%   indices for unknown elements are 0, with unchanged L(idx == 0) = x(idx == 0)
%
% [L,IDX] = PARSELIST(x,LIST,'element','in my list of ...') - customize error message to:
%       "Unknown element(s): 'a','b',.. 'in my list of ...'"
%
% See also: VALIDATESTRING, RENAMEFIELDS, SHORTLISTSTR, NTHINGS

    [opt,varargin] = getflagoptions(varargin,{'-regexp','-regexpx','-matchcase','-soft'});
    
    if ~isempty(varargin)
        assert(iscellstr(varargin) && numel(varargin) <= 2,'Unrecognized arguments'); %#ok<ISCLSTR>
        varargin(end+1:2) = {''};
        [stuffname,sufix] = deal(varargin{:});
    else, stuffname = 'element'; sufix = '';
    end
        
    if isstruct(LIST)
        tgt = fieldnames(LIST);
        src = struct2cell(LIST);
    elseif iscell(LIST) && size(LIST,2) == 2 && iscellstr(LIST(:,1))
        tgt = LIST(:,1);
        src = LIST(:,2);
    elseif iscellstr(LIST) || isstring(LIST)
        src = LIST;
        tgt = LIST;
        if opt.regexp, opt.regexp = false; opt.regexpx = true; end
    end

    try
        tgt = cellstr(tgt);
        try
            src = cellstr(src);
            ir = (1:numel(src))';
        catch
            src = cellfun(@(x) reshape(cellstr(x),1,[]),src,'unif',0);
            ir = repelem(1:numel(src),cellfun(@numel,src))';
            src = cat(2,src{:})';
        end
    catch
        error('Invalid LIST or DICT');
    end
    
    singleinput = ischar(x);
    
    if ischar(x) && strcmpi(x,'all')
        ix = (1:numel(tgt))';
        x = tgt;
        return;
    elseif isempty(x), ix = []; return;
    elseif ischar(x), x = {x};
    end

    assert(iscellstr(x) || isstring(x),'Expecting numeric, char, or cellstr %s list',stuffname);
    
    inverted = strcmp(x{1},'^');
    if inverted, x(1) = []; end

    if isempty(src)
        if opt.soft
            ix = zeros(size(x));
            return;
        else
            found = false(size(x)); % let it crash below
        end
    else
        [L,src] = meshgrid(x,src);
        if opt.regexpx  
            if ~opt.matchcase, M = ~cellfun(@isempty,regexpi(src,L,'start'));
            else, M = ~cellfun(@isempty,regexp(src,L,'start'));
            end
        elseif opt.regexp  
            if ~opt.matchcase, M = ~cellfun(@isempty,regexpi(L,src,'start'));
            else, M = ~cellfun(@isempty,regexp(L,src,'start'));
            end
        else
            if ~opt.matchcase, M = strcmpi(L,src);
            else, M = strcmp(L,src);
            end
        end
        M = groupsummary(M,ir,'sum')' > 0; % group columns that point to same target

        ambiguous = sum(M,2) > 1;
        if ~opt.soft && any(ambiguous)
            error('parselist:ambiguous',['More than one LIST element match for ', ...
               shortliststr(x(ambiguous),'field')]);
        end
        found = any(M,2);
    end
    
    if ~opt.soft && ~all(found)
        ERR = MException('parselist:unknown',strjoin({shortliststr(x(~found), ...
                                  ['Unknown ',stuffname],'quotes','"','colon',':'),sufix},' '));
        throwAsCaller(ERR);
    end
    
    if ~inverted
        if ~any(ambiguous)
            ix = zeros(numel(x),1);
            [ib,~] = find(M');
            ix(found) = ib;
            x(found) = tgt(ib)';
        else
            ix = M;
            x(found) = arrayfun(@(j) tgt(M(j,:)>0),find(found),'unif',0);
        end
        if singleinput, x = x{:}; end
    else
        ix = find(~any(M,1));
        x = tgt(ix);
    end
end