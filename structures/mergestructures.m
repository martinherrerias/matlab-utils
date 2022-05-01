function X = mergestructures(varargin)
% M = MERGESTRUCTURES(X,IDX,N) - outer join of structures X, treated as tables. Returns a merged
%   structure M such that (*) X{j} = FILTERSTRUCTURE(M,IDX{j},N).
%   X can be an cell array of dissimilar structures, or a regular structure array.
%   Each field X{j}.(k) that has N rows is placed into the IDX{j} rows of the matching M.(k)
%
% X = MERGESTRUCTURES(X,x,IDX,N,'-fillmissing') - fill any gaps in FILTER/IDX with NaN
%
% X = MERGESTRUCTURES(X,FILTER) - works the same, with logical indices FILTER{j}.
%
% S = MERGESTRUCTURES(s,..) - if s is a scalar structure, FILTER and IDX need not be cell arrays.
%   The result will be the inverse (*) of s = FILTERSTRUCTURE(S,..).
%
% X = MERGESTRUCTURES(A,B,C,..) treats all fields of structures A, B, C as non-separable (*).
%
% (*) NOTE: Unlike with FILTERSTRUCTURE, where non-separable fields are just cloned to the children
% structure, MERGESTRUCTURES needs to find a way to cope with these fields to avoid loss of info:
%
%   1. Fields of x that don't exist in X are cloned: if ~isfield(X,f), X.f <- x.f
%   2. Matching fields are preserved: if X.f == x.f, well, don't do anything
%   3. If x.f is empty, X.f will be left undisturbed
%   4. Classes should overload the method (e.g. STOPWATCH, SHADINGRESULTS)
%   5. Anything else will be piled up into cell-arrays, X.f = a, x.f = b -> X.f = {a,b}
%   ... and tested for uniqueness, i.e. X.f = {a,b}, x.f = {a,c} -> X.f = {a,b,c}
%
% NOTE that the criteria to decide which fields to filter (size(x,1) == n) is relatively weak, so 
% it might be a good idea to pass a partial structure with only the fields to be merged.
%
% TODO: handle tabular objects: call to outerjoin?
%
% See also: FILTERSTRUCTURE, NESTEDFIELDNAMES, AVGDOWNSAMPLE, REVERTFILTER

    if nargin > 1 && (isstruct(varargin{2}) || isobject(varargin{2})), varargin = {varargin}; end
    [varargin{end+1:4}] = deal([]);
    [S,idx,N,fillmissing] = deal(varargin{:}); clear varargin;

    assert(iscell(S) || isstruct(S) || isobject(S),'Expecting array of structures/objects');
    m = numel(S);

    if isempty(idx), idx = {}; end
    if isempty(N), N = NaN; end
    assert(isempty(idx) || (isnumeric(N) && isscalar(N)),'Expecting numeric scalar N');
    
    if ~isempty(fillmissing) && (strcmpi(fillmissing,'-fillmissing') || isequal(fillmissing,true))
        fillmissing = true;
    else
        fillmissing = false; 
    end
        
    if ~isempty(idx)
        if ~iscell(idx), idx = {idx}; end
        assert(numel(idx) == m,'Number of indices must match number of parts');
        if ~isequal(size(idx),size(S)), idx = idx(:); S = S(:); end

        % Parse all filters/indices, check for uniqueness and completeness
        [idx,N,missing] = parseindices(idx,N);
        if ~fillmissing, missing(:) = false; end
    else
        missing = [];
    end

    % Actual merging can start with a cell-array or a struct-like object, but both become
    % recursive (and might call each other) if contents require it.
    if iscell(S)
        X = mergearray(S,idx,N,missing);
    else
        X = mergefields(S,idx,N,missing); 
    end
end

function V = mergearray(v,idx,N,missing,V)

    if nargin < 5, V = []; end
    m = numel(v);
    
    divisible = ~isempty(idx) && all(arrayfun(@(x,n) size(x{1},1) == n,v,cellfun(@numel,idx)));
    similar = allsimilar(v);

    if divisible && similar && ~ischar(v)
    % Matching (i.e. divisible) field
    
        if isempty(V)
        % Expand singleton/empty dimensions of V to avoid crash
            V = eval([class(v{1}) '.empty']);
            if isdatetime(V),V.TimeZone = v{1}.TimeZone; end
            s = size(v{1}); s(1) = 0;
            V = reshape(V,s);
        end

        % Replace the affected rows, and put everything back into place.
        for i = 1:m, V(idx{i},:) = v{i}(:,:); end

        if any(missing)
            V(missing,:) = filling(v{1});
        elseif size(V,1) < N
            V(end+1:N,:) = filling(v{1});  % In any case ensure size(V,1) == N
        end
        
    elseif similar && (isstruct(v{1}) || (isobject(v{1}) && hasinternalmerge(v{1}))) && ...
                all(cellfun(@(x) isequal(size(x),size(v{1})),v))
    % [Equal-sized array of] structure(s): recursive merge subfields

        sz = size(v{1});
        if prod(sz) > 1, v = cellfun(@(v) v(:),v,'unif',0); end
        v = cat(2,v{:});
        
        if isobject(v(1)) && hasinternalmerge(v(1))
           V = mergestructures(v,idx,N,any(missing));
           return;
        end
        
        if isempty(V)
            V = filling(v(1));
            V(prod(sz)) = V; 
            V = reshape(V,sz); 
        end

        for j = 1:size(v,1)
            V(j) = mergefields(v(j,:)',idx,N,missing);
        end
    else
    % Otherwise concatenate in a cell-array (removing duplicates)

        if iscellstr(v{1}) && similar
            v = cat(1,v{:});
        end
        V = uniquecell(v,'stable');
        if numel(V) == 1, V = V{:}; end
    end
end

function X = mergefields(S,idx,N,fillmissing)

    fields = fieldnames(S);
    isthere = @(S,x) isfield(S,x) || isprop(S,x);
    
    % For non-structure objects, take care not to try to get/set something we shouldn't
    if ~isstruct(S), fields = filterfields(fields,metaclass(S)); end
    
    % Allocate an empty merging structure with the same fields as S(1)
    X = filling(S);
    
    for k = 1:numel(fields)
        % get field k from X, if existing*
        if isthere(X,fields{k}), V = X.(fields{k}); else, V = []; end
        
        % get corresponding fields k from S
        v = {S.(fields{k})}';
        
        X.(fields{k}) = mergearray(v,idx,N,fillmissing,V);
    end
end

function a = allsimilar(x)

    if isempty(x), a = true; return; end

    c = class(x{1});
    a = all(cellfun(@(y) isequal(class(y),c),x));
    if ~a, return; end
    
    if ~isscalar(x{1})
        s = size(x{1});
        a = all(cellfun(@(y) isequal(s(2:end),arrayfun(@(j) size(y,j),2:numel(s))),x));
        if ~a, return; end
    end
    
    if isstruct(x{1})
        f = fieldnames(x{1});
        a = all(cellfun(@(y) isempty(setxor(f,fieldnames(y))),x));
        if ~a, return; end
    end
end

function a = hasinternalmerge(x)
    MC = metaclass(x);
    a = any(strcmp('mergestructures',{MC.MethodList.Name}));
end

function x0 = filling(x)
    if isnumeric(x), x0 = NaN;
    elseif islogical(x), x0 = false;
    elseif isstruct(x)
        f = fieldnames(x);
        x0 = cell2struct(repmat({[]},numel(f),1),f);
    elseif iscell(x), x0 = {};
    elseif isdatetime(x)
        x0 = NaT(1,'TimeZone',x.TimeZone);
    else % struct / obj
       x0 = eval(class(x)); % call class constructor
    end
end

function fields = filterfields(fields,mc)
% FIELDS = FILTERFIELDS(NESTEDFIELDNAMES(X),METACLASS(X))
% Take a cellstr of NESTEDFIELDNAMES(X), and remove any that cannot be accessed or set, according
% to the PropertyList in METACLASS(X). 

    p = mc.PropertyList;

    unsettable = [p.NonCopyable] | [p.Abstract] | [p.Constant] | [p.Dependent];
    restricted = ~strcmp({p.GetAccess},'public') | ~strcmp({p.SetAccess},'public');
    restricted = restricted & ~unsettable;
    p = {p.Name};

    if any(unsettable)
        warning('mergein:unsettable',... 
            'Abstract/Constant/Dependent/NonCopyable %s will not be merged',...
            shortliststr(p(unsettable),{'property','properties'},'quotes',''''));
    end
    if any(restricted)
        warning('mergein:restricted',... 
            'Set/GetAccess is restricted for %s, and will not be merged',...
            shortliststr(p(restricted),{'property','properties'},'quotes',''''));
    end
    p = p(~unsettable & ~restricted);
    if isempty(p), return; end

    q = regexprep(fields,'([^.]*).*','$1');  % property to which nested fields belong
    fields(~ismember(q,p)) = [];
end

function [idx,N,missing] = parseindices(idx,N)
% Parse all filters/indices (turn logical to numeric), get/check N, check for uniqueness and
% completeness.

    m = numel(idx);
    N = repmat(N,1,m);
    for j = 1:numel(idx)
        if islogical(idx{j})
            if isnan(N(j)), N(j) = numel(idx{j});
            else
                assert(numel(idx{j}) == N(j),'Inconsistent N and filter-size',j);
            end
            idx{j} = find(full(idx{j}(:)));
        else
            if isnan(N(j)), N(j) = max(idx{j});
            else
                assert(N(j) >= max(idx{j}),'Inconsistent N(%d) and IDX',j);
            end
            assert(all(idx{j} > 0) && all(mod(idx{j},1) == 0),'Bad index (%d)',j);
            idx{j} = idx{j}(:);
        end
    end
    N = max(N,[],'omitnan');
    
    used = accumarray(cat(1,idx{:}),1,[N,1]);
    assert(all(used <= 1),'Indices are not unique');
    
    missing = ~used;
end
