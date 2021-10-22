function B = avgdownsample(A,m,varargin)
% B = AVGDOWNSAMPLE(A,m) - Take an array A, and average every m rows to get a floor(N/m)-array,
%   each element B(j,k) = mean(A(m·(j-1) + 1:m·j,k)).
%
% B = AVGDOWNSAMPLE(..,FLAG) for FLAG in {'omitnan','includenan','double','native','default'}
%   passes the corresponding FLAG(s) to the MEAN function. e.g. by default, elements B(j,k) for
%   which any A(m·(j-1) + 1:m·j,k) is NaN  will be NaN. Use 'omitnan' to override this behavior.
%
% B = AVGDOWNSAMPLE(A,m,W) - where W is a size(A,1) vector (numeric or logical), will return 
%  weighted means B(j,k) = sum(A(J,k).*w,FLAG)/sum(w,FLAG), where w = W(J).*isfinite(A(J,k)) 
%  and J = m·(j-1) + 1:m·j  
%
% R = AVGDOWNSAMPLE(S,m,N) - where S is a structure, applies AVGDOWNSAMPLE(A,m) to any fields
%   and sub-fields A of S which are N-row numeric arrays. In this case R will be a copy 
%   of the structure S, with (some) fields reduced to ~(N/m) rows.
%
% R = AVGDOWNSAMPLE(S,m,W) - where S is a structure, applies AVGDOWNSAMPLE(A,m,W) to any fields
%   and sub-fields A of S which are numel(W)-row numeric arrays.

% B = AVGDOWNSAMPLE(..,'full') return a ceil(N/m) array, where the last element is a partial mean
%   of mod(N,m) elements. Sets the 'omitnan' by default, unless 'includenan' is passed explicitly. 
%
% B = AVGDOWNSAMPLE(..,'offset',X) - applies a [fractional] row offset X before downsampling, i.e. 
%   shifts the rows of A [and W] by -fix(X), setting the first/last fix(X) elements to NaN, then 
%   gets the set means normally.
%   Fractional offsets F = rem(X,1) ~= 0, are resolved by averaging over m+1 elements (i.e. adding
%   an extra element to each set) and weighting the extreme values by F and (1-F), respectively.
%
% B = AVGDOWNSAMPLE(..,'logical','all'/'any'/X) - allow downsampling of logical arrays/fields, 
%   equivalent to L = AVGDOWNSAMPLE(double(L),..) > X for a numeric 0 < X < 1. Keywords 'all'
%   and 'any' are equivalent to X = 1-eps and X = 0, respectively. NaNs will result in False.
%   X = [] (default) will not filter logical fields!.
%
% TODO: add a 'smart' NaN FLAG that estimates the error of incomplete sample means (e.g. using
%   random bootstrapping of complete samples) and sets a minimal availability threshold to comply
%   with a given tolerance.
%
% TODO: should this just be a special case of RESAMPLESTRUCTURE?
%
% See also: NESTEDFIELDNAMES, FILTERSTRUCTURE, RESAMPLESTRUCTURE

    narginchk(2,7);
    FCNFLAGS = {'omitnan','includenan','double','native','default'};
    
    [opt,varargin] = getflagoptions(varargin,[FCNFLAGS,{'full'}]);
    opt.offset = 0;
    opt.logical = [];
    [opt,varargin] = getpairedoptions(varargin,opt);
    flags = FCNFLAGS(cellfun(@(f) opt.(f),FCNFLAGS));
    
    if ischar(opt.logical) && strcmpi(opt.logical,'all'), opt.logical = 1-eps(); end
    if ischar(opt.logical) && strcmpi(opt.logical,'any'), opt.logical = 0; end
    
    if isempty(opt.logical)
        isnumeric_ish = @(x) isnumeric(x) || isa(x,'datetime');
    else
        LT = opt.logical;
        assert(isnumeric(LT) && isscalar(LT) && LT >= 0 && LT < 1,'Bad logical threshold');
        isnumeric_ish = @(x) isnumeric(x) || islogical(x) || isa(x,'datetime');
    end
    
    if ~isempty(varargin)
        assert(numel(varargin) == 1,'Unrecognized arguments');
        W = varargin{1};
        if isscalar(W) && mod(W,1) == 0, N = W; W = true;  % AVGDOWNSAMPLE(S,m,N)
        else
            assert(isvector(W) && isnumeric(W) || islogical(W),...
                'Expecting scalar N or numeric/boolean vector W/F');
            N = numel(W);
            if isnumeric_ish(A), assert(N == size(A,1),'W/F must be a size(A,1) vector'); end
        end
    else
        if ~isnumeric_ish(A)
            if isstruct(A)
                error('AVGDOWNSAMPLE requires explicit N/W/F for structure input');
            elseif islogical(A)
                error('Logical arrays require explicit ''logical'',X argument pair');
            else
                error('Unrecognized argument A');
            end
        end
        W = true;
        N = size(A,1);
    end
    
    if opt.full && ~opt.includenan, opt.omitnan = true; end
    assert(~(opt.double && opt.native && opt.default) && ~(opt.includenan && opt.omitnan),...
        'Inconsistent flags');
    
    assert(abs(opt.offset) < N,'Offset cannot exceed the total number of rows');
    
    assert(isnumeric(m) && mod(m,1) == 0,'Expecting integer m');
    if m == 1 && opt.offset == 0, B = A; return; end
    
    % Numeric arrays:
    if isnumeric_ish(A)  
        B = arrayavg(A,m,W,opt.full,flags,opt.offset,opt.logical);
        return;
    end
    
    % Structures:
    assert(isstruct(A) && isscalar(A),'AVGDOWNSAMPLE takes only matrices or scalar structures');

    fields = nestedfieldnames(A);
    for f = 1:numel(fields)
        ff = strsplit(fields{f},'.');
        v = getfield(A,{1},ff{:});

        % Only for numeric fields (and, when restricted, N-row fields)...
        if ~isnumeric_ish(v), continue; end
        if size(v,1) ~= N, continue; end
        
        % AVGDOWNSAMPLE each sub-field
        A = setfield(A,{1},ff{:},arrayavg(v,m,W,opt.full,flags,opt.offset,opt.logical));
    end
    B = A;
end

function B = arrayavg(A,m,w,full,flags,offset,logic_threshold)

    waslogical = ~isempty(logic_threshold) && islogical(A);
    
    N = size(A,1);
    if ~ismatrix(A)
        sz = size(A);
        A = reshape(A,N,[]);
    else
        sz = [];
    end
    
    if abs(offset) >= 1
        A = circshift(A,-fix(offset),1);
        w = circshift(w,-fix(offset),1);
        if offset > 0
            A(end-fix(offset)+1:end,:) = NaN;
        else
            A(1:-fix(offset),:) = NaN;
        end
        offset = offset - fix(offset);
    end

    if islogical(w) && offset == 0
    % Set any ~w to NaN, let flags do the rest
        if ~all(w), A(~w,:) = NaN; end
        
        if mod(N,m) > 0
            if full, A(end+1:ceil(N/m)*m,:) = NaN;
            else, A(floor(N/m)*m+1:N,:) = [];
            end
        end

        A = reshape(A,m,[],size(A,2));
        B = mean(A,1,flags{:});
    else
    % Weighted mean
    
        if isscalar(w), w = repmat(w,N,1); end
        if mod(N,m) > 0
            if full
                A(end+1:ceil(N/m)*m,:) = NaN;
                w(end+1:ceil(N/m)*m) = 0;
            else
                A(floor(N/m)*m+1:N,:) = [];
                w(floor(N/m)*m+1:N) = [];
            end
        end

        w = ~isnan(A).*w;
        w = reshape(w,m,[],size(A,2));
        A = reshape(A,m,[],size(A,2));
        
        if offset > 0
        % add end+1 element as new row (extrapolate for last)
            A(m+1,:) = [A(m+1:m:end),NaN]';
            w(m+1,:) = [w(m+1:m:end),0]'*offset;  % add OFFSET weight to new row
            w(1,:,:) = w(1,:,:)*(1-offset);       % remove OFFSET weight from 1st
        elseif offset < 0
        % add -1 element as new row
            A(m+1,:) = [NaN,A(m:m:end-1)]';
            w(m+1,:) = [0,w(m:m:end-1)]'*(-offset);  % add -OFFSET weight to new row
            w(m,:,:) = w(m,:,:)*(1+offset);          % remove -OFFSET weight from last
        end
        
        B = sum(A.*w,1,flags{:})./sum(w,1,flags{:});
    end

    B = shiftdim(B,1);

    if ~isempty(sz)
        B = reshape(B,[size(B,1),sz(2:end)]);
    end
    
    if waslogical
        B = B > logic_threshold;
    end
end
