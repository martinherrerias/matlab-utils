function A = revertfilter(B,filter,dim,missing)
% A = REVERTFILTER(B,FILTER,DIM,MISSING) - Revert logical-indexing B = A(FILTER) upon dimension
%   DIM, filling any ~FILTER slices with MISSING.
%
%   E.g. if B = C(:,:,FILTER,:), A = REVERTFILTER(B,FILTER,3,NaN) will return a size(C) array,
%   such that all C(:,:,FILTER,:) == A(:,:,FILTER,:), and for which all A(:,:,~FILTER,:) == NaN.
%
% A = REVERTFILTER(B,FILTER,[],MISSING) - Revert logical-indexing B = A(FILTER), setting the size
%   of A from FILTER. For this to work DIM must be missing/empty, and nnz(FILTER) == numel(B).
%
% A = REVERTFILTER(B,IDX,SZ,MISSING) - Revert linear indexing B = A(IDX), where size(A) = SZ 
%
%   E.g. C = rand(1,2,3); FILTER = C > 0.5; B = C(FILTER); A = REVERTFILTER(B,FILTER)
%
% See also: MERGESTRUCTURES

    if nargin < 3, dim = []; end
    if nargin < 4
        try
            missing = filling(B);
        catch ERR
            B = num2cell(B);
            missing = filling(B);
            warning(ERR.identifyer,'%s',ERR.message)
        end
    end
    
    if isnumeric(B), missing = double(missing); 
    elseif islogical(B), missing = missing > 0;
    else
      assert(isequal(class(B),class(missing)),'Class of MISSING must match B');
    end

    if numel(dim) > 1 && isnumeric(filter) && all(filter > 0 & mod(filter,1) == 0)
    % A = REVERTFILTER(B,IDX,SZ,MISSING)
        A = repmat(missing,dim);
        A(filter) = B;
        return;
    end
    
    if isnumeric(filter) && all(filter == 0 | filter == 1), filter = filter == 1; end
    assert(islogical(filter),'REVERSEFILTER only works for logical indexing');

    if nnz(filter) == numel(B) && isempty(dim)
    % A = REVERTFILTER(B,FILTER,[],MISSING)
    
        A = repmat(missing,size(filter));
        A(filter) = B;
        
    elseif isvector(filter)
    % A = REVERTFILTER(B,FILTER,DIM,MISSING)
        
        if isempty(dim)
            if size(filter,2) > 1, dim = 2; else, dim = 1; end 
        end 
        assert(size(B,dim) == nnz(filter),'nnz(FILTER) ~= size(B,DIM)');
        sz = size(B);
        sz(dim) = numel(filter);
        A = repmat(missing,sz);
        
        % Place enough colons before and after dim for A(..,filter,..) = B to work
        args = [repmat({':'},1,dim-1),{filter},repmat({':'},1,numel(sz)-dim)];
        A(args{:}) = B;
    else
       error('Unrecognized syntax/arguments');
    end

end

function x0 = filling(x)
    if isnumeric(x), x0 = NaN(1,'like',x);
    elseif islogical(x), x0 = false;
    elseif isstruct(x)
        f = fieldnames(x);
        x0 = cell2struct(repmat({[]},numel(f),1),f);
    elseif iscell(x)
        try
            x0 = {eval([class(x{end}) '.empty'])};
        catch
            x0 = {[]};
        end
    elseif isdatetime(x)
        x0 = NaT(1,'TimeZone',x.TimeZone);
    else % struct / obj
        try
            x0 = eval(class(x)); % call class constructor
            assert(isscalar(x0));
        catch
            ERR = MException('revertfilter:missing',...
                'No default missing value for class: %s, provide one explicitly',class(x));
            throwAsCaller(ERR);
        end
    end
end