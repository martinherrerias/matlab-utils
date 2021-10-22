function B = filterstructure(A,filter,varargin)
% B = FILTERSTRUCTURE(A,FILTER) - Take stucture A, and apply a logical FILTER to any and all 
% fields and  subfields which have N rows, where N = numel(FILTER). In other words, for each 
% field k in NESTEDFIELDNAMES(A): if size(A.k,1) == N, B.k = A.k(FILTER), otherwise B.k = A.k.
% NOTE that A.k above can stand for a potentially nested field, e.g. A.b.c.d...
%
% B = FILTERSTRUCTURE(A,IDX,N) - Do the same, but using numerical indexing IDX and explicitly
% defining number of rows N. (see MERGESTRUCTURES for inverse-indexing).
%
% B = FILTERSTRUCTURE(A,M) - where M is a matrix, uses matrix-multiplication B.k = M·A.k for any
%   fields A.k with size(A.k,1) == size(M,2), thus allowing for resampling, averaging, etc.
%   See INTERPMATRIX for generation of interpolation matrices.
%
% B = FILTERSTRUCTURE(..,'dim',K) - Apply filter along the K'th dimension. (Experimental) 
%
% Y = FILTERSTRUCTURE(X,...) - Works also for TABLE, TIMETABLE, and DATETIME objects, as well
%   as for simple numeric and logical arrays. For other types, function should be overloaded 
%   inside the class (implementing a valid struct call if required).
%
% NOTE that the criteria to decide which fields to filter (size(x,1) == N) is relatively weak, so 
% it might be a good idea to pass a partial structure with only the fields to be filtered.
%
% See also: MERGESTRUCTURES, NESTEDFIELDNAMES, AVGDOWNSAMPLE, INTERPMATRIX

    if ~isempty(varargin)
        [opt,varargin] = getpairedoptions(varargin,{'dim','n'},{1,[]});
        if ~isempty(varargin)
            assert(numel(varargin) == 1,'Unrecognized argument(s)');
            opt.n = varargin{:};
        end
        DIM = opt.dim;
        N = opt.n;
        assert(mod(DIM,1) == 0 && DIM > 0,'Expecting integer dimension > 0');
        assert(isempty(N) || mod(N,1) == 0 && N > 0,'Expecting integer N > 0');
    else
        N = [];
        DIM = 1;
    end
    matrixfilter = size(filter,1) > 1 && size(filter,2) > 1;
    
    if matrixfilter
    % B = FILTERSTRUCTURE(A,M)
        if isempty(N), N = size(filter,2);
        elseif N ~= size(filter,2), filter = ''; % crash below
        end
        assert(DIM == 1,'Matrix filters do not yet work with DIM > 1');
        assert(isnumeric(filter) && ismatrix(filter),'Expecting N-vector or M·N matrix filter');
    elseif isnumeric(filter)
    % B = FILTERSTRUCTURE(A,IDX,N)
        if isempty(N), N = max(filter); end
        assert(max(filter) <= N && all(mod(filter,1)==0) && min(filter) > 0,'Index out of bounds')
        % assert(all(mod(filter,1)==0) && min(filter) > 0,'Index out of bounds')
    elseif islogical(filter)
    % B = FILTERSTRUCTURE(A,FILTER)
        if isempty(N), N = numel(filter); end
        assert(isvector(filter) && numel(filter) == N,'Expecting logical FILTER to be N,1 vector');
    else
        error('Expecting logical vector, array of indices, or M·N matrix filter');
    end
    
    if isnumeric(A) || isa(A,'tabular') || isdatetime(A)
        assert(size(A,DIM) == N,'filter does not match A');
    end
    
    if isnumeric(A) || isdatetime(A)
        B = filterarray(A,filter,DIM);
    elseif ~matrixfilter && ~takesrecursivefilter(A) && size(A,DIM) == N
        B = filterarray(A,filter,DIM); % simple indexing
    else
        if istable(A)
            P = A.Properties;
            A = table2struct(A,'toscalar',true);
        elseif istimetable(A)
            P = A.Properties;
            A = timetable2struct(A);
        else
            P = [];
            assert(isstruct(A),'FILTERSTRUCT(OBJ) not defined for %s',metaclass(A));
        end
            
        B = recursivefilter(A,filter,N,DIM);
        
        if isempty(P), return; end
        
        switch class(P)
        case 'matlab.tabular.TimetableProperties', B = struct2timetable(B,P.DimensionNames{1});
        case 'matlab.tabular.TableProperties', B = struct2table(B);
        end
        
        for f = {'Description','UserData','DimensionNames','VariableNames',...
                'VariableDescriptions','VariableUnits','VariableContinuity','CustomProperties'}
            B.Properties.(f{1}) = P.(f{1});
        end
    end
end

function x = filterarray(v,filter,DIM)

    matrixfilter = size(filter,1) > 1 && size(filter,2) > 1;
    
    if matrixfilter
        z = size(filter,1); % New size(1)
    end  
    
    s = size(v);
    if ~matrixfilter
        idxargs = repmat({':'},1,numel(s));
        idxargs{DIM} = filter;
        x = v(idxargs{:});
    else
        istime = isdatetime(v);
        if istime
            t0 = v(1); % preserves timezone!
            v = days(v-t0); 
        end
        if issparse(filter) && ~isa(v,'double')
            oldclass = class(v);
            v = filter*reshape(double(v),s(1),[]);
            switch oldclass
                case 'logical', v(~isfinite(v)) = false;
                % others?
            end
            v = cast(v,oldclass);
        else
            v = filter*reshape(v,s(1),[]);
        end
        if istime, v = t0 + days(v); end
        x = reshape(v,[z,s(2:end)]);
    end
end

function B = recursivefilter(B,filter,N,DIM)
    
    fields = fieldnames(B);
    for k = 1:numel(B)
        for f = 1:numel(fields)
            if takesrecursivefilter(B(k).(fields{f}))
            % Go down one level, filter whatever is on B(k).(fields{f})
                B(k).(fields{f}) = recursivefilter(B(k).(fields{f}),filter,N,DIM);
                
            elseif size(B(k).(fields{f}),DIM) == N
            % Apply filter if field has length N
                B(k).(fields{f}) = filterarray(B(k).(fields{f}),filter,DIM);
            end
        end
    end
end

function a = takesrecursivefilter(x)
    % if isstruct(x), a = true; return; end
    a = isstruct(x) || isa(x,'ShadingResults');
    
    % % For some reason this takes ages:
    % MC = metaclass(x);
    % methods = {MC.MethodList.Name};
    % a = any(strcmp('filterstructure',methods));
end
