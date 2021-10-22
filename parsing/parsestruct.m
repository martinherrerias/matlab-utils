function varargout = parsestruct(S,fields,varargin)
% PARSESTRUCT(S,FIELDS,..) - assert that S is a structure (or structure-like class), and that is
%   has fields/propperties FIELDS, with additional optional constraints (see below). Throw an error
%   (or issue just a warning*) otherwise.
%
% PARSESTRUCT(C,FIELDS,..) - with cell array C is (roughly) equivalent to a CELLFUN call to
%   VALIDATEATTRIBUTES, with the additional options described below. FIELDS are used to label 
%   arguments in warnings / errors. Use VALIDATEATTRIBUTES directly if possible!
%
% PARSESTRUCT(S,FIELDS,..,'X',Y,..) specify constraint(s), common to all FIELDS. Most arguments
%   are passed directly to VALIDATEATTRIBUTES. Exceptions are:
%
%  shorhands for common attributes:  -f = finite, -v = vector, -m = 2d, -s = scalar
%     -Z = nonzero, -p = nonnegative, -i = integer, -l = binary, -r = real.
%
%     'class',X - specify VALIDATEATTRIBUTES class(es). Default is a list of common 'primitives'
%                 {'numeric','logical','char','cell'}
%     -n is a shorthand for 'class','numeric'
%
%     -e, equalsize: check that the size of all fields is the same
%     -c, -compatible: check that fields are singleton-expansion-compatible with each other
%     'compatible',SZ: check that fields are compatible with each other AND with size SZ 
%      @f(x) - custom validation function handle(s), that return a scalar boolean
%
% PARSESTRUCT(..,'compatible',SZ) chech that fields are singleton-expansion compatible
%   with an array of size SZ (see COMPATIBLESIZE).
%
% PARSESTRUCT(..,'opt',OPTFIELDS,...) check that optional fields OPTFIELDS, if available,
%   comply with any requirements set by further flags/argument-pairs.
%
% P = PARSESTRUCT(..,'-soft') - instead of throwing an error, just issue a warning informing
%   which fields do not comply with their constraints. Return a structure P with those fields
%   removed.
%
% P = PARSESTRUCT(..,'-reduced') - By default, any fields in S that are not parsed as invalid
%   will be copied onto P. With '-reduced' any field that is not in FIELDS or OPTFIELDS will
%   not show on P.
%
% See also: INPUTPARSER, ARGUMENTS, VALIDATEATTRIBUTES, GETPAIREDOPTIONS, ISNESTEDFIELD
         
    narginchk(2,Inf);
    
    ALIAS = {'finite','-f';'numeric','-n';'vector','-v';'2d','-m';'scalar','-s';... 
             'nonzero','-Z';'nonnegative','-p';'integer','-i';'real','-r';...
             'binary','-l';'equalsize','-e';'compatible','c';'-compatible','-c'};
             
    FLAGS = {'-soft','-reduced','-compatible','equalsize','numeric'};
         
    narginchk(2,Inf);
    
    keys = cellfun(@ischar,varargin);
    [ia,ib] = ismember(varargin(keys),ALIAS(:,2));
    keys(keys) = ia;
    varargin(keys) = ALIAS(ib(ia),1);
    
    [opt,varargin] = getflagoptions(varargin,FLAGS);    
    opt.opt = {};
    if opt.numeric, opt.class = 'numeric'; else, opt.class = {'numeric','logical','char','cell'}; end
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    isfcn = cellfun(@(x) isa(x,'function_handle'),varargin);
    custom_fcn = varargin(isfcn);
    varargin(isfcn) = []; % remaining arguments for validateattributes
    
    if opt.soft, assert(nargout ~= 0,'-soft flag requires an output argument'); end
    if nargout < 1, opt.reduced = false; end
    
    if ~isempty(fields)
        if ischar(fields), fields = {fields}; 
        else
            assert(iscellstr(fields),'Expecting string/cell-array of strings FIELDS'); %#ok<ISCLSTR>
            fields = fields(:)';
        end   
    end

    if iscell(S)
        assert(numel(S) == numel(fields),'For cell input, numel(FIELDS) must match numel(C)');
        assert(isempty(opt.opt),'Cannot parse optional fields with cell input');
        opt.reduced = false;
        S = cell2struct(S(:),fields(:));
        report = @(msg) shortliststr(msg,[],'-newlines');
        FLD = 'argument';
        hasfields = @(f) true(size(f));
    else
        if isstruct(S)
            FLD = 'field';
            OBJ = 'structure';
            hasfields = @(f) isfield(S,f);
        else
            assert(isobject(S) && ~isempty(getfield(metaclass(S),{1},'PropertyList')),...
                'Expecting cell, structure or struct-like object');
            FLD = 'property';
            OBJ = 'object';
            hasfields = @(f) cellfun(@(f) isprop(S,f),f);
        end
        
        varname = inputname(1);
        if isempty(varname)
            report = @(msg) shortliststr(msg,[],'-newlines');
        else
            report = @(msg) shortliststr([{['In ' OBJ ' ' varname ':']}; msg],[],'-newlines');
        end
    end

    msg = {};
    if ~isempty(fields)
        missing = ~hasfields(fields);
        if any(missing)
            msg = {['Missing ' shortliststr(fields(missing),'field',Inf,'quotes','''')]};
            assert(opt.soft,report(msg));
            fields = fields(~missing);
        end
    end
    
    if ~isempty(opt.opt)
        assert(iscellstr(opt.opt),'Expecting string/cell-array of strings OPTFIELDS');
        opt.opt = opt.opt(:)';
            
        fields = [fields,opt.opt(hasfields(opt.opt))]; % id.
    end

    if opt.reduced
        assert(isstruct(S),'-reduced flag requires a structure');
        S = rmfield(S,setdiff(fieldnames(S),fields));
    end
    if isempty(fields)
        if ~isempty(msg), warning(report(msg)); end
        return; 
    end
    
    m = numel(custom_fcn);   
    n = numel(fields);
    bad = false(n,1);
    fldmsg = cell(n,1);
    
    for j = 1:n
        x = S.(fields{j});
        % x = getnestedfield(S,fields{j});
        try
            validateattributes(x,opt.class,varargin,'',fields{j});
            for k = 1:m
                assert(custom_fcn{k}(x),['failed to satisfy ' func2str(custom_fcn{k})]);
            end
        catch ERR
            bad(j) = true;
            fldmsg{j} = ERR.message;
        end
    end
    fldmsg = fldmsg(bad);
    
    givensize = isnumeric(opt.compatible);
    if givensize ||  opt.compatible || opt.equalsize
        
        sz = cellfun(@(f) size(S.(f)),fields,'unif',0);
        
        if givensize, sz{n+1} = opt.compatible; end
        d = max(cellfun(@numel,sz));
        for j = 1:numel(sz), sz{j}(end+1:d) = 1; end
        SZ = cat(1,sz{:});
        sz = max(SZ,[],1);

        if opt.equalsize
            failed = ~all(SZ == sz,'all');
        else
            failed = ~all(all(SZ == sz | SZ == 1,2) | all(SZ <= 1,2));
        end
        
        if failed
            bad(:) = true;
            fldmsg{end+1} = 'arguments have incompatible sizes';
        end
    end
    
    if ~any(bad)
        if ~isempty(msg), warning(report(msg)); end
        if nargout > 0, varargout{1} = S; end
        return; 
    end
        
    if opt.reduced, S = rmfield(S,fields(bad)); end
    if nargout > 0, varargout{1} = S; end

    msg{end+1} = [nthings(nnz(bad),FLD) ' failed to satisfy constraints:'];    
    msg =[ msg(:); fldmsg(:) ];

    if ~opt.soft
        err = MException('parsestruct:fail',report(msg));
        throwAsCaller(err);
    else
        warning(report(msg));
    end
end

function test()
%%
    S.x = magic(3);
    S.y = rand(10,1) > 0;
    S.z = @(x) 2;
    
    S.x2 = magic(4);
    S.x3 = magic(5);
    S.x4 = magic(6);
    
    tic()
        assert(isfield(S,'x'));
        assert(isfield(S,'y'));
        assert(isfield(S,'z'));
        validateattributes(S.x,{'numeric'},{'2d','real','finite','integer','positive'});
        validateattributes(S.y,{'numeric','logical'},{'vector','binary','size',[10,1]});
        validateattributes(S.z,{'function_handle'},{});
    toc()
    
    tic()
        parsestruct(S,{'x'},'numeric','2d','real','finite','integer','positive');
        parsestruct(S,{'y'},'class',{'numeric','logical'},'vector','binary','size',[10,1]);
        parsestruct(S,{'z'},'class','function_handle');
    toc()
    
    tic()
        assert(all(isfield(S,{'x','x2','x3','x4'})));
        for f = {'x','x2','x3','x4'}
            validateattributes(S.(f{1}),{'numeric'},{'2d','real','finite','integer','positive'});
        end
    toc()
    
    tic()
        parsestruct(S,{'x','x2','x3','x4'},'numeric','2d','real','finite','integer','positive');
    toc()

end
