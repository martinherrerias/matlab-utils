function [R,varargout] = runparallel(fcn,args,sidx,varargin)
% R = RUNPARALLEL(FCN,ARGS,SIDX,[M]) - calculate the result of the single-output function
%   FCN(ARGS{:}) by splitting a subset SIDX of its arguments, i.e. ARGS(SIDX) into M partial copies
%   each, packing them (along with a copy of any non-separable arguments) into self-contained -mat 
%   files, running partial computations in batch mode (using PACKEDFCN), and compiling the results 
%   with MERGERESULTS.
%
% R = RUNPARALLEL(FCN,ARGS,SIDX,[M,P]) - split the function into M partial copies, but run only
%   P batches at a time.
%
% INPUT:
%   FCN - function handle, for the function to be calculated. Expected to produce a single output
%       that is row-wise independent for each row of the arguments ARGS(SIDX).
%   ARGS - cell array of arguments for FCN. The equivalent is R = FCN(ARGS{:})
%     Separable arguments (indexed by SIDX) are to be row-split by FILTERSTRUCTURE, i.e. recognized
%     types are numerical arrays, (nested) structures, tables, and any user-defined class (like
%     SHADINGRESULTS) that overloads the FILTERSTRUCTURE method. 
%     Constant arguments (any ARGS no indexed by SIDX) are passed as identical copies to all groups.
%   SIDX - logical or integer index for separable arguments (subset of ARGS).
%   M - number of partial copies into which the arguments should be split. The default is set by
%       ceil(N/SimOptions.chunksize), where N is the number of rows of separable arguments.
%   P - number of processes to run simultaneously, it defaults to M, unless limited by SimOptions 
%       or available memmory. See GETNWORKERS for details.  
%
% R = RUNPARALLEL(...,'-simoptions') - includes a copy of the global SimOptions in each partial
%   file, so that any call to GETSIMOPTION within FCN is consistent* with a single-thread call of
%   FCN(ARGS{:}) in the current session. Minor exceptions (logs, version, etc.) are listed in
%   PARALLELCONFIG.option_overrides.
%
% R = RUNPARALLEL(..,'-backup') - include autosave-crash-resume functionality. For this to work, 
%   FCN must recognize the name-value pair ..,'backup',X. Specifically the calls:
%
%       FCN(args{:},'backup',FILE~) - where FILE~ is an autogenerated backup file name
%       FCN('backup',S) - where S = load(FILE~,'-mat') contains the contents of an interrupted run
%
%   On a fresh run (first syntax), FCN is informed that is should start dumping any provisional 
%   results into a MAT file named FILE~. On a recovery run (i.e. if FILE~ exists), FCN is expected
%   to recover and resume its task from the contents of S.
%
% R = RUNPARALLEL(..,'N',N) - for complex input arguments, it might not always be straightforward
%   to decide what a 'row' is, and how to divide the arguments. Providing an explicit N will allow
%   better parsing of structures with uneven fields.
%
% TODO: + take 'dim',D to slice accross other dimensions
%       + use varargout to handle multiple outputs
%
% OUTPUT:
%   R - result compiled by MERGERESULTS, and ultimately by MERGESTRUCTURES. NOTE that there is no
%   unequivocal way to handle non-separable fields & properties of many classes, so for complex
%   result structures it might be a good idea to define a results class and overload this last
%   method (see e.g. SHADINGRESULTS.MERGERESULTS).
% 
% [R,CLEANER] = RUNPARALLEL(...) - When RUNPARALLEL runs without problems, all provisional
%   files and batch jobs are deleted after results are merged. To delay this deletion, e.g. until
%   merged results have been saved to disk, function handle CLEANER can be passed around and 
%   called only once the results are safe.
%
% TEMPORARY FILES: (See naming conventions in PARALLELCONFIG).
%
% See also: PACKEDFCN, PARALLELCONFIG, MERGERESULTS
    
    if debugging() && nargin == 0
        test(); 
        return; 
    end
    narginchk(3,9)
    
    script = which('runparallel.sh');
    assert(~isempty(script),'Failed to find runparallel.sh');

    if ischar(fcn) || isstring(fcn), fcn = str2func(fcn); end
    assert(isa(fcn,'function_handle'),'FCN must be a function handle');
    
    assert(iscell(args) && isvector(args),'Expecting cell-vector of arguments ARGS');
    
    [opt,varargin] = getflagoptions(varargin,{'-simoptions','-backup'});
    opt.N = [];
    [opt,varargin] = getpairedoptions(varargin,opt);
    if ~isempty(varargin)
        if isscalar(varargin), varargin{2} = NaN; end
        [varargin{cellfun(@isempty,varargin)}] = deal(NaN);
        assert(numel(varargin) == 2 && all(cellfun(@(x) isscalar(x) && isnumeric(x) && ...
            (isnan(x) || mod(x,1) == 0),varargin)),'Unrecognized argument(s)');
        [M,P] = deal(varargin{1:2});
       if isempty(P), P = NaN; end       
    else
        P = NaN; M = NaN;
    end

    narg = numel(args);
    separable = false(1,narg);
    try separable(sidx) = true;
    catch
        error('Expecting logical or integer "separable" index SIDX for %d argument array',narg);
    end

    % Actual work: divide arguments into M PACKEDFCN-friendly files
    [POPT,partfiles,~,P] = divideargs(fcn,args,separable,M,P,opt);
    clear args
    
    subdir = fullfile(POPT.path,POPT.subdir);
    idxfile = fullfile(subdir,[POPT.prefix POPT.index]);
    
    % Change to .par directory, fall back to current directory on crash
    basedir = pwd(); lastwill = onCleanup(@() cd(basedir)); cd(subdir);

    % savepath('pathdef.m'); % doesn't do what it's supposed to on Hazelhen,
                             % runparallel.sh must call matlabsetup.sh to fix the issue
                             
    % Check completion
    finished = cellfun(@(f) ismember('res',who('-file',f)),partfiles);
    if any(finished)
        fprintf('\t(Skipping %d completed jobs)\n',nnz(finished));
    end
    
    % Skip finished (e.g. resuming a crashed simulation)
    finished = cellfun(@(f) ismember('res',who('-file',f)),partfiles);
    partfiles = partfiles(~finished);
    
    if ~all(finished)
        % call to RUNPARALLEL script
        % Run from shell: $ ./runparallel.sh -p P part1 part2 ...
        bashcall = strjoin(cellfun(@(s) ['"' s '"'],partfiles,'unif',0),' ');
        bashcall = [script sprintf(' -p %d ',P) bashcall];
        system(bashcall);
    end
    
    finished = cellfun(@(f) ismember('res',who('-file',f)),partfiles);
    assert(all(finished),'%d of %d parts terminated on error',nnz(~finished),numel(partfiles));
    
    varname = @(x) inputname(1);  % (clearer var tracing, avoids NASGU/ASGLU warnings)
    clear(varname(lastwill));

    [R,cleaner] = mergeresults(idxfile);
    if nargout < 2, cleaner(); 
    else, varargout = {cleaner}; 
    end
end

function [POPT,parts,idx,P] = divideargs(fcn,args,separable,M,P,opt)
% [POPT,PARTS,IDX,P] = DIVIDEARGS(FCN,ARGS,SEPARABLE,M,P,OPT) - take the separable
%   arguments ARGS(SEPARABLE) and split them row-wise(*) into M partial copies; then put these
%   together with any remaining ARGS(~SEPARABLE) into M cell-arrays, and save them as 'args' 
%   into -mat files, along with an 'idx' specifying which rows of the original arguments went
%   into the -mat file, a function handle 'fcn' for FCN, and if OPT.simoptions, a (modified) copy 
%   of the global SimOptions.
%
%   If OPT.backup, the list of ARGS is appended with ..,'backup',BACKUP_FILE_NAME, to let FCN know,
%   upon runtime, where it is expected to dump its backups.
%
% OUTPUT:
%   POPT - PARALLELCONFIG output, based on OPT.partfilename (if available) or 'prjname_fcn'. It
%       defines PATH, SUBDIR, PREFIX, and extensions PART and IDX used below.
%       e.g. the index-file is F = fullfile(POPT.path,POPT.subdir,[POPT.prefix POPT.index])
%
%   PARTS - M cellstr of -mat file names where partial-argument copies were stored. File names
%       follow the convention ./PREFIX_XX.PARTS. Note(!) PARTS are relative paths from PATH/SUBDIR.
%
%   IDX - M cell array of vectors, indicating which rows of the original separable arguments 
%       t1,t2,..tN are stored on each partial copy.
%
%   P - parsed number of workers, limited by memmory and SimOptions.
%
% WRITTEN FILES (See naming conventions in PARALLELCONFIG).
%   PATH/SUBDIR/PREFIX_XX.PARTS (for XX in 1..M) containing a cell-array of partial arguments and
%       a vector of indices indicating which rows of the original separable variables are included
%       in the given file.
%   PATH/SUBDIR/PREFIX.IDX - a single file containing this function's outputs SUBS, IDX.

    MAX_P = 36; % avoid 1000 MATLAB instances, even if user sets maxthreads = 1000

    global SimOptions;
    if isempty(SimOptions)
       SimOptions = struct();
       if opt.simoptions
          warning('No SimOptions structure to store');
          opt.simoptions = false;
       end
    end
    
    Nt = checkseparable(args(separable),opt.N);
    
    fcn_str = matlab.lang.makeValidName(func2str(fcn));
    if  isfield(SimOptions,'prjname') && ~isempty(SimOptions.prjname)
        POPT = parallelconfig([SimOptions.prjname '_' fcn_str]);
    else
        POPT = parallelconfig(fcn_str); % naming conventions and other Parallel-Options
    end
    
    % Override options included as last argument
    if opt.simoptions
        opts = completestruct(POPT.option_overrides,SimOptions);
    else
        opts = struct();
    end
    subpath = fullfile(POPT.path,POPT.subdir);
    
    if isempty(dir(subpath)), mkdir(subpath); end
    
    % Change to project directory, fall back to current directory on crash
    basedir = pwd(); cleaner = onCleanup(@() cd(basedir)); cd(subpath);
    
    idxfile = fullfile([POPT.prefix POPT.index]);
    
    if ~isempty(dir(idxfile))
    % If IDXFILE already exists, try to stick to the previously defined indexing, in case the last
    % simulation crashed and [partial] results are included in PARTS or PARTS~ files...
        [fcn0,parts0,idx0] = load2vars(idxfile,{'fcn','parts','index'});
        if ~isequal(func2str(fcn0),func2str(fcn)) || numel(unique([idx0{:}])) ~= Nt
        % Check that fcn0 == fcn and Nt @ idx0 == Nt
        % (handles for identical functions differ, so have to check func2str) :/
            idx0 = {};
        end
    else
    % Alternatively, try to rebuild index from existing SUBS files
        parts0 = pickfile([POPT.prefix '_*' POPT.parts],Inf);        
        if ~isempty(parts0)
            idx0 = cellfun(@(s) textscan(s,['./' POPT.prefix '_%d' POPT.parts]),parts0);
            parts0(cellfun(@isempty,idx0)) = []; % refine search: PREFIX_N.PARTS
        end
        if ~isempty(parts0)
            [fcn0,idx0] = cellfun(@(f) load2vars(f,{'fcn','idx'}),parts0,'unif',0);

            % Keep SUBS files with fcn0 == fcn
            matching = cellfun(@(x) isequal(func2str(x),func2str(fcn)),fcn0);
            parts0 = parts0(matching);
            idx0 = idx0(matching);

            if numel(unique([idx0{:}])) ~= Nt, idx0 = {}; end
        else
           idx0 = {}; 
        end
    end
    
    % Set M from existing index/part-files, whenever possible
    if ~isempty(idx0)
        M0 = numel(parts0);
        if isfinite(M) && M0 ~= M
        % M0 is allright, but doesn't match provided M: override
            warning('Overriding M (%d) to match existing part-files (%d)',M,M0);
            M = M0;
        elseif ~isfinite(M)
        % Use M0 as default
            M = M0;
        end
    end
    
    % Put all non-separable arguments into a base structure
    prj0 = cell(size(args)); 
    prj0(~separable) = args(~separable);
    
    if opt.backup
        prj0(end+1:end+2) = {'backup','foo'}; % foo will be replaced (§)
        separable(end+1:end+2) = false;
    end
    
    % Get a memmory-bound number of parallel threads
    % (the whole of varargin is already loaded on memmory, so it shouldn't get worse)
    [B,Np] = getnworkers(P,'vars',whos('prj0'));
    P = min([M,P,B,Nt]);
    assert(isfinite(P) && mod(P,1) == 0 && P < MAX_P, '%d is not a reasonable number of threads',P);
        
    if ~isfinite(M), M = max(P,ceil(Nt/SimOptions.chunksize)); end
    M = min([M,Nt]);
    
    if M <= 1
        warning('runparallel:alone',['Using M = 1 (single-MATLAB instance), ',...
                 'check SimOptions.runparallel/maxthreads and memmory requirements']); 
    else
        fprintf('\tUsing %s on %s (%s)...\n',...
            nthings(M,'chunk'),nthings(P,'worker'),nthings(Np,'processor'));
    end

    if ~isempty(idx0) && numel(parts0) == M
    % Matching (existing) index
        fprintf('\t(Resuming from existing part-files and index: %s)\n',relativepath(idxfile));
        parts = parts0; idx = idx0;
        
        if isempty(dir(idxfile))
            savewithnames(idxfile,{'fcn','parts','index'},{fcn,parts,idx});
        end
        return;
    end
    
    % Non-matching (existing) index/part-files: 
    
    formatstr = ['%s_%0',num2str(floor(log10(M))+1),'d%s']; % %s_%0Xd%s, e.g. file001.ext
    parts = arrayfun(@(j) sprintf(formatstr,POPT.prefix,j,POPT.parts),1:M,'unif',0);
    logs = arrayfun(@(j) sprintf(formatstr,POPT.prefix,j,'.log'),1:M,'unif',0);
    backups = cellfun(@(s) [s POPT.backup],parts,'unif',0);
    
    % delete any conflicting files
    f = [idxfile;parts(:);logs(:);backups(:)];
    f = f(cellfun(@(f) ~isempty(dir(f)),f));
    if ~isempty(f)
       warning('Deleting inconsistent part-files and index: %s',idxfile);
       cellfun(@delete,f);
    end
    
    % Generate indices from scratch
    % NOTE: not divided in chunks, but alternated and shuffled! (for more uniform timing)
    shuffle = @(x) x(randperm(numel(x))');
    idx = arrayfun(@(j) shuffle(j:M:Nt),1:M,'unif',0);

    fprintf('\tDividing calculation into Workers...\n');
    
    % Generate and write part-files
    for j = 1:M
        subprj = prj0; % start with copy of base project
        for k = find(separable)
            % subprj{k} = splitarg(args{k},idx{j},Nt); 
            subprj{k} = filterstructure(args{k},idx{j},Nt); 
        end
        if opt.backup, subprj(end) = backups(j); end % (§)
        savewithnames(parts{j},{'fcn','args','idx','opt'},{fcn,subprj,idx{j},opts});
    end
    % Save index last, only if things went smoothly
    savewithnames(idxfile,{'fcn','parts','index'},{fcn,parts,idx});

    fprintf('\t...Data split into %d pieces.\n',M);
end
    
function N0 = checkseparable(X,N0)
    N = cell(numel(X),1);
    for j = 1:numel(N)
        if isstruct(X{j})
            N{j} = cellfun(@(f) size(X{j}.(f),1),fieldnames(X{j}));
            if all(N{j}==N{j}(1)), N{j} = N{j}(1); end
        elseif isnumeric(X{j}) || isdatetime(X{j}) || isa(X{j},'tabular')
            N{j} = size(X{j},1);
        elseif isa(X{j},'ShadingResults')
            N{j} = X{j}.Nt;
        else
            N{j} = 1;
        end
    end
    if isempty(N0) || ~isfinite(N0)
        n = cellfun(@numel,N);
        unclear = n > 1;
        if any(unclear)
            N0 = cat(1,N{~unclear});
            N(unclear) = cellfun(@(x) intersect(x,N0),N(unclear),'unif',0);
        end
        N0 = mode(cat(1,N{:}));
    end
    bad = ~cellfun(@(n) any(n == N0),N);
    assert(~any(bad),'%d Incoherent/unrecognized separable-arguments',nnz(bad));
end

% function x = splitarg(X,idx,N)
% 
%     if isstruct(X) || isa(X,'ShadingResults')
%         x = filterstructure(X,idx,N); % SHADINGRESULTS uses own overloaded method
%     elseif isnumeric(X)
%         s = size(X); 
%         s(1) = nnz(idx);
%         x = reshape(X(idx,:),s);
%     else
%         error('Cannot split argument of class %s',class(X))
%     end
% end

function savewithnames(filename,varnames,vars)
% Save cell-array of variables VARS, with names VARNAMES in file FILENAME
    for j = 1:numel(varnames)
       F.(varnames{j}) = vars{j};
    end
    if exist(filename,'file')
        save(filename,'-append','-struct','F');
    else
        save(filename,'-struct','F');
    end
end

function varargout = load2vars(filename,varnames)
% [a,b,c] = load2vars('filename',{'a','b','c'})

    warning_resetter = naptime('MATLAB:load:variableNotFound'); %#ok<NASGU>
    
    f = load(filename,varnames{:},'-mat');
    assert(all(isfield(f,varnames)),'Fields not found');
    varargout = cellfun(@(s) f.(s),varnames,'unif',0);
end

function test()
    
    lastwill = restoreSimOptions(); %#ok<NASGU>
    
    setSimOption('prjname','test');
    setSimOption('runparallel',true);
    setSimOption('maxthreads',31416);
       
    X = (1:20)';
    Y = X;
    M = 5;
    P = 2;
    
    Z =  dummy(X,false,Y);
    P = runparallel(@dummy,{X,true,Y},[1,3],M,P,'-backup','-simoptions');
    assert(isequal(Z,P));
    fprintf('Imagine that, it worked!\n');
end

function z = dummy(varargin)
% z = dummy(x,printstuff,y,['backup',X])
% z = dummy('backup',S)
% Multiply X,Y, and if printstuff, print some output, if SimOptions.prjname doesn't match the
% secret password, crash with a backup.

    [opt,varargin] = getpairedoptions(varargin,{'backup'},{''});
    switch numel(varargin)
    case 0
    % recover from crash
        x = opt.backup.x;
        y = opt.backup.y;
        printstuff = opt.backup.printing;
        opt.backup = opt.backup.opt.backup;
    case 3
    % fresh run
        [x,printstuff,y] = deal(varargin{:});
    otherwise
        error('Unrecognized arguments');
    end

    if printstuff, fprintf('\n\ndummy started on %s\n',datestr(now)); end

    if ~isequal(getSimOption('maxthreads'),31416)
        if ~isempty(opt.backup), save(opt.backup,'x','y','opt'); end
        error('DUMMY is running with wrong SimOptions');
    end

    z = x.*y;
    if printstuff
        for j = 1:numel(z)
            fprintf('%g x %g = %g\n',x(j),y(j),z(j));
            pause(1);
        end
    end
end