function [R,janitor] = mergeresults(varargin)    
% R = MERGERESULTS() - merge the partial results of an array of PACKEDFCN calls generated by 
%   RUNPARALLEL. Without arguments, MERGERESULTS will look for .par/*.idx or .par/*.parts files
%   in the current and/or project directories.
% 
% R = MERGERESULTS(IDX) - take the index-file-name IDX, with the full list of PART files.
% R = MERGERESULTS(SUBS) - use directly the PART files SUBS (cellstr of file names)
% 
% NOTE that in all cases, MERGERESULTS will attempt to find the underlying index-file, to make sure
% that all (and only) the partial files for a particular RUNPARALLEL function call are merged.
%
% [R,CLEANER] = MERGERESULTS(...) - When MERGERESULTS runs without problems, all of the relevant
%   files and jobs (i.e. SUBS, IDX, JOBS) are deleted after merging. To delay this deletion,
%   Function handle CLEANER can be passed around and called e.g. only after results have been 
%   saved to disk by the main calling function.
%
% See also: RUNPARALLEL, PACKEDFCN

    POPT = parallelconfig(); % Load file naming conventions and other options
    
    % Get a full list of the relevant jobs and/or subs-files
    [subfiles,idxfile] = findparts(varargin{:});
    backups = cellfun(@(f) [f POPT.backup],subfiles,'unif',0);
    logs = strrep(subfiles,POPT.parts,'.log');
    
    % Change to project directory, fall back to current directory on crash
    basedir = pwd(); lastwill = onCleanup(@() cd(basedir)); cd(fileparts(idxfile));
    
    warning_resetter = naptime(); %#ok<NASGU>
    warning('error','MATLAB:load:variableNotFound'); %#ok<CTPCT> undocumented :/

    % Read indices from files, and check that they make up a complete Nt set    
    idx = cellfun(@(f) load2vars(f,{'idx'}),subfiles,'unif',0);
    allidx = sort([idx{:}]);
    Nt = max(allidx);
    assert(isequal(allidx,1:Nt),'Missing or duplicate indices in provided jobs');
    
    % (try to) load the partial-results
    R = cell(numel(subfiles),1);
    finished = false(numel(subfiles),1);
    for j = 1:numel(subfiles)
        try
            R{j} = load2vars(subfiles{j},{'res'});
            finished(j) = true;
        catch ERR
            if ~strcmp(ERR.identifier,'MATLAB:load:variableNotFound'), rethrow(ERR); end
        end
    end

    if ~all(finished)
        % script = which('runparallel.sh');
        % assert(~isempty(script),'Failed to find runparallel.sh to complete jobs');
        % args = strjoin(cellfun(@(s) ['"' s '"'],subfiles,'unif',0),' ');
        % system([script ' ' args]);
        
        args = strjoin(cellfun(@(s) ['"' s '"'],subfiles(~finished),'unif',0),' ');
        error(['Not all files to be merged contain results. ',...
               'RUNPARALLEL might have terminated on error. \n',...
               'Run PACKEDFCN on each file, or from a console: \n\nrunparallel.sh %s'],args);
    end
    
    % Merge results into a single structure R
    R = mergestructures(R,idx,Nt);
    
    varname = @(x) inputname(1);  % (clearer var tracing, avoids NASGU/ASGLU warnings)
    clear(varname(lastwill));
    
    if nargout > 1
        janitor = @() cleanpartfiles(subfiles,idxfile,backups,logs);
    else
        cleanpartfiles(subfiles,idxfile,backups,logs);
    end
    
    function cleanpartfiles(partfiles,idxfile,backups,logs)
    % Clean-up if things will likely be ok
        warning_resetter_2 = naptime('MATLAB:DELETE:FileNotFound'); %#ok<NASGU>
        
        cd(fileparts(idxfile)); 
        cellfun(@delete,[partfiles(:);backups(:);logs(:);idxfile]);
        delete('parallel.log');
        delete('SimOption.xml');
        % cellfun(@delete,{'pathdef.m','runparallel.sh'});
              
        if numel(dir()) < 3, cd('..'); rmdir(fileparts(idxfile)); end % remove empty .par directory
    end
end

function [subs,idxfile] = findparts(arg)
% Parsing function for MERGERESULTS

    POPT = parallelconfig(); % Load file naming conventions and other options

    % R = MERGERESULTS() - look for .par/*.idx or .par/*.parts files
    if nargin < 1 || isempty(arg)
        basedir = pwd(); lastwill = onCleanup(@() cd(basedir));
        if ~isempty(dir(POPT.subdir))
            cd(POPT.subdir);
        elseif ~isempty(dir(fullfile(getSimOption('prjname'),POPT.subdir)))
            cd(fullfile(getSimoption('prjname'),POPT.subdir));
        end
        idxfile = pickfile(['*' POPT.index],Inf);
        if isempty(idxfile)
            subs = pickfile(['*' POPT.parts],Inf);
            assert(~isempty(subs),'Nothing to merge');
        else
           idxfile = idxfile{1}; 
        end
    elseif ischar(arg)
        idxfile = arg; subs = {};
    else
        subs = arg; idxfile = {};
    end
        
    if isempty(idxfile)
    % R = MERGERESULTS(SUBS) - do the same but from the input-output files SUBS (cellstr of file names)

        for j = 1:numel(subs)
        % Attempt to find index-file based on the names of sub-files

            % Guess prefix, removing digits and extension from end of subs-file name
            prefix = regexp(subs{j},['^(.*?)\d+' POPT.parts '$'],'tokens');
            while iscell(prefix) && ~isempty(prefix), prefix = prefix{:}; end
            if isempty(prefix), continue; end

            % Find files that mach the found prefix and POPT.index
            idxfile = arrayfun(@(f) fullfile(f.folder,f.name),dir(fileparts(prefix)),'unif',0); 
            idxfile(~contains(idxfile,prefix)) = [];
            idxfile(~endsWith(idxfile,POPT.index)) = [];

            if isscalar(idxfile)
               idxfile = idxfile{1};
               subs0 = load2vars(idxfile,{'parts'});
               if isempty(setdiff(subs,subs0)), break; else, idxfile = {}; end
            end
        end
        if isempty(idxfile)
            error('Failed to find valid index-file for (all) part-files: %s*',prefix);
        else
            [subs,idxfile] = findparts(idxfile); % recursive call, lands below
            return;
        end
    end
    assert(ischar(idxfile),'Unrecognized argument');
            
    % R = MERGERESULTS(IDX) - take the index-file-name IDX, with the full list of subs.
    idxfile = pickfile(idxfile,1,'Pick index file','fullpath',true);
    [subs,idx] = load2vars(idxfile,{'parts','index'});
    basedir = pwd(); lastwill = onCleanup(@() cd(basedir)); cd(fileparts(idxfile));

    % Compare index-file with indices in files  
    idxf = cellfun(@(f) load2vars(f,{'idx'}),subs,'unif',0);
    assert(isequal(idxf,idx),'Index-file idx does not match subs-file-indices');

    missing = cellfun(@(x) isempty(dir(x)),subs);
    assert(~any(missing),shortliststr(subs(missing),'Missing file'));
end

function varargout = load2vars(filename,varnames)
% [a,b,c] = load2vars('filename',{'a','b','c'})

    try
        f = load(filename,varnames{:},'-mat');
    catch ERR
        throwAsCaller(ERR);
    end
    %assert(all(isfield(f,varnames)),'Fields not found');
    varargout = cellfun(@(s) f.(s),varnames,'unif',0);
end


